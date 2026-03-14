"""
Main Calculator class — unified entry point for the SPARC-X framework.

Mirrors the design of the SPARC-X-API calculator: attach a Calculator to
a galaxy profile (analogous to an ASE Atoms object), then call methods to
compute energies, forces (accelerations), and derived quantities.

Example
-------
>>> from sparc_x import Calculator
>>> from sparc_x.profiles import exponential_disk, GalaxyProfile
>>> import numpy as np
>>>
>>> r = np.linspace(0.5, 30, 60)
>>> v_disk = exponential_disk(r, M_disk=5e10, R_d=3.0)
>>> profile = GalaxyProfile(
...     name="synthetic", r_kpc=r, v_obs=v_disk,
...     v_gas=np.zeros_like(r), v_disk=v_disk,
...     v_bul=np.zeros_like(r), e_v_obs=np.ones_like(r),
... )
>>> calc = Calculator(profile, interpolation="rar")
>>> v_pred = calc.get_rotation_curve()
>>> rho_dm = calc.get_dark_matter_density()
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from sparc_x.constants import A0, G, KPC, KMS
from sparc_x.profiles import GalaxyProfile, predict_mond_rc, sparc_prefactor
from sparc_x.mond import a_obs as mond_a_obs, sigma_static
from sparc_x.kuramoto import KuramotoSolver, KuramotoState, FixedPointResult
from sparc_x.lyapunov import stability_analysis, LyapunovResult


@dataclass
class CalculatorResults:
    """Cache of computed quantities."""

    v_mond: np.ndarray | None = None
    v_kuramoto: np.ndarray | None = None
    rho_dark: np.ndarray | None = None
    kuramoto_state: KuramotoState | None = None
    lyapunov: LyapunovResult | None = None
    prefactors: np.ndarray | None = None


class Calculator:
    """Unified calculator for the synchronisation-gravity framework.

    Parameters
    ----------
    profile : GalaxyProfile
        Galaxy profile to analyse.
    a0 : float
        MOND acceleration scale [m s^-2].
    interpolation : {"simple", "standard", "rar"}
        Which MOND interpolating function to use.
    K_coupling : float
        Kuramoto coupling strength (in units of a0).
    damping : float
        Under-relaxation for the Kuramoto fixed-point iteration.
    """

    def __init__(
        self,
        profile: GalaxyProfile,
        *,
        a0: float = A0,
        interpolation: str = "rar",
        K_coupling: float = 1.0,
        damping: float = 0.5,
    ) -> None:
        self.profile = profile
        self.a0 = a0
        self.interpolation = interpolation
        self.K_coupling = K_coupling
        self.damping = damping
        self._results = CalculatorResults()

    # ------------------------------------------------------------------
    # MOND-level predictions
    # ------------------------------------------------------------------

    def get_rotation_curve(self, *, method: str = "mond") -> np.ndarray:
        """Predicted circular velocity [km/s].

        Parameters
        ----------
        method : {"mond", "kuramoto"}
            "mond"     — algebraic MOND boost.
            "kuramoto" — full Kuramoto fixed-point solve.
        """
        if method == "mond":
            if self._results.v_mond is None:
                self._results.v_mond = predict_mond_rc(
                    self.profile, a0=self.a0,
                    interpolation=self.interpolation,
                )
            return self._results.v_mond

        if method == "kuramoto":
            if self._results.v_kuramoto is None:
                self._solve_kuramoto()
            return self._results.v_kuramoto  # type: ignore[return-value]

        raise ValueError(f"Unknown method {method!r}")

    def get_accelerations(self) -> dict[str, np.ndarray]:
        """Return a_N, a_obs, and a_predicted [m s^-2]."""
        a_pred = mond_a_obs(
            self.profile.a_N, a0=self.a0,
            interpolation=self.interpolation,
        )
        return {
            "a_N": self.profile.a_N,
            "a_obs": self.profile.a_obs,
            "a_pred": a_pred,
        }

    def get_dark_matter_density(self) -> np.ndarray:
        r"""Synchronisation deficit  rho_dark(r)  [kg m^-3].

        From the Hamiltonian constraint:
            rho_dark = (1 / 16 pi G) (R + K^2 - K_ij K^ij) - rho_bary

        Approximated in the Newtonian limit as:
            rho_dark ~ (a_obs - a_N) / (4 pi G r)
        """
        if self._results.rho_dark is None:
            a_total = mond_a_obs(
                self.profile.a_N, a0=self.a0,
                interpolation=self.interpolation,
            )
            deficit = np.maximum(a_total - self.profile.a_N, 0.0)
            r_m = self.profile.r_kpc * KPC
            r_m = np.where(r_m > 0, r_m, 1e-10)
            self._results.rho_dark = deficit / (4.0 * np.pi * G * r_m)
        return self._results.rho_dark

    # ------------------------------------------------------------------
    # Kuramoto-level predictions
    # ------------------------------------------------------------------

    def _solve_kuramoto(self) -> None:
        """Run the full Kuramoto fixed-point solver."""
        solver = KuramotoSolver(
            r_grid=self.profile.r_kpc,
            omega=self.profile.omega,
            K_coupling=self.K_coupling,
            damping=self.damping,
        )
        state = solver.solve()
        self._results.kuramoto_state = state

        # Convert coherence to predicted acceleration and velocity
        # a ~ c^2 grad(r)  in the Newtonian limit; here we use
        # a_pred = a_N / mu(x)  with mu approximated by r(x).
        coherence = np.maximum(state.coherence, 1e-10)
        a_pred = self.profile.a_N / coherence
        r_m = self.profile.r_kpc * KPC
        r_m = np.where(r_m > 0, r_m, 1e-10)
        self._results.v_kuramoto = np.sqrt(a_pred * r_m) / KMS

        # dark matter from Kuramoto
        deficit = np.maximum(a_pred - self.profile.a_N, 0.0)
        self._results.rho_dark = deficit / (4.0 * np.pi * G * r_m)

        # store solver for stability analysis
        self._solver = solver

    def get_kuramoto_state(self) -> KuramotoState:
        """Return the converged Kuramoto state."""
        if self._results.kuramoto_state is None:
            self._solve_kuramoto()
        return self._results.kuramoto_state  # type: ignore[return-value]

    def get_lyapunov(self) -> LyapunovResult:
        """Lyapunov stability analysis at the Kuramoto fixed point."""
        if self._results.lyapunov is None:
            if self._results.kuramoto_state is None:
                self._solve_kuramoto()
            self._results.lyapunov = stability_analysis(
                self._solver, self._results.kuramoto_state  # type: ignore
            )
        return self._results.lyapunov

    # ------------------------------------------------------------------
    # Prefactor / audit helpers
    # ------------------------------------------------------------------

    def get_prefactors(self) -> np.ndarray:
        r"""SPARC prefactor  eta = a_obs / sqrt(a_N * a0).

        eta -> 1 in deep-MOND.  Deviations reveal interpolation structure.
        """
        if self._results.prefactors is None:
            self._results.prefactors = sparc_prefactor(self.profile)
        return self._results.prefactors

    def get_sigma(self, v: np.ndarray | None = None) -> np.ndarray:
        r"""Surface-density function  Sigma(v).

        If *v* is None, uses the observed rotation velocities.
        """
        if v is None:
            v = self.profile.v_obs * KMS  # convert to m/s
        return sigma_static(v, a0=self.a0)

    def get_regime_mask(self) -> dict[str, np.ndarray]:
        """Boolean masks for Newtonian (x > 1) and deep-MOND (x < 1) radii."""
        x = self.profile.x
        return {
            "newtonian": x > 1.0,
            "deep_mond": x < 1.0,
            "transition": (x > 0.3) & (x < 3.0),
        }

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def summary(self) -> dict[str, Any]:
        """Return a summary dict suitable for inspection or serialisation."""
        masks = self.get_regime_mask()
        return {
            "name": self.profile.name,
            "n_points": len(self.profile.r_kpc),
            "r_range_kpc": (float(self.profile.r_kpc[0]),
                            float(self.profile.r_kpc[-1])),
            "v_obs_range_kms": (float(self.profile.v_obs.min()),
                                float(self.profile.v_obs.max())),
            "frac_newtonian": float(masks["newtonian"].mean()),
            "frac_deep_mond": float(masks["deep_mond"].mean()),
            "a0": self.a0,
            "interpolation": self.interpolation,
        }
