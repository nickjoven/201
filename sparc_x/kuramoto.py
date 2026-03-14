"""
Continuum Kuramoto model on a radial manifold.

Implements the self-consistency equation for the order parameter r(x)
and the coupled fixed-point iteration that determines the synchronization
deficit (dark-matter density) from the baryonic natural-frequency profile.

Reference: kuramoto_einstein_mapping.md §1-§4
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy import integrate, interpolate


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class KuramotoState:
    """Snapshot of the Kuramoto fields on a 1-D radial grid."""

    r_grid: np.ndarray          # radial positions  [arbitrary units]
    omega: np.ndarray           # natural frequency  omega(r) = sqrt(4 pi G rho_bary)
    coherence: np.ndarray       # order parameter  r(x) in [0, 1]
    mean_phase: np.ndarray      # psi(x) — mean phase field
    K_eff: np.ndarray           # effective coupling at each radius
    n_iter: int = 0             # iterations to converge


@dataclass
class FixedPointResult:
    """Result of the self-consistency fixed-point iteration."""

    state: KuramotoState
    converged: bool
    residual: float
    rho_dark: np.ndarray        # synchronization deficit -> dark matter density
    v_circ: np.ndarray          # predicted circular velocity


# ---------------------------------------------------------------------------
# Coupling kernel
# ---------------------------------------------------------------------------

def gravitational_kernel(r1: np.ndarray, r2: np.ndarray) -> np.ndarray:
    """Green's function of the radial Laplacian (Newtonian potential).

    K(r1, r2) ~ 1 / |r1 - r2|  for r1 != r2, regularised at r1 == r2.
    Returns the full 2-D kernel matrix  K[i, j].
    """
    R1, R2 = np.meshgrid(r1, r2, indexing="ij")
    sep = np.abs(R1 - R2)
    eps = 0.5 * np.min(np.diff(r1)) if len(r1) > 1 else 1.0
    return 1.0 / np.sqrt(sep ** 2 + eps ** 2)


def enclosed_mass_kernel(r_grid: np.ndarray,
                         rho: np.ndarray) -> np.ndarray:
    """Cumulative (enclosed) mass profile from a density profile.

    M(<r) = 4 pi int_0^r rho(r') r'^2 dr'
    """
    integrand = 4.0 * np.pi * rho * r_grid ** 2
    return integrate.cumulative_trapezoid(integrand, r_grid, initial=0.0)


# ---------------------------------------------------------------------------
# Self-consistency solver
# ---------------------------------------------------------------------------

class KuramotoSolver:
    """Solve the Kuramoto self-consistency equation on a radial grid.

    At steady state the coherence field satisfies:

        r(x) = F[omega(x), K[gamma[r]]]

    We iterate this map starting from an initial guess until convergence.

    Parameters
    ----------
    r_grid : array, shape (N,)
        Radial positions [kpc or natural units].
    omega : array, shape (N,)
        Natural frequency at each radius (sqrt(4 pi G rho_bary)).
    K_coupling : float
        Global coupling strength (related to a0).
    kernel : callable, optional
        Coupling kernel K(r1, r2).  Defaults to gravitational_kernel.
    damping : float
        Under-relaxation factor in (0, 1].  Lower = more stable.
    """

    def __init__(
        self,
        r_grid: np.ndarray,
        omega: np.ndarray,
        K_coupling: float = 1.0,
        kernel: Callable | None = None,
        damping: float = 0.5,
    ) -> None:
        self.r_grid = np.asarray(r_grid, dtype=float)
        self.omega = np.asarray(omega, dtype=float)
        self.K_coupling = K_coupling
        self.damping = damping

        if kernel is None:
            self._K = gravitational_kernel(self.r_grid, self.r_grid)
        else:
            self._K = kernel(self.r_grid, self.r_grid)

        self._dr = np.gradient(self.r_grid)
        # volume element  dmu = r^2 dr  (spherical)
        self._dmu = self.r_grid ** 2 * self._dr

    # ------------------------------------------------------------------
    def _mean_field(self, coherence: np.ndarray) -> np.ndarray:
        """Compute effective coupling  K_eff(x) = K * int K(x,x') r(x') dmu."""
        weighted = self._K @ (coherence * self._dmu)
        norm = self._K @ self._dmu
        norm = np.where(norm > 0, norm, 1.0)
        return self.K_coupling * weighted / norm

    # ------------------------------------------------------------------
    def _self_consistency_step(
        self, coherence: np.ndarray
    ) -> np.ndarray:
        """One Kuramoto self-consistency update.

        For the standard Kuramoto model with Lorentzian g(omega):
            r_new = sqrt(1 - K_c / K_eff)   when K_eff > K_c
            r_new = 0                        otherwise
        where K_c = 2 / (pi * g(0))  is the critical coupling.
        We use a smooth approximation amenable to non-uniform omega(x).
        """
        K_eff = self._mean_field(coherence)
        # Local critical coupling ~ 2 * omega (width of frequency dist.)
        K_c_local = 2.0 * self.omega
        ratio = K_eff / np.maximum(K_c_local, 1e-30)
        # Square-root onset  r ~ sqrt(K/Kc - 1)  for K > Kc
        safe_ratio = np.maximum(ratio, 1.0)
        r_new = np.where(ratio > 1.0,
                         np.sqrt(1.0 - 1.0 / safe_ratio),
                         0.0)
        return np.clip(r_new, 0.0, 1.0)

    # ------------------------------------------------------------------
    def solve(
        self,
        r0: np.ndarray | None = None,
        max_iter: int = 500,
        tol: float = 1e-8,
    ) -> KuramotoState:
        """Iterate the self-consistency map to convergence.

        Parameters
        ----------
        r0 : array, optional
            Initial guess for coherence.  Defaults to 0.5 everywhere.
        max_iter : int
            Maximum number of iterations.
        tol : float
            Convergence tolerance on max |r_new - r_old|.

        Returns
        -------
        KuramotoState
        """
        r = np.full_like(self.r_grid, 0.5) if r0 is None else np.array(r0)
        alpha = self.damping

        for n in range(1, max_iter + 1):
            r_new = self._self_consistency_step(r)
            r_update = alpha * r_new + (1.0 - alpha) * r
            residual = np.max(np.abs(r_update - r))
            r = r_update
            if residual < tol:
                break

        K_eff = self._mean_field(r)
        psi = np.zeros_like(r)  # mean phase (trivial in steady state)

        return KuramotoState(
            r_grid=self.r_grid,
            omega=self.omega,
            coherence=r,
            mean_phase=psi,
            K_eff=K_eff,
            n_iter=n,
        )


# ---------------------------------------------------------------------------
# Finite-N fluctuation analysis  (audit item 3 — deep-MOND prefactors)
# ---------------------------------------------------------------------------

def finite_N_fluctuations(N: int, K: float, omega_spread: float,
                          n_samples: int = 2000,
                          rng: np.random.Generator | None = None
                          ) -> dict:
    """Monte-Carlo estimate of order-parameter fluctuations for finite N.

    In the thermodynamic limit (N -> inf) the Kuramoto transition is sharp.
    For finite N the order parameter fluctuates around the mean-field value
    with variance ~ 1/N.  In the deep-MOND regime (K ~ K_c) these
    fluctuations set the prefactor of the sqrt(a_N * a0) scaling.

    Parameters
    ----------
    N : int
        Number of discrete oscillators.
    K : float
        Coupling strength.
    omega_spread : float
        Half-width of the uniform frequency distribution.
    n_samples : int
        Monte-Carlo realisations.

    Returns
    -------
    dict with keys: r_mean, r_std, r_samples, K_c, x (= K/K_c).
    """
    if rng is None:
        rng = np.random.default_rng(42)

    K_c = 2.0 * omega_spread / np.pi  # critical coupling (Lorentzian approx)
    r_samples = np.empty(n_samples)

    for s in range(n_samples):
        omegas = rng.uniform(-omega_spread, omega_spread, size=N)
        # evolve phases with simple Euler for a few relaxation times
        theta = rng.uniform(0, 2 * np.pi, size=N)
        dt = 0.05 / max(K, 1e-6)
        for _ in range(300):
            # mean-field coupling
            z = np.mean(np.exp(1j * theta))
            r_inst = np.abs(z)
            psi_inst = np.angle(z)
            dtheta = omegas + K * r_inst * np.sin(psi_inst - theta)
            theta += dt * dtheta
            theta %= 2 * np.pi
        z = np.mean(np.exp(1j * theta))
        r_samples[s] = np.abs(z)

    return {
        "r_mean": float(np.mean(r_samples)),
        "r_std": float(np.std(r_samples)),
        "r_samples": r_samples,
        "K_c": K_c,
        "x": K / K_c if K_c > 0 else np.inf,
    }
