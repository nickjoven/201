"""
SPARC galaxy profile handling and rotation-curve utilities.

Loads observed baryonic profiles (surface brightness -> mass density),
computes Newtonian accelerations, and converts to the framework's
natural-frequency field omega(r).

Reference: renzos_rule_derivation.md §7 item 1 & 4
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import interpolate, integrate

from sparc_x.constants import G, A0, M_SUN, KPC, KMS


# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------

@dataclass
class GalaxyProfile:
    """Observed SPARC galaxy profile."""

    name: str
    r_kpc: np.ndarray           # radial positions  [kpc]
    v_obs: np.ndarray           # observed rotation velocity  [km/s]
    v_gas: np.ndarray           # gas contribution  [km/s]
    v_disk: np.ndarray          # stellar disk contribution  [km/s]
    v_bul: np.ndarray           # bulge contribution  [km/s]
    e_v_obs: np.ndarray         # error on v_obs  [km/s]
    distance: float = 0.0       # distance  [Mpc]
    luminosity: float = 0.0     # total luminosity  [L_sun]
    inc: float = 0.0            # inclination  [deg]
    quality: int = 0            # quality flag (1 = best)

    # derived ---
    @property
    def a_N(self) -> np.ndarray:
        """Newtonian (baryonic) centripetal acceleration [m s^-2]."""
        v_bar_sq = self.v_gas ** 2 + self.v_disk ** 2 + self.v_bul ** 2
        r_m = self.r_kpc * KPC
        r_m = np.where(r_m > 0, r_m, 1e-10)
        return np.abs(v_bar_sq) * KMS ** 2 / r_m

    @property
    def a_obs(self) -> np.ndarray:
        """Observed centripetal acceleration [m s^-2]."""
        r_m = self.r_kpc * KPC
        r_m = np.where(r_m > 0, r_m, 1e-10)
        return (self.v_obs * KMS) ** 2 / r_m

    @property
    def x(self) -> np.ndarray:
        """Dimensionless acceleration  x = a_N / a0."""
        return self.a_N / A0

    @property
    def omega(self) -> np.ndarray:
        """Natural frequency field  omega(r) = sqrt(4 pi G rho_bary).

        Approximated from the Newtonian acceleration profile via
        rho ~ a_N / (4 pi G r).
        """
        r_m = self.r_kpc * KPC
        r_m = np.where(r_m > 0, r_m, 1e-10)
        rho_approx = self.a_N / (4.0 * np.pi * G * r_m)
        rho_approx = np.maximum(rho_approx, 0.0)
        return np.sqrt(4.0 * np.pi * G * rho_approx)


# ---------------------------------------------------------------------------
# I/O — load a single SPARC rotmod file
# ---------------------------------------------------------------------------

def load_rotmod(path: str | Path) -> GalaxyProfile:
    """Load a SPARC *_rotmod.dat file.

    Expected columns (whitespace-separated):
        Rad  Vobs  errV  Vgas  Vdisk  Vbul  SBdisk  SBbul

    Only the first six columns are required.
    """
    path = Path(path)
    data = np.loadtxt(path, comments="#")
    ncol = data.shape[1]

    r   = data[:, 0]
    vo  = data[:, 1]
    ev  = data[:, 2]
    vg  = data[:, 3]
    vd  = data[:, 4]
    vb  = data[:, 5] if ncol > 5 else np.zeros_like(r)

    name = path.stem.replace("_rotmod", "")
    return GalaxyProfile(
        name=name, r_kpc=r, v_obs=vo, v_gas=vg,
        v_disk=vd, v_bul=vb, e_v_obs=ev,
    )


def load_sparc_table(path: str | Path) -> dict[str, GalaxyProfile]:
    """Load all *_rotmod.dat files in a directory.

    Returns a dict  {galaxy_name: GalaxyProfile}.
    """
    root = Path(path)
    profiles: dict[str, GalaxyProfile] = {}
    for f in sorted(root.glob("*_rotmod.dat")):
        gp = load_rotmod(f)
        profiles[gp.name] = gp
    return profiles


# ---------------------------------------------------------------------------
# Synthetic profile builders
# ---------------------------------------------------------------------------

def exponential_disk(r_kpc: np.ndarray, *,
                     M_disk: float = 5e10,
                     R_d: float = 3.0) -> np.ndarray:
    """Circular velocity from a thin exponential disk (Freeman 1970).

    Parameters
    ----------
    r_kpc : array
        Radial positions [kpc].
    M_disk : float
        Total disk mass [solar masses].
    R_d : float
        Disk scale length [kpc].

    Returns
    -------
    v_disk : ndarray
        Disk contribution to circular velocity [km/s].
    """
    from scipy.special import i0, i1, k0, k1

    y = 0.5 * r_kpc / R_d
    y = np.maximum(y, 1e-10)
    # Freeman formula:  v^2 = 4 pi G Sigma_0 R_d y^2 [I0 K0 - I1 K1]
    factor = i0(y) * k0(y) - i1(y) * k1(y)
    M_kg = M_disk * M_SUN
    R_m = R_d * KPC
    v2 = 4.0 * np.pi * G * (M_kg / (2.0 * np.pi * R_m)) * (y ** 2) * factor
    # convert to km/s
    return np.sqrt(np.maximum(v2, 0.0)) / KMS


def point_mass(r_kpc: np.ndarray, M: float) -> np.ndarray:
    """Keplerian velocity from a point mass [km/s]."""
    r_m = r_kpc * KPC
    r_m = np.where(r_m > 0, r_m, 1e-10)
    v = np.sqrt(G * M * M_SUN / r_m)
    return v / KMS


# ---------------------------------------------------------------------------
# Rotation-curve prediction using MOND
# ---------------------------------------------------------------------------

def predict_mond_rc(profile: GalaxyProfile, *,
                    a0: float = A0,
                    interpolation: str = "rar") -> np.ndarray:
    """Predict the total rotation curve from the baryonic profile using MOND.

    Returns v_pred [km/s].
    """
    from sparc_x.mond import a_obs as _a_obs

    a_total = _a_obs(profile.a_N, a0=a0, interpolation=interpolation)
    r_m = profile.r_kpc * KPC
    r_m = np.where(r_m > 0, r_m, 1e-10)
    v_m = np.sqrt(a_total * r_m)
    return v_m / KMS


# ---------------------------------------------------------------------------
# Prefactor analysis  (audit item 3)
# ---------------------------------------------------------------------------

def sparc_prefactor(profile: GalaxyProfile) -> np.ndarray:
    r"""Compute the SPARC prefactor  eta(r) = a_obs / sqrt(a_N * a0).

    In the deep-MOND regime (x << 1)  eta -> 1 exactly.
    In the Newtonian regime (x >> 1)  eta -> sqrt(x).
    Deviations from these limits reveal the interpolation structure.
    """
    a_n = profile.a_N
    a_o = profile.a_obs
    safe = np.maximum(a_n * A0, 1e-60)
    return a_o / np.sqrt(safe)
