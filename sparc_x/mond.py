"""
MOND interpolation functions and the Stribeck-MOND correspondence.

The interpolating function mu(x) bridges Newtonian (x >> 1) and deep-MOND
(x << 1) regimes, where x = a / a0.  Each variant has a direct reading as
a Stribeck friction curve with velocity-weakening below the threshold.

Reference: joven_unifying_framework.md §2.1, §6
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from sparc_x.constants import A0


# ---------------------------------------------------------------------------
# Interpolating functions  mu(x)  where x = a_N / a0
# ---------------------------------------------------------------------------

def mu_simple(x: ArrayLike) -> np.ndarray:
    """Simple interpolating function: mu = x / (1 + x)."""
    x = np.asarray(x, dtype=float)
    return x / (1.0 + x)


def mu_standard(x: ArrayLike) -> np.ndarray:
    """Standard interpolating function: mu = x / sqrt(1 + x^2)."""
    x = np.asarray(x, dtype=float)
    return x / np.sqrt(1.0 + x * x)


def mu_rar(x: ArrayLike) -> np.ndarray:
    """RAR (McGaugh+2016) interpolating function: mu = 1 - exp(-sqrt(x))."""
    x = np.asarray(x, dtype=float)
    return 1.0 - np.exp(-np.sqrt(np.maximum(x, 0.0)))


# ---------------------------------------------------------------------------
# Inverse:  nu(y)  such that  a_obs = nu(a_N / a0) * a_N
# ---------------------------------------------------------------------------

def nu_simple(y: ArrayLike) -> np.ndarray:
    """Inverse (boost) for the simple interpolation."""
    y = np.asarray(y, dtype=float)
    return 0.5 + np.sqrt(0.25 + 1.0 / np.maximum(y, 1e-30))


def nu_standard(y: ArrayLike) -> np.ndarray:
    """Inverse (boost) for the standard interpolation."""
    y = np.asarray(y, dtype=float)
    return 1.0 / np.sqrt(1.0 - np.exp(-np.sqrt(np.maximum(y, 0.0))))


def nu_rar(y: ArrayLike) -> np.ndarray:
    """Inverse (boost) for the RAR interpolation."""
    y = np.asarray(y, dtype=float)
    e = np.exp(-np.sqrt(np.maximum(y, 0.0)))
    return 1.0 / (1.0 - e)


# ---------------------------------------------------------------------------
# Observed acceleration from Newtonian
# ---------------------------------------------------------------------------

def a_obs(a_N: ArrayLike, *, a0: float = A0,
          interpolation: str = "rar") -> np.ndarray:
    """Return observed (total) acceleration from Newtonian acceleration.

    Parameters
    ----------
    a_N : array_like
        Newtonian (baryonic) acceleration [m s^-2].
    a0 : float
        MOND acceleration scale.
    interpolation : {"simple", "standard", "rar"}
        Which interpolating function to use.

    Returns
    -------
    a_total : ndarray
        Predicted total acceleration [m s^-2].
    """
    _nu = {"simple": nu_simple, "standard": nu_standard, "rar": nu_rar}
    if interpolation not in _nu:
        raise ValueError(f"Unknown interpolation {interpolation!r}")
    y = np.asarray(a_N, dtype=float) / a0
    return np.asarray(a_N, dtype=float) * _nu[interpolation](y)


# ---------------------------------------------------------------------------
# Surface-density function  Sigma(v)  (§7.1 audit item 1)
# ---------------------------------------------------------------------------

def sigma_static(v: ArrayLike, *, Sigma_0: float = 1.0,
                 a0: float = A0) -> np.ndarray:
    r"""Static form of the surface-density function Sigma(v).

    In the synchronization framework the effective surface density that
    enters the Tully-Fisher-like relation is:

        Sigma(v) = v^4 / (G * a0)            (deep-MOND, x << 1)
        Sigma(v) = v^2 / (G * r)             (Newtonian,  x >> 1)

    This helper returns the deep-MOND form normalised by *Sigma_0*.
    """
    from sparc_x.constants import G as _G
    v = np.asarray(v, dtype=float)
    return Sigma_0 * v ** 4 / (_G * a0)
