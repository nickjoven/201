"""
Stribeck friction curve and its mapping to MOND interpolation.

The Stribeck curve describes how the friction coefficient mu_f depends on
sliding velocity v_s.  The velocity-weakening branch below the critical
velocity maps structurally onto the MOND interpolating function with
x = a / a0.

Reference: joven_unifying_framework.md §2, §6
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def stribeck_curve(v: ArrayLike, *,
                   mu_s: float = 1.0,
                   mu_k: float = 0.4,
                   v_s: float = 1.0,
                   delta: float = 2.0) -> np.ndarray:
    r"""Classical Stribeck friction coefficient.

    .. math::
        \mu_f(v) = \mu_k + (\mu_s - \mu_k)\,\exp\!\bigl[-(v / v_s)^\delta\bigr]

    Parameters
    ----------
    v : array_like
        Sliding velocity (or acceleration proxy).
    mu_s : float
        Static friction coefficient (peak at v = 0).
    mu_k : float
        Kinetic friction coefficient (plateau at high v).
    v_s : float
        Stribeck velocity (transition scale, analogous to a0).
    delta : float
        Stribeck exponent controlling sharpness.

    Returns
    -------
    mu_f : ndarray
        Friction coefficient at each velocity.
    """
    v = np.asarray(v, dtype=float)
    return mu_k + (mu_s - mu_k) * np.exp(-np.abs(v / v_s) ** delta)


def stribeck_to_mond(v: ArrayLike, *,
                     v_s: float = 1.0,
                     delta: float = 2.0) -> np.ndarray:
    r"""Map Stribeck velocity to MOND-like coupling excess.

    Returns the *excess coupling* below the threshold:

    .. math::
        \xi(v) = \exp\!\bigl[-(v / v_s)^\delta\bigr]

    This is structurally parallel to 1 - mu(x) in the MOND framework.
    When xi -> 0  (v >> v_s) the system is Newtonian (kinetic plateau).
    When xi -> 1  (v << v_s) the system shows maximum excess coupling
    (deep-MOND / static-friction regime).
    """
    v = np.asarray(v, dtype=float)
    return np.exp(-np.abs(v / v_s) ** delta)


def mond_from_stribeck(x: ArrayLike, *, delta: float = 2.0) -> np.ndarray:
    r"""MOND interpolating function derived from Stribeck structure.

    .. math::
        \mu(x) = 1 - \exp(-x^\delta)

    With delta = 0.5 this recovers the RAR form  mu = 1 - exp(-sqrt(x)).
    """
    x = np.asarray(x, dtype=float)
    return 1.0 - np.exp(-np.abs(x) ** delta)
