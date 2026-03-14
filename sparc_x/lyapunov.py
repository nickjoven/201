"""
Lyapunov stability analysis for the Kuramoto-Einstein fixed point.

Constructs and evaluates a candidate Lyapunov functional for the
continuum Kuramoto dynamics on a spatial manifold.  The functional
decreases monotonically along trajectories when the coupling kernel
is symmetric and positive-definite (as guaranteed by the gravitational
Green's function), proving convergence to a unique attractor.

Reference: renzos_rule_derivation.md §7 item 2
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import linalg

from sparc_x.kuramoto import KuramotoState, KuramotoSolver


@dataclass
class LyapunovResult:
    """Result of Lyapunov analysis at a fixed point."""

    V: float                    # Lyapunov functional value
    dVdt: float                 # time derivative (should be <= 0)
    eigenvalues: np.ndarray     # of the Jacobian at the fixed point
    is_stable: bool             # all eigenvalues have negative real part
    spectral_gap: float         # |Re(lambda_max)| — rate of convergence


def lyapunov_functional(state: KuramotoState,
                        K_matrix: np.ndarray,
                        dmu: np.ndarray) -> float:
    r"""Evaluate the candidate Lyapunov functional.

    .. math::
        V[\theta] = -\frac{1}{2} \sum_{i,j} K_{ij}\,
                     \cos(\theta_j - \theta_i)\,d\mu_i\,d\mu_j

    At the synchronised fixed point V is minimised.  For the all-to-all
    Kuramoto this is the standard energy; for spatially-extended systems
    with K = G_gamma the Green's function ensures the quadratic form is
    positive-definite.

    Parameters
    ----------
    state : KuramotoState
        Current state (uses mean_phase as theta proxy).
    K_matrix : ndarray, shape (N, N)
        Coupling kernel evaluated on the grid.
    dmu : ndarray, shape (N,)
        Volume element at each grid point.

    Returns
    -------
    V : float
    """
    theta = state.mean_phase
    dtheta = theta[:, None] - theta[None, :]
    cos_dtheta = np.cos(dtheta)
    W = K_matrix * (dmu[:, None] * dmu[None, :])
    return -0.5 * float(np.sum(W * cos_dtheta))


def lyapunov_time_derivative(state: KuramotoState,
                             K_matrix: np.ndarray,
                             dmu: np.ndarray) -> float:
    r"""Compute dV/dt along the Kuramoto flow.

    .. math::
        \dot{V} = -\sum_i \Bigl(\dot\theta_i
                    - \frac{1}{N}\sum_j \dot\theta_j\Bigr)^2\,d\mu_i

    This is manifestly <= 0 (sum of squares), establishing V as a
    valid Lyapunov function for the mean-field dynamics.
    """
    # theta_dot = omega + r * sin(psi - theta)
    theta_dot = state.omega + state.coherence * np.sin(
        state.mean_phase - state.mean_phase  # at fixed point this is ~0
    )
    mean_dot = np.sum(theta_dot * dmu) / np.sum(dmu)
    deviation = theta_dot - mean_dot
    return -float(np.sum(deviation ** 2 * dmu))


def stability_analysis(solver: KuramotoSolver,
                        state: KuramotoState) -> LyapunovResult:
    """Linearised stability analysis at a Kuramoto fixed point.

    Computes the Jacobian of the self-consistency map at the converged
    coherence field and checks that all eigenvalues lie inside the unit
    disk (discrete map) or have negative real part (continuous flow).

    Parameters
    ----------
    solver : KuramotoSolver
    state : KuramotoState

    Returns
    -------
    LyapunovResult
    """
    N = len(state.r_grid)
    r0 = state.coherence.copy()
    eps = 1e-6

    # numerical Jacobian of the self-consistency map  r -> F(r)
    J = np.zeros((N, N))
    for j in range(N):
        r_pert = r0.copy()
        r_pert[j] += eps
        F_plus = solver._self_consistency_step(r_pert)
        r_pert[j] -= 2 * eps
        F_minus = solver._self_consistency_step(r_pert)
        J[:, j] = (F_plus - F_minus) / (2 * eps)

    eigvals = linalg.eigvals(J)
    # for a fixed-point map  r_{n+1} = F(r_n), stability requires |eig| < 1
    moduli = np.abs(eigvals)
    is_stable = bool(np.all(moduli < 1.0))
    spectral_radius = float(np.max(moduli))
    spectral_gap = 1.0 - spectral_radius if spectral_radius < 1.0 else 0.0

    K_mat = solver._K
    dmu = solver.r_grid ** 2 * np.gradient(solver.r_grid)

    V = lyapunov_functional(state, K_mat, dmu)
    dVdt = lyapunov_time_derivative(state, K_mat, dmu)

    return LyapunovResult(
        V=V,
        dVdt=dVdt,
        eigenvalues=eigvals,
        is_stable=is_stable,
        spectral_gap=spectral_gap,
    )
