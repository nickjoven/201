"""Tests for the Kuramoto solver."""

import numpy as np
import pytest

from sparc_x.kuramoto import (
    KuramotoSolver, gravitational_kernel,
    enclosed_mass_kernel, finite_N_fluctuations,
)


class TestGravitationalKernel:
    def test_symmetric(self):
        r = np.linspace(1, 10, 20)
        K = gravitational_kernel(r, r)
        np.testing.assert_allclose(K, K.T, atol=1e-15)

    def test_positive(self):
        r = np.linspace(1, 10, 20)
        K = gravitational_kernel(r, r)
        assert np.all(K > 0)


class TestEnclosedMass:
    def test_uniform_density(self):
        r = np.linspace(0.01, 1.0, 200)
        rho = np.ones_like(r)  # uniform density = 1
        M = enclosed_mass_kernel(r, rho)
        # M(<R) = (4/3) pi R^3 for uniform rho=1
        expected = (4.0 / 3.0) * np.pi * r ** 3
        # skip first few points where trapezoidal rule has large relative error
        np.testing.assert_allclose(M[5:], expected[5:], rtol=0.02)


class TestKuramotoSolver:
    def test_converges(self):
        r = np.linspace(1, 20, 40)
        omega = 1.0 / r  # decreasing natural frequency
        solver = KuramotoSolver(r, omega, K_coupling=5.0, damping=0.3)
        state = solver.solve()
        assert state.n_iter < 500
        assert np.all(state.coherence >= 0)
        assert np.all(state.coherence <= 1)

    def test_strong_coupling_gives_high_coherence(self):
        r = np.linspace(1, 10, 30)
        omega = 0.1 * np.ones_like(r)
        solver = KuramotoSolver(r, omega, K_coupling=100.0, damping=0.3)
        state = solver.solve()
        # strong coupling => nearly full synchronisation
        assert np.mean(state.coherence) > 0.5

    def test_weak_coupling_gives_low_coherence(self):
        r = np.linspace(1, 10, 30)
        omega = 10.0 * np.ones_like(r)
        solver = KuramotoSolver(r, omega, K_coupling=0.01, damping=0.5)
        state = solver.solve()
        assert np.mean(state.coherence) < 0.5


class TestFiniteNFluctuations:
    def test_above_threshold(self):
        result = finite_N_fluctuations(N=100, K=5.0, omega_spread=1.0,
                                       n_samples=50)
        assert result["x"] > 1.0
        assert result["r_mean"] > 0.1

    def test_below_threshold(self):
        result = finite_N_fluctuations(N=100, K=0.1, omega_spread=1.0,
                                       n_samples=50)
        assert result["x"] < 1.0
        # below threshold, r should be small (but finite-N gives fluctuations)
        assert result["r_mean"] < 0.5
