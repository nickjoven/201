"""Tests for MOND interpolation functions."""

import numpy as np
import pytest

from sparc_x.mond import (
    mu_simple, mu_standard, mu_rar,
    nu_simple, nu_rar,
    a_obs, sigma_static,
)


class TestInterpolatingFunctions:
    def test_mu_simple_limits(self):
        assert mu_simple(0.0) == pytest.approx(0.0, abs=1e-10)
        assert mu_simple(1e6) == pytest.approx(1.0, abs=1e-4)

    def test_mu_standard_limits(self):
        assert mu_standard(0.0) == pytest.approx(0.0, abs=1e-10)
        assert mu_standard(1e6) == pytest.approx(1.0, abs=1e-4)

    def test_mu_rar_limits(self):
        assert mu_rar(0.0) == pytest.approx(0.0, abs=1e-10)
        assert mu_rar(100.0) == pytest.approx(1.0, abs=1e-3)

    def test_mu_vectorized(self):
        x = np.array([0.01, 0.1, 1.0, 10.0])
        result = mu_simple(x)
        assert result.shape == (4,)
        assert np.all(result >= 0)
        assert np.all(result <= 1)

    def test_nu_simple_is_boost(self):
        """nu(y) * y should give x such that mu(x) * x = y."""
        x = np.array([0.1, 1.0, 10.0])
        # For the simple form: a_obs = nu(a_N/a0) * a_N
        # and nu is defined so that a_obs > a_N in MOND regime
        nu = nu_simple(x)
        assert np.all(nu >= 1.0)  # boost factor always >= 1


class TestAobs:
    def test_deep_mond(self):
        """In deep-MOND: a_obs ~ sqrt(a_N * a0)."""
        a_N = 1e-13  # well below a0
        result = a_obs(a_N, interpolation="simple")
        expected = np.sqrt(a_N * 1.2e-10)
        assert result == pytest.approx(expected, rel=0.5)

    def test_newtonian(self):
        """In Newtonian regime: a_obs ~ a_N."""
        a_N = 1e-7  # well above a0
        result = a_obs(a_N, interpolation="simple")
        assert result == pytest.approx(a_N, rel=0.01)

    def test_bad_interpolation(self):
        with pytest.raises(ValueError):
            a_obs(1.0, interpolation="nonexistent")


class TestSigmaStatic:
    def test_scaling(self):
        v = np.array([100.0, 200.0])
        s = sigma_static(v)
        # Sigma ~ v^4, so ratio should be 2^4 = 16
        assert s[1] / s[0] == pytest.approx(16.0, rel=1e-6)
