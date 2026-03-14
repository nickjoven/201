"""Tests for Stribeck friction curve."""

import numpy as np
import pytest

from sparc_x.stribeck import stribeck_curve, stribeck_to_mond, mond_from_stribeck


class TestStribeckCurve:
    def test_limits(self):
        # At v=0: mu_f = mu_s
        assert stribeck_curve(0.0) == pytest.approx(1.0)
        # At v >> v_s: mu_f -> mu_k
        assert stribeck_curve(100.0) == pytest.approx(0.4, abs=0.01)

    def test_monotone_decreasing(self):
        v = np.linspace(0, 10, 50)
        mu = stribeck_curve(v)
        assert np.all(np.diff(mu) <= 1e-10)


class TestMondFromStribeck:
    def test_recovers_rar(self):
        """delta=0.5 gives mu = 1 - exp(-sqrt(x)), the RAR form."""
        x = np.array([0.01, 0.1, 1.0, 10.0])
        mu = mond_from_stribeck(x, delta=0.5)
        expected = 1.0 - np.exp(-np.sqrt(x))
        np.testing.assert_allclose(mu, expected, rtol=1e-12)

    def test_limits(self):
        assert mond_from_stribeck(0.0, delta=0.5) == pytest.approx(0.0, abs=1e-8)
        assert mond_from_stribeck(100.0, delta=0.5) == pytest.approx(1.0, abs=1e-3)
