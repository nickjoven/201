"""Tests for the main Calculator class."""

import numpy as np
import pytest

from sparc_x.calculator import Calculator
from sparc_x.profiles import GalaxyProfile, exponential_disk


@pytest.fixture
def synthetic_profile():
    r = np.linspace(0.5, 30, 60)
    v_disk = exponential_disk(r, M_disk=5e10, R_d=3.0)
    return GalaxyProfile(
        name="test_galaxy",
        r_kpc=r,
        v_obs=v_disk * 1.3,  # observed > baryonic (needs dark matter)
        v_gas=np.zeros_like(r),
        v_disk=v_disk,
        v_bul=np.zeros_like(r),
        e_v_obs=np.ones_like(r) * 5.0,
    )


class TestCalculator:
    def test_mond_rotation_curve(self, synthetic_profile):
        calc = Calculator(synthetic_profile)
        v = calc.get_rotation_curve(method="mond")
        assert v.shape == synthetic_profile.r_kpc.shape
        assert np.all(np.isfinite(v))
        assert np.all(v > 0)

    def test_accelerations(self, synthetic_profile):
        calc = Calculator(synthetic_profile)
        acc = calc.get_accelerations()
        assert "a_N" in acc
        assert "a_obs" in acc
        assert "a_pred" in acc
        assert np.all(acc["a_pred"] >= acc["a_N"])

    def test_dark_matter_density(self, synthetic_profile):
        calc = Calculator(synthetic_profile)
        rho = calc.get_dark_matter_density()
        assert rho.shape == synthetic_profile.r_kpc.shape
        assert np.all(rho >= 0)

    def test_prefactors(self, synthetic_profile):
        calc = Calculator(synthetic_profile)
        eta = calc.get_prefactors()
        assert eta.shape == synthetic_profile.r_kpc.shape
        assert np.all(np.isfinite(eta))

    def test_regime_mask(self, synthetic_profile):
        calc = Calculator(synthetic_profile)
        masks = calc.get_regime_mask()
        assert "newtonian" in masks
        assert "deep_mond" in masks
        assert "transition" in masks

    def test_summary(self, synthetic_profile):
        calc = Calculator(synthetic_profile)
        s = calc.summary()
        assert s["name"] == "test_galaxy"
        assert s["n_points"] == 60

    def test_sigma(self, synthetic_profile):
        calc = Calculator(synthetic_profile)
        sigma = calc.get_sigma()
        assert sigma.shape == synthetic_profile.r_kpc.shape
        assert np.all(sigma >= 0)

    def test_kuramoto_method(self, synthetic_profile):
        calc = Calculator(synthetic_profile, K_coupling=5.0)
        v = calc.get_rotation_curve(method="kuramoto")
        assert v.shape == synthetic_profile.r_kpc.shape
        assert np.all(np.isfinite(v))

    def test_lyapunov(self, synthetic_profile):
        calc = Calculator(synthetic_profile, K_coupling=5.0)
        ly = calc.get_lyapunov()
        assert ly.dVdt <= 0  # Lyapunov function must decrease
        assert len(ly.eigenvalues) == len(synthetic_profile.r_kpc)
