"""Unit tests for the first-lobe Bessel / diffusion plasma profile."""

import unittest

import numpy as np
import scipy.special as sp

from plasmeep.lib import J01, Bessel_profile, beta_from_h, h_from_beta


N0 = 1.0e18
R_P = 0.015
PRIOR_H = (0.05, 0.2, 0.5, 0.8)


class TestBesselProfile(unittest.TestCase):
    def test_axis_density_is_n0(self):
        for h in PRIOR_H:
            self.assertAlmostEqual(
                Bessel_profile(0.0, N0, R_P, h=h), N0, places=9
            )
        for beta in (0.5, 1.8, J01):
            self.assertAlmostEqual(
                Bessel_profile(0.0, N0, R_P, beta=beta), N0, places=9
            )

    def test_edge_density_is_h_times_n0(self):
        for h in (0.0,) + PRIOR_H:
            n_edge = Bessel_profile(R_P, N0, R_P, h=h)
            self.assertAlmostEqual(n_edge, h * N0, delta=1e-9 * N0)

    def test_beta_h_round_trip(self):
        betas = np.array([0.05, 0.4, 1.2, 1.8, 2.3, J01])
        np.testing.assert_allclose(
            beta_from_h(h_from_beta(betas)), betas, rtol=1e-10, atol=1e-12
        )
        hs = np.array([0.0, 0.05, 0.2, 0.5, 0.8, 0.999])
        np.testing.assert_allclose(
            h_from_beta(beta_from_h(hs)), hs, rtol=1e-10, atol=1e-12
        )

    def test_profile_nonnegative_and_decreasing_on_column(self):
        r = np.linspace(0.0, R_P, 201)
        for h in (0.0,) + PRIOR_H:
            n = Bessel_profile(r, N0, R_P, h=h)
            self.assertTrue(np.all(n >= -1e-15 * N0))
            self.assertTrue(np.all(np.diff(n) <= 1e-12 * N0))

    def test_h_is_j0_beta_not_beta_itself(self):
        beta = 1.8
        h = h_from_beta(beta)
        self.assertAlmostEqual(h, float(sp.j0(beta)), places=12)
        self.assertNotAlmostEqual(h, beta)
        self.assertGreater(h, 0.0)
        self.assertLess(h, 1.0)

    def test_schottky_and_near_flat_limits(self):
        self.assertAlmostEqual(beta_from_h(0.0), J01, places=12)
        self.assertEqual(h_from_beta(J01), 0.0)
        self.assertAlmostEqual(
            Bessel_profile(R_P, N0, R_P, beta=J01), 0.0, delta=1e-12 * N0
        )
        self.assertLess(beta_from_h(0.8), beta_from_h(0.05))

    def test_invalid_h(self):
        for h in (-0.01, 1.0, 1.2):
            with self.assertRaises(ValueError):
                Bessel_profile(0.0, N0, R_P, h=h)
            with self.assertRaises(ValueError):
                beta_from_h(h)

    def test_invalid_beta(self):
        for beta in (0.0, -1.0, J01 + 1e-6):
            with self.assertRaises(ValueError):
                Bessel_profile(0.0, N0, R_P, beta=beta)
            with self.assertRaises(ValueError):
                h_from_beta(beta)

    def test_invalid_radius_and_density(self):
        with self.assertRaises(ValueError):
            Bessel_profile(0.0, N0, 0.0, h=0.5)
        with self.assertRaises(ValueError):
            Bessel_profile(0.0, N0, -R_P, h=0.5)
        with self.assertRaises(ValueError):
            Bessel_profile(0.0, -N0, R_P, h=0.5)

    def test_invalid_radial_coordinates(self):
        with self.assertRaises(ValueError):
            Bessel_profile(-1e-9, N0, R_P, h=0.5)
        with self.assertRaises(ValueError):
            Bessel_profile(R_P + 1e-9, N0, R_P, h=0.5)
        with self.assertRaises(ValueError):
            Bessel_profile(np.array([0.0, R_P, 1.1 * R_P]), N0, R_P, h=0.5)

    def test_beta_and_h_exclusivity(self):
        with self.assertRaises(ValueError):
            Bessel_profile(0.0, N0, R_P)
        with self.assertRaises(ValueError):
            Bessel_profile(0.0, N0, R_P, beta=1.0, h=0.5)


if __name__ == "__main__":
    unittest.main()
