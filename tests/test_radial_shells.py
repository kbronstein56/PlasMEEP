"""Unit tests for concentric Drude-shell plasma columns."""

import unittest

import numpy as np

from plasmeep.lib import WP, Bessel_profile, Plasmeep


A = 0.01
R_P = 1.5
N0 = 1.0e18
CENTER = [0.0, 0.0, 0.0]


def _uniform(n):
    return lambda r: n


def _bessel(n0=N0, R_p=R_P, h=0.5):
    return lambda r: Bessel_profile(r, n0, R_p, h=h)


def _drude_freq(cylinder):
    susc = cylinder.material.E_susceptibilities
    if not susc:
        return 0.0
    return float(susc[0].frequency)


class TestAddRodRadialShells(unittest.TestCase):
    def setUp(self):
        self.model = Plasmeep(A, 10, 0.0, 4, 4, verbose=False)

    def test_requested_shell_count(self):
        shells = self.model.Add_Rod_radial_shells(
            R_P, CENTER, _uniform(N0), n_shells=12
        )
        self.assertEqual(len(shells), 12)
        self.assertEqual(len(self.model.geometry), 12)

    def test_default_shell_count_is_in_guide_range(self):
        shells = self.model.Add_Rod_radial_shells(R_P, CENTER, _uniform(N0))
        self.assertEqual(len(shells), 12)
        self.assertGreaterEqual(len(shells), 10)
        self.assertLessEqual(len(shells), 15)

    def test_radial_boundaries_cover_column(self):
        n_shells = 10
        shells = self.model.Add_Rod_radial_shells(
            R_P, CENTER, _uniform(N0), n_shells=n_shells
        )
        self.assertAlmostEqual(shells[0]["r_inner"], 0.0)
        self.assertAlmostEqual(shells[-1]["r_outer"], R_P)
        for i, shell in enumerate(shells):
            self.assertAlmostEqual(shell["r_inner"], i * R_P / n_shells)
            self.assertAlmostEqual(shell["r_outer"], (i + 1) * R_P / n_shells)
            self.assertLess(shell["r_inner"], shell["r_outer"])
        for left, right in zip(shells, shells[1:]):
            self.assertAlmostEqual(left["r_outer"], right["r_inner"])

    def test_profile_sampled_at_midpoints(self):
        recorded = []

        def profile(r):
            recorded.append(float(r))
            return N0

        n_shells = 8
        shells = self.model.Add_Rod_radial_shells(
            R_P, CENTER, profile, n_shells=n_shells
        )
        expected_mids = [(i + 0.5) * R_P / n_shells for i in range(n_shells)]
        self.assertEqual(len(recorded), n_shells)
        np.testing.assert_allclose(recorded, expected_mids, rtol=0.0, atol=1e-12)
        np.testing.assert_allclose(
            [s["r_mid"] for s in shells], expected_mids, rtol=0.0, atol=1e-12
        )

    def test_density_to_meep_wp_includes_2pi(self):
        wp_correct = self.model.density_to_meep_wp(N0)
        f_hz = WP(N0) / (2.0 * np.pi)
        self.assertAlmostEqual(
            wp_correct, self.model.Nondimensionalize_Freq(f_hz), places=12
        )
        wp_missing_2pi = self.model.Nondimensionalize_Freq(WP(N0))
        self.assertGreater(wp_missing_2pi / wp_correct, 6.0)

        shells = self.model.Add_Rod_radial_shells(
            R_P, CENTER, _uniform(N0), n_shells=4
        )
        for shell in shells:
            self.assertAlmostEqual(shell["wp"], wp_correct, places=12)
        # Get_Med is unchanged: Drude frequency is the pre-converted a-unit value.
        self.assertAlmostEqual(_drude_freq(self.model.geometry[0]), wp_correct)

    def test_bessel_profile_center_to_edge_trend(self):
        shells = self.model.Add_Rod_radial_shells(
            R_P, CENTER, _bessel(h=0.2), n_shells=12
        )
        n_vals = np.array([s["n_e"] for s in shells])
        wp_vals = np.array([s["wp"] for s in shells])
        self.assertTrue(np.all(np.diff(n_vals) < 0.0))
        self.assertTrue(np.all(np.diff(wp_vals) < 0.0))
        self.assertGreater(n_vals[0], n_vals[-1])
        self.assertAlmostEqual(
            shells[0]["n_e"], Bessel_profile(shells[0]["r_mid"], N0, R_P, h=0.2)
        )
        self.assertAlmostEqual(
            shells[-1]["n_e"],
            Bessel_profile(shells[-1]["r_mid"], N0, R_P, h=0.2),
        )

    def test_shell_ordering_later_cylinder_is_the_core(self):
        n_shells = 6
        shells = self.model.Add_Rod_radial_shells(
            R_P, CENTER, _bessel(h=0.3), n_shells=n_shells
        )
        geom = self.model.geometry
        self.assertAlmostEqual(geom[0].radius, R_P)
        self.assertAlmostEqual(geom[-1].radius, R_P / n_shells)
        radii = [cyl.radius for cyl in geom]
        self.assertTrue(np.all(np.diff(radii) < 0.0))
        # Last-added (smallest) cylinder carries the innermost midpoint density.
        self.assertAlmostEqual(_drude_freq(geom[-1]), shells[0]["wp"], places=12)
        self.assertAlmostEqual(_drude_freq(geom[0]), shells[-1]["wp"], places=12)
        self.assertGreater(_drude_freq(geom[-1]), _drude_freq(geom[0]))

    def test_zero_density_shell_has_no_drude_term(self):
        self.model.Add_Rod_radial_shells(
            R_P, CENTER, _uniform(0.0), n_shells=3
        )
        for cyl in self.model.geometry:
            self.assertEqual(len(cyl.material.E_susceptibilities), 0)
            self.assertEqual(_drude_freq(cyl), 0.0)

    def test_gamma_uses_add_rod_a_units_convention(self):
        gamma_a = self.model.Nondimensionalize_Freq(1.0e9)
        self.model.Add_Rod_radial_shells(
            R_P, CENTER, _uniform(N0), n_shells=2, gamma=gamma_a
        )
        susc = self.model.geometry[0].material.E_susceptibilities[0]
        self.assertAlmostEqual(float(susc.gamma), gamma_a)

    def test_invalid_radius_shell_count_and_profile(self):
        with self.assertRaises(ValueError):
            self.model.Add_Rod_radial_shells(0.0, CENTER, _uniform(N0))
        with self.assertRaises(ValueError):
            self.model.Add_Rod_radial_shells(-R_P, CENTER, _uniform(N0))
        with self.assertRaises(ValueError):
            self.model.Add_Rod_radial_shells(R_P, CENTER, _uniform(N0), n_shells=0)
        with self.assertRaises(ValueError):
            self.model.Add_Rod_radial_shells(R_P, CENTER, "not-callable")
        with self.assertRaises(ValueError):
            self.model.Add_Rod_radial_shells(
                R_P, CENTER, lambda r: -1.0, n_shells=2
            )
        with self.assertRaises(ValueError):
            self.model.Add_Rod_radial_shells(
                R_P, CENTER, lambda r: np.array([N0, N0]), n_shells=2
            )

    def test_tiny_geometry_constructs_without_timestepping(self):
        self.model.Add_Rod_radial_shells(
            R_P, CENTER, _bessel(h=0.5), n_shells=4
        )
        sim = self.model.Get_Sim()
        self.assertEqual(len(sim.geometry), 4)
        sim.init_sim()
        sim.reset_meep()


if __name__ == "__main__":
    unittest.main()
