"""Unit tests for reusable 2D dielectric horn / feed geometry."""

import math
import unittest

import meep as mp
import numpy as np

from plasmeep.lib import Plasmeep


A = 0.01
EPS = 4.0
FEED_W = 0.40
FEED_L = 1.20
FLARE_L = 0.80
APERTURE_W = 1.00


def _new_model(nx=8, ny=8):
    return Plasmeep(A, 10, 0.5, nx, ny, verbose=False)


def _add_straight(model, center=(0.0, 0.0), theta=0.0, **kwargs):
    return model.Add_Horn(
        center, theta, FEED_W, FEED_L, EPS, **kwargs
    )


def _add_flared(model, center=(0.0, 0.0), theta=0.0, **kwargs):
    return model.Add_Horn(
        center, theta, FEED_W, FEED_L, EPS,
        flare_length=FLARE_L, aperture_width=APERTURE_W, **kwargs
    )


def _prism_xyz(prism):
    return np.array([[v.x, v.y, v.z] for v in prism.vertices], dtype=float)


def _eps(geom):
    return float(geom.material.epsilon_diag.x)


def _project(points, origin, axis, normal):
    d = points - origin
    return d @ axis, d @ normal


class TestAddHorn(unittest.TestCase):
    def test_straight_feed_dimensions(self):
        model = _new_model()
        meta = _add_straight(model, center=(0.0, 0.0), theta=0.0)

        self.assertEqual(len(model.geometry), 1)
        self.assertIsInstance(model.geometry[0], mp.Prism)
        self.assertAlmostEqual(_eps(model.geometry[0]), EPS)
        self.assertIsNone(meta["flare_vertices"])
        self.assertAlmostEqual(meta["aperture_width"], FEED_W)
        self.assertAlmostEqual(meta["flare_length"], 0.0)

        verts = _prism_xyz(model.geometry[0])
        np.testing.assert_allclose(verts, meta["feed_vertices"], atol=1e-12)
        along, across = _project(
            verts, meta["center"], meta["axis"], meta["normal"]
        )
        self.assertAlmostEqual(along.max() - along.min(), FEED_L)
        self.assertAlmostEqual(across.max() - across.min(), FEED_W)
        np.testing.assert_allclose(along.max(), 0.5 * FEED_L, atol=1e-12)
        np.testing.assert_allclose(along.min(), -0.5 * FEED_L, atol=1e-12)

    def test_flare_geometry_dimensions(self):
        model = _new_model()
        meta = _add_flared(model, center=(0.0, 0.0), theta=0.0)

        self.assertEqual(len(model.geometry), 2)
        feed = model.geometry[0]
        flare = model.geometry[1]
        self.assertAlmostEqual(_eps(feed), EPS)
        self.assertAlmostEqual(_eps(flare), EPS)

        feed_along, feed_across = _project(
            _prism_xyz(feed), meta["center"], meta["axis"], meta["normal"]
        )
        self.assertAlmostEqual(feed_along.max() - feed_along.min(), FEED_L)
        self.assertAlmostEqual(feed_across.max() - feed_across.min(), FEED_W)

        flare_verts = _prism_xyz(flare)
        np.testing.assert_allclose(
            flare_verts, meta["flare_vertices"], atol=1e-12
        )
        along, across = _project(
            flare_verts, meta["feed_front"], meta["axis"], meta["normal"]
        )
        throat = np.isclose(along, 0.0, atol=1e-12)
        mouth = np.isclose(along, FLARE_L, atol=1e-12)
        self.assertEqual(int(np.count_nonzero(throat)), 2)
        self.assertEqual(int(np.count_nonzero(mouth)), 2)
        self.assertAlmostEqual(across[throat].max() - across[throat].min(), FEED_W)
        self.assertAlmostEqual(across[mouth].max() - across[mouth].min(), APERTURE_W)
        self.assertAlmostEqual(along.max() - along.min(), FLARE_L)

        expected_aperture = meta["feed_front"] + FLARE_L * meta["axis"]
        np.testing.assert_allclose(meta["aperture_center"], expected_aperture)

    def test_rotations_0_90_and_arbitrary(self):
        angles = (0.0, 0.5 * math.pi, 0.37)
        for theta in angles:
            with self.subTest(theta=theta):
                model = _new_model()
                meta = _add_flared(model, center=(0.0, 0.0), theta=theta)
                u = np.array([math.cos(theta), math.sin(theta), 0.0])
                n = np.array([-math.sin(theta), math.cos(theta), 0.0])
                np.testing.assert_allclose(meta["axis"], u, atol=1e-12)
                np.testing.assert_allclose(meta["normal"], n, atol=1e-12)
                self.assertAlmostEqual(meta["theta"], theta)

                feed_along, feed_across = _project(
                    meta["feed_vertices"], meta["center"], u, n
                )
                self.assertAlmostEqual(feed_along.max() - feed_along.min(), FEED_L)
                self.assertAlmostEqual(
                    feed_across.max() - feed_across.min(), FEED_W
                )
                np.testing.assert_allclose(
                    meta["feed_front"] - meta["feed_back"], FEED_L * u
                )

                mouth = meta["aperture_center"]
                np.testing.assert_allclose(
                    mouth - meta["feed_front"], FLARE_L * u, atol=1e-12
                )
                meep_feed = _prism_xyz(model.geometry[0])
                np.testing.assert_allclose(meep_feed, meta["feed_vertices"])

    def test_translated_centers(self):
        origin_model = _new_model()
        origin_meta = _add_flared(origin_model, center=(0.0, 0.0), theta=0.7)
        shift = np.array([1.5, -0.8, 0.0])
        model = _new_model()
        meta = _add_flared(model, center=shift[:2], theta=0.7)

        np.testing.assert_allclose(meta["center"], shift)
        np.testing.assert_allclose(meta["port_center"], shift)
        np.testing.assert_allclose(
            meta["feed_vertices"], origin_meta["feed_vertices"] + shift
        )
        np.testing.assert_allclose(
            meta["flare_vertices"], origin_meta["flare_vertices"] + shift
        )
        np.testing.assert_allclose(
            meta["port_p1"], origin_meta["port_p1"] + shift
        )
        np.testing.assert_allclose(
            _prism_xyz(model.geometry[0]),
            _prism_xyz(origin_model.geometry[0]) + shift,
        )

    def test_aperture_expands_outward(self):
        theta = 2.1
        center = np.array([-0.4, 0.6, 0.0])
        model = _new_model()
        meta = _add_flared(model, center=center[:2], theta=theta)
        u = meta["axis"]

        launch = (meta["aperture_center"] - meta["center"]) @ u
        self.assertGreater(launch, 0.5 * FEED_L)
        self.assertAlmostEqual(launch, 0.5 * FEED_L + FLARE_L)
        self.assertGreater(
            (meta["feed_front"] - meta["feed_back"]) @ u, 0.0
        )

        along, across = _project(
            meta["flare_vertices"], meta["feed_front"], u, meta["normal"]
        )
        width_at = {}
        for s, t in zip(along, across):
            width_at.setdefault(round(float(s), 12), []).append(t)
        throat_w = max(width_at[0.0]) - min(width_at[0.0])
        mouth_w = max(width_at[round(FLARE_L, 12)]) - min(
            width_at[round(FLARE_L, 12)]
        )
        self.assertAlmostEqual(throat_w, FEED_W)
        self.assertAlmostEqual(mouth_w, APERTURE_W)
        self.assertGreater(mouth_w, throat_w)

    def test_port_plane_in_uniform_feed(self):
        theta = 1.1
        center = np.array([0.25, -0.15, 0.0])
        model = _new_model()
        meta = _add_flared(model, center=center[:2], theta=theta)

        np.testing.assert_allclose(meta["port_center"], meta["center"])
        span = meta["port_p2"] - meta["port_p1"]
        self.assertAlmostEqual(float(span @ meta["axis"]), 0.0)
        self.assertAlmostEqual(np.linalg.norm(span), FEED_W)
        mid = 0.5 * (meta["port_p1"] + meta["port_p2"])
        np.testing.assert_allclose(mid, meta["center"])

        s_port = float((meta["port_center"] - meta["feed_back"]) @ meta["axis"])
        s_front = float((meta["feed_front"] - meta["feed_back"]) @ meta["axis"])
        s_aperture = float(
            (meta["aperture_center"] - meta["feed_back"]) @ meta["axis"]
        )
        self.assertAlmostEqual(s_port, 0.5 * FEED_L)
        self.assertAlmostEqual(s_front, FEED_L)
        self.assertGreater(s_front, s_port)
        self.assertGreater(s_aperture, s_front)
        self.assertLess(s_port, FEED_L)
        self.assertGreater(s_port, 0.0)

        # Cross-section endpoints sit on the uniform-feed sidewalls.
        for p in (meta["port_p1"], meta["port_p2"]):
            s = float((p - meta["center"]) @ meta["axis"])
            t = float((p - meta["center"]) @ meta["normal"])
            self.assertAlmostEqual(s, 0.0)
            self.assertAlmostEqual(abs(t), 0.5 * FEED_W)

    def test_invalid_dimensions_fail_clearly(self):
        model = _new_model()
        center = (0.0, 0.0)
        with self.assertRaises(ValueError):
            model.Add_Horn(center, 0.0, 0.0, FEED_L, EPS)
        with self.assertRaises(ValueError):
            model.Add_Horn(center, 0.0, -FEED_W, FEED_L, EPS)
        with self.assertRaises(ValueError):
            model.Add_Horn(center, 0.0, FEED_W, 0.0, EPS)
        with self.assertRaises(ValueError):
            model.Add_Horn(center, 0.0, FEED_W, FEED_L, 0.0)
        with self.assertRaises(ValueError):
            model.Add_Horn(center, 0.0, FEED_W, FEED_L, -1000.0)
        with self.assertRaises(ValueError):
            model.Add_Horn(center, 0.0, FEED_W, FEED_L, EPS, flare_length=-0.1)
        with self.assertRaises(ValueError):
            model.Add_Horn(
                center, 0.0, FEED_W, FEED_L, EPS, flare_length=FLARE_L
            )
        with self.assertRaises(ValueError):
            model.Add_Horn(
                center, 0.0, FEED_W, FEED_L, EPS,
                flare_length=FLARE_L, aperture_width=FEED_W,
            )
        with self.assertRaises(ValueError):
            model.Add_Horn(
                center, 0.0, FEED_W, FEED_L, EPS,
                flare_length=0.0, aperture_width=APERTURE_W,
            )
        with self.assertRaises(ValueError):
            model.Add_Horn((0.0,), 0.0, FEED_W, FEED_L, EPS)
        with self.assertRaises(TypeError):
            model.Add_Horn(center, 0.0, FEED_W, FEED_L)
        self.assertEqual(len(model.geometry), 0)

    def test_geometry_constructs_in_meep(self):
        model = _new_model()
        _add_flared(model, center=(0.0, 0.0), theta=math.pi / 5)
        sim = model.Get_Sim()
        self.assertEqual(len(sim.geometry), 2)
        sim.init_sim()
        sim.reset_meep()

    def test_helper_is_not_hardcoded_to_five_horns(self):
        model = _new_model()
        _add_straight(model, center=(1.5, 0.0), theta=0.0)
        _add_straight(model, center=(0.0, 1.5), theta=0.5 * math.pi)
        self.assertEqual(len(model.geometry), 2)


if __name__ == "__main__":
    unittest.main()
