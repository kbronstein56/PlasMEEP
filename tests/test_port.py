"""Unit tests for reusable MEEP horn-port helpers."""

import math
import unittest

import meep as mp
import meep.adjoint as mpa
import numpy as np

from plasmeep.lib import (
    PORT_AWAY_FROM_APERTURE,
    PORT_TOWARD_APERTURE,
    Plasmeep,
)


A = 0.01
EPS = 4.0
FEED_W = 0.40
FEED_L = 1.20
FLARE_L = 0.80
APERTURE_W = 1.00
FCEN = 0.20


def _new_model(nx=8, ny=8):
    return Plasmeep(A, 10, 0.5, nx, ny, verbose=False)


def _add_flared_horn(model, center=(0.0, 0.0), theta=0.0):
    return model.Add_Horn(
        center, theta, FEED_W, FEED_L, EPS,
        flare_length=FLARE_L, aperture_width=APERTURE_W,
    )


def _xyz(vec):
    return np.array([vec.x, vec.y, vec.z], dtype=float)


class TestAddPort(unittest.TestCase):
    def test_source_and_monitor_center_match_horn_port(self):
        model = _new_model()
        horn = _add_flared_horn(model, center=(0.25, -0.15), theta=0.0)
        port = model.Add_Port(horn, PORT_TOWARD_APERTURE, frequency=FCEN)
        source = model.Add_Port_Source(port)
        sim = model.Get_Sim()
        monitor = model.Add_Port_Monitor(port, sim)

        np.testing.assert_allclose(port["port_center"], horn["port_center"])
        np.testing.assert_allclose(_xyz(source.center), horn["port_center"])
        np.testing.assert_allclose(
            _xyz(monitor.volume.center), horn["port_center"]
        )
        np.testing.assert_allclose(_xyz(port["volume"].center), horn["port_center"])

    def test_transverse_extent_equals_feed_width(self):
        model = _new_model()
        horn = _add_flared_horn(model, center=(0.0, 0.0), theta=0.0)
        port = model.Add_Port(horn, PORT_TOWARD_APERTURE, frequency=FCEN)
        source = model.Add_Port_Source(port)
        sim = model.Get_Sim()
        monitor = model.Add_Port_Monitor(port, sim)

        np.testing.assert_allclose(port["volume_size"], [0.0, FEED_W, 0.0])
        np.testing.assert_allclose(_xyz(source.size), [0.0, FEED_W, 0.0])
        np.testing.assert_allclose(_xyz(monitor.volume.size), [0.0, FEED_W, 0.0])
        self.assertAlmostEqual(np.linalg.norm(port["volume_size"]), FEED_W)
        self.assertAlmostEqual(
            np.linalg.norm(horn["port_p2"] - horn["port_p1"]), FEED_W
        )

    def test_arbitrary_horn_rotations(self):
        angles = (0.0, 0.5 * math.pi, 0.37, 2.0 * math.pi / 5.0)
        for theta in angles:
            with self.subTest(theta=theta):
                model = _new_model()
                horn = _add_flared_horn(model, center=(0.3, -0.2), theta=theta)
                port = model.Add_Port(
                    horn, PORT_TOWARD_APERTURE, frequency=FCEN
                )
                source = model.Add_Port_Source(port)
                np.testing.assert_allclose(
                    _xyz(source.center), horn["port_center"]
                )
                np.testing.assert_allclose(port["axis"], horn["axis"])
                np.testing.assert_allclose(_xyz(source.eig_kpoint), horn["axis"])
                self.assertAlmostEqual(
                    np.linalg.norm(port["volume_size"]), FEED_W
                )
                self.assertEqual(
                    int(np.count_nonzero(port["volume_size"])), 1
                )
                n = horn["normal"]
                if abs(n[0]) >= abs(n[1]):
                    self.assertAlmostEqual(port["volume_size"][0], FEED_W)
                    self.assertAlmostEqual(port["volume_size"][1], 0.0)
                else:
                    self.assertAlmostEqual(port["volume_size"][1], FEED_W)
                    self.assertAlmostEqual(port["volume_size"][0], 0.0)

    def test_toward_and_away_map_to_horn_axis(self):
        model = _new_model()
        horn = _add_flared_horn(model, center=(-0.4, 0.6), theta=2.1)
        toward = model.Add_Port(horn, PORT_TOWARD_APERTURE, frequency=FCEN)
        away = model.Add_Port(horn, PORT_AWAY_FROM_APERTURE, frequency=FCEN)
        src_toward = model.Add_Port_Source(toward)
        src_away = model.Add_Port_Source(away)

        np.testing.assert_allclose(toward["kpoint"], horn["axis"])
        np.testing.assert_allclose(away["kpoint"], -horn["axis"])
        np.testing.assert_allclose(
            toward["kpoint_toward_aperture"], horn["axis"]
        )
        np.testing.assert_allclose(
            toward["kpoint_away_from_aperture"], -horn["axis"]
        )
        np.testing.assert_allclose(_xyz(src_toward.eig_kpoint), horn["axis"])
        np.testing.assert_allclose(_xyz(src_away.eig_kpoint), -horn["axis"])
        self.assertEqual(src_toward.direction, mp.NO_DIRECTION)
        self.assertEqual(src_away.direction, mp.NO_DIRECTION)

        sim = model.Get_Sim()
        mon_toward = model.Add_Port_Monitor(toward, sim, PORT_TOWARD_APERTURE)
        mon_away = model.Add_Port_Monitor(toward, sim, PORT_AWAY_FROM_APERTURE)
        k_t = _xyz(mon_toward.kpoint_func(FCEN, 1))
        k_a = _xyz(mon_away.kpoint_func(FCEN, 1))
        np.testing.assert_allclose(k_t, horn["axis"])
        np.testing.assert_allclose(k_a, -horn["axis"])
        self.assertEqual(mon_toward.kpoint_func_overlap_idx, 0)
        self.assertEqual(mon_away.kpoint_func_overlap_idx, 0)

    def test_eigenmode_band_and_parity_preserved(self):
        model = _new_model()
        horn = _add_flared_horn(model)
        port = model.Add_Port(
            horn, PORT_TOWARD_APERTURE,
            frequency=FCEN, eig_band=2, eig_parity=mp.ODD_Z,
        )
        source = model.Add_Port_Source(port)
        sim = model.Get_Sim()
        monitor = model.Add_Port_Monitor(port, sim)

        self.assertEqual(port["eig_band"], 2)
        self.assertEqual(port["eig_parity"], mp.ODD_Z)
        self.assertEqual(source.eig_band, 2)
        self.assertEqual(source.eig_parity, mp.ODD_Z)
        self.assertEqual(monitor.mode, 2)
        self.assertEqual(monitor.eigenmode_kwargs["eig_parity"], mp.ODD_Z)

        default_port = model.Add_Port(horn, PORT_TOWARD_APERTURE, frequency=FCEN)
        self.assertEqual(default_port["eig_parity"], mp.ODD_Z)
        self.assertEqual(default_port["eig_band"], 1)

    def test_invalid_direction_band_and_frequency(self):
        model = _new_model()
        horn = _add_flared_horn(model)
        with self.assertRaises(ValueError):
            model.Add_Port(horn, "forward")
        with self.assertRaises(ValueError):
            model.Add_Port(horn, "+x")
        with self.assertRaises(ValueError):
            model.Add_Port(horn, PORT_TOWARD_APERTURE, eig_band=0)
        with self.assertRaises(ValueError):
            model.Add_Port(horn, PORT_TOWARD_APERTURE, eig_band=-1)
        with self.assertRaises(ValueError):
            model.Add_Port(horn, PORT_TOWARD_APERTURE, frequency=0.0)
        with self.assertRaises(ValueError):
            model.Add_Port(horn, PORT_TOWARD_APERTURE, frequency=-0.2)
        with self.assertRaises(ValueError):
            model.Add_Port(
                horn, PORT_TOWARD_APERTURE, frequencies=[0.2, -0.1]
            )
        with self.assertRaises(ValueError):
            model.Add_Port(horn, PORT_TOWARD_APERTURE, fwidth=0.0)
        with self.assertRaises(ValueError):
            model.Add_Port({}, PORT_TOWARD_APERTURE)
        port = model.Add_Port(horn, PORT_TOWARD_APERTURE)
        with self.assertRaises(ValueError):
            model.Add_Port_Source(port)
        with self.assertRaises(ValueError):
            model.Add_Port_Monitor(port, None)
        with self.assertRaises(ValueError):
            model.Add_Port_Monitor(port, object(), sense="incoming")

    def test_multiple_rotated_horns(self):
        model = _new_model()
        thetas = (0.0, 2.0 * math.pi / 5.0, 4.0 * math.pi / 5.0)
        ports = []
        for i, theta in enumerate(thetas):
            center = (1.5 * math.cos(theta), 1.5 * math.sin(theta))
            horn = _add_flared_horn(model, center=center, theta=theta)
            port = model.Add_Port(
                horn, PORT_TOWARD_APERTURE,
                frequency=FCEN, eig_band=1, eig_parity=mp.ODD_Z,
            )
            source = model.Add_Port_Source(port)
            ports.append((horn, port, source))
            np.testing.assert_allclose(_xyz(source.center), horn["port_center"])
            np.testing.assert_allclose(_xyz(source.eig_kpoint), horn["axis"])
            self.assertNotEqual(
                id(source), id(ports[0][2]) if i else id(object())
            )

        self.assertEqual(len(model.sources), 3)
        self.assertEqual(len(model.ports), 3)
        sim = model.Get_Sim()
        self.assertEqual(len(sim.sources), 3)
        for horn, port, source in ports:
            monitor = model.Add_Port_Monitor(
                port, sim, sense=PORT_AWAY_FROM_APERTURE
            )
            np.testing.assert_allclose(
                _xyz(monitor.volume.center), _xyz(source.center)
            )
            np.testing.assert_allclose(
                _xyz(monitor.volume.size), _xyz(source.size)
            )
            np.testing.assert_allclose(
                _xyz(monitor.kpoint_func(FCEN, 1)), -horn["axis"]
            )

    def test_meep_objects_construct_without_timestepping(self):
        model = _new_model()
        horn = _add_flared_horn(model, center=(0.0, 0.0), theta=math.pi / 5)
        port = model.Add_Port(
            horn, PORT_TOWARD_APERTURE,
            frequency=FCEN, fwidth=0.1 * FCEN, eig_band=1,
        )
        source = model.Add_Port_Source(port)
        self.assertIsInstance(source, mp.EigenModeSource)
        self.assertEqual(source.direction, mp.NO_DIRECTION)
        sim = model.Get_Sim()
        self.assertEqual(len(sim.sources), 1)
        self.assertIs(sim.sources[0], source)
        monitor = model.Add_Port_Monitor(port, sim)
        self.assertIsInstance(monitor, mpa.EigenmodeCoefficient)
        self.assertEqual(len(sim.geometry), 2)
        # Geometry-only init: no EigenModeSource MPB solve, no FDTD steps.
        geom_model = _new_model()
        _add_flared_horn(geom_model, center=(0.0, 0.0), theta=0.0)
        geom_sim = geom_model.Get_Sim()
        geom_sim.init_sim()
        geom_sim.reset_meep()

    def test_source_is_not_placed_in_the_flare(self):
        model = _new_model()
        horn = _add_flared_horn(model, theta=0.8)
        port = model.Add_Port(horn, PORT_TOWARD_APERTURE, frequency=FCEN)
        source = model.Add_Port_Source(port)
        s_src = float(
            (_xyz(source.center) - horn["feed_back"]) @ horn["axis"]
        )
        s_front = float((horn["feed_front"] - horn["feed_back"]) @ horn["axis"])
        s_ap = float(
            (horn["aperture_center"] - horn["feed_back"]) @ horn["axis"]
        )
        self.assertGreater(s_src, 0.0)
        self.assertLess(s_src, s_front)
        self.assertLess(s_front, s_ap)
        self.assertAlmostEqual(s_src, 0.5 * FEED_L)


if __name__ == "__main__":
    unittest.main()
