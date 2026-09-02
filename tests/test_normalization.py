"""Unit tests for stored complex S-parameter normalization."""

import unittest

import meep as mp
import numpy as np

from plasmeep.normalization import (
    IncidentNormalization,
    NormalizationTable,
)


FREQS = np.array([0.20, 0.25, 0.30])
C_PLUS = np.array([1.0 + 0.5j, 2.0 - 0.25j, -0.5 + 1.25j], dtype=complex)
PARITY = int(mp.ODD_Z)


def _norm(**kwargs):
    defaults = dict(
        coefficients=C_PLUS,
        frequencies=FREQS,
        source_id="port_a",
        eig_band=1,
        eig_parity=PARITY,
    )
    defaults.update(kwargs)
    return IncidentNormalization(**defaults)


class TestIncidentNormalization(unittest.TestCase):
    def test_is_complex_amplitude_ratio_not_power(self):
        norm = _norm()
        outgoing = np.array(
            [0.5 + 0.25j, 1.0 - 0.125j, -0.25 + 0.625j], dtype=complex
        )
        S = norm.normalize(outgoing, source_id="port_a")
        expected = outgoing / C_PLUS
        np.testing.assert_allclose(S.real, expected.real)
        np.testing.assert_allclose(S.imag, expected.imag)
        power_ratio = (np.abs(outgoing) ** 2) / (np.abs(C_PLUS) ** 2)
        np.testing.assert_allclose(S, 0.5 + 0j)
        self.assertFalse(np.allclose(np.abs(S), power_ratio))

    def test_frequency_and_source_indexing(self):
        norm = _norm()
        outgoing = 2.0 * C_PLUS
        S = norm.normalize(
            outgoing, source_id="port_a", frequencies=FREQS
        )
        np.testing.assert_allclose(S, 2.0 + 0j)
        S_one = _norm(
            coefficients=[C_PLUS[1]], frequencies=[FREQS[1]]
        ).normalize([3.0 * C_PLUS[1]], source_id="port_a")
        self.assertAlmostEqual(S_one[0].real, 3.0)
        self.assertAlmostEqual(S_one[0].imag, 0.0)

    def test_incompatible_source_frequency_and_mode_fail(self):
        norm = _norm()
        outgoing = C_PLUS
        with self.assertRaises(ValueError):
            norm.normalize(outgoing, source_id="port_b")
        with self.assertRaises(ValueError):
            norm.normalize(
                outgoing, source_id="port_a", frequencies=[0.20, 0.25, 0.31]
            )
        with self.assertRaises(ValueError):
            norm.normalize(
                outgoing, source_id="port_a", frequencies=[0.20, 0.25]
            )
        with self.assertRaises(ValueError):
            norm.normalize(outgoing, source_id="port_a", eig_band=2)
        with self.assertRaises(ValueError):
            norm.normalize(
                outgoing, source_id="port_a", eig_parity=int(mp.NO_PARITY)
            )

    def test_zero_denominator_fails_clearly(self):
        with self.assertRaises(ValueError):
            _norm(coefficients=[0.0], frequencies=[0.25])
        near_zero = np.array(C_PLUS, copy=True)
        near_zero[1] = 1e-20
        with self.assertRaises(ValueError):
            _norm(coefficients=near_zero)

    def test_reflection_and_transmission_share_denominator(self):
        norm = _norm()
        c_refl = np.array([0.1j, -0.2, 0.05 - 0.05j], dtype=complex)
        c_tran = np.array([0.9, 1.8 + 0.1j, -0.4 + 1.0j], dtype=complex)
        S11 = norm.s_parameter(c_refl, source_id="port_a")
        S21 = norm.s_parameter(c_tran, source_id="port_a")
        np.testing.assert_allclose(S11, c_refl / C_PLUS)
        np.testing.assert_allclose(S21, c_tran / C_PLUS)
        np.testing.assert_allclose(S11 * C_PLUS, c_refl)
        np.testing.assert_allclose(S21 * C_PLUS, c_tran)
        self.assertTrue(np.allclose(S11 * C_PLUS / c_refl, 1.0))
        self.assertTrue(np.allclose(S21 * C_PLUS / c_tran, 1.0))

    def test_from_port_and_table(self):
        port = {
            "frequency": 0.25,
            "frequencies": None,
            "eig_band": 1,
            "eig_parity": PARITY,
            "direction": "toward_aperture",
        }
        stored = IncidentNormalization.from_port(
            port, [2.0 + 1.0j], source_id="in"
        )
        S = stored.normalize([1.0 + 0.5j], source_id="in")
        np.testing.assert_allclose(S, 0.5 + 0j)

        table = NormalizationTable()
        table.add(stored)
        with self.assertRaises(ValueError):
            table.add(stored)
        with self.assertRaises(ValueError):
            table.normalize([1.0], source_id="other")
        S_table = table.normalize([1.0 + 0.5j], source_id="in")
        np.testing.assert_allclose(S_table, 0.5 + 0j)
        self.assertIn("in", table)
        self.assertEqual(len(table), 1)


if __name__ == "__main__":
    unittest.main()
