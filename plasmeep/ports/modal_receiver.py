"""
Modal overlap receiver for numerical port modes.
"""

from __future__ import annotations

from typing import Any, Dict, List

import meep as mp
import numpy as np

from plasmeep.ports.mode_profile import complex_overlap
from plasmeep.ports.numerical_mode import NumericalPortMode


def add_modal_overlap_monitor(
    sim: mp.Simulation,
    center_xy: np.ndarray,
    tangent_xy: np.ndarray,
    mode: NumericalPortMode,
    *,
    frequency: float,
) -> Dict[str, Any]:
    """
  Register point DFT Hz monitors along the port tangent at the monitor plane.

  Returns a monitor bundle consumed by extract_modal_coefficient().
  """
    center_xy = np.asarray(center_xy, dtype=float).reshape(2)
    tangent = np.asarray(tangent_xy, dtype=float).reshape(2)
    tangent = tangent / max(np.linalg.norm(tangent), 1e-30)
    offsets = np.asarray(mode.offsets_a, dtype=float)
    dft_objs: List = []
    for s in offsets:
        xy = center_xy + s * tangent
        dft_objs.append(
            sim.add_dft_fields(
                [mp.Hz],
                frequency,
                0,
                1,
                center=mp.Vector3(xy[0], xy[1], 0),
                size=mp.Vector3(0, 0, 0),
            )
        )
    return {
        "dft_objs": dft_objs,
        "center": center_xy,
        "tangent": tangent,
        "offsets": offsets,
        "mode": mode,
        "frequency": frequency,
    }


def sample_hz_line(sim: mp.Simulation, mon_info: Dict[str, Any]) -> np.ndarray:
    hz = [
        complex(np.squeeze(sim.get_dft_array(dft, mp.Hz, 0)))
        for dft in mon_info["dft_objs"]
    ]
    return np.asarray(hz, dtype=complex)


def extract_modal_coefficient(sim: mp.Simulation, mon_info: Dict[str, Any]) -> complex:
    """
  Complex modal amplitude a_n = ⟨φ_n | H⟩.

  Both φ_n and H are L2-normalized on the same discrete offset grid before
  inner product (see mode_profile.complex_overlap).

  Transmitted power proxy: |a_n|^2 when incident launch is normalized to
  unit modal power via the reference normalization run.
  """
    field = sample_hz_line(sim, mon_info)
    mode: NumericalPortMode = mon_info["mode"]
    return mode.receiver_coefficient(field)


def extract_modal_power(sim: mp.Simulation, mon_info: Dict[str, Any]) -> float:
    """Return |⟨φ_n|H⟩|^2."""
    coeff = extract_modal_coefficient(sim, mon_info)
    return float(np.abs(coeff) ** 2)


def modal_coefficient_metrics(
    coeff_ab: complex,
    coeff_ba: complex,
) -> Dict[str, float]:
    """Reciprocity diagnostics on complex modal amplitudes."""
    eps = 1e-30
    amp_ratio = abs(coeff_ab) / max(abs(coeff_ba), eps)
    amp_err_dB = 20.0 * np.log10(amp_ratio) if amp_ratio > 0 else float("inf")
    phase_diff_deg = float(np.degrees(np.angle(coeff_ab) - np.angle(coeff_ba)))
    mean_amp = 0.5 * (abs(coeff_ab) + abs(coeff_ba))
    sym_err = abs(coeff_ab - coeff_ba) / max(mean_amp, eps)
    return {
        "amp_ratio": float(amp_ratio),
        "amp_err_dB": float(amp_err_dB),
        "phase_diff_deg": phase_diff_deg,
        "complex_sym_err": float(sym_err),
    }
