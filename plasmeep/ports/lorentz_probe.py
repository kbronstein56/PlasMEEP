"""
Direct Lorentz reciprocity probes (no port power normalization).

For 2D TE (Hz, Ex, Ey) in a passive linear geometry at B=0, compare complex
Hz transfer when a localized Hz current is placed at port A and measured at B
versus the swapped configuration.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

import meep as mp
import numpy as np


@dataclass(frozen=True)
class ProbeSites:
    """Launch and receive coordinates for one port (a-units)."""

    port_index: int
    source_xy: np.ndarray
    monitor_xy: np.ndarray
    outward_dir: np.ndarray

    @property
    def tangent(self) -> np.ndarray:
        u = self.outward_dir / np.linalg.norm(self.outward_dir)
        return np.array([-u[1], u[0]])


@dataclass(frozen=True)
class LorentzTransfer:
    """Complex Hz DFT at receiver when Hz source is at launcher."""

    launcher: int
    receiver: int
    hz_complex: complex
    source_xy: Tuple[float, float]
    monitor_xy: Tuple[float, float]


@dataclass(frozen=True)
class LorentzPairResult:
    """Swapped transfers for one port pair."""

    port_a: int
    port_b: int
    a_to_b: LorentzTransfer
    b_to_a: LorentzTransfer
    amp_ratio: float
    amp_err_dB: float
    phase_diff_deg: float
    complex_sym_err: float

    def as_dict(self) -> Dict[str, Any]:
        return {
            "port_a": self.port_a,
            "port_b": self.port_b,
            "a_to_b_hz": complex(self.a_to_b.hz_complex),
            "b_to_a_hz": complex(self.b_to_a.hz_complex),
            "amp_ratio": self.amp_ratio,
            "amp_err_dB": self.amp_err_dB,
            "phase_diff_deg": self.phase_diff_deg,
            "complex_sym_err": self.complex_sym_err,
        }


def make_point_hz_source(
    xy: np.ndarray,
    *,
    frequency: float,
    fwidth: float,
    amplitude: float = 1.0,
) -> List[mp.Source]:
    xy = np.asarray(xy, dtype=float).reshape(2)
    return [
        mp.Source(
            src=mp.GaussianSource(frequency=frequency, fwidth=fwidth),
            component=mp.Hz,
            center=mp.Vector3(xy[0], xy[1], 0),
            amplitude=float(amplitude),
        )
    ]


def hz_dft_at_point(
    sim: mp.Simulation,
    xy: np.ndarray,
    *,
    frequency: float,
) -> complex:
    """Return complex Hz DFT coefficient at a single grid point."""
    xy = np.asarray(xy, dtype=float).reshape(2)
    cell = sim.fields.get_dft_cell_volume(
        mp.Vector3(xy[0], xy[1], 0), mp.Vector3(1, 1, 1)
    )
    dft_obj = sim.add_dft_fields([mp.Hz], frequency, 0, 1, where=cell)
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Hz, xy, 1e-4))
    arr = sim.get_dft_array(dft_obj, mp.Hz, 0)
    return complex(np.squeeze(arr))


def run_hz_transfer(
    p_device,
    *,
    res: int,
    source_xy: np.ndarray,
    monitor_xy: np.ndarray,
    frequency: float,
    fwidth: float,
    run_time: float,
) -> complex:
    """
    One-shot Hz point-source excitation; DFT Hz at monitor point.

    Uses until_after_sources=run_time (caller supplies decay-safe duration).
    """
    p_device.sources = make_point_hz_source(
        source_xy, frequency=frequency, fwidth=fwidth, amplitude=1.0
    )
    sim = p_device.Get_Sim()
    xy = np.asarray(monitor_xy, dtype=float).reshape(2)
    dft_obj = sim.add_dft_fields(
        [mp.Hz],
        frequency,
        0,
        1,
        center=mp.Vector3(xy[0], xy[1], 0),
        size=mp.Vector3(0, 0, 0),
    )
    sim.run(until_after_sources=run_time)
    arr = sim.get_dft_array(dft_obj, mp.Hz, 0)
    return complex(np.squeeze(arr))


def lorentz_pair_metrics(h_ab: complex, h_ba: complex) -> Tuple[float, float, float, float]:
    """Amplitude ratio, amplitude error (dB), phase difference (deg), |ΔH|/mean|H|."""
    eps = 1e-30
    amp_ratio = abs(h_ab) / max(abs(h_ba), eps)
    amp_err_dB = 20.0 * np.log10(amp_ratio) if amp_ratio > 0 else float("inf")
    phase_diff_deg = float(np.degrees(np.angle(h_ab) - np.angle(h_ba)))
    mean_amp = 0.5 * (abs(h_ab) + abs(h_ba))
    sym_err = abs(h_ab - h_ba) / max(mean_amp, eps)
    return float(amp_ratio), float(amp_err_dB), phase_diff_deg, float(sym_err)


def evaluate_lorentz_pair(
    p_device,
    sites_a: ProbeSites,
    sites_b: ProbeSites,
    *,
    res: int,
    frequency: float,
    fwidth: float,
    run_time: float,
) -> LorentzPairResult:
    """Swap localized Hz sources; compare complex Hz at the paired monitor."""
    h_ab = run_hz_transfer(
        p_device,
        res=res,
        source_xy=sites_a.source_xy,
        monitor_xy=sites_b.monitor_xy,
        frequency=frequency,
        fwidth=fwidth,
        run_time=run_time,
    )
    h_ba = run_hz_transfer(
        p_device,
        res=res,
        source_xy=sites_b.source_xy,
        monitor_xy=sites_a.monitor_xy,
        frequency=frequency,
        fwidth=fwidth,
        run_time=run_time,
    )
    amp_ratio, amp_err_dB, phase_diff_deg, sym_err = lorentz_pair_metrics(h_ab, h_ba)
    return LorentzPairResult(
        port_a=sites_a.port_index,
        port_b=sites_b.port_index,
        a_to_b=LorentzTransfer(
            sites_a.port_index,
            sites_b.port_index,
            h_ab,
            tuple(sites_a.source_xy),
            tuple(sites_b.monitor_xy),
        ),
        b_to_a=LorentzTransfer(
            sites_b.port_index,
            sites_a.port_index,
            h_ba,
            tuple(sites_b.source_xy),
            tuple(sites_a.monitor_xy),
        ),
        amp_ratio=amp_ratio,
        amp_err_dB=amp_err_dB,
        phase_diff_deg=phase_diff_deg,
        complex_sym_err=sym_err,
    )
