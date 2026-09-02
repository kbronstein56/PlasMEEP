"""
Direct Lorentz reciprocity probes (no port power normalization).

For 2D TE (Hz, Ex, Ey) in a passive linear geometry at B=0, compare complex
Hz transfer when a localized Hz current is placed at port A and measured at B
versus the swapped configuration.

Also provides a matched discrete overlap test: compact source/receiver
distributions on the same Hz Yee grid DOFs with identical weights in both roles.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Literal, Optional, Tuple

import meep as mp
import numpy as np

WeightKind = Literal["uniform", "gaussian"]


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
        h_ab = complex(self.a_to_b.hz_complex)
        h_ba = complex(self.b_to_a.hz_complex)
        return {
            "port_a": self.port_a,
            "port_b": self.port_b,
            "a_to_b_hz": h_ab,
            "b_to_a_hz": h_ba,
            "abs_a_to_b": float(abs(h_ab)),
            "abs_b_to_a": float(abs(h_ba)),
            "amp_ratio": self.amp_ratio,
            "amp_err_dB": self.amp_err_dB,
            "phase_diff_deg": self.phase_diff_deg,
            "complex_sym_err": self.complex_sym_err,
            "source_a_xy": list(self.a_to_b.source_xy),
            "monitor_b_xy": list(self.a_to_b.monitor_xy),
            "source_b_xy": list(self.b_to_a.source_xy),
            "monitor_a_xy": list(self.b_to_a.monitor_xy),
        }


def lorentz_test_specification() -> Dict[str, Any]:
    """Document the direct Lorentz reciprocity observable (Phase 1 audit)."""
    return {
        "observable": "complex Hz DFT at monitor point",
        "source_type": "mp.Source with mp.GaussianSource, component=mp.Hz",
        "source_amplitude": "1.0 (real, positive)",
        "source_phase": "0 rad (Meep default)",
        "receiver_component": "mp.Hz",
        "reciprocity_relation": (
            "For passive linear reciprocal 2D TE at B=0: H_AB should equal H_BA "
            "(same complex transfer when source/monitor roles swap). "
            "NOT H_AB = -H_BA (that would be odd symmetry)."
        ),
        "swap_convention": "A excites at port_a source_center, measured at port_b monitor_center; "
        "then roles swap.",
        "dft_frequency": "fs_a (single frequency, index 0)",
        "timing": "sim.run(until_after_sources=run_time) — fixed duration after source turn-off",
        "normalization": "none (raw complex Hz DFT coefficient)",
        "small_field_warning": (
            "If |H| < 1e-3, dB amplitude ratios are numerically fragile; "
            "check complex_sym_err and raw magnitudes."
        ),
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


@dataclass(frozen=True)
class HzYeeSite:
    """One Hz field DOF on the Yee grid (2D cell center)."""

    requested_xy: Tuple[float, float]
    grid_xy: Tuple[float, float]
    cell_index: Tuple[int, int]
    fractional_offset_cells: Tuple[float, float]

    def as_dict(self) -> Dict[str, Any]:
        return {
            "requested_xy": list(self.requested_xy),
            "grid_xy": list(self.grid_xy),
            "cell_index": list(self.cell_index),
            "fractional_offset_cells": list(self.fractional_offset_cells),
            "component": "Hz (Yee cell center)",
        }


@dataclass(frozen=True)
class HzGridPatch:
    """Compact Hz distribution on exact Yee grid points."""

    label: str
    requested_center: Tuple[float, float]
    coords: np.ndarray
    weights: np.ndarray
    cell_indices: np.ndarray
    half_width_cells: int
    res: int

    def as_dict(self) -> Dict[str, Any]:
        return {
            "label": self.label,
            "requested_center": list(self.requested_center),
            "half_width_cells": self.half_width_cells,
            "res": self.res,
            "n_dofs": int(self.coords.shape[0]),
            "coords": self.coords.tolist(),
            "weights": self.weights.tolist(),
            "cell_indices": self.cell_indices.tolist(),
        }


@dataclass(frozen=True)
class DiscreteTransfer:
    launcher_label: str
    receiver_label: str
    response: complex

    def as_dict(self) -> Dict[str, Any]:
        return {
            "launcher": self.launcher_label,
            "receiver": self.receiver_label,
            "response": complex(self.response),
            "abs_response": float(abs(self.response)),
        }


@dataclass(frozen=True)
class DiscreteReciprocityResult:
    patch_a: HzGridPatch
    patch_b: HzGridPatch
    a_to_b: DiscreteTransfer
    b_to_a: DiscreteTransfer
    amp_ratio: float
    amp_err_dB: float
    phase_diff_deg: float
    complex_sym_err: float

    def as_dict(self) -> Dict[str, Any]:
        g_ab = complex(self.a_to_b.response)
        g_ba = complex(self.b_to_a.response)
        return {
            "patch_a": self.patch_a.as_dict(),
            "patch_b": self.patch_b.as_dict(),
            "g_a_to_b": g_ab,
            "g_b_to_a": g_ba,
            "abs_g_a_to_b": float(abs(g_ab)),
            "abs_g_b_to_a": float(abs(g_ba)),
            "amp_ratio": self.amp_ratio,
            "amp_err_dB": self.amp_err_dB,
            "phase_diff_deg": self.phase_diff_deg,
            "complex_sym_err": self.complex_sym_err,
            "fields_large_enough": abs(g_ab) > 1e-6 and abs(g_ba) > 1e-6,
        }


def hz_yee_site(requested_xy: np.ndarray, *, res: int) -> HzYeeSite:
    """Map a physical coordinate to the nearest Hz Yee cell-center DOF."""
    xy = np.asarray(requested_xy, dtype=float).reshape(2)
    ix = int(np.round(xy[0] * res - 0.5))
    iy = int(np.round(xy[1] * res - 0.5))
    grid = np.array([(ix + 0.5) / res, (iy + 0.5) / res], dtype=float)
    frac = (xy - grid) * res
    return HzYeeSite(
        requested_xy=(float(xy[0]), float(xy[1])),
        grid_xy=(float(grid[0]), float(grid[1])),
        cell_index=(ix, iy),
        fractional_offset_cells=(float(frac[0]), float(frac[1])),
    )


def build_hz_grid_patch(
    center_xy: np.ndarray,
    *,
    res: int,
    half_width_cells: int = 1,
    weight_kind: WeightKind = "uniform",
    label: str = "",
) -> HzGridPatch:
    """
    Build a compact Hz distribution on exact Yee cell centers.

    Weights are identical when the patch acts as source or receiver overlap.
    """
    center_xy = np.asarray(center_xy, dtype=float).reshape(2)
    site = hz_yee_site(center_xy, res=res)
    ix0, iy0 = site.cell_index
    hw = int(half_width_cells)
    coords: List[np.ndarray] = []
    indices: List[Tuple[int, int]] = []
    for di in range(-hw, hw + 1):
        for dj in range(-hw, hw + 1):
            ix = ix0 + di
            iy = iy0 + dj
            indices.append((ix, iy))
            coords.append(np.array([(ix + 0.5) / res, (iy + 0.5) / res], dtype=float))
    coords_arr = np.stack(coords, axis=0)
    idx_arr = np.asarray(indices, dtype=int)
    if weight_kind == "uniform":
        weights = np.ones(coords_arr.shape[0], dtype=float)
    elif weight_kind == "gaussian":
        sigma = max(hw / 2.0, 0.5)
        r2 = di_grid = None
        di_grid = np.array([i[0] - ix0 for i in indices], dtype=float)
        dj_grid = np.array([i[1] - iy0 for i in indices], dtype=float)
        r2 = di_grid**2 + dj_grid**2
        weights = np.exp(-0.5 * r2 / sigma**2)
    else:
        raise ValueError(f"unknown weight_kind: {weight_kind}")
    weights = weights / np.sum(weights)
    return HzGridPatch(
        label=label,
        requested_center=(float(center_xy[0]), float(center_xy[1])),
        coords=coords_arr,
        weights=weights,
        cell_indices=idx_arr,
        half_width_cells=hw,
        res=int(res),
    )


def _patch_weight_grid(
    patch: HzGridPatch,
    xs: np.ndarray,
    ys: np.ndarray,
) -> np.ndarray:
    """Map patch DOF weights onto a Meep DFT block grid (x-major indexing)."""
    wmap = np.zeros((len(xs), len(ys)), dtype=float)
    for xy, weight in zip(patch.coords, patch.weights):
        ix = int(np.argmin(np.abs(xs - xy[0])))
        iy = int(np.argmin(np.abs(ys - xy[1])))
        wmap[ix, iy] += float(weight)
    return wmap


def discrete_reciprocity_specification() -> Dict[str, Any]:
    """Document the matched discrete overlap reciprocity observable."""
    return {
        "observable": "weighted complex Hz DFT overlap on Yee grid DOFs",
        "source_type": "mp.Source (GaussianSource) per patch DOF, amplitude=weight_k",
        "receiver": "sum_k weight_k * Hz_DFT(x_k) on receiver patch DOFs",
        "reciprocity_relation": "G_AB = G_BA for passive linear reciprocal media at B=0",
        "grid_convention": "Hz at Yee cell centers (i+0.5)/res, (j+0.5)/res",
        "interpolation": "none — sources and receivers use identical discrete coordinates",
        "timing": "sim.run(until_after_sources=run_time)",
        "normalization": "weights sum to 1 per patch",
    }


def make_patch_hz_sources(
    patch: HzGridPatch,
    *,
    frequency: float,
    fwidth: float,
) -> List[mp.Source]:
    sources: List[mp.Source] = []
    for xy, weight in zip(patch.coords, patch.weights):
        if abs(weight) < 1e-15:
            continue
        sources.append(
            mp.Source(
                src=mp.GaussianSource(frequency=frequency, fwidth=fwidth),
                component=mp.Hz,
                center=mp.Vector3(float(xy[0]), float(xy[1]), 0),
                amplitude=float(weight),
            )
        )
    return sources


def measure_weighted_hz_overlap(
    sim: mp.Simulation,
    receiver_patch: HzGridPatch,
    *,
    frequency: float,
    run_time: float,
) -> complex:
    """Weighted sum of Hz DFT coefficients on exact receiver grid DOFs."""
    site = hz_yee_site(receiver_patch.requested_center, res=receiver_patch.res)
    hw = receiver_patch.half_width_cells
    span = (2 * hw + 1) / receiver_patch.res
    center = mp.Vector3(site.grid_xy[0], site.grid_xy[1], 0)
    size = mp.Vector3(span, span, 0)
    dft_obj = sim.add_dft_fields(
        [mp.Hz],
        frequency,
        0,
        1,
        center=center,
        size=size,
    )
    sim.run(until_after_sources=run_time)
    arr = np.squeeze(sim.get_dft_array(dft_obj, mp.Hz, 0))
    xs, ys, _, _ = sim.get_array_metadata(center=center, size=size)
    wmap = _patch_weight_grid(receiver_patch, np.asarray(xs), np.asarray(ys))
    if arr.shape != wmap.shape:
        raise ValueError(f"DFT grid {arr.shape} != weight grid {wmap.shape}")
    return complex(np.sum(arr * wmap))


def run_discrete_hz_transfer(
    p_device,
    source_patch: HzGridPatch,
    receiver_patch: HzGridPatch,
    *,
    frequency: float,
    fwidth: float,
    run_time: float,
) -> complex:
    """Excite source_patch; return weighted Hz overlap with receiver_patch."""
    p_device.sources = make_patch_hz_sources(
        source_patch, frequency=frequency, fwidth=fwidth
    )
    sim = p_device.Get_Sim()
    return measure_weighted_hz_overlap(
        sim, receiver_patch, frequency=frequency, run_time=run_time
    )


def evaluate_discrete_reciprocity_pair(
    p_device,
    patch_a: HzGridPatch,
    patch_b: HzGridPatch,
    *,
    frequency: float,
    fwidth: float,
    run_time: float,
) -> DiscreteReciprocityResult:
    g_ab = run_discrete_hz_transfer(
        p_device,
        patch_a,
        patch_b,
        frequency=frequency,
        fwidth=fwidth,
        run_time=run_time,
    )
    g_ba = run_discrete_hz_transfer(
        p_device,
        patch_b,
        patch_a,
        frequency=frequency,
        fwidth=fwidth,
        run_time=run_time,
    )
    amp_ratio, amp_err_dB, phase_diff_deg, sym_err = lorentz_pair_metrics(g_ab, g_ba)
    return DiscreteReciprocityResult(
        patch_a=patch_a,
        patch_b=patch_b,
        a_to_b=DiscreteTransfer(patch_a.label, patch_b.label, g_ab),
        b_to_a=DiscreteTransfer(patch_b.label, patch_a.label, g_ba),
        amp_ratio=amp_ratio,
        amp_err_dB=amp_err_dB,
        phase_diff_deg=phase_diff_deg,
        complex_sym_err=sym_err,
    )


def probe_sites_grid_audit(sites: ProbeSites, *, res: int) -> Dict[str, Any]:
    """Report Yee-grid registration for point-test source and monitor coordinates."""
    source_site = hz_yee_site(sites.source_xy, res=res)
    monitor_site = hz_yee_site(sites.monitor_xy, res=res)
    return {
        "port": sites.port_index,
        "source": source_site.as_dict(),
        "monitor": monitor_site.as_dict(),
    }


def evaluate_reciprocity_pair(
    p_device,
    sites_a: ProbeSites,
    sites_b: ProbeSites,
    *,
    res: int,
    frequency: float,
    fwidth: float,
    run_time: float,
    observable: Literal["point", "discrete", "both"] = "discrete",
    half_width_cells: int = 1,
    weight_kind: WeightKind = "uniform",
) -> Dict[str, Any]:
    """
    Canonical reciprocity evaluation.

    Prefer ``observable='discrete'`` for full-device PMM arrays: matched Yee-grid
    source/receiver patches avoid off-grid interpolation artifacts that inflate the
    legacy point-probe metric in strongly scattering geometries.
    """
    out: Dict[str, Any] = {"port_a": sites_a.port_index, "port_b": sites_b.port_index}
    if observable in ("point", "both"):
        out["point"] = evaluate_lorentz_pair(
            p_device,
            sites_a,
            sites_b,
            res=res,
            frequency=frequency,
            fwidth=fwidth,
            run_time=run_time,
        ).as_dict()
    if observable in ("discrete", "both"):
        patch_a = build_hz_grid_patch(
            sites_a.source_xy,
            res=res,
            half_width_cells=half_width_cells,
            weight_kind=weight_kind,
            label=f"P{sites_a.port_index + 1}_source",
        )
        patch_b = build_hz_grid_patch(
            sites_b.source_xy,
            res=res,
            half_width_cells=half_width_cells,
            weight_kind=weight_kind,
            label=f"P{sites_b.port_index + 1}_source",
        )
        out["discrete"] = evaluate_discrete_reciprocity_pair(
            p_device,
            patch_a,
            patch_b,
            frequency=frequency,
            fwidth=fwidth,
            run_time=run_time,
        ).as_dict()
    out["grid_audit"] = {
        "port_a": probe_sites_grid_audit(sites_a, res=res),
        "port_b": probe_sites_grid_audit(sites_b, res=res),
    }
    return out


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
