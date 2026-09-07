"""
Exact-DOF grid-index modal launch/receive helpers for Hz port diagnostics.

This is a targeted reciprocity diagnostic: source and receiver use the same
snapped Hz Yee DOFs and the same discrete complex mode vector on that line.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import meep as mp
import numpy as np

from plasmeep.ports.grid_index_port import GridIndexPortLine, build_grid_index_port_line
from plasmeep.ports.lorentz_probe import hz_yee_site
from plasmeep.ports.numerical_mode import NumericalPortMode


@dataclass(frozen=True)
class GridModalProfile:
    """Discrete port mode collapsed onto one exact Hz Yee cross-section."""

    line: GridIndexPortLine
    phi: np.ndarray  # unit-L2 normalized complex mode on line.grid_xy order

    def as_dict(self) -> Dict[str, Any]:
        return {
            "line": self.line.as_dict(),
            "phi_abs": np.abs(self.phi).tolist(),
            "phi_phase_deg": np.degrees(np.angle(self.phi)).tolist(),
        }

    def source_amplitudes(self, target_power: float = 1.0) -> np.ndarray:
        """Adjoint-matched source weights on the exact same DOFs."""
        phi = np.asarray(self.phi, dtype=complex).reshape(-1)
        wsum = float(np.sum(np.abs(phi) ** 2))
        if wsum <= 0:
            return np.ones_like(phi) / max(len(phi), 1)
        scale = np.sqrt(target_power / wsum)
        return np.conj(phi) * scale

    def receiver_coefficient(self, field: np.ndarray) -> complex:
        """Adjoint modal coefficient <phi|H> on the same exact DOFs."""
        field = np.asarray(field, dtype=complex).reshape(-1)
        phi = np.asarray(self.phi, dtype=complex).reshape(-1)
        if field.shape != phi.shape:
            raise ValueError(
                f"field/mode shape mismatch: field={field.shape} mode={phi.shape}"
            )
        return complex(np.vdot(phi, field))


def build_grid_modal_profile(
    mode: NumericalPortMode,
    *,
    center_xy: np.ndarray,
    outward_dir: np.ndarray,
    res: int,
    span_a: Optional[float] = None,
    n_target: Optional[int] = None,
) -> GridModalProfile:
    """
    Collapse a continuous cached mode profile onto a unique snapped Hz grid line.

    Each cached mode sample is snapped to the nearest Hz Yee site at the target
    cross-section; duplicate snaps are merged by summing complex amplitudes.
    The final vector is normalized to unit L2 norm on the unique line DOFs.
    """
    center_xy = np.asarray(center_xy, dtype=float).reshape(2)
    u = np.asarray(outward_dir, dtype=float).reshape(2)
    u = u / max(np.linalg.norm(u), 1e-30)
    tangent = np.array([-u[1], u[0]])
    offsets = np.asarray(mode.offsets_a, dtype=float)
    field = np.asarray(mode.normalized_field(), dtype=complex)
    span = float(mode.span_a if span_a is None else span_a)
    line = build_grid_index_port_line(
        mode.port_index,
        res=res,
        center_xy=center_xy,
        outward_dir=u,
        span_a=span,
        n_target=int(len(offsets) if n_target is None else n_target),
    )

    merged: Dict[tuple[int, int], complex] = {}
    for s, phi_k in zip(offsets, field):
        site = hz_yee_site(center_xy + float(s) * tangent, res=res)
        merged[site.cell_index] = merged.get(site.cell_index, 0.0 + 0.0j) + complex(phi_k)

    vals: List[complex] = []
    ix0, iy0 = line.anchor_cell
    for dix, diy in line.relative_cells:
        key = (ix0 + dix, iy0 + diy)
        vals.append(merged.get(key, 0.0 + 0.0j))
    phi = np.asarray(vals, dtype=complex)
    nrm = float(np.sqrt(np.sum(np.abs(phi) ** 2)))
    if nrm > 0:
        phi = phi / nrm
    return GridModalProfile(line=line, phi=phi)


def make_grid_modal_hz_sources(
    profile: GridModalProfile,
    *,
    frequency: float,
    fwidth: float,
    target_power: float = 1.0,
) -> List[mp.Source]:
    """Launch the exact discrete modal profile on its snapped Hz Yee DOFs."""
    amps = profile.source_amplitudes(target_power=target_power)
    sources: List[mp.Source] = []
    for xy, amp in zip(profile.line.grid_xy, amps):
        if abs(amp) < 1e-30:
            continue
        sources.append(
            mp.Source(
                src=mp.GaussianSource(frequency=frequency, fwidth=fwidth),
                component=mp.Hz,
                center=mp.Vector3(float(xy[0]), float(xy[1]), 0),
                amplitude=complex(amp),
            )
        )
    return sources


def add_grid_modal_overlap_monitor(
    sim: mp.Simulation,
    profile: GridModalProfile,
    *,
    frequency: float,
) -> Dict[str, Any]:
    """Register Hz DFT monitors on the exact line.grid_xy DOFs."""
    dft_objs: List[Any] = []
    for x, y in profile.line.grid_xy:
        dft_objs.append(
            sim.add_dft_fields(
                [mp.Hz],
                frequency,
                0,
                1,
                center=mp.Vector3(float(x), float(y), 0),
                size=mp.Vector3(),
            )
        )
    return {"dft_objs": dft_objs, "profile": profile, "frequency": frequency}


def sample_grid_modal_field(sim: mp.Simulation, mon_info: Dict[str, Any]) -> np.ndarray:
    hz = [
        complex(np.squeeze(sim.get_dft_array(dft, mp.Hz, 0)))
        for dft in mon_info["dft_objs"]
    ]
    return np.asarray(hz, dtype=complex)


def extract_grid_modal_coefficient(sim: mp.Simulation, mon_info: Dict[str, Any]) -> complex:
    field = sample_grid_modal_field(sim, mon_info)
    profile: GridModalProfile = mon_info["profile"]
    return profile.receiver_coefficient(field)
