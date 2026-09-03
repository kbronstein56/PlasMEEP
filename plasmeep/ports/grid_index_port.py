"""
Grid-index-first port cross-sections on the Hz Yee lattice.

Physical horn geometry is unchanged; source/monitor sample DOFs are defined by
integer cell indices relative to a snapped anchor, then mapped to a-units.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np

from plasmeep.ports.lorentz_probe import hz_yee_site


@dataclass(frozen=True)
class GridIndexPortLine:
    """Hz Yee DOFs along a port cross-section."""

    port_index: int
    res: int
    anchor_cell: Tuple[int, int]
    relative_cells: List[Tuple[int, int]]
    grid_xy: List[Tuple[float, float]]
  # tangent/outward in a-units
    tangent: Tuple[float, float]
    outward_dir: Tuple[float, float]
    requested_center: Tuple[float, float]
    span_a: float

    @property
    def n_samples(self) -> int:
        return len(self.grid_xy)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "port_index": self.port_index,
            "res": self.res,
            "anchor_cell": list(self.anchor_cell),
            "relative_cells": [list(rc) for rc in self.relative_cells],
            "grid_xy": [list(xy) for xy in self.grid_xy],
            "tangent": list(self.tangent),
            "outward_dir": list(self.outward_dir),
            "requested_center": list(self.requested_center),
            "span_a": float(self.span_a),
            "n_samples": self.n_samples,
        }


def _unit(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float).reshape(2)
    n = float(np.linalg.norm(v))
    if n <= 0:
        raise ValueError("zero tangent/outward vector")
    return v / n


def cell_to_grid_xy(ix: int, iy: int, *, res: int) -> Tuple[float, float]:
    return (float((ix + 0.5) / res), float((iy + 0.5) / res))


def build_grid_index_port_line(
    port_index: int,
    *,
    res: int,
    center_xy: np.ndarray,
    outward_dir: np.ndarray,
    span_a: float,
    n_target: int = 31,
) -> GridIndexPortLine:
    """
    Build a cross-section by snapping a continuous chord to unique Hz Yee DOFs.

    The anchor is the snapped monitor/source center; samples walk along the
    tangent with sub-cell spacing, merging duplicate snaps (grid-first DOFs).
    """
    center_xy = np.asarray(center_xy, dtype=float).reshape(2)
    u = _unit(outward_dir)
    tangent = np.array([-u[1], u[0]])
    anchor_site = hz_yee_site(center_xy, res=res)
    ix0, iy0 = anchor_site.cell_index

    if n_target <= 1:
        offsets = np.array([0.0])
    else:
        offsets = np.linspace(-span_a / 2, span_a / 2, n_target)

    merged: Dict[Tuple[int, int], Tuple[float, float]] = {}
    order: List[Tuple[int, int]] = []
    for s in offsets:
        xy = center_xy + s * tangent
        site = hz_yee_site(xy, res=res)
        key = site.cell_index
        if key not in merged:
            merged[key] = (float(site.grid_xy[0]), float(site.grid_xy[1]))
            order.append(key)

    rel = [(ix - ix0, iy - iy0) for ix, iy in order]
    grid = [merged[k] for k in order]
    return GridIndexPortLine(
        port_index=int(port_index),
        res=int(res),
        anchor_cell=(int(ix0), int(iy0)),
        relative_cells=rel,
        grid_xy=grid,
        tangent=(float(tangent[0]), float(tangent[1])),
        outward_dir=(float(u[0]), float(u[1])),
        requested_center=(float(center_xy[0]), float(center_xy[1])),
        span_a=float(span_a),
    )


def audit_port_line_vs_continuous(
    line: GridIndexPortLine,
    *,
    continuous_offsets: Sequence[float],
) -> Dict[str, Any]:
    """Compare grid-index line to explicit continuous-coordinate snaps."""
    center = np.array(line.requested_center, dtype=float)
    tangent = np.array(line.tangent, dtype=float)
    res = line.res
    cont_cells: List[Tuple[int, int]] = []
    for s in continuous_offsets:
        site = hz_yee_site(center + s * tangent, res=res)
        cont_cells.append(site.cell_index)
    grid_set = set(line.relative_cells)
    rel_from_anchor = [
        (ix - line.anchor_cell[0], iy - line.anchor_cell[1])
        for ix, iy in cont_cells
    ]
    matched = sum(1 for rc in rel_from_anchor if rc in grid_set)
    return {
        "n_grid_index": line.n_samples,
        "n_continuous_snaps": len(cont_cells),
        "n_unique_continuous": len(set(cont_cells)),
        "matched_relative_cells": int(matched),
        "anchor_cell": list(line.anchor_cell),
    }
