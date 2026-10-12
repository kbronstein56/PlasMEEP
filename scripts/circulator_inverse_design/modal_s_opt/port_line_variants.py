#!/usr/bin/env python3
"""Monitor/source line variants for modal-S extraction audits."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sixport_common as sc  # noqa: E402


def make_port_lines(
    n_points: int = 41,
    monitor_shift_outward: float = 0.0,
    span_factor: float = 0.92,
) -> list:
    """Build six port lines; shift monitors along +outward (away from device).

    Positive shift moves monitor toward the source (more outward).
    Keep monitors between source and flare: nominal mon is 1.2 inward of src,
    so shift ∈ [-0.4, +0.4] stays in the uniform feed for production geometry.
    """
    lines = []
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    span = span_factor * sc.clear_width
    offsets = np.linspace(-span / 2.0, span / 2.0, n_points)
    ds = span / max(n_points - 1, 1)
    weights = np.full(n_points, ds, float)
    if n_points > 1:
        weights[0] *= 0.5
        weights[-1] *= 0.5
    else:
        weights[:] = span
    for port in range(6):
        horn = sc.full_horns[port]
        n_hat = np.asarray(sc.port_dirs[port], float)
        n_hat /= np.linalg.norm(n_hat)
        tang = np.array([-n_hat[1], n_hat[0]])
        src = np.asarray(horn["source_center"], float) + shift
        mon = np.asarray(horn["monitor_center"], float) + shift
        mon = mon + float(monitor_shift_outward) * n_hat
        lines.append({
            "port": port,
            "outward": n_hat,
            "tangent": tang,
            "source_xy": np.vstack([src + float(s) * tang for s in offsets]),
            "monitor_xy": np.vstack([mon + float(s) * tang for s in offsets]),
            "weights": weights.copy(),
            "offsets": offsets.copy(),
            "monitor_shift_outward": float(monitor_shift_outward),
            "n_points": int(n_points),
            "span": float(span),
        })
    return lines


def te1_profile(offsets: np.ndarray, W: float) -> np.ndarray:
    """L2-normalized odd TE1 (Neumann n=1): φ ∝ sin(π s / W) on [-W/2,W/2]."""
    phi = np.sin(np.pi * offsets / W)
    # renormalize with trapezoid later by caller; return shape only
    return phi.astype(float)


def mode_beta(k0: float, n: int, W: float) -> complex:
    kc = n * np.pi / W
    disc = k0**2 - kc**2
    if disc >= 0:
        return complex(np.sqrt(disc), 0.0)
    return complex(0.0, np.sqrt(-disc))
