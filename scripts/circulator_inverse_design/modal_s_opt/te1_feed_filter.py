#!/usr/bin/env python3
"""Add a slight TE1-cutoff constriction in each feed between monitor and flare.

Local clear width W_filter < π/k0 ≈ 1.947 at 3.85 GHz makes TE1 evanescent
before the monitor, so fundamental-only S more closely matches a closed
single-mode port. Constriction is symmetric about the feed axis.
"""
from __future__ import annotations

import numpy as np

import sixport_common as sc


def te1_filter_wall_mask(
    points: np.ndarray,
    tris: np.ndarray,
    *,
    W_filter: float = 1.90,
    # Along outward: monitor at mon·n; flare is more inward (smaller n·x).
    # Place filter between monitor and cavity: ζ in [ζ_in, ζ_out] with
    # ζ = x·n̂, ζ_out = mon·n̂ - 0.15, ζ_in = mon·n̂ - 0.55.
    depth_outward_from_mon: float = -0.15,
    depth_inward_span: float = 0.45,
) -> np.ndarray:
    """Return triangle mask to OR into walls (ρ→0)."""
    cents = points[tris].mean(axis=1)
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    W0 = float(sc.clear_width)
    half_gap = 0.5 * W_filter
    # wall occupies |s| > half_gap where s is transverse offset
    mask = np.zeros(len(tris), dtype=bool)
    for port in range(6):
        horn = sc.full_horns[port]
        n_hat = np.asarray(sc.port_dirs[port], float)
        n_hat /= np.linalg.norm(n_hat)
        tang = np.array([-n_hat[1], n_hat[0]])
        mon = np.asarray(horn["monitor_center"], float) + shift
        zeta_mon = float(mon @ n_hat)
        # more inward = smaller ζ
        z_hi = zeta_mon + depth_outward_from_mon  # e.g. mon - 0.15
        z_lo = z_hi - depth_inward_span
        zeta = cents @ n_hat
        s = (cents - mon) @ tang  # transverse about monitor center line
        in_slab = (zeta >= z_lo) & (zeta <= z_hi)
        # add metal where |s| > half_gap but |s| < W0/2 + wall (inside feed)
        in_side = (np.abs(s) >= half_gap) & (np.abs(s) <= 0.5 * W0 + sc.wall_thickness)
        mask |= in_slab & in_side
    return mask
