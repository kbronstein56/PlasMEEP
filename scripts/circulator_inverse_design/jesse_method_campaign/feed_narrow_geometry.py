#!/usr/bin/env python3
"""Narrowed-feed geometry helpers for forward-only TE1-suppression tests.

Parallel-plate convention (matches modes.py / production DtN):
  clear aperture width W (lattice units), physical W_m = W * a
  TE_m cutoff: f_c,m = m * c / (2 W_m)   (m = 1,2,...)
  TE0 has no cutoff.

Production feed: W = clear_width = 2.0 lu = 40 mm → f_c,TE1 ≈ 3.747 GHz
  so TE1 propagates throughout the 4–6 GHz advisor band.

Narrowed feed target: f_c,TE1 = 6.5 GHz > Band-C upper edge (5.5 GHz)
  → W_m = c/(2*6.5e9) ≈ 23.061 mm → W ≈ 1.153 lu

Stage-6 NOT STARTED. Forward-only use.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from dtn_port_geometry import port_outer_center  # noqa: E402
from modes import verify_cutoffs_at_f  # noqa: E402

C0 = 299792458.0

# Production clear width (lattice units / meters)
W_ORIG_LU = float(sc.clear_width)
W_ORIG_M = W_ORIG_LU * float(sc.a)

# TE1 cutoff safely above Band C (4.5–5.5 GHz) and Band B (≤5.25 GHz)
FC_TE1_TARGET_HZ = 6.5e9
W_NARROW_M = C0 / (2.0 * FC_TE1_TARGET_HZ)
W_NARROW_LU = W_NARROW_M / float(sc.a)


def te_cutoffs(W_lu: float) -> dict:
    """Analytic cutoffs for parallel-plate PEC feed used by the solver."""
    W_m = float(W_lu) * float(sc.a)
    rows = []
    for m in range(0, 4):
        if m == 0:
            fc = 0.0
        else:
            fc = m * C0 / (2.0 * W_m)
        rows.append(
            {
                "m": m,
                "label": f"TE{m}" if m > 0 else "TE0",
                "f_cut_Hz": fc,
                "f_cut_GHz": fc / 1e9,
                "kc_lu": m * np.pi / float(W_lu),
            }
        )
    return {
        "W_lu": float(W_lu),
        "W_m": W_m,
        "W_mm": W_m * 1e3,
        "a_m": float(sc.a),
        "formula": "f_c,m = m * c / (2 W_m), W_m = W_lu * a; TE0 cutoff-free",
        "modes": rows,
        "TE1_cutoff_GHz": rows[1]["f_cut_GHz"],
    }


def feed_geometry_document() -> dict:
    orig = te_cutoffs(W_ORIG_LU)
    narrow = te_cutoffs(W_NARROW_LU)
    # sanity: mode_constants agree
    check_orig = verify_cutoffs_at_f(5.0e9, W_ORIG_LU)
    check_nar = verify_cutoffs_at_f(5.0e9, W_NARROW_LU)
    return {
        "stage6_status": "NOT STARTED / AWAITING AUTHORIZATION",
        "original_feed": orig,
        "narrowed_feed": narrow,
        "target": {
            "TE1_cutoff_target_GHz": FC_TE1_TARGET_HZ / 1e9,
            "rationale": (
                "Place TE1 cutoff at 6.5 GHz so TE1 is evanescent over Band B "
                "(≤5.25 GHz) and Band C (≤5.5 GHz) with margin."
            ),
        },
        "solver_check_at_5GHz": {
            "original_TE1_propagating": check_orig["modes"][1]["propagating"],
            "narrowed_TE1_propagating": check_nar["modes"][1]["propagating"],
            "original": check_orig,
            "narrowed": check_nar,
        },
    }


def narrowed_feed_wall_mask(
    points: np.ndarray,
    tris: np.ndarray,
    *,
    W_filter: float = W_NARROW_LU,
) -> np.ndarray:
    """OR-mask: metalize |transverse| > W_filter/2 along each straight feed.

    Covers the uniform feed from just outside the flare throat to the outer
    DtN plane, matching the production TE1-filter idea but over the full
    straight section so the port plane sees the narrowed aperture.
    """
    cents = points[tris].mean(axis=1)
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    W0 = float(sc.clear_width)
    half = 0.5 * float(W_filter)
    mask = np.zeros(len(tris), dtype=bool)
    for port in range(6):
        horn = sc.full_horns[port]
        n_hat = np.asarray(sc.port_dirs[port], float)
        n_hat /= np.linalg.norm(n_hat)
        tang = np.array([-n_hat[1], n_hat[0]])
        throat = np.asarray(horn["throat_center"], float) + shift
        outer = port_outer_center(horn, n_hat) + shift
        z_lo = float(throat @ n_hat) + 0.05  # just into straight feed
        z_hi = float(outer @ n_hat) + 0.05  # past outer plane slightly
        zeta = cents @ n_hat
        # transverse about throat axis (same as feed centerline)
        s = (cents - throat) @ tang
        in_feed = (zeta >= z_lo) & (zeta <= z_hi)
        in_side = (np.abs(s) >= half) & (np.abs(s) <= 0.5 * W0 + sc.wall_thickness + 1e-9)
        mask |= in_feed & in_side
    return mask


def write_geometry_doc(path: Path):
    doc = feed_geometry_document()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, indent=2) + "\n")
    return doc


if __name__ == "__main__":
    doc = feed_geometry_document()
    print(json.dumps(doc, indent=2))
