#!/usr/bin/env python3
"""Numerically check whether C6/C3 rotationally related sources match predicted S.

For a C6-symmetric design, rotating the TE0 response of source 0 by 60° should
reproduce sources 1..5. For C3, two inequivalent classes may be needed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "jesse_method_campaign"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from eval_forward_only import eval_q_forward  # noqa: E402
from multimode_opt_core import B_PLUS, F_HZ, orbits_and_masks  # noqa: E402
from threefold_orbits import build_orbits_c3, expand_tied as expand_c3, map_c6_q_to_c3  # noqa: E402
from sixfold_orbits import expand_tied as expand_c6  # noqa: E402
import sixport_common as sc  # noqa: E402
from full91_geometry import bulb_centers_device  # noqa: E402


def rotate_port_index(p, k):
    return (p + k) % 6


def rotate_S_te0_columns(S, k):
    """Permute 12×12 S under 60°*k rotation of ports (TE0/TE1 pairs)."""
    Sp = np.zeros_like(S)
    for j in range(6):
        j2 = rotate_port_index(j, k)
        for p in range(6):
            p2 = rotate_port_index(p, k)
            for m in (0, 1):
                for n in (0, 1):
                    Sp[2 * p2 + m, 2 * j2 + n] = S[2 * p + m, 2 * j + n]
    return Sp


def te0_column_err(S, S_ref_rot, j):
    c = 2 * j
    return float(np.linalg.norm(S[:, c] - S_ref_rot[:, c]) / max(np.linalg.norm(S_ref_rot[:, c]), 1e-30))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--design", default="opt_uniform_C_warm_M")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    orbits6, _, _ = orbits_and_masks()
    p = ROOT / "outputs/circulator_inverse_design/multimode_optimization" / f"{args.design}.json"
    d = json.loads(p.read_text())
    q6 = np.asarray(d.get("q_best", d.get("q", d.get("q0"))), float)

    bun, _ = eval_q_forward(args.grade, q6, orbits6, f_hz=F_HZ, b_tesla=B_PLUS)
    S = bun["S"]
    c6_errs = []
    for k in range(1, 6):
        Sr = rotate_S_te0_columns(S, k)
        # Compare column k of S to rotated column 0
        err = te0_column_err(S, Sr, k)
        # Actually: S[:, 2k] should equal rotate_k(S)[:, 2k] built from source 0
        # Build predicted column k from source 0 by rotating ports
        pred = np.zeros(12, dtype=complex)
        for p0 in range(6):
            for m in (0, 1):
                pred[2 * rotate_port_index(p0, k) + m] = S[2 * p0 + m, 0]
        act = S[:, 2 * k]
        err = float(np.linalg.norm(act - pred) / max(np.linalg.norm(act), 1e-30))
        c6_errs.append({"k": k, "rel_err": err})

    # C3-mapped design: may break perfect C6; check 120° and two source classes
    centers = bulb_centers_device(sc)
    orbits3, _ = build_orbits_c3(centers)
    q3 = map_c6_q_to_c3(q6, orbits6, orbits3)
    # For a C6 design mapped to C3, s should match → same S. Use free C3 only if we had one.
    # Check 120° on the C6-symmetric S:
    c3_errs = []
    for k in (2, 4):  # 120°, 240°
        pred = np.zeros(12, dtype=complex)
        for p0 in range(6):
            for m in (0, 1):
                pred[2 * rotate_port_index(p0, k) + m] = S[2 * p0 + m, 0]
        act = S[:, 2 * k]
        err = float(np.linalg.norm(act - pred) / max(np.linalg.norm(act), 1e-30))
        c3_errs.append({"k": k, "rel_err": err})

    # Class check: source 0 vs source 1 under C3 (should differ if only C3 symmetry)
    # For C6 design they should still match via 60° — report both.
    report = {
        "design": args.design,
        "grade": args.grade,
        "C6_rotation_from_source0": c6_errs,
        "C6_max_rel_err": max(e["rel_err"] for e in c6_errs),
        "C6_pass_1e-6": max(e["rel_err"] for e in c6_errs) < 1e-6,
        "C3_subset_120deg_from_source0": c3_errs,
        "note": (
            "If C6_pass, one representative TE0 source determines all six by rotation "
            "for a C6-symmetric design (LU still dominates; RHS savings small). "
            "C3-only designs would need two inequivalent source classes {0} and {1}."
        ),
    }
    (OUT / "rotation_source_check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
