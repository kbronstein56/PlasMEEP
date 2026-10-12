#!/usr/bin/env python3
"""Phase 3: modal-S adjoint FD at nonuniform Pref designs (full quotient rule)."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "modal_s_opt"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from circulator_core import F_HZ, B_PLUS, load_device, solve_six  # noqa: E402
from modal_s_adjoint import dJ_ds_selected, modal_J  # noqa: E402
from modal_s_core import dump_json  # noqa: E402
from sixfold_orbits import build_orbits  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design" / "modal_s_optimization"
DESIGNS = ROOT / "outputs" / "circulator_inverse_design"

# Predeclared gate (same as uniform ~1e-4 class; allow 5e-4 for nonuniform)
REL_ERR_MAX = 5e-4


def load_s(tag: str) -> np.ndarray:
    d = json.loads((DESIGNS / f"{tag}.json").read_text())
    return np.asarray(d["best"]["s"], float)


def pick_rods(centers: np.ndarray) -> dict:
    """Select representative rods: center, inner, outer, near-horn, symmetry pair."""
    c = np.asarray(centers, float)
    mid = c.mean(axis=0)
    r = np.linalg.norm(c - mid, axis=1)
    # near horn 0: outward +x
    horn0 = np.array([1.0, 0.0])
    along = c @ horn0
    rods = {
        "center": int(np.argmin(r)),
        "inner_ring": int(np.argsort(r)[len(r) // 4]),
        "outer_ring": int(np.argmax(r)),
        "near_horn0": int(np.argmax(along)),
    }
    # symmetry-related pair: C6 of near_horn0
    orbits, _origin = build_orbits(c)
    nh = rods["near_horn0"]
    for o in orbits:
        if nh in o["rods"] and len(o["rods"]) >= 2:
            rods["sym_pair_a"] = int(o["rods"][0])
            rods["sym_pair_b"] = int(o["rods"][1])
            break
    else:
        rods["sym_pair_a"] = nh
        rods["sym_pair_b"] = int((nh + 1) % 91)
    return rods


def fd_sweep(dev, s0, f_hz, B, rod: int, hs: list[float]) -> dict:
    bun0 = solve_six(dev, s0, f_hz, B, P_ref=1.0)
    adj = dJ_ds_selected(dev, bun0, [rod], freeze_a=False)
    g = adj["grads"][rod]
    rows = []
    best = None
    for h in hs:
        sp = s0.copy()
        sm = s0.copy()
        sp[rod] = float(np.clip(s0[rod] + h, 0.05, 1.80))
        sm[rod] = float(np.clip(s0[rod] - h, 0.05, 1.80))
        Jp = modal_J(dev, solve_six(dev, sp, f_hz, B, P_ref=1.0))
        Jm = modal_J(dev, solve_six(dev, sm, f_hz, B, P_ref=1.0))
        denom = sp[rod] - sm[rod]
        fd = (Jp - Jm) / denom if abs(denom) > 1e-18 else np.nan
        rel = abs(g - fd) / max(abs(fd), 1e-30)
        rec = {"h": h, "adjoint": g, "fd": float(fd), "rel_err": float(rel), "Jp": Jp, "Jm": Jm}
        rows.append(rec)
        if best is None or rel < best["rel_err"]:
            best = rec
    return {"rod": rod, "adjoint": g, "h_sweep": rows, "best": best}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--f-hz", type=float, default=F_HZ)
    ap.add_argument("--B", type=float, default=B_PLUS)
    ap.add_argument("--designs", default="opt_single_uniform_M_Bp0.0500,opt_broadband_C_Bp0.0500")
    ap.add_argument("--hs", default="1e-3,1e-4,1e-5,3e-6")
    args = ap.parse_args()
    hs = [float(x) for x in args.hs.split(",")]
    OUT.mkdir(parents=True, exist_ok=True)

    print(f"Loading grade={args.grade}", flush=True)
    dev = load_device(args.grade)
    rod_map = pick_rods(dev.centers)
    print("rods", rod_map, flush=True)

    summary = {
        "grade": args.grade,
        "REL_ERR_MAX": REL_ERR_MAX,
        "hs": hs,
        "rod_map": rod_map,
        "designs": {},
        "PASS": True,
    }

    for tag in args.designs.split(","):
        tag = tag.strip()
        print(f"=== design {tag} ===", flush=True)
        s0 = load_s(tag)
        t0 = time.perf_counter()
        # unique rods
        rods = sorted(set(rod_map.values()))
        rod_results = {}
        max_rel = 0.0
        for rid in rods:
            print(f"  rod {rid} h-sweep...", flush=True)
            rr = fd_sweep(dev, s0, args.f_hz, args.B, rid, hs)
            rod_results[str(rid)] = rr
            max_rel = max(max_rel, rr["best"]["rel_err"])
            print(f"    best rel={rr['best']['rel_err']:.3e} h={rr['best']['h']} "
                  f"adj={rr['best']['adjoint']:.6f} fd={rr['best']['fd']:.6f}", flush=True)
        # symmetry pair consistency of adjoint (same orbit should match under C6 design)
        sym = None
        if "sym_pair_a" in rod_map and "sym_pair_b" in rod_map:
            a, b = rod_map["sym_pair_a"], rod_map["sym_pair_b"]
            ga = rod_results[str(a)]["adjoint"]
            gb = rod_results[str(b)]["adjoint"]
            # For C6-tied peak design, orbit mates should have equal s and similar |grad|
            sym = {
                "a": a, "b": b,
                "s_a": float(s0[a]), "s_b": float(s0[b]),
                "g_a": ga, "g_b": gb,
                "rel_grad_diff": abs(ga - gb) / max(abs(ga), abs(gb), 1e-30),
            }
        passed = max_rel <= REL_ERR_MAX
        summary["designs"][tag] = {
            "PASS": bool(passed),
            "max_best_rel_err": max_rel,
            "rods": rod_results,
            "symmetry": sym,
            "wall_s": time.perf_counter() - t0,
            "s_range": [float(s0.min()), float(s0.max())],
        }
        summary["PASS"] = summary["PASS"] and passed
        dump_json(OUT / f"grad_valid_{tag}_{args.grade}.json", summary["designs"][tag])
        print(f"  design PASS={passed} max_rel={max_rel:.3e}", flush=True)

    dump_json(OUT / f"grad_valid_summary_{args.grade}.json", summary)
    print("OVERALL", "PASS" if summary["PASS"] else "FAIL", flush=True)
    return 0 if summary["PASS"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
