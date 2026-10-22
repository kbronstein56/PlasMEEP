#!/usr/bin/env python3
"""Evaluate Objectives A and B on existing designs (forward-only, single freq).

Writes MODE_OBJECTIVE_WEIGHTS_FROZEN.json before any comparison.
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
from mode_objectives import OBJECTIVE_FREEZE, write_freeze  # noqa: E402
from multimode_opt_core import B_PLUS, F_HZ, expand_tied, orbits_and_masks  # noqa: E402
from sixfold_orbits import build_orbits  # noqa: E402
import sixport_common as sc  # noqa: E402
from full91_geometry import bulb_centers_device  # noqa: E402


def load_design(name, orbits):
    if name == "uniform":
        return np.ones(len(orbits), float), "uniform"
    from design_io import load_design_q

    q = load_design_q(name)
    return q, name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument(
        "--designs",
        default="uniform,opt_uniform_C_warm_M,opt_uniform_C",
    )
    ap.add_argument("--f-ghz", type=float, default=F_HZ / 1e9)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    write_freeze(OUT)

    orbits, _, _ = orbits_and_masks()
    rows = []
    for name in args.designs.split(","):
        name = name.strip()
        try:
            q, tag = load_design(name, orbits)
        except FileNotFoundError:
            print("skip missing", name, flush=True)
            continue
        if len(q) != len(orbits):
            print(f"skip {name}: len(q)={len(q)} != {len(orbits)}", flush=True)
            continue
        print(f"eval {tag}...", flush=True)
        bun, s = eval_q_forward(args.grade, q, orbits, f_hz=args.f_ghz * 1e9, b_tesla=B_PLUS)
        mr = bun["mode_resolved"]
        rows.append(
            {
                "design": tag,
                "f_GHz": args.f_ghz,
                "J_A": mr["J_A"],
                "J_B": mr["J_B"],
                "des0": mr["des0_avg"],
                "des1": mr["des1_avg"],
                "des_tot": mr["des_tot_avg"],
                "rev0": mr["rev0_avg"],
                "rev1": mr["rev1_avg"],
                "rev_tot": mr["rev_tot_avg"],
                "iso_te0_dB": mr["iso_te0_avg_dB"],
                "iso_tot_dB": mr["iso_tot_avg_dB"],
                "TE1_frac_of_des": mr["des1_avg"] / max(mr["des_tot_avg"], 1e-30),
            }
        )
        print(
            f"  JA={mr['J_A']:.4f} JB={mr['J_B']:.4f} des0={mr['des0_avg']:.4f} "
            f"des1={mr['des1_avg']:.4f} isoTE0={mr['iso_te0_avg_dB']:.2f}",
            flush=True,
        )

    report = {"weights_frozen": OBJECTIVE_FREEZE, "grade": args.grade, "rows": rows}
    (OUT / "mode_baselines.json").write_text(json.dumps(report, indent=2) + "\n")
    print("wrote mode_baselines.json", flush=True)


if __name__ == "__main__":
    main()
