#!/usr/bin/env python3
"""Forward-only dense local audit of a feasibility-band optimum.

Uses 3.125 MHz (and optional 1.5625 MHz around A2<=-5 minima) within the
optimized band. Stage-6 NOT STARTED.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
FEAS = OUT / "feasibility"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "jesse_method_campaign"),
]

from audit_bb_stage5_dense import freq_grid, parallel_scan  # noqa: E402
from multimode_opt_core import B_PLUS  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--opt-json", required=True)
    ap.add_argument("--grade", default="C")
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--omp-threads", type=int, default=1)
    ap.add_argument("--df-MHz", type=float, default=3.125)
    ap.add_argument("--b-tesla", type=float, default=None)
    args = ap.parse_args()

    for k, v in (
        ("OMP_NUM_THREADS", str(args.omp_threads)),
        ("OPENBLAS_NUM_THREADS", str(args.omp_threads)),
        ("MKL_NUM_THREADS", str(args.omp_threads)),
        ("PYTHONUNBUFFERED", "1"),
    ):
        os.environ.setdefault(k, v)

    opt_path = Path(args.opt_json)
    d = json.loads(opt_path.read_text())
    q = np.asarray(d["q_best"], float)
    f0, f1 = float(d["f0_GHz"]), float(d["f1_GHz"])
    df = args.df_MHz * 1e-3
    freqs = freq_grid(f0, f1, df)
    tag = f"audit_{opt_path.stem}_{args.df_MHz:g}MHz"
    resume = FEAS / f"{tag}.json"
    print(
        f"FEASIBILITY AUDIT {opt_path.name} band=[{f0},{f1}] df={args.df_MHz} MHz "
        f"n={len(freqs)} — Stage-6 NOT STARTED",
        flush=True,
    )
    # Note: parallel_scan uses B_PLUS via eval_q_forward; b override not wired here.
    # For B≠0.05 audits, use a dedicated path later.
    if args.b_tesla is not None and abs(args.b_tesla - B_PLUS) > 1e-12:
        raise SystemExit(
            "feasibility_band_audit currently supports production B_PLUS only; "
            "extend eval_forward_only path for other B."
        )
    out = parallel_scan(
        args.grade,
        q,
        freqs,
        tag=tag,
        resume_path=resume,
        seed_paths=[],
        n_workers=args.n_workers,
        omp_threads=args.omp_threads,
    )
    rows = out["rows"]
    worst = min(rows, key=lambda r: r["J_A2"])
    summary = {
        "opt": str(opt_path),
        "J_BB_best": d.get("J_BB_best"),
        "band": [f0, f1],
        "df_MHz": args.df_MHz,
        "n_freq": len(rows),
        "worst_A2": worst,
        "n_A2_le_m5": sum(1 for r in rows if r["J_A2"] <= -5),
        "n_A2_le_m7": sum(1 for r in rows if r["J_A2"] <= -7),
        "stage6_status": "NOT STARTED / AWAITING AUTHORIZATION",
    }
    (FEAS / f"{tag}_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        f"TRUE WORST A2={worst['J_A2']:.4f} @{worst['f_GHz']:.5f} "
        f"IL={worst['IL_dB']:.2f} iso={worst['isolation_dB']:.2f}",
        flush=True,
    )
    print("STAGE-6: NOT STARTED / AWAITING AUTHORIZATION", flush=True)


if __name__ == "__main__":
    main()
