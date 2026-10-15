#!/usr/bin/env python3
"""Phase 2: evaluate three starts under frozen direct 12×12 DtN S."""
from __future__ import annotations

import argparse
import time

import numpy as np

from multimode_opt_core import (
    B_PLUS,
    F_HZ,
    OUT,
    WEIGHTS,
    dump_json,
    eval_design,
    load_start,
    orbits_and_masks,
)


STARTS = [
    ("uniform", "uniform"),
    ("pref_peak", "opt_single_uniform_M_Bp0.0500"),
    ("pref_broadband", "opt_broadband_C_Bp0.0500"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--B", type=float, default=B_PLUS)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    orbits, _, _ = orbits_and_masks()
    rows = []
    t0 = time.time()
    for label, tag in STARTS:
        q, s, _ = load_start(tag, orbits)
        print(f"=== {label} grade={args.grade} ===", flush=True)
        bun = eval_design(args.grade, s, F_HZ, args.B, with_absorption=True)
        met = bun["metrics"]
        np.savez_compressed(
            OUT / f"S12_baseline_{label}_{args.grade}_Bp.npz",
            S=bun["S"],
            q=q,
            s=s,
            label=label,
            design_tag=tag,
        )
        rec = {
            "label": label,
            "design_tag": tag,
            "grade": args.grade,
            "B": args.B,
            "f_Hz": F_HZ,
            "q": q.tolist(),
            "s": s.tolist(),
            "q_min": float(np.min(q)),
            "q_max": float(np.max(q)),
            "weights": WEIGHTS,
            "metrics": met,
        }
        dump_json(OUT / f"baseline_{label}_{args.grade}.json", rec)
        rows.append(rec)
        print(
            f"  J={met['J']:.6f} Tdes={met['T_des_avg']:.4f} IL={met['IL_avg_dB']:.2f} "
            f"iso={met['isolation_avg_dB']:.2f} RL={met['RL_avg_dB']:.2f} "
            f"TE1={met['TE1_avg']:.4f} Pmat={met.get('P_material_avg', float('nan')):.4f}",
            flush=True,
        )
    summary = {
        "grade": args.grade,
        "B": args.B,
        "f_Hz": F_HZ,
        "weights": WEIGHTS,
        "starts": [
            {
                "label": r["label"],
                "design_tag": r["design_tag"],
                "J": r["metrics"]["J"],
                "T_des_avg": r["metrics"]["T_des_avg"],
                "T_des_min": r["metrics"]["T_des_min"],
                "IL_avg_dB": r["metrics"]["IL_avg_dB"],
                "IL_worst_dB": r["metrics"]["IL_worst_dB"],
                "isolation_avg_dB": r["metrics"]["isolation_avg_dB"],
                "isolation_min_dB": r["metrics"]["isolation_min_dB"],
                "RL_avg_dB": r["metrics"]["RL_avg_dB"],
                "RL_min_dB": r["metrics"]["RL_min_dB"],
                "T_thr_avg": r["metrics"]["T_thr_avg"],
                "T_other_avg": r["metrics"]["T_other_avg"],
                "R_avg": r["metrics"]["R_avg"],
                "TE1_avg": r["metrics"]["TE1_avg"],
                "P_material_avg": r["metrics"].get("P_material_avg"),
                "P_out_avg": r["metrics"]["P_out_avg"],
                "q_min": r["q_min"],
                "q_max": r["q_max"],
            }
            for r in rows
        ],
        "elapsed_s": time.time() - t0,
    }
    dump_json(OUT / f"baselines_summary_{args.grade}.json", summary)

    lines = [
        "# MULTIMODE_BASELINE_REPORT",
        "",
        f"Frozen direct DtN 12×12, measure=projected, f=3.85 GHz, B=+{args.B} T, grade={args.grade}.",
        f"Weights: {WEIGHTS}.",
        "",
        "## Aggregate (mean over 6 TE0 sources)",
        "",
        "| start | J | T_des | IL (dB) | iso (dB) | RL (dB) | thr | other | R | TE1 | P_mat |",
        "|-------|--:|------:|--------:|---------:|--------:|----:|------:|--:|----:|------:|",
    ]
    for r in summary["starts"]:
        lines.append(
            f"| {r['label']} | {r['J']:.4f} | {r['T_des_avg']:.4f} | {r['IL_avg_dB']:.2f} | "
            f"{r['isolation_avg_dB']:.2f} | {r['RL_avg_dB']:.2f} | {r['T_thr_avg']:.4f} | "
            f"{r['T_other_avg']:.4f} | {r['R_avg']:.4f} | {r['TE1_avg']:.4f} | "
            f"{r['P_material_avg']:.4f} |"
        )
    lines += [
        "",
        "Pref designs are warm starts only under the multimode S; they are **not** the final design.",
        "",
        f"Elapsed: {summary['elapsed_s']:.1f} s.",
        "",
    ]
    (OUT / "MULTIMODE_BASELINE_REPORT.md").write_text("\n".join(lines) + "\n")
    print("wrote", OUT / "MULTIMODE_BASELINE_REPORT.md", flush=True)


if __name__ == "__main__":
    main()
