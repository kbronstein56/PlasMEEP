#!/usr/bin/env python3
"""Phase 7: re-evaluate best C6 design on a finer mesh (FULL91-F or fallback)."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from multimode_opt_core import F_HZ, B_PLUS, OUT, dump_json, eval_design


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--opt-json", required=True)
    ap.add_argument("--from-grade", default="M")
    ap.add_argument("--to-grade", default="F")
    ap.add_argument("--label", default="best")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    opt = json.loads(Path(args.opt_json).read_text())
    s = np.asarray(opt["best"]["s"], float)
    t0 = time.time()

    print(f"eval {args.from_grade} ...", flush=True)
    bun_lo = eval_design(args.from_grade, s, F_HZ, B_PLUS, with_absorption=True)
    print(f"eval {args.to_grade} ...", flush=True)
    try:
        bun_hi = eval_design(args.to_grade, s, F_HZ, B_PLUS, with_absorption=True)
        hi_ok = True
        hi_err = None
    except Exception as e:
        bun_hi = None
        hi_ok = False
        hi_err = repr(e)
        print(f"FINE MESH FAILED: {e}", flush=True)

    m_lo = bun_lo["metrics"]
    report = {
        "opt_json": str(args.opt_json),
        "label": args.label,
        "from_grade": args.from_grade,
        "to_grade": args.to_grade,
        "metrics_from": m_lo,
        "fine_ok": hi_ok,
        "fine_error": hi_err,
        "elapsed_s": time.time() - t0,
    }
    if hi_ok:
        m_hi = bun_hi["metrics"]
        report["metrics_to"] = m_hi
        keys = [
            "J",
            "T_des_avg",
            "IL_avg_dB",
            "isolation_avg_dB",
            "RL_avg_dB",
            "TE1_avg",
            "T_thr_avg",
            "T_other_avg",
            "R_avg",
            "P_material_avg",
        ]
        diffs = {}
        for k in keys:
            if k in m_lo and k in m_hi:
                diffs[k] = float(m_hi[k] - m_lo[k])
        report["deltas_to_minus_from"] = diffs
        np.savez_compressed(
            OUT / f"S12_finemesh_{args.label}_{args.to_grade}_Bp.npz",
            S=bun_hi["S"],
            s=s,
        )
    dump_json(OUT / f"fine_mesh_{args.label}_{args.from_grade}_to_{args.to_grade}.json", report)

    lines = [
        "# MULTIMODE_FINE_MESH_REPORT",
        "",
        f"Design `{args.opt_json}` transferred {args.from_grade} → {args.to_grade}.",
        "",
    ]
    if not hi_ok:
        lines += [
            f"**FINE MESH NOT REACHED:** {hi_err}",
            "",
            f"Production metrics remain those of grade {args.from_grade}; "
            "do not label them mesh-converged.",
            "",
            f"**MULTIMODE_FINE_MESH: NOT_REACHED**",
            "",
        ]
    else:
        lines += [
            "| metric | " + args.from_grade + " | " + args.to_grade + " | Δ |",
            "|--------|---:|---:|--:|",
        ]
        for k, d in report["deltas_to_minus_from"].items():
            lines.append(
                f"| {k} | {m_lo[k]:.6g} | {report['metrics_to'][k]:.6g} | {d:.4g} |"
            )
        lines += ["", f"**MULTIMODE_FINE_MESH: PASS** (evaluation completed)", ""]
    (OUT / "MULTIMODE_FINE_MESH_REPORT.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({k: report[k] for k in ("fine_ok", "elapsed_s") if k in report}, indent=2))


if __name__ == "__main__":
    main()
