#!/usr/bin/env python3
"""Phase 5: evaluate best M design on F and apply frozen MESH_TRANSFER_CRITERION."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "multimode_optimization"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from efficient_adjoint import eval_and_grad_efficient  # noqa: E402
from multimode_opt_core import B_PLUS, F_HZ, dump_json, orbits_and_masks, te0_absorption_mean  # noqa: E402

# Frozen BEFORE judging — see MESH_TRANSFER_CRITERION.md
BARS = {
    "J": 0.05,
    "T_des_avg": 0.04,
    "IL_avg_dB": 0.40,
    "isolation_avg_dB": 1.5,
    "RL_avg_dB": 1.5,
    "TE1_avg": 0.015,
    "T_thr_avg": 0.010,
    "T_other_avg": 0.030,
}


def metrics_bundle(grade, q, orbits):
    bun, g, s = eval_and_grad_efficient(grade, q, orbits, f_hz=F_HZ, b_tesla=B_PLUS)
    met = dict(bun["metrics"])
    met["P_material_avg"] = te0_absorption_mean(bun["st"], bun["lu"])
    met["rem_avg"] = float(1.0 - met["P_out_avg"] - met["P_material_avg"])
    return {
        "J": float(bun["J"]),
        "metrics": met,
        "s": np.asarray(s, float).tolist(),
        "q": np.asarray(q, float).tolist(),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--opt-json", required=True, help="M-optimized design JSON")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    d = json.loads(Path(args.opt_json).read_text())
    q = np.asarray(d["best"]["q"], float)
    orbits, _, _ = orbits_and_masks()

    t0 = time.time()
    print("eval M ...", flush=True)
    m_rec = metrics_bundle("M", q, orbits)
    print(f"  M J={m_rec['J']:.4f} Tdes={m_rec['metrics']['T_des_avg']:.4f}", flush=True)
    print("eval F ...", flush=True)
    f_rec = metrics_bundle("F", q, orbits)
    print(f"  F J={f_rec['J']:.4f} Tdes={f_rec['metrics']['T_des_avg']:.4f}", flush=True)

    keys = list(BARS.keys())
    deltas = {}
    pass_flags = {}
    for k in keys:
        if k == "J":
            mv, fv = m_rec["J"], f_rec["J"]
        else:
            mv, fv = m_rec["metrics"][k], f_rec["metrics"][k]
        deltas[k] = float(fv - mv)
        pass_flags[k] = abs(deltas[k]) <= BARS[k]

    transfer_ok = all(pass_flags.values())
    report = {
        "opt_json": args.opt_json,
        "M": m_rec,
        "F": f_rec,
        "deltas_F_minus_M": deltas,
        "bars": BARS,
        "pass_per_metric": pass_flags,
        "MESH_TRANSFER": "PASS" if transfer_ok else "FAIL",
        "MESH_CONVERGENCE": "PASS" if transfer_ok else "FAIL",
        "elapsed_s": time.time() - t0,
    }
    dump_json(OUT / "m_to_f_transfer_report.json", report)

    lines = [
        "# M→F transfer (M-optimized design)",
        "",
        f"Design: `{args.opt_json}`",
        f"Criterion: `MESH_TRANSFER_CRITERION.md` (frozen before judgment).",
        "",
        f"**MESH_TRANSFER: {report['MESH_TRANSFER']}**",
        f"**MESH_CONVERGENCE: {report['MESH_CONVERGENCE']}**",
        "",
        "| metric | M | F | Δ(F−M) | bar | pass |",
        "|--------|--:|--:|-------:|----:|:----:|",
    ]
    for k in keys:
        if k == "J":
            mv, fv = m_rec["J"], f_rec["J"]
        else:
            mv, fv = m_rec["metrics"][k], f_rec["metrics"][k]
        lines.append(
            f"| {k} | {mv:.4f} | {fv:.4f} | {deltas[k]:+.4f} | {BARS[k]} | "
            f"{'Y' if pass_flags[k] else 'N'} |"
        )
    if not transfer_ok:
        lines += [
            "",
            "M→F remains large → do **not** call mesh-converged.",
            "Investigate local interface refinement / finer reference rather than",
            "only more optimizer iterations.",
        ]
    (OUT / "M_TO_F_TRANSFER_REPORT.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"MESH_TRANSFER": report["MESH_TRANSFER"], "deltas": deltas}, indent=2))


if __name__ == "__main__":
    main()
