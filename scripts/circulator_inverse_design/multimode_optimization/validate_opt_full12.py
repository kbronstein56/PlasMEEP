#!/usr/bin/env python3
"""Phase 6: full 12×12 Onsager + power balance on best optimized design."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

from multimode_opt_core import (
    BALANCE_BAR,
    B_PLUS,
    F_HZ,
    ONSAGER_BAR,
    OUT,
    _channel_basis,
    dump_json,
    eval_design,
    te0_metrics_from_S,
)

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
]
from adjoint_direct_dtn import build_system  # noqa: E402
from boundary_ports import add_incident_rhs_from_c, modal_c_vector  # noqa: E402
from direct_dtn_s import onsager_metrics  # noqa: E402
from modes import power_scale  # noqa: E402
from validate_modal_s import material_absorption  # noqa: E402


def power_balance_all(st, lu, S):
    points, tris, rho, ports, modes, k0 = (
        st["points"],
        st["tris"],
        st["rho"],
        st["ports"],
        st["modes"],
        st["k0"],
    )
    rows = []
    ok = True
    for n in range(12):
        p, m = divmod(n, 2)
        idx, c = modal_c_vector(points, ports[p], m)
        scale = power_scale(modes[m]["beta"], k0)
        a_field = 1.0 / max(abs(scale), 1e-30)
        b = np.zeros(len(points), dtype=np.complex128)
        b = add_incident_rhs_from_c(b, idx, c, modes[m]["beta"], a_field)
        x = lu.solve(b)
        P_mat = float(material_absorption(points, tris, x, rho, k0))
        P_out = float(np.sum(np.abs(S[:, n]) ** 2))
        eng = abs(1.0 - P_out - P_mat)
        pass_n = eng <= BALANCE_BAR and P_mat >= -0.02
        ok = ok and pass_n
        rows.append(
            {
                "channel": n,
                "port": p,
                "mode": m,
                "P_out": P_out,
                "P_material": P_mat,
                "eng_residual": float(eng),
                "PASS": bool(pass_n),
            }
        )
    return rows, ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--opt-json", required=True)
    ap.add_argument("--grade", default=None)
    ap.add_argument("--label", default="best")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    opt = json.loads(Path(args.opt_json).read_text())
    grade = args.grade or opt["grade"]
    s = np.asarray(opt["best"]["s"], float)
    t0 = time.time()
    print(f"validating {args.opt_json} grade={grade}", flush=True)

    bun_p = eval_design(grade, s, F_HZ, +B_PLUS, with_absorption=True)
    bun_m = eval_design(grade, s, F_HZ, -B_PLUS, with_absorption=False)
    bun_0 = eval_design(grade, s, F_HZ, 0.0, with_absorption=False)
    Sp, Sm, S0 = bun_p["S"], bun_m["S"], bun_0["S"]
    ons = onsager_metrics(Sp, Sm)
    ons_pass = float(ons["r_F"]) <= ONSAGER_BAR

    st = build_system(grade, s, F_HZ, +B_PLUS, "projected")
    lu, cols, C, U, V, betas, S_chk = _channel_basis(st)
    bal_rows, bal_ok = power_balance_all(st, lu, Sp)

    np.savez_compressed(OUT / f"S12_validated_{args.label}_{grade}_Bp.npz", S=Sp, s=s)
    np.savez_compressed(OUT / f"S12_validated_{args.label}_{grade}_Bm.npz", S=Sm, s=s)
    np.savez_compressed(OUT / f"S12_validated_{args.label}_{grade}_B0.npz", S=S0, s=s)

    met_p = te0_metrics_from_S(Sp)
    met_m = te0_metrics_from_S(Sm)
    met_0 = te0_metrics_from_S(S0)
    report = {
        "opt_json": str(args.opt_json),
        "label": args.label,
        "grade": grade,
        "onsager": ons,
        "ONSAGER_BAR": ONSAGER_BAR,
        "OPTIMIZED_FULL12_ONSAGER": "PASS" if ons_pass else "FAIL",
        "power_balance": {
            "bar": BALANCE_BAR,
            "per_channel": bal_rows,
            "max_residual": float(max(r["eng_residual"] for r in bal_rows)),
            "OPTIMIZED_POWER_BALANCE": "PASS" if bal_ok else "FAIL",
        },
        "metrics_Bp": met_p,
        "metrics_Bm": met_m,
        "metrics_B0": met_0,
        "B_reversal": {
            "Bp_Tdes": met_p["T_des_avg"],
            "Bp_Trev": met_p["T_rev_avg"],
            "Bm_Tdes": met_m["T_des_avg"],
            "Bm_Trev": met_m["T_rev_avg"],
            "note": "At -B, CW/CCW swap: compare Bp Tdes ≈ Bm Trev",
            "swap_abs_err": float(abs(met_p["T_des_avg"] - met_m["T_rev_avg"])),
        },
        "elapsed_s": time.time() - t0,
    }
    dump_json(OUT / f"full12_validation_{args.label}_{grade}.json", report)

    lines = [
        "# MULTIMODE_OPT_FULL12_VALIDATION",
        "",
        f"Design: `{args.opt_json}` label={args.label} grade={grade}.",
        "",
        f"- Onsager r_F = {ons['r_F']:.3e} (bar {ONSAGER_BAR}) → "
        f"**{report['OPTIMIZED_FULL12_ONSAGER']}**",
        f"- Power balance max residual = "
        f"{report['power_balance']['max_residual']:.4f} (bar {BALANCE_BAR}) → "
        f"**{report['power_balance']['OPTIMIZED_POWER_BALANCE']}**",
        "",
        "## TE0 metrics",
        "",
        "| B | T_des | T_rev | IL | iso | RL | TE1 |",
        "|---|------:|------:|---:|----:|---:|----:|",
    ]
    for name, met in [("+0.05", met_p), ("−0.05", met_m), ("0", met_0)]:
        lines.append(
            f"| {name} | {met['T_des_avg']:.4f} | {met['T_rev_avg']:.4f} | "
            f"{met['IL_avg_dB']:.2f} | {met['isolation_avg_dB']:.2f} | "
            f"{met['RL_avg_dB']:.2f} | {met['TE1_avg']:.4f} |"
        )
    lines += [
        "",
        f"B-reversal |T_des(+B) − T_rev(−B)| = {report['B_reversal']['swap_abs_err']:.4e}",
        "",
        "## Power balance per channel",
        "",
        "| ch | port | mode | P_out | P_mat | residual | PASS |",
        "|---:|-----:|-----:|------:|------:|---------:|:----:|",
    ]
    for r in bal_rows:
        lines.append(
            f"| {r['channel']} | {r['port']} | {r['mode']} | {r['P_out']:.4f} | "
            f"{r['P_material']:.4f} | {r['eng_residual']:.4f} | {r['PASS']} |"
        )
    lines.append("")
    (OUT / "MULTIMODE_OPT_FULL12_VALIDATION.md").write_text("\n".join(lines) + "\n")
    print(
        json.dumps(
            {
                "OPTIMIZED_FULL12_ONSAGER": report["OPTIMIZED_FULL12_ONSAGER"],
                "r_F": report["onsager"]["r_F"],
                "OPTIMIZED_POWER_BALANCE": report["power_balance"]["OPTIMIZED_POWER_BALANCE"],
                "max_residual": report["power_balance"]["max_residual"],
                "elapsed_s": report["elapsed_s"],
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
