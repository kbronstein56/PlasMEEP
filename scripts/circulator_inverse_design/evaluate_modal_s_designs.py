#!/usr/bin/env python3
"""Re-evaluate existing circulator designs with true modal S-parameters.

Does NOT re-optimize. Compares conventional |S|^2 metrics against the
frozen-P_ref campaign quantities (t = P/P_ref), and corrects the misuse of
the word \"insertion\" for 10 log10(P/P_ref).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from circulator_core import F_HZ, load_device, metrics_from_power, solve_six  # noqa: E402
from modal_s_core import (  # noqa: E402
    circulator_s_metrics,
    complex_matrix_to_json,
    dump_json,
    modal_port_powers,
    onsager_error,
    reciprocity_error,
    scattering_matrix,
)
from validate_modal_s import material_absorption  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design" / "modal_s"
DESIGNS = ROOT / "outputs" / "circulator_inverse_design"


def load_design(path: Path) -> dict:
    d = json.loads(path.read_text())
    best = d.get("best", d)
    s = np.asarray(best["s"], float)
    return {
        "tag": path.stem,
        "path": str(path),
        "s": s,
        "P_ref_frozen": float(d.get("P_ref", best.get("metrics", {}).get("P_ref", 0.0))),
        "pref_metrics": best.get("metrics", {}),
        "pref_J": float(best.get("J", best.get("metrics", {}).get("J", np.nan))),
    }


def eval_one(dev, s, f_hz, B, P_ref_frozen):
    bun = solve_six(dev, s, f_hz, B, P_ref=P_ref_frozen)
    A, Bmat, details = modal_port_powers(dev.sampler, bun.fields, bun.rho, bun.k0, dev.lines)
    S = scattering_matrix(A, Bmat)
    sm = circulator_s_metrics(S, A, Bmat)
    # Material absorption per source
    abs_mat = [
        material_absorption(dev.points, dev.tris, bun.fields[j], bun.rho, bun.k0) for j in range(6)
    ]
    sm["P_absorbed_material_avg"] = float(np.mean(abs_mat))
    sm["P_absorbed_material_per_source"] = abs_mat
    pref = bun.metrics
    return {
        "B_T": B,
        "f_Hz": f_hz,
        "k0": bun.k0,
        "factor_s": bun.factor_s,
        "solve_s": bun.solve_s,
        "S": complex_matrix_to_json(S),
        "A_diag": [[A[j, j].real, A[j, j].imag] for j in range(6)],
        "modal_metrics": sm,
        "pref_metrics": {
            "desired_avg": pref["desired_avg"],
            "reverse_avg": pref["reverse_avg"],
            "through_avg": pref["through_avg"],
            "leak_avg": pref["leak_avg"],
            "accepted_avg": pref["accepted_avg"],
            "isolation_avg_dB": pref["isolation_avg_dB"],
            "pref_gain_dB_misnamed_insertion": pref["insertion_avg_dB"],
            "J": pref["J"],
            "P_ref": pref["P_ref"],
            "note": (
                "pref_gain_dB_misnamed_insertion = 10 log10(t_des) with t=P/P_ref. "
                "This is NOT conventional insertion loss."
            ),
        },
        "comparison": {
            "T_des_modal": sm["T_des_avg"],
            "t_des_pref": pref["desired_avg"],
            "IL_dB_conventional": sm["IL_dB_avg"],
            "pref_gain_dB": pref["insertion_avg_dB"],
            "isolation_dB_modal": sm["isolation_dB_avg"],
            "isolation_dB_pref": pref["isolation_avg_dB"],
            "return_loss_dB": sm["return_loss_dB_avg"],
            "T_jj": sm["T_jj_avg"],
            "T_thr": sm["T_thr_avg"],
            "T_leak": sm["T_leak_avg"],
            "P_abs_modal_avg": sm["P_absorbed_modal_avg"],
            "P_abs_material_avg": sm["P_absorbed_material_avg"],
            "P_inc_avg": sm["P_inc_avg"],
            "accepted_over_Pref": pref["accepted_avg"],
        },
        "port_details": details,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="M")
    ap.add_argument("--f-hz", type=float, default=F_HZ)
    ap.add_argument("--B", type=float, default=0.05)
    ap.add_argument(
        "--designs",
        default="opt_single_uniform_M_Bp0.0500,opt_broadband_C_Bp0.0500,uniform",
        help="Comma tags: json stems under outputs/... or 'uniform'",
    )
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    print(f"Loading device grade={args.grade}", flush=True)
    t0 = time.perf_counter()
    dev = load_device(args.grade)
    print(f"  device ready in {time.perf_counter()-t0:.1f}s", flush=True)

    summary = {"grade": args.grade, "f_Hz": args.f_hz, "B_T": args.B, "designs": {}}

    for tag in args.designs.split(","):
        tag = tag.strip()
        if tag == "uniform":
            design = {
                "tag": "uniform_s1",
                "s": np.ones(91),
                "P_ref_frozen": None,
                "pref_metrics": {},
                "pref_J": None,
            }
            # freeze Pref from uniform accepted power
            bun0 = solve_six(dev, design["s"], args.f_hz, args.B, P_ref=None)
            design["P_ref_frozen"] = bun0.P_ref
        else:
            path = DESIGNS / f"{tag}.json"
            design = load_design(path)

        print(f"=== Design {design['tag']} ===", flush=True)
        rec = {"tag": design["tag"], "P_ref_frozen": design["P_ref_frozen"], "bias": {}}

        for B, label in [(args.B, "plusB"), (-args.B, "minusB"), (0.0, "B0")]:
            print(f"  solve B={B} ...", flush=True)
            t1 = time.perf_counter()
            ev = eval_one(dev, design["s"], args.f_hz, B, design["P_ref_frozen"])
            ev["wall_s"] = time.perf_counter() - t1
            rec["bias"][label] = ev
            c = ev["comparison"]
            print(
                f"    T_des={c['T_des_modal']:.4f} IL={c['IL_dB_conventional']:.2f} dB "
                f"iso={c['isolation_dB_modal']:.2f} dB RL={c['return_loss_dB']:.2f} dB "
                f"Tjj={c['T_jj']:.4f} thr={c['T_thr']:.4f} leak={c['T_leak']:.4f} "
                f"Pabs={c['P_abs_modal_avg']:.4e} (pref t_des={c['t_des_pref']:.4f} "
                f"pref_gain={c['pref_gain_dB']:.2f} dB)",
                flush=True,
            )

        Sp = np.array([[complex(*z) for z in row] for row in rec["bias"]["plusB"]["S"]])
        Sm = np.array([[complex(*z) for z in row] for row in rec["bias"]["minusB"]["S"]])
        S0 = np.array([[complex(*z) for z in row] for row in rec["bias"]["B0"]["S"]])
        rec["onsager_rel_F"] = onsager_error(Sp, Sm)
        rec["reciprocity_B0_rel_F"] = reciprocity_error(S0)
        # Circulator strength judgment under true S
        c = rec["bias"]["plusB"]["comparison"]
        strong = (
            c["T_des_modal"] > 0.5
            and c["isolation_dB_modal"] > 10.0
            and c["T_thr"] < 0.1
            and c["T_jj"] < 0.3
        )
        rec["circulator_judgment"] = {
            "strong_under_modal_S": bool(strong),
            "criteria": "T_des>0.5, iso>10 dB, T_thr<0.1, T_jj<0.3",
            "plusB": c,
        }
        dump_json(OUT / f"eval_{design['tag']}_{args.grade}.json", rec)
        summary["designs"][design["tag"]] = {
            "onsager_rel_F": rec["onsager_rel_F"],
            "reciprocity_B0_rel_F": rec["reciprocity_B0_rel_F"],
            "circulator_judgment": rec["circulator_judgment"],
            "plusB_comparison": c,
        }
        print(
            f"  Onsager={rec['onsager_rel_F']:.3e} B0-recip={rec['reciprocity_B0_rel_F']:.3e} "
            f"strong={strong}",
            flush=True,
        )

    dump_json(OUT / f"DESIGN_EVAL_SUMMARY_{args.grade}.json", summary)
    print("Wrote", OUT / f"DESIGN_EVAL_SUMMARY_{args.grade}.json", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
