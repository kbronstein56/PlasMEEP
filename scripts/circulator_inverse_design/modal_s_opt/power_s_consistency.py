#!/usr/bin/env python3
"""Phase 2: modal S power balance with material vs PML loss separated."""
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
from modal_s_core import dump_json, modal_port_powers, port_roles, scattering_matrix  # noqa: E402
from validate_modal_s import material_absorption  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design" / "modal_s_optimization"
DESIGNS = ROOT / "outputs" / "circulator_inverse_design"

# Predeclared residual gate: |1 - Σ|S|² - P_mat/P_inc - P_PML/P_inc| ≤ 0.05
# with P_PML := P_inc - Σ|b|² - P_mat  (nonnegative expected)
BALANCE_RESIDUAL_MAX = 0.05


def load_s(tag: str) -> np.ndarray:
    if tag == "uniform":
        return np.ones(91)
    d = json.loads((DESIGNS / f"{tag}.json").read_text())
    return np.asarray(d["best"]["s"], float)


def analyze(dev, s, f_hz, B, tag: str) -> dict:
    bun = solve_six(dev, s, f_hz, B, P_ref=1.0)
    A, Bmat, details = modal_port_powers(dev.sampler, bun.fields, bun.rho, bun.k0, dev.lines)
    S = scattering_matrix(A, Bmat)
    rows = []
    ok = True
    for j in range(6):
        r = port_roles(j)
        P_inc = float(abs(A[j, j]) ** 2)
        T = {i: float(abs(S[i, j]) ** 2) for i in range(6)}
        sum_T = float(sum(T.values()))
        P_out = float(np.sum(np.abs(Bmat[:, j]) ** 2))
        P_mat = material_absorption(dev.points, dev.tris, bun.fields[j], bun.rho, bun.k0)
        # Modal outgoing may miss TE1+ higher → count in PML/other
        P_pml = P_inc - P_out - P_mat
        residual = abs(1.0 - sum_T - P_mat / max(P_inc, 1e-30) - max(P_pml, 0.0) / max(P_inc, 1e-30))
        # Identity: sum_i |S_ij|^2 + P_mat/P_inc + P_other/P_inc = 1
        # with P_other = P_inc - P_out_fund - P_mat, and sum|S|^2 = P_out_fund/P_inc
        # so residual of (sum|S|^2 + P_mat/Pinc + P_other/Pinc) should be ~0 by algebra
        # if P_other uses same P_out. Check algebraic identity + passivity.
        identity = sum_T + P_mat / max(P_inc, 1e-30) + P_pml / max(P_inc, 1e-30)
        pass_j = (
            abs(identity - 1.0) < 1e-9  # algebraic
            and P_mat >= -0.01 * P_inc
            and P_out <= P_inc * 1.05
            and sum_T <= 1.05
        )
        # Engineering residual: unaccounted if we force P_pml>=0
        P_pml_clip = max(P_pml, 0.0)
        eng_identity = sum_T + P_mat / max(P_inc, 1e-30) + P_pml_clip / max(P_inc, 1e-30)
        eng_res = abs(1.0 - eng_identity)
        if P_pml < -0.02 * P_inc:
            pass_j = False
        if eng_res > BALANCE_RESIDUAL_MAX:
            pass_j = False
        ok = ok and pass_j
        rows.append({
            "source": j,
            "P_inc": P_inc,
            "R": T[j],
            "T_des": T[r["desired"]],
            "T_rev": T[r["reverse"]],
            "T_thr": T[r["through"]],
            "T_leak": T[r["leak_a"]] + T[r["leak_b"]],
            "sum_T": sum_T,
            "P_out_fundamental": P_out,
            "P_material": P_mat,
            "P_PML_and_higher_modes": P_pml,
            "identity_raw": float(identity),
            "eng_residual": float(eng_res),
            "te1_frac_driven": details[j][j].get("te1_power_frac", 0.0),
            "PASS": bool(pass_j),
        })
    return {
        "tag": tag,
        "B_T": B,
        "f_Hz": f_hz,
        "per_source": rows,
        "PASS": bool(ok),
        "threshold": BALANCE_RESIDUAL_MAX,
        "note": (
            "P_PML_and_higher_modes = P_inc - Σ|b_fund|² - P_material. "
            "Includes PML absorption and power in TE1+ not captured by fundamental S."
        ),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="M")
    ap.add_argument("--f-hz", type=float, default=F_HZ)
    ap.add_argument("--B", type=float, default=B_PLUS)
    ap.add_argument(
        "--designs",
        default="uniform,opt_single_uniform_M_Bp0.0500,opt_broadband_C_Bp0.0500",
    )
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Loading {args.grade}", flush=True)
    dev = load_device(args.grade)
    summary = {"grade": args.grade, "designs": {}}
    all_ok = True
    for tag in args.designs.split(","):
        tag = tag.strip()
        print(f"=== {tag} ===", flush=True)
        s = load_s(tag)
        t0 = time.perf_counter()
        rec = analyze(dev, s, args.f_hz, args.B, tag)
        rec["wall_s"] = time.perf_counter() - t0
        dump_json(OUT / f"power_balance_{tag}_{args.grade}.json", rec)
        summary["designs"][tag] = {
            "PASS": rec["PASS"],
            "mean_eng_residual": float(np.mean([r["eng_residual"] for r in rec["per_source"]])),
            "mean_sum_T": float(np.mean([r["sum_T"] for r in rec["per_source"]])),
            "mean_P_mat_frac": float(np.mean([r["P_material"] / max(r["P_inc"], 1e-30) for r in rec["per_source"]])),
            "mean_P_pml_frac": float(
                np.mean([r["P_PML_and_higher_modes"] / max(r["P_inc"], 1e-30) for r in rec["per_source"]])
            ),
            "mean_T_des": float(np.mean([r["T_des"] for r in rec["per_source"]])),
        }
        all_ok = all_ok and rec["PASS"]
        print(f"  PASS={rec['PASS']} {summary['designs'][tag]}", flush=True)
    summary["PASS"] = bool(all_ok)
    dump_json(OUT / f"power_balance_summary_{args.grade}.json", summary)
    print("OVERALL", summary["PASS"], flush=True)
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
