#!/usr/bin/env python3
"""Phases 5–7: FULL12 Onsager + power balance + TE0 metrics for existing designs."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "multimode_ports"
DES = ROOT / "outputs" / "circulator_inverse_design"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from full91_dtn_solve import full12_s, onsager_residual  # noqa: E402
from validate_modal_s import material_absorption  # noqa: E402

ONSAGER_BAR = 5e-3
BALANCE_BAR = 0.05


def load_s(tag: str) -> np.ndarray:
    if tag == "uniform":
        return np.ones(91)
    d = json.loads((DES / f"{tag}.json").read_text())
    return np.asarray(d["best"]["s"], float)


def port_roles(j: int) -> dict:
    return {
        "desired": (j + 1) % 6,
        "reverse": (j - 1) % 6,
        "through": (j + 3) % 6,
        "leak_a": (j + 2) % 6,
        "leak_b": (j + 4) % 6,
    }


def te0_metrics(S: np.ndarray) -> dict:
    """S is 12×12; TE0 channels are even indices 2p."""
    rows = []
    for j in range(6):
        col = 2 * j  # TE0 drive
        r = port_roles(j)
        T = {i: float(abs(S[2 * i, col]) ** 2) for i in range(6)}
        te1_all = float(sum(abs(S[2 * i + 1, col]) ** 2 for i in range(6)))
        T_des = T[r["desired"]]
        rows.append({
            "source": j,
            "T_des_TE0": T_des,
            "T_rev_TE0": T[r["reverse"]],
            "T_thr_TE0": T[r["through"]],
            "T_leak_TE0": T[r["leak_a"]] + T[r["leak_b"]],
            "R_TE0": T[j],
            "TE1_total": te1_all,
            "IL_TE0_dB": float(-10 * np.log10(max(T_des, 1e-30))),
            "isolation_TE0_dB": float(10 * np.log10(max(T_des, 1e-30) / max(T[r["reverse"]], 1e-30))),
            "RL_TE0_dB": float(-10 * np.log10(max(T[j], 1e-30))),
            "undesired_modal": float(
                T[r["reverse"]] + T[r["through"]] + T[r["leak_a"]] + T[r["leak_b"]] + T[j] + te1_all
            ),
        })
    avg = {k: float(np.mean([row[k] for row in rows])) for k in rows[0] if k != "source"}
    return {"per_source": rows, "mean": avg}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument(
        "--designs",
        default="uniform,opt_single_uniform_M_Bp0.0500,opt_broadband_C_Bp0.0500",
    )
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    report = {"grade": args.grade, "onsager_bar": ONSAGER_BAR, "designs": {}}
    t0 = time.time()
    for tag in args.designs.split(","):
        tag = tag.strip()
        s = load_s(tag)
        print(f"=== {tag} ===", flush=True)
        Sp, meta = full12_s(args.grade, s, b_tesla=0.05)
        Sm, _ = full12_s(args.grade, s, b_tesla=-0.05)
        S0, _ = full12_s(args.grade, s, b_tesla=0.0)
        np.savez_compressed(OUT / f"S12_{tag}_Bp.npz", S=Sp)
        r_F = onsager_residual(Sp, Sm)
        r_0 = onsager_residual(S0, S0)
        mets = te0_metrics(Sp)
        # Power: for each TE0 column, sum |S|^2 over 12 channels (unit incident via S def)
        bal_rows = []
        for j in range(6):
            col = 2 * j
            P_out = float(np.sum(np.abs(Sp[:, col]) ** 2))
            P_te0 = float(sum(abs(Sp[2 * i, col]) ** 2 for i in range(6)))
            P_te1 = float(sum(abs(Sp[2 * i + 1, col]) ** 2 for i in range(6)))
            # Material absorption requires fields — skip full field dump here;
            # remaining = 1 - P_out is non-port / absorption lump for S-normalized.
            rem = 1.0 - P_out
            bal_rows.append({
                "source": j,
                "P_out_prop": P_out,
                "P_TE0": P_te0,
                "P_TE1": P_te1,
                "P_remaining_abs_or_nonport": rem,
            })
        mean_rem = float(np.mean([r["P_remaining_abs_or_nonport"] for r in bal_rows]))
        design_rec = {
            "r_F_plus": r_F,
            "r_F_zero": r_0,
            "FULL12_ONSAGER": "PASS" if r_F <= ONSAGER_BAR else "FAIL",
            "te0_metrics": mets,
            "power_S_normalized": {
                "per_TE0_source": bal_rows,
                "mean_remaining": mean_rem,
                "note": (
                    "With power-normalized S, incident |a|=1 so Σ|b|² + remaining = 1. "
                    "remaining lumps material absorption + non-port radiation. "
                    "TE1 is reported separately, not as absorption."
                ),
            },
            "meta": meta,
        }
        report["designs"][tag] = design_rec
        (OUT / f"eval_{tag}_{args.grade}.json").write_text(json.dumps(design_rec, indent=2) + "\n")
        print(
            f"  r_F={r_F:.4e} T_des={mets['mean']['T_des_TE0']:.4f} "
            f"IL={mets['mean']['IL_TE0_dB']:.2f}dB iso={mets['mean']['isolation_TE0_dB']:.2f}dB "
            f"TE1={mets['mean']['TE1_total']:.4f} rem={mean_rem:.4f}",
            flush=True,
        )

    report["elapsed_s"] = time.time() - t0
    (OUT / "EXISTING_DESIGNS_FULL12_REPORT.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"elapsed_s": report["elapsed_s"], "tags": list(report["designs"])}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
