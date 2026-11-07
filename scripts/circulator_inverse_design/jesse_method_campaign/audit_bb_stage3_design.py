#!/usr/bin/env python3
"""Dense spectral audit of completed Stage-3 A2 design + Stage-2 comparison.

Loads q_best from opt_broadband_stage3_A2_C.json.
Reports Objective A2 alongside engineering metrics.
Does NOT re-optimize.
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
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from eval_forward_only import eval_q_forward  # noqa: E402
from multimode_opt_core import B_PLUS, orbits_and_masks  # noqa: E402
from objective_A2 import objective_A2_from_S  # noqa: E402
from spectral_recon import ensure_summary_fields, summarize  # noqa: E402

HOLE_FREQS = [5.75, 5.825]
FOCUS_FREQS = [4.1, 4.3, 5.65, 5.75, 5.825]


def load_q(path: Path):
    d = json.loads(path.read_text())
    return np.asarray(d["q_best"], float), d


def scan(grade, q, orbits, freqs, tag, resume_path=None):
    rows = []
    done = {}
    if resume_path and Path(resume_path).exists():
        prev = json.loads(Path(resume_path).read_text())
        for r in prev.get("rows", []):
            done[round(r["f_GHz"], 6)] = r
            rows.append(r)
    for f_ghz in freqs:
        key = round(float(f_ghz), 6)
        if key in done:
            print(f"skip {f_ghz}", flush=True)
            continue
        bun, _ = eval_q_forward(grade, q, orbits, f_hz=f_ghz * 1e9, b_tesla=B_PLUS)
        m = bun["metrics"]
        J_A2 = float(objective_A2_from_S(bun["S"]))
        row = {
            "f_GHz": float(f_ghz),
            "J_A2": J_A2,
            "J_A": bun["J_A"],
            "J_B": bun["J_B"],
            "T_des": m["T_des_avg"],
            "T_rev": m["T_rev_avg"],
            "IL_dB": m["IL_avg_dB"],
            "isolation_dB": m["isolation_avg_dB"],
            "RL_dB": m["RL_avg_dB"],
            "R": m["R_avg"],
            "TE1": m["TE1_avg"],
        }
        rows.append(row)
        print(
            f"{tag} f={f_ghz:.3f} A2={J_A2:.4f} IL={row['IL_dB']:.2f} "
            f"iso={row['isolation_dB']:.2f} Tdes={row['T_des']:.4e}",
            flush=True,
        )
        if resume_path:
            Path(resume_path).write_text(
                json.dumps({"tag": tag, "grade": grade, "rows": rows}, indent=2) + "\n"
            )
    rows = sorted(rows, key=lambda r: r["f_GHz"])
    summary = summarize(tag, grade, rows)
    summary = ensure_summary_fields(summary, grade=grade)
    out = {"tag": tag, "grade": grade, "rows": rows, **summary}
    return out


def row_at(rows, f_ghz, tol=1e-6):
    for r in rows:
        if abs(r["f_GHz"] - f_ghz) <= tol:
            return r
    # nearest
    return min(rows, key=lambda r: abs(r["f_GHz"] - f_ghz))


def hole_report(label, rows):
    out = {}
    for f in FOCUS_FREQS:
        r = row_at(rows, f)
        out[f"{f:g}"] = {
            "f_GHz": r["f_GHz"],
            "J_A2": r["J_A2"],
            "IL_dB": r["IL_dB"],
            "isolation_dB": r["isolation_dB"],
            "T_des": r["T_des"],
        }
    worst = min(rows, key=lambda r: r["J_A2"])
    out["worst_A2"] = {
        "f_GHz": worst["f_GHz"],
        "J_A2": worst["J_A2"],
        "IL_dB": worst["IL_dB"],
        "isolation_dB": worst["isolation_dB"],
    }
    out["label"] = label
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--stage", choices=("50", "25", "compare", "all"), default="all")
    args = ap.parse_args()

    q3, d3 = load_q(OUT / "opt_broadband_stage3_A2_C.json")
    orbits, _, _ = orbits_and_masks()
    print(
        f"Stage-3 audit q_best J_BB={d3.get('J_BB_best')} "
        f"accepted_nit={d3.get('accepted_nit_total')}",
        flush=True,
    )

    report = {
        "stage3_opt": {
            "J_BB_best": d3.get("J_BB_best"),
            "accepted_nit": d3.get("accepted_nit_total") or d3.get("accepted_nit"),
            "elapsed_s": d3.get("elapsed_s"),
            "status": d3.get("status"),
        }
    }

    a50 = None
    if args.stage in ("50", "all"):
        freqs50 = [4.0 + 0.05 * i for i in range(41)]
        a50 = scan(
            args.grade,
            q3,
            orbits,
            freqs50,
            "bb_stage3_audit_50MHz",
            resume_path=OUT / "bb_stage3_audit_50MHz.json",
        )
        (OUT / "bb_stage3_audit_50MHz.json").write_text(json.dumps(a50, indent=2) + "\n")
        report["audit_50MHz"] = hole_report("stage3_50MHz", a50["rows"])

    if args.stage in ("25", "all"):
        freqs25 = []
        f = 5.50
        while f <= 6.00 + 1e-12:
            freqs25.append(round(f, 6))
            f += 0.025
        a25 = scan(
            args.grade,
            q3,
            orbits,
            freqs25,
            "bb_stage3_audit_25MHz_local",
            resume_path=OUT / "bb_stage3_audit_25MHz_local.json",
        )
        (OUT / "bb_stage3_audit_25MHz_local.json").write_text(
            json.dumps(a25, indent=2) + "\n"
        )
        report["audit_25MHz_local"] = hole_report("stage3_25MHz", a25["rows"])

    # Compare vs Stage-2 audit if present
    s2_path = OUT / "bb_stage2_audit_50MHz.json"
    if a50 is not None and s2_path.exists():
        s2 = json.loads(s2_path.read_text())
        s2_rows = s2.get("rows", [])
        cmp = {
            "stage2": hole_report("stage2_50MHz", s2_rows),
            "stage3": hole_report("stage3_50MHz", a50["rows"]),
            "holes": {},
        }
        for f in HOLE_FREQS:
            r2 = row_at(s2_rows, f)
            r3 = row_at(a50["rows"], f)
            repaired = (r3["J_A2"] > r2["J_A2"] + 1.0) and (
                r3["isolation_dB"] > r2["isolation_dB"] + 5.0
            )
            improved = r3["J_A2"] > r2["J_A2"]
            cmp["holes"][f"{f:g}"] = {
                "stage2": {
                    "J_A2": r2["J_A2"],
                    "IL_dB": r2["IL_dB"],
                    "isolation_dB": r2["isolation_dB"],
                },
                "stage3": {
                    "J_A2": r3["J_A2"],
                    "IL_dB": r3["IL_dB"],
                    "isolation_dB": r3["isolation_dB"],
                },
                "delta_A2": r3["J_A2"] - r2["J_A2"],
                "delta_iso_dB": r3["isolation_dB"] - r2["isolation_dB"],
                "improved": bool(improved),
                "repaired_strong": bool(repaired),
            }
        # 4.1 preservation
        r2 = row_at(s2_rows, 4.1)
        r3 = row_at(a50["rows"], 4.1)
        cmp["at_4p1"] = {
            "stage2": {"T_des": r2["T_des"], "IL_dB": r2["IL_dB"], "isolation_dB": r2["isolation_dB"], "J_A2": r2["J_A2"]},
            "stage3": {"T_des": r3["T_des"], "IL_dB": r3["IL_dB"], "isolation_dB": r3["isolation_dB"], "J_A2": r3["J_A2"]},
            "preserved": bool(r3["T_des"] >= 0.5 * r2["T_des"] and r3["IL_dB"] <= r2["IL_dB"] + 3.0),
        }
        report["stage2_vs_stage3"] = cmp

    (OUT / "BB_STAGE3_AUDIT_REPORT.json").write_text(json.dumps(report, indent=2) + "\n")
    print("wrote", OUT / "BB_STAGE3_AUDIT_REPORT.json", flush=True)
    print(json.dumps(report.get("stage2_vs_stage3", report.get("audit_50MHz")), indent=2), flush=True)


if __name__ == "__main__":
    main()
