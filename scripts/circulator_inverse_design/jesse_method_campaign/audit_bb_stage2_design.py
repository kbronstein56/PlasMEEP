#!/usr/bin/env python3
"""Dense spectral audit of completed Stage-2 A2 design.

Loads q_best from opt_broadband_stage2_A2_C.json.
Reports Objective A2 alongside engineering metrics.
Does NOT re-optimize. Does NOT start Stage-3.
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
from spectral_recon import (  # noqa: E402
    GRID_BARS,
    compare_grids,
    ensure_summary_fields,
    summarize,
)

# Stage-2 frozen objective grid
F_OBJ = [
    4.0,
    4.1,
    4.2,
    4.3,
    4.4,
    4.5,
    4.55,
    4.6,
    4.7,
    4.8,
    4.85,
    4.9,
    5.0,
    5.1,
    5.2,
    5.3,
    5.35,
    5.4,
    5.5,
    5.6,
    5.625,
    5.65,
    5.7,
    5.8,
    5.85,
    5.9,
    5.95,
    6.0,
]


def load_q_best():
    d = json.loads((OUT / "opt_broadband_stage2_A2_C.json").read_text())
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
        mr = bun["mode_resolved"]
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
            "TE1": m["TE1_avg"],
            "R": m["R_avg"],
            "des0": mr["des0_avg"],
            "des1": mr["des1_avg"],
            "des_tot": mr["des_tot_avg"],
            "rev0": mr["rev0_avg"],
            "rev1": mr["rev1_avg"],
            "rev_tot": mr["rev_tot_avg"],
            "iso_tot_dB": mr["iso_tot_avg_dB"],
            "timings": bun["timings"],
        }
        rows.append(row)
        print(
            f"{tag} f={f_ghz:.3f} A2={row['J_A2']:.4f} iso={row['isolation_dB']:.2f} "
            f"IL={row['IL_dB']:.2f} TE1={row['TE1']:.4f} Tdes={row['T_des']:.4f}",
            flush=True,
        )
        out = {
            "tag": tag,
            "grade": grade,
            "GRID_BARS": GRID_BARS,
            "rows": sorted(rows, key=lambda r: r["f_GHz"]),
        }
        (OUT / f"{tag}.json").write_text(json.dumps(out, indent=2) + "\n")
    # summarize expects J_A for worst_J_A — also record worst_J_A2 manually
    summary = summarize(tag, grade, sorted(rows, key=lambda r: r["f_GHz"]))
    a2 = [r["J_A2"] for r in rows if r.get("J_A2") is not None]
    if a2:
        i = int(np.argmin(a2))
        summary["worst_J_A2"] = {"value": float(a2[i]), "f_GHz": rows[i]["f_GHz"]}
        summary["mean_J_A2"] = float(np.mean(a2))
    (OUT / f"{tag}.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def between_grid_failures(audit_rows, obj_freqs):
    obj_set = {round(f, 6) for f in obj_freqs}

    def worst(rows, key, sense="min"):
        vals = [r[key] for r in rows]
        i = int(np.argmin(vals) if sense == "min" else np.argmax(vals))
        return rows[i]

    on_obj = [r for r in audit_rows if round(r["f_GHz"], 6) in obj_set]
    off_obj = [r for r in audit_rows if round(r["f_GHz"], 6) not in obj_set]

    w_iso_all = worst(audit_rows, "isolation_dB", "min")
    w_iso_obj = worst(on_obj, "isolation_dB", "min")
    w_il_all = worst(audit_rows, "IL_dB", "max")
    w_il_obj = worst(on_obj, "IL_dB", "max")
    key_j = "J_A2" if "J_A2" in audit_rows[0] else "J_A"
    w_j_all = worst(audit_rows, key_j, "min")
    w_j_obj = worst(on_obj, key_j, "min")

    iso_floor = w_iso_all["isolation_dB"] + 0.5
    il_ceil = w_il_all["IL_dB"] - 0.2
    j_floor = w_j_all[key_j] + 0.05

    add = []
    for r in off_obj:
        reasons = []
        if r["isolation_dB"] <= iso_floor:
            reasons.append("near_worst_iso")
        if r["IL_dB"] >= il_ceil:
            reasons.append("near_worst_IL")
        if r[key_j] <= j_floor:
            reasons.append("near_worst_J")
        if reasons:
            add.append({"f_GHz": r["f_GHz"], "reasons": reasons, key_j: r[key_j],
                        "IL_dB": r["IL_dB"], "isolation_dB": r["isolation_dB"]})

    must = []
    if round(w_j_all["f_GHz"], 6) not in obj_set and w_j_all[key_j] < w_j_obj[key_j] - 1e-9:
        must.append({"f_GHz": w_j_all["f_GHz"], "why": "off_grid_worst_J", key_j: w_j_all[key_j]})
    if round(w_iso_all["f_GHz"], 6) not in obj_set and w_iso_all["isolation_dB"] < w_iso_obj["isolation_dB"] - 1e-9:
        must.append({"f_GHz": w_iso_all["f_GHz"], "why": "off_grid_worst_iso", "isolation_dB": w_iso_all["isolation_dB"]})
    if round(w_il_all["f_GHz"], 6) not in obj_set and w_il_all["IL_dB"] > w_il_obj["IL_dB"] + 1e-9:
        must.append({"f_GHz": w_il_all["f_GHz"], "why": "off_grid_worst_IL", "IL_dB": w_il_all["IL_dB"]})

    return {
        "J_key": key_j,
        "worst_iso_all": {"f_GHz": w_iso_all["f_GHz"], "value": w_iso_all["isolation_dB"]},
        "worst_iso_obj": {"f_GHz": w_iso_obj["f_GHz"], "value": w_iso_obj["isolation_dB"]},
        "worst_IL_all": {"f_GHz": w_il_all["f_GHz"], "value": w_il_all["IL_dB"]},
        "worst_IL_obj": {"f_GHz": w_il_obj["f_GHz"], "value": w_il_obj["IL_dB"]},
        "worst_J_all": {"f_GHz": w_j_all["f_GHz"], "value": w_j_all[key_j]},
        "worst_J_obj": {"f_GHz": w_j_obj["f_GHz"], "value": w_j_obj[key_j]},
        "pass_iso_bar": abs(w_iso_all["isolation_dB"] - w_iso_obj["isolation_dB"]) < GRID_BARS["worst_isolation_dB_change"]
        and abs(w_iso_all["f_GHz"] - w_iso_obj["f_GHz"]) <= 0.1,
        "pass_IL_bar": abs(w_il_all["IL_dB"] - w_il_obj["IL_dB"]) < GRID_BARS["worst_IL_dB_change"]
        and abs(w_il_all["f_GHz"] - w_il_obj["f_GHz"]) <= 0.1,
        "must_add": must,
        "candidate_add": sorted(add, key=lambda x: x.get(key_j, 0)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--stage", choices=("50", "25", "compare", "all"), default="all")
    ap.add_argument("--f-lo", type=float, default=5.5)
    ap.add_argument("--f-hi", type=float, default=5.8)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    q, opt = load_q_best()
    orbits, _, _ = orbits_and_masks()
    print(
        f"Stage-2 audit design J_BB_best={opt.get('J_BB_best')} "
        f"status={opt.get('status')}",
        flush=True,
    )

    if args.stage in ("50", "all"):
        freqs50 = [4.0 + 0.05 * i for i in range(41)]  # 4.0..6.0
        scan(
            args.grade,
            q,
            orbits,
            freqs50,
            "bb_stage2_audit_50MHz",
            resume_path=OUT / "bb_stage2_audit_50MHz.json",
        )

    if args.stage in ("25", "all"):
        # targeted local midpoints — caller may restrict after 50 MHz
        freqs25 = []
        f = args.f_lo
        while f <= args.f_hi + 1e-12:
            freqs25.append(round(f, 6))
            f += 0.025
        scan(
            args.grade,
            q,
            orbits,
            freqs25,
            "bb_stage2_audit_25MHz_local",
            resume_path=OUT / "bb_stage2_audit_25MHz_local.json",
        )

    if args.stage in ("compare", "all", "50", "25"):
        aud_path = OUT / "bb_stage2_audit_50MHz.json"
        if not aud_path.exists():
            print("no 50 MHz Stage-2 audit yet", flush=True)
            return
        aud = json.loads(aud_path.read_text())
        report = between_grid_failures(aud["rows"], F_OBJ)
        p25 = OUT / "bb_stage2_audit_25MHz_local.json"
        if p25.exists():
            a25 = ensure_summary_fields(json.loads(p25.read_text()), grade=args.grade)
            f_lo = min(r["f_GHz"] for r in a25["rows"])
            f_hi = max(r["f_GHz"] for r in a25["rows"])
            rows50 = [r for r in aud["rows"] if f_lo - 1e-9 <= r["f_GHz"] <= f_hi + 1e-9]
            s50 = summarize("stage2_50_local", args.grade, rows50)
            report["50_vs_25_local"] = compare_grids(s50, a25)
        report["enrichment_needed"] = (
            not report["pass_iso_bar"]
            or not report["pass_IL_bar"]
            or bool(report["must_add"])
        )
        to_add = sorted({round(m["f_GHz"], 6) for m in report["must_add"]})
        for c in report["candidate_add"]:
            f = round(c["f_GHz"], 6)
            if f not in to_add and (
                "near_worst_iso" in c["reasons"] or "near_worst_J" in c["reasons"]
            ):
                to_add.append(f)
        to_add = sorted(set(to_add))
        report["frequencies_to_add_GHz"] = to_add
        report["recommended_new_objective_grid_GHz"] = sorted(set(F_OBJ) | set(to_add))
        (OUT / "bb_stage2_between_grid_report.json").write_text(
            json.dumps(report, indent=2) + "\n"
        )
        print(json.dumps({k: report[k] for k in report if k != "candidate_add"}, indent=2))
        print("candidates:", json.dumps(report["candidate_add"][:15], indent=2))


if __name__ == "__main__":
    main()
