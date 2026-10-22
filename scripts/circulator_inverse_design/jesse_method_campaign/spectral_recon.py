#!/usr/bin/env python3
"""Staged 4–6 GHz spectral reconnaissance of a fixed design (forward-only).

Grid-sufficiency bars are FROZEN before interpreting results.
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

# Frozen BEFORE seeing results (Jesse campaign Phase 7)
GRID_BARS = {
    "worst_isolation_dB_change": 0.5,
    "worst_IL_dB_change": 0.2,
    "worst_RL_dB_change": 0.2,
    "note": (
        "Refining the audit grid is considered converged for sampling purposes "
        "when worst isolation changes by <0.5 dB and worst IL/RL by <0.2 dB, "
        "and the frequency location of the worst isolation does not jump by "
        ">100 MHz. Justification: experimental VNA uncertainty and continuum "
        "SF isolation uncertainty (~0.13 dB) make sub-0.5 dB grid jitter "
        "sub-dominant to continuum/mesh error."
    ),
}


def load_q(design: str, orbits):
    from design_io import load_design_q

    return load_design_q(design)


def scan(grade, q, orbits, freqs_ghz, tag, resume_path=None):
    rows = []
    done = {}
    if resume_path and Path(resume_path).exists():
        prev = json.loads(Path(resume_path).read_text())
        for r in prev.get("rows", []):
            done[round(r["f_GHz"], 6)] = r
            rows.append(r)
    for f_ghz in freqs_ghz:
        key = round(float(f_ghz), 6)
        if key in done:
            print(f"skip {f_ghz}", flush=True)
            continue
        bun, s = eval_q_forward(grade, q, orbits, f_hz=f_ghz * 1e9, b_tesla=B_PLUS)
        m = bun["metrics"]
        mr = bun["mode_resolved"]
        row = {
            "f_GHz": float(f_ghz),
            "J_A": bun["J_A"],
            "J_B": bun["J_B"],
            "T_des": m["T_des_avg"],
            "IL_dB": m["IL_avg_dB"],
            "isolation_dB": m["isolation_avg_dB"],
            "RL_dB": m["RL_avg_dB"],
            "TE1": m["TE1_avg"],
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
        done[key] = row
        print(
            f"{tag} f={f_ghz:.3f} JA={row['J_A']:.4f} iso={row['isolation_dB']:.2f} "
            f"IL={row['IL_dB']:.2f} TE1={row['TE1']:.4f}",
            flush=True,
        )
        # checkpoint after each freq
        out = {
            "tag": tag,
            "grade": grade,
            "GRID_BARS": GRID_BARS,
            "rows": sorted(rows, key=lambda r: r["f_GHz"]),
        }
        (OUT / f"spectral_{tag}.json").write_text(json.dumps(out, indent=2) + "\n")
    rows = sorted(rows, key=lambda r: r["f_GHz"])
    return summarize(tag, grade, rows)


def summarize(tag, grade, rows):
    def worst(key, sense="min"):
        vals = [r[key] for r in rows]
        if sense == "min":
            i = int(np.argmin(vals))
        else:
            i = int(np.argmax(vals))
        return {"value": vals[i], "f_GHz": rows[i]["f_GHz"]}

    summary = {
        "tag": tag,
        "grade": grade,
        "n_freq": len(rows),
        "spacing_GHz": None
        if len(rows) < 2
        else float(np.median(np.diff([r["f_GHz"] for r in rows]))),
        "worst_isolation_dB": worst("isolation_dB", "min"),
        "worst_IL_dB": worst("IL_dB", "max"),
        "worst_RL_dB": worst("RL_dB", "min"),
        "worst_J_A": worst("J_A", "min"),
        "worst_des_tot": worst("des_tot", "min"),
        "max_rev_tot": worst("rev_tot", "max"),
        "rows": rows,
        "GRID_BARS": GRID_BARS,
    }
    (OUT / f"spectral_{tag}.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def compare_grids(coarse, fine):
    """Apply frozen bars: does fine change worst metrics vs coarse?"""
    def delta(a, b):
        return abs(a["value"] - b["value"])

    d_iso = delta(coarse["worst_isolation_dB"], fine["worst_isolation_dB"])
    d_il = delta(coarse["worst_IL_dB"], fine["worst_IL_dB"])
    d_rl = delta(coarse["worst_RL_dB"], fine["worst_RL_dB"])
    df = abs(coarse["worst_isolation_dB"]["f_GHz"] - fine["worst_isolation_dB"]["f_GHz"])
    return {
        "d_worst_iso_dB": d_iso,
        "d_worst_IL_dB": d_il,
        "d_worst_RL_dB": d_rl,
        "d_worst_iso_location_GHz": df,
        "pass_iso": d_iso < GRID_BARS["worst_isolation_dB_change"],
        "pass_IL": d_il < GRID_BARS["worst_IL_dB_change"],
        "pass_RL": d_rl < GRID_BARS["worst_RL_dB_change"],
        "pass_location": df <= 0.1,
        "pass_all": (
            d_iso < GRID_BARS["worst_isolation_dB_change"]
            and d_il < GRID_BARS["worst_IL_dB_change"]
            and d_rl < GRID_BARS["worst_RL_dB_change"]
            and df <= 0.1
        ),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--design", default="opt_uniform_C_warm_M")
    ap.add_argument(
        "--stage",
        choices=["200", "100", "50", "compare", "all"],
        default="200",
    )
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "GRID_SUFFICIENCY_BARS_FROZEN.json").write_text(json.dumps(GRID_BARS, indent=2) + "\n")

    orbits, _, _ = orbits_and_masks()
    q = load_q(args.design, orbits)

    if args.stage in ("200", "all"):
        freqs = [4.0 + 0.2 * i for i in range(11)]
        scan(args.grade, q, orbits, freqs, "200MHz")
    if args.stage in ("100", "all"):
        freqs = [4.0 + 0.1 * i for i in range(21)]
        scan(
            args.grade,
            q,
            orbits,
            freqs,
            "100MHz",
            resume_path=OUT / "spectral_100MHz.json",
        )
    if args.stage in ("50", "all"):
        # selective: start with full 50 MHz (41 pts) — expensive; allow resume
        freqs = [4.0 + 0.05 * i for i in range(41)]
        scan(
            args.grade,
            q,
            orbits,
            freqs,
            "50MHz",
            resume_path=OUT / "spectral_50MHz.json",
        )
    if args.stage in ("compare", "all"):
        paths = {
            "200": OUT / "spectral_200MHz.json",
            "100": OUT / "spectral_100MHz.json",
            "50": OUT / "spectral_50MHz.json",
        }
        loaded = {k: json.loads(p.read_text()) for k, p in paths.items() if p.exists()}
        comps = {}
        if "200" in loaded and "100" in loaded:
            comps["200_vs_100"] = compare_grids(loaded["200"], loaded["100"])
        if "100" in loaded and "50" in loaded:
            comps["100_vs_50"] = compare_grids(loaded["100"], loaded["50"])
        (OUT / "spectral_grid_comparison.json").write_text(json.dumps(comps, indent=2) + "\n")
        print(json.dumps(comps, indent=2))


if __name__ == "__main__":
    main()
