#!/usr/bin/env python3
"""Short localization audit of Stage-3 off-grid notches (forward-only).

Resolves the 5.8125 GHz and 4.775 GHz features with 6.25 MHz (and optional
3.125 MHz) sampling. No optimization.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "jesse_method_campaign"),
]

from audit_bb_stage3_dense import freq_grid, parallel_scan  # noqa: E402


def analyze_notch(rows: list[dict], label: str) -> dict:
    rs = sorted(rows, key=lambda r: r["f_GHz"])
    worst = min(rs, key=lambda r: r["J_A2"])
    # Contiguous region where A2 is within 3 of the worst (notch core)
    thr = float(worst["J_A2"]) + 3.0
    in_core = [r for r in rs if r["J_A2"] <= thr]
    # Also isolation-based core: iso < 0 and near worst
    iso_bad = [r for r in rs if r["isolation_dB"] < 0.0]
    # Half-depth width: A2 below midpoint between worst and best neighbor floor
    best = max(rs, key=lambda r: r["J_A2"])
    half = 0.5 * (float(worst["J_A2"]) + float(best["J_A2"]))
    # More stable: points where A2 < -4 or iso < 0 around the worst
    f0 = float(worst["f_GHz"])
    # Walk left/right from worst until A2 recovers above -3 or iso>5
    idx = next(i for i, r in enumerate(rs) if abs(r["f_GHz"] - f0) < 1e-9)

    def recovered(r):
        return r["J_A2"] > -3.0 and r["isolation_dB"] > 5.0

    iL = idx
    while iL > 0 and not recovered(rs[iL - 1]):
        iL -= 1
    iR = idx
    while iR < len(rs) - 1 and not recovered(rs[iR + 1]):
        iR += 1
    # Narrower: A2 within 50% of depth from local baseline
    # Use adjacent local maxima as baseline
    left_max = max(rs[: idx + 1], key=lambda r: r["J_A2"])
    right_max = max(rs[idx:], key=lambda r: r["J_A2"])
    baseline = 0.5 * (left_max["J_A2"] + right_max["J_A2"])
    depth = baseline - worst["J_A2"]
    half_level = worst["J_A2"] + 0.5 * depth
    half_pts = [r for r in rs if r["J_A2"] <= half_level]
    if len(half_pts) >= 2:
        f_half_lo = half_pts[0]["f_GHz"]
        f_half_hi = half_pts[-1]["f_GHz"]
        width_half = f_half_hi - f_half_lo
    else:
        f_half_lo = f_half_hi = worst["f_GHz"]
        width_half = 0.0

    # Count local A2 minima (interior)
    minima = []
    for i in range(1, len(rs) - 1):
        if rs[i]["J_A2"] < rs[i - 1]["J_A2"] and rs[i]["J_A2"] < rs[i + 1]["J_A2"]:
            minima.append(rs[i])

    table = []
    for r in rs:
        wc = r.get("worst_competing") or {}
        table.append(
            {
                "f_GHz": r["f_GHz"],
                "J_A2": r["J_A2"],
                "T_des": r["T_des"],
                "IL_dB": r["IL_dB"],
                "isolation_dB": r["isolation_dB"],
                "worst_competing_port": wc.get("competing_port"),
                "worst_source_port": wc.get("source_port"),
                "competing_T_rev": wc.get("T_rev"),
            }
        )

    return {
        "label": label,
        "n_freq": len(rs),
        "spacing_GHz": round(rs[1]["f_GHz"] - rs[0]["f_GHz"], 6) if len(rs) > 1 else None,
        "notch_center_GHz": worst["f_GHz"],
        "notch_A2": worst["J_A2"],
        "notch_IL_dB": worst["IL_dB"],
        "notch_isolation_dB": worst["isolation_dB"],
        "notch_T_des": worst["T_des"],
        "width_half_depth_GHz": width_half,
        "width_half_depth_MHz": width_half * 1000.0,
        "half_depth_span_GHz": [f_half_lo, f_half_hi],
        "recovery_span_GHz": [rs[iL]["f_GHz"], rs[iR]["f_GHz"]],
        "recovery_width_MHz": (rs[iR]["f_GHz"] - rs[iL]["f_GHz"]) * 1000.0,
        "left_shoulder": {
            "f_GHz": left_max["f_GHz"],
            "J_A2": left_max["J_A2"],
            "isolation_dB": left_max["isolation_dB"],
        },
        "right_shoulder": {
            "f_GHz": right_max["f_GHz"],
            "J_A2": right_max["J_A2"],
            "isolation_dB": right_max["isolation_dB"],
        },
        "n_local_minima": len(minima),
        "local_minima": [
            {
                "f_GHz": m["f_GHz"],
                "J_A2": m["J_A2"],
                "isolation_dB": m["isolation_dB"],
                "IL_dB": m["IL_dB"],
            }
            for m in minima
        ],
        "table": table,
    }


def needs_finer(analysis: dict, current_df_mhz: float) -> bool:
    """Refine if half-depth width is < ~2 sample spacings or single-point floor."""
    w = analysis["width_half_depth_MHz"]
    if w <= 0:
        return True
    if w < 2.0 * current_df_mhz:
        return True
    # If only one point is near the floor (much worse than both neighbors)
    tab = analysis["table"]
    center = analysis["notch_center_GHz"]
    i = next(i for i, r in enumerate(tab) if abs(r["f_GHz"] - center) < 1e-9)
    if 0 < i < len(tab) - 1:
        dL = tab[i]["J_A2"] - tab[i - 1]["J_A2"]
        dR = tab[i]["J_A2"] - tab[i + 1]["J_A2"]
        # Very sharp single-sample spike
        if dL < -5 and dR < -5 and w <= current_df_mhz + 1e-6:
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--omp-threads-per-worker", type=int, default=1)
    ap.add_argument("--grade", default="C")
    args = ap.parse_args()

    for k, v in (
        ("OMP_NUM_THREADS", "1"),
        ("OPENBLAS_NUM_THREADS", "1"),
        ("MKL_NUM_THREADS", "1"),
    ):
        os.environ.setdefault(k, v)

    design = json.loads((OUT / "opt_broadband_stage3_A2_C.json").read_text())
    q = np.asarray(design["q_best"], float)
    print(
        f"LOCALIZE Stage-3 q_best J_BB={design.get('J_BB_best')} — no optimization",
        flush=True,
    )

    seeds = [
        OUT / "bb_stage3_dense_25MHz.json",
        OUT / "bb_stage3_dense_12p5MHz_local.json",
        OUT / "bb_stage3_localize_5p8_6p25MHz.json",
        OUT / "bb_stage3_localize_4p775_6p25MHz.json",
    ]

    # Pass 1: 5.775–5.850 @ 6.25 MHz
    freqs_58 = freq_grid(5.775, 5.850, 0.00625)
    print(f"PASS1 5.8-region 6.25 MHz n={len(freqs_58)}: {freqs_58}", flush=True)
    a58 = parallel_scan(
        args.grade,
        q,
        freqs_58,
        tag="bb_stage3_localize_5p8_6p25MHz",
        resume_path=OUT / "bb_stage3_localize_5p8_6p25MHz.json",
        seed_paths=seeds,
        n_workers=args.n_workers,
        omp_threads=args.omp_threads_per_worker,
    )
    an58 = analyze_notch(a58["rows"], "5.8_region_6p25MHz")
    print(
        f"5.8 notch center={an58['notch_center_GHz']} "
        f"half-width={an58['width_half_depth_MHz']:.2f} MHz "
        f"n_minima={an58['n_local_minima']}",
        flush=True,
    )

    # Optional 3.125 MHz refine around notch
    an58_fine = None
    if needs_finer(an58, 6.25):
        c = float(an58["notch_center_GHz"])
        # ±18.75 MHz window at 3.125 MHz
        freqs_fine = freq_grid(c - 0.01875, c + 0.01875, 0.003125)
        print(
            f"PASS1b 5.8 refine 3.125 MHz around {c}: n={len(freqs_fine)} {freqs_fine}",
            flush=True,
        )
        a58f = parallel_scan(
            args.grade,
            q,
            freqs_fine,
            tag="bb_stage3_localize_5p8_3p125MHz",
            resume_path=OUT / "bb_stage3_localize_5p8_3p125MHz.json",
            seed_paths=seeds + [OUT / "bb_stage3_localize_5p8_6p25MHz.json"],
            n_workers=args.n_workers,
            omp_threads=args.omp_threads_per_worker,
        )
        an58_fine = analyze_notch(a58f["rows"], "5.8_region_3p125MHz")
        print(
            f"5.8 fine center={an58_fine['notch_center_GHz']} "
            f"half-width={an58_fine['width_half_depth_MHz']:.2f} MHz "
            f"n_minima={an58_fine['n_local_minima']}",
            flush=True,
        )
    else:
        print("5.8 feature adequately sampled at 6.25 MHz — skip 3.125", flush=True)

    # Pass 2: 4.725–4.825 @ 6.25 MHz
    freqs_47 = freq_grid(4.725, 4.825, 0.00625)
    print(f"PASS2 4.775-region 6.25 MHz n={len(freqs_47)}: {freqs_47}", flush=True)
    a47 = parallel_scan(
        args.grade,
        q,
        freqs_47,
        tag="bb_stage3_localize_4p775_6p25MHz",
        resume_path=OUT / "bb_stage3_localize_4p775_6p25MHz.json",
        seed_paths=seeds
        + [
            OUT / "bb_stage3_localize_5p8_6p25MHz.json",
            OUT / "bb_stage3_localize_5p8_3p125MHz.json",
        ],
        n_workers=args.n_workers,
        omp_threads=args.omp_threads_per_worker,
    )
    an47 = analyze_notch(a47["rows"], "4.775_region_6p25MHz")
    print(
        f"4.775 notch center={an47['notch_center_GHz']} "
        f"half-width={an47['width_half_depth_MHz']:.2f} MHz "
        f"n_minima={an47['n_local_minima']}",
        flush=True,
    )

    an47_fine = None
    if needs_finer(an47, 6.25):
        c = float(an47["notch_center_GHz"])
        freqs_fine = freq_grid(c - 0.01875, c + 0.01875, 0.003125)
        print(
            f"PASS2b 4.775 refine 3.125 MHz around {c}: n={len(freqs_fine)} {freqs_fine}",
            flush=True,
        )
        a47f = parallel_scan(
            args.grade,
            q,
            freqs_fine,
            tag="bb_stage3_localize_4p775_3p125MHz",
            resume_path=OUT / "bb_stage3_localize_4p775_3p125MHz.json",
            seed_paths=seeds
            + [
                OUT / "bb_stage3_localize_4p775_6p25MHz.json",
                OUT / "bb_stage3_localize_5p8_6p25MHz.json",
            ],
            n_workers=args.n_workers,
            omp_threads=args.omp_threads_per_worker,
        )
        an47_fine = analyze_notch(a47f["rows"], "4.775_region_3p125MHz")
        print(
            f"4.775 fine center={an47_fine['notch_center_GHz']} "
            f"half-width={an47_fine['width_half_depth_MHz']:.2f} MHz",
            flush=True,
        )
    else:
        print("4.775 feature adequately sampled at 6.25 MHz — skip 3.125", flush=True)

    report = {
        "design": "opt_broadband_stage3_A2_C.json",
        "J_BB_best": design.get("J_BB_best"),
        "note": "Localization only; Stage-4 NOT STARTED",
        "region_5p8": {
            "pass_6p25MHz": an58,
            "pass_3p125MHz": an58_fine,
        },
        "region_4p775": {
            "pass_6p25MHz": an47,
            "pass_3p125MHz": an47_fine,
        },
    }
    outp = OUT / "BB_STAGE3_NOTCH_LOCALIZATION.json"
    outp.write_text(json.dumps(report, indent=2) + "\n")
    print("wrote", outp, flush=True)
    print("LOCALIZE_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
