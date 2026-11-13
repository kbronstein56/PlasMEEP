#!/usr/bin/env python3
"""Localize Stage-4 dense-audit holes at ~5.08125 and ~4.725 GHz (forward-only).

Uses fixed Stage-4 q_best. Seeds from existing Stage-4 dense 25/12.5/6.25 audits
so only missing nodes are solved. Optional 3.125 MHz core refine when 6.25 MHz
does not resolve the notch width. No optimization. Stage-5 NOT started.
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

from audit_bb_stage4_dense import freq_grid, parallel_scan  # noqa: E402


def analyze_notch(rows: list[dict], label: str) -> dict:
    rs = sorted(rows, key=lambda r: r["f_GHz"])
    worst = min(rs, key=lambda r: r["J_A2"])
    f0 = float(worst["f_GHz"])
    idx = next(i for i, r in enumerate(rs) if abs(r["f_GHz"] - f0) < 1e-9)

    def recovered(r):
        return r["J_A2"] > -3.0 and r["isolation_dB"] > 5.0

    iL = idx
    while iL > 0 and not recovered(rs[iL - 1]):
        iL -= 1
    iR = idx
    while iR < len(rs) - 1 and not recovered(rs[iR + 1]):
        iR += 1

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

    deep = [r for r in rs if r["J_A2"] < -7.0]
    bad = [r for r in rs if r["J_A2"] < -5.0]
    deep_span = (
        [deep[0]["f_GHz"], deep[-1]["f_GHz"]] if deep else [f0, f0]
    )
    bad_span = [bad[0]["f_GHz"], bad[-1]["f_GHz"]] if bad else [f0, f0]

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
        "deep_core_A2_lt_m7_span_GHz": deep_span,
        "deep_core_width_MHz": (deep_span[1] - deep_span[0]) * 1000.0,
        "bad_core_A2_lt_m5_span_GHz": bad_span,
        "bad_core_width_MHz": (bad_span[1] - bad_span[0]) * 1000.0,
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
    w = analysis["width_half_depth_MHz"]
    if w <= 0:
        return True
    if w < 2.0 * current_df_mhz:
        return True
    tab = analysis["table"]
    center = analysis["notch_center_GHz"]
    i = next(i for i, r in enumerate(tab) if abs(r["f_GHz"] - center) < 1e-9)
    if 0 < i < len(tab) - 1:
        dL = tab[i]["J_A2"] - tab[i - 1]["J_A2"]
        dR = tab[i]["J_A2"] - tab[i + 1]["J_A2"]
        if dL < -3 and dR < -3 and w <= current_df_mhz + 1e-6:
            return True
    # Deep core narrower than ~2 samples
    if analysis["deep_core_width_MHz"] <= current_df_mhz + 1e-6:
        return True
    return False


def load_stage3_5p8125_ref() -> dict | None:
    p = OUT / "BB_STAGE3_NOTCH_LOCALIZATION.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    block = (d.get("region_5p8") or {}).get("pass_6p25MHz") or {}
    return {
        "source": str(p),
        "notch_center_GHz": block.get("notch_center_GHz"),
        "notch_A2": block.get("notch_A2"),
        "width_half_depth_MHz": block.get("width_half_depth_MHz"),
        "deep_core_hint": "Stage-3 report: deep core A2<-7 spans 5.8125–5.81875 (~6.25 MHz)",
        "bad_core_hint": "Stage-3 report: bad core A2<-5 spans 5.80625–5.81875 (~12.5 MHz)",
        "raw": {
            k: block.get(k)
            for k in (
                "notch_center_GHz",
                "notch_A2",
                "notch_IL_dB",
                "notch_isolation_dB",
                "width_half_depth_MHz",
                "n_local_minima",
            )
        },
    }


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

    design = json.loads((OUT / "opt_broadband_stage4_A2_C.json").read_text())
    q = np.asarray(design["q_best"], float)
    print(
        f"LOCALIZE Stage-4 q_best J_BB={design.get('J_BB_best')} — no optimization",
        flush=True,
    )

    seeds = [
        OUT / "bb_stage4_dense_25MHz.json",
        OUT / "bb_stage4_dense_12p5MHz_local.json",
        OUT / "bb_stage4_dense_6p25MHz_local.json",
    ]

    # Pass 1: 5.05–5.15 @ 6.25 MHz (covers true worst + secondary lobe)
    freqs_51 = freq_grid(5.050, 5.150, 0.00625)
    print(f"PASS1 5.08-region 6.25 MHz n={len(freqs_51)}: {freqs_51}", flush=True)
    a51 = parallel_scan(
        args.grade,
        q,
        freqs_51,
        tag="bb_stage4_localize_5p08_6p25MHz",
        resume_path=OUT / "bb_stage4_localize_5p08_6p25MHz.json",
        seed_paths=seeds,
        n_workers=args.n_workers,
        omp_threads=args.omp_threads_per_worker,
    )
    an51 = analyze_notch(a51["rows"], "5p08_region_6p25MHz")
    print(
        f"5.08 notch center={an51['notch_center_GHz']} "
        f"A2={an51['notch_A2']:.3f} "
        f"half-width={an51['width_half_depth_MHz']:.2f} MHz "
        f"deep={an51['deep_core_width_MHz']:.2f} MHz "
        f"bad={an51['bad_core_width_MHz']:.2f} MHz "
        f"n_minima={an51['n_local_minima']}",
        flush=True,
    )

    an51_fine = None
    if needs_finer(an51, 6.25):
        c = float(an51["notch_center_GHz"])
        freqs_fine = freq_grid(c - 0.01875, c + 0.01875, 0.003125)
        print(
            f"PASS1b 5.08 refine 3.125 MHz around {c}: n={len(freqs_fine)} {freqs_fine}",
            flush=True,
        )
        a51f = parallel_scan(
            args.grade,
            q,
            freqs_fine,
            tag="bb_stage4_localize_5p08_3p125MHz",
            resume_path=OUT / "bb_stage4_localize_5p08_3p125MHz.json",
            seed_paths=seeds + [OUT / "bb_stage4_localize_5p08_6p25MHz.json"],
            n_workers=args.n_workers,
            omp_threads=args.omp_threads_per_worker,
        )
        an51_fine = analyze_notch(a51f["rows"], "5p08_region_3p125MHz")
        print(
            f"5.08 fine center={an51_fine['notch_center_GHz']} "
            f"A2={an51_fine['notch_A2']:.3f} "
            f"half-width={an51_fine['width_half_depth_MHz']:.2f} MHz "
            f"deep={an51_fine['deep_core_width_MHz']:.2f} MHz",
            flush=True,
        )
    else:
        print("5.08 feature adequately sampled at 6.25 MHz — skip 3.125", flush=True)

    # Pass 2: 4.70–4.78 @ 6.25 MHz (covers 4.725–4.73125 corpse + shoulder to 4.775)
    freqs_47 = freq_grid(4.700, 4.780, 0.00625)
    print(f"PASS2 4.73-region 6.25 MHz n={len(freqs_47)}: {freqs_47}", flush=True)
    a47 = parallel_scan(
        args.grade,
        q,
        freqs_47,
        tag="bb_stage4_localize_4p73_6p25MHz",
        resume_path=OUT / "bb_stage4_localize_4p73_6p25MHz.json",
        seed_paths=seeds
        + [
            OUT / "bb_stage4_localize_5p08_6p25MHz.json",
            OUT / "bb_stage4_localize_5p08_3p125MHz.json",
        ],
        n_workers=args.n_workers,
        omp_threads=args.omp_threads_per_worker,
    )
    an47 = analyze_notch(a47["rows"], "4p73_region_6p25MHz")
    print(
        f"4.73 notch center={an47['notch_center_GHz']} "
        f"A2={an47['notch_A2']:.3f} "
        f"half-width={an47['width_half_depth_MHz']:.2f} MHz "
        f"deep={an47['deep_core_width_MHz']:.2f} MHz "
        f"bad={an47['bad_core_width_MHz']:.2f} MHz "
        f"n_minima={an47['n_local_minima']}",
        flush=True,
    )

    an47_fine = None
    if needs_finer(an47, 6.25):
        c = float(an47["notch_center_GHz"])
        freqs_fine = freq_grid(c - 0.01875, c + 0.01875, 0.003125)
        print(
            f"PASS2b 4.73 refine 3.125 MHz around {c}: n={len(freqs_fine)} {freqs_fine}",
            flush=True,
        )
        a47f = parallel_scan(
            args.grade,
            q,
            freqs_fine,
            tag="bb_stage4_localize_4p73_3p125MHz",
            resume_path=OUT / "bb_stage4_localize_4p73_3p125MHz.json",
            seed_paths=seeds
            + [
                OUT / "bb_stage4_localize_4p73_6p25MHz.json",
                OUT / "bb_stage4_localize_5p08_6p25MHz.json",
            ],
            n_workers=args.n_workers,
            omp_threads=args.omp_threads_per_worker,
        )
        an47_fine = analyze_notch(a47f["rows"], "4p73_region_3p125MHz")
        print(
            f"4.73 fine center={an47_fine['notch_center_GHz']} "
            f"A2={an47_fine['notch_A2']:.3f} "
            f"half-width={an47_fine['width_half_depth_MHz']:.2f} MHz "
            f"deep={an47_fine['deep_core_width_MHz']:.2f} MHz",
            flush=True,
        )
    else:
        print("4.73 feature adequately sampled at 6.25 MHz — skip 3.125", flush=True)

    # Best available analysis per region
    best_51 = an51_fine or an51
    best_47 = an47_fine or an47
    ref_581 = load_stage3_5p8125_ref()

    widths = {
        "stage3_5p8125_MHz": {
            "deep_core_hint_MHz": 6.25,
            "bad_core_hint_MHz": 12.5,
            "localization": ref_581,
        },
        "stage4_5p08_MHz": {
            "center_GHz": best_51["notch_center_GHz"],
            "half_depth_MHz": best_51["width_half_depth_MHz"],
            "deep_core_MHz": best_51["deep_core_width_MHz"],
            "bad_core_MHz": best_51["bad_core_width_MHz"],
            "A2": best_51["notch_A2"],
        },
        "stage4_4p73_MHz": {
            "center_GHz": best_47["notch_center_GHz"],
            "half_depth_MHz": best_47["width_half_depth_MHz"],
            "deep_core_MHz": best_47["deep_core_width_MHz"],
            "bad_core_MHz": best_47["bad_core_width_MHz"],
            "A2": best_47["notch_A2"],
        },
    }
    observed = [
        widths["stage4_5p08_MHz"]["deep_core_MHz"],
        widths["stage4_4p73_MHz"]["deep_core_MHz"],
        6.25,  # Stage-3 5.8125 deep core
    ]
    # Treat single-sample deep core (0 width) as one spacing of the resolving grid
    resolved = []
    for w, spacing in (
        (widths["stage4_5p08_MHz"]["deep_core_MHz"], 3.125 if an51_fine else 6.25),
        (widths["stage4_4p73_MHz"]["deep_core_MHz"], 3.125 if an47_fine else 6.25),
        (6.25, 6.25),
    ):
        resolved.append(spacing if w <= 1e-9 else w)
    narrowest = min(resolved)

    report = {
        "design": "opt_broadband_stage4_A2_C.json",
        "J_BB_best": design.get("J_BB_best"),
        "note": "Localization only; Stage-5 NOT STARTED / AWAITING AUTHORIZATION",
        "region_5p08": {
            "pass_6p25MHz": an51,
            "pass_3p125MHz": an51_fine,
        },
        "region_4p73": {
            "pass_6p25MHz": an47,
            "pass_3p125MHz": an47_fine,
        },
        "width_comparison": widths,
        "narrowest_reliably_observed_bad_feature_width_MHz": narrowest,
        "sampling_implication": (
            f"Narrowest deep-core width ≈ {narrowest:.3f} MHz. "
            "A fixed grid coarser than this can park a corpse between nodes. "
            "12.5 MHz alone is insufficient in sharp regions; 6.25 MHz is the "
            "minimum protective spacing for known sharp notches, with 3.125 MHz "
            "only where 6.25 leaves a single-sample floor."
        ),
    }
    outp = OUT / "BB_STAGE4_NOTCH_LOCALIZATION.json"
    outp.write_text(json.dumps(report, indent=2) + "\n")
    print("wrote", outp, flush=True)
    print(
        f"NARROWEST_DEEP_CORE_MHz={narrowest:.3f}",
        flush=True,
    )
    print("LOCALIZE_COMPLETE", flush=True)
    print("STAGE-5 STATUS: NOT STARTED / AWAITING AUTHORIZATION", flush=True)


if __name__ == "__main__":
    main()
