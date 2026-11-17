#!/usr/bin/env python3
"""Stage-5 q_best local notch-width characterization (forward-only).

Deterministic windows at 1.5625 MHz around Stage-5 dense-audit residuals;
optional 0.78125 MHz core refine when the floor is single-sample or the
deep core remains unresolved. No optimization. Stage-6 NOT STARTED.
"""
from __future__ import annotations

import argparse
import csv
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

from audit_bb_stage5_dense import (  # noqa: E402
    freq_grid,
    load_seed_rows,
    local_minima,
    parallel_scan,
)

# Centers from Stage-5 3.125 MHz dense audit (A2 ≲ -8)
CENTERS = {
    "4p734": 4.734375,
    "4p228": 4.228125,
    "5p084": 5.084375,
    "5p572": 5.571875,
}

HALF_WINDOW_GHZ = 0.040  # ±40 MHz captures shoulders for observed few-MHz cores
DF_FINE = 0.0015625  # 1.5625 MHz
DF_ULTRA = 0.00078125  # 0.78125 MHz
ULTRA_HALF_GHZ = 0.008  # ±8 MHz around unresolved floor


def analyze_notch(rows: list[dict], label: str, spacing_GHz: float) -> dict:
    rs = sorted(rows, key=lambda r: r["f_GHz"])
    worst = min(rs, key=lambda r: r["J_A2"])
    f0 = float(worst["f_GHz"])
    idx = next(i for i, r in enumerate(rs) if abs(r["f_GHz"] - f0) < 1e-12)

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
        f_half_lo, f_half_hi = half_pts[0]["f_GHz"], half_pts[-1]["f_GHz"]
        width_half = f_half_hi - f_half_lo
    else:
        f_half_lo = f_half_hi = f0
        width_half = 0.0

    def contiguous_span(pred):
        """Largest contiguous run containing the worst that satisfies pred."""
        if not pred(worst):
            return [f0, f0], 0.0, 0
        lo = hi = idx
        while lo > 0 and pred(rs[lo - 1]):
            lo -= 1
        while hi < len(rs) - 1 and pred(rs[hi + 1]):
            hi += 1
        span = [rs[lo]["f_GHz"], rs[hi]["f_GHz"]]
        # Inclusive width of sampled span; report both span and n_samples
        width = span[1] - span[0]
        return span, width, hi - lo + 1

    deep_span, deep_w, n_deep = contiguous_span(lambda r: r["J_A2"] < -7.0)
    bad_span, bad_w, n_bad = contiguous_span(lambda r: r["J_A2"] < -5.0)

    mins = local_minima(rs, key="J_A2", pad=1)
    important_mins = [m for m in mins if m["J_A2"] <= -5.0]

    # Recovery rate: ΔA2 / Δf on each side using first recovered neighbor or window edge
    def side_rate(side: str):
        if side == "left":
            pts = list(reversed(rs[: idx + 1]))
        else:
            pts = rs[idx:]
        if len(pts) < 2:
            return None
        # slope from worst to first point with A2 > -3 or last available
        target = pts[0]
        for p in pts[1:]:
            target = p
            if p["J_A2"] > -3.0:
                break
        df = abs(target["f_GHz"] - f0)
        if df < 1e-12:
            return None
        return {
            "df_GHz": df,
            "dA2": float(target["J_A2"] - worst["J_A2"]),
            "dA2_per_MHz": float((target["J_A2"] - worst["J_A2"]) / (df * 1e3)),
            "end_f_GHz": target["f_GHz"],
            "end_A2": target["J_A2"],
        }

    # Single-sample floor: worst is a local min and both neighbors are substantially better
    neighbors = []
    if idx > 0:
        neighbors.append(rs[idx - 1])
    if idx < len(rs) - 1:
        neighbors.append(rs[idx + 1])
    single_sample_floor = False
    if neighbors:
        gap = min(r["J_A2"] - worst["J_A2"] for r in neighbors)
        single_sample_floor = gap > 1.0  # >1 A2 unit jump to nearest sample

    needs_ultra = (
        single_sample_floor
        or n_deep <= 1
        or (n_deep == 2 and deep_w <= spacing_GHz + 1e-12)
    )

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
            }
        )

    return {
        "label": label,
        "spacing_GHz": spacing_GHz,
        "spacing_MHz": spacing_GHz * 1e3,
        "n_freq": len(rs),
        "notch_center_GHz": f0,
        "notch_A2": worst["J_A2"],
        "notch_IL_dB": worst["IL_dB"],
        "notch_isolation_dB": worst["isolation_dB"],
        "notch_T_des": worst["T_des"],
        "worst_competing": worst.get("worst_competing"),
        "width_half_depth_GHz": width_half,
        "width_half_depth_MHz": width_half * 1e3,
        "half_depth_span_GHz": [f_half_lo, f_half_hi],
        "deep_core_A2_lt_m7_span_GHz": deep_span,
        "deep_core_width_MHz": deep_w * 1e3,
        "deep_core_n_samples": n_deep,
        "bad_core_A2_lt_m5_span_GHz": bad_span,
        "bad_core_width_MHz": bad_w * 1e3,
        "bad_core_n_samples": n_bad,
        "recovery_span_GHz": [rs[iL]["f_GHz"], rs[iR]["f_GHz"]],
        "recovery_width_MHz": (rs[iR]["f_GHz"] - rs[iL]["f_GHz"]) * 1e3,
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
        "recovery_rate_left": side_rate("left"),
        "recovery_rate_right": side_rate("right"),
        "n_local_minima": len(mins),
        "n_important_local_minima_A2_le_m5": len(important_mins),
        "local_minima": [
            {
                "f_GHz": m["f_GHz"],
                "J_A2": m["J_A2"],
                "IL_dB": m["IL_dB"],
                "isolation_dB": m["isolation_dB"],
            }
            for m in mins
        ],
        "important_local_minima": [
            {
                "f_GHz": m["f_GHz"],
                "J_A2": m["J_A2"],
                "IL_dB": m["IL_dB"],
                "isolation_dB": m["isolation_dB"],
            }
            for m in important_mins
        ],
        "single_sample_floor": single_sample_floor,
        "needs_0p78125_refine": needs_ultra,
        "feature_smooth_deterministic_note": (
            "Measured FEM samples only; no Lorentzian assumed. "
            "Widths are contiguous-sample spans at the stated A2 thresholds."
        ),
        "table": table,
    }


def window_freqs(center: float, half: float, df: float) -> list[float]:
    """Symmetric grid about *center* so the notch seed lands on a sample."""
    n = int(round(half / df))
    return [round(center + i * df, 6) for i in range(-n, n + 1)]


def write_region_csv(path: Path, rows: list[dict]):
    fields = [
        "f_GHz",
        "J_A2",
        "T_des",
        "IL_dB",
        "isolation_dB",
        "worst_competing_port",
        "worst_source_port",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            wc = r.get("worst_competing") or {}
            w.writerow(
                {
                    "f_GHz": r["f_GHz"],
                    "J_A2": r["J_A2"],
                    "T_des": r["T_des"],
                    "IL_dB": r["IL_dB"],
                    "isolation_dB": r["isolation_dB"],
                    "worst_competing_port": wc.get("competing_port"),
                    "worst_source_port": wc.get("source_port"),
                }
            )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--omp-threads", type=int, default=1)
    ap.add_argument("--skip-ultra", action="store_true")
    args = ap.parse_args()

    for k, v in (
        ("OMP_NUM_THREADS", str(args.omp_threads)),
        ("OPENBLAS_NUM_THREADS", str(args.omp_threads)),
        ("MKL_NUM_THREADS", str(args.omp_threads)),
        ("PYTHONUNBUFFERED", "1"),
    ):
        os.environ.setdefault(k, v)

    design = OUT / "opt_broadband_stage5_A2_C.json"
    d = json.loads(design.read_text())
    q = np.asarray(d["q_best"], float)
    print(
        f"Stage-5 notch-width audit design={design.name} "
        f"J_BB={d.get('J_BB_best'):.6f} n_workers={args.n_workers} "
        f"— Stage-6 NOT STARTED",
        flush=True,
    )

    seed_paths = [
        OUT / "bb_stage5_dense_6p25MHz.json",
        OUT / "bb_stage5_dense_3p125MHz_local.json",
        OUT / "bb_stage5_opt_grid_seed.json",
    ]

    regions = {}
    for tag, center in CENTERS.items():
        freqs = window_freqs(center, HALF_WINDOW_GHZ, DF_FINE)
        resume = OUT / f"bb_stage5_notch_{tag}_1p5625MHz.json"
        print(
            f"\n=== {tag} center={center:.6f} GHz  "
            f"window=[{center-HALF_WINDOW_GHZ:.6f},{center+HALF_WINDOW_GHZ:.6f}] "
            f"df=1.5625 MHz n={len(freqs)} ===",
            flush=True,
        )
        out = parallel_scan(
            args.grade,
            q,
            freqs,
            tag=f"bb_stage5_notch_{tag}_1p5625MHz",
            resume_path=resume,
            seed_paths=seed_paths,
            n_workers=args.n_workers,
            omp_threads=args.omp_threads,
        )
        rows = out["rows"]
        write_region_csv(resume.with_suffix(".csv"), rows)
        analysis = analyze_notch(rows, f"{tag}_1p5625MHz", DF_FINE)
        regions[tag] = {"pass_1p5625MHz": analysis, "pass_0p78125MHz": None}

        if (not args.skip_ultra) and analysis["needs_0p78125_refine"]:
            fmin = analysis["notch_center_GHz"]
            ufreqs = window_freqs(fmin, ULTRA_HALF_GHZ, DF_ULTRA)
            uresume = OUT / f"bb_stage5_notch_{tag}_0p78125MHz.json"
            print(
                f"--- {tag} ULTRA refine around {fmin:.6f} "
                f"df=0.78125 MHz n={len(ufreqs)} ---",
                flush=True,
            )
            uout = parallel_scan(
                args.grade,
                q,
                ufreqs,
                tag=f"bb_stage5_notch_{tag}_0p78125MHz",
                resume_path=uresume,
                seed_paths=seed_paths + [resume],
                n_workers=args.n_workers,
                omp_threads=args.omp_threads,
            )
            urows = uout["rows"]
            write_region_csv(uresume.with_suffix(".csv"), urows)
            uanalysis = analyze_notch(urows, f"{tag}_0p78125MHz", DF_ULTRA)
            regions[tag]["pass_0p78125MHz"] = uanalysis

    # Compact width table (authoritative pass = ultra if present else fine)
    compact_rows = []
    for tag, blob in regions.items():
        fine = blob["pass_1p5625MHz"]
        ultra = blob["pass_0p78125MHz"]
        auth = ultra or fine
        compact_rows.append(
            {
                "region": tag,
                "authoritative_spacing_MHz": auth["spacing_MHz"],
                "min_f_GHz": auth["notch_center_GHz"],
                "min_A2": auth["notch_A2"],
                "IL_dB": auth["notch_IL_dB"],
                "isolation_dB": auth["notch_isolation_dB"],
                "T_des": auth["notch_T_des"],
                "A2_lt_m5_width_MHz": auth["bad_core_width_MHz"],
                "A2_lt_m5_n_samples": auth["bad_core_n_samples"],
                "A2_lt_m7_width_MHz": auth["deep_core_width_MHz"],
                "A2_lt_m7_n_samples": auth["deep_core_n_samples"],
                "half_depth_width_MHz": auth["width_half_depth_MHz"],
                "n_important_minima": auth["n_important_local_minima_A2_le_m5"],
                "single_sample_floor_at_1p5625": fine["single_sample_floor"],
                "ultra_refined": ultra is not None,
            }
        )

    # Historical comparison stubs from prior localization JSONs
    s3 = json.loads((OUT / "BB_STAGE3_NOTCH_LOCALIZATION.json").read_text())
    s4 = json.loads((OUT / "BB_STAGE4_NOTCH_LOCALIZATION.json").read_text())

    def hist_row(stage, name, a2, deep_mhz, bad_mhz, note):
        return {
            "stage": stage,
            "feature": name,
            "min_A2": a2,
            "deep_core_A2_lt_m7_MHz": deep_mhz,
            "bad_core_A2_lt_m5_MHz": bad_mhz,
            "note": note,
        }

    history = [
        hist_row(
            3,
            "5.8125",
            s3["region_5p8"]["pass_6p25MHz"]["notch_A2"],
            6.25,
            12.5,
            "Stage-3 deep-core hint from 6.25 MHz localization",
        ),
        hist_row(
            4,
            "5.08125",
            s4["region_5p08"]["pass_3p125MHz"]["notch_A2"]
            if s4["region_5p08"].get("pass_3p125MHz")
            else s4["region_5p08"]["pass_6p25MHz"]["notch_A2"],
            s4["width_comparison"]["stage4_5p08_MHz"]["deep_core_MHz"],
            s4["width_comparison"]["stage4_5p08_MHz"]["bad_core_MHz"],
            "Stage-4 3.125 MHz localization where available",
        ),
        hist_row(
            4,
            "~4.73125",
            s4["region_4p73"]["pass_6p25MHz"]["notch_A2"],
            s4["width_comparison"]["stage4_4p73_MHz"]["deep_core_MHz"],
            s4["width_comparison"]["stage4_4p73_MHz"]["bad_core_MHz"],
            "Stage-4 6.25 MHz localization",
        ),
    ]
    for row in compact_rows:
        history.append(
            hist_row(
                5,
                f"{row['region']}@{row['min_f_GHz']:.6f}",
                row["min_A2"],
                row["A2_lt_m7_width_MHz"],
                row["A2_lt_m5_width_MHz"],
                f"Stage-5 q_best; spacing={row['authoritative_spacing_MHz']} MHz",
            )
        )

    # Narrowest reliably resolved important notch:
    # among Stage-5 features with A2<=-7 and deep_core_n_samples>=2,
    # take the minimum deep_core_width; if all single-sample even at ultra,
    # report ultra spacing as unresolved upper bound.
    resolved = []
    unresolved = []
    for row in compact_rows:
        if row["min_A2"] > -7:
            continue
        if row["A2_lt_m7_n_samples"] >= 2 and row["A2_lt_m7_width_MHz"] > 0:
            resolved.append(row)
        else:
            unresolved.append(row)

    if resolved:
        narrowest = min(resolved, key=lambda r: r["A2_lt_m7_width_MHz"])
        narrowest_statement = {
            "quantity": "narrowest_reliably_resolved_important_notch_width_MHz",
            "value_MHz": narrowest["A2_lt_m7_width_MHz"],
            "definition": (
                "Minimum contiguous A2<-7 sample-span width among Stage-5 "
                "important notches with >=2 samples below -7 (measured FEM)."
            ),
            "feature": narrowest,
        }
    else:
        # All deep cores single-sample → width unresolved; report sample spacing
        auth_spacings = [r["authoritative_spacing_MHz"] for r in compact_rows]
        narrowest_statement = {
            "quantity": "narrowest_reliably_resolved_important_notch_width_MHz",
            "value_MHz": None,
            "upper_bound_unresolved_core_MHz": min(auth_spacings) if auth_spacings else None,
            "definition": (
                "Deep cores remain single-sample even at finest spacing used; "
                "true A2<-7 width is unresolved and may be < that spacing."
            ),
            "unresolved_features": unresolved,
        }

    summary = {
        "design": str(design),
        "J_BB_best": d.get("J_BB_best"),
        "stage6_status": "NOT STARTED / AWAITING AUTHORIZATION",
        "note": (
            "Forward-only notch-width characterization on FINAL Stage-5 q_best. "
            "No optimization. Measured FEM samples are authoritative."
        ),
        "centers_GHz": CENTERS,
        "pass_1p5625_window_half_GHz": HALF_WINDOW_GHZ,
        "regions": {
            tag: {
                "pass_1p5625MHz": {
                    k: v
                    for k, v in blob["pass_1p5625MHz"].items()
                    if k != "table"
                },
                "pass_0p78125MHz": (
                    {k: v for k, v in blob["pass_0p78125MHz"].items() if k != "table"}
                    if blob["pass_0p78125MHz"]
                    else None
                ),
            }
            for tag, blob in regions.items()
        },
        "compact_width_table": compact_rows,
        "stages_3_to_5_width_comparison": history,
        "narrowest_reliably_resolved_important_notch": narrowest_statement,
        "prior_stage4_narrowest_MHz": s4.get(
            "narrowest_reliably_observed_bad_feature_width_MHz"
        ),
    }

    # Full tables separately (keep summary readable)
    full = {
        **summary,
        "region_tables": {
            tag: {
                "pass_1p5625MHz_table": blob["pass_1p5625MHz"]["table"],
                "pass_0p78125MHz_table": (
                    blob["pass_0p78125MHz"]["table"] if blob["pass_0p78125MHz"] else None
                ),
            }
            for tag, blob in regions.items()
        },
    }

    out_json = OUT / "BB_STAGE5_NOTCH_WIDTHS.json"
    out_json.write_text(json.dumps(full, indent=2) + "\n")
    compact_csv = OUT / "BB_STAGE5_NOTCH_WIDTHS.csv"
    with compact_csv.open("w", newline="") as f:
        fields = list(compact_rows[0].keys()) if compact_rows else []
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in compact_rows:
            w.writerow(r)

    hist_csv = OUT / "BB_STAGE5_NOTCH_WIDTHS_HISTORY.csv"
    with hist_csv.open("w", newline="") as f:
        fields = list(history[0].keys()) if history else []
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in history:
            w.writerow(r)

    print("\n=== COMPACT WIDTH TABLE ===", flush=True)
    for r in compact_rows:
        print(
            f"{r['region']}: f={r['min_f_GHz']:.6f} A2={r['min_A2']:.3f} "
            f"A2<-5={r['A2_lt_m5_width_MHz']:.3f} MHz (n={r['A2_lt_m5_n_samples']}) "
            f"A2<-7={r['A2_lt_m7_width_MHz']:.3f} MHz (n={r['A2_lt_m7_n_samples']}) "
            f"ultra={r['ultra_refined']}",
            flush=True,
        )
    print("narrowest:", json.dumps(narrowest_statement, indent=2), flush=True)
    print("wrote", out_json, flush=True)
    print("wrote", compact_csv, flush=True)
    print("STAGE-6: NOT STARTED / AWAITING AUTHORIZATION", flush=True)


if __name__ == "__main__":
    main()
