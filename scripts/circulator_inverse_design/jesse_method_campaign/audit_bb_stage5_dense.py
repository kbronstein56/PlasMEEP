#!/usr/bin/env python3
"""Post-Stage-5 spectral audit (forward-only).

Baseline: full-band 6.25 MHz (same as Stage-5 objective grid).
Then refine suspicious local minima (and historical failure neighborhoods)
to 3.125 MHz.

Does not change physics. Does not launch Stage-6.
"""
from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
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

from spectral_recon import ensure_summary_fields, summarize  # noqa: E402
from stage5_frequency_set import stage5_freqs_GHz  # noqa: E402

HISTORICAL = [4.73125, 4.775, 5.08125, 5.5875, 5.75, 5.8125, 5.825, 5.875]


def _init_worker(omp_threads: int):
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = str(int(omp_threads))


def _worker_eval(payload: dict) -> dict:
    _init_worker(int(payload["omp_threads"]))
    from eval_forward_only import eval_q_forward
    from multimode_opt_core import B_PLUS, orbits_and_masks
    from objective_A2 import objective_A2_from_S

    q = np.asarray(payload["q"], float)
    orbits, _, _ = orbits_and_masks()
    f_ghz = float(payload["f_GHz"])
    bun, _ = eval_q_forward(
        payload["grade"], q, orbits, f_hz=f_ghz * 1e9, b_tesla=B_PLUS
    )
    return row_from_bundle(f_ghz, bun)


def row_from_bundle(f_ghz: float, bun: dict) -> dict:
    from objective_A2 import objective_A2_from_S

    m = bun["metrics"]
    J_A2 = float(objective_A2_from_S(bun["S"]))
    per = m.get("per_source") or []
    if per:
        worst_iso_src = min(per, key=lambda r: r["isolation_dB"])
        j = int(worst_iso_src["source"])
        worst_competing = {
            "source_port": j + 1,
            "desired_port": (j + 1) % 6 + 1,
            "competing_port": (j - 1) % 6 + 1,
            "T_des": float(worst_iso_src["T_des"]),
            "T_rev": float(worst_iso_src["T_rev"]),
            "isolation_dB": float(worst_iso_src["isolation_dB"]),
            "IL_dB": float(worst_iso_src["IL_dB"]),
        }
    else:
        worst_competing = None
    return {
        "f_GHz": float(f_ghz),
        "J_A2": J_A2,
        "J_A": float(bun["J_A"]),
        "J_B": float(bun["J_B"]),
        "T_des": float(m["T_des_avg"]),
        "T_rev": float(m["T_rev_avg"]),
        "IL_dB": float(m["IL_avg_dB"]),
        "IL_worst_dB": float(m.get("IL_worst_dB", m["IL_avg_dB"])),
        "isolation_dB": float(m["isolation_avg_dB"]),
        "isolation_min_dB": float(m.get("isolation_min_dB", m["isolation_avg_dB"])),
        "RL_dB": float(m["RL_avg_dB"]),
        "R": float(m["R_avg"]),
        "TE1": float(m["TE1_avg"]),
        "worst_competing": worst_competing,
    }


def freq_grid(f0: float, f1: float, df: float) -> list[float]:
    n = int(round((f1 - f0) / df)) + 1
    return [round(f0 + i * df, 6) for i in range(n)]


def load_seed_rows(paths: list[Path]) -> dict[float, dict]:
    done: dict[float, dict] = {}
    for p in paths:
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        for r in d.get("rows", []):
            if "J_A2" not in r:
                continue
            key = round(float(r["f_GHz"]), 6)
            if key not in done or (
                done[key].get("worst_competing") is None
                and r.get("worst_competing") is not None
            ):
                done[key] = r
    return done


def write_checkpoint(tag: str, grade: str, rows: list[dict], path: Path):
    rows = sorted(rows, key=lambda r: r["f_GHz"])
    summary = summarize(tag, grade, rows)
    summary = ensure_summary_fields(summary, grade=grade)
    worst_a2 = min(rows, key=lambda r: r["J_A2"])
    worst_il = max(rows, key=lambda r: r["IL_dB"])
    worst_iso = min(rows, key=lambda r: r["isolation_dB"])
    out = {
        "tag": tag,
        "grade": grade,
        "rows": rows,
        **summary,
        "worst_A2": {
            "f_GHz": worst_a2["f_GHz"],
            "J_A2": worst_a2["J_A2"],
            "IL_dB": worst_a2["IL_dB"],
            "isolation_dB": worst_a2["isolation_dB"],
            "T_des": worst_a2["T_des"],
        },
        "worst_IL": {
            "f_GHz": worst_il["f_GHz"],
            "IL_dB": worst_il["IL_dB"],
            "J_A2": worst_il["J_A2"],
            "isolation_dB": worst_il["isolation_dB"],
        },
        "worst_isolation": {
            "f_GHz": worst_iso["f_GHz"],
            "isolation_dB": worst_iso["isolation_dB"],
            "J_A2": worst_iso["J_A2"],
            "IL_dB": worst_iso["IL_dB"],
        },
    }
    path.write_text(json.dumps(out, indent=2) + "\n")
    return out


def write_csv(rows: list[dict], path: Path):
    fields = [
        "f_GHz", "J_A2", "J_A", "J_B", "T_des", "T_rev", "IL_dB", "IL_worst_dB",
        "isolation_dB", "isolation_min_dB", "RL_dB", "R", "TE1",
    ]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in sorted(rows, key=lambda x: x["f_GHz"]):
            w.writerow({k: r.get(k) for k in fields})


def local_minima(rows: list[dict], key: str = "J_A2", pad: int = 1) -> list[dict]:
    rs = sorted(rows, key=lambda r: r["f_GHz"])
    out = []
    for i in range(pad, len(rs) - pad):
        v = rs[i][key]
        if all(v < rs[i + d][key] for d in range(-pad, pad + 1) if d != 0):
            out.append(rs[i])
    return out


def parallel_scan(grade, q, freqs, *, tag, resume_path, seed_paths, n_workers, omp_threads):
    done = load_seed_rows(seed_paths + ([resume_path] if resume_path.exists() else []))
    want = {round(float(f), 6) for f in freqs}
    rows = [done[k] for k in sorted(want) if k in done]
    missing = [f for f in freqs if round(float(f), 6) not in done]
    print(
        f"{tag}: {len(freqs)} requested, {len(rows)} seeded, {len(missing)} to solve "
        f"(n_workers={n_workers}, omp={omp_threads})",
        flush=True,
    )
    if not missing:
        out = write_checkpoint(tag, grade, rows, resume_path)
        write_csv(rows, resume_path.with_suffix(".csv"))
        return out

    thread_keys = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
    saved = {k: os.environ.get(k) for k in thread_keys}
    for k in thread_keys:
        os.environ[k] = str(int(omp_threads))
    t0 = time.perf_counter()
    payloads = [
        {"grade": grade, "q": np.asarray(q, float).tolist(), "f_GHz": float(f),
         "omp_threads": int(omp_threads)}
        for f in missing
    ]
    ctx = mp.get_context("spawn")
    try:
        with ProcessPoolExecutor(
            max_workers=max(1, int(n_workers)),
            mp_context=ctx,
            initializer=_init_worker,
            initargs=(int(omp_threads),),
        ) as ex:
            futs = {ex.submit(_worker_eval, p): p["f_GHz"] for p in payloads}
            n_done = 0
            for fut in as_completed(futs):
                row = fut.result()
                rows.append(row)
                n_done += 1
                print(
                    f"{tag} [{n_done}/{len(missing)}] f={row['f_GHz']:.5f} "
                    f"A2={row['J_A2']:.4f} IL={row['IL_dB']:.2f} "
                    f"iso={row['isolation_dB']:.2f} Tdes={row['T_des']:.4e}",
                    flush=True,
                )
                if n_done % 12 == 0 or n_done == len(missing):
                    write_checkpoint(tag, grade, rows, resume_path)
                    write_csv(rows, resume_path.with_suffix(".csv"))
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    print(f"{tag}: solved {len(missing)} freqs in {time.perf_counter()-t0:.1f}s", flush=True)
    out = write_checkpoint(tag, grade, rows, resume_path)
    write_csv(rows, resume_path.with_suffix(".csv"))
    return out


def seed_from_opt_rec(rec_best: dict | None) -> list[dict]:
    """Optional seed from Stage-5 objective-grid metrics (same frequencies)."""
    if not rec_best:
        return []
    freqs = rec_best.get("freqs_GHz") or []
    Js = rec_best.get("J_per_freq") or []
    mets = rec_best.get("metrics_per_freq") or []
    if not (freqs and Js and mets) or not (len(freqs) == len(Js) == len(mets)):
        return []
    rows = []
    for f, J, m in zip(freqs, Js, mets):
        rows.append(
            {
                "f_GHz": float(f),
                "J_A2": float(J),
                "J_A": float("nan"),
                "J_B": float("nan"),
                "T_des": float(m.get("T_des_avg", m.get("T_des", float("nan")))),
                "T_rev": float(m.get("T_rev_avg", m.get("T_rev", float("nan")))),
                "IL_dB": float(m.get("IL_avg_dB", m.get("IL_dB", float("nan")))),
                "IL_worst_dB": float(m.get("IL_worst_dB", m.get("IL_avg_dB", float("nan")))),
                "isolation_dB": float(
                    m.get("isolation_avg_dB", m.get("isolation_dB", float("nan")))
                ),
                "isolation_min_dB": float(
                    m.get("isolation_min_dB", m.get("isolation_avg_dB", float("nan")))
                ),
                "RL_dB": float(m.get("RL_avg_dB", m.get("RL_dB", float("nan")))),
                "R": float(m.get("R_avg", m.get("R", float("nan")))),
                "TE1": float(m.get("TE1_avg", m.get("TE1", float("nan")))),
                "worst_competing": None,
                "seeded_from_opt": True,
            }
        )
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--omp-threads-per-worker", type=int, default=1)
    ap.add_argument(
        "--stage",
        choices=("baseline625", "refine3125", "all"),
        default="all",
    )
    ap.add_argument("--refine-half-width-ghz", type=float, default=0.025)
    ap.add_argument("--a2-suspect-threshold", type=float, default=-3.0)
    ap.add_argument("--force-rescan-baseline", action="store_true")
    args = ap.parse_args()

    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(k, "1")

    design_path = OUT / "opt_broadband_stage5_A2_C.json"
    d5 = json.loads(design_path.read_text())
    q5 = np.asarray(d5["q_best"], float)
    print(
        f"Stage-5 dense audit design={design_path.name} "
        f"J_BB={d5.get('J_BB_best')} nit={d5.get('accepted_nit_total', d5.get('accepted_nit'))} "
        f"n_eval={d5.get('n_eval')}",
        flush=True,
    )

    report = {
        "design": str(design_path),
        "J_BB_best": d5.get("J_BB_best"),
        "accepted_nit": d5.get("accepted_nit_total", d5.get("accepted_nit")),
        "n_eval": d5.get("n_eval"),
        "n_workers": args.n_workers,
        "stage5_status": "audit_only — Stage-6 NOT STARTED",
    }

    a625 = None
    baseline_path = OUT / "bb_stage5_dense_6p25MHz.json"
    if args.stage in ("baseline625", "all"):
        freqs = stage5_freqs_GHz()
        assert len(freqs) == 321
        # Seed from opt rec_best when available (same grid); still re-solve if forced
        seed_paths = []
        if not args.force_rescan_baseline:
            opt_seed = OUT / "bb_stage5_opt_grid_seed.json"
            rows_seed = seed_from_opt_rec(d5.get("rec_best"))
            if rows_seed:
                write_checkpoint("bb_stage5_opt_grid_seed", args.grade, rows_seed, opt_seed)
                seed_paths.append(opt_seed)
                print(f"Seeded {len(rows_seed)} objective-grid rows from rec_best", flush=True)
        print(f"BASELINE 6.25 MHz: n={len(freqs)} {freqs[0]}–{freqs[-1]} GHz", flush=True)
        a625 = parallel_scan(
            args.grade,
            q5,
            freqs,
            tag="bb_stage5_dense_6p25MHz",
            resume_path=baseline_path,
            seed_paths=seed_paths,
            n_workers=args.n_workers,
            omp_threads=args.omp_threads_per_worker,
        )
        mins = local_minima(a625["rows"])
        report["dense_6p25MHz"] = {
            "n_freq": a625["n_freq"],
            "spacing_MHz": 6.25,
            "worst_A2": a625["worst_A2"],
            "worst_IL": a625["worst_IL"],
            "worst_isolation": a625["worst_isolation"],
            "local_A2_minima": [
                {
                    "f_GHz": r["f_GHz"],
                    "J_A2": r["J_A2"],
                    "IL_dB": r["IL_dB"],
                    "isolation_dB": r["isolation_dB"],
                }
                for r in mins
            ],
        }

    if args.stage in ("refine3125", "all"):
        if a625 is None:
            a625 = json.loads(baseline_path.read_text())
        centers = set(HISTORICAL)
        for r in local_minima(a625["rows"]):
            if r["J_A2"] <= args.a2_suspect_threshold:
                centers.add(float(r["f_GHz"]))
        # Always include true worst
        centers.add(float(a625["worst_A2"]["f_GHz"]))
        hw = float(args.refine_half_width_ghz)
        freqs3125 = set()
        for c in sorted(centers):
            for f in freq_grid(c - hw, c + hw, 0.003125):
                if 4.0 - 1e-12 <= f <= 6.0 + 1e-12:
                    freqs3125.add(round(f, 6))
        freqs3125 = sorted(freqs3125)
        print(
            f"REFINE 3.125 MHz: {len(centers)} centers → {len(freqs3125)} freqs "
            f"centers={sorted(centers)}",
            flush=True,
        )
        a3125 = parallel_scan(
            args.grade,
            q5,
            freqs3125,
            tag="bb_stage5_dense_3p125MHz_local",
            resume_path=OUT / "bb_stage5_dense_3p125MHz_local.json",
            seed_paths=[baseline_path],
            n_workers=args.n_workers,
            omp_threads=args.omp_threads_per_worker,
        )
        report["dense_3p125MHz_local"] = {
            "n_freq": a3125["n_freq"],
            "spacing_MHz": 3.125,
            "centers_GHz": sorted(centers),
            "worst_A2": a3125["worst_A2"],
            "worst_IL": a3125["worst_IL"],
            "worst_isolation": a3125["worst_isolation"],
            "local_A2_minima": [
                {
                    "f_GHz": r["f_GHz"],
                    "J_A2": r["J_A2"],
                    "IL_dB": r["IL_dB"],
                    "isolation_dB": r["isolation_dB"],
                }
                for r in local_minima(a3125["rows"])
            ],
        }

    # Historical watch table from baseline
    if a625 is not None:
        byf = {round(r["f_GHz"], 6): r for r in a625["rows"]}
        report["historical_watch"] = {}
        for f in HISTORICAL:
            r = byf.get(round(f, 6))
            if r:
                report["historical_watch"][str(f)] = {
                    "J_A2": r["J_A2"],
                    "IL_dB": r["IL_dB"],
                    "isolation_dB": r["isolation_dB"],
                    "T_des": r["T_des"],
                }

    # Compare vs Stage-3/4 dense audits if present
    comparisons = {}
    for label, path in (
        ("stage3", OUT / "BB_STAGE3_DENSE_AUDIT_SUMMARY.json"),
        ("stage4", OUT / "BB_STAGE4_DENSE_AUDIT_SUMMARY.json"),
    ):
        if path.exists():
            comparisons[label] = json.loads(path.read_text())
    report["prior_summaries_present"] = list(comparisons.keys())

    # Conclusion scaffold (filled/refined by results writer)
    true_worst = None
    if "dense_3p125MHz_local" in report:
        true_worst = report["dense_3p125MHz_local"]["worst_A2"]
    elif "dense_6p25MHz" in report:
        true_worst = report["dense_6p25MHz"]["worst_A2"]
    report["true_worst"] = true_worst

    out_path = OUT / "BB_STAGE5_DENSE_AUDIT_SUMMARY.json"
    out_path.write_text(json.dumps(report, indent=2) + "\n")
    print("wrote", out_path, flush=True)
    if true_worst:
        print(
            f"TRUE WORST A2={true_worst['J_A2']:.4f} @ {true_worst['f_GHz']:.5f} GHz "
            f"IL={true_worst['IL_dB']:.2f} iso={true_worst['isolation_dB']:.2f}",
            flush=True,
        )
    print("STAGE-6 STATUS: NOT STARTED", flush=True)


if __name__ == "__main__":
    main()
