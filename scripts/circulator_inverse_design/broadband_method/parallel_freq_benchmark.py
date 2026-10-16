#!/usr/bin/env python3
"""Phase 6: benchmark concurrent frequency eval+grad under RAM limits."""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "broadband_method"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]


def _worker(payload):
    # Limit BLAS/OpenMP oversubscription inside each worker
    os.environ["OMP_NUM_THREADS"] = payload.get("omp", "1")
    os.environ["OPENBLAS_NUM_THREADS"] = payload.get("omp", "1")
    os.environ["MKL_NUM_THREADS"] = payload.get("omp", "1")
    from efficient_adjoint import eval_and_grad_efficient
    from multimode_opt_core import orbits_and_masks

    grade = payload["grade"]
    q = np.asarray(payload["q"], float)
    f = payload["f_Hz"]
    orbits, _, _ = orbits_and_masks()
    t0 = time.perf_counter()
    bun, g, s = eval_and_grad_efficient(grade, q, orbits, f_hz=f)
    return {
        "f_Hz": f,
        "J": bun["J"],
        "elapsed_s": time.perf_counter() - t0,
        "timings": bun.get("timings"),
    }


def rss_mb():
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024.0
    except Exception:
        return float("nan")
    return float("nan")


def run_batch(grade, q, freqs, workers: int):
    payloads = [
        {"grade": grade, "q": q.tolist(), "f_Hz": float(f), "omp": "1"}
        for f in freqs
    ]
    t0 = time.perf_counter()
    results = []
    if workers <= 1:
        for p in payloads:
            results.append(_worker(p))
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(_worker, p) for p in payloads]
            for fut in as_completed(futs):
                results.append(fut.result())
    wall = time.perf_counter() - t0
    return results, wall


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument(
        "--design",
        default="outputs/circulator_inverse_design/multimode_optimization/opt_uniform_C.json",
    )
    ap.add_argument("--n-freq", type=int, default=4, help="subset of objective grid")
    ap.add_argument("--workers", default="1,2", help="comma list; 4 may OOM on C")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    d = json.loads(Path(args.design).read_text())
    q = np.asarray(d["best"]["q"], float)
    # sample first n-freq of 11-point grid
    freqs = [ (4.0 + 0.2 * i) * 1e9 for i in range(args.n_freq) ]
    worker_list = [int(x) for x in args.workers.split(",")]

    rows = []
    for w in worker_list:
        print(f"workers={w} freqs={len(freqs)} ...", flush=True)
        # RSS of parent only; child peak not fully captured — note limitation
        rss0 = rss_mb()
        results, wall = run_batch(args.grade, q, freqs, w)
        rss1 = rss_mb()
        per = float(np.mean([r["elapsed_s"] for r in results]))
        rows.append(
            {
                "workers": w,
                "n_freq": len(freqs),
                "wall_s": wall,
                "mean_per_freq_s": per,
                "parent_rss_mb_before": rss0,
                "parent_rss_mb_after": rss1,
                "Js": [r["J"] for r in results],
            }
        )
        print(f"  wall={wall:.1f}s mean/freq={per:.1f}s", flush=True)

    report = {
        "grade": args.grade,
        "design": args.design,
        "freqs_Hz": freqs,
        "note": (
            "Each frequency: 1 LU + 6 TE0 fwd + 6 adj (efficient path). "
            "OMP_NUM_THREADS=1 per worker to avoid BLAS oversubscription. "
            "Peak aggregate RAM for multi-process is ~workers × single-process footprint; "
            "parent VmRSS does not capture children."
        ),
        "runs": rows,
    }
    (OUT / "parallel_freq_benchmark.json").write_text(json.dumps(report, indent=2) + "\n")

    lines = [
        "# PARALLEL_FREQUENCY_STRATEGY",
        "",
        report["note"],
        "",
        "| workers | N_f | wall (s) | mean s/freq |",
        "|--------:|----:|---------:|------------:|",
    ]
    for r in rows:
        lines.append(
            f"| {r['workers']} | {r['n_freq']} | {r['wall_s']:.1f} | {r['mean_per_freq_s']:.1f} |"
        )
    best = min(rows, key=lambda r: r["wall_s"])
    lines += [
        "",
        f"**Recommended concurrency (this machine / grade {args.grade}):** "
        f"{best['workers']} worker(s) for the tested subset.",
        "",
        "For full 11-frequency objective+grad, expect roughly "
        f"`11 × (efficient single-freq time) / workers` wall seconds per L-BFGS step, "
        "plus orbit dA cost already inside each worker.",
        "",
    ]
    (OUT / "PARALLEL_FREQUENCY_STRATEGY.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"best_workers": best["workers"], "best_wall_s": best["wall_s"]}, indent=2))


if __name__ == "__main__":
    main()
