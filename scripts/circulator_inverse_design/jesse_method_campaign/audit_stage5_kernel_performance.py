#!/usr/bin/env python3
"""Stage-5 kernel performance campaign: profile, equivalence, SF + 48-freq benches.

Does NOT launch Stage-5. Physics/objective unchanged.
Writes STAGE5_KERNEL_PERFORMANCE_AUDIT.{md,json} companions under outputs/.../
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
sys.path[:0] = [
    str(Path(__file__).resolve().parent),
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = "1"

from a2_efficient_core import eval_and_grad_A2, eval_bb_A2, rss_gib  # noqa: E402
from fem_assembly_fast import clear_caches  # noqa: E402
from multimode_opt_core import orbits_and_masks  # noqa: E402


def peak_rss_gib() -> float:
    # Linux: ru_maxrss is KiB
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0**2)


def load_q():
    return np.asarray(
        json.loads((OUT / "opt_broadband_stage4_A2_C.json").read_text())["q_best"],
        float,
    )


def probe_48():
    # Same representative set as BROADBAND_PARALLEL_SCALING_48FREQ.json
    return [4.0 + i * (2.0 / 47.0) for i in range(48)]


def compare_sf(grade, q, orbits, f_ghz, warmup=True):
    """Compare reference vs fast. Report warm fast time (production-relevant)."""
    f_hz = f_ghz * 1e9

    # Reference (no durable cache)
    clear_caches()
    if warmup:
        eval_and_grad_A2(grade, q, orbits, f_hz, use_fast_kernel=False, profile=False)
    t0 = time.perf_counter()
    bun_r, g_r, _ = eval_and_grad_A2(
        grade, q, orbits, f_hz, use_fast_kernel=False, profile=True, return_fields=True
    )
    wall_r = time.perf_counter() - t0
    rss_r = rss_gib()
    tim_r = bun_r["timings"]

    # Fast: one cold build (cache fill) then timed warm eval
    clear_caches()
    t_cold0 = time.perf_counter()
    bun_cold, _, _ = eval_and_grad_A2(
        grade, q, orbits, f_hz, use_fast_kernel=True, profile=True
    )
    wall_cold = time.perf_counter() - t_cold0
    t0 = time.perf_counter()
    bun_f, g_f, _ = eval_and_grad_A2(
        grade, q, orbits, f_hz, use_fast_kernel=True, profile=True, return_fields=True
    )
    wall_f = time.perf_counter() - t0
    rss_f = rss_gib()
    tim_f = bun_f["timings"]

    S_r, S_f = bun_r["S"], bun_f["S"]
    U_r, U_f = bun_r["U_te0"], bun_f["U_te0"]
    dJ = abs(bun_f["J"] - bun_r["J"])
    relJ = dJ / max(abs(bun_r["J"]), 1e-30)
    dg = np.linalg.norm(g_f - g_r)
    relg = dg / max(np.linalg.norm(g_r), 1e-30)
    dS = float(np.max(np.abs(S_f - S_r)))
    dU = float(np.max(np.abs(U_f - U_r)))
    return {
        "f_GHz": f_ghz,
        "J_ref": float(bun_r["J"]),
        "J_fast": float(bun_f["J"]),
        "abs_dJ": float(dJ),
        "rel_dJ": float(relJ),
        "grad_norm_ref": float(np.linalg.norm(g_r)),
        "abs_dgrad": float(dg),
        "rel_dgrad": float(relg),
        "max_abs_dS": dS,
        "max_abs_dU": dU,
        "wall_ref_s": wall_r,
        "wall_fast_cold_s": wall_cold,
        "wall_fast_warm_s": wall_f,
        "wall_fast_s": wall_f,  # production-relevant = warm
        "speedup": wall_r / max(wall_f, 1e-9),
        "speedup_vs_cold": wall_r / max(wall_cold, 1e-9),
        "rss_end_ref_GiB": rss_r,
        "rss_end_fast_GiB": rss_f,
        "timings_ref": tim_r,
        "timings_fast_warm": tim_f,
        "timings_fast_cold": bun_cold["timings"],
        "pass": bool(relJ < 1e-10 and relg < 1e-8 and dS < 1e-8),
    }


def bench_48(grade, q, orbits, use_fast: bool, n_workers=12):
    clear_caches()
    freqs = [f * 1e9 for f in probe_48()]
    # one-freq warmup in-process so first parallel worker cold-start is comparable
    eval_and_grad_A2(
        grade, q, orbits, freqs[0], use_fast_kernel=use_fast, profile=False
    )
    t0 = time.perf_counter()
    out = eval_bb_A2(
        grade,
        q,
        orbits,
        freqs,
        tau=0.05,
        alpha=0.25,
        n_workers=n_workers,
        omp_threads_per_worker=1,
        profile=True,
        use_fast_kernel=use_fast,
    )
    wall = time.perf_counter() - t0
    return {
        "use_fast_kernel": use_fast,
        "n_workers": n_workers,
        "n_freq": len(freqs),
        "wall_s": float(wall),
        "wall_reported_s": float(out["wall_s"]),
        "J_BB": float(out["J_BB"]),
        "grad_norm": float(out["grad_norm"]),
        "grad_q": np.asarray(out["grad_q"], float).tolist(),
        "J_per_freq": out["J_per_freq"],
        "peak_rss_parent_GiB": peak_rss_gib(),
        "rss_now_GiB": rss_gib(),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--skip-48", action="store_true")
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument(
        "--freqs",
        default="4.5,4.73125,5.08125,5.8125",
        help="single-freq GHz list",
    )
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    q = load_q()
    orbits, _, _ = orbits_and_masks()
    freqs = [float(x) for x in args.freqs.split(",")]

    report = {
        "status": "RUNNING",
        "stage5": "NOT STARTED / AWAITING AUTHORIZATION",
        "grade": args.grade,
        "omp": 1,
        "q_source": "opt_broadband_stage4_A2_C.json",
        "single_freq": [],
        "parallel_48": {},
    }

    print("=== single-frequency reference vs fast ===", flush=True)
    for i, f in enumerate(freqs):
        print(f"-- {f} GHz --", flush=True)
        rec = compare_sf(args.grade, q, orbits, f, warmup=(i == 0))
        print(
            f"  pass={rec['pass']} speedup={rec['speedup']:.2f}x "
            f"ref={rec['wall_ref_s']:.1f}s fast_warm={rec['wall_fast_warm_s']:.1f}s "
            f"fast_cold={rec['wall_fast_cold_s']:.1f}s "
            f"relJ={rec['rel_dJ']:.3e} relg={rec['rel_dgrad']:.3e} dS={rec['max_abs_dS']:.3e}",
            flush=True,
        )
        # Drop bulky U from stored JSON (timings only)
        report["single_freq"].append(rec)

    sf_pass = all(r["pass"] for r in report["single_freq"])
    report["single_freq_all_pass"] = sf_pass

    if not args.skip_48:
        if not sf_pass:
            print("SKIP 48-freq: single-freq equivalence FAILED", flush=True)
            report["parallel_48"] = {"skipped": True, "reason": "sf_fail"}
        else:
            print("=== 48-freq parallel n_workers=12 reference ===", flush=True)
            ref48 = bench_48(args.grade, q, orbits, use_fast=False, n_workers=args.n_workers)
            print(
                f"  ref wall={ref48['wall_s']:.1f}s J={ref48['J_BB']:.8f} "
                f"||g||={ref48['grad_norm']:.6e}",
                flush=True,
            )
            print("=== 48-freq parallel n_workers=12 fast ===", flush=True)
            fast48 = bench_48(args.grade, q, orbits, use_fast=True, n_workers=args.n_workers)
            print(
                f"  fast wall={fast48['wall_s']:.1f}s J={fast48['J_BB']:.8f} "
                f"||g||={fast48['grad_norm']:.6e}",
                flush=True,
            )
            dJ = abs(fast48["J_BB"] - ref48["J_BB"])
            relJ = dJ / max(abs(ref48["J_BB"]), 1e-30)
            g_r = np.asarray(ref48["grad_q"], float)
            g_f = np.asarray(fast48["grad_q"], float)
            dg = float(np.linalg.norm(g_f - g_r))
            relg = dg / max(np.linalg.norm(g_r), 1e-30)
            speedup = ref48["wall_s"] / max(fast48["wall_s"], 1e-9)
            # Stage-5 cost from measured 48-freq wall at 12 workers
            t_batch = fast48["wall_s"] / max(int(np.ceil(48 / args.n_workers)), 1)
            n_freq_s5 = 321
            n_batches = int(np.ceil(n_freq_s5 / args.n_workers))
            t_eval = n_batches * (fast48["wall_s"] / max(int(np.ceil(48 / args.n_workers)), 1))
            # More direct: scale by batch count ratio
            batches_48 = int(np.ceil(48 / args.n_workers))
            batches_321 = int(np.ceil(321 / args.n_workers))
            t_per_eval = fast48["wall_s"] * (batches_321 / batches_48)
            report["parallel_48"] = {
                "ref": ref48,
                "fast": fast48,
                "abs_dJ_BB": float(dJ),
                "rel_dJ_BB": float(relJ),
                "abs_dgrad": dg,
                "rel_dgrad": relg,
                "speedup": float(speedup),
                "pass": bool(relJ < 1e-10 and relg < 1e-8),
                "stage5_estimate": {
                    "n_freq": 321,
                    "n_workers": args.n_workers,
                    "batches_per_eval": batches_321,
                    "wall_per_eval_s": float(t_per_eval),
                    "wall_per_eval_h": float(t_per_eval / 3600.0),
                    "total_h_50": float(50 * t_per_eval / 3600.0),
                    "total_h_60": float(60 * t_per_eval / 3600.0),
                    "total_h_74": float(74 * t_per_eval / 3600.0),
                    "method": "scale measured 48-freq/12w wall by ceil(321/12)/ceil(48/12)",
                },
            }
            print(
                f"  48-eq pass={report['parallel_48']['pass']} speedup={speedup:.2f}x "
                f"Stage5 ~{t_per_eval/3600:.2f} h/eval",
                flush=True,
            )

    # Decision
    mean_speedup = float(np.mean([r["speedup"] for r in report["single_freq"]]))
    p48 = report.get("parallel_48") or {}
    end_to_end = float(p48.get("speedup") or 0.0)
    if not report.get("single_freq_all_pass", False) or (
        p48 and not p48.get("skipped") and not p48.get("pass", False)
    ):
        decision = "C"
        decision_text = (
            "optimization attempt risks changing scientific behavior — "
            "revert and retain validated implementation"
        )
    elif end_to_end >= 1.5 or mean_speedup >= 1.5:
        decision = "A"
        decision_text = (
            "meaningful implementation speedup achieved — use optimized kernel for Stage-5"
        )
    else:
        decision = "B"
        decision_text = (
            "only modest speedup — Stage-5 remains approximately multi-day and should be "
            "considered as such"
        )
    report["decision"] = decision
    report["decision_text"] = decision_text
    report["mean_sf_speedup"] = mean_speedup
    report["status"] = "COMPLETE"

    json_path = OUT / "STAGE5_KERNEL_PERFORMANCE_AUDIT.json"
    json_path.write_text(json.dumps(report, indent=2))
    print("wrote", json_path, flush=True)
    print("DECISION", decision, decision_text, flush=True)


if __name__ == "__main__":
    main()
