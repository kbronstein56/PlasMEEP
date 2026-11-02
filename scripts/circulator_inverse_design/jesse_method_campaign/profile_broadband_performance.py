#!/usr/bin/env python3
"""Broadband A2 performance audit for Stage-2 / Stage-3 planning.

Measures serial vs process-parallel frequency evals, OpenMP/BLAS thread
sensitivity, SuperLU path, CPU utilization proxies, and RAM.

Does NOT start Stage-3. Does NOT change physics.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
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

from a2_efficient_core import (  # noqa: E402
    eval_and_grad_A2,
    eval_bb_A2,
    rss_gib,
)
from multimode_opt_core import orbits_and_masks  # noqa: E402

# Stage-3 proposed frequency set (30)
F_STAGE3_GHZ = [
    4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.55, 4.6, 4.7, 4.8, 4.85, 4.9,
    5.0, 5.1, 5.2, 5.3, 5.35, 5.4, 5.5, 5.6, 5.625, 5.65, 5.7,
    5.75, 5.8, 5.825, 5.85, 5.9, 5.95, 6.0,
]
F_STAGE2_GHZ = [f for f in F_STAGE3_GHZ if f not in (5.75, 5.825)]


def set_threads(n: int):
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = str(n)


def cpu_sample(pid: int, dt: float = 1.0) -> dict:
    """Sample /proc/<pid>/stat for user+system jiffies over dt seconds."""
    def read():
        with open(f"/proc/{pid}/stat") as f:
            parts = f.read().split()
        # utime 14, stime 15 (1-indexed) → indices 13,14
        return int(parts[13]) + int(parts[14])

    try:
        clk = os.sysconf(os.sysconf_names["SC_CLK_TCK"])
        a = read()
        t0 = time.perf_counter()
        time.sleep(dt)
        b = read()
        wall = time.perf_counter() - t0
        cpu_s = (b - a) / float(clk)
        return {
            "sample_wall_s": wall,
            "cpu_seconds": cpu_s,
            "cpu_utilization_equiv_cores": cpu_s / wall if wall > 0 else float("nan"),
        }
    except Exception as e:
        return {"error": str(e)}


def load_q():
    p = OUT / "opt_broadband_stage2_A2_C.json"
    d = json.loads(p.read_text())
    return np.asarray(d["q_best"], float), d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--n-freq-probe", type=int, default=4,
                    help="subset size for parallel/thread sweeps")
    ap.add_argument("--full-serial-profile", action="store_true",
                    help="profile full Stage-2 28-freq serial once (expensive)")
    ap.add_argument("--skip-parallel", action="store_true")
    ap.add_argument("--skip-threads", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    q, opt = load_q()
    orbits, _, _ = orbits_and_masks()
    # Representative mid-band frequencies for short probes
    probe_ghz = [4.1, 5.0, 5.65, 5.75][: args.n_freq_probe]
    probe_hz = [f * 1e9 for f in probe_ghz]

    report: dict = {
        "machine": {
            "nproc": int(subprocess.check_output(["nproc"]).decode().strip()),
            "mem_GiB_total": None,
        },
        "design": "opt_broadband_stage2_A2_C.json q_best",
        "grade": args.grade,
        "stage2_historical": {
            "elapsed_s": opt.get("elapsed_s"),
            "n_eval": opt.get("n_eval") or opt.get("accepted_nit_total"),
            "accepted_nit": 30,
            "n_freq": 28,
            "note": "From completed Stage-2 run; ~n_eval BB evals × 28 serial freqs",
        },
        "solver": {
            "path": "scipy.sparse.linalg.splu → SuperLU",
            "umfpack": False,
            "pypardiso": False,
            "lu_per_frequency": 1,
            "forward_rhs_per_frequency": 6,
            "adjoint_rhs_per_frequency": 6,
            "frequencies_currently": "SERIAL in eval_bb_A2 for-loop",
        },
        "thread_env_at_start": {
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
            "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
            "MKL_NUM_THREADS": os.environ.get("MKL_NUM_THREADS"),
        },
    }
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    report["machine"]["mem_GiB_total"] = int(line.split()[1]) / (1024 ** 2)
                    break
    except Exception:
        pass

    # --- Single-frequency detailed profile (default OMP=8 as historically used) ---
    print("=== SF profile @ 5.65 GHz with OMP=8 ===", flush=True)
    set_threads(8)
    # Warm assemble path once
    bun, g, _ = eval_and_grad_A2(args.grade, q, orbits, 5.65e9, profile=True)
    # Fresh timed
    t0 = time.perf_counter()
    bun, g, _ = eval_and_grad_A2(args.grade, q, orbits, 5.65e9, profile=True)
    wall = time.perf_counter() - t0
    cpu = cpu_sample(os.getpid(), dt=0.5)  # after; also sample during next
    # Sample during LU-heavy call
    set_threads(8)
    import threading

    samples = []

    def _sampler():
        samples.append(cpu_sample(os.getpid(), dt=2.0))

    th = threading.Thread(target=_sampler)
    th.start()
    bun2, g2, _ = eval_and_grad_A2(args.grade, q, orbits, 5.65e9, profile=True)
    th.join()
    report["single_freq_profile_OMP8"] = {
        "f_GHz": 5.65,
        "wall_s": wall,
        "timings": bun["timings"],
        "grad_norm": float(np.linalg.norm(g)),
        "J": bun["J"],
        "cpu_sample_during_second_call": samples[0] if samples else None,
        "note": (
            "OMP/OPENBLAS/MKL affect dense BLAS inside SuperLU/SciPy where used; "
            "SuperLU factorization itself is largely serial. Measure equiv_cores "
            "during the call — not just env vars."
        ),
    }
    print(json.dumps(report["single_freq_profile_OMP8"]["timings"], indent=2), flush=True)
    print("CPU sample:", samples, flush=True)

    # --- Thread sweep on single frequency ---
    thread_rows = []
    if not args.skip_threads:
        for nth in [1, 2, 4, 8]:
            print(f"=== thread sweep OMP={nth} @ 5.0 GHz ===", flush=True)
            set_threads(nth)
            # one warmup
            eval_and_grad_A2(args.grade, q, orbits, 5.0e9, profile=False)
            samples = []
            th = threading.Thread(target=lambda: samples.append(cpu_sample(os.getpid(), 3.0)))
            th.start()
            t0 = time.perf_counter()
            bun, g, _ = eval_and_grad_A2(args.grade, q, orbits, 5.0e9, profile=True)
            wall = time.perf_counter() - t0
            th.join()
            row = {
                "omp_threads": nth,
                "wall_s": wall,
                "timings": bun["timings"],
                "cpu_sample": samples[0] if samples else None,
                "rss_GiB": bun["timings"].get("rss_GiB_end"),
            }
            thread_rows.append(row)
            print(
                f"  wall={wall:.1f}s factor={bun['timings'].get('factor_lu'):.1f} "
                f"cores≈{(samples[0] or {}).get('cpu_utilization_equiv_cores')}",
                flush=True,
            )
    report["thread_sweep_single_freq"] = thread_rows

    # --- Process-parallel probe ---
    parallel_rows = []
    if not args.skip_parallel:
        set_threads(1)  # parent; workers set their own
        # Serial baseline on probe set
        print(f"=== serial BB probe n_freq={len(probe_hz)} ===", flush=True)
        ser = eval_bb_A2(
            args.grade, q, orbits, probe_hz, n_workers=1, profile=True
        )
        parallel_rows.append(
            {
                "n_workers": 1,
                "omp_per_worker": "parent_threads",
                "wall_s": ser["wall_s"],
                "J_BB": ser["J_BB"],
                "grad_norm": ser["grad_norm"],
                "J_per_freq": ser["J_per_freq"],
                "mode": ser["mode"],
                "rss_parent_GiB": rss_gib(),
            }
        )
        print(f"  serial wall={ser['wall_s']:.1f}s J_BB={ser['J_BB']:.6f}", flush=True)

        for nw in [2, 4]:
            print(f"=== process-parallel workers={nw} omp/worker=1 ===", flush=True)
            par = eval_bb_A2(
                args.grade,
                q,
                orbits,
                probe_hz,
                n_workers=nw,
                omp_threads_per_worker=1,
                profile=True,
            )
            # Equivalence vs serial
            j_err = abs(par["J_BB"] - ser["J_BB"]) / max(abs(ser["J_BB"]), 1e-12)
            g_err = float(
                np.linalg.norm(np.asarray(par["grad_q"]) - np.asarray(ser["grad_q"]))
                / max(np.linalg.norm(ser["grad_q"]), 1e-12)
            )
            jp_err = float(
                np.max(np.abs(np.asarray(par["J_per_freq"]) - np.asarray(ser["J_per_freq"])))
            )
            row = {
                "n_workers": nw,
                "omp_per_worker": 1,
                "wall_s": par["wall_s"],
                "J_BB": par["J_BB"],
                "grad_norm": par["grad_norm"],
                "J_per_freq": par["J_per_freq"],
                "speedup_vs_serial": ser["wall_s"] / par["wall_s"] if par["wall_s"] else None,
                "rel_err_J_BB": j_err,
                "rel_err_grad": g_err,
                "max_abs_err_J_per_freq": jp_err,
                "pass_equivalence": bool(j_err < 1e-10 and g_err < 1e-9 and jp_err < 1e-10),
                "mode": par["mode"],
                "rss_parent_GiB": rss_gib(),
            }
            parallel_rows.append(row)
            print(
                f"  wall={par['wall_s']:.1f}s speedup={row['speedup_vs_serial']:.2f}x "
                f"eq={row['pass_equivalence']} j_err={j_err:.2e} g_err={g_err:.2e}",
                flush=True,
            )

        # Optional 6 workers if RAM looks safe (~rss_after_factor * workers)
        sf_rss = (report.get("single_freq_profile_OMP8") or {}).get("timings", {}).get(
            "rss_GiB_after_factor"
        )
        if sf_rss and sf_rss * 6 < 0.4 * (report["machine"].get("mem_GiB_total") or 200):
            nw = 6
            print(f"=== process-parallel workers={nw} omp/worker=1 ===", flush=True)
            par = eval_bb_A2(
                args.grade, q, orbits, probe_hz, n_workers=nw, omp_threads_per_worker=1
            )
            j_err = abs(par["J_BB"] - ser["J_BB"]) / max(abs(ser["J_BB"]), 1e-12)
            g_err = float(
                np.linalg.norm(np.asarray(par["grad_q"]) - np.asarray(ser["grad_q"]))
                / max(np.linalg.norm(ser["grad_q"]), 1e-12)
            )
            parallel_rows.append(
                {
                    "n_workers": nw,
                    "omp_per_worker": 1,
                    "wall_s": par["wall_s"],
                    "J_BB": par["J_BB"],
                    "speedup_vs_serial": ser["wall_s"] / par["wall_s"],
                    "rel_err_J_BB": j_err,
                    "rel_err_grad": g_err,
                    "pass_equivalence": bool(j_err < 1e-10 and g_err < 1e-9),
                    "mode": par["mode"],
                }
            )
            print(
                f"  wall={par['wall_s']:.1f}s speedup={ser['wall_s']/par['wall_s']:.2f}x eq={j_err<1e-10}",
                flush=True,
            )
    report["parallel_probe"] = {
        "freqs_GHz": probe_ghz,
        "rows": parallel_rows,
    }

    # --- Optional full Stage-2 serial profile ---
    if args.full_serial_profile:
        print("=== FULL Stage-2 28-freq serial profile (expensive) ===", flush=True)
        set_threads(8)
        freqs = [f * 1e9 for f in F_STAGE2_GHZ]
        full = eval_bb_A2(args.grade, q, orbits, freqs, n_workers=1, profile=True)
        report["full_stage2_serial_profile"] = {
            "n_freq": 28,
            "wall_s": full["wall_s"],
            "s_per_freq": full["wall_s"] / 28,
            "J_BB": full["J_BB"],
            "per_freq_timings": full.get("per_freq_timings"),
        }
        print(f"  full wall={full['wall_s']:.1f}s ({full['wall_s']/28:.1f}s/freq)", flush=True)

    # --- Stage-3 estimates ---
    # Prefer measured SF time from OMP=8 profile; scale
    sf = report["single_freq_profile_OMP8"]["timings"]["total"]
    # Historical Stage-2: 281238s / 83 evals / 28 freqs ≈ 121 s/freq
    hist_s_per_freq = None
    if opt.get("elapsed_s") and opt.get("n_eval"):
        hist_s_per_freq = float(opt["elapsed_s"]) / float(opt["n_eval"]) / 28.0
    elif opt.get("elapsed_s"):
        hist_s_per_freq = float(opt["elapsed_s"]) / 83.0 / 28.0  # from checkpoint knowledge

    best_par = None
    for r in parallel_rows:
        if r.get("n_workers", 1) > 1 and r.get("pass_equivalence"):
            if best_par is None or r["wall_s"] < best_par["wall_s"]:
                best_par = r
    speedup = (best_par["speedup_vs_serial"] if best_par else 1.0) or 1.0

    n_freq_s3 = 30
    # Conservatively assume ~2.5 evals per accepted iter (line search) × 30 iters + start
    # Stage-2 had 83 evals / 30 accepted ≈ 2.77
    evals_per_accepted = 83 / 30.0
    n_evals_s3 = 1 + 30 * evals_per_accepted

    s_per_freq = hist_s_per_freq or sf
    wall_serial_s3 = n_evals_s3 * n_freq_s3 * s_per_freq
    wall_fast_s3 = wall_serial_s3 / speedup

    report["stage3_estimate"] = {
        "n_freq": n_freq_s3,
        "freqs_GHz": F_STAGE3_GHZ,
        "s_per_freq_basis": s_per_freq,
        "s_per_freq_source": "stage2_historical" if hist_s_per_freq else "sf_profile",
        "assumed_n_evals": n_evals_s3,
        "evals_per_accepted_from_stage2": evals_per_accepted,
        "wall_hours_current_serial": wall_serial_s3 / 3600.0,
        "best_verified_parallel": best_par,
        "wall_hours_with_best_parallel": wall_fast_s3 / 3600.0,
        "speedup_applied": speedup,
        "note": "Estimate only; line-search eval count varies. Stage-3 NOT launched.",
    }

    # Classification
    report["work_classification"] = {
        "frequencies": "independent — SAFE to process-parallel",
        "within_frequency_assemble": "serial (mesh/materials)",
        "within_frequency_SuperLU_factor": "largely serial (SuperLU default)",
        "within_frequency_forward_6RHS": "serial backsolves after factor; could batch",
        "within_frequency_adjoint_6RHS": "serial backsolves after factor; could batch",
        "within_frequency_dA_orbits": "serial over 16 orbits; threaded BLAS may help mats",
        "OMP_OPENBLAS_MKL": (
            "Affect dense kernels / some sparse internals; do NOT treat as "
            "frequency parallelism. High OMP with process workers causes oversubscription."
        ),
    }

    report["recommendation"] = {
        "production_default": "serial (n_workers=1) until operator enables parallel",
        "safe_fast_config": {
            "n_workers": (best_par or {}).get("n_workers", 1),
            "omp_threads_per_worker": 1,
            "parent_OMP": 1,
            "enabled_via": "--n-workers N on production_broadband_opt_stage2.py",
        },
        "do_not": "blindly set OMP=64 or n_workers=64",
        "stage3_status": "NOT STARTED / AWAITING AUTHORIZATION",
    }

    out_json = OUT / "BROADBAND_PERFORMANCE_AUDIT.json"
    out_json.write_text(json.dumps(report, indent=2) + "\n")
    print("wrote", out_json, flush=True)
    print("STAGE3 estimate hours serial/fast:",
          report["stage3_estimate"]["wall_hours_current_serial"],
          report["stage3_estimate"]["wall_hours_with_best_parallel"], flush=True)


if __name__ == "__main__":
    main()
