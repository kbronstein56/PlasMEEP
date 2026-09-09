#!/usr/bin/env python3
"""
MPI scaling benchmark: one device-port excitation of the six-port geometry.

Times a single source-port Meep run (after a one-port normalization) at fixed
res/run_time while varying MPI ranks. Uses OMP_NUM_THREADS=1 by default.

Does not change production physics. Prefer extending this harness over adding
new one-off MPI scripts.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VAL = os.path.dirname(os.path.abspath(__file__))
if VAL not in sys.path:
    sys.path.insert(0, VAL)
if os.path.join(ROOT, "scripts") not in sys.path:
    sys.path.insert(0, os.path.join(ROOT, "scripts"))


WORKER = r'''
import json, os, sys, time
sys.path.insert(0, os.environ["SIXPORT_VAL"])
sys.path.insert(0, os.environ["SIXPORT_SCRIPTS"])
import meep as mp
from sixport_common import (
    ensure_normalizations,
    default_uniform_rho,
    build_circulator_device,
    simulate_circulator,
    physical_resolution_report,
)
from port_formulations import get_formulation

res = int(os.environ["BENCH_RES"])
run_time = float(os.environ["BENCH_RUNTIME"])
port = int(os.environ["BENCH_PORT"])
formul = os.environ.get("BENCH_FORMULATION", "te1_hz_line")

t0 = time.time()
cache = ensure_normalizations(
    res,
    run_time=run_time,
    force=True,
    ports=[port],
    verbose=False,
    formulation=formul,
)
t_norm = time.time() - t0

rho = default_uniform_rho()
t1 = time.time()
# Single-port device excitation via the same path as reciprocity studies
result = simulate_circulator(
    rho,
    [0.0, 0.0, 0.0],
    res=res,
    run_time=run_time,
    verbose=False,
    incident_cache=cache,
    ports=[port],
    formulation=formul,
)
t_dev = time.time() - t1
flux = float(abs(result["power_matrix"][0, 0]))

out = {
    "rank_env": int(os.environ.get("OMPI_COMM_WORLD_SIZE", os.environ.get("PMI_SIZE", "1"))),
    "t_norm_s": t_norm,
    "t_device_s": t_dev,
    "flux": flux,
    "res": res,
    "run_time": run_time,
    "port": port,
    "formulation": formul,
    "physical_resolution": physical_resolution_report(res),
}
path = os.environ["BENCH_OUT"]
rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", os.environ.get("PMI_RANK", "0")))
if rank == 0:
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out))
'''


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, default=48)
    parser.add_argument("--run-time", type=float, default=30.0)
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--ranks", type=str, default="4,8,16,32")
    parser.add_argument("--port-formulation", type=str, default="te1_hz_line")
    parser.add_argument("--omp-threads", type=int, default=1)
    parser.add_argument("--json-out", type=str, default="")
    parser.add_argument(
        "--mpirun",
        type=str,
        default="/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/mpirun",
    )
    parser.add_argument(
        "--python",
        type=str,
        default="/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/python",
    )
    parser.add_argument(
        "--mode",
        type=str,
        default="oneport",
        choices=["oneport", "reciprocity_pair", "throughput"],
        help="oneport=device-port worker; reciprocity_pair=P1P2 study; throughput=parallel pairs",
    )
    parser.add_argument("--pair-ports", type=str, default="0,1")
    parser.add_argument(
        "--device-mode",
        type=str,
        default="full",
        choices=["full", "horns_only"],
        help="Device build for reciprocity_pair / throughput modes",
    )
    args = parser.parse_args()

    if args.mode == "oneport":
        return _run_oneport_benchmark(args)

    if args.mode == "reciprocity_pair":
        return _run_reciprocity_pair_benchmark(args)

    return _run_throughput_benchmark(args)


def _mpi_env_base(omp_threads: int) -> dict:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(omp_threads)
    env["MKL_NUM_THREADS"] = str(omp_threads)
    env["OPENBLAS_NUM_THREADS"] = str(omp_threads)
    env.setdefault("FI_PROVIDER", "tcp")
    env.setdefault("MPICH_CH4_NETMOD", "ofi")
    env.setdefault("UCX_TLS", "tcp,self")
    env["PYTHONPATH"] = f"{VAL}:{os.path.join(ROOT, 'scripts')}:" + env.get(
        "PYTHONPATH", ""
    )
    return env


def _run_reciprocity_pair_benchmark(args) -> int:
    """Time a full P1↔P2 reciprocity_b0_study (norm + device) at given res."""
    out_dir = os.path.join(ROOT, "outputs", "validation", "mpi_scale")
    os.makedirs(out_dir, exist_ok=True)
    ranks_list = [int(x) for x in args.ranks.split(",") if x.strip()]
    json_out = args.json_out or os.path.join(
        out_dir,
        f"reciprocity_pair_res{args.res}_rt{args.run_time:g}_{args.port_formulation}.json",
    )
    results = []
    study = os.path.join(VAL, "reciprocity_b0_study.py")
    for n in ranks_list:
        tag = f"mpi_pair_r{n}_res{args.res}"
        out_json = os.path.join(out_dir, f"{tag}.json")
        log_path = os.path.join(out_dir, f"{tag}.log")
        cmd = [
            args.mpirun,
            "-np",
            str(n),
            args.python,
            study,
            "--res",
            str(args.res),
            "--run-time",
            str(args.run_time),
            "--ports",
            args.pair_ports,
            "--port-formulation",
            args.port_formulation,
            "--device-mode",
            args.device_mode,
            "--label",
            tag,
            "--json-out",
            out_json,
            "--mpi-note",
            f"np={n}",
            "--skip-norm-if-cached",
        ]
        print(f"\n=== reciprocity_pair ranks={n} res={args.res} ===", flush=True)
        t0 = time.time()
        env = _mpi_env_base(args.omp_threads)
        with open(log_path, "w", encoding="utf-8") as log:
            log.write(f"# {' '.join(cmd)}\n")
            proc = subprocess.run(cmd, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        wall = time.time() - t0
        row = {"ranks": n, "wall_s": wall, "returncode": proc.returncode}
        if os.path.isfile(out_json):
            with open(out_json, encoding="utf-8") as f:
                data = json.load(f)
            row["total_s"] = data.get("timings_s", {}).get("total")
            row["P12_dB"] = data.get("reciprocity", {}).get("max_abs_diff_dB")
        results.append(row)
        print(f"  wall={wall:.1f}s total={row.get('total_s')}", flush=True)

    summary = {
        "mode": "reciprocity_pair",
        "res": args.res,
        "run_time": args.run_time,
        "formulation": args.port_formulation,
        "results": results,
    }
    with open(json_out, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"Wrote {json_out}")
    return 0


def _run_throughput_benchmark(args) -> int:
    """Two independent reciprocity jobs in parallel (e.g. 2×np=32)."""
    import concurrent.futures

    out_dir = os.path.join(ROOT, "outputs", "validation", "mpi_scale")
    os.makedirs(out_dir, exist_ok=True)
    ranks = int(args.ranks.split(",")[0])
    study = os.path.join(VAL, "reciprocity_b0_study.py")

    def _one_job(job_id: int, ports: str) -> dict:
        tag = f"throughput_j{job_id}_r{ranks}_res{args.res}"
        out_json = os.path.join(out_dir, f"{tag}.json")
        log_path = os.path.join(out_dir, f"{tag}.log")
        cmd = [
            args.mpirun,
            "-np",
            str(ranks),
            args.python,
            study,
            "--res",
            str(args.res),
            "--run-time",
            str(args.run_time),
            "--ports",
            ports,
            "--port-formulation",
            args.port_formulation,
            "--device-mode",
            args.device_mode,
            "--label",
            tag,
            "--json-out",
            out_json,
            "--skip-norm-if-cached",
        ]
        env = _mpi_env_base(args.omp_threads)
        t0 = time.time()
        with open(log_path, "w", encoding="utf-8") as log:
            proc = subprocess.run(cmd, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        return {
            "job": job_id,
            "ports": ports,
            "wall_s": time.time() - t0,
            "returncode": proc.returncode,
            "path": out_json,
        }

    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        futs = [
            ex.submit(_one_job, 0, "0,1"),
            ex.submit(_one_job, 1, "1,2"),
        ]
        jobs = [f.result() for f in futs]
    wall = time.time() - t0
    summary = {
        "mode": "throughput",
        "ranks_per_job": ranks,
        "res": args.res,
        "run_time": args.run_time,
        "wall_s_parallel": wall,
        "jobs": jobs,
    }
    json_out = args.json_out or os.path.join(
        out_dir, f"throughput_2xnp{ranks}_res{args.res}.json"
    )
    with open(json_out, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"parallel wall={wall:.1f}s  jobs={jobs}")
    print(f"Wrote {json_out}")
    return 0


def _run_oneport_benchmark(args) -> int:
    out_dir = os.path.join(ROOT, "outputs", "validation", "mpi_scale")
    os.makedirs(out_dir, exist_ok=True)
    worker_path = os.path.join(out_dir, "_mpi_worker.py")
    with open(worker_path, "w", encoding="utf-8") as handle:
        handle.write(WORKER)

    ranks = [int(x) for x in args.ranks.split(",") if x.strip()]
    results = []
    env_base = os.environ.copy()
    env_base["SIXPORT_VAL"] = VAL
    env_base["SIXPORT_SCRIPTS"] = os.path.join(ROOT, "scripts")
    env_base["BENCH_RES"] = str(args.res)
    env_base["BENCH_RUNTIME"] = str(args.run_time)
    env_base["BENCH_PORT"] = str(args.port)
    env_base["BENCH_FORMULATION"] = args.port_formulation
    env_base["OMP_NUM_THREADS"] = str(args.omp_threads)
    env_base["MKL_NUM_THREADS"] = str(args.omp_threads)
    env_base["OPENBLAS_NUM_THREADS"] = str(args.omp_threads)
    # Prefer OFI/TCP on WSL: default UCX shared-memory hits tiny /dev/shm limits
    # at high rank counts unless the host remounts shm.
    env_base.setdefault("FI_PROVIDER", "tcp")
    env_base.setdefault("MPICH_CH4_NETMOD", "ofi")
    env_base.setdefault("UCX_TLS", "tcp,self")

    for n in ranks:
        out_path = os.path.join(
            out_dir,
            f"scale_r{n}_res{args.res}_rt{args.run_time:g}_{args.port_formulation}.json",
        )
        env = env_base.copy()
        env["BENCH_OUT"] = out_path
        cmd = [
            args.mpirun,
            "-np",
            str(n),
            args.python,
            worker_path,
        ]
        print(f"\n=== MPI ranks={n} OMP={args.omp_threads} ===", flush=True)
        print(" ".join(cmd), flush=True)
        t0 = time.time()
        proc = subprocess.run(cmd, env=env, cwd=ROOT)
        wall = time.time() - t0
        row = {"ranks": n, "wall_s": wall, "returncode": proc.returncode}
        if os.path.exists(out_path):
            with open(out_path, encoding="utf-8") as handle:
                row.update(json.load(handle))
        results.append(row)
        print(
            f"wall={wall:.1f}s device={row.get('t_device_s', 'NA')} "
            f"flux={row.get('flux', 'NA')}",
            flush=True,
        )

    # Speeds relative to the first successful rank count
    base_dev = None
    for row in results:
        if row.get("returncode", 1) == 0 and row.get("t_device_s"):
            base_dev = float(row["t_device_s"])
            base_ranks = int(row["ranks"])
            break

    enriched = []
    for row in results:
        t = row.get("t_device_s")
        entry = dict(row)
        if base_dev and t:
            speedup = base_dev / float(t)
            entry["speedup_vs_first"] = speedup
            entry["parallel_efficiency_vs_first"] = speedup / (
                int(row["ranks"]) / base_ranks
            )
        enriched.append(entry)

    summary = {
        "res": args.res,
        "run_time": args.run_time,
        "port": args.port,
        "formulation": args.port_formulation,
        "omp_num_threads": args.omp_threads,
        "results": enriched,
    }
    json_out = args.json_out or os.path.join(
        out_dir,
        f"mpi_scale_summary_res{args.res}_{args.port_formulation}.json",
    )
    with open(json_out, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print("\n## MPI scaling summary")
    for row in enriched:
        print(
            f"  ranks={row['ranks']}: device={row.get('t_device_s')}s "
            f"wall={row['wall_s']:.1f}s "
            f"speedup={row.get('speedup_vs_first', float('nan')):.2f} "
            f"eff={row.get('parallel_efficiency_vs_first', float('nan')):.2%} "
            f"flux={row.get('flux')}"
        )
    print(f"Wrote {json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
