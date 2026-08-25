#!/usr/bin/env python3
"""
MPI scaling benchmark: one device-port excitation of the six-port geometry.

Times a single source-port Meep run (after a one-port normalization) at fixed
res/run_time while varying MPI ranks. Does not change production physics.
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
    make_port_source,
    make_flux_region,
    port_monitor_centers,
    port_dirs,
    fs_a,
)

res = int(os.environ["BENCH_RES"])
run_time = float(os.environ["BENCH_RUNTIME"])
port = int(os.environ["BENCH_PORT"])

t0 = time.time()
cache = ensure_normalizations(res, run_time=run_time, force=True, ports=[port], verbose=False)
t_norm = time.time() - t0

rho = default_uniform_rho()
t1 = time.time()
_pmm, P_device, _wp = build_circulator_device(rho, [0, 0, 0], res=res)
P_device.sources = make_port_source(port)
sim = P_device.Get_Sim()
region, sign = make_flux_region(port_monitor_centers[port], port_dirs[port])
mon = sim.add_flux(fs_a, 0, 1, region)
sim.load_minus_flux_data(mon, cache["incident_flux_data_by_port"][port])
sim.run(until_after_sources=run_time)
flux = abs(sign * mp.get_fluxes(mon)[0])
t_dev = time.time() - t1

out = {
    "rank_env": int(os.environ.get("OMPI_COMM_WORLD_SIZE", os.environ.get("PMI_SIZE", "1"))),
    "t_norm_s": t_norm,
    "t_device_s": t_dev,
    "flux": float(flux),
    "res": res,
    "run_time": run_time,
    "port": port,
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
    parser.add_argument("--res", type=int, default=32)
    parser.add_argument("--run-time", type=float, default=40.0)
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--ranks", type=str, default="1,2,4,8")
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
    args = parser.parse_args()

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

    for n in ranks:
        out_path = os.path.join(out_dir, f"scale_r{n}_res{args.res}.json")
        env = env_base.copy()
        env["BENCH_OUT"] = out_path
        cmd = [
            args.mpirun,
            "-np",
            str(n),
            args.python,
            worker_path,
        ]
        print(f"\n=== MPI ranks={n} ===", flush=True)
        print(" ".join(cmd), flush=True)
        t0 = time.time()
        proc = subprocess.run(cmd, env=env, cwd=ROOT)
        wall = time.time() - t0
        row = {"ranks": n, "wall_s": wall, "returncode": proc.returncode}
        if os.path.exists(out_path):
            with open(out_path, encoding="utf-8") as handle:
                row.update(json.load(handle))
        results.append(row)
        print(f"wall={wall:.1f}s device={row.get('t_device_s', 'NA')}", flush=True)

    summary = {
        "res": args.res,
        "run_time": args.run_time,
        "port": args.port,
        "results": results,
    }
    json_out = args.json_out or os.path.join(out_dir, f"mpi_scale_summary_res{args.res}.json")
    with open(json_out, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    print("\n## MPI scaling summary")
    base = None
    for row in results:
        t = row.get("t_device_s", row.get("wall_s"))
        if base is None and t:
            base = t
        speedup = (base / t) if (base and t) else float("nan")
        print(
            f"  ranks={row['ranks']}: device={row.get('t_device_s')}s "
            f"wall={row['wall_s']:.1f}s speedup={speedup:.2f}"
        )
    print(f"Wrote {json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
