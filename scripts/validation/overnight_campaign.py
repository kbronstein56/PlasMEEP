#!/usr/bin/env python3
"""
Overnight validation campaign: timing, MPI, run-time convergence, reciprocity,
mode profiles, caching policy. Saves incremental JSON after each phase.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
VAL = ROOT / "scripts" / "validation"
OUT = ROOT / "outputs" / "validation" / "overnight"
PYTHON = "/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/python"
MPIRUN = "/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/mpirun"


def _save(name: str, data: Dict[str, Any]) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
        f.write("\n")
    print(f"[saved] {path}", flush=True)
    return path


def _load_or_empty(name: str) -> Dict[str, Any]:
    path = OUT / name
    if path.is_file():
        return json.loads(path.read_text(encoding="utf-8"))
    return {}


def _mpi_env() -> dict:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env.setdefault("FI_PROVIDER", "tcp")
    env.setdefault("MPICH_CH4_NETMOD", "ofi")
    env["PYTHONPATH"] = f"{VAL}:{ROOT / 'scripts'}:" + env.get("PYTHONPATH", "")
    return env


def phase_timing(res: int, run_time: float) -> Dict[str, Any]:
    existing = _load_or_empty(f"01_timing_breakdown_res{res}.json")
    if existing.get("breakdown"):
        print("Phase 1: timing already done, skipping", flush=True)
        return existing

    cmd = [
        PYTHON,
        str(VAL / "profile_horns_only.py"),
        "--res",
        str(res),
        "--run-time",
        str(run_time),
        "--json-out",
        str(OUT / f"01_timing_breakdown_res{res}.json"),
    ]
    subprocess.run(cmd, cwd=ROOT, env=_mpi_env(), check=True)
    return json.loads((OUT / f"01_timing_breakdown_res{res}.json").read_text())


def phase_mpi_cheap(res: int, run_time: float, ranks: List[int]) -> Dict[str, Any]:
    worker = r'''
import json, os, sys, time
sys.path.insert(0, os.environ["VAL"])
sys.path.insert(0, os.environ["SCRIPTS"])
from plasmeep.ports.lorentz_probe import run_hz_transfer
from sixport_common import (
    build_circulator_device, default_uniform_rho, fs_a, source_df,
    horn_for_port, monitor_center_for_port, set_geometry_context,
)

res = int(os.environ["BENCH_RES"])
run_time = float(os.environ["BENCH_RT"])
port_a = int(os.environ["PORT_A"])
port_b = int(os.environ["PORT_B"])

set_geometry_context(res=res, horn_walls="prism")
rho = default_uniform_rho()
B = [0.0, 0.0, 0.0]
t0 = time.perf_counter()
_pmm, dev, _ = build_circulator_device(rho, B, res=res, device_mode="horns_only")
t_build = time.perf_counter() - t0

src = horn_for_port(port_a, res)["source_center"]
mon = monitor_center_for_port(port_b, res)
t1 = time.perf_counter()
hz = run_hz_transfer(
    dev, res=res, source_xy=src, monitor_xy=mon,
    frequency=fs_a, fwidth=source_df, run_time=run_time,
)
t_xfer = time.perf_counter() - t1

out = {
    "ranks": int(os.environ.get("OMPI_COMM_WORLD_SIZE", os.environ.get("PMI_SIZE", "1"))),
    "t_build_s": t_build,
    "t_transfer_s": t_xfer,
    "t_total_s": t_build + t_xfer,
    "hz_abs": abs(hz),
}
rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", os.environ.get("PMI_RANK", "0")))
if rank == 0:
    with open(os.environ["OUT"], "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
        f.write("\n")
'''
    worker_path = OUT / "_mpi_lorentz_worker.py"
    worker_path.write_text(worker, encoding="utf-8")

    results = []
    for n in ranks:
        out_json = OUT / f"02_mpi_r{n}_res{res}.json"
        if out_json.is_file():
            row = json.loads(out_json.read_text())
            results.append(row)
            continue
        env = _mpi_env()
        env["VAL"] = str(VAL)
        env["SCRIPTS"] = str(ROOT / "scripts")
        env["BENCH_RES"] = str(res)
        env["BENCH_RT"] = str(run_time)
        env["PORT_A"] = "0"
        env["PORT_B"] = "1"
        env["OUT"] = str(out_json)
        cmd = [MPIRUN, "-np", str(n), PYTHON, str(worker_path)]
        print(f"Phase 2: MPI np={n}", flush=True)
        t0 = time.perf_counter()
        proc = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True)
        wall = time.perf_counter() - t0
        row = {"ranks": n, "wall_s": wall, "returncode": proc.returncode}
        if out_json.is_file():
            try:
                row.update(json.loads(out_json.read_text()))
            except json.JSONDecodeError:
                out_json.unlink(missing_ok=True)
                raise
        results.append(row)
        _save(f"02_mpi_partial.json", {"results": results})

    best = min(results, key=lambda r: r.get("wall_s", 1e99))
    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "case": "horns_only lorentz single transfer P1->P2",
        "res": res,
        "run_time": run_time,
        "results": results,
        "optimal_ranks": best["ranks"],
        "optimal_wall_s": best.get("wall_s"),
        "note": "Use optimal_ranks for cheap res32 diagnostics; keep np=32 for res96.",
    }
    _save(f"02_mpi_summary_res{res}.json", summary)
    return summary


def phase_runtime_convergence(res: int) -> Dict[str, Any]:
    existing = _load_or_empty(f"03_runtime_convergence_res{res}.json")
    if existing.get("sweep"):
        print("Phase 3: runtime convergence already done", flush=True)
        return existing

    sys.path.insert(0, str(VAL))
    sys.path.insert(0, str(ROOT / "scripts"))
    from plasmeep.ports.lorentz_probe import (  # noqa: E402
        ProbeSites,
        evaluate_lorentz_pair,
    )
    from sixport_common import (  # noqa: E402
        build_circulator_device,
        default_uniform_rho,
        effective_port_dir,
        fs_a,
        horn_for_port,
        monitor_center_for_port,
        set_geometry_context,
        source_df,
    )

    set_geometry_context(res=res, horn_walls="prism")
    rho = default_uniform_rho()
    B = np.zeros(3)
    _pmm, dev, _ = build_circulator_device(rho, B, res=res, device_mode="horns_only")

    def sites(pa: int, pb: int) -> tuple[ProbeSites, ProbeSites]:
        return (
            ProbeSites(
                pa,
                np.asarray(horn_for_port(pa, res)["source_center"], float),
                np.asarray(monitor_center_for_port(pa, res), float),
                np.asarray(effective_port_dir(pa), float),
            ),
            ProbeSites(
                pb,
                np.asarray(horn_for_port(pb, res)["source_center"], float),
                np.asarray(monitor_center_for_port(pb, res), float),
                np.asarray(effective_port_dir(pb), float),
            ),
        )

    sa, sb = sites(0, 1)
    ref_rt = 80.0
    ref = evaluate_lorentz_pair(
        dev, sa, sb, res=res, frequency=fs_a, fwidth=source_df, run_time=ref_rt
    )

    run_times = [5, 10, 15, 20, 25, 30, 40, 60, 80]
    sweep = []
    gate_db = 0.01  # much tighter than 0.2 dB reciprocity gate

    for rt in run_times:
        pair = evaluate_lorentz_pair(
            dev, sa, sb, res=res, frequency=fs_a, fwidth=source_df, run_time=rt
        )
        h_err = abs(pair.a_to_b.hz_complex - ref.a_to_b.hz_complex) / max(
            abs(ref.a_to_b.hz_complex), 1e-30
        )
        sweep.append(
            {
                "run_time": rt,
                "amp_err_dB_reciprocity": pair.amp_err_dB,
                "phase_diff_deg": pair.phase_diff_deg,
                "h01_rel_err_vs_rt80": float(h_err),
            }
        )
        partial = {"reference_run_time": ref_rt, "sweep": sweep}
        _save(f"03_runtime_convergence_res{res}.json", partial)

    # find shortest rt where reciprocity amp err < gate and H stable vs ref
    safe = None
    for row in sweep:
        if (
            row["amp_err_dB_reciprocity"] < gate_db
            and row["h01_rel_err_vs_rt80"] < 0.02
        ):
            safe = row["run_time"]
            break

    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "res": res,
        "reference_run_time": ref_rt,
        "reciprocity_gate_dB": 0.2,
        "convergence_gate_dB": gate_db,
        "sweep": sweep,
        "recommended_run_time": safe if safe is not None else 40.0,
        "note": "recommended_run_time: shortest rt with |P12-P21|<0.01 dB and H01 within 2% of rt=80",
    }
    _save(f"03_runtime_convergence_res{res}.json", summary)
    return summary


def phase_direct_reciprocity(res: int, run_time: float, mpi_ranks: int) -> Dict[str, Any]:
    out_json = OUT / f"04_lorentz_direct_res{res}_rt{run_time:g}.json"
    if out_json.is_file():
        return json.loads(out_json.read_text())

    cmd = [
        MPIRUN,
        "-np",
        str(mpi_ranks),
        PYTHON,
        str(VAL / "direct_reciprocity_test.py"),
        "--res",
        str(res),
        "--run-time",
        str(run_time),
        "--pairs",
        "0,1;1,2",
    ]
    env = _mpi_env()
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)
    src = ROOT / "outputs" / "validation" / "lorentz_direct" / f"lorentz_direct_res{res}.json"
    data = json.loads(src.read_text())
    data["overnight_mpi_ranks"] = mpi_ranks
    data["overnight_run_time"] = run_time
    _save(out_json.name, data)
    return data


def phase_mode_profiles(res: int, run_time: float, mpi_ranks: int) -> Dict[str, Any]:
    out_json = ROOT / "outputs" / "validation" / "mode_profiles" / f"mode_profiles_res{res}.json"
    if out_json.is_file():
        return json.loads(out_json.read_text())

    cmd = [
        MPIRUN,
        "-np",
        str(mpi_ranks),
        PYTHON,
        str(VAL / "mode_profile_study.py"),
        "--res",
        str(res),
        "--run-time",
        str(run_time),
    ]
    subprocess.run(cmd, cwd=ROOT, env=_mpi_env(), check=True)
    return json.loads(out_json.read_text())


def write_cache_policy() -> Path:
    text = """# Validation cache policy

## Safe to cache (invalidate on listed key changes)

| Artifact | Cache key fields | Reuse scope |
|---|---|---|
| Source normalization pickle | res, run_time, formulation, grid_offset, monitor_offset, horn_walls, coord_rotation, ports | Same geometry + source type |
| Numerical port mode JSON | res, frequency, horn_walls, grid_offset, port symmetry class | Launch/measure for matching horns |
| Lorentz direct reciprocity JSON | res, run_time, horn_walls, grid_offset | Diagnostic replay / plotting |
| Mode profile JSON | res, run_time, horn_walls, ports | Overlap analysis without re-sim |
| MPI timing summary | res, case type | Rank selection only |

## Must recalculate when changing

- `res` (grid spacing)
- `fs_a` / source frequency or `source_df`
- `horn_walls`, `grid_offset_cells`, `coord_rotation_deg`
- `device_mode` (horns_only vs full)
- plasma `rho`, `B`, or bulb geometry (full device)
- `run_time` below convergence-safe minimum

## Symmetry reuse

- Ports P2/P3 and P5/P6 (+60°/−60° pairs): numerical mode profiles may share one reference per |angle|.
- P2↔P3 reciprocity control: run once; do not duplicate for symmetric copies.

## Do not cache

- Full field HDF5 arrays (unless explicitly needed)
- Meep `Simulation` objects across processes
- Normalizations computed with `force=True` mid-sweep without updating cache key

## Overnight cheap defaults (after profiling)

- horns_only res32: use MPI ranks from `02_mpi_summary_res32.json`
- run_time from `03_runtime_convergence_res32.json` `recommended_run_time`
- Skip completed JSON phases in `overnight_campaign.py`
"""
    path = OUT / "CACHE_POLICY.md"
    path.write_text(text, encoding="utf-8")
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, default=32)
    parser.add_argument("--run-time", type=float, default=40.0)
    parser.add_argument("--phases", type=str, default="all")
    args = parser.parse_args()

    phases = (
        ["timing", "mpi", "runtime", "lorentz", "modes", "policy"]
        if args.phases == "all"
        else args.phases.split(",")
    )

    report: Dict[str, Any] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "phases": {},
    }

    mpi_ranks = 4
    safe_rt = args.run_time

    if "timing" in phases:
        report["phases"]["timing"] = phase_timing(args.res, args.run_time)

    if "mpi" in phases:
        mpi_summary = phase_mpi_cheap(args.res, args.run_time, [4, 8, 16, 32])
        report["phases"]["mpi"] = mpi_summary
        mpi_ranks = int(mpi_summary["optimal_ranks"])

    if "runtime" in phases:
        rt_summary = phase_runtime_convergence(args.res)
        report["phases"]["runtime"] = rt_summary
        safe_rt = float(rt_summary["recommended_run_time"])

    if "lorentz" in phases:
        report["phases"]["lorentz"] = phase_direct_reciprocity(
            args.res, safe_rt, mpi_ranks
        )

    if "modes" in phases:
        report["phases"]["modes"] = phase_mode_profiles(args.res, safe_rt, mpi_ranks)

    if "policy" in phases:
        report["phases"]["cache_policy"] = str(write_cache_policy())

    _save("overnight_summary.json", report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
