#!/usr/bin/env python3
"""
Wall-clock breakdown for one representative horns_only reciprocity case.

Separates Python geometry, Meep init/rasterization, normalization, FDTD,
flux extraction, and I/O.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from typing import Any, Dict

import meep as mp
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VAL = os.path.dirname(os.path.abspath(__file__))
for p in (VAL, os.path.join(ROOT, "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from port_formulations import (  # noqa: E402
    add_flux_monitor,
    get_formulation,
    make_flux_region_for_formulation,
    port_measure_center,
)
from sixport_common import (  # noqa: E402
    build_circulator_device,
    default_uniform_rho,
    effective_port_dir,
    ensure_normalizations,
    physical_resolution_report,
    reciprocity_metrics,
    set_geometry_context,
    simulate_circulator,
)

OUT = os.path.join(ROOT, "outputs", "validation", "overnight")


def _pct(part: float, total: float) -> float:
    return 100.0 * part / total if total > 0 else 0.0


def profile_reciprocity_pair(
    *,
    res: int,
    run_time: float,
    formulation: str,
    ports: tuple[int, int],
) -> Dict[str, Any]:
    set_geometry_context(res=res, horn_walls="prism")
    form = get_formulation(formulation)
    port_list = list(ports)

    timings: Dict[str, float] = {}
    t_all = time.perf_counter()

    # --- geometry construction (Python) ---
    t0 = time.perf_counter()
    rho = default_uniform_rho()
    B = np.zeros(3)
    _pmm, p_device, _ = build_circulator_device(
        rho, B, res=res, device_mode="horns_only"
    )
    timings["python_geometry_build_s"] = time.perf_counter() - t0

    # --- source normalization (reference straight guides) ---
    t0 = time.perf_counter()
    norm = ensure_normalizations(
        res,
        run_time=run_time,
        force=True,
        ports=port_list,
        verbose=False,
        formulation=formulation,
    )
    timings["source_normalization_s"] = time.perf_counter() - t0

    # --- one device excitation with init vs FDTD split ---
    source_port = port_list[0]
    p_device.sources = form.make_sources(source_port)
    t0 = time.perf_counter()
    sim = p_device.Get_Sim()
    timings["meep_sim_construct_s"] = time.perf_counter() - t0

    measure_xy = port_measure_center(form.name, port_list[1])
    regions, sign = make_flux_region_for_formulation(
        form.name, measure_xy, effective_port_dir(port_list[1])
    )
    monitor = add_flux_monitor(sim, regions)

    t0 = time.perf_counter()
    sim.run(until_after_sources=run_time)
    timings["fdtd_timestepping_s"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    raw_flux = mp.get_fluxes(monitor)[0]
    signed = sign * raw_flux
    power = abs(signed) / norm["incident_power_by_port"][source_port]
    timings["flux_extraction_s"] = time.perf_counter() - t0
    timings["_sample_transmitted_power"] = float(power)

    # --- full reciprocity pair via harness (norm already done) ---
    t0 = time.perf_counter()
    result = simulate_circulator(
        rho,
        B,
        res=res,
        run_time=run_time,
        verbose=False,
        incident_cache=norm,
        ports=port_list,
        formulation=formulation,
        device_mode="horns_only",
    )
    timings["full_pair_simulate_circulator_s"] = time.perf_counter() - t0

    timings["total_wall_s"] = time.perf_counter() - t_all

    # Meep init is embedded in each Get_Sim+first run inside simulate_circulator;
    # approximate FDTD-only from single-run measurement scaled to pair count.
    n_excitations = len(port_list)
    timings["estimated_pair_fdtd_s"] = timings["fdtd_timestepping_s"] * n_excitations
    timings["estimated_pair_total_s"] = (
        timings["python_geometry_build_s"]
        + timings["source_normalization_s"]
        + timings["meep_sim_construct_s"] * n_excitations
        + timings["fdtd_timestepping_s"] * n_excitations
        + timings["flux_extraction_s"] * n_excitations
    )

    metrics = reciprocity_metrics(result["power_matrix_dB"], port_ids=result["ports"])

    breakdown = {
        "python_geometry_build": {
            "seconds": timings["python_geometry_build_s"],
            "percent": _pct(
                timings["python_geometry_build_s"], timings["total_wall_s"]
            ),
        },
        "source_normalization": {
            "seconds": timings["source_normalization_s"],
            "percent": _pct(timings["source_normalization_s"], timings["total_wall_s"]),
        },
        "meep_sim_construct_one_excitation": {
            "seconds": timings["meep_sim_construct_s"],
            "note": "Get_Sim() for one excitation; repeated per port run",
        },
        "fdtd_timestepping_one_excitation": {
            "seconds": timings["fdtd_timestepping_s"],
            "percent_of_single_run": 100.0
            * timings["fdtd_timestepping_s"]
            / max(
                timings["meep_sim_construct_s"] + timings["fdtd_timestepping_s"],
                1e-30,
            ),
        },
        "flux_extraction_one_excitation": {
            "seconds": timings["flux_extraction_s"],
        },
        "full_pair_simulate_circulator": {
            "seconds": timings["full_pair_simulate_circulator_s"],
            "percent": _pct(
                timings["full_pair_simulate_circulator_s"], timings["total_wall_s"]
            ),
        },
        "total_wall": {"seconds": timings["total_wall_s"]},
    }

    conclusion = (
        "Python horn geometry is a negligible fraction of wall time; "
        "normalization + FDTD dominate."
        if timings["python_geometry_build_s"] < 0.05 * timings["total_wall_s"]
        else "Python geometry build is non-negligible; investigate further."
    )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "case": {
            "device_mode": "horns_only",
            "res": res,
            "run_time": run_time,
            "formulation": formulation,
            "ports": list(ports),
        },
        "physical_resolution": physical_resolution_report(res),
        "reciprocity_max_abs_diff_dB": metrics["max_abs_diff_dB"],
        "timings_s": timings,
        "breakdown": breakdown,
        "conclusion": conclusion,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, default=32)
    parser.add_argument("--run-time", type=float, default=40.0)
    parser.add_argument("--formulation", type=str, default="te1_hz_line")
    parser.add_argument("--ports", type=str, default="0,1")
    parser.add_argument("--json-out", type=str, default="")
    args = parser.parse_args()

    ports = tuple(int(x) for x in args.ports.split(","))
    os.makedirs(OUT, exist_ok=True)
    payload = profile_reciprocity_pair(
        res=args.res,
        run_time=args.run_time,
        formulation=args.formulation,
        ports=ports,  # type: ignore[arg-type]
    )

    out = args.json_out or os.path.join(
        OUT, f"timing_breakdown_res{args.res}_rt{args.run_time:g}.json"
    )
    with open(out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
        f.write("\n")

    print(json.dumps(payload["breakdown"], indent=2))
    print(f"\nConclusion: {payload['conclusion']}")
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
