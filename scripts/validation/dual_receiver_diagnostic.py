#!/usr/bin/env python3
"""
Dual-receiver B=0 port diagnostic: flux/subtraction vs direct modal overlap.

In one FDTD run per source port, measure BOTH:
  1. guide-normal flux (optional load_minus_flux_data on source monitor)
  2. numerical-mode overlap on raw DFT Hz (NO flux subtraction)

This answers whether load_minus_flux_data / flux integration causes S21!=S12
without doubling wall time.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

import meep as mp
import numpy as np

from physical_units import add_resolution_arguments, resolve_simulation_resolution
from sixport_common import (
    a,
    build_circulator_device,
    default_uniform_rho,
    effective_port_dir,
    ensure_normalizations,
    fs_Hz,
    geometry_summary,
    monitor_center_for_port,
    reciprocity_metrics,
    set_geometry_context,
    fs_a,
    source_df,
)


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, complex):
        return {"real": float(obj.real), "imag": float(obj.imag), "abs": float(abs(obj))}
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if obj is None or isinstance(obj, (str, bool, int, float)):
        return obj
    return str(obj)


def _parse_ports(s: str) -> List[int]:
    ports = [int(p.strip()) for p in s.split(",") if p.strip() != ""]
    for p in ports:
        if p < 0 or p > 5:
            raise ValueError(f"port {p} out of range")
    return ports


def _mpi_info() -> Dict[str, Any]:
    try:
        from mpi_runner import mpi_info

        return mpi_info()
    except Exception:
        return {"ranks": 1, "rank": 0}


def _pair_err_dB(M_dB: np.ndarray, i: int, j: int) -> float:
    aij = M_dB[i, j]
    aji = M_dB[j, i]
    if not (np.isfinite(aij) and np.isfinite(aji)):
        return float("nan")
    return float(abs(aij - aji))


def run_dual(
    *,
    ports: List[int],
    res: int,
    run_time: float,
    device_mode: str,
    skip_flux_subtraction: bool,
    monitor_outward_cells: float,
) -> Dict[str, Any]:
    from port_formulations import (
        add_flux_monitor,
        add_modal_overlap_monitor_for_port,
        extract_flux_powers,
        extract_modal_coefficients,
        get_formulation,
        make_flux_region_for_formulation,
        make_numerical_mode_sources,
    )
    from plasmeep.ports.modal_receiver import modal_coefficient_metrics

    set_geometry_context(
        res=res,
        monitor_offset_cells=(monitor_outward_cells, 0.0),
    )
    form_flux = get_formulation("num_mode_guide_normal")

    # Incident for flux method (with flux data for optional subtraction).
    # Match reciprocity_b0_study cache naming so full-device ppc50 norms reuse.
    cache_dir = os.path.join(os.path.dirname(__file__), ".cache")
    os.makedirs(cache_dir, exist_ok=True)
    from physical_units import points_per_cm_from_meep_resolution

    ppc = points_per_cm_from_meep_resolution(res, a_m=a)
    off_tag = f"g0_0_m{monitor_outward_cells:g}_0"
    port_tag = "p" + "".join(str(p) for p in ports)
    flux_cache = os.path.join(
        cache_dir,
        f"norm_num_mode_guide_normal_{device_mode}_prism_rot0_pec_{off_tag}_ppc{ppc:g}_res{res}_rt{run_time:g}_{port_tag}.pkl",
    )
    flux_norm = ensure_normalizations(
        res=res,
        run_time=run_time,
        ports=ports,
        cache_path=flux_cache,
        force=False,
        verbose=True,
        formulation="num_mode_guide_normal",
    )

    # Incident for modal (no flux data needed)
    modal_cache = os.path.join(
        cache_dir,
        f"norm_num_mode_modal_{device_mode}_prism_rot0_pec_{off_tag}_ppc{ppc:g}_res{res}_rt{run_time:g}_{port_tag}.pkl",
    )
    modal_norm = ensure_normalizations(
        res=res,
        run_time=run_time,
        ports=ports,
        cache_path=modal_cache,
        force=False,
        verbose=True,
        formulation="num_mode_modal",
    )

    rho = default_uniform_rho()
    B = np.zeros(3)
    _pmm, P_device, _wp = build_circulator_device(
        rho, B, res=res, device_mode=device_mode
    )

    n = len(ports)
    local = {p: k for k, p in enumerate(ports)}
    flux_raw = np.full((n, n), np.nan)
    flux_normed = np.full((n, n), np.nan)
    modal_power = np.full((n, n), np.nan)
    modal_coeff: Dict[str, complex] = {}
    hz_profiles: Dict[str, Any] = {}

    for source_port in ports:
        print(f"\n=== dual receive, source P{source_port + 1} ===", flush=True)
        P_device.sources = make_numerical_mode_sources(source_port)
        sim = P_device.Get_Sim()

        flux_mons = []
        flux_signs = []
        modal_infos = []
        for out_port in ports:
            center = monitor_center_for_port(out_port, res)
            regions, sign = make_flux_region_for_formulation(
                form_flux.name, center, effective_port_dir(out_port), out_port
            )
            flux_mons.append(add_flux_monitor(sim, regions))
            flux_signs.append(sign)
            modal_infos.append(add_modal_overlap_monitor_for_port(sim, out_port))

        src_local = local[source_port]
        if not skip_flux_subtraction:
            sim.load_minus_flux_data(
                flux_mons[src_local],
                flux_norm["incident_flux_data_by_port"][source_port],
            )

        sim.run(until_after_sources=run_time)

        f_raw = extract_flux_powers(sim, flux_mons, flux_signs)
        coeffs = extract_modal_coefficients(sim, modal_infos)
        m_pow = np.array([float(np.real(c * np.conj(c))) for c in coeffs])

        from plasmeep.ports.modal_receiver import sample_hz_line
        from plasmeep.ports.mode_registry import get_numerical_mode

        for out_port, info, c in zip(ports, modal_infos, coeffs):
            hz = sample_hz_line(sim, info)
            mode = get_numerical_mode(out_port, res=res, frequency_a=fs_a)
            key = f"P{out_port + 1}_from_P{source_port + 1}"
            hz_profiles[key] = {
                "hz_abs": np.abs(hz).tolist(),
                "hz_phase_deg": np.degrees(np.angle(hz)).tolist(),
                "offsets_a": mode.offsets_a.tolist(),
                "mode_overlap_abs": float(abs(c)),
                "mode_power": float(abs(c) ** 2),
                "overlap_with_cached_mode": float(abs(mode.receiver_coefficient(hz))),
            }

        flux_raw[:, src_local] = f_raw
        flux_normed[:, src_local] = f_raw / flux_norm["incident_power_by_port"][source_port]
        modal_power[:, src_local] = m_pow / modal_norm["incident_power_by_port"][source_port]
        for out_port, c in zip(ports, coeffs):
            modal_coeff[f"P{out_port + 1}_from_P{source_port + 1}"] = complex(c)

    def to_dB(M: np.ndarray) -> np.ndarray:
        out = np.full_like(M, np.nan, dtype=float)
        pos = M > 0
        out[pos] = 10.0 * np.log10(M[pos])
        return out

    flux_dB = to_dB(flux_normed)
    modal_dB = to_dB(modal_power)

    # Complex modal reciprocity for first pair if 2-port
    modal_complex_metrics = None
    if n == 2:
        pa, pb = ports
        cab = modal_coeff[f"P{pb + 1}_from_P{pa + 1}"]
        cba = modal_coeff[f"P{pa + 1}_from_P{pb + 1}"]
        modal_complex_metrics = modal_coefficient_metrics(cab, cba)

    # Discrete control
    discrete = None
    if n == 2:
        from plasmeep.ports.lorentz_probe import ProbeSites, evaluate_reciprocity_pair
        from sixport_common import horn_for_port

        pa, pb = ports
        sites_a = ProbeSites(
            pa,
            np.asarray(horn_for_port(pa, res)["source_center"], dtype=float),
            np.asarray(monitor_center_for_port(pa, res), dtype=float),
            np.asarray(effective_port_dir(pa), dtype=float),
        )
        sites_b = ProbeSites(
            pb,
            np.asarray(horn_for_port(pb, res)["source_center"], dtype=float),
            np.asarray(monitor_center_for_port(pb, res), dtype=float),
            np.asarray(effective_port_dir(pb), dtype=float),
        )
        disc = evaluate_reciprocity_pair(
            P_device,
            sites_a,
            sites_b,
            res=res,
            frequency=fs_a,
            fwidth=source_df,
            run_time=run_time,
            observable="discrete",
        )
        discrete = disc.get("discrete", {})

    # Flux subtraction health on normalized flux matrix
    diag_pos = 0
    diag_vals = {}
    for k, p in enumerate(ports):
        v = float(flux_normed[k, k])
        diag_vals[f"P{p + 1}"] = v
        if np.isfinite(v) and v > 0:
            diag_pos += 1

    pair_summaries = []
    for i in range(n):
        for j in range(i + 1, n):
            pair_summaries.append(
                {
                    "ports": [ports[i] + 1, ports[j] + 1],
                    "flux_diff_dB": _pair_err_dB(flux_dB, i, j),
                    "modal_power_diff_dB": _pair_err_dB(modal_dB, i, j),
                    "flux_Tij_dB": float(flux_dB[j, i]) if np.isfinite(flux_dB[j, i]) else None,
                    "flux_Tji_dB": float(flux_dB[i, j]) if np.isfinite(flux_dB[i, j]) else None,
                    "modal_Tij_dB": float(modal_dB[j, i]) if np.isfinite(modal_dB[j, i]) else None,
                    "modal_Tji_dB": float(modal_dB[i, j]) if np.isfinite(modal_dB[i, j]) else None,
                }
            )

    return {
        "ports": ports,
        "skip_flux_subtraction": skip_flux_subtraction,
        "monitor_outward_cells": monitor_outward_cells,
        "flux_incident": {str(k): float(v) for k, v in flux_norm["incident_power_by_port"].items()},
        "modal_incident": {str(k): float(v) for k, v in modal_norm["incident_power_by_port"].items()},
        "flux_raw_matrix": flux_raw,
        "flux_normalized": flux_normed,
        "flux_normalized_dB": flux_dB,
        "modal_power_normalized": modal_power,
        "modal_power_normalized_dB": modal_dB,
        "modal_coefficients": modal_coeff,
        "hz_profiles": hz_profiles,
        "modal_complex_metrics": modal_complex_metrics,
        "flux_subtraction_health": {
            "diagonal_normalized_power": diag_vals,
            "positive_reflection_count": diag_pos,
            "healthy": diag_pos == 0,
            "note": "load_minus_flux_data only affects source-port diagonal; transmission off-diagonals are always raw flux/incident",
        },
        "pair_summaries": pair_summaries,
        "discrete_control": discrete,
        "flux_reciprocity": reciprocity_metrics(flux_dB, port_ids=ports),
        "modal_reciprocity": reciprocity_metrics(modal_dB, port_ids=ports),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_resolution_arguments(parser)
    parser.add_argument("--run-time", type=float, default=20.0)
    parser.add_argument("--ports", type=str, default="0,1")
    parser.add_argument("--device-mode", choices=["full", "horns_only"], default="horns_only")
    parser.add_argument("--skip-flux-subtraction", action="store_true")
    parser.add_argument(
        "--monitor-outward-cells",
        type=float,
        default=0.0,
        help="Shift monitors farther into the feed (+outward, in cells).",
    )
    parser.add_argument("--label", type=str, default="dual")
    parser.add_argument("--json-out", type=str, default="")
    parser.add_argument("--discrete-control", action="store_true")
    args = parser.parse_args()

    sim_res, ppc, res_report = resolve_simulation_resolution(
        args, a_m=a, default_points_per_cm=50.0
    )
    ports = _parse_ports(args.ports)
    print(geometry_summary())
    print(f"ppc={ppc} res={sim_res} device={args.device_mode} ports={ports}")
    print(f"skip_flux_subtraction={args.skip_flux_subtraction}")
    print(f"monitor_outward_cells={args.monitor_outward_cells}")

    t0 = time.perf_counter()
    result = run_dual(
        ports=ports,
        res=sim_res,
        run_time=args.run_time,
        device_mode=args.device_mode,
        skip_flux_subtraction=args.skip_flux_subtraction,
        monitor_outward_cells=args.monitor_outward_cells,
    )
    # discrete already inside if 2-port; flag kept for CLI compatibility
    _ = args.discrete_control
    elapsed = time.perf_counter() - t0

    payload = {
        "label": args.label,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "mpi": _mpi_info(),
        "settings": {
            "points_per_cm": ppc,
            "res": sim_res,
            "run_time": args.run_time,
            "device_mode": args.device_mode,
            "ports": ports,
            "a": a,
            "fs_Hz": fs_Hz,
            "skip_flux_subtraction": args.skip_flux_subtraction,
            "monitor_outward_cells": args.monitor_outward_cells,
            "source": "num_mode",
            "receivers": ["guide_normal_flux", "numerical_mode_overlap"],
        },
        "physical_resolution": res_report,
        "timings_s": {"total": elapsed},
        "results": result,
    }

    try:
        from mpi4py import MPI

        rank = int(MPI.COMM_WORLD.Get_rank())
    except Exception:
        rank = 0

    json_out = args.json_out
    if not json_out:
        json_out = f"dual_receiver_{args.label}.json"

    if rank == 0:
        with open(json_out, "w", encoding="utf-8") as f:
            json.dump(_json_safe(payload), f, indent=2)
            f.write("\n")
        print(f"\nWrote {json_out}")
        for row in result["pair_summaries"]:
            print(
                f"  P{row['ports'][0]}<->P{row['ports'][1]}: "
                f"flux={row['flux_diff_dB']:.4g} dB  "
                f"modal={row['modal_power_diff_dB']:.4g} dB"
            )
        if result.get("discrete_control"):
            print(f"  discrete={result['discrete_control'].get('amp_err_dB')} dB")


if __name__ == "__main__":
    main()
