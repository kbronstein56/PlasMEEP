#!/usr/bin/env python3
"""
Per-direction port power decomposition for P1↔P2 axis diagnosis.

Logs incident power, raw/subtracted flux, modal overlap transmission, and
Yee-line Hz samples for one source → one receiver excitation.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import meep as mp
import numpy as np

VAL = Path(__file__).resolve().parent
ROOT = VAL.parents[1]
for p in (str(VAL), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from port_formulations import (  # noqa: E402
    add_flux_monitor,
    add_modal_overlap_monitor_for_port,
    make_flux_region_for_formulation,
    make_numerical_mode_sources,
)
from plasmeep.ports.grid_index_port import build_grid_index_port_line  # noqa: E402
from plasmeep.ports.lorentz_probe import hz_yee_site  # noqa: E402
from plasmeep.ports.modal_receiver import (  # noqa: E402
    extract_modal_coefficient,
    sample_hz_line,
)
from sixport_common import (  # noqa: E402
    build_circulator_device,
    clear_width,
    default_uniform_rho,
    effective_port_dir,
    ensure_normalizations,
    fs_a,
    horn_for_port,
    monitor_center_for_port,
    normalize_port,
    set_geometry_context,
    source_df,
)


OUT = ROOT / "outputs" / "validation" / "port_axis"


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, complex):
        return {"real": obj.real, "imag": obj.imag}
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def _site_dict(xy: np.ndarray, *, res: int) -> Dict[str, Any]:
    s = hz_yee_site(xy, res=res)
    return s.as_dict()


def _read_hz_from_existing_dfts(
    sim: mp.Simulation,
    mon_info: Dict[str, Any],
) -> np.ndarray:
    return sample_hz_line(sim, mon_info)


def _flux_signed(
    sim: mp.Simulation,
    monitor,
    sign: float,
) -> float:
    return float(sign * mp.get_fluxes(monitor)[0])


def run_direction_diagnostic(
    *,
    source_port: int,
    receive_port: int,
    res: int,
    run_time: float,
    device_mode: str,
    formulation: str = "num_mode_guide_normal",
    grid_offset: Tuple[float, float] = (0.0, 0.0),
    verbose: bool = True,
) -> Dict[str, Any]:
    """One excitation: source_port in, measure receive_port (+ source reflection)."""
    set_geometry_context(res=res, grid_offset_cells=grid_offset)
    rho = default_uniform_rho()
    B = np.zeros(3)

    incident_power, incident_flux_data = normalize_port(
        source_port,
        res=res,
        run_time=run_time,
        verbose=verbose,
        formulation=formulation,
    )

    src_center = np.asarray(horn_for_port(source_port, res)["source_center"], dtype=float)
    src_mon = np.asarray(monitor_center_for_port(source_port, res), dtype=float)
    rcv_mon = np.asarray(monitor_center_for_port(receive_port, res), dtype=float)
    u_src = effective_port_dir(source_port)
    u_rcv = effective_port_dir(receive_port)
    span = 0.96 * clear_width

    def _one_device_run(*, subtract: bool) -> Dict[str, Any]:
        _pmm, p_dev, _ = build_circulator_device(
            rho,
            B,
            res=res,
            device_mode=device_mode,
        )
        p_dev.sources = make_numerical_mode_sources(source_port)
        sim = p_dev.Get_Sim()

        src_regions, src_sign = make_flux_region_for_formulation(
            formulation, src_mon, u_src, source_port
        )
        rcv_regions, rcv_sign = make_flux_region_for_formulation(
            formulation, rcv_mon, u_rcv, receive_port
        )
        mon_src = add_flux_monitor(sim, src_regions)
        mon_rcv = add_flux_monitor(sim, rcv_regions)
        modal_rcv = add_modal_overlap_monitor_for_port(sim, receive_port)
        modal_src = add_modal_overlap_monitor_for_port(sim, source_port)

        if subtract and incident_flux_data is not None:
            sim.load_minus_flux_data(mon_src, incident_flux_data)

        t0 = time.perf_counter()
        sim.run(until_after_sources=run_time)
        wall = time.perf_counter() - t0

        flux_src = _flux_signed(sim, mon_src, src_sign)
        flux_rcv = _flux_signed(sim, mon_rcv, rcv_sign)
        coeff_rcv = extract_modal_coefficient(sim, modal_rcv)
        coeff_src = extract_modal_coefficient(sim, modal_src)
        hz_rcv = _read_hz_from_existing_dfts(sim, modal_rcv)

        return {
            "wall_s": wall,
            "flux_source_signed": flux_src,
            "flux_receive_signed": flux_rcv,
            "flux_source_norm": flux_src / incident_power,
            "flux_receive_norm": flux_rcv / incident_power,
            "modal_coeff_receive": complex(coeff_rcv),
            "modal_coeff_source": complex(coeff_src),
            "modal_power_receive": float(np.abs(coeff_rcv) ** 2),
            "modal_power_source": float(np.abs(coeff_src) ** 2),
            "modal_power_receive_norm": float(np.abs(coeff_rcv) ** 2) / incident_power,
            "hz_rcv_rms": float(np.sqrt(np.mean(np.abs(hz_rcv) ** 2))),
        }

    raw = _one_device_run(subtract=False)
    net = _one_device_run(subtract=True)

    grid_line_src = build_grid_index_port_line(
        source_port,
        res=res,
        center_xy=src_mon,
        outward_dir=u_src,
        span_a=span,
    )
    grid_line_rcv = build_grid_index_port_line(
        receive_port,
        res=res,
        center_xy=rcv_mon,
        outward_dir=u_rcv,
        span_a=span,
    )

    return {
        "source_port": source_port,
        "receive_port": receive_port,
        "direction": f"P{source_port + 1}->P{receive_port + 1}",
        "res": res,
        "run_time": run_time,
        "device_mode": device_mode,
        "formulation": formulation,
        "grid_offset_cells": list(grid_offset),
        "incident_power": float(incident_power),
        "geometry": {
            "source_center": _site_dict(src_center, res=res),
            "source_monitor": _site_dict(src_mon, res=res),
            "receive_monitor": _site_dict(rcv_mon, res=res),
            "grid_line_source": grid_line_src.as_dict(),
            "grid_line_receive": grid_line_rcv.as_dict(),
        },
        "raw_device": raw,
        "after_subtraction": net,
        "subtraction_delta": {
            "flux_source": net["flux_source_signed"] - raw["flux_source_signed"],
            "flux_receive": net["flux_receive_signed"] - raw["flux_receive_signed"],
        },
        "reflection_positive_after_sub": bool(net["flux_source_norm"] > 0),
    }


def run_pair_diagnostic(
    *,
    res: int,
    run_time: float,
    device_mode: str,
    formulation: str,
    grid_offset: Tuple[float, float],
) -> Dict[str, Any]:
    """P1→P2 and P2→P1 with reciprocity metrics on flux and modal overlap."""
    d12 = run_direction_diagnostic(
        source_port=0,
        receive_port=1,
        res=res,
        run_time=run_time,
        device_mode=device_mode,
        formulation=formulation,
        grid_offset=grid_offset,
        verbose=False,
    )
    d21 = run_direction_diagnostic(
        source_port=1,
        receive_port=0,
        res=res,
        run_time=run_time,
        device_mode=device_mode,
        formulation=formulation,
        grid_offset=grid_offset,
        verbose=False,
    )

    t12 = d12["after_subtraction"]["flux_receive_norm"]
    t21 = d21["after_subtraction"]["flux_receive_norm"]
    m12 = d12["after_subtraction"]["modal_power_receive_norm"]
    m21 = d21["after_subtraction"]["modal_power_receive_norm"]

    def _recip_err(a: float, b: float) -> float:
        if a <= 0 or b <= 0:
            return float("nan")
        return abs(20 * np.log10(a / b))

    return {
        "res": res,
        "device_mode": device_mode,
        "formulation": formulation,
        "grid_offset_cells": list(grid_offset),
        "run_time": run_time,
        "P1_to_P2": d12,
        "P2_to_P1": d21,
        "reciprocity": {
            "flux_T12_norm": float(t12),
            "flux_T21_norm": float(t21),
            "flux_recip_err_dB": _recip_err(t12, t21),
            "modal_T12_norm": float(m12),
            "modal_T21_norm": float(m21),
            "modal_recip_err_dB": _recip_err(m12, m21),
            "P1_reflection_positive": d12["reflection_positive_after_sub"],
            "P2_reflection_positive": d21["reflection_positive_after_sub"],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, default=32)
    parser.add_argument("--run-time", type=float, default=5.0)
    parser.add_argument("--device-mode", type=str, default="horns_only")
    parser.add_argument("--formulation", type=str, default="num_mode_guide_normal")
    parser.add_argument("--grid-offset", type=str, default="0,0")
    parser.add_argument("--label", type=str, default="")
    parser.add_argument("--json-out", type=str, default="")
    args = parser.parse_args()

    parts = [p.strip() for p in args.grid_offset.split(",")]
    grid_offset = (float(parts[0]), float(parts[1]))

    result = run_pair_diagnostic(
        res=args.res,
        run_time=args.run_time,
        device_mode=args.device_mode,
        formulation=args.formulation,
        grid_offset=grid_offset,
    )
    result["timestamp_utc"] = datetime.now(timezone.utc).isoformat()
    result["label"] = args.label or (
        f"axis_diag_{args.device_mode}_{args.formulation}_res{args.res}_rt{args.run_time:g}"
    )

    OUT.mkdir(parents=True, exist_ok=True)
    json_out = Path(args.json_out) if args.json_out else OUT / f"{result['label']}.json"

    try:
        from mpi4py import MPI

        rank = int(MPI.COMM_WORLD.Get_rank())
    except Exception:
        rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", "0"))

    if rank == 0:
        json_out.write_text(json.dumps(_json_safe(result), indent=2) + "\n", encoding="utf-8")
        r = result["reciprocity"]
        print(f"Wrote {json_out}")
        print(
            f"flux reciprocity {r['flux_recip_err_dB']:.3f} dB  "
            f"modal {r['modal_recip_err_dB']:.3f} dB  "
            f"T12={r['flux_T12_norm']:.4e} T21={r['flux_T21_norm']:.4e}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
