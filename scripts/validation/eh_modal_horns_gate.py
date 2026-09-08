#!/usr/bin/env python3
"""
Horns-only gate for a full electromagnetic E/H modal receiver at 50 ppc.

This script:
1. extracts full Ex/Ey/Hz reference modes for selected ports from straight-guide runs,
2. normalizes each reference mode to unit modal power using a Meep-style bilinear form,
3. runs reciprocal horns-only cases, and
4. reports forward/backward modal amplitudes and S21/S12 from |alpha|^2.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
SCRIPTS_DIR = str(Path(__file__).resolve().parents[1])
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from plasmeep.ports.em_mode import (
    EMLineMode,
    add_em_line_monitor,
    build_em_mode_from_samples,
    sample_em_line_fields,
)
from plasmeep.ports.mode_registry import get_numerical_mode

import meep as mp
from mpi4py import MPI

from PMMCirculatorInverse import PMMI
from mpi_runner import mpi_info
from physical_units import add_resolution_arguments, resolve_simulation_resolution
from port_formulations import make_yee_snapped_numerical_mode_sources
from sixport_common import (
    a,
    build_circulator_device,
    clear_width,
    default_uniform_rho,
    effective_port_dir,
    feed_center_for_port,
    fs_a,
    get_coord_rotation_deg,
    get_grid_offset_cells,
    get_horn_walls,
    horns_for_device,
    nx_ports,
    ny_ports,
    dpml_ports,
    rotated_wall,
    set_geometry_context,
    wall_thickness,
)

COMM = MPI.COMM_WORLD
RANK = int(COMM.Get_rank())


def _port_tangent(outward: np.ndarray) -> np.ndarray:
    u = np.asarray(outward, dtype=float)
    u = u / max(np.linalg.norm(u), 1e-30)
    return np.array([-u[1], u[0]], dtype=float)


def _complex_dict(z: complex) -> Dict[str, float]:
    return {"real": float(np.real(z)), "imag": float(np.imag(z))}


def _ensure_finite(name: str, arr: np.ndarray) -> None:
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name}: non-finite values detected")
    if np.max(np.abs(arr)) > 1e12:
        raise ValueError(f"{name}: absurdly large magnitude detected")


def _is_root() -> bool:
    return RANK == 0


def _barrier() -> None:
    try:
        COMM.Barrier()
    except Exception:
        pass


def _straight_guide_sim_for_port(port_index: int, *, res: int) -> mp.Simulation:
    u = np.asarray(effective_port_dir(port_index), dtype=float)
    u = u / max(np.linalg.norm(u), 1e-30)
    n = np.array([-u[1], u[0]], dtype=float)

    pmm_ref = PMMI(
        a=a,
        res=res,
        nx=nx_ports,
        ny=ny_ports,
        dpml=dpml_ports,
        B=np.array([0.0, 0.0, 0.0]),
    )
    pref = pmm_ref.Build_Sim()

    reference_center = np.asarray(
        horns_for_device(res)[port_index]["source_center"], dtype=float
    )
    guide_length = 2.0 * np.hypot(nx_ports, ny_ports)
    wall_center_offset = clear_width / 2 + wall_thickness / 2
    wall_a = rotated_wall(
        centerline=reference_center,
        axis=u,
        normal=n,
        normal_offset=+wall_center_offset,
        length=guide_length,
        thickness=wall_thickness,
    )
    wall_b = rotated_wall(
        centerline=reference_center,
        axis=u,
        normal=n,
        normal_offset=-wall_center_offset,
        length=guide_length,
        thickness=wall_thickness,
    )
    pref.Add_Prism(vertices=wall_a, axis=np.array([0, 0, 1]), PEC=True)
    pref.Add_Prism(vertices=wall_b, axis=np.array([0, 0, 1]), PEC=True)
    pref.sources = make_yee_snapped_numerical_mode_sources(port_index)
    return pref.Get_Sim()


def _reference_mode_path(res: int, port_index: int) -> Path:
    key = (
        f"res{res}_f{fs_a:.6f}_{get_horn_walls()}"
        f"_g{get_grid_offset_cells()[0]:g}_{get_grid_offset_cells()[1]:g}"
        f"_rot{get_coord_rotation_deg():g}"
    )
    return (
        ROOT
        / "outputs"
        / "validation"
        / "mode_profiles"
        / "em_modes"
        / f"res{res}"
        / key
        / f"eh_mode_P{port_index + 1}.json"
    )


def _reference_meta_path(res: int, port_index: int) -> Path:
    return _reference_mode_path(res, port_index).with_name(
        f"eh_mode_P{port_index + 1}_meta.json"
    )


def extract_reference_mode(
    port_index: int,
    *,
    res: int,
    run_time: float,
    force: bool,
) -> Dict[str, Any]:
    path = _reference_mode_path(res, port_index)
    meta_path = _reference_meta_path(res, port_index)
    if path.is_file() and meta_path.is_file() and not force:
        try:
            mode = EMLineMode.load(path)
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            return {
                "mode": mode,
                "cache_path": str(path),
                "cached": True,
                "raw_overlap": meta["raw_overlap"],
            }
        except Exception:
            pass

    num_mode = get_numerical_mode(
        port_index,
        res=res,
        frequency_a=fs_a,
        horn_walls=get_horn_walls(),
        grid_offset_cells=get_grid_offset_cells(),
        coord_rotation_deg=get_coord_rotation_deg(),
    )
    outward = np.asarray(effective_port_dir(port_index), dtype=float)
    tangent = _port_tangent(outward)
    center = feed_center_for_port(port_index, res, 0.50)

    sim = _straight_guide_sim_for_port(port_index, res=res)
    mon = add_em_line_monitor(
        sim,
        center_xy=center,
        tangent_xy=tangent,
        offsets_a=num_mode.offsets_a,
        frequency=fs_a,
    )
    sim.run(until_after_sources=run_time)
    fields = sample_em_line_fields(sim, mon)
    _ensure_finite(f"P{port_index+1} ref Ex", fields["ex"])
    _ensure_finite(f"P{port_index+1} ref Ey", fields["ey"])
    _ensure_finite(f"P{port_index+1} ref Hz", fields["hz"])

    raw_mode = build_em_mode_from_samples(
        port_index=port_index,
        label=f"P{port_index + 1} raw reference at clean plane",
        center_xy=center,
        tangent_xy=tangent,
        outward_dir=outward,
        offsets_a=num_mode.offsets_a,
        ex=fields["ex"],
        ey=fields["ey"],
        hz=fields["hz"],
    )

    # The reference guide run observes the mode traveling inward toward the device.
    # Flip Hz so the stored canonical mode points outward from the device.
    if raw_mode.mode_power < 0:
        raw_mode = build_em_mode_from_samples(
            port_index=port_index,
            label=f"P{port_index + 1} outward reference at clean plane",
            center_xy=center,
            tangent_xy=tangent,
            outward_dir=outward,
            offsets_a=num_mode.offsets_a,
            ex=fields["ex"],
            ey=fields["ey"],
            hz=-fields["hz"],
        )

    mode = raw_mode.normalized_to_unit_power()
    mode.assert_finite()
    if _is_root():
        path.parent.mkdir(parents=True, exist_ok=True)
        mode.save(path)

    raw_overlap = mode.overlap_amplitudes(
        ex_field=fields["ex"],
        ey_field=fields["ey"],
        hz_field=fields["hz"],
    )
    raw_overlap_json = {
        "alpha_plus": _complex_dict(complex(raw_overlap["alpha_plus"])),
        "alpha_minus": _complex_dict(complex(raw_overlap["alpha_minus"])),
        "power_plus": float(raw_overlap["power_plus"]),
        "power_minus": float(raw_overlap["power_minus"]),
        "i1": _complex_dict(complex(raw_overlap["i1"])),
        "i2": _complex_dict(complex(raw_overlap["i2"])),
    }
    if _is_root():
        meta_path.write_text(
            json.dumps({"raw_overlap": raw_overlap_json}, indent=2) + "\n",
            encoding="utf-8",
        )
    _barrier()
    return {
        "mode": mode,
        "cache_path": str(path),
        "cached": False,
        "raw_overlap": raw_overlap_json,
    }


@dataclass
class PortObservation:
    alpha_plus: complex
    alpha_minus: complex
    power_plus: float
    power_minus: float

    def as_dict(self) -> Dict[str, Any]:
        return {
            "alpha_plus": _complex_dict(self.alpha_plus),
            "alpha_minus": _complex_dict(self.alpha_minus),
            "power_plus": self.power_plus,
            "power_minus": self.power_minus,
        }


def run_case(
    source_port: int,
    receive_port: int,
    *,
    res: int,
    run_time: float,
    rho: np.ndarray,
    mode_by_port: Dict[int, EMLineMode],
    device_mode: str,
) -> Dict[str, Any]:
    set_geometry_context(res=res)
    _pmm, device, _wp = build_circulator_device(
        rho,
        np.array([0.0, 0.0, 0.0]),
        res=res,
        wall_pec=True,
        device_mode=device_mode,
    )
    device.sources = make_yee_snapped_numerical_mode_sources(source_port)
    sim = device.Get_Sim()

    monitors = {}
    for port_index, mode in mode_by_port.items():
        monitors[port_index] = add_em_line_monitor(
            sim,
            center_xy=np.asarray(mode.center_xy, dtype=float),
            tangent_xy=np.asarray(mode.tangent_xy, dtype=float),
            offsets_a=mode.offsets_a,
            frequency=fs_a,
        )

    sim.run(until_after_sources=run_time)
    obs = {}
    for port_index, mode in mode_by_port.items():
        fields = sample_em_line_fields(sim, monitors[port_index])
        _ensure_finite(f"P{port_index+1} case Ex", fields["ex"])
        _ensure_finite(f"P{port_index+1} case Ey", fields["ey"])
        _ensure_finite(f"P{port_index+1} case Hz", fields["hz"])
        overlap = mode.overlap_amplitudes(
            ex_field=fields["ex"], ey_field=fields["ey"], hz_field=fields["hz"]
        )
        obs[port_index] = PortObservation(
            alpha_plus=complex(overlap["alpha_plus"]),
            alpha_minus=complex(overlap["alpha_minus"]),
            power_plus=float(overlap["power_plus"]),
            power_minus=float(overlap["power_minus"]),
        )

    return {
        "device_mode": device_mode,
        "source_port": source_port,
        "receive_port": receive_port,
        "observations": {f"P{k+1}": v.as_dict() for k, v in obs.items()},
        "transmitted_power": float(obs[receive_port].power_plus),
        "source_forward_power": float(obs[source_port].power_plus),
        "source_backward_power": float(obs[source_port].power_minus),
    }


def db_diff(x: float, y: float) -> float:
    eps = 1e-300
    return float(10.0 * np.log10(max(x, eps) / max(y, eps)))


def main() -> None:
    ap = argparse.ArgumentParser()
    add_resolution_arguments(ap)
    ap.add_argument("--run-time", type=float, default=20.0)
    ap.add_argument("--force-modes", action="store_true")
    ap.add_argument(
        "--device-mode",
        choices=("horns_only", "full"),
        default="horns_only",
    )
    ap.add_argument(
        "--output",
        type=str,
        default="",
    )
    args = ap.parse_args()

    sim_res, _ppc, res_report = resolve_simulation_resolution(args)
    set_geometry_context(res=sim_res)
    info = mpi_info()
    if int(info.get("ranks", 1)) != 32:
        raise RuntimeError(f"Expected 32 MPI ranks, got {info}")

    t0 = time.time()
    mode_meta = {}
    mode_by_port: Dict[int, EMLineMode] = {}
    for port in (0, 1):
        meta = extract_reference_mode(
            port, res=sim_res, run_time=args.run_time, force=args.force_modes
        )
        mode = meta.pop("mode")
        mode_meta[f"P{port+1}"] = {
            **meta,
            "mode_power": mode.mode_power,
            "n_samples": int(mode.offsets_a.size),
            "center_xy": list(mode.center_xy),
            "tangent_xy": list(mode.tangent_xy),
            "outward_dir": list(mode.outward_dir),
        }
        # Incident wave for the straight-guide reference is inward, i.e. alpha_minus.
        raw = meta.get("raw_overlap")
        if raw is not None:
            mode_meta[f"P{port+1}"]["reference_alpha_plus"] = raw["alpha_plus"]
            mode_meta[f"P{port+1}"]["reference_alpha_minus"] = raw["alpha_minus"]
            mode_meta[f"P{port+1}"]["reference_power_plus"] = float(raw["power_plus"])
            mode_meta[f"P{port+1}"]["reference_power_minus"] = float(raw["power_minus"])
        mode_by_port[port] = mode

    # Straight-guide normalization convention: incident power is the inward coefficient.
    incident_p1 = float(mode_meta["P1"].get("reference_power_minus", 1.0))
    incident_p2 = float(mode_meta["P2"].get("reference_power_minus", 1.0))

    output_path = args.output
    if not output_path:
        stem = f"{args.device_mode}_P1P2_ppc50_rt20.json"
        output_path = str(ROOT / "outputs" / "validation" / "eh_modal" / stem)

    rho = default_uniform_rho()
    p1_to_p2 = run_case(
        0,
        1,
        res=sim_res,
        run_time=args.run_time,
        rho=rho,
        mode_by_port=mode_by_port,
        device_mode=args.device_mode,
    )
    p2_to_p1 = run_case(
        1,
        0,
        res=sim_res,
        run_time=args.run_time,
        rho=rho,
        mode_by_port=mode_by_port,
        device_mode=args.device_mode,
    )

    s21 = float(p1_to_p2["transmitted_power"] / max(incident_p1, 1e-300))
    s12 = float(p2_to_p1["transmitted_power"] / max(incident_p2, 1e-300))
    diff_db = abs(db_diff(s21, s12))

    out = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "mpi": info,
        "device_mode": args.device_mode,
        "resolution_report": res_report,
        "run_time": args.run_time,
        "mode_formula": {
            "mode_power": "Re integral(E_t* conj(H_z) dl) = 1",
            "alpha_plus": "0.5 * integral(conj(E_t_mode) * H_z + E_t * conj(H_z_mode)) dl",
            "alpha_minus": "0.5 * integral(conj(E_t_mode) * H_z - E_t * conj(H_z_mode)) dl",
        },
        "reference_modes": mode_meta,
        "cases": {
            "P1_to_P2": p1_to_p2,
            "P2_to_P1": p2_to_p1,
        },
        "incident_modal_powers": {"P1": incident_p1, "P2": incident_p2},
        "transmitted_modal_powers": {
            "S21_numer": p1_to_p2["transmitted_power"],
            "S12_numer": p2_to_p1["transmitted_power"],
        },
        "normalized_transmission": {
            "S21": s21,
            "S12": s12,
            "difference_db": diff_db,
            "known_direct_poynting_difference_db": 0.0861,
        },
        "elapsed_s": time.time() - t0,
    }

    if _is_root():
        out_path = Path(output_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(out["normalized_transmission"], indent=2))


if __name__ == "__main__":
    main()
