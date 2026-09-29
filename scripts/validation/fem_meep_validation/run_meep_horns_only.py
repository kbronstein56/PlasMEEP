#!/usr/bin/env python3
"""Trusted Meep six-horn geometry, no quartz and no plasma. Source P1, receivers P1–P6."""
from __future__ import annotations

import json
import os
import pickle
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
VAL = ROOT / "scripts" / "validation"
sys.path[:0] = [str(VAL), str(ROOT / "scripts")]

import meep as mp  # noqa: E402
import sixport_common as sc  # noqa: E402
from physical_units import meep_resolution_from_points_per_cm  # noqa: E402
from port_formulations import (  # noqa: E402
    add_flux_monitor,
    extract_flux_powers,
    get_formulation,
    make_flux_region_for_formulation,
    port_measure_center,
)
from sixport_common import build_circulator_device, default_uniform_rho, set_geometry_context  # noqa: E402

OUT = Path(os.environ.get("MEEP_OUT", str(ROOT / "outputs" / "validation" / "fem_meep_validation" / "horn_localization")))
FORMULATION = "num_mode_guide_normal"


def _rank() -> int:
    try:
        from mpi4py import MPI

        return int(MPI.COMM_WORLD.Get_rank())
    except Exception:
        return int(os.environ.get("PMI_RANK", os.environ.get("OMPI_COMM_WORLD_RANK", "0")))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    points_per_cm = float(os.environ.get("HORN_PPC", "25"))
    run_time = float(os.environ.get("HORN_RUN_TIME", "20"))
    res = int(meep_resolution_from_points_per_cm(points_per_cm, a_m=sc.a))
    set_geometry_context(
        res=res,
        horn_walls="prism",
        grid_offset_cells=(0.0, 0.0),
        monitor_offset_cells=(0.0, 0.0),
        coord_rotation_deg=0.0,
    )
    cache_dir = VAL / ".cache"
    stem = f"norm_{FORMULATION}_full_prism_rot0_pec_g0_0_m0_0_ppc{points_per_cm:g}_res{res}_rt20_p0"
    cands = sorted(cache_dir.glob(stem + "*.pkl"))
    exact = cache_dir / f"{stem}.pkl"
    norm_path = exact if exact.is_file() else (cands[0] if cands else None)
    if norm_path is None:
        raise FileNotFoundError(f"no incident-power cache for ppc={points_per_cm:g} res={res}")
    with open(norm_path, "rb") as f:
        cache = pickle.load(f)
    P_inc = float(cache["incident_power_by_port"][0])
    form = get_formulation(FORMULATION)
    from sixport_common import default_uniform_rho

    material = os.environ.get("MEEP_MATERIAL", "horns")
    if material == "quartz":
        rho = np.zeros(91)
        device_mode, plasma_fill = "full", "geometry_only"
    elif material == "plasma":
        rho = default_uniform_rho()
        device_mode, plasma_fill = "full", "active"
    elif material == "horns":
        rho = np.zeros(91)
        device_mode, plasma_fill = "horns_only", "active"
    else:
        raise ValueError(f"unknown MEEP_MATERIAL={material}")
    _pmm, P_device, _wp = build_circulator_device(
        rho,
        np.zeros(3),
        res=res,
        device_mode=device_mode,
        wall_pec=True,
        plasma_fill=plasma_fill,
    )
    if material == "plasma":
        if abs(float(_wp[0]) - float(sc.fp_a)) > 1e-9:
            raise RuntimeError(f"plasma frequency {_wp[0]} is not fp_a={sc.fp_a}")
        # Meep prints the instantaneous epsilon; the Drude frequency must be fp.
        n_drude = 0
        for geom in P_device.geometry:
            sus = list(geom.material.E_susceptibilities or [])
            if not sus:
                continue
            freq = float(sus[0].frequency)
            if abs(freq - float(sc.fp_a)) > 1e-6:
                raise RuntimeError(f"Drude frequency {freq} is not fp_a={sc.fp_a}")
            n_drude += 1
        if n_drude != 91:
            raise RuntimeError(f"expected 91 Drude bulbs, found {n_drude}")
        if _rank() == 0:
            print(f"PLASMA_CHECK fp_a={sc.fp_a:.6f} gamma_a={sc.gamma_a:.6e} n_drude={n_drude}", flush=True)
    P_device.sources = form.make_sources(0)
    sim = mp.Simulation(
        cell_size=P_device.cell,
        boundary_layers=P_device.pml,
        geometry=P_device.geometry,
        sources=P_device.sources,
        resolution=P_device.res,
        default_material=mp.Medium(epsilon=1, mu=1),
        eps_averaging=os.environ.get("HORN_EPS_AVG", "1") != "0",
    )
    monitors = []
    signs = []
    centers = {}
    for out_p in range(6):
        measure_xy = port_measure_center(form.name, out_p)
        centers[out_p] = np.asarray(measure_xy, dtype=float).tolist()
        regions, sign = make_flux_region_for_formulation(
            form.name, measure_xy, sc.effective_port_dir(out_p), out_p
        )
        monitors.append(add_flux_monitor(sim, regions))
        signs.append(sign)
    sim.load_minus_flux_data(monitors[0], cache["incident_flux_data_by_port"][0])
    skip_fields = os.environ.get("HORN_SKIP_FIELDS", "0") == "1" or abs(points_per_cm - 25.0) > 1e-6

    fs = float(sc.fs_a)
    dfts = {}
    if not skip_fields:
        span = float(0.96 * sc.clear_width)
        h0 = sc.horn_for_port(0, res)
        n_hat = np.asarray(sc.effective_port_dir(0), dtype=float)
        n_hat = n_hat / np.linalg.norm(n_hat)
        lines = {
            "source": np.asarray(h0["source_center"], dtype=float),
            "throat": np.asarray(h0["throat_center"], dtype=float),
            "monitor": np.asarray(h0["monitor_center"], dtype=float),
            "center": np.zeros(2),
        }
        for name, xy in lines.items():
            if name == "center":
                vol = mp.Volume(center=mp.Vector3(0, 0, 0), size=mp.Vector3(8.0, 0, 0))
            elif abs(n_hat[0]) >= abs(n_hat[1]):
                vol = mp.Volume(center=mp.Vector3(float(xy[0]), float(xy[1]), 0), size=mp.Vector3(0, span, 0))
            else:
                vol = mp.Volume(center=mp.Vector3(float(xy[0]), float(xy[1]), 0), size=mp.Vector3(span, 0, 0))
            dfts[name] = sim.add_dft_fields([mp.Hz, mp.Ex, mp.Ey], fs, 0, 1, where=vol)

    if _rank() == 0:
        print(
            f"Meep horns_only ppc={points_per_cm:g} res={res} rt={run_time} "
            f"P_inc={P_inc:.6e} norm={norm_path.name}",
            flush=True,
        )
    t0 = time.perf_counter()
    sim.run(until_after_sources=run_time)
    wall = time.perf_counter() - t0
    flux = extract_flux_powers(sim, monitors, signs)
    norm = flux / P_inc
    if _rank() == 0:
        print("POWERS", norm.tolist(), flush=True)
    # get_dft_array is collective: every rank must call it.
    fields = {}
    for name, dft in dfts.items():
        fields[name] = {
            "Hz": np.array(sim.get_dft_array(dft, mp.Hz, 0)),
            "Ex": np.array(sim.get_dft_array(dft, mp.Ex, 0)),
            "Ey": np.array(sim.get_dft_array(dft, mp.Ey, 0)),
        }
    if _rank() != 0:
        return 0
    tag = f"meep_{material}"
    if abs(points_per_cm - 25.0) > 1e-6 or material != "horns":
        tag += f"_ppc{points_per_cm:g}"
    if abs(run_time - 20.0) > 1e-9:
        tag += f"_rt{run_time:g}"
    if os.environ.get("HORN_EPS_AVG", "1") == "0":
        tag += "_noavg"
    if fields:
        np.savez_compressed(
            OUT / f"{tag}_fields.npz",
            **{f"{k}_{c}": v for k, d in fields.items() for c, v in d.items()},
        )
    dx_a = 1.0 / float(res)
    out = {
        "device_mode": device_mode,
        "plasma_fill": plasma_fill,
        "material": material,
        "formulation": FORMULATION,
        "res": res,
        "points_per_cm": points_per_cm,
        "dx_a": dx_a,
        "dx_mm": dx_a * float(sc.a) * 1e3,
        "run_time": run_time,
        "eps_averaging": os.environ.get("HORN_EPS_AVG", "1") != "0",
        "norm_cache": norm_path.name,
        "P_inc": P_inc,
        "raw_flux": flux.tolist(),
        "normalized_power": norm.tolist(),
        "power_dB": [float(10 * np.log10(v)) if v > 0 else None for v in norm],
        "monitor_centers_a": centers,
        "wall_s": wall,
        "note": "P11 uses load_minus_flux_data. Cached FluxData has been all zeros before, so P11 may be unsubtracted.",
    }
    (OUT / f"{tag}.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in ("raw_flux", "normalized_power", "power_dB", "wall_s")}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
