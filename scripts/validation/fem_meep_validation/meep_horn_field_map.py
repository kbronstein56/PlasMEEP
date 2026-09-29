#!/usr/bin/env python3
"""Dump harmonic Hz on the horns-only Meep cell for the field-difference map."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "scripts" / "validation"), str(ROOT / "scripts")]

import meep as mp  # noqa: E402
import sixport_common as sc  # noqa: E402
from physical_units import meep_resolution_from_points_per_cm  # noqa: E402
from port_formulations import get_formulation  # noqa: E402
from sixport_common import build_circulator_device, set_geometry_context  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "boundary_pml"


def rank() -> int:
    try:
        from mpi4py import MPI

        return int(MPI.COMM_WORLD.Get_rank())
    except Exception:
        return int(os.environ.get("PMI_RANK", "0"))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    res = int(meep_resolution_from_points_per_cm(25.0, a_m=sc.a))
    set_geometry_context(res=res, horn_walls="prism", grid_offset_cells=(0.0, 0.0), coord_rotation_deg=0.0)
    _pmm, dev, _wp = build_circulator_device(
        np.zeros(91), np.zeros(3), res=res, device_mode="horns_only", wall_pec=True
    )
    dev.sources = get_formulation("num_mode_guide_normal").make_sources(0)
    sim = dev.Get_Sim()
    vol = mp.Volume(center=mp.Vector3(0, 0, 0), size=mp.Vector3(sc.nx_ports, sc.ny_ports, 0))
    dft = sim.add_dft_fields([mp.Hz], float(sc.fs_a), 0, 1, where=vol)
    if rank() == 0:
        print("running meep field map", flush=True)
    sim.run(until_after_sources=20)
    hz = np.array(sim.get_dft_array(dft, mp.Hz, 0))
    if rank() != 0:
        return 0
    # Meep 2d arrays are indexed (x, y). Coordinates span the centered cell.
    nx, ny = hz.shape
    xs = np.linspace(-sc.nx_ports / 2.0, sc.nx_ports / 2.0, nx)
    ys = np.linspace(-sc.ny_ports / 2.0, sc.ny_ports / 2.0, ny)
    np.savez_compressed(OUT / "meep_hz_grid.npz", Hz=hz, x=xs, y=ys)
    print("saved", hz.shape, "max", float(np.max(np.abs(hz))), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
