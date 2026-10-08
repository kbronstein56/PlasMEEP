#!/usr/bin/env python3
"""Tiny CW / solve_cw compatibility probe (benchmark-only)."""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

import meep as mp
import numpy as np

OUT = Path(__file__).resolve().parent


def try_case(name: str, medium_factory, use_gyro: bool = False) -> dict:
    result = {"name": name, "ok": False, "error": None, "wall_s": None, "solve_cw_return": None}
    try:
        cell = mp.Vector3(4, 4)
        geom = [
            mp.Cylinder(radius=0.5, center=mp.Vector3(), material=medium_factory()),
        ]
        sources = [
            mp.Source(
                mp.ContinuousSource(frequency=0.25, is_integrated=True),
                component=mp.Hz,
                center=mp.Vector3(-1.2, 0),
            )
        ]
        sim = mp.Simulation(
            cell_size=cell,
            geometry=geom,
            sources=sources,
            resolution=20,
            boundary_layers=[mp.PML(0.5)],
        )
        # flux at fs
        fr = mp.FluxRegion(center=mp.Vector3(1.2, 0), size=mp.Vector3(0, 1.0))
        flux = sim.add_flux(0.25, 0, 1, fr)
        t0 = time.perf_counter()
        sim.init_sim()
        # solve_cw after fields exist
        ret = sim.solve_cw(tol=1e-6, maxiters=2000, L=10)
        wall = time.perf_counter() - t0
        fluxes = mp.get_fluxes(flux)
        result.update(
            ok=True,
            wall_s=wall,
            solve_cw_return=str(ret)[:200],
            flux=float(fluxes[0]) if fluxes else None,
            use_gyro=use_gyro,
        )
    except Exception as e:
        result["error"] = f"{type(e).__name__}: {e}"
        result["traceback"] = traceback.format_exc()[-1500:]
    return result


def main() -> int:
    cases = []
    # vacuum control
    cases.append(try_case("vacuum", lambda: mp.Medium(epsilon=1)))
    # ordinary Drude (B=0 plasma-like)
    cases.append(
        try_case(
            "drude",
            lambda: mp.Medium(
                epsilon=1,
                E_susceptibilities=[
                    mp.DrudeSusceptibility(frequency=0.5, gamma=1e-4, sigma=1.0)
                ],
            ),
        )
    )
    # gyrotropic Drude (B!=0)
    cases.append(
        try_case(
            "gyrotropic_drude",
            lambda: mp.Medium(
                epsilon=1,
                E_susceptibilities=[
                    mp.GyrotropicDrudeSusceptibility(
                        frequency=0.5, gamma=1e-4, sigma=1.0, bias=mp.Vector3(0, 0, 0.1)
                    )
                ],
            ),
            use_gyro=True,
        )
    )
    out = OUT / "cw_compat.json"
    out.write_text(json.dumps(cases, indent=2) + "\n")
    print(json.dumps(cases, indent=2))
    print("Wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
