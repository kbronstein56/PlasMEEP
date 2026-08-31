#!/usr/bin/env python3
"""
Controlled EigenModeSource test on a straight PEC parallel-plate reference guide.

Documents why full PEC-horn EigenModeSource failed and whether a simplified
geometry can supply a numerical mode.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict

import meep as mp
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VAL = os.path.dirname(os.path.abspath(__file__))
for p in (VAL, os.path.join(ROOT, "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from sixport_common import a, clear_width, fs_a, source_df  # noqa: E402

OUT = Path(ROOT) / "outputs" / "validation" / "eigenmode_reference"


def build_straight_guide_sim(res: int, guide_length_a: float = 3.0) -> mp.Simulation:
    """Infinite PEC walls along x, clear aperture = clear_width."""
    dpml = 0.5
    sx = guide_length_a + 2 * dpml
    sy = clear_width + 1.0
    cell = mp.Vector3(sx, sy, 0)

    # PEC top/bottom plates (high-eps blocks)
    wall_eps = 1e20
    wall_h = 0.5
    geometry = [
        mp.Block(
            center=mp.Vector3(0, clear_width / 2 + wall_h / 2, 0),
            size=mp.Vector3(sx, wall_h, mp.inf),
            material=mp.Medium(epsilon=wall_eps),
        ),
        mp.Block(
            center=mp.Vector3(0, -clear_width / 2 - wall_h / 2, 0),
            size=mp.Vector3(sx, wall_h, mp.inf),
            material=mp.Medium(epsilon=wall_eps),
        ),
    ]

    pml_layers = [mp.PML(thickness=dpml)]
    return mp.Simulation(
        cell_size=cell,
        resolution=res,
        boundary_layers=pml_layers,
        geometry=geometry,
        default_material=mp.air,
    )


def run_eigenmode_launch(
    sim: mp.Simulation,
    *,
    frequency: float,
    fwidth: float,
    run_time: float,
) -> Dict[str, Any]:
    """EigenModeSource at x=-0.8a, flux monitors upstream/downstream."""
    src_x = -0.8
    mon_up_x = -1.2
    mon_dn_x = 0.8
    span_y = 0.9 * clear_width

    src = mp.EigenModeSource(
        src=mp.GaussianSource(frequency=frequency, fwidth=fwidth),
        center=mp.Vector3(src_x, 0, 0),
        size=mp.Vector3(0, span_y, 0),
        direction=mp.X,
        eig_band=1,
        eig_parity=mp.ODD_Z,
        eig_match_freq=True,
        eig_kpoint=mp.Vector3(1, 0, 0),
    )
    sim.sources = [src]

    flux_up = sim.add_flux(
        frequency, 0, 1, mp.FluxRegion(center=mp.Vector3(mon_up_x, 0, 0), size=mp.Vector3(0, span_y, 0))
    )
    flux_dn = sim.add_flux(
        frequency, 0, 1, mp.FluxRegion(center=mp.Vector3(mon_dn_x, 0, 0), size=mp.Vector3(0, span_y, 0))
    )

    # Line DFT for profile at downstream monitor
    n_pts = 41
    offsets = np.linspace(-span_y / 2, span_y / 2, n_pts)
    dft_objs = []
    for oy in offsets:
        dft_objs.append(
            sim.add_dft_fields(
                [mp.Hz],
                frequency,
                0,
                1,
                center=mp.Vector3(mon_dn_x, oy, 0),
                size=mp.Vector3(0, 0, 0),
            )
        )

    t0 = time.perf_counter()
    sim.run(until_after_sources=run_time)
    elapsed = time.perf_counter() - t0

    up = abs(complex(np.squeeze(mp.get_fluxes(flux_up)[0])))
    dn = abs(complex(np.squeeze(mp.get_fluxes(flux_dn)[0])))
    hz = [
        complex(np.squeeze(sim.get_dft_array(d, mp.Hz, 0))) for d in dft_objs
    ]

    return {
        "flux_upstream": up,
        "flux_downstream": dn,
        "transmission_ratio": dn / max(up, 1e-30),
        "hz_profile": {
            "offsets_a": offsets.tolist(),
            "hz_complex": [{"real": z.real, "imag": z.imag} for z in hz],
        },
        "wall_time_s": elapsed,
    }


def run_pec_point_source_baseline(
    sim: mp.Simulation,
    *,
    frequency: float,
    fwidth: float,
    run_time: float,
) -> Dict[str, Any]:
    """Time-domain TE1 cos line source (no MPB) for comparison."""
    span_y = 0.9 * clear_width
    offsets = np.linspace(-span_y / 2, span_y / 2, 31)
    weights = np.cos(np.pi * offsets / span_y).clip(min=0.0)
    weights = weights / max(np.sum(weights), 1e-30)

    sources = []
    for oy, w in zip(offsets, weights):
        sources.append(
            mp.Source(
                mp.GaussianSource(frequency=frequency, fwidth=fwidth),
                component=mp.Hz,
                center=mp.Vector3(-0.8, oy, 0),
                amplitude=float(w),
            )
        )
    sim.sources = sources

    mon_dn_x = 0.8
    flux_dn = sim.add_flux(
        frequency, 0, 1, mp.FluxRegion(center=mp.Vector3(mon_dn_x, 0, 0), size=mp.Vector3(0, span_y, 0))
    )
    t0 = time.perf_counter()
    sim.run(until_after_sources=run_time)
    elapsed = time.perf_counter() - t0
    dn = abs(complex(np.squeeze(mp.get_fluxes(flux_dn)[0])))
    return {"flux_downstream": dn, "wall_time_s": elapsed}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, default=32)
    parser.add_argument("--run-time", type=float, default=30.0)
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)

    eigen_sim = build_straight_guide_sim(args.res)
    eigen_result = run_eigenmode_launch(
        eigen_sim,
        frequency=fs_a,
        fwidth=source_df,
        run_time=args.run_time,
    )

    td_sim = build_straight_guide_sim(args.res)
    td_result = run_pec_point_source_baseline(
        td_sim,
        frequency=fs_a,
        fwidth=source_df,
        run_time=args.run_time,
    )

    eigen_ok = eigen_result["flux_downstream"] > 1e-6 and eigen_result["transmission_ratio"] > 0.1

    payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "settings": {"res": args.res, "run_time": args.run_time, "a_m": a},
        "full_horn_eigenmode_failure": {
            "reasons": [
                "PEC invalid for MPB (perfect conductor not representable in MPB dielectric solve)",
                "near-zero launched flux on diagonal ports in 6-port staircased horns",
                "unphysical MPB |k| for f~0.36 in angled NO_DIRECTION sources",
                "eigenmode coefficients near zero on 60° feeds",
            ],
            "source": "outputs/validation/ports/eigenmode_verdict.json",
        },
        "straight_guide_eigenmode": eigen_result,
        "straight_guide_te1_cos_td": td_result,
        "verdict": {
            "eigenmode_viable_on_straight_guide": bool(eigen_ok),
            "recommendation": (
                "Use time-domain reference-guide field extraction for angled PEC horns; "
                "EigenModeSource may work only on simplified high-eps or straight guides."
                if eigen_ok
                else "Fall back to time-domain reference-guide extraction; EigenModeSource "
                "failed even on straight PEC guide at this resolution."
            ),
        },
    }

    out_path = OUT / f"eigenmode_reference_res{args.res}.json"
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
        f.write("\n")

    print(f"Eigen downstream flux: {eigen_result['flux_downstream']:.6g}")
    print(f"TD te1 cos downstream flux: {td_result['flux_downstream']:.6g}")
    print(f"Eigen viable: {eigen_ok}")
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
