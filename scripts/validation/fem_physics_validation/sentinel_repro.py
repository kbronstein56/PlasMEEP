#!/usr/bin/env python3
"""Small sentinel reruns of committed B=0 checks.

Writes outputs/validation/fem_physics_validation/sentinel_repro.json.
Does not overwrite the stored validation JSON. Select one group:

    mms, planar, oblique, recip, power, guide, bare, three, seven, coated
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation" / "sentinel_repro.json"


def _load():
    if OUT.exists():
        return json.loads(OUT.read_text())
    return {}


def _save(blob):
    OUT.write_text(json.dumps(blob, indent=2) + "\n")


def run_mms():
    from mms import CASES, run_case

    rec = run_case("D_full_complex", CASES["D_full_complex"], [64])
    return rec["levels"][-1]


def run_planar():
    from planar_fem import driven_slabs, pec_guide

    guide = pec_guide([0.025])[0]
    slabs = driven_slabs([0.02])
    plasma = [r for r in slabs if r["stack"] == "plasma_0.40" and abs(r["f_GHz"] - 3.85) < 1e-6][0]
    return {"pec_guide_h0025": guide, "plasma_slab_0.40_h002": plasma}


def run_oblique():
    from oblique_pml_fit import solve_one

    rows = []
    for ang in (0.0, 25.0, 60.0):
        rows.append(solve_one(ang, 1.2, 1.0, 2.4, 1.0, 0.02))
    return rows


def run_recip():
    from conservation import reciprocity

    return reciprocity()


def run_power():
    from conservation import absorption_sign

    return absorption_sign()


def run_guide():
    from guide_eh_close import guide_close

    return guide_close(0.00064)


def _scatter(name, centers, radius, coated, levels, fs_scale):
    import sixport_common as sc
    from analytic_sweeps import plasma
    from fem_scatterers import solve_case

    saved = float(sc.fs_a)
    try:
        sc.fs_a = saved * fs_scale
        eps = plasma(sc.fs_a)
        rows = solve_case(name, centers, radius, eps, coated, levels)
    finally:
        sc.fs_a = saved
    rec = rows[0]
    # Drop nothing large; probes are small.
    return rec


def run_bare():
    return _scatter("sentinel_bare", np.zeros((1, 2)), 0.230, False, [(0.04, 0.012)], 1.0)


def run_three():
    from analytic_sweeps import _three

    return _scatter("sentinel_three", _three(), 0.230, True, [(0.04, 0.012)], 1.0)


def run_seven():
    from analytic_sweeps import _seven

    return _scatter("sentinel_seven", _seven(), 0.230, True, [(0.02, 0.0045)], 1.0)


def run_coated():
    return _scatter("sentinel_coated_530", np.zeros((1, 2)), 0.230, True, [(0.02, 0.0022, 0.0032, 3.2)], 5.30 / 3.85)


def main():
    group = sys.argv[1]
    fn = {
        "mms": run_mms,
        "planar": run_planar,
        "oblique": run_oblique,
        "recip": run_recip,
        "power": run_power,
        "guide": run_guide,
        "bare": run_bare,
        "three": run_three,
        "seven": run_seven,
        "coated": run_coated,
    }[group]
    blob = _load()
    blob[group] = fn()
    _save(blob)
    print("SENTINEL_SAVED", group, flush=True)


if __name__ == "__main__":
    main()
