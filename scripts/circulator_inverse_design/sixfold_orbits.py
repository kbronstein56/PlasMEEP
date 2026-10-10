#!/usr/bin/env python3
"""Build and verify C6 rotational orbits of the 91-bulb lattice."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sixport_common as sc  # noqa: E402
from full91_geometry import bulb_centers_device  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design"


def rotate60(xy, center, k: int):
    th = k * np.pi / 3.0
    c, s = np.cos(th), np.sin(th)
    R = np.array([[c, -s], [s, c]])
    return (xy - center) @ R.T + center


def build_orbits(centers: np.ndarray, tol: float = 1e-6):
    centers = np.asarray(centers, float)
    origin = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    n = len(centers)
    used = np.zeros(n, dtype=bool)
    orbits = []
    for i in range(n):
        if used[i]:
            continue
        members = []
        for k in range(6):
            target = rotate60(centers[i], origin, k)
            d = np.linalg.norm(centers - target, axis=1)
            j = int(np.argmin(d))
            if d[j] > tol:
                raise RuntimeError(f"orbit break at rod {i} rot {k}: min dist {d[j]}")
            members.append(j)
        members = sorted(set(members))
        for m in members:
            used[m] = True
        r = float(np.linalg.norm(centers[members[0]] - origin))
        orbits.append({"rods": members, "radius_a": r, "size": len(members)})
    orbits.sort(key=lambda o: (o["radius_a"], o["rods"][0]))
    return orbits, origin


def expand_tied(q: np.ndarray, orbits) -> np.ndarray:
    n = sum(o["size"] for o in orbits)
    # n should be 91; build from max rod index
    n = 1 + max(max(o["rods"]) for o in orbits)
    s = np.zeros(n, float)
    for m, o in enumerate(orbits):
        for i in o["rods"]:
            s[i] = q[m]
    return s


def reduce_gradient(g_s: np.ndarray, orbits) -> np.ndarray:
    g_q = np.zeros(len(orbits), float)
    for m, o in enumerate(orbits):
        g_q[m] = float(sum(g_s[i] for i in o["rods"]))
    return g_q


def main():
    centers = bulb_centers_device(sc)
    orbits, origin = build_orbits(centers)
    sizes = [o["size"] for o in orbits]
    assert sum(sizes) == 91, sizes
    n_center = sum(1 for o in orbits if o["size"] == 1)
    n_six = sum(1 for o in orbits if o["size"] == 6)
    assert n_center == 1 and n_six == 15, (n_center, n_six, sizes)
    # geometric check: rotating any non-center rod stays in its orbit
    ok = True
    for o in orbits:
        for i in o["rods"]:
            for k in range(6):
                t = rotate60(centers[i], origin, k)
                d = np.linalg.norm(centers[o["rods"]] - t, axis=1).min()
                if d > 1e-6:
                    ok = False
    report = {
        "n_rods": 91,
        "n_orbits": len(orbits),
        "n_center_orbits": n_center,
        "n_six_orbits": n_six,
        "origin": origin.tolist(),
        "orbits": orbits,
        "geometry_ok": bool(ok),
        "pass": bool(ok and n_center == 1 and n_six == 15),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "sixfold_orbits.json").write_text(json.dumps(report, indent=2) + "\n")
    md = [
        "# Sixfold design orbits",
        "",
        f"Verified geometrically: **{len(orbits)}** orbits = 1 center + 15×6.",
        f"Domain center used for C6: `({origin[0]:.3f}, {origin[1]:.3f})` a.",
        "",
        "| orbit | size | radius (a) | rod indices |",
        "|------:|-----:|-----------:|-------------|",
    ]
    for m, o in enumerate(orbits):
        md.append(f"| {m} | {o['size']} | {o['radius_a']:.6f} | {o['rods']} |")
    md += [
        "",
        "Tied variables: `s_i = q_m` for `i ∈ orbit m`.",
        "Gradient: `∂J/∂q_m = ∑_{i∈m} ∂J/∂s_i`.",
        "",
        f"**SIXFOLD_ORBITS: {'PASS' if report['pass'] else 'FAIL'}**",
        "",
    ]
    (OUT / "SIXFOLD_DESIGN_ORBITS.md").write_text("\n".join(md))
    print("orbits", len(orbits), "pass", report["pass"], flush=True)


if __name__ == "__main__":
    main()
