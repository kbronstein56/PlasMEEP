#!/usr/bin/env python3
"""Build and verify C3 (120°) rotational orbits of the 91-bulb lattice.

Expected: 1 center + 30 threefold orbits = 31 design variables.
Does not modify production C6 sixfold_orbits.py.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
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

OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"


def rotate120(xy, center, k: int):
    th = k * 2.0 * np.pi / 3.0
    c, s = np.cos(th), np.sin(th)
    R = np.array([[c, -s], [s, c]])
    return (xy - center) @ R.T + center


def build_orbits_c3(centers: np.ndarray, tol: float = 1e-6):
    centers = np.asarray(centers, float)
    origin = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    n = len(centers)
    used = np.zeros(n, dtype=bool)
    orbits = []
    for i in range(n):
        if used[i]:
            continue
        members = []
        for k in range(3):
            target = rotate120(centers[i], origin, k)
            d = np.linalg.norm(centers - target, axis=1)
            j = int(np.argmin(d))
            if d[j] > tol:
                raise RuntimeError(f"C3 orbit break at rod {i} rot {k}: min dist {d[j]}")
            members.append(j)
        members = sorted(set(members))
        for m in members:
            used[m] = True
        r = float(np.linalg.norm(centers[members[0]] - origin))
        orbits.append({"rods": members, "radius_a": r, "size": len(members)})
    orbits.sort(key=lambda o: (o["radius_a"], o["rods"][0]))
    return orbits, origin


def expand_tied(q: np.ndarray, orbits) -> np.ndarray:
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


def map_c6_q_to_c3(q_c6, orbits_c6, orbits_c3) -> np.ndarray:
    """Neutral warm start: expand C6 design to rods, then average within C3 orbits."""
    from sixfold_orbits import expand_tied as expand_c6

    s = expand_c6(np.asarray(q_c6, float), orbits_c6)
    q3 = np.zeros(len(orbits_c3), float)
    for m, o in enumerate(orbits_c3):
        q3[m] = float(np.mean([s[i] for i in o["rods"]]))
    return q3


def main():
    centers = bulb_centers_device(sc)
    orbits, origin = build_orbits_c3(centers)
    sizes = [o["size"] for o in orbits]
    assert sum(sizes) == 91, sizes
    n_center = sum(1 for o in orbits if o["size"] == 1)
    n_three = sum(1 for o in orbits if o["size"] == 3)
    # C3 may also have size-1 only at center; all others should be 3
    ok_sizes = n_center == 1 and n_three == 30 and len(orbits) == 31
    ok_geo = True
    for o in orbits:
        for i in o["rods"]:
            for k in range(3):
                t = rotate120(centers[i], origin, k)
                d = np.linalg.norm(centers[o["rods"]] - t, axis=1).min()
                if d > 1e-6:
                    ok_geo = False
    report = {
        "n_rods": 91,
        "n_orbits": len(orbits),
        "n_center_orbits": n_center,
        "n_three_orbits": n_three,
        "n_design_vars": len(orbits),
        "origin": origin.tolist(),
        "orbits": orbits,
        "geometry_ok": bool(ok_geo),
        "pass": bool(ok_sizes and ok_geo),
        "port_classes": {
            "odd_ports": [0, 2, 4],
            "even_ports": [1, 3, 5],
            "note": "Under exact C3, {1,3,5} and {2,4,6} are inequivalent rotational classes.",
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "threefold_orbits.json").write_text(json.dumps(report, indent=2) + "\n")
    md = [
        "# Threefold (C3) design orbits",
        "",
        f"Verified: **{len(orbits)}** orbits = 1 center + 30×3 → **{len(orbits)}** design variables.",
        f"Domain center: `({origin[0]:.3f}, {origin[1]:.3f})` a.",
        "",
        "Port rotational classes under C3:",
        "- `{1,3,5}` (odd)",
        "- `{2,4,6}` (even)",
        "",
        f"**C3_ORBITS: {'PASS' if report['pass'] else 'FAIL'}**",
        "",
        "| orbit | size | radius (a) | rod indices |",
        "|------:|-----:|-----------:|-------------|",
    ]
    for m, o in enumerate(orbits):
        md.append(f"| {m} | {o['size']} | {o['radius_a']:.6f} | {o['rods']} |")
    (OUT / "THREEFOLD_DESIGN_ORBITS.md").write_text("\n".join(md) + "\n")
    print("C3 orbits", len(orbits), "pass", report["pass"], flush=True)


if __name__ == "__main__":
    main()
