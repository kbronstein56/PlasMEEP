#!/usr/bin/env python3
"""Parameter sweeps of the independent analytic references.

No FEM and no Meep. Truncation, radii, materials, probe sets, conditioning,
and a contour power balance.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import hankel1

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from analytic_maxwell import cluster_field, mie_T, solve_clusters, stack_response  # noqa: E402
from faraday_benchmark import gyrotropic_drude_eps_eta  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"
SRC = np.array([-4.5, 0.0])
R_CORE = 0.230
R_GAP = float(sc.r_bulb_inner)
R_SHELL = float(sc.r_bulb_outer)


def plasma(f_a: float) -> complex:
    eps, _ = gyrotropic_drude_eps_eta(float(f_a), float(sc.fp_a), float(sc.gamma_a), 0.0)
    return complex(eps)


def probes(outer: float) -> dict[str, tuple[float, float]]:
    c60 = np.cos(np.deg2rad(60.0))
    s60 = np.sin(np.deg2rad(60.0))
    return {
        "forward": (2.8, 0.0),
        "backward": (-2.8, 0.0),
        "side_p60": (2.8 * c60, 2.8 * s60),
        "side_m60": (2.8 * c60, -2.8 * s60),
        "side_90": (0.0, 2.8),
        "gap_air": (0.55, 0.0),
        "outside_shell": (outer + 0.06, 0.0),
    }


def ratio_map(centers, k0, radius, eps, nmax, coated=None):
    xy = np.array(list(probes(radius if coated is None else coated[0][-1]).values()), float)
    names = list(probes(radius if coated is None else coated[0][-1]))
    b = solve_clusters(centers, k0, radius, eps, SRC, nmax, coated=coated)
    outer = radius if coated is None else coated[0][-1]
    field = cluster_field(xy, centers, k0, outer, eps, SRC, b, nmax)
    vac = hankel1(0, k0 * np.linalg.norm(xy - SRC, axis=1))
    cond = float(getattr(solve_clusters, "last_cond", np.nan))
    out = {}
    for name, z, v in zip(names, field, vac):
        if not np.isfinite(z):
            out[name] = None
        else:
            r = z / v
            out[name] = [float(r.real), float(r.imag)]
    return out, cond


def truncation_table():
    k0 = 2 * np.pi * float(sc.fs_a)
    eps = plasma(float(sc.fs_a))
    coated = ([R_CORE, R_GAP, R_SHELL], [eps, 1.0 + 0j, 3.8 + 0j])
    geometries = {
        "bare_1": (np.zeros((1, 2)), None, R_CORE, eps),
        "bare_3": (_three(), None, R_CORE, eps),
        "bare_7": (_seven(), None, R_CORE, eps),
        "coated_1": (np.zeros((1, 2)), coated, R_CORE, eps),
        "coated_3": (_three(), coated, R_CORE, eps),
        "coated_7": (_seven(), coated, R_CORE, eps),
    }
    rows = []
    prev = {}
    for nmax in (2, 3, 4, 5, 6, 8, 10, 12, 14):
        for name, (centers, coat, radius, ee) in geometries.items():
            vals, cond = ratio_map(centers, k0, radius, ee, nmax, coated=coat)
            rec = {"case": name, "nmax": nmax, "cond": cond, "forward": vals["forward"], "side_90": vals["side_90"], "gap_air": vals["gap_air"]}
            if name in prev and prev[name]["forward"] and vals["forward"]:
                a = complex(*prev[name]["forward"])
                b = complex(*vals["forward"])
                rec["d_forward"] = float(abs(b - a))
                rec["d_side"] = float(abs(complex(*vals["side_90"]) - complex(*prev[name]["side_90"])))
            prev[name] = vals
            rows.append(rec)
            print("trunc", name, nmax, "cond", f"{cond:.3e}", "dF", rec.get("d_forward"), flush=True)
    return rows


def _three():
    path = ROOT / "outputs/validation/fem_meep_validation/plasma_ring_campaign/fem_cluster3.json"
    return np.array(json.loads(path.read_text())[-1]["centers_a"], float)


def _seven():
    # Same production hex as the saved FEM file when it is available.
    path = ROOT / "outputs/validation/fem_meep_validation/plasma_ring_campaign/fem_cluster7.json"
    if path.exists():
        return np.array(json.loads(path.read_text())[-1]["centers_a"], float)
    ang = np.arange(6) * np.pi / 3.0
    ring = np.column_stack([np.cos(ang), np.sin(ang)])
    return np.vstack([np.zeros(2), ring])


def mie_matrix():
    f_a = float(sc.fs_a)
    k0 = 2 * np.pi * f_a
    eps_p = plasma(f_a)
    materials = {
        "dielectric_3.8": 3.8 + 0j,
        "neg_1.5": -1.5 + 0j,
        "neg_3.3": -3.3 + 0j,
        "neg_8": -8.0 + 0.01j,
        "plasma_fs": eps_p,
    }
    radii = {"small": 0.12, "production": R_CORE, "large": 0.40}
    rows = []
    for rname, radius in radii.items():
        for mname, eps in materials.items():
            centers = np.zeros((1, 2))
            vals, cond = ratio_map(centers, k0, radius, complex(eps), 12, coated=None)
            T = mie_T(k0, complex(eps), radius, 12)
            rows.append({
                "radius_name": rname,
                "radius": radius,
                "material": mname,
                "eps": [complex(eps).real, complex(eps).imag],
                "max_abs_T": float(np.max(np.abs(T))),
                "cond": cond,
                "probes": vals,
            })
            print("mie", rname, mname, "Tmax", rows[-1]["max_abs_T"], flush=True)
    # Frequency sweep at three radii for the production plasma.
    freq_rows = []
    for f_GHz in np.linspace(3.2, 6.5, 34):
        fa = f_GHz * 1e9 * 0.02 / 2.99792458e8
        kk = 2 * np.pi * fa
        ee = plasma(fa)
        for radius in (0.12, R_CORE, 0.40):
            T = mie_T(kk, ee, radius, 10)
            freq_rows.append({"f_GHz": float(f_GHz), "radius": radius, "max_abs_T": float(np.max(np.abs(T))), "eps_re": ee.real})
    return rows, freq_rows


def coated_frequency():
    rows = []
    center = np.zeros((1, 2))
    for f_GHz in (3.50, 3.85, 4.20, 5.30):
        fa = f_GHz * 1e9 * 0.02 / 2.99792458e8
        kk = 2 * np.pi * fa
        ee = plasma(fa)
        coated = ([R_CORE, R_GAP, R_SHELL], [ee, 1.0 + 0j, 3.8 + 0j])
        vals, cond = ratio_map(center, kk, R_CORE, ee, 12, coated=coated)
        rows.append({"f_GHz": f_GHz, "eps": [ee.real, ee.imag], "cond": cond, "probes": vals})
        print("coated f", f_GHz, vals["forward"], flush=True)
    return rows


def cluster_probe_table():
    k0 = 2 * np.pi * float(sc.fs_a)
    eps = plasma(float(sc.fs_a))
    coated = ([R_CORE, R_GAP, R_SHELL], [eps, 1.0 + 0j, 3.8 + 0j])
    ang = np.deg2rad(30.0)
    rot = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    pair = np.array([[0.0, 0.0], [1.0, 0.0]])
    cases = {
        "bare_1": (np.zeros((1, 2)), None),
        "bare_2_axis": (pair, None),
        "bare_2_rot30": (pair @ rot.T, None),
        "bare_3": (_three(), None),
        "bare_7": (_seven(), None),
        "coated_1": (np.zeros((1, 2)), coated),
        "coated_2_axis": (pair, coated),
        "coated_3": (_three(), coated),
        "coated_7": (_seven(), coated),
    }
    out = {}
    for name, (centers, coat) in cases.items():
        vals, cond = ratio_map(centers, k0, R_CORE, eps, 12, coated=coat)
        out[name] = {"cond": cond, "probes": vals, "n": int(len(centers))}
        print("probes", name, vals["forward"], "cond", f"{cond:.3e}", flush=True)
    return out


def contour_power():
    """Net outward power on a circle that encloses the scatterer and not the source."""
    k0 = 2 * np.pi * float(sc.fs_a)
    eps = plasma(float(sc.fs_a))
    rows = []
    specs = [
        ("dielectric", 3.8 + 0j, None),
        ("plasma", eps, None),
        ("coated", eps, ([R_CORE, R_GAP, R_SHELL], [eps, 1.0 + 0j, 3.8 + 0j])),
    ]
    ang = np.linspace(0, 2 * np.pi, 721, endpoint=False)
    radius_c = 1.2
    xy = np.column_stack([radius_c * np.cos(ang), radius_c * np.sin(ang)])
    dtheta = ang[1] - ang[0]
    dl = radius_c * dtheta
    normal = xy / radius_c

    def poynting(field_fn):
        # Central differences on the analytic field, step 1e-5.
        step = 1e-5
        hz = field_fn(xy)
        dx = field_fn(xy + np.array([step, 0.0])) - field_fn(xy - np.array([step, 0.0]))
        dy = field_fn(xy + np.array([0.0, step])) - field_fn(xy - np.array([0.0, step]))
        dhz_dx = dx / (2 * step)
        dhz_dy = dy / (2 * step)
        # Air outside. Ex,Ey = (i/k0) (dHz/dy, -dHz/dx)
        ex = (1j / k0) * dhz_dy
        ey = (1j / k0) * (-dhz_dx)
        sx = 0.5 * np.real(ey * np.conj(hz))
        sy = -0.5 * np.real(ex * np.conj(hz))
        return float(np.sum((sx * normal[:, 0] + sy * normal[:, 1]) * dl))

    for name, ee, coat in specs:
        def fn(pts, ee=ee, coat=coat):
            b = solve_clusters(np.zeros((1, 2)), k0, R_CORE, complex(ee), SRC, 10, coated=coat)
            outer = R_CORE if coat is None else coat[0][-1]
            return cluster_field(pts, np.zeros((1, 2)), k0, outer, complex(ee), SRC, b, 10)

        p_obj = poynting(fn)
        def fn_vac(pts):
            return hankel1(0, k0 * np.linalg.norm(pts - SRC, axis=1))

        p_vac = poynting(fn_vac)
        rows.append({"case": name, "power_object": p_obj, "power_vacuum": p_vac, "absorption_proxy": p_vac - p_obj})
        print("power", rows[-1], flush=True)
    # Lossless quartz slab power from the transfer matrix.
    resp = stack_response(k0, [(0.40, 3.8 + 0j)])
    rows.append({"case": "quartz_slab_0.40", "R": abs(resp["r"]) ** 2, "T": abs(resp["t_over_vacuum"]) ** 2, "sum": abs(resp["r"]) ** 2 + abs(resp["t_over_vacuum"]) ** 2})
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mie_rows, freq_rows = mie_matrix()
    out = {
        "source": SRC.tolist(),
        "truncation": truncation_table(),
        "mie_matrix": mie_rows,
        "mie_frequency": freq_rows,
        "coated_frequency": coated_frequency(),
        "cluster_probes": cluster_probe_table(),
        "contour_power": contour_power(),
    }
    (OUT / "analytic_sweeps.json").write_text(json.dumps(out, indent=2) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    for case in ("bare_1", "bare_7", "coated_7"):
        xs, ys = [], []
        for rec in out["truncation"]:
            if rec["case"] == case and "d_forward" in rec:
                xs.append(rec["nmax"])
                ys.append(rec["d_forward"])
        ax.semilogy(xs, ys, "o-", label=case)
    ax.set_xlabel("m_max")
    ax.set_ylabel("change in forward ratio")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "tmatrix_truncation.png", dpi=140)
    print("ANALYTIC_SWEEPS_DONE", flush=True)


if __name__ == "__main__":
    main()
