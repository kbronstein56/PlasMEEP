#!/usr/bin/env python3
"""Detailed per-frequency kernel profile for Stage-5 planning (no optimization).

Breaks build_system + A2 eval/grad into wall-clock phases on Stage-4 q_best.
Does not change physics. Writes STAGE5_KERNEL_BASELINE_PROFILE.json.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
sys.path[:0] = [
    str(Path(__file__).resolve().parent),
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(k, "1")

import sixport_common as sc  # noqa: E402
from adjoint_direct_dtn import build_system  # noqa: E402
from a2_efficient_core import eval_and_grad_A2, rss_gib  # noqa: E402
from boundary_ports import build_port_from_boundary  # noqa: E402
from diff_fem_core import assemble_dA_rod, per_rod_plasma_masks, rho_from_scales  # noqa: E402
from direct_dtn_s import assemble_port_operator, build_Q_columns  # noqa: E402
from dtn_port_geometry import port_outer_center  # noqa: E402
from efficient_adjoint import S_from_U, _pack_C, _solve_multi  # noqa: E402
from fem_validated_solver import assemble_anisotropic  # noqa: E402
from full91_dtn_solve import exterior_beyond_ports_mask  # noqa: E402
from full91_geometry import bulb_centers_device  # noqa: E402
from full91_solve import load_mesh, material_masks  # noqa: E402
from modes import mode_constants  # noqa: E402
from multimode_opt_core import B_PLUS, expand_tied, orbits_and_masks, te0_metrics_from_S  # noqa: E402
from objective_A2 import EPS, objective_A2_from_S, scalar_adjoint_rhs_A2  # noqa: E402


def timed(fn, *a, **kw):
    t0 = time.perf_counter()
    out = fn(*a, **kw)
    return out, time.perf_counter() - t0


def profile_build_system(grade, s_vec, f_hz, b_tesla=B_PLUS, measure="projected"):
    tim = {}
    t_all = time.perf_counter()

    mesh, tim["load_mesh"] = timed(load_mesh, grade)
    points, tris = mesh
    masks, tim["material_masks"] = timed(material_masks, points, tris, sc)
    plasma, quartz, walls, _ = masks
    centers = bulb_centers_device(sc)
    r_p = 4.6 * sc.r_bulb_inner / 6.5
    rod_masks, tim["rod_masks"] = timed(per_rod_plasma_masks, points, tris, centers, r_p)
    rho_pack, tim["rho_from_scales"] = timed(
        rho_from_scales,
        len(tris),
        rod_masks,
        s_vec,
        quartz,
        walls,
        f_hz,
        b_tesla,
        sc.a,
        sc.fp_Hz,
        sc.gamma_Hz,
    )
    rho, f_ord, gamma, fc = rho_pack
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    rough = []
    for p in range(6):
        n_hat = np.asarray(sc.port_dirs[p], float)
        n_hat /= np.linalg.norm(n_hat)
        c = port_outer_center(sc.full_horns[p], n_hat) + shift
        rough.append({"center": c, "outward": n_hat})
    ext, tim["exterior_mask"] = timed(exterior_beyond_ports_mask, points, tris, rough)
    for arr in rho:
        arr[ext] = 0.0
    k0 = 2.0 * np.pi * (f_hz * sc.a / 299792458.0)
    pec0 = np.zeros(len(points), dtype=bool)
    A0, tim["assemble_bulk"] = timed(assemble_anisotropic, points, tris, *rho, k0, pec0)
    ports = []
    t_ports = 0.0
    for p in range(6):
        n_hat = np.asarray(sc.port_dirs[p], float)
        n_hat /= np.linalg.norm(n_hat)
        tang = np.array([-n_hat[1], n_hat[0]])
        c = port_outer_center(sc.full_horns[p], n_hat) + shift
        t0 = time.perf_counter()
        ps = build_port_from_boundary(
            points, tris, rho, c, n_hat, tang, W=float(sc.clear_width), plane_tol=0.12, port=p
        )
        t_ports += time.perf_counter() - t0
        ps["measure"] = measure
        ports.append(ps)
    tim["build_ports"] = t_ports
    modes = mode_constants(k0, float(sc.clear_width), n_modes=4)
    A, tim["assemble_port_operator"] = timed(
        assemble_port_operator, A0, points, ports, modes, n_prop=2, n_evan=2
    )
    tim["build_system_total"] = time.perf_counter() - t_all
    st = {
        "points": points,
        "tris": tris,
        "rho": rho,
        "A": A,
        "k0": k0,
        "ports": ports,
        "modes": modes,
        "rod_masks": rod_masks,
        "f_ord": f_ord,
        "gamma": gamma,
        "fc": fc,
        "fp_ref_ord": sc.fp_Hz * sc.a / 299792458.0,
    }
    return st, tim


def profile_eval(grade, q, orbits, f_hz):
    tim = {}
    rss0 = rss_gib()
    t_all = time.perf_counter()
    s = expand_tied(np.asarray(q, float), orbits)
    st, tim_bs = profile_build_system(grade, s, f_hz)
    for k, v in tim_bs.items():
        tim[f"bs_{k}"] = v

    t0 = time.perf_counter()
    _, _, o_masks = orbits_and_masks(rod_masks=st["rod_masks"])
    points, A = st["points"], st["A"]
    cols = build_Q_columns(points, st["ports"], st["modes"], n_prop=2)
    C, betas = _pack_C(points, cols)
    te0_cols = [2 * j for j in range(6)]
    tim["pack_ports"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    lu = splu(A.tocsc())
    tim["factor_lu"] = time.perf_counter() - t0
    tim["rss_after_factor"] = rss_gib()

    t0 = time.perf_counter()
    U = _solve_multi(lu, C[:, te0_cols], trans="N")
    tim["forward_solve"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    S = S_from_U(C, U, betas, te0_cols)
    J = objective_A2_from_S(S, eps=EPS)
    met = te0_metrics_from_S(S)
    tim["S_and_objective"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    R = np.column_stack(
        [scalar_adjoint_rhs_A2(C, betas, S, j, eps=EPS) for j in range(6)]
    )
    Lam = _solve_multi(lu, R, trans="T")
    tim["adjoint_solve"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    g = np.zeros(len(orbits), float)
    dA_times = []
    contract_times = []
    for m, o in enumerate(orbits):
        s_m = float(s[o["rods"][0]])
        t1 = time.perf_counter()
        dA = assemble_dA_rod(
            points,
            st["tris"],
            o_masks[m],
            st["f_ord"],
            st["fp_ref_ord"],
            st["gamma"],
            st["fc"],
            st["k0"],
            s=s_m,
        )
        dA_times.append(time.perf_counter() - t1)
        t1 = time.perf_counter()
        acc = 0.0
        for j_loc in range(6):
            dAu = dA @ U[:, j_loc]
            acc += 2.0 * np.real(np.dot(Lam[:, j_loc], dAu))
        g[m] = float(acc / 6.0)
        contract_times.append(time.perf_counter() - t1)
    tim["dA_assemble_sum"] = float(sum(dA_times))
    tim["dA_matvec_sum"] = float(sum(contract_times))
    tim["dA_contraction_total"] = time.perf_counter() - t0
    tim["dA_assemble_mean"] = float(np.mean(dA_times))
    tim["dA_matvec_mean"] = float(np.mean(contract_times))
    tim["total"] = time.perf_counter() - t_all
    tim["rss_start"] = rss0
    tim["rss_end"] = rss_gib()
    tim["n_dof"] = int(A.shape[0])
    tim["n_tris"] = int(len(st["tris"]))
    tim["orbit_mask_sizes"] = [int(m.sum()) for m in o_masks]
    return {"J": float(J), "g": np.asarray(g, float), "metrics": met, "timings": tim}


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument(
        "--freqs",
        default="4.5,4.73125,5.08125,5.8125",
        help="GHz list",
    )
    args = ap.parse_args()
    q = np.asarray(
        json.loads((OUT / "opt_broadband_stage4_A2_C.json").read_text())["q_best"],
        float,
    )
    orbits, _, _ = orbits_and_masks()
    freqs = [float(x) for x in args.freqs.split(",")]
    report = {
        "grade": args.grade,
        "q_source": "opt_broadband_stage4_A2_C.json",
        "omp": os.environ.get("OMP_NUM_THREADS"),
        "freqs_GHz": freqs,
        "profiles": [],
    }
    # Warmup once on the REFERENCE kernel (this script profiles the legacy path)
    print("warmup @ 5.0 GHz via reference eval_and_grad_A2", flush=True)
    eval_and_grad_A2(args.grade, q, orbits, 5.0e9, profile=True, use_fast_kernel=False)
    for f in freqs:
        print(f"=== detailed profile @ {f} GHz ===", flush=True)
        rec = profile_eval(args.grade, q, orbits, f * 1e9)
        tim = rec["timings"]
        total = tim["total"]
        shares = {
            k: (100.0 * v / total if isinstance(v, float) and k.startswith(("bs_", "dA_", "factor", "forward", "adjoint", "pack", "S_")) else None)
            for k, v in tim.items()
        }
        print(
            f"  total={total:.2f}s  assemble_bulk={tim['bs_assemble_bulk']:.2f}s "
            f"dA_asm={tim['dA_assemble_sum']:.2f}s factor={tim['factor_lu']:.2f}s "
            f"ports={tim['bs_build_ports']:.2f}s",
            flush=True,
        )
        report["profiles"].append(
            {
                "f_GHz": f,
                "J": rec["J"],
                "grad_norm": float(np.linalg.norm(rec["g"])),
                "timings": tim,
                "share_pct_of_total": {
                    "assemble_bulk": 100.0 * tim["bs_assemble_bulk"] / total,
                    "build_ports": 100.0 * tim["bs_build_ports"] / total,
                    "assemble_port_op": 100.0 * tim["bs_assemble_port_operator"] / total,
                    "load_masks_rho": 100.0
                    * (
                        tim["bs_load_mesh"]
                        + tim["bs_material_masks"]
                        + tim["bs_rod_masks"]
                        + tim["bs_rho_from_scales"]
                        + tim["bs_exterior_mask"]
                    )
                    / total,
                    "factor_lu": 100.0 * tim["factor_lu"] / total,
                    "forward": 100.0 * tim["forward_solve"] / total,
                    "adjoint": 100.0 * tim["adjoint_solve"] / total,
                    "dA_assemble": 100.0 * tim["dA_assemble_sum"] / total,
                    "dA_matvec": 100.0 * tim["dA_matvec_sum"] / total,
                    "S_obj": 100.0 * tim["S_and_objective"] / total,
                },
            }
        )
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "STAGE5_KERNEL_BASELINE_PROFILE.json"
    path.write_text(json.dumps(report, indent=2))
    print("wrote", path, flush=True)


if __name__ == "__main__":
    main()
