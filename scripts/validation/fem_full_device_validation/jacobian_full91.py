#!/usr/bin/env python3
"""Full-91 material Jacobian: selected-rod FD checks + optional full 91-column port J.

One factorization of A at production s=1. Sensitivity columns reuse that LU.
Finite differences reassemble/refactor only for the independent check.

Usage:
  jacobian_full91.py [grade] [--fd-rods] [--full-J] [--B 0.0]
Defaults: grade=F, B=0, both --fd-rods and a reduced full-J for one source.
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sixport_common as sc  # noqa: E402
from diff_fem_core import (  # noqa: E402
    monitor_amp,
    per_rod_plasma_masks,
    sensitivity_column,
    solve_forward,
)
from fem_validated_solver import ElementSampler, assemble_anisotropic, guide_normal_flux  # noqa: E402
from full91_geometry import OUT, bulb_centers_device, production_constants  # noqa: E402
from full91_solve import inject, load_mesh, material_masks, port_lines, free_gib  # noqa: E402

HS = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4)
REL_TOL = 1e-4


def peak_rss_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0**2)


def select_rods(centers: np.ndarray) -> dict:
    """Representative rods: center, inner, outer, near horn 0, symmetry pair."""
    c0 = centers - centers.mean(0)
    r = np.linalg.norm(c0, axis=1)
    # center: smallest r
    i_center = int(np.argmin(r))
    # inner ring ~ first shell
    shells = np.unique(np.round(r, 3))
    shells = shells[shells > 0.05]
    r_inner = shells[0] if len(shells) else None
    r_outer = shells[-1] if len(shells) else None
    i_inner = int(np.argmin(np.abs(r - r_inner))) if r_inner is not None else i_center
    i_outer = int(np.argmin(np.abs(r - r_outer))) if r_outer is not None else i_center
    # near port-0 horn: maximize x (port 0 is +x in device coords typically)
    # use distance to port-0 monitor
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    mon0 = np.asarray(sc.full_horns[0]["monitor_center"], float) + shift
    d_horn = np.linalg.norm(centers - mon0, axis=1)
    i_horn = int(np.argmin(d_horn))
    # symmetry pair about x-axis: rod with y>0 and its mirror
    # pick an outer-ish rod with y > 0.2
    cand = np.where((c0[:, 1] > 0.3) & (r > 0.5 * r.max()))[0]
    if len(cand) == 0:
        cand = np.where(c0[:, 1] > 0.2)[0]
    i_sym_a = int(cand[int(np.argmin(np.abs(r[cand] - np.median(r[cand]))))]) if len(cand) else i_outer
    # mirror: closest to (x, -y)
    target = centers[i_sym_a] * np.array([1.0, -1.0]) + np.array([0.0, 2 * centers.mean(0)[1]])
    # device y mirror about domain center
    cy = sc.ny_ports / 2.0
    mirror = np.array([centers[i_sym_a, 0], 2 * cy - centers[i_sym_a, 1]])
    i_sym_b = int(np.argmin(np.linalg.norm(centers - mirror, axis=1)))
    return {
        "center": i_center,
        "inner": i_inner,
        "outer": i_outer,
        "near_horn0": i_horn,
        "sym_a": i_sym_a,
        "sym_b": i_sym_b,
    }


def gate(errs):
    imin = int(np.argmin(errs))
    minimum = errs[imin]
    larger = max(errs[:imin]) if imin > 0 else (max(errs[1:]) if len(errs) > 1 else errs[0])
    return imin, minimum, larger, bool(minimum <= REL_TOL and larger >= 10.0 * minimum)


def rho_walls_quartz_only(n, quartz, walls):
    """Helper unused — scales go through solve_forward."""
    return None


def fd_rod_check(grade, b_tesla, rod_ids, labels, source_port=0):
    _sc, const = production_constants()
    points, tris = load_mesh(grade)
    masks_all = material_masks(points, tris, sc)
    plasma, quartz, walls, centers = masks_all
    centers = np.asarray(centers, float)
    r_p = const["r_plasma_material_a"]
    rod_masks = per_rod_plasma_masks(points, tris, centers, r_p)
    union = np.zeros(len(tris), dtype=bool)
    for m in rod_masks:
        union |= m
    print(
        f"grade={grade} B={b_tesla} nodes={len(points)} rods={len(rod_masks)} "
        f"plasma_union={int(union.sum())} plasma_mask={int(plasma.sum())} free={free_gib():.0f}GiB",
        flush=True,
    )
    sampler = ElementSampler(points, tris)
    lines = port_lines()
    bvec = inject(points, tris, sampler, lines[source_port]["source_xy"], lines[source_port]["weights"])
    s0 = np.ones(len(rod_masks))
    t0 = time.perf_counter()
    state = solve_forward(
        points, tris, rod_masks, quartz, walls, s0, bvec,
        sc.fs_Hz, b_tesla, sc.a, sc.fp_Hz, sc.gamma_Hz,
    )
    t_fac = time.perf_counter() - t0
    print(f"base factor+solve {t_fac:.1f}s peak_rss={peak_rss_gib():.1f}GiB", flush=True)

    # receiving ports: all except source
    recv = [i for i in range(6) if i != source_port]
    # direct port-amp derivatives for selected rods
    results = []
    for lab, k in zip(labels, rod_ids):
        t1 = time.perf_counter()
        dx = sensitivity_column(state, k)
        t_rhs = time.perf_counter() - t1
        da = {}
        for i in recv:
            da[i] = monitor_amp(
                sampler, dx, state.rho, state.k0, lines[i]["monitor_xy"], lines[i]["weights"]
            )
        # also d of power on one receive port via finite flux derivative from dx
        # (rho vacuum at monitor)
        i_pow = recv[0]
        p0 = guide_normal_flux(
            points, state.x, None, None, i_pow, sampler=sampler, tris=tris, rho=state.rho, omega=state.k0
        )
        # FD sweep
        sweep = []
        for h in HS:
            sp = s0.copy(); sp[k] = 1.0 + h
            sm = s0.copy(); sm[k] = 1.0 - h
            t2 = time.perf_counter()
            st_p = solve_forward(
                points, tris, rod_masks, quartz, walls, sp, bvec,
                sc.fs_Hz, b_tesla, sc.a, sc.fp_Hz, sc.gamma_Hz,
            )
            st_m = solve_forward(
                points, tris, rod_masks, quartz, walls, sm, bvec,
                sc.fs_Hz, b_tesla, sc.a, sc.fp_Hz, sc.gamma_Hz,
            )
            t_fd = time.perf_counter() - t2
            dx_fd = (st_p.x - st_m.x) / (2.0 * h)
            field_rel = float(np.linalg.norm(dx_fd - dx) / max(np.linalg.norm(dx), 1e-30))
            amp_rels = []
            for i in recv:
                ap = monitor_amp(sampler, st_p.x, st_p.rho, st_p.k0, lines[i]["monitor_xy"], lines[i]["weights"])
                am = monitor_amp(sampler, st_m.x, st_m.rho, st_m.k0, lines[i]["monitor_xy"], lines[i]["weights"])
                da_fd = (ap - am) / (2.0 * h)
                amp_rels.append(abs(da_fd - da[i]) / max(abs(da[i]), 1e-30))
            sweep.append({
                "h": h,
                "field_rel": field_rel,
                "amp_rel_max": float(max(amp_rels)),
                "amp_rel_by_port": {str(i): float(r) for i, r in zip(recv, amp_rels)},
                "fd_pair_s": t_fd,
            })
            print(
                f"rod {lab}({k}) h={h:.1e} field {field_rel:.3e} amp_max {max(amp_rels):.3e} "
                f"fd_pair={t_fd:.1f}s",
                flush=True,
            )
            # free LU memory from FD states
            del st_p, st_m
        ferrs = [r["field_rel"] for r in sweep]
        aerrs = [r["amp_rel_max"] for r in sweep]
        fi, fmin, flarger, fpass = gate(ferrs)
        ai, amin, alarger, apass = gate(aerrs)
        row = {
            "label": lab,
            "rod_index": int(k),
            "center_xy": centers[k].tolist(),
            "plasma_tris": int(rod_masks[k].sum()),
            "rhs_s": t_rhs,
            "p0_recv": float(p0),
            "da_abs": {str(i): abs(da[i]) for i in recv},
            "minimum_field_rel": fmin,
            "minimum_amp_rel": amin,
            "h_field": sweep[fi]["h"],
            "h_amp": sweep[ai]["h"],
            "field_pass": bool(fpass),
            "amp_pass": bool(apass),
            "pass": bool(fpass and apass),
            "sweep": sweep,
        }
        results.append(row)
        print("ROD", lab, "PASS" if row["pass"] else "FAIL", "field", fmin, "amp", amin, flush=True)
    return {
        "grade": grade,
        "B_T": b_tesla,
        "source_port": source_port,
        "nodes": int(len(points)),
        "n_rods": len(rod_masks),
        "base_factor_solve_s": t_fac,
        "peak_rss_gib": peak_rss_gib(),
        "rods": results,
        "pass": all(r["pass"] for r in results),
    }


def full_port_jacobian(grade, b_tesla, source_port=0):
    """Build J (n_recv_ports × 91) for complex monitor amplitudes. One LU."""
    _sc, const = production_constants()
    points, tris = load_mesh(grade)
    plasma, quartz, walls, centers = material_masks(points, tris, sc)
    centers = np.asarray(centers, float)
    rod_masks = per_rod_plasma_masks(points, tris, centers, const["r_plasma_material_a"])
    sampler = ElementSampler(points, tris)
    lines = port_lines()
    bvec = inject(points, tris, sampler, lines[source_port]["source_xy"], lines[source_port]["weights"])
    s0 = np.ones(len(rod_masks))
    t0 = time.perf_counter()
    state = solve_forward(
        points, tris, rod_masks, quartz, walls, s0, bvec,
        sc.fs_Hz, b_tesla, sc.a, sc.fp_Hz, sc.gamma_Hz,
    )
    t_fac = time.perf_counter() - t0
    recv = [i for i in range(6) if i != source_port]
    J = np.zeros((len(recv), len(rod_masks)), np.complex128)
    t_rhs = 0.0
    for k in range(len(rod_masks)):
        t1 = time.perf_counter()
        dx = sensitivity_column(state, k)
        t_rhs += time.perf_counter() - t1
        for m, i in enumerate(recv):
            J[m, k] = monitor_amp(
                sampler, dx, state.rho, state.k0, lines[i]["monitor_xy"], lines[i]["weights"]
            )
        if (k + 1) % 10 == 0 or k == 0:
            print(f"J col {k+1}/{len(rod_masks)} cum_rhs={t_rhs:.1f}s rss={peak_rss_gib():.1f}GiB", flush=True)
    out = {
        "grade": grade,
        "B_T": b_tesla,
        "source_port": source_port,
        "recv_ports": recv,
        "shape": [int(J.shape[0]), int(J.shape[1])],
        "factor_solve_s": t_fac,
        "all_rhs_s": t_rhs,
        "s_per_column": t_rhs / max(len(rod_masks), 1),
        "total_jacobian_s": t_fac + t_rhs,
        "peak_rss_gib": peak_rss_gib(),
        "J_abs_mean": float(np.mean(np.abs(J))),
        "J_abs_max": float(np.max(np.abs(J))),
        # store compact: real/imag as lists
        "J_re": J.real.tolist(),
        "J_im": J.imag.tolist(),
    }
    tag = f"full91_jacobian_{grade}_B{b_tesla:+.4f}_src{source_port}".replace("+", "p").replace("-", "m")
    (OUT / f"{tag}.json").write_text(json.dumps(out) + "\n")
    print("WROTE", tag, "total_s", out["total_jacobian_s"], flush=True)
    return out


def main():
    argv = sys.argv[1:]
    grade = "F"
    b_tesla = 0.0
    do_fd = "--fd-rods" in argv or "--full-J" not in argv
    do_full = "--full-J" in argv
    if not argv:
        do_fd = True
        do_full = True
    for a in argv:
        if a.startswith("--B"):
            b_tesla = float(a.split("=")[1]) if "=" in a else float(argv[argv.index(a) + 1])
        elif a in ("F", "M", "C", "A", "Q") or (not a.startswith("--") and a[0].isalpha()):
            if a not in ("--fd-rods", "--full-J"):
                grade = a
    # positional grade
    if argv and not argv[0].startswith("--"):
        grade = argv[0]

    points, tris = load_mesh(grade)
    _, _, _, centers = material_masks(points, tris, sc)
    centers = np.asarray(centers, float)
    sel = select_rods(centers)
    print("selected rods", sel, flush=True)
    # unique ordered selection for FD
    labels, ids = [], []
    for lab in ("center", "inner", "outer", "near_horn0", "sym_a", "sym_b"):
        k = sel[lab]
        if k not in ids:
            labels.append(lab)
            ids.append(k)
        else:
            labels.append(f"{lab}_dup{k}")
            # still record but skip duplicate FD
    # dedupe for FD work
    seen = set()
    labs_u, ids_u = [], []
    for lab, k in zip(
        ("center", "inner", "outer", "near_horn0", "sym_a", "sym_b"),
        (sel["center"], sel["inner"], sel["outer"], sel["near_horn0"], sel["sym_a"], sel["sym_b"]),
    ):
        if k in seen:
            continue
        seen.add(k)
        labs_u.append(lab)
        ids_u.append(k)

    summary = {"selection": {k: int(v) for k, v in sel.items()}, "grade": grade, "B_T": b_tesla}
    tag = f"{grade}_B{b_tesla:+.4f}".replace("+", "p").replace("-", "m")
    if do_fd:
        fd = fd_rod_check(grade, b_tesla, ids_u, labs_u)
        summary["fd_rods"] = fd
        (OUT / f"full91_jacobian_fd_selected_{tag}.json").write_text(json.dumps(fd, indent=2) + "\n")
        (OUT / "full91_jacobian_fd_selected.json").write_text(json.dumps(fd, indent=2) + "\n")
        print("FD overall", fd["pass"], flush=True)
    if do_full:
        fj = full_port_jacobian(grade, b_tesla)
        summary["full_J"] = {
            "shape": fj["shape"],
            "factor_solve_s": fj["factor_solve_s"],
            "all_rhs_s": fj["all_rhs_s"],
            "s_per_column": fj["s_per_column"],
            "total_jacobian_s": fj["total_jacobian_s"],
            "peak_rss_gib": fj["peak_rss_gib"],
            "J_abs_mean": fj["J_abs_mean"],
            "J_abs_max": fj["J_abs_max"],
        }
    (OUT / f"full91_jacobian_summary_{tag}.json").write_text(json.dumps(summary, indent=2) + "\n")
    (OUT / "full91_jacobian_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("WROTE full91_jacobian_summary.json", flush=True)


if __name__ == "__main__":
    main()
