#!/usr/bin/env python3
"""Full91 adjoint spot-check: one scalar FoM, a few rods, vs forward+FD.

J = |monitor_amplitude on port 3|^2 for source port 0.
Grade default M for speed; pass F via argv.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sixport_common as sc  # noqa: E402
from diff_fem_core import assemble_dA_rod, per_rod_plasma_masks, sensitivity_column, solve_forward  # noqa: E402
from fem_validated_solver import ElementSampler  # noqa: E402
from full91_geometry import production_constants  # noqa: E402
from full91_solve import inject, load_mesh, material_masks, port_lines  # noqa: E402
from jacobian_full91 import select_rods, gate, HS, REL_TOL  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_full_device_validation"


def monitor_L(sampler, tris, xy, weights, n_dof):
    """Assemble L such that L x ≈ weighted Hz samples (barycentric)."""
    L = sparse.lil_matrix((1, n_dof), dtype=np.complex128)
    tid = sampler.locate(xy)
    for n in range(len(xy)):
        t = int(tid[n])
        if t < 0:
            continue
        w = sampler._bary(t, xy[n])
        for k in range(3):
            L[0, int(tris[t, k])] += complex(w[k]) * float(weights[n])
    return L.tocsr()


def main():
    grade = sys.argv[1] if len(sys.argv) > 1 else "M"
    b_tesla = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    _sc, const = production_constants()
    points, tris = load_mesh(grade)
    plasma, quartz, walls, centers = material_masks(points, tris, sc)
    centers = np.asarray(centers, float)
    rod_masks = per_rod_plasma_masks(points, tris, centers, const["r_plasma_material_a"])
    sel = select_rods(centers)
    rods = [sel["center"], sel["near_horn0"], sel["outer"]]
    labels = ["center", "near_horn0", "outer"]
    sampler = ElementSampler(points, tris)
    lines = port_lines()
    bvec = inject(points, tris, sampler, lines[0]["source_xy"], lines[0]["weights"])
    L = monitor_L(sampler, tris, lines[3]["monitor_xy"], lines[3]["weights"], len(points))
    s0 = np.ones(len(rod_masks))
    t0 = time.perf_counter()
    state = solve_forward(
        points, tris, rod_masks, quartz, walls, s0, bvec,
        sc.fs_Hz, b_tesla, sc.a, sc.fp_Hz, sc.gamma_Hz,
    )
    print(f"factor {time.perf_counter()-t0:.1f}s grade={grade} B={b_tesla}", flush=True)
    y = complex((L @ state.x)[0])
    rhs = np.asarray(L.conj().T @ np.array([y])).ravel()
    try:
        lam = state.lu.solve(rhs, trans="H")
    except TypeError:
        from scipy.sparse.linalg import splu
        lam = splu(state.A.conj().T.tocsc()).solve(rhs)

    rows = []
    for lab, k in zip(labels, rods):
        dA = assemble_dA_rod(
            points, tris, rod_masks[k], state.f_ord, state.fp_ref_ord,
            state.gamma, state.fc, state.k0, 1.0,
        )
        dx = state.lu.solve(-(dA @ state.x))
        dJ_fwd = 2.0 * np.real(np.conj(y) * complex((L @ dx)[0]))
        dJ_adj = -2.0 * np.real(np.vdot(lam, dA @ state.x))
        fwd_adj = abs(dJ_fwd - dJ_adj) / max(abs(dJ_fwd), 1e-30)
        sweep = []
        for h in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3):
            sp = s0.copy(); sp[k] = 1.0 + h
            sm = s0.copy(); sm[k] = 1.0 - h
            st_p = solve_forward(
                points, tris, rod_masks, quartz, walls, sp, bvec,
                sc.fs_Hz, b_tesla, sc.a, sc.fp_Hz, sc.gamma_Hz,
            )
            st_m = solve_forward(
                points, tris, rod_masks, quartz, walls, sm, bvec,
                sc.fs_Hz, b_tesla, sc.a, sc.fp_Hz, sc.gamma_Hz,
            )
            Jp = abs(complex((L @ st_p.x)[0])) ** 2
            Jm = abs(complex((L @ st_m.x)[0])) ** 2
            dJ_fd = (Jp - Jm) / (2.0 * h)
            rel = abs(dJ_fd - dJ_adj) / max(abs(dJ_adj), 1e-30)
            sweep.append({"h": h, "fd_vs_adj": float(rel), "dJ_fd": float(dJ_fd)})
            print(f"{lab} h={h:.1e} fd-adj {rel:.3e} fwd-adj {fwd_adj:.3e}", flush=True)
            del st_p, st_m
        errs = [r["fd_vs_adj"] for r in sweep]
        _, minimum, larger, passed = gate(errs)
        rows.append({
            "label": lab,
            "rod": int(k),
            "dJ_fwd": float(dJ_fwd),
            "dJ_adj": float(dJ_adj),
            "fwd_vs_adj": float(fwd_adj),
            "minimum_fd_vs_adj": minimum,
            "pass": bool(passed and fwd_adj <= REL_TOL),
            "sweep": sweep,
        })
        print("ROD", lab, "PASS" if rows[-1]["pass"] else "FAIL", flush=True)
    out = {
        "grade": grade,
        "B_T": b_tesla,
        "y_abs": abs(y),
        "rods": rows,
        "overall_pass": all(r["pass"] for r in rows),
    }
    tag = f"adjoint_full91_spot_{grade}_B{b_tesla:+.4f}".replace("+", "p").replace("-", "m")
    (OUT / f"{tag}.json").write_text(json.dumps(out, indent=2) + "\n")
    print("WROTE", tag, out["overall_pass"], flush=True)


if __name__ == "__main__":
    main()
