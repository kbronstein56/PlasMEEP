#!/usr/bin/env python3
"""Three-way gradient check: FD, forward Jacobian chain, and discrete adjoint.

Real objective J = |Hz(probe)|² on the one-cylinder mesh.
Adjoint: A^H λ = L^H y  with y = L x,  dJ/ds = -2 Re(λ^H (dA/ds) x).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sixport_common as sc  # noqa: E402
from fem_validated_solver import ElementSampler, assemble_anisotropic  # noqa: E402
from gyrotropic_tensor import ordinary_from_si, rho_xy, tensor_ordinary  # noqa: E402
from plasma_sensitivity import element_drho  # noqa: E402
from jacobian_one_cylinder import build_mesh, rho_elements, solve_state  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_full_device_validation"
HS = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5)
REL_TOL = 1e-4
C = 299792458.0


def probe_row(sampler, tris, probe, n_dof):
    """Sparse row L such that L x = Hz(probe) via barycentric interpolation."""
    t = int(sampler.locate(probe[None, :])[0])
    w = sampler._bary(t, probe)
    L = sparse.lil_matrix((1, n_dof), dtype=np.complex128)
    for k in range(3):
        L[0, int(tris[t, k])] = complex(w[k])
    return L.tocsr()


def gate(errs):
    imin = int(np.argmin(errs))
    minimum = errs[imin]
    larger = max(errs[:imin]) if imin > 0 else (max(errs[1:]) if len(errs) > 1 else errs[0])
    return imin, minimum, larger, bool(minimum <= REL_TOL and larger >= 10.0 * minimum)


def run_case(b_tesla):
    pts, tris, inside, center = build_mesh()
    sc.fs_a = sc.fs_Hz * sc.a / C
    k0 = 2 * np.pi * sc.fs_a
    sampler = ElementSampler(pts, tris)
    src = center + np.array([-1.5, 0.2])
    bvec = np.zeros(len(pts), np.complex128)
    t_src = int(sampler.locate(src[None, :])[0])
    w = sampler._bary(t_src, src)
    for k in range(3):
        bvec[int(tris[t_src, k])] += complex(w[k])
    probe = center + np.array([1.3, 0.4])
    L = probe_row(sampler, tris, probe, len(pts))

    x, A, ords, rho = solve_state(pts, tris, inside, 1.0, b_tesla, k0, bvec)
    f, fp, gamma, fc = ords
    fp_ref = sc.fp_Hz * sc.a / C
    dr = element_drho(inside, f, fp_ref, gamma, fc, 1.0)
    dA = assemble_anisotropic(
        pts, tris, *dr, k0, np.zeros(len(pts), dtype=bool), mass_scale=0.0, pin_empty=False
    )
    lu = splu(A.tocsc())
    # forward sensitivity
    dx = lu.solve(-(dA @ x))
    y = complex((L @ x)[0])
    dy_fwd = complex((L @ dx)[0])
    dJ_fwd = 2.0 * np.real(np.conj(y) * dy_fwd)

    # adjoint: A^H λ = L^H y
    # SuperLU of A: solve A^H λ = rhs via lu.solve(rhs, trans='H') if available
    rhs_adj = np.asarray(L.conj().T @ np.array([y])).ravel()
    try:
        lam = lu.solve(rhs_adj, trans="H")
    except TypeError:
        # fallback: factor A^H explicitly
        Ah = A.conj().T.tocsc()
        lam = splu(Ah).solve(rhs_adj)
    dJ_adj = -2.0 * np.real(np.vdot(lam, dA @ x))

    sweep = []
    for h in HS:
        xp, _, _, _ = solve_state(pts, tris, inside, 1.0 + h, b_tesla, k0, bvec)
        xm, _, _, _ = solve_state(pts, tris, inside, 1.0 - h, b_tesla, k0, bvec)
        yp = complex((L @ xp)[0])
        ym = complex((L @ xm)[0])
        dJ_fd = (abs(yp) ** 2 - abs(ym) ** 2) / (2.0 * h)
        rel_fwd = abs(dJ_fd - dJ_fwd) / max(abs(dJ_fwd), 1e-30)
        rel_adj = abs(dJ_fd - dJ_adj) / max(abs(dJ_adj), 1e-30)
        rel_fa = abs(dJ_fwd - dJ_adj) / max(abs(dJ_fwd), 1e-30)
        sweep.append({
            "h": h,
            "fd_vs_fwd": float(rel_fwd),
            "fd_vs_adj": float(rel_adj),
            "fwd_vs_adj": float(rel_fa),
            "dJ_fd": float(dJ_fd),
        })
        print(
            f"B={b_tesla} h={h:.1e} fd-fwd {rel_fwd:.3e} fd-adj {rel_adj:.3e} fwd-adj {rel_fa:.3e}",
            flush=True,
        )

    e_fwd = [r["fd_vs_fwd"] for r in sweep]
    e_adj = [r["fd_vs_adj"] for r in sweep]
    _, m_fwd, l_fwd, p_fwd = gate(e_fwd)
    _, m_adj, l_adj, p_adj = gate(e_adj)
    fwd_adj = float(abs(dJ_fwd - dJ_adj) / max(abs(dJ_fwd), 1e-30))
    return {
        "B_T": b_tesla,
        "y_abs": abs(y),
        "dJ_fwd": float(dJ_fwd),
        "dJ_adj": float(dJ_adj),
        "fwd_vs_adj_rel": fwd_adj,
        "minimum_fd_vs_fwd": m_fwd,
        "minimum_fd_vs_adj": m_adj,
        "pass": bool(p_fwd and p_adj and fwd_adj <= REL_TOL),
        "sweep": sweep,
    }


def main():
    rows = []
    for b in (0.0, 0.05, -0.05):
        row = run_case(b)
        rows.append(row)
        print("CASE", b, "PASS" if row["pass"] else "FAIL", row["fwd_vs_adj_rel"], flush=True)
    out = {"criteria_rel": REL_TOL, "cases": rows, "overall_pass": all(r["pass"] for r in rows)}
    (OUT / "adjoint_validation.json").write_text(json.dumps(out, indent=2) + "\n")
    print("WROTE adjoint_validation.json", out["overall_pass"], flush=True)


if __name__ == "__main__":
    main()
