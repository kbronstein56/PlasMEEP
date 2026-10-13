#!/usr/bin/env python3
"""Phase 8: direct-DtN multimode objective gradient vs centered FD.

Objective (TE0 incidence only, average over 6 source ports):

  J = W_des * T_des_TE0
    - W_rev * T_rev_TE0
    - W_thr * T_thr_TE0
    - W_leak * T_other_TE0
    - W_refl * R_TE0
    - W_te1 * Σ_ports |S_TE1←TE0|²

S = I + 2i Q, Q_mn = √(β_m β_n) c_mᵀ A⁻¹ c_n.
Port loads c independent of rod scales s; ∂A/∂s from plasma only.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "multimode_ports"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from boundary_ports import modal_c_vector  # noqa: E402
from diff_fem_core import assemble_dA_rod, per_rod_plasma_masks, rho_from_scales  # noqa: E402
from direct_dtn_s import assemble_port_operator, build_Q_columns, direct_S_from_factor  # noqa: E402
from dtn_port_geometry import port_outer_center  # noqa: E402
from fem_validated_solver import assemble_anisotropic  # noqa: E402
from full91_direct_dtn import load_s  # noqa: E402
from full91_dtn_solve import exterior_beyond_ports_mask  # noqa: E402
from full91_geometry import bulb_centers_device  # noqa: E402
from full91_solve import load_mesh, material_masks  # noqa: E402
from modes import mode_constants  # noqa: E402

W_DES, W_REV, W_THR, W_LEAK, W_REFL, W_TE1 = 1.0, 2.0, 1.5, 1.0, 0.25, 1.0
FD_REL_BAR = 1e-3
# Match circulator_core design box (s may exceed 1)
S_LO, S_HI = 0.05, 1.80


def build_system(grade, s_vec, f_hz, b_tesla, measure="projected"):
    points, tris = load_mesh(grade)
    masks = material_masks(points, tris, sc)
    plasma, quartz, walls, _ = masks
    centers = bulb_centers_device(sc)
    r_p = 4.6 * sc.r_bulb_inner / 6.5
    rod_masks = per_rod_plasma_masks(points, tris, centers, r_p)
    rho, f_ord, gamma, fc = rho_from_scales(
        len(tris), rod_masks, s_vec, quartz, walls, f_hz, b_tesla, sc.a, sc.fp_Hz, sc.gamma_Hz
    )
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    rough = []
    for p in range(6):
        n_hat = np.asarray(sc.port_dirs[p], float)
        n_hat /= np.linalg.norm(n_hat)
        c = port_outer_center(sc.full_horns[p], n_hat) + shift
        rough.append({"center": c, "outward": n_hat})
    ext = exterior_beyond_ports_mask(points, tris, rough)
    for arr in rho:
        arr[ext] = 0.0
    k0 = 2.0 * np.pi * (f_hz * sc.a / 299792458.0)
    A0 = assemble_anisotropic(points, tris, *rho, k0, np.zeros(len(points), dtype=bool))
    ports = []
    for p in range(6):
        n_hat = np.asarray(sc.port_dirs[p], float)
        n_hat /= np.linalg.norm(n_hat)
        tang = np.array([-n_hat[1], n_hat[0]])
        c = port_outer_center(sc.full_horns[p], n_hat) + shift
        from boundary_ports import build_port_from_boundary

        ps = build_port_from_boundary(
            points, tris, rho, c, n_hat, tang, W=float(sc.clear_width), plane_tol=0.12, port=p
        )
        ps["measure"] = measure
        ports.append(ps)
    modes = mode_constants(k0, float(sc.clear_width), n_modes=4)
    A = assemble_port_operator(A0, points, ports, modes, n_prop=2, n_evan=2)
    return {
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


def objective_from_S(S: np.ndarray) -> float:
    """TE0-drive multimode circulator objective (mean over 6 sources)."""
    js = []
    for j in range(6):
        col = 2 * j
        des, rev, thr = (j + 1) % 6, (j - 1) % 6, (j + 3) % 6
        Tdes = abs(S[2 * des, col]) ** 2
        Trev = abs(S[2 * rev, col]) ** 2
        Tthr = abs(S[2 * thr, col]) ** 2
        R = abs(S[col, col]) ** 2
        other = 0.0
        for p in range(6):
            if p in (j, des, rev, thr):
                continue
            other += abs(S[2 * p, col]) ** 2
        te1 = float(sum(abs(S[2 * i + 1, col]) ** 2 for i in range(6)))
        js.append(
            W_DES * Tdes
            - W_REV * Trev
            - W_THR * Tthr
            - W_LEAK * other
            - W_REFL * R
            - W_TE1 * te1
        )
    return float(np.mean(js))


def channel_weights(j_src: int) -> np.ndarray:
    """w[m] for |S_m, col|² terms; col = 2*j_src (TE0 drive)."""
    w = np.zeros(12, float)
    des, rev, thr = (j_src + 1) % 6, (j_src - 1) % 6, (j_src + 3) % 6
    w[2 * des] = W_DES
    w[2 * rev] = -W_REV
    w[2 * thr] = -W_THR
    w[2 * j_src] = -W_REFL
    for p in range(6):
        if p not in (j_src, des, rev, thr):
            w[2 * p] = -W_LEAK
        w[2 * p + 1] = -W_TE1
    return w


def _channel_basis(st):
    points, A, ports, modes = st["points"], st["A"], st["ports"], st["modes"]
    lu = splu(A.tocsc())
    cols = build_Q_columns(points, ports, modes, n_prop=2)
    n_dof = len(points)
    n_ch = 12
    C = np.zeros((n_dof, n_ch), dtype=np.complex128)
    betas = np.zeros(n_ch, dtype=np.complex128)
    for j, (idx, c, beta) in enumerate(cols):
        betas[j] = beta
        if len(idx):
            C[idx, j] = c
    U = np.zeros_like(C)
    for j in range(n_ch):
        U[:, j] = lu.solve(C[:, j])
    luT = splu(A.T.tocsc())
    V = np.zeros_like(C)
    for m in range(n_ch):
        V[:, m] = luT.solve(C[:, m])
    S = direct_S_from_factor(lu, n_dof, cols)
    return lu, cols, C, U, V, betas, S, objective_from_S(S)


def _grad_from_dS(S, betas, dS_fn) -> float:
    """Mean over TE0 sources of Σ_m w_m · 2 Re(conj(S) dS)."""
    g = 0.0
    n_ch = 12
    for j_src in range(6):
        n = 2 * j_src
        w = channel_weights(j_src)
        for m in range(n_ch):
            if abs(w[m]) < 1e-30:
                continue
            dS = dS_fn(m, n)
            g += w[m] * 2.0 * np.real(np.conj(S[m, n]) * dS)
    return float(g / 6.0)


def adjoint_gradient(st, s_vec, rod_indices):
    """∂J/∂s_k for selected rods via Aᵀ adjoints.

    Also returns forward-sensitivity grads (A du = −dA u) for cross-check.
    """
    points, A = st["points"], st["A"]
    lu, cols, C, U, V, betas, S, J0 = _channel_basis(st)
    n_ch = 12
    grads = {}
    grads_fwd = {}
    for k in rod_indices:
        dA = assemble_dA_rod(
            points,
            st["tris"],
            st["rod_masks"][k],
            st["f_ord"],
            st["fp_ref_ord"],
            st["gamma"],
            st["fc"],
            st["k0"],
            s=float(s_vec[k]),
        )
        # Precompute dA u_n for TE0 drives
        dAu = {}
        dU = {}
        for j_src in range(6):
            n = 2 * j_src
            dAu[n] = dA @ U[:, n]
            dU[n] = lu.solve(-dAu[n])

        def dS_adj(m, n, _dAu=dAu, _V=V, _b=betas):
            sbeta = np.sqrt(_b[m] * _b[n])
            return -2j * sbeta * complex(np.dot(_V[:, m], _dAu[n]))

        def dS_fwd(m, n, _dU=dU, _C=C, _b=betas):
            sbeta = np.sqrt(_b[m] * _b[n])
            return 2j * sbeta * complex(np.dot(_C[:, m], _dU[n]))

        grads[int(k)] = _grad_from_dS(S, betas, dS_adj)
        grads_fwd[int(k)] = _grad_from_dS(S, betas, dS_fwd)
    return J0, S, grads, grads_fwd


def fd_gradient(grade, s0, rod, h, b_tesla, measure):
    sp = s0.copy()
    sm = s0.copy()
    # Clip only to the design box [S_LO, S_HI], never to 1.0 (many rods have s>1).
    sp[rod] = min(S_HI, sp[rod] + h)
    sm[rod] = max(S_LO, sm[rod] - h)
    hp = float(sp[rod] - s0[rod])
    hm = float(s0[rod] - sm[rod])
    denom = hp + hm
    if denom < 1e-18:
        raise RuntimeError(f"FD step collapsed for rod {rod} at s={s0[rod]}")
    st_p = build_system(grade, sp, 3.85e9, b_tesla, measure)
    st_m = build_system(grade, sm, 3.85e9, b_tesla, measure)
    Sp = direct_S_from_factor(
        splu(st_p["A"].tocsc()), len(st_p["points"]), build_Q_columns(st_p["points"], st_p["ports"], st_p["modes"], 2)
    )
    Sm = direct_S_from_factor(
        splu(st_m["A"].tocsc()), len(st_m["points"]), build_Q_columns(st_m["points"], st_m["ports"], st_m["modes"], 2)
    )
    Jp, Jm = objective_from_S(Sp), objective_from_S(Sm)
    return float((Jp - Jm) / denom), Jp, Jm, denom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--design", default="opt_single_uniform_M_Bp0.0500")
    ap.add_argument("--rods", default="0,45,30")
    ap.add_argument("--hs", default="1e-3,3e-4")
    ap.add_argument("--B", type=float, default=0.05)
    ap.add_argument("--measure", default="projected")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    s = load_s(args.design)
    rods = [int(x) for x in args.rods.split(",")]
    hs = [float(x) for x in args.hs.split(",")]
    t0 = time.time()
    st = build_system(args.grade, s, 3.85e9, args.B, args.measure)
    J0, S, g_adj, g_fwd = adjoint_gradient(st, s, rods)
    fwd_ok = True
    fwd_rows = []
    for rod in rods:
        ga, gf = g_adj[rod], g_fwd[rod]
        rel_f = abs(ga - gf) / max(abs(ga), abs(gf), 1e-30)
        ok_f = rel_f <= 1e-10
        fwd_ok = fwd_ok and ok_f
        fwd_rows.append({"rod": rod, "adjoint": ga, "forward": gf, "rel_err": float(rel_f), "PASS": ok_f})
        print(f"rod {rod} forward≡adjoint rel={rel_f:.3e} PASS={ok_f}", flush=True)
    rows = []
    fd_ok = True
    for rod in rods:
        best = None
        sweeps = []
        for h in hs:
            fd, Jp, Jm, h_eff = fd_gradient(args.grade, s, rod, h, args.B, args.measure)
            g = g_adj[rod]
            rel = abs(fd - g) / max(abs(fd), abs(g), 1e-30)
            rec = {"h": h, "h_eff": h_eff, "adjoint": g, "fd": fd, "rel_err": float(rel), "Jp": Jp, "Jm": Jm}
            sweeps.append(rec)
            if best is None or rel < best["rel_err"]:
                best = rec
            print(f"rod {rod} h={h:g} adj={g:.6e} fd={fd:.6e} rel={rel:.3e}", flush=True)
        ok = best["rel_err"] <= FD_REL_BAR
        fd_ok = fd_ok and ok
        rows.append({"rod": rod, "adjoint": g_adj[rod], "best": best, "h_sweep": sweeps, "PASS": ok})
    all_ok = fwd_ok and fd_ok
    report = {
        "grade": args.grade,
        "design": args.design,
        "B": args.B,
        "measure": args.measure,
        "J0": J0,
        "weights": {
            "des": W_DES, "rev": W_REV, "thr": W_THR,
            "leak": W_LEAK, "refl": W_REFL, "te1": W_TE1,
        },
        "fd_rel_bar": FD_REL_BAR,
        "forward_vs_adjoint": fwd_rows,
        "FORWARD_EQ_ADJOINT": "PASS" if fwd_ok else "FAIL",
        "rods": rows,
        "DIRECT_DTN_ADJOINT_FD": "PASS" if all_ok else "FAIL",
        "elapsed_s": time.time() - t0,
    }
    (OUT / f"adjoint_direct_dtn_{args.grade}.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = [
        "# Direct DtN multimode adjoint vs FD",
        "",
        f"Grade {args.grade}, design `{args.design}`, B={args.B}, measure={args.measure}.",
        f"J0 = {J0:.6e}. FD rel bar = {FD_REL_BAR}.",
        "",
        "## Forward ≡ adjoint (machine check)",
        "",
        "| rod | g_adj | g_fwd | rel | PASS |",
        "|-----|------:|------:|----:|:----:|",
    ]
    for r in fwd_rows:
        lines.append(
            f"| {r['rod']} | {r['adjoint']:.4e} | {r['forward']:.4e} | {r['rel_err']:.3e} | {r['PASS']} |"
        )
    lines += [
        "",
        f"**FORWARD_EQ_ADJOINT: {report['FORWARD_EQ_ADJOINT']}**",
        "",
        "## Centered FD",
        "",
        "| rod | g_adj | best FD | rel | PASS |",
        "|-----|------:|--------:|----:|:----:|",
    ]
    for r in rows:
        b = r["best"]
        lines.append(
            f"| {r['rod']} | {b['adjoint']:.4e} | {b['fd']:.4e} | {b['rel_err']:.3e} | {r['PASS']} |"
        )
    lines += ["", f"**DIRECT_DTN_ADJOINT_FD: {report['DIRECT_DTN_ADJOINT_FD']}**", ""]
    (OUT / "DIRECT_DTN_ADJOINT_FD_VALIDATION.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({
        "DIRECT_DTN_ADJOINT_FD": report["DIRECT_DTN_ADJOINT_FD"],
        "FORWARD_EQ_ADJOINT": report["FORWARD_EQ_ADJOINT"],
        "J0": J0,
        "elapsed_s": report["elapsed_s"],
    }, indent=2))
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
