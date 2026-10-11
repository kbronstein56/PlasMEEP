#!/usr/bin/env python3
"""Adjoint gradient of modal |S_ij|² / circulator S-metrics (post-validation).

For a fixed driven port j with incident â_j = A[j,j] treated as independent of
the design to first order when the feed is vacuum (source and monitor outside
the plasma), S_ij = b_i / a_j and

    ∂(|S_ij|²)/∂s = (2/|a_j|²) Re( b_i* ∂b_i/∂s )   if a_j frozen,

or with full quotient rule if a_j depends on s (generally weak for vacuum feeds).

Modal b at port i is linear in the nodal field x through the extraction:

    b = α† x ,   α assembled from φ-weighted (½ scale) Hz and ∂n terms.

Then ∂b/∂s_k = α† ∂x/∂s_k with A ∂x = −(∂A/∂s_k) x, so the adjoint solve
A† λ = α (or the metric-weighted combination of α's) yields

    ∂b/∂s_k = − λ† (∂A/∂s_k) x .

This module implements that gradient and a finite-difference spot check.
It is intended to run only after VALIDATION_SUMMARY.json reports all_PASS.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
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

from circulator_core import F_HZ, B_PLUS, load_device, solve_six  # noqa: E402
from diff_fem_core import assemble_dA_rod  # noqa: E402
from modal_s_core import dump_json, extract_ab, modal_port_powers, scattering_matrix  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design" / "modal_s"


def modal_b_functional(sampler, tris, line, k0: float) -> np.ndarray:
    """Complex vector α such that b ≈ α† x (power-normalized outgoing).

    From extract_ab:
      c_H = Σ w φ Hz = Σ_nodes (Σ_s w_s φ_s ℓ_{s,node}) x_node
      c_d = Σ w φ ∂n = Σ_nodes (Σ_s w_s φ_s g_n,node) x_node
      b_f = ½ (c_H + c_d/(iβ))
      b   = b_f √(β/(2ω))
    So α_node = scale * ½ * Σ_s w_s φ_s (ℓ_s + g_n/(iβ))
    and b = α† x with the conjugate absorbed into the definition
    Re(q† dx) style — here we keep b = α · x (linear form, not conjugate).
    """
    n = len(sampler.points)
    alpha = np.zeros(n, dtype=np.complex128)
    omega = float(k0)
    beta = float(k0)
    scale = np.sqrt(beta / (2.0 * omega))
    xy = np.asarray(line["monitor_xy"], float)
    n_hat = np.asarray(line["outward"], float)
    n_hat = n_hat / np.linalg.norm(n_hat)
    weights = np.asarray(line["weights"], float)
    off = np.asarray(line.get("offsets", []), float)
    span = float(off[-1] - off[0]) if len(off) > 1 else float(np.sum(weights))
    phi = 1.0 / np.sqrt(max(span, 1e-30))
    tid = sampler.locate(xy)
    for s in range(len(xy)):
        t = int(tid[s])
        if t < 0:
            continue
        nodes = tris[t]
        ell = np.asarray(sampler._bary(t, xy[s]), float)
        gn = n_hat[0] * sampler.bx[t] + n_hat[1] * sampler.by[t]
        ws = float(weights[s]) * phi
        coeff = 0.5 * scale * ws
        for k in range(3):
            nd = int(nodes[k])
            # b_f = ½(c_H + c_d/(iβ)) → contrib ℓ and gn/(iβ)
            alpha[nd] += coeff * (ell[k] + gn[k] / (1j * beta))
    return alpha


def modal_a_functional(sampler, tris, line, k0: float) -> np.ndarray:
    """α_a such that a = α_a · x (power-normalized incident)."""
    n = len(sampler.points)
    alpha = np.zeros(n, dtype=np.complex128)
    omega = float(k0)
    beta = float(k0)
    scale = np.sqrt(beta / (2.0 * omega))
    xy = np.asarray(line["monitor_xy"], float)
    n_hat = np.asarray(line["outward"], float)
    n_hat = n_hat / np.linalg.norm(n_hat)
    weights = np.asarray(line["weights"], float)
    off = np.asarray(line.get("offsets", []), float)
    span = float(off[-1] - off[0]) if len(off) > 1 else float(np.sum(weights))
    phi = 1.0 / np.sqrt(max(span, 1e-30))
    tid = sampler.locate(xy)
    for s in range(len(xy)):
        t = int(tid[s])
        if t < 0:
            continue
        nodes = tris[t]
        ell = np.asarray(sampler._bary(t, xy[s]), float)
        gn = n_hat[0] * sampler.bx[t] + n_hat[1] * sampler.by[t]
        ws = float(weights[s]) * phi
        coeff = 0.5 * scale * ws
        for k in range(3):
            nd = int(nodes[k])
            # a_f = ½(c_H - c_d/(iβ))
            alpha[nd] += coeff * (ell[k] - gn[k] / (1j * beta))
    return alpha


def dJ_ds_selected(
    dev,
    bun,
    j_src: int,
    rod_indices: list[int],
    *,
    freeze_a: bool = True,
) -> dict:
    """Gradient of J_S = T_des - 2 T_rev - 1.5 T_thr - T_leak - 0.25 |S_jj|²
    (modal analogue of the Pref objective weights) for selected rods.

    Uses quotient rule with optional freeze of a_j (vacuum-feed approximation).
    """
    from circulator_core import W_ACC, W_DES, W_LEAK, W_REV, W_THR, port_roles

    roles = port_roles(j_src)
    x = bun.fields[j_src]
    k0 = bun.k0
    A = bun.A
    # Functionals for each port's b and driven a
    alpha_b = [modal_b_functional(dev.sampler, dev.tris, line, k0) for line in dev.lines]
    alpha_a_j = modal_a_functional(dev.sampler, dev.tris, dev.lines[j_src], k0)
    b = np.array([alpha_b[i] @ x for i in range(6)], dtype=np.complex128)
    a_j = complex(alpha_a_j @ x)
    # Verify against extract_ab
    m_chk = extract_ab(dev.sampler, x, bun.rho, k0, dev.lines[j_src])
    # Metric weights on |S_pj|² = |b_p/a_j|²
    w_port = np.zeros(6, float)
    w_port[roles["desired"]] = W_DES
    w_port[roles["reverse"]] = -W_REV
    w_port[roles["through"]] = -W_THR
    w_port[roles["leak_a"]] = -W_LEAK
    w_port[roles["leak_b"]] = -W_LEAK
    w_port[roles["source"]] = -W_ACC  # penalize |S_jj|² like return

    # d(|b/a|²) = (2/|a|²) Re(b* db) - (2|b|²/|a|⁴) Re(a* da)
    # Adjoint RHS for the combined linear form on dx:
    # g† dx contributes to Σ w_p [ (2/|a|²) Re(b_p* (α_b_p·dx)) - ... ]
    # For real objective J = Σ w_p |S_p|², dJ = 2 Re(λ† ... ) style:
    # ∂J/∂x direction: for each p, contrib (w_p * 2 / |a|²) conj(b_p) α_b_p
    # minus (if not freeze) Σ_p w_p * 2 |b_p|²/|a|⁴ * conj(a) α_a
    inv_a2 = 1.0 / max(abs(a_j) ** 2, 1e-30)
    rhs = np.zeros(len(x), dtype=np.complex128)
    for p in range(6):
        if abs(w_port[p]) < 1e-30:
            continue
        rhs += (w_port[p] * 2.0 * inv_a2 * np.conj(b[p])) * alpha_b[p]
        if not freeze_a:
            rhs -= (w_port[p] * 2.0 * (abs(b[p]) ** 2) * (inv_a2**2) * np.conj(a_j)) * alpha_a_j

    # Solve A† λ = rhs  (for dJ = 2 Re(λ† (-dA x))? )
    # J is real; dJ = Re(rhs† dx) if rhs is built with the 2conj factors...
    # We built rhs so that dJ = Re(rhs† dx) with the factor 2 absorbed in rhs? 
    # |z|² differential: d(|z|²) = 2 Re(z* dz). We put (2 w /|a|²) conj(b) α into rhs,
    # so dJ = Re(rhs† dx) only if we use  Re(conj(coeff)* ...) carefully.
    # With rhs = Σ c_p α_p and c_p = (2w/|a|²) conj(b_p), then
    # rhs† dx = Σ conj(c_p) (α·dx)? No — we want Σ 2w/|a|² Re(b* db) = Re( Σ (2w/|a|²) conj(b) db )? 
    # 2 Re(b* db) = 2 Re(conj(b) db) wait: b* means conj(b).
    # 2 Re(conj(b) db) = b* db + b (db)* ...
    # Standard: dJ = Re(q† dx) with q = Σ (2w/|a|²) b α  if db=α·dx and we use
    #   2 Re(conj(b) (α·dx)) = 2 Re( conj(b α†? ) )
    # If db = α·dx (row vector), then 2 Re(conj(b) db) = 2 Re(conj(b) α · dx).
    # For dJ = Re(q† dx), need q† dx = 2 conj(b) α · dx + c.c. handling.
    # Equivalently q = 2 b α  (then q† dx = 2 conj(b) conj(α)·...? messy with complex α).
    #
    # Practical approach matching circulator_core: use forward FD on selected rods
    # for the check, and adjoint via A^H λ = v with v chosen so
    # dJ = −2 Re(λ† dA x) works for holomorphic |·|².
    #
    # Here: v := Σ_p (w_p / |a|²) conj(b_p) α_b_p   [and a-terms],
    # solve A^H λ = v, then ∂J/∂s = −2 Re(λ† (∂A/∂s) x).

    v = np.zeros(len(x), dtype=np.complex128)
    for p in range(6):
        if abs(w_port[p]) < 1e-30:
            continue
        v += (w_port[p] * inv_a2 * np.conj(b[p])) * alpha_b[p]
        if not freeze_a:
            v -= (w_port[p] * (abs(b[p]) ** 2) * (inv_a2**2) * np.conj(a_j)) * alpha_a_j

    # SuperLU on A^H: solve A^H λ = v  <=>  λ = (A^H)^{-1} v
    # splu(A).solve(v, trans='H')
    lu = bun.lu
    lam = lu.solve(v, trans="H")

    grads = {}
    for k in rod_indices:
        dA = assemble_dA_rod(
            dev.points,
            dev.tris,
            bun.rho,
            bun.k0,
            bun.f_ord,
            bun.gamma,
            bun.fc,
            bun.fp_ref_ord,
            bun.s_vec[k],
            dev.rod_masks[k],
        )
        # ∂J/∂s_k = −2 Re(λ† dA x)
        grads[int(k)] = float(-2.0 * np.real(np.vdot(lam, dA @ x)))

    return {
        "j_src": j_src,
        "a_j": [a_j.real, a_j.imag],
        "a_j_extract": [m_chk.a.real, m_chk.a.imag],
        "a_functional_err": abs(a_j - m_chk.a) / max(abs(m_chk.a), 1e-30),
        "b": [[z.real, z.imag] for z in b],
        "freeze_a": freeze_a,
        "grads": grads,
        "w_port": w_port.tolist(),
    }


def fd_check(dev, s0, f_hz, B, j_src, rod_indices, eps=1e-4, freeze_a=True):
    from circulator_core import W_ACC, W_DES, W_LEAK, W_REV, W_THR, port_roles

    def J_of(s):
        bun = solve_six(dev, s, f_hz, B, P_ref=1.0)  # Pref irrelevant for modal J
        A, Bmat, _ = modal_port_powers(dev.sampler, bun.fields, bun.rho, bun.k0, dev.lines)
        S = scattering_matrix(A, Bmat)
        r = port_roles(j_src)
        return (
            W_DES * abs(S[r["desired"], j_src]) ** 2
            - W_REV * abs(S[r["reverse"], j_src]) ** 2
            - W_THR * abs(S[r["through"], j_src]) ** 2
            - W_LEAK * (abs(S[r["leak_a"], j_src]) ** 2 + abs(S[r["leak_b"], j_src]) ** 2)
            - W_ACC * abs(S[r["source"], j_src]) ** 2
        )

    bun0 = solve_six(dev, s0, f_hz, B, P_ref=1.0)
    adj = dJ_ds_selected(dev, bun0, j_src, rod_indices, freeze_a=freeze_a)
    rows = []
    for k in rod_indices:
        s_p = s0.copy()
        s_m = s0.copy()
        s_p[k] += eps
        s_m[k] -= eps
        # Keep bounds soft
        s_p[k] = float(np.clip(s_p[k], 0.05, 1.80))
        s_m[k] = float(np.clip(s_m[k], 0.05, 1.80))
        Jp = J_of(s_p)
        Jm = J_of(s_m)
        fd = (Jp - Jm) / (s_p[k] - s_m[k])
        g = adj["grads"][int(k)]
        rel = abs(g - fd) / max(abs(fd), 1e-30)
        rows.append({"rod": int(k), "adjoint": g, "fd": float(fd), "rel_err": float(rel), "Jp": Jp, "Jm": Jm})
    return adj, rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--f-hz", type=float, default=F_HZ)
    ap.add_argument("--B", type=float, default=B_PLUS)
    ap.add_argument("--rods", default="0,7,14")
    ap.add_argument("--eps", type=float, default=1e-4)
    ap.add_argument("--require-validation", action="store_true", default=True)
    ap.add_argument("--no-require-validation", action="store_false", dest="require_validation")
    args = ap.parse_args()

    summary_path = OUT / "VALIDATION_SUMMARY.json"
    if args.require_validation:
        if not summary_path.is_file():
            print("FAIL: validation summary missing; refuse adjoint", flush=True)
            return 2
        summ = json.loads(summary_path.read_text())
        if not summ.get("all_PASS"):
            print("FAIL: modal S validation not all_PASS; refuse adjoint", flush=True)
            return 2

    rods = [int(x) for x in args.rods.split(",")]
    print(f"Loading grade={args.grade}", flush=True)
    dev = load_device(args.grade)
    s0 = np.ones(91)
    print(f"FD check rods={rods}", flush=True)
    t0 = time.perf_counter()
    adj, rows = fd_check(dev, s0, args.f_hz, args.B, j_src=0, rod_indices=rods, eps=args.eps)
    wall = time.perf_counter() - t0
    max_rel = max(r["rel_err"] for r in rows)
    out = {
        "grade": args.grade,
        "f_Hz": args.f_hz,
        "B_T": args.B,
        "adjoint_meta": adj,
        "fd_rows": rows,
        "max_rel_err": max_rel,
        "PASS": bool(max_rel < 0.15),
        "threshold_rel": 0.15,
        "wall_s": wall,
        "note": "No new 91-variable optimization launched.",
    }
    dump_json(OUT / f"modal_s_adjoint_fd_{args.grade}.json", out)
    print(json.dumps({k: out[k] for k in ("max_rel_err", "PASS", "fd_rows", "wall_s")}, indent=2), flush=True)
    return 0 if out["PASS"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
