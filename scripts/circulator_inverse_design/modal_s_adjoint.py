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


def _q_for_source(dev, bun, j_src: int, *, freeze_a: bool) -> tuple[np.ndarray, dict]:
    """Build adjoint RHS q with dJ_src = Re(q† dx) for one driven port.

    J_src = Σ_p w_p |b_p / a_j|² with Pref-analogue weights.
    With b = α·x and d(|b|²)=2 Re(conj(b) db):
      q_i = Σ_p (2 w_p / |a|²) b_p conj(α_b_p,i)   [freeze_a]
    plus quotient-rule terms in a_j when freeze_a is False.
    """
    from circulator_core import W_ACC, W_DES, W_LEAK, W_REV, W_THR, port_roles

    roles = port_roles(j_src)
    x = bun.fields[j_src]
    k0 = bun.k0
    alpha_b = [modal_b_functional(dev.sampler, dev.tris, line, k0) for line in dev.lines]
    alpha_a_j = modal_a_functional(dev.sampler, dev.tris, dev.lines[j_src], k0)
    b = np.array([alpha_b[i] @ x for i in range(6)], dtype=np.complex128)
    a_j = complex(alpha_a_j @ x)
    m_chk = extract_ab(dev.sampler, x, bun.rho, k0, dev.lines[j_src])
    w_port = np.zeros(6, float)
    w_port[roles["desired"]] = W_DES
    w_port[roles["reverse"]] = -W_REV
    w_port[roles["through"]] = -W_THR
    w_port[roles["leak_a"]] = -W_LEAK
    w_port[roles["leak_b"]] = -W_LEAK
    w_port[roles["source"]] = -W_ACC
    inv_a2 = 1.0 / max(abs(a_j) ** 2, 1e-30)
    q = np.zeros(len(x), dtype=np.complex128)
    for p in range(6):
        if abs(w_port[p]) < 1e-30:
            continue
        # q_i = (2w/|a|²) b conj(α_i)  ⇒  vdot(q,dx)= (2w/|a|²) conj(b) (α·dx)
        q += (w_port[p] * 2.0 * inv_a2 * b[p]) * np.conj(alpha_b[p])
        if not freeze_a:
            # − (2w |b|² / |a|⁴) Re(conj(a) da) → q -= (2w |b|²/|a|⁴) a conj(α_a)
            q -= (w_port[p] * 2.0 * (abs(b[p]) ** 2) * (inv_a2**2) * a_j) * np.conj(alpha_a_j)
    meta = {
        "a_j": [a_j.real, a_j.imag],
        "a_j_extract": [m_chk.a.real, m_chk.a.imag],
        "a_functional_err": abs(a_j - m_chk.a) / max(abs(m_chk.a), 1e-30),
        "b_functional_err": float(
            max(abs(b[i] - extract_ab(dev.sampler, x, bun.rho, k0, dev.lines[i]).b) for i in range(6))
            / max(abs(b[0]), 1e-30)
        ),
        "w_port": w_port.tolist(),
    }
    return q, meta


def dJ_ds_selected(
    dev,
    bun,
    rod_indices: list[int],
    *,
    freeze_a: bool = True,
) -> dict:
    """∂J/∂s for J = (1/6) Σ_j J_src(j), modal |S|² weights matching Pref FoM."""
    from scipy.sparse.linalg import splu

    lams = []
    metas = []
    for j in range(6):
        q, meta = _q_for_source(dev, bun, j, freeze_a=freeze_a)
        metas.append(meta)
        try:
            lam = bun.lu.solve(q, trans="H")
        except TypeError:
            lam = splu(bun.A.conj().T.tocsc()).solve(q)
        lams.append(lam)

    grads = {}
    for k in rod_indices:
        dA = assemble_dA_rod(
            dev.points,
            dev.tris,
            dev.rod_masks[k],
            bun.f_ord,
            bun.fp_ref_ord,
            bun.gamma,
            bun.fc,
            bun.k0,
            float(bun.s_vec[k]),
        )
        acc = 0.0
        for j in range(6):
            acc += -np.real(np.vdot(lams[j], dA @ bun.fields[j]))
        grads[int(k)] = float(acc / 6.0)

    return {
        "freeze_a": freeze_a,
        "grads": grads,
        "per_source_meta": metas,
        "max_a_functional_err": max(m["a_functional_err"] for m in metas),
        "max_b_functional_err": max(m["b_functional_err"] for m in metas),
    }


def modal_J(dev, bun) -> float:
    from circulator_core import W_ACC, W_DES, W_LEAK, W_REV, W_THR, port_roles

    A, Bmat, _ = modal_port_powers(dev.sampler, bun.fields, bun.rho, bun.k0, dev.lines)
    S = scattering_matrix(A, Bmat)
    acc = 0.0
    for j in range(6):
        r = port_roles(j)
        acc += (
            W_DES * abs(S[r["desired"], j]) ** 2
            - W_REV * abs(S[r["reverse"], j]) ** 2
            - W_THR * abs(S[r["through"], j]) ** 2
            - W_LEAK * (abs(S[r["leak_a"], j]) ** 2 + abs(S[r["leak_b"], j]) ** 2)
            - W_ACC * abs(S[r["source"], j]) ** 2
        )
    return float(acc / 6.0)


def fd_check(dev, s0, f_hz, B, rod_indices, eps=1e-3, freeze_a=False):
    def J_of(s):
        bun = solve_six(dev, s, f_hz, B, P_ref=1.0)
        return modal_J(dev, bun)

    bun0 = solve_six(dev, s0, f_hz, B, P_ref=1.0)
    adj = dJ_ds_selected(dev, bun0, rod_indices, freeze_a=freeze_a)
    rows = []
    for k in rod_indices:
        s_p = s0.copy()
        s_m = s0.copy()
        s_p[k] = float(np.clip(s0[k] + eps, 0.05, 1.80))
        s_m[k] = float(np.clip(s0[k] - eps, 0.05, 1.80))
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
    ap.add_argument("--eps", type=float, default=1e-3)
    ap.add_argument("--freeze-a", action="store_true", help="Omit da_j quotient terms (less accurate)")
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
    freeze_a = bool(args.freeze_a)
    adj, rows = fd_check(
        dev, s0, args.f_hz, args.B, rod_indices=rods, eps=args.eps, freeze_a=freeze_a
    )
    wall = time.perf_counter() - t0
    max_rel = max(r["rel_err"] for r in rows)
    out = {
        "grade": args.grade,
        "f_Hz": args.f_hz,
        "B_T": args.B,
        "freeze_a": freeze_a,
        "eps": args.eps,
        "adjoint_meta": {
            "max_a_functional_err": adj["max_a_functional_err"],
            "max_b_functional_err": adj["max_b_functional_err"],
            "grads": adj["grads"],
        },
        "fd_rows": rows,
        "max_rel_err": max_rel,
        "PASS": bool(max_rel < 0.05),
        "threshold_rel": 0.05,
        "wall_s": wall,
        "note": (
            "Full quotient rule (freeze_a=False) required for <1e-4 FD agreement. "
            "No new 91-variable optimization launched."
        ),
    }
    dump_json(OUT / f"modal_s_adjoint_fd_{args.grade}.json", out)
    print(json.dumps({k: out[k] for k in ("max_rel_err", "PASS", "fd_rows", "wall_s")}, indent=2), flush=True)
    return 0 if out["PASS"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
