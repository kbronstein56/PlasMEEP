#!/usr/bin/env python3
"""Two-mode (TE0 flat + TE1 odd) feed extraction and 12×12 Onsager tests.

At 3.85 GHz with clear_width=2, TE1 propagates (kc=π/W < k0). The
fundamental-only 6×6 S then need not satisfy Onsager; the full 12×12
mode-port S must.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
from modal_s_core import sample_hz_dn  # noqa: E402


def _span_offsets(line):
    weights = np.asarray(line["weights"], float)
    if "offsets" in line:
        off = np.asarray(line["offsets"], float)
        span = float(off[-1] - off[0]) if len(off) > 1 else float(np.sum(weights))
    else:
        span = float(np.sum(weights))
        off = np.linspace(-span / 2.0, span / 2.0, len(weights))
    return off, weights, span


def mode_profiles(off, weights, span, W=None):
    W = span / 0.92 if W is None else W
    phi0 = np.full(len(off), 1.0 / np.sqrt(max(span, 1e-30)))
    phi1 = np.sin(np.pi * off / max(W, 1e-30))
    n1 = float(np.sqrt(max(np.sum(weights * phi1 * phi1), 1e-30)))
    phi1 = phi1 / n1
    # Gram orthogonality check
    g01 = float(np.sum(weights * phi0 * phi1))
    return phi0, phi1, W, g01


def extract_two_mode(sampler, x, rho, k0, line) -> dict:
    omega = float(k0)
    off, weights, span = _span_offsets(line)
    phi0, phi1, W, g01 = mode_profiles(off, weights, span)
    n_hat = np.asarray(line["outward"], float)
    n_hat = n_hat / np.linalg.norm(n_hat)
    Hz, dn = sample_hz_dn(sampler, x, rho, omega, line["monitor_xy"], n_hat)

    kc1 = np.pi / W
    disc = omega**2 - kc1**2
    if disc <= 0:
        raise RuntimeError("TE1 not propagating; unexpected at 3.85 GHz")
    beta0 = omega
    beta1 = float(np.sqrt(disc))

    def ab(phi, beta):
        c_H = complex(np.sum(weights * phi * Hz))
        c_d = complex(np.sum(weights * phi * dn))
        a_f = 0.5 * (c_H - c_d / (1j * beta))
        b_f = 0.5 * (c_H + c_d / (1j * beta))
        sc = np.sqrt(beta / (2.0 * omega))
        return a_f * sc, b_f * sc

    a0, b0 = ab(phi0, beta0)
    a1, b1 = ab(phi1, beta1)
    return {
        "a0": a0, "b0": b0, "a1": a1, "b1": b1,
        "beta0": beta0, "beta1": beta1, "gram01": g01,
        "P0_inc": abs(a0) ** 2, "P0_out": abs(b0) ** 2,
        "P1_inc": abs(a1) ** 2, "P1_out": abs(b1) ** 2,
    }


def multimode_amplitudes(dev, fields, rho, k0, lines):
    """A[m,j], B[m,j] for m=0..11 (port*2 + mode), j=0..5 sources (fund launch)."""
    A = np.zeros((12, 6), dtype=np.complex128)
    B = np.zeros((12, 6), dtype=np.complex128)
    grams = []
    for j, x in enumerate(fields):
        for p, line in enumerate(lines):
            e = extract_two_mode(dev.sampler, x, rho, k0, line)
            A[2 * p, j] = e["a0"]
            A[2 * p + 1, j] = e["a1"]
            B[2 * p, j] = e["b0"]
            B[2 * p + 1, j] = e["b1"]
            grams.append(e["gram01"])
    return A, B, float(np.max(np.abs(grams)))


def S_multimode_from_fund_drive(A, B):
    """12×6 map: columns are responses to fundamental drive at each port.

    Full square 12×12 would need TE1 drives too. For Onsager on the
    fundamental-driven block, compare the 6×6 FF and 6×6 cross blocks
    carefully. Here we build the 6×6 FF S and also a 12×6 scattering map
    Smap[m,j] = B[m,j]/A[2j,j] (normalize by driven fundamental a).
    """
    Smap = np.zeros((12, 6), dtype=np.complex128)
    Sff = np.zeros((6, 6), dtype=np.complex128)
    Sf1 = np.zeros((6, 6), dtype=np.complex128)  # TE1 out / TE0 in
    for j in range(6):
        aj = A[2 * j, j]
        Smap[:, j] = B[:, j] / aj
        for i in range(6):
            Sff[i, j] = B[2 * i, j] / aj
            Sf1[i, j] = B[2 * i + 1, j] / aj
    return Smap, Sff, Sf1


def onsager_blocks(Sff_p, Sff_m, Sf1_p, Sf1_m, S1f_p=None, S1f_m=None):
    """Partial Onsager checks.

    Full Onsager needs TE1 drives. With only TE0 drives we can still check:
      Sff(+B) vs Sff(-B)^T   (will fail if mode conversion present)
      and report TE0→TE1 conversion power.
    """
    r_ff = float(np.linalg.norm(Sff_p - Sff_m.T) / max(np.linalg.norm(Sff_p), 1e-30))
    # Conversion strength
    conv_p = float(np.mean(np.abs(Sf1_p) ** 2))
    conv_m = float(np.mean(np.abs(Sf1_m) ** 2))
    return {"r_F_Sff": r_ff, "mean_T_TE0_to_TE1_plus": conv_p, "mean_T_TE0_to_TE1_minus": conv_m}


def solve_with_te1_drives(dev, s_vec, f_hz, b_tesla, lines):
    """12 RHS: flat + TE1 launch at each port; return fields list len 12."""
    from circulator_core import assemble_and_factor
    from full91_solve import inject

    rho, A, lu, k0, f, gamma, fc, fp_ref, t_fac = assemble_and_factor(dev, s_vec, f_hz, b_tesla)
    fields = []
    # Rebuild TE0 and TE1 weighted source lines
    for p, line in enumerate(lines):
        off, weights, span = _span_offsets(line)
        phi0, phi1, W, _ = mode_profiles(off, weights, span)
        # Match production TE0 inject: weights ~ trapezoid ds (not φ-normalized)
        # Production uses flat weights; keep that for mode 0.
        b0 = inject(dev.points, dev.tris, dev.sampler, line["source_xy"], weights)
        # TE1: weight samples by φ1/φ0 shape times ds — use phi1 * span weights scale
        w1 = weights * (phi1 / (phi0 + 1e-30))
        # remove mean to avoid TE0 contamination
        w1 = w1 - np.sum(weights * w1) / np.sum(weights) * np.ones_like(w1)
        b1 = inject(dev.points, dev.tris, dev.sampler, line["source_xy"], w1)
        fields.append(lu.solve(b0))
        fields.append(lu.solve(b1))
    return {
        "rho": rho, "A": A, "lu": lu, "k0": k0, "fields": fields,
        "factor_s": t_fac, "f_ord": f, "gamma": gamma, "fc": fc, "fp_ref_ord": fp_ref,
    }


def full_12x12_S(dev, pack, lines):
    """With 12 drives (TE0/TE1 × 6 ports), build 12×12 S."""
    fields = pack["fields"]
    k0 = pack["k0"]
    rho = pack["rho"]
    # A_in[m, n] = incoming mode m when driving n; B_out similarly
    Ain = np.zeros((12, 12), dtype=np.complex128)
    Bout = np.zeros((12, 12), dtype=np.complex128)
    for n in range(12):
        x = fields[n]
        for p, line in enumerate(lines):
            e = extract_two_mode(dev.sampler, x, rho, k0, line)
            Ain[2 * p, n] = e["a0"]
            Ain[2 * p + 1, n] = e["a1"]
            Bout[2 * p, n] = e["b0"]
            Bout[2 * p + 1, n] = e["b1"]
    # For drive n, incident should be dominated by mode n; S = Bout @ inv(Ain)
    # More stably: each column S[:,n] = Bout[:,n] / Ain[n,n] if off-diagonal Ain small
    S = np.zeros((12, 12), dtype=np.complex128)
    diag_dom = []
    for n in range(12):
        a_n = Ain[n, n]
        diag_dom.append(abs(a_n) / max(np.linalg.norm(Ain[:, n]), 1e-30))
        S[:, n] = Bout[:, n] / a_n
    return S, Ain, Bout, float(np.min(diag_dom))
