#!/usr/bin/env python3
"""True modal six-port scattering amplitudes for the FULL91 FEM.

Convention (matches production horn geometry)
--------------------------------------------
- Port outward normal n̂ points *away* from the device.
- Source is outside the monitor (src·n̂ > mon·n̂); launch is into the device (−n̂).
- At a monitor plane: incident a travels −n̂; outgoing b travels +n̂.
- Fundamental feed mode is the flat Hz PEC mode used by full91_solve.port_lines.
- Propagation constant in the vacuum feed: β = k0 (TEM / cutoff-0).

Bidirectional decomposition (phase reference at the monitor)
-----------------------------------------------------------
With L2-normalized flat profile φ (∫|φ|² ds = 1):

    Hz      = (a_f + b_f) φ
    ∂_n Hz  = i β (b_f − a_f) φ

Power-normalized amplitudes (vacuum feed, e^{-iωt}):

    â = a_f √(β / (2 ω)),   ˆb = b_f √(β / (2 ω))

so that |â|² and |ˆb|² equal incident and outgoing modal powers
(∫ S·n̂ = |ˆb|² − |â|² for a pure bidirectional fundamental).

Scattering matrix for a unit incident at driven port j (a_i=0, i≠j):

    S_ij = ˆb_i / â_j

Net Poynting flux is *not* used as a proxy for reflection; a and b are
separated by the (Hz, ∂n Hz) modal projection.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def flat_mode_weights(n_points: int, span: float) -> tuple[np.ndarray, np.ndarray, float]:
    """Trapezoid quadrature and L2-normalized constant mode φ = 1/√span."""
    offsets = np.linspace(-span / 2.0, span / 2.0, n_points)
    if n_points == 1:
        w = np.array([span], float)
    else:
        ds = span / (n_points - 1)
        w = np.full(n_points, ds, float)
        w[0] *= 0.5
        w[-1] *= 0.5
    phi = np.full(n_points, 1.0 / np.sqrt(span), dtype=float)
    return offsets, w, phi


def sample_hz_dn(
    sampler,
    x: np.ndarray,
    rho,
    omega: float,
    xy: np.ndarray,
    n_hat: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Hz and ∂n Hz on sample points via element barycentric gradients."""
    xy = np.atleast_2d(np.asarray(xy, float))
    n_hat = np.asarray(n_hat, float)
    n_hat = n_hat / np.linalg.norm(n_hat)
    tid = sampler.locate(xy)
    Hz = np.zeros(len(xy), dtype=np.complex128)
    dn = np.zeros(len(xy), dtype=np.complex128)
    found = tid >= 0
    if np.any(found):
        t = tid[found]
        nodes = sampler.tris[t]
        w = np.zeros((len(t), 3), dtype=float)
        for i, ti in enumerate(t):
            w[i] = sampler._bary(int(ti), xy[found][i])
        un = x[nodes]
        Hz[found] = np.sum(w * un, axis=1)
        gx = np.sum(sampler.bx[t] * un, axis=1)
        gy = np.sum(sampler.by[t] * un, axis=1)
        dn[found] = n_hat[0] * gx + n_hat[1] * gy
    if np.any(~found):
        for i in np.flatnonzero(~found):
            j = int(np.argmin(np.sum((sampler.points - xy[i]) ** 2, axis=1)))
            Hz[i] = x[j]
    return Hz, dn


@dataclass
class ModalAmplitudes:
    a: complex  # power-normalized incident (−n̂)
    b: complex  # power-normalized outgoing (+n̂)
    a_field: complex
    b_field: complex
    P_inc: float
    P_out: float
    P_net_modal: float
    P_net_poynting: float
    overlap_quality: float  # |projected Hz|² / ∫|Hz|² ; ~1 for pure fundamental
    te1_power_frac: float = 0.0  # |c_H1|² / ∫|Hz|² for odd TE1
    te1_P_out: float = 0.0  # outgoing TE1 modal power (if propagating)


def extract_ab(
    sampler,
    x: np.ndarray,
    rho,
    k0: float,
    line: dict,
    beta: Optional[float] = None,
    report_te1: bool = True,
) -> ModalAmplitudes:
    """Extract power-normalized a,b at one port monitor line."""
    omega = float(k0)  # c=1 units: ω = k0
    beta = float(k0 if beta is None else beta)
    xy = np.asarray(line["monitor_xy"], float)
    n_hat = np.asarray(line["outward"], float)
    n_hat = n_hat / np.linalg.norm(n_hat)
    weights = np.asarray(line["weights"], float)
    span = float(np.sum(weights))  # trapezoid total length ≈ span
    # Prefer explicit span from offsets if present
    if "offsets" in line:
        off = np.asarray(line["offsets"], float)
        span = float(off[-1] - off[0]) if len(off) > 1 else float(np.sum(weights))
    else:
        off = np.linspace(-span / 2.0, span / 2.0, len(weights))
    phi = np.full(len(weights), 1.0 / np.sqrt(max(span, 1e-30)), dtype=float)

    Hz, dn = sample_hz_dn(sampler, x, rho, omega, xy, n_hat)
    # Overlaps: c = ∫ φ* f ds ; N = ∫|φ|² = 1
    c_H = complex(np.sum(weights * phi * Hz))
    c_d = complex(np.sum(weights * phi * dn))
    # a_f + b_f = c_H;  iβ (b_f - a_f) = c_d
    sum_ab = c_H
    diff_ba = c_d / (1j * beta) if abs(beta) > 1e-30 else 0.0
    a_f = 0.5 * (sum_ab - diff_ba)
    b_f = 0.5 * (sum_ab + diff_ba)

    scale = np.sqrt(beta / (2.0 * omega)) if omega > 0 else 0.0
    a = a_f * scale
    b = b_f * scale
    P_inc = float(abs(a) ** 2)
    P_out = float(abs(b) ** 2)

    # Poynting cross-check on the same samples (vacuum ρ on feed)
    # S·n = ½ Re[(-i/ω) (∂n Hz) Hz*]
    alpha = -1j / omega
    P_poynt = float(np.sum(weights * 0.5 * np.real(alpha * dn * np.conj(Hz))))

    hz_pow = float(np.sum(weights * np.abs(Hz) ** 2))
    proj_pow = float(abs(c_H) ** 2)  # since ∫|φ|²=1, Parseval for one mode
    qual = proj_pow / max(hz_pow, 1e-30)

    te1_frac = 0.0
    te1_P_out = 0.0
    if report_te1:
        # Odd Neumann TE1 on clear width W (use production clear_width via span/0.92)
        W = span / 0.92 if span > 0 else span
        phi1 = np.sin(np.pi * off / max(W, 1e-30))
        nrm = float(np.sqrt(max(np.sum(weights * phi1 * phi1), 1e-30)))
        phi1 = phi1 / nrm
        c_H1 = complex(np.sum(weights * phi1 * Hz))
        c_d1 = complex(np.sum(weights * phi1 * dn))
        te1_frac = float(abs(c_H1) ** 2 / max(hz_pow, 1e-30))
        kc1 = np.pi / max(W, 1e-30)
        disc = omega**2 - kc1**2
        if disc > 0:
            beta1 = float(np.sqrt(disc))
            b1_f = 0.5 * (c_H1 + c_d1 / (1j * beta1))
            te1_P_out = float(abs(b1_f * np.sqrt(beta1 / (2.0 * omega))) ** 2)

    return ModalAmplitudes(
        a=complex(a),
        b=complex(b),
        a_field=complex(a_f),
        b_field=complex(b_f),
        P_inc=P_inc,
        P_out=P_out,
        P_net_modal=P_out - P_inc,
        P_net_poynting=P_poynt,
        overlap_quality=qual,
        te1_power_frac=te1_frac,
        te1_P_out=te1_P_out,
    )


def modal_port_powers(
    sampler,
    fields: Sequence[np.ndarray],
    rho,
    k0: float,
    lines: Sequence[dict],
) -> tuple[np.ndarray, np.ndarray, list]:
    """For each source j and port i, return a[i,j], b[i,j] and detail list."""
    n = len(lines)
    A = np.zeros((n, n), dtype=np.complex128)
    B = np.zeros((n, n), dtype=np.complex128)
    details = []
    for j, x in enumerate(fields):
        row = []
        for i, line in enumerate(lines):
            m = extract_ab(sampler, x, rho, k0, line)
            A[i, j] = m.a
            B[i, j] = m.b
            row.append({
                "a": [m.a.real, m.a.imag],
                "b": [m.b.real, m.b.imag],
                "P_inc": m.P_inc,
                "P_out": m.P_out,
                "P_net_modal": m.P_net_modal,
                "P_net_poynting": m.P_net_poynting,
                "overlap_quality": m.overlap_quality,
                "te1_power_frac": m.te1_power_frac,
                "te1_P_out": m.te1_P_out,
            })
        details.append(row)
    return A, B, details


def scattering_matrix(A: np.ndarray, B: np.ndarray, eps: float = 1e-30) -> np.ndarray:
    """S_ij = b_i / a_j using the driven-port incident amplitude a_j = A[j,j]."""
    n = A.shape[0]
    S = np.zeros((n, n), dtype=np.complex128)
    for j in range(n):
        aj = A[j, j]
        if abs(aj) < eps:
            S[:, j] = np.nan
        else:
            S[:, j] = B[:, j] / aj
    return S


def port_roles(j: int) -> dict:
    return {
        "desired": (j + 1) % 6,
        "reverse": (j - 1) % 6,
        "through": (j + 3) % 6,
        "leak_a": (j + 2) % 6,
        "leak_b": (j + 4) % 6,
        "source": j,
    }


def circulator_s_metrics(S: np.ndarray, A: np.ndarray, B: np.ndarray) -> dict:
    """Conventional S-parameter circulator metrics (CCW desired at +B)."""
    rows = []
    for j in range(6):
        r = port_roles(j)
        S_des = S[r["desired"], j]
        S_rev = S[r["reverse"], j]
        S_thr = S[r["through"], j]
        S_jj = S[j, j]
        leak = abs(S[r["leak_a"], j]) ** 2 + abs(S[r["leak_b"], j]) ** 2
        P_inc = float(abs(A[j, j]) ** 2)
        P_out_sum = float(np.sum(np.abs(B[:, j]) ** 2))
        P_abs = P_inc - P_out_sum
        rows.append({
            "source": j,
            "S_des": [S_des.real, S_des.imag],
            "S_rev": [S_rev.real, S_rev.imag],
            "S_thr": [S_thr.real, S_thr.imag],
            "S_jj": [S_jj.real, S_jj.imag],
            "T_des": float(abs(S_des) ** 2),
            "T_rev": float(abs(S_rev) ** 2),
            "T_thr": float(abs(S_thr) ** 2),
            "T_jj": float(abs(S_jj) ** 2),
            "T_leak": float(leak),
            "IL_dB": float(-10 * np.log10(max(abs(S_des) ** 2, 1e-30))),
            "isolation_dB": float(
                10 * np.log10(max(abs(S_des) ** 2, 1e-30) / max(abs(S_rev) ** 2, 1e-30))
            ),
            "return_loss_dB": float(-10 * np.log10(max(abs(S_jj) ** 2, 1e-30))),
            "P_inc": P_inc,
            "P_out_sum": P_out_sum,
            "P_absorbed_modal": P_abs,
        })

    def avg(key):
        return float(np.mean([r[key] for r in rows]))

    return {
        "per_source": rows,
        "T_des_avg": avg("T_des"),
        "T_rev_avg": avg("T_rev"),
        "T_thr_avg": avg("T_thr"),
        "T_jj_avg": avg("T_jj"),
        "T_leak_avg": avg("T_leak"),
        "IL_dB_avg": avg("IL_dB"),
        "isolation_dB_avg": avg("isolation_dB"),
        "return_loss_dB_avg": avg("return_loss_dB"),
        "P_absorbed_modal_avg": avg("P_absorbed_modal"),
        "P_inc_avg": avg("P_inc"),
        "note": (
            "IL_dB = -10 log10(|S_des|^2) is conventional insertion loss. "
            "It is NOT 10 log10(P/P_ref) from the frozen Pref campaign."
        ),
    }


def complex_matrix_to_json(M: np.ndarray) -> list:
    return [[[float(z.real), float(z.imag)] for z in row] for row in M]


def frobenius_rel(A: np.ndarray, B: np.ndarray) -> float:
    return float(np.linalg.norm(A - B) / max(np.linalg.norm(A), 1e-30))


def reciprocity_error(S: np.ndarray) -> float:
    """‖S − Sᵀ‖_F / ‖S‖_F"""
    return frobenius_rel(S, S.T)


def onsager_error(S_plus: np.ndarray, S_minus: np.ndarray) -> float:
    """‖S(+B) − Sᵀ(−B)‖_F / ‖S(+B)‖_F"""
    return frobenius_rel(S_plus, S_minus.T)


def dump_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n")
