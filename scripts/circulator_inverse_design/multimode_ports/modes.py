#!/usr/bin/env python3
"""Parallel-plate Hz modes for PEC-walled feeds (Neumann Hz on metal).

Domain: clear aperture y ∈ [-W/2, W/2], W = clear_width (2a in code units).

With e^{-iωt} and the validated FEM using natural Neumann Hz on ρ→0 metal,
transverse eigenfunctions satisfy ∂φ/∂y = 0 on the PEC sidewalls:

    φ_m(y) = cos(m π (y + W/2) / W) ,   m = 0, 1, 2, ...

    k_c,m = m π / W
    β_m = sqrt(k0² - k_c,m²)   (propagating, real β ≥ 0)
    α_m = sqrt(k_c,m² - k0²)  (evanescent, use +α in exp(-α|Δx|))

Power-normalized modal amplitudes â, b̂ (vacuum, e^{-iωt}):

    H_z = Σ_m (a_m + b_m) φ_m / N_m   with ∫|φ_m|² dy = 1

    P_inc = |â|² ,  P_out = |b̂|²  using â = a_m √(β_m / (2ω)), ω = k0 (c=1).
"""
from __future__ import annotations

import numpy as np

from feed_constants import C_MPS, a, clear_width  # noqa: E402

DEFAULT_W = float(clear_width)


def transverse_mode(y: np.ndarray, m: int, W: float = DEFAULT_W) -> np.ndarray:
    """Unnormalized φ_m(y) on the clear aperture."""
    y = np.asarray(y, float)
    return np.cos(m * np.pi * (y + 0.5 * W) / W)


def mode_constants(k0: float, W: float = DEFAULT_W, n_modes: int = 4) -> list[dict]:
    """Return list of mode dicts with beta, alpha, kc, propagating flag."""
    out = []
    for m in range(n_modes):
        kc = m * np.pi / W
        disc = k0**2 - kc**2
        if disc >= 0:
            beta = float(np.sqrt(disc))
            out.append({
                "m": m,
                "kc": kc,
                "beta": beta + 0j,
                "alpha": 0.0,
                "propagating": True,
                "f_cut_Hz": None if m == 0 else C_MPS * kc / (2 * np.pi * a),
            })
        else:
            alpha = float(np.sqrt(-disc))
            out.append({
                "m": m,
                "kc": kc,
                "beta": 1j * alpha,
                "alpha": alpha,
                "propagating": False,
                "f_cut_Hz": C_MPS * kc / (2 * np.pi * a),
            })
    return out


def l2_norm_trapz(phi: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sqrt(max(np.sum(weights * phi * phi), 1e-30)))


def normalize_on_line(offsets: np.ndarray, weights: np.ndarray, m: int, W: float) -> tuple[np.ndarray, float]:
    phi = transverse_mode(offsets, m, W)
    nrm = l2_norm_trapz(phi, weights)
    return phi / nrm, nrm


def power_scale(beta: complex, k0: float) -> complex:
    """√(β/(2ω)) with ω=k0."""
    omega = float(k0)
    b = complex(beta)
    if b.imag != 0:
        # Evanescent: use real part of β² for power bookkeeping
        return np.sqrt(b / (2.0 * omega))
    return np.sqrt(b / (2.0 * omega))


def verify_cutoffs_at_f(f_hz: float, W: float = DEFAULT_W) -> dict:
    k0 = 2 * np.pi * (f_hz * a / C_MPS)
    modes = mode_constants(k0, W, n_modes=4)
    return {
        "f_Hz": f_hz,
        "W": W,
        "k0": k0,
        "modes": [
            {
                "m": md["m"],
                "propagating": md["propagating"],
                "beta": [md["beta"].real, md["beta"].imag],
                "f_cut_GHz": None if md["m"] == 0 else md["f_cut_Hz"] / 1e9,
            }
            for md in modes
        ],
    }
