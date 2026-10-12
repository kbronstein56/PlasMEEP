#!/usr/bin/env python3
"""Two-plane fundamental modal extraction (no ∂n).

At monitor offsets ζ1, ζ2 along outward (local coordinate increasing outward):

    c(ζ) = ∫ φ* Hz ds = a e^{+i β ζ_in} + b e^{+i β ζ_out}

With ζ measured along outward from a reference, incident (~−n̂) ~ e^{+iβ ζ}
if ζ increases outward (incident phase increases as we go outward, toward source).

Convention consistent with modal_s_core:
  incident a travels −n̂; at position ζ (outward from ref):
    Hz_inc ~ a exp(+i β ζ)   # because moving +outward is against propagation → +phase for e^{-iωt}? 

Careful with e^{-iωt}: wave traveling in +x is exp(+i k x).
Traveling −n̂ (inward): field ~ exp(−i β ζ) if ζ along outward.
Traveling +n̂ (outward): field ~ exp(+i β ζ).

So c(ζ) = a exp(−i β ζ) + b exp(+i β ζ).

Solve 2×2 at ζ1, ζ2.
"""
from __future__ import annotations

import numpy as np

from modal_s_core import sample_hz_dn


def _project_c(sampler, x, rho, k0, line, zeta_shift: float):
    """Project flat mode at monitor shifted by zeta_shift along outward."""
    n_hat = np.asarray(line["outward"], float)
    n_hat = n_hat / np.linalg.norm(n_hat)
    tang = np.asarray(line["tangent"], float)
    off = np.asarray(line["offsets"], float)
    weights = np.asarray(line["weights"], float)
    span = float(off[-1] - off[0])
    phi = np.full(len(off), 1.0 / np.sqrt(span))
    # base monitor center = mean of monitor_xy
    xy0 = np.asarray(line["monitor_xy"], float)
    center = xy0.mean(axis=0)
    xy = np.vstack([center + zeta_shift * n_hat + float(s) * tang for s in off])
    Hz, _ = sample_hz_dn(sampler, x, rho, float(k0), xy, n_hat)
    return complex(np.sum(weights * phi * Hz))


def extract_ab_two_plane(sampler, x, rho, k0, line, d_plane: float = 0.25) -> dict:
    """Extract a,b using planes at ζ=−d/2 and +d/2 relative to nominal monitor."""
    beta = float(k0)
    omega = float(k0)
    z1, z2 = -0.5 * d_plane, +0.5 * d_plane
    c1 = _project_c(sampler, x, rho, k0, line, z1)
    c2 = _project_c(sampler, x, rho, k0, line, z2)
    # c = a e^{-iβζ} + b e^{+iβζ}
    # [e^{-iβz1} e^{iβz1}] [a] = [c1]
    # [e^{-iβz2} e^{iβz2}] [b]   [c2]
    M = np.array([
        [np.exp(-1j * beta * z1), np.exp(+1j * beta * z1)],
        [np.exp(-1j * beta * z2), np.exp(+1j * beta * z2)],
    ], dtype=np.complex128)
    ab = np.linalg.solve(M, np.array([c1, c2], dtype=np.complex128))
    a_f, b_f = complex(ab[0]), complex(ab[1])
    scale = np.sqrt(beta / (2.0 * omega))
    a, b = a_f * scale, b_f * scale
    return {
        "a": a, "b": b,
        "P_inc": float(abs(a) ** 2),
        "P_out": float(abs(b) ** 2),
        "d_plane": d_plane,
        "c1": [c1.real, c1.imag],
        "c2": [c2.real, c2.imag],
    }


def scattering_two_plane(dev, fields, rho, k0, lines, d_plane=0.25):
    A = np.zeros((6, 6), dtype=np.complex128)
    B = np.zeros((6, 6), dtype=np.complex128)
    for j, x in enumerate(fields):
        for i, line in enumerate(lines):
            e = extract_ab_two_plane(dev.sampler, x, rho, k0, line, d_plane=d_plane)
            A[i, j] = e["a"]
            B[i, j] = e["b"]
    S = np.zeros((6, 6), dtype=np.complex128)
    for j in range(6):
        S[:, j] = B[:, j] / A[j, j]
    return S, A, B
