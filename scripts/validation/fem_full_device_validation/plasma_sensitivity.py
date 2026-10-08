#!/usr/bin/env python3
"""Analytic dε/ds, dρ/ds, and dA/ds for a density scale on plasma elements.

s = 1 is the reference plasma frequency. fp^2(s) = s * fp_ref^2.
The mass term of A does not depend on s. dA is the stiffness of dρ only.
"""
from __future__ import annotations

import numpy as np

from gyrotropic_tensor import rho_xy, tensor_ordinary


def deps_ds(f: float, fp_ref: float, gamma: float, fc: float, s: float = 1.0):
    """Derivative of the in-plane permittivity with respect to s."""
    fp = fp_ref * np.sqrt(max(s, 0.0))
    # d(eps)/d(fp^2) * d(fp^2)/ds, and fp^2 = s * fp_ref^2.
    u = f + 1j * gamma
    denom = u * u - fc * fc
    de_perp = -(fp_ref ** 2) * u / (f * denom)
    deta = (fp_ref ** 2) * fc / (f * denom)
    return de_perp, 1j * deta, -1j * deta, de_perp, fp


def drho_ds(f: float, fp_ref: float, gamma: float, fc: float, s: float = 1.0):
    fp = fp_ref * np.sqrt(s)
    exx, exy, eyx, eyy, _, _ = tensor_ordinary(f, fp, gamma, fc)
    rxx, rxy, ryx, ryy, _ = rho_xy(exx, exy, eyx, eyy)
    dexx, dexy, deyx, deyy, _ = deps_ds(f, fp_ref, gamma, fc, s)
    # dρ = -ρ (dε) ρ
    rho = np.array([[rxx, rxy], [ryx, ryy]], np.complex128)
    deps = np.array([[dexx, dexy], [deyx, deyy]], np.complex128)
    dr = -rho @ deps @ rho
    return dr[0, 0], dr[0, 1], dr[1, 0], dr[1, 1]


def element_drho(mask: np.ndarray, f, fp_ref, gamma, fc, s: float = 1.0):
    dxx, dxy, dyx, dyy = drho_ds(f, fp_ref, gamma, fc, s)
    z = np.zeros(len(mask), np.complex128)
    out = [z.copy(), z.copy(), z.copy(), z.copy()]
    out[0][mask] = dxx
    out[1][mask] = dxy
    out[2][mask] = dyx
    out[3][mask] = dyy
    return out
