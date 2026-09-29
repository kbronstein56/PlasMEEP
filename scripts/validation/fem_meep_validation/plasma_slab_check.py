#!/usr/bin/env python3
"""Normal-incidence plasma slab. Analytic transfer matrix versus the FEM weak form."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "overnight_b0"


def analytic(k0, eps, L):
    kp = k0 * np.sqrt(eps)
    # Unknowns R, A, B, T.
    # x=0: 1+R = A+B
    #       ik0(1-R) = (ikp/eps)(A-B)
    # x=L: A e^{ikp L}+B e^{-ikp L} = T
    #       (ikp/eps)(A e^{ikp L}-B e^{-ikp L}) = ik0 T
    e = np.exp(1j * kp * L)
    em = np.exp(-1j * kp * L)
    M = np.array(
        [
            [1, -1, -1, 0],
            [-1j * k0, -1j * kp / eps, 1j * kp / eps, 0],
            [0, e, em, -1],
            [0, 1j * kp / eps * e, -1j * kp / eps * em, -1j * k0],
        ],
        dtype=np.complex128,
    )
    rhs = np.array([-1, -1j * k0, 0, 0], dtype=np.complex128)
    R, A, B, T = np.linalg.solve(M, rhs)
    return {"R": complex(R), "T": complex(T), "kp": complex(kp), "abs_T2": float(abs(T) ** 2), "abs_R2": float(abs(R) ** 2)}


def fem_slab(k0, eps, L, h=0.02):
    from scipy.spatial import Delaunay

    import sixport_common as sc
    from fem_validated_solver import ElementSampler, assemble_anisotropic, solve_system

    # Long guide so the PML does not sit on the slab. Neumann on y is the 1D reduction.
    nx, ny = 24.0, 1.0
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = nx, ny, 2.0
    xs = np.arange(0.0, nx + 1e-9, h)
    ys = np.arange(0.0, ny + 1e-9, h)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    tris = Delaunay(pts).simplices.astype(int)
    rho_p = 1.0 / eps
    cents = pts[tris].mean(1)
    x0, x1 = 8.0, 8.0 + L
    inside = (cents[:, 0] >= x0) & (cents[:, 0] <= x1)
    T = len(tris)
    rho = [
        np.ones(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.ones(T, dtype=np.complex128),
    ]
    rho[0][inside] = rho_p
    rho[3][inside] = rho_p
    pec = np.zeros(len(pts), dtype=bool)
    A = assemble_anisotropic(pts, tris, *rho, k0, pec)
    sampler = ElementSampler(pts, tris)
    # Uniform-in-y source at x=4, left of the slab.
    src = np.column_stack([np.full(8, 4.0), np.linspace(0.05, ny - 0.05, 8)])
    b = np.zeros(len(pts), dtype=np.complex128)
    for p in src:
        t = int(sampler.locate(p[None, :])[0])
        w = sampler._bary(t, p)
        for k in range(3):
            b[int(tris[t, k])] += complex(w[k])
    x, *_rest = solve_system(A, b, pec)
    # Sample a vertical line left of the slab and right of it, away from PML.
    def line(xv):
        xy = np.column_stack([np.full(20, xv), np.linspace(0.1, ny - 0.1, 20)])
        hz, _, _ = sampler.fields(x, *rho, k0, xy)
        return complex(np.mean(hz))

    left = [line(xv) for xv in (5.0, 5.5, 6.0, 6.5)]
    right = [line(xv) for xv in (12.0, 12.5, 13.0, 13.5)]
    # Incident+reflected on the left: fit Hz = I e^{ikx} + R e^{-ikx} with I known up to scale.
    xl = np.array([5.0, 5.5, 6.0, 6.5])
    zl = np.array(left)
    # Two-wave fit.
    M = np.column_stack([np.exp(1j * k0 * xl), np.exp(-1j * k0 * xl)])
    coef, *_ = np.linalg.lstsq(M, zl, rcond=None)
    I, Rr = coef
    xr = np.array([12.0, 12.5, 13.0, 13.5])
    zr = np.array(right)
    # Transmitted wave referenced to the slab exit x=x1, plus a small left-going residual.
    M2 = np.column_stack([np.exp(1j * k0 * (xr - x1)), np.exp(-1j * k0 * (xr - x1))])
    Tt, leftgoing = np.linalg.lstsq(M2, zr, rcond=None)[0]
    # Analytic T is transmitted / incident amplitude.
    ratio = Tt / I
    return {
        "T_over_I": [ratio.real, ratio.imag],
        "abs_T2": float(abs(ratio) ** 2),
        "leftgoing_over_I": float(abs(leftgoing / I)),
        "R_over_I": [complex(Rr / I).real, complex(Rr / I).imag],
    }


def main():
    import sixport_common as sc
    from fem_validated_solver import eps_tensor_at_bias

    k0 = 2 * np.pi * float(sc.fs_a)
    eps = complex(eps_tensor_at_bias(0.0)[0])
    L = 0.46  # one plasma diameter
    ref = analytic(k0, eps, L)
    num = fem_slab(k0, eps, L)
    rec = {
        "k0": k0,
        "eps": [eps.real, eps.imag],
        "L_a": L,
        "analytic_abs_T2": ref["abs_T2"],
        "analytic_abs_R2": ref["abs_R2"],
        "analytic_T": [ref["T"].real, ref["T"].imag],
        "analytic_R": [ref["R"].real, ref["R"].imag],
        "fem": num,
        "T_abs_dB": float(10 * np.log10(num["abs_T2"] / ref["abs_T2"])),
    }
    print(json.dumps(rec, indent=2))
    (OUT / "plasma_slab.json").write_text(json.dumps(rec, indent=2) + "\n")


if __name__ == "__main__":
    main()
