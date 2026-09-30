#!/usr/bin/env python3
"""Method of manufactured solutions for the P1 Hz weak form.

The assembled operator is

    A = K - k0^2 M

with

    K_ij = ∫ (ρ ∇φ_j) · ∇φ_i dA
    M_ij = ∫ φ_j φ_i dA

when the PML stretch is identically 1. That is the Galerkin form of

    -div(ρ ∇Hz) - k0^2 Hz = f.

The homogeneous problem is the same statement as

    div(ρ ∇Hz) + k0^2 Hz = 0.

The manufactured load uses the first, signed, equation. ρ is element-wise
constant, evaluated at the element centroid. Off-diagonal entries are included
so the ρxy / ρyx index order is tested.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from fem_validated_solver import assemble_anisotropic  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"
K0 = 1.7

# Degree-5 triangle rule. Weights sum to 1; multiply by the element area.
_QUAD = [
    (0.333333333333333, 0.333333333333333, 0.333333333333333, 0.225000000000000),
    (0.059715871789770, 0.470142064105115, 0.470142064105115, 0.132394152788506),
    (0.470142064105115, 0.059715871789770, 0.470142064105115, 0.132394152788506),
    (0.470142064105115, 0.470142064105115, 0.059715871789770, 0.132394152788506),
    (0.797426985353087, 0.101286507323456, 0.101286507323456, 0.125939180544827),
    (0.101286507323456, 0.797426985353087, 0.101286507323456, 0.125939180544827),
    (0.101286507323456, 0.101286507323456, 0.797426985353087, 0.125939180544827),
]


def unit_mesh(n: int, origin: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Right triangles on the unit square, translated into the PML-free interior."""
    xs = origin[0] + np.linspace(0.0, 1.0, n + 1)
    ys = origin[1] + np.linspace(0.0, 1.0, n + 1)
    pts = np.array([(x, y) for y in ys for x in xs], float)
    tris = []
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i
            b = a + 1
            c = a + (n + 1)
            d = c + 1
            tris.append((a, b, d))
            tris.append((a, d, c))
    return pts, np.asarray(tris, int)


def manufactured(xy: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """u and derivatives on the unit-square coordinates (X, Y) = xy - origin."""
    # Caller passes physical xy. The function below is in local coordinates.
    raise RuntimeError("use local_fields")


def local_fields(x: np.ndarray, y: np.ndarray):
    """u = sin(pi x) sin(pi y) + (0.25+0.15j) sin(2 pi x) sin(pi y), zero on the unit square."""
    s1, c1 = np.sin(np.pi * x), np.cos(np.pi * x)
    s2, c2 = np.sin(2 * np.pi * x), np.cos(2 * np.pi * x)
    sy, cy = np.sin(np.pi * y), np.cos(np.pi * y)
    a = 0.25 + 0.15j
    u = s1 * sy + a * s2 * sy
    ux = np.pi * c1 * sy + a * 2 * np.pi * c2 * sy
    uy = np.pi * s1 * cy + a * np.pi * s2 * cy
    uxx = -(np.pi**2) * s1 * sy + a * (-(2 * np.pi) ** 2) * s2 * sy
    uyy = -(np.pi**2) * s1 * sy + a * (-(np.pi**2)) * s2 * sy
    uxy = (np.pi**2) * c1 * cy + a * 2 * (np.pi**2) * c2 * cy
    return u, ux, uy, uxx, uyy, uxy


def forcing(xx, yy, uxx, uyy, uxy, u, rho):
    """f in -div(ρ ∇u) - k0^2 u = f."""
    rxx, rxy, ryx, ryy = rho
    div = rxx * uxx + (rxy + ryx) * uxy + ryy * uyy
    return -div - (K0**2) * u


CASES = {
    "A_real_scalar": (2.0 + 0j, 0j, 0j, 2.0 + 0j),
    "B_complex_scalar": (1.5 + 0.2j, 0j, 0j, 1.5 + 0.2j),
    "C_diagonal_complex": (1.2 + 0.1j, 0j, 0j, 0.8 - 0.05j),
    "D_full_complex": (1.1 + 0.05j, 0.3 + 0.1j, 0.15 - 0.05j, 0.9 + 0.2j),
}


def apply_dirichlet(A, b, bound, values):
    A = A.tolil()
    for i in np.flatnonzero(bound):
        A.rows[i] = [int(i)]
        A.data[i] = [1.0 + 0j]
        b[i] = values[i]
    return A.tocsr(), b


def run_case(name: str, rho, ns: list[int]) -> dict:
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 4.0, 4.0, 1.0
    origin = np.array([1.25, 1.25])
    rows = []
    for n in ns:
        pts, tris = unit_mesh(n, origin)
        local = pts - origin
        bound = (
            (np.abs(local[:, 0]) < 1e-12)
            | (np.abs(local[:, 0] - 1.0) < 1e-12)
            | (np.abs(local[:, 1]) < 1e-12)
            | (np.abs(local[:, 1] - 1.0) < 1e-12)
        )
        T = len(tris)
        rho_e = [np.full(T, rho[k], dtype=np.complex128) for k in range(4)]
        A = assemble_anisotropic(pts, tris, *rho_e, K0, np.zeros(len(pts), dtype=bool))
        # Confirm the manufactured box sits where the PML stretch is 1.
        from fem_validated_solver import pml_sx_sy

        sx, sy = pml_sx_sy(pts[tris].mean(1))
        if np.max(np.abs(sx - 1)) > 1e-12 or np.max(np.abs(sy - 1)) > 1e-12:
            raise RuntimeError(f"PML stretch is active inside the MMS box: {np.max(np.abs(sx-1))}")
        b = np.zeros(len(pts), dtype=np.complex128)
        for t, tri in enumerate(tris):
            xy = pts[tri]
            area = 0.5 * np.abs(
                (xy[1, 0] - xy[0, 0]) * (xy[2, 1] - xy[0, 1]) - (xy[2, 0] - xy[0, 0]) * (xy[1, 1] - xy[0, 1])
            )
            for l1, l2, l3, w in _QUAD:
                p = l1 * xy[0] + l2 * xy[1] + l3 * xy[2]
                loc = p - origin
                u, ux, uy, uxx, uyy, uxy = local_fields(loc[0], loc[1])
                f = forcing(loc[0], loc[1], uxx, uyy, uxy, u, rho)
                for a, lam in enumerate((l1, l2, l3)):
                    b[tri[a]] += area * w * f * lam
        u_exact, *_ = local_fields(local[:, 0], local[:, 1])
        A, b = apply_dirichlet(A, b, bound, u_exact)
        lu = sparse.linalg.splu(A.tocsc())
        uh = lu.solve(b)
        err = uh - u_exact
        # L2 and H1 via element centroid sampling of the P1 interpolant error.
        l2_acc = 0.0
        h1_acc = 0.0
        for tri in tris:
            xy = pts[tri]
            area = 0.5 * np.abs(
                (xy[1, 0] - xy[0, 0]) * (xy[2, 1] - xy[0, 1]) - (xy[2, 0] - xy[0, 0]) * (xy[1, 1] - xy[0, 1])
            )
            acc_u = 0.0
            acc_g = 0.0
            for l1, l2, l3, w in _QUAD:
                p = l1 * xy[0] + l2 * xy[1] + l3 * xy[2]
                loc = p - origin
                u, ux, uy, *_ = local_fields(loc[0], loc[1])
                uh_q = l1 * uh[tri[0]] + l2 * uh[tri[1]] + l3 * uh[tri[2]]
                # P1 gradient is constant on the element.
                twice = (xy[1, 0] - xy[0, 0]) * (xy[2, 1] - xy[0, 1]) - (xy[2, 0] - xy[0, 0]) * (xy[1, 1] - xy[0, 1])
                bx = np.array([xy[1, 1] - xy[2, 1], xy[2, 1] - xy[0, 1], xy[0, 1] - xy[1, 1]]) / twice
                by = np.array([xy[2, 0] - xy[1, 0], xy[0, 0] - xy[2, 0], xy[1, 0] - xy[0, 0]]) / twice
                ghx = bx @ uh[tri]
                ghy = by @ uh[tri]
                du = abs(complex(uh_q) - complex(u))
                acc_u += w * du * du
                acc_g += w * (abs(complex(ghx) - complex(ux)) ** 2 + abs(complex(ghy) - complex(uy)) ** 2)
            l2_acc += area * acc_u
            h1_acc += area * acc_g
        h = 1.0 / n
        rec = {
            "n": n,
            "h": h,
            "dofs": int(len(pts)),
            "L2": float(np.sqrt(l2_acc)),
            "H1": float(np.sqrt(h1_acc)),
            "max_nodal": float(np.max(np.abs(err))),
        }
        if rows:
            rec["L2_order"] = float(np.log(rows[-1]["L2"] / rec["L2"]) / np.log(rows[-1]["h"] / rec["h"]))
            rec["H1_order"] = float(np.log(rows[-1]["H1"] / rec["H1"]) / np.log(rows[-1]["h"] / rec["h"]))
        rows.append(rec)
        print(name, rec, flush=True)
    return {"levels": rows, "rho_reim": [[z.real, z.imag] for z in rho]}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ns = [8, 16, 32, 64]
    out = {}
    for name, rho in CASES.items():
        out[name] = run_case(name, rho, ns)
    (OUT / "mms.json").write_text(json.dumps(out, indent=2) + "\n")
    # Convergence plot.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.8))
    for name, rec in out.items():
        h = [r["h"] for r in rec["levels"]]
        ax[0].loglog(h, [r["L2"] for r in rec["levels"]], "o-", label=name)
        ax[1].loglog(h, [r["H1"] for r in rec["levels"]], "o-", label=name)
    href = np.array([1 / 8, 1 / 64])
    ax[0].loglog(href, 0.15 * href**2, "k--", lw=0.8, label="h^2")
    ax[1].loglog(href, 0.8 * href, "k--", lw=0.8, label="h")
    ax[0].set_xlabel("h")
    ax[1].set_xlabel("h")
    ax[0].set_ylabel("L2 error")
    ax[1].set_ylabel("H1 seminorm error")
    ax[0].invert_xaxis()
    ax[1].invert_xaxis()
    ax[0].legend(fontsize=7)
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "mms_convergence.png", dpi=140)
    print("MMS_DONE", flush=True)


if __name__ == "__main__":
    main()
