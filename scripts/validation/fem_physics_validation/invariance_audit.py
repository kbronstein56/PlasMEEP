#!/usr/bin/env python3
"""Numerically irrelevant changes. Observables should match to roundoff or to the mesh error."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from analytic_maxwell import cluster_field, kz_of, solve_clusters  # noqa: E402
from analytic_sweeps import _seven, plasma  # noqa: E402
from fem_validated_solver import assemble_anisotropic  # noqa: E402
from planar_fem import assign_rho, dirichlet, rect_mesh  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation" / "invariance_audit.json"


def guide_problem(pts, tris, phase=1.0 + 0j):
    width, k0 = 1.0, 5.0
    ky = np.pi / width
    beta = kz_of(k0, 1.0, ky)
    x0 = float(pts[:, 0].min())
    y0 = float(pts[:, 1].min())
    length = float(pts[:, 0].max() - x0)
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 8.0, 8.0, 1.0

    def ufn(x, y):
        return phase * np.cos(ky * (y - y0)) * np.exp(1j * beta * (x - x0))

    bound = (np.abs(pts[:, 0] - x0) < 1e-10) | (np.abs(pts[:, 0] - (x0 + length)) < 1e-10)
    rho = assign_rho(pts[tris].mean(1), lambda x, y: 1.0 + 0j)
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    b = np.zeros(len(pts), np.complex128)
    A, b = dirichlet(A, b, bound, ufn(pts[:, 0], pts[:, 1]))
    return A, b, ufn, k0, beta


def rel(a, b):
    return float(np.linalg.norm(a - b) / (np.linalg.norm(b) + 1e-30))


def main():
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 8.0, 8.0, 1.0
    pts, tris, _, _ = rect_mesh(2.2, 3.4, 2.2, 3.2, 0.05)
    A, b, ufn, k0, beta = guide_problem(pts, tris)
    uh_col = splu(A.tocsc(), permc_spec="COLAMD").solve(b)
    uh_nat = splu(A.tocsc(), permc_spec="NATURAL").solve(b)
    exact = ufn(pts[:, 0], pts[:, 1])

    # Opposite diagonal in each quad: (a,b,d)+(a,d,c) versus (a,b,c)+(b,d,c).
    nx = int(round(1.2 / 0.05))
    alt = []
    for j in range(int(round(1.0 / 0.05))):
        for i in range(nx):
            a = j * (nx + 1) + i
            bb = a + 1
            c = a + (nx + 1)
            d = c + 1
            alt.append((a, bb, c))
            alt.append((bb, d, c))
    alt = np.asarray(alt, int)
    A2, b2, _, _, _ = guide_problem(pts, alt)
    uh_alt = splu(A2.tocsc()).solve(b2)
    exact_alt = ufn(pts[:, 0], pts[:, 1])

    # Reverse the last two local nodes. Orientation flips; node ids do not.
    flipped = tris[:, [0, 2, 1]]
    A3, b3, _, _, _ = guide_problem(pts, flipped)
    uh_flip = splu(A3.tocsc()).solve(b3)

    phase = np.exp(1j * 0.7)
    A4, b4, ufn4, _, _ = guide_problem(pts, tris, phase=phase)
    uh_phase = splu(A4.tocsc()).solve(b4)

    # Translate the geometry and the exact mode together. PML box stays large enough.
    shift = np.array([0.37, -0.21])
    pts_s = pts + shift
    A5, b5, ufn5, _, _ = guide_problem(pts_s, tris)
    uh_shift = splu(A5.tocsc()).solve(b5)

    # Mirror through the guide centerline. cos(ky (y-y0)) = cos(ky (y_mirror-y0))
    # because the mode is symmetric and the width is one half-wave. Compare the
    # nodal values after reflecting the y coordinate of the sample.
    y0, y1 = 2.2, 3.2
    pts_m = pts.copy()
    pts_m[:, 1] = y0 + y1 - pts[:, 1]
    A6, b6, _, _, _ = guide_problem(pts_m, tris)
    uh_mirror = splu(A6.tocsc()).solve(b6)

    # Monitor in a homogeneous guide: two x stations of the element-free nodal field.
    x_a, x_b = 2.55, 2.95
    ia = int(np.argmin(np.abs(pts[:, 0] - x_a) + np.abs(pts[:, 1] - 2.45)))
    ib = int(np.argmin(np.abs(pts[:, 0] - x_b) + np.abs(pts[:, 1] - 2.45)))
    ratio_nodes = uh_col[ib] / uh_col[ia]
    ratio_exact = exact[ib] / exact[ia]

    # 60-degree rotation of the production hex. Hz is a scalar, so rotating the
    # source and the sample together leaves Hz unchanged. nmax=6 is enough for
    # a roundoff identity; it is not a new accuracy claim.
    centers = _seven()
    eps = plasma(float(sc.fs_a))
    k0p = 2 * np.pi * float(sc.fs_a)
    src = np.array([-4.5, 0.0])
    sample = np.array([[2.8, 0.4]])
    ang = np.pi / 3.0
    rot = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    coef = solve_clusters(centers, k0p, 0.230, complex(eps), src, 6, coated=([0.230, 0.325, 0.375], [complex(eps), 1.0, 3.8]))
    hz0 = cluster_field(sample, centers, k0p, 0.375, complex(eps), src, coef, 6)
    coef_r = solve_clusters(centers @ rot.T, k0p, 0.230, complex(eps), src @ rot.T, 6, coated=([0.230, 0.325, 0.375], [complex(eps), 1.0, 3.8]))
    hz_r = cluster_field(sample @ rot.T, centers @ rot.T, k0p, 0.375, complex(eps), src @ rot.T, coef_r, 6)
    # Centers of a hex are already invariant, so rotating only the centers is a reordering.
    coef_c = solve_clusters((centers @ rot.T), k0p, 0.230, complex(eps), src, 6, coated=([0.230, 0.325, 0.375], [complex(eps), 1.0, 3.8]))
    hz_same = cluster_field(sample, centers @ rot.T, k0p, 0.375, complex(eps), src, coef_c, 6)

    rec = {
        "guide_dofs": int(len(pts)),
        "colamd_vs_exact": rel(uh_col, exact),
        "natural_vs_colamd": rel(uh_nat, uh_col),
        "opposite_diagonal_vs_exact": rel(uh_alt, exact_alt),
        "opposite_diagonal_vs_first": rel(uh_alt, uh_col),
        "reversed_local_nodes_vs_first": rel(uh_flip, uh_col),
        "source_phase_scale_error": rel(uh_phase, phase * uh_col),
        "translation_vs_unshifted": rel(uh_shift, uh_col),
        "mirror_sign_flip_error": rel(uh_mirror, -uh_col),
        "monitor_point_ratio_abs_diff": abs((uh_col[ia] / exact[ia]) - (uh_col[ib] / exact[ib])),
        "monitor_node_rel_a": abs(uh_col[ia] / exact[ia] - 1.0),
        "monitor_node_rel_b": abs(uh_col[ib] / exact[ib] - 1.0),
        "hex_rotate_source_and_sample_abs_diff": abs(hz_r[0] - hz0[0]),
        "hex_center_reorder_abs_diff": abs(hz_same[0] - hz0[0]),
        "beta": [beta.real, beta.imag],
    }
    OUT.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=2), flush=True)


if __name__ == "__main__":
    main()
