#!/usr/bin/env python3
"""Phase 1: matrix-level Onsager audit of FULL91 + modal DtN operator."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "multimode_ports"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from dtn_assembly import add_modal_robin, port_line_integrals  # noqa: E402
from full91_dtn_solve import (  # noqa: E402
    add_all_robin,
    assemble_base,
    build_port_structs,
)
from modes import mode_constants, normalize_on_line  # noqa: E402


def fro_rel(A: sparse.spmatrix, B: sparse.spmatrix) -> float:
    D = (A - B).tocoo()
    num = float(np.sqrt(np.sum(np.abs(D.data) ** 2)))
    den = float(np.sqrt(np.sum(np.abs(A.tocoo().data) ** 2)))
    return num / max(den, 1e-30)


def complex_sym_rel(A: sparse.spmatrix) -> float:
    """||A - A^T|| / ||A|| (complex transpose, not Hermitian)."""
    return fro_rel(A, A.T)


def port_blocks(port_structs, modes, n_nodes: int, n_prop=2, n_evan=2):
    """Build separate sparse contributions: TE0, TE1, evanescent."""
    W = float(sc.clear_width)
    A0 = sparse.csr_matrix((n_nodes, n_nodes), dtype=np.complex128)
    A1 = sparse.csr_matrix((n_nodes, n_nodes), dtype=np.complex128)
    Aev = sparse.csr_matrix((n_nodes, n_nodes), dtype=np.complex128)
    for ps in port_structs:
        if ps["n_nodes"] < 3:
            continue
        for md in modes[: n_prop + n_evan]:
            m = md["m"]
            phi, _ = normalize_on_line(ps["y"], ps["w"], m, W)
            if m == 0:
                A0 = add_modal_robin(A0, ps["nodes"], ps["y"], ps["w"], phi, md["beta"])
            elif m == 1:
                A1 = add_modal_robin(A1, ps["nodes"], ps["y"], ps["w"], phi, md["beta"])
            else:
                Aev = add_modal_robin(Aev, ps["nodes"], ps["y"], ps["w"], phi, md["beta"])
    return A0, A1, Aev


def audit_port_geometry(port_structs, points) -> dict:
    """Check whether mesh_port_nodes collects a node cloud vs a 1-D boundary."""
    rows = []
    for ps in port_structs:
        idx = np.asarray(ps["nodes"], int)
        y = np.asarray(ps["y"], float)
        n = np.asarray(ps["outward"], float)
        n = n / np.linalg.norm(n)
        c = np.asarray(ps["center"], float)
        along = (points[idx] - c) @ n
        # unique-ish along bins
        n_unique_y = len(np.unique(np.round(y, 8)))
        span_y = float(y[-1] - y[0]) if len(y) > 1 else 0.0
        rows.append({
            "port": int(ps["port"]),
            "n_nodes": int(len(idx)),
            "n_unique_y_8dec": int(n_unique_y),
            "y_span": span_y,
            "W_clear": float(sc.clear_width),
            "along_min": float(along.min()) if len(along) else None,
            "along_max": float(along.max()) if len(along) else None,
            "along_std": float(along.std()) if len(along) else None,
            "n_distinct_along_layers": int(len(np.unique(np.round(along, 4)))),
            "weight_sum": float(np.sum(ps["w"])),
        })
    return {"ports": rows, "plane_tol_used": 0.06}


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    s = np.ones(91)
    report = {"grade": args.grade, "f_Hz": 3.85e9}

    # Geometry audit (B-independent)
    pts, tris, rho, A_bulk_p, k0, _, meta = assemble_base(args.grade, 3.85e9, 0.05, s)
    ps = meta["port_structs"]
    report["geometry"] = audit_port_geometry(ps, pts)
    report["n_exterior_tris"] = meta["n_exterior_tris"]
    report["k0"] = float(k0)
    report["n_nodes"] = int(len(pts))

    modes = mode_constants(k0, float(sc.clear_width), n_modes=4)
    A_te0, A_te1, A_ev = port_blocks(ps, modes, len(pts))

    report["port_symmetry"] = {
        "TE0": complex_sym_rel(A_te0),
        "TE1": complex_sym_rel(A_te1),
        "evanescent": complex_sym_rel(A_ev),
        "note": "||A-A^T||/||A||; reciprocal Robin must be complex-symmetric",
    }

    # Bulk at ±B, 0
    _, _, _, A_bulk_m, _, _, _ = assemble_base(args.grade, 3.85e9, -0.05, s)
    _, _, _, A_bulk_0, _, _, _ = assemble_base(args.grade, 3.85e9, 0.0, s)

    report["bulk"] = {
        "r_A_plus_minus": fro_rel(A_bulk_p, A_bulk_m.T),
        "r_A_zero": complex_sym_rel(A_bulk_0),
    }

    # Full A = bulk + all Robin
    A_full_p = add_all_robin(A_bulk_p.copy(), ps, modes)
    A_full_m = add_all_robin(A_bulk_m.copy(), ps, modes)
    A_full_0 = add_all_robin(A_bulk_0.copy(), ps, modes)

    report["full_with_dtn"] = {
        "r_A_plus_minus": fro_rel(A_full_p, A_full_m.T),
        "r_A_zero": complex_sym_rel(A_full_0),
        "r_port_only_zero": complex_sym_rel(A_te0 + A_te1 + A_ev),
    }

    # Discrete c vectors: check real
    c_imag = []
    W = float(sc.clear_width)
    for pstruct in ps:
        for m in (0, 1):
            phi, _ = normalize_on_line(pstruct["y"], pstruct["w"], m, W)
            c = port_line_integrals(pstruct["nodes"], pstruct["y"], pstruct["w"], phi)
            c_imag.append(float(np.max(np.abs(c.imag))))
    report["max_abs_imag_c"] = float(max(c_imag) if c_imag else 0.0)

    bar = 1e-10  # matrix-level should be near machine / assembly noise
    report["DTN_MATRIX_ONSAGER"] = (
        "PASS"
        if (
            report["full_with_dtn"]["r_A_plus_minus"] < 1e-8
            and report["full_with_dtn"]["r_A_zero"] < 1e-8
            and report["port_symmetry"]["TE0"] < 1e-12
            and report["port_symmetry"]["TE1"] < 1e-12
        )
        else "FAIL"
    )
    report["bars"] = {"full_A_rel": 1e-8, "port_sym": 1e-12}

    (OUT / "dtn_matrix_onsager_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["DTN_MATRIX_ONSAGER"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
