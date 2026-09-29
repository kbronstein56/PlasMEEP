#!/usr/bin/env python3
"""
B=0 diagnosis: compare cosine vs numerical-mode source, and |Hz|² proxy vs Poynting.
Stops if we cannot get within ~1 dB of Meep on important ports.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from fem_validated_solver import (  # noqa: E402
    OUT,
    assemble_anisotropic,
    assign_rho,
    fields_from_hz,
    guide_normal_flux,
    inject_numerical_mode_rhs,
    load_mesh,
    load_numerical_mode,
    log,
    origin_to_mesh,
    pec_mask_horns,
    solve_system,
)
from sixport_common import set_geometry_context, monitor_center_for_port, effective_port_dir  # noqa: E402

MEEP = {
    0: 0.16826699269846146,
    1: 0.0137296820,
    2: 0.00220668062,
    3: 0.318803961,
    4: 0.00220668062,
    5: 0.0137296820,
}


def cosine_rhs(points, pec, port=0):
    set_geometry_context(res=50, horn_walls="prism")
    center = origin_to_mesh(sc.horn_for_port(port, 50)["source_center"])
    u = np.asarray(effective_port_dir(port), dtype=float)
    u /= np.linalg.norm(u)
    tang = np.array([-u[1], u[0]])
    span = 0.96 * sc.clear_width
    b = np.zeros(len(points), dtype=np.complex128)
    for s in np.linspace(-span / 2, span / 2, 41):
        xy = center + s * tang
        d2 = np.sum((points - xy) ** 2, axis=1)
        d2 = d2.copy()
        d2[pec] = np.inf
        i = int(np.argmin(d2))
        amp = max(np.cos(np.pi * s / span), 0.0)
        b[i] += amp
    b /= np.linalg.norm(b) + 1e-30
    return b


def hz2_proxy(points, u, port):
    set_geometry_context(res=50, horn_walls="prism")
    center = origin_to_mesh(monitor_center_for_port(port, 50))
    n_hat = np.asarray(effective_port_dir(port), dtype=float)
    n_hat /= np.linalg.norm(n_hat)
    tang = np.array([-n_hat[1], n_hat[0]])
    span = 0.96 * sc.clear_width
    acc = 0.0
    wsum = 0.0
    for s in np.linspace(-span / 2, span / 2, 41):
        xy = center + s * tang
        i = int(np.argmin(np.sum((points - xy) ** 2, axis=1)))
        w = max(np.cos(np.pi * s / span), 0.0)
        acc += w * float(np.abs(u[i]) ** 2)
        wsum += w
    return float(acc / (wsum + 1e-30))


def db_err(fem, meep):
    if fem <= 0 or meep <= 0:
        return float("nan")
    return 10 * np.log10(fem) - 10 * np.log10(meep)


def main():
    grade = "FEM-M"  # cheaper for diagnosis
    out_dir = OUT / "phase7"
    out_dir.mkdir(parents=True, exist_ok=True)

    points, tris = load_mesh(grade)
    rho = assign_rho(points, tris, 0.0)
    pec = pec_mask_horns(points)
    k0 = 2 * np.pi * sc.fs_a
    log("Assembling...")
    t0 = time.perf_counter()
    A = assemble_anisotropic(points, tris, *rho, k0, pec)
    log(f"  assemble {time.perf_counter()-t0:.2f}s")
    x0, lu, t_fac, _, resid = solve_system(A, cosine_rhs(points, pec), pec)
    log(f"  factor {t_fac:.2f}s resid={resid:.2e}")

    src = load_numerical_mode(OUT / "phase2" / "numerical_mode_P1_res50.json")
    b_num = inject_numerical_mode_rhs(points, src, pec)
    x_num = lu.solve(b_num.copy())
    # zero pec already in solve path for cosine; for num need pec zero
    b_num[pec] = 0
    x_num = lu.solve(b_num)

    cases = {}
    for name, x in (("cosine", x0), ("numerical", x_num)):
        Ex, Ey = fields_from_hz(points, tris, x, *rho, 2 * np.pi * sc.fs_a)
        raw_poynt = {p: guide_normal_flux(points, x, Ex, Ey, p) for p in range(6)}
        raw_hz2 = {p: hz2_proxy(points, x, p) for p in range(6)}
        # Normalize transmissions by |P0_raw| for poynting (entering), and by port0 for hz2
        p0 = raw_poynt[0]
        norm_poynt_enter = {p: raw_poynt[p] / abs(p0) for p in range(6)}
        # For poynting with incident-style: use abs entering as P_inc proxy for transmission only
        # Reflection estimate: not valid without field subtraction
        hz2_0 = raw_hz2[0]
        norm_hz2 = {p: raw_hz2[p] / hz2_0 for p in range(6)}
        cases[name] = {
            "raw_poynting": {str(k): float(v) for k, v in raw_poynt.items()},
            "norm_poynting_by_abs_P0": {str(k): float(v) for k, v in norm_poynt_enter.items()},
            "norm_hz2_by_P0": {str(k): float(v) for k, v in norm_hz2.items()},
            "dB_err_hz2_vs_meep": {
                str(p): float(db_err(norm_hz2[p], MEEP[p])) for p in range(1, 6)
            },
            "dB_err_poynt_absP0_vs_meep": {
                str(p): float(db_err(norm_poynt_enter[p], MEEP[p])) for p in range(1, 6)
            },
        }
        log(f"=== {name} ===")
        log(f"  hz2  P12={norm_hz2[1]:.4f} P14={norm_hz2[3]:.4f}  dBerr12={db_err(norm_hz2[1],MEEP[1]):.2f} dBerr14={db_err(norm_hz2[3],MEEP[3]):.2f}")
        log(f"  poynt P12={norm_poynt_enter[1]:.4f} P14={norm_poynt_enter[3]:.4f} dBerr12={db_err(norm_poynt_enter[1],MEEP[1]):.2f} dBerr14={db_err(norm_poynt_enter[3],MEEP[3]):.2f}")

    # Also try Poynting normalized by FEM vacuum P_inc if cache exists
    inc_path = OUT / "phase3" / "incident_FEM-M.json"
    # build quick incident if missing
    if not inc_path.is_file():
        log("Building FEM-M incident...")
        from fem_validated_solver import fem_incident_reference
        meta = fem_incident_reference(grade, src)
        inc_path.parent.mkdir(parents=True, exist_ok=True)
        inc_path.write_text(json.dumps(meta, indent=2) + "\n")
    else:
        meta = json.loads(inc_path.read_text())

    P_inc = float(meta["P_inc"])
    P_raw_ref = float(meta["P_raw_ref"])
    for name, x in (("cosine", x0), ("numerical", x_num)):
        Ex, Ey = fields_from_hz(points, tris, x, *rho, 2 * np.pi * sc.fs_a)
        raw = {p: guide_normal_flux(points, x, Ex, Ey, p) for p in range(6)}
        # Field-subtraction-free transmission
        T = {p: raw[p] / P_inc for p in range(6)}
        # Scalar subtraction for P11 (known imperfect)
        T[0] = (raw[0] - P_raw_ref) / P_inc
        cases[name]["norm_poynting_by_Pinc"] = {str(k): float(v) for k, v in T.items()}
        cases[name]["dB_err_Pinc_vs_meep"] = {
            str(p): float(db_err(T[p], MEEP[p])) for p in range(6) if T[p] > 0 and MEEP[p] > 0
        }
        log(f"=== {name} / P_inc === P11={T[0]:.4f} P12={T[1]:.4f} P14={T[3]:.4f}")

    summary = {
        "grade": grade,
        "meep": MEEP,
        "incident": meta,
        "cases": cases,
        "interpretation": [
            "If cosine+|Hz2| recovers Meep-like P14, physics geometry is OK and issue is source/Poynting/norm.",
            "If numerical+|Hz2| matches cosine, source profile is not the main issue.",
            "If Poynting ratios differ from |Hz2|, E-field reconstruction or S formula is suspect.",
        ],
    }
    path = out_dir / "b0_source_flux_diagnosis.json"
    path.write_text(json.dumps(summary, indent=2) + "\n")
    log(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
