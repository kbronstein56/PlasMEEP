#!/usr/bin/env python3
"""Export real complex Hz fields for presentation clips 02/03.

Single forward solves of EXISTING designs only — no optimization.
Writes compact regular-grid NPZ + meta.json under presentation_fields/.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
FIELD_ROOT = OUT / "presentation_fields"
sys_paths = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "jesse_method_campaign"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sys

sys.path[:0] = sys_paths

from adjoint_direct_dtn import build_system  # noqa: E402
from a2_efficient_core import rss_gib  # noqa: E402
from direct_dtn_s import build_Q_columns  # noqa: E402
from efficient_adjoint import S_from_U, _pack_C, _solve_multi  # noqa: E402
from multimode_opt_core import (  # noqa: E402
    B_PLUS,
    expand_tied,
    orbits_and_masks,
    te0_metrics_from_S,
)
from matplotlib.tri import LinearTriInterpolator, Triangulation  # noqa: E402


def _design_q(tag: str) -> tuple[np.ndarray, dict]:
    orbits, _, _ = orbits_and_masks()
    if tag == "uniform":
        q = np.ones(len(orbits), float)
        meta = {"design_id": "uniform_q=1", "source_file": None}
        return q, meta
    if tag == "stage2_best":
        p = OUT / "opt_broadband_stage2_A2_C.json"
        d = json.loads(p.read_text())
        q = np.asarray(d["q_best"], float)
        meta = {
            "design_id": "stage2_A2_C_q_best",
            "source_file": str(p),
            "J_BB": d.get("J_BB_best", d.get("J_BB")),
            "worst_f_GHz": d.get("worst_f_GHz"),
        }
        return q, meta
    if tag == "stage3_best":
        p = OUT / "opt_broadband_stage3_A2_C.json"
        d = json.loads(p.read_text())
        q = np.asarray(d["q_best"], float)
        meta = {
            "design_id": "stage3_A2_C_q_best",
            "source_file": str(p),
            "J_BB": d.get("J_BB_best", d.get("J_BB")),
            "accepted_nit": d.get("accepted_nit"),
        }
        return q, meta
    if tag == "stage4_best":
        p = OUT / "opt_broadband_stage4_A2_C.json"
        d = json.loads(p.read_text())
        q = np.asarray(d["q_best"], float)
        meta = {
            "design_id": "stage4_A2_C_q_best",
            "source_file": str(p),
            "J_BB": d.get("J_BB_best", d.get("J_BB")),
            "accepted_nit": d.get("accepted_nit_total", d.get("accepted_nit")),
        }
        return q, meta
    raise ValueError(tag)


def interpolate_to_grid(points, tris, Hz, nx=480, ny=448):
    """Linear FEM nodal → regular Cartesian grid for imshow movies."""
    points = np.asarray(points, float)
    tris = np.asarray(tris, int)
    Hz = np.asarray(Hz, complex)
    xmin, xmax = float(points[:, 0].min()), float(points[:, 0].max())
    ymin, ymax = float(points[:, 1].min()), float(points[:, 1].max())
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    X, Y = np.meshgrid(x, y)
    tri = Triangulation(points[:, 0], points[:, 1], tris)
    re_i = LinearTriInterpolator(tri, np.real(Hz))
    im_i = LinearTriInterpolator(tri, np.imag(Hz))
    Re = np.ma.filled(re_i(X, Y), np.nan)
    Im = np.ma.filled(im_i(X, Y), np.nan)
    return x, y, Re.astype(np.float32), Im.astype(np.float32)


def export_one(
    grade: str,
    design_tag: str,
    f_ghz: float,
    source_port: int = 1,
    nx: int = 480,
    ny: int = 448,
    b_tesla: float = B_PLUS,
):
    """Forward-only TE0 solve for one source port; write viz-ready artifacts.

    Ports are 1-indexed in metadata (P1..P6); internally 0-indexed.
    Convention: animate Re[Hz * exp(-i ω t)], ω=2πf.
    """
    assert 1 <= source_port <= 6
    q, dmeta = _design_q(design_tag)
    orbits, _, _ = orbits_and_masks()
    f_hz = float(f_ghz) * 1e9
    t0 = time.perf_counter()
    rss0 = rss_gib()

    s = expand_tied(q, orbits)
    st = build_system(grade, s, f_hz, b_tesla, "projected")
    points, A = st["points"], st["A"]
    cols = build_Q_columns(points, st["ports"], st["modes"], n_prop=2)
    C, betas = _pack_C(points, cols)
    te0_cols = [2 * j for j in range(6)]
    j_src = source_port - 1  # P1 → column 0 of TE0 drives

    t_lu0 = time.perf_counter()
    lu = splu(A.tocsc())
    t_factor = time.perf_counter() - t_lu0
    # Solve all 6 TE0 for S metrics, but only export source-port field
    t_s0 = time.perf_counter()
    U = _solve_multi(lu, C[:, te0_cols], trans="N")
    t_solve = time.perf_counter() - t_s0
    S = S_from_U(C, U, betas, te0_cols)
    met = te0_metrics_from_S(S)
    Hz = np.asarray(U[:, j_src], complex)

    x, y, Re, Im = interpolate_to_grid(points, st["tris"], Hz, nx=nx, ny=ny)
    wall = time.perf_counter() - t0

    # Use 4 decimals so 5.8125 does not collide with 5.812 under :.3f rounding
    stem = f"{design_tag}_P{source_port}_{f_ghz:.4f}GHz_Bplus".replace(".", "p")
    out_dir = FIELD_ROOT / stem
    out_dir.mkdir(parents=True, exist_ok=True)

    np.save(out_dir / "x.npy", x)
    np.save(out_dir / "y.npy", y)
    np.save(out_dir / "Hz_real.npy", Re)
    np.save(out_dir / "Hz_imag.npy", Im)
    # Compact archive + nodal subset for reproducibility audits
    np.savez_compressed(
        out_dir / "field_grid.npz",
        x=x,
        y=y,
        Hz_real=Re,
        Hz_imag=Im,
    )
    # Save a downsampled nodal dump (every Nth node) for provenance, not full 823k
    stride = max(1, len(points) // 50000)
    np.savez_compressed(
        out_dir / "field_nodal_stride.npz",
        points=np.asarray(points[::stride], float),
        Hz_real=np.real(Hz[::stride]).astype(np.float32),
        Hz_imag=np.imag(Hz[::stride]).astype(np.float32),
        stride=stride,
        n_dof_full=len(points),
    )

    meta = {
        "design_tag": design_tag,
        "design": dmeta,
        "f_GHz": float(f_ghz),
        "f_Hz": float(f_hz),
        "omega_rad_s": float(2.0 * np.pi * f_hz),
        "source_port": int(source_port),
        "source_port_label": f"P{source_port}",
        "port_indexing": "P1..P6 = indices 0..5; desired for B_PLUS is (j+1)%6",
        "B_tesla": float(b_tesla),
        "B_sign": "B_PLUS",
        "mesh_grade": grade,
        "n_dof": int(A.shape[0]),
        "n_tris": int(len(st["tris"])),
        "grid_nx": int(nx),
        "grid_ny": int(ny),
        "xy_units": "cm (device domain as in production FEM)",
        "field_quantity": "Hz (out-of-plane H, complex phasor at nodes → linear interp to grid)",
        "normalization": (
            "As returned by production Direct-DtN TE0 drive column; "
            "absolute scale is solver-consistent within a movie, not absolute SI calib."
        ),
        "phase_time_convention": (
            "Phasor convention for animation: frame(t) = Re[Hz * exp(-1j * omega * t)], "
            "omega=2*pi*f. DO NOT animate |Hz| as oscillating."
        ),
        "solver": "scipy.sparse.linalg.splu (SuperLU) via production adjoint_direct_dtn",
        "metrics_te0": {k: (float(v) if np.ndim(v) == 0 else v) for k, v in met.items()
                        if k in ("T_des_avg", "T_rev_avg", "IL_avg_dB", "isolation_avg_dB",
                                 "RL_avg_dB", "R_avg", "TE1_avg")},
        "timing_s": {
            "wall_total": wall,
            "factor_lu": t_factor,
            "forward_6rhs": t_solve,
            "rss_GiB_start": rss0,
            "rss_GiB_end": rss_gib(),
        },
        "geometry_reference": (
            "Same FULL91 / grade-C mesh and port layout as Stage-2 production "
            "(orbits_and_masks / build_system projected)."
        ),
    }
    # JSON-safe metrics
    for k, v in list(meta["metrics_te0"].items()):
        try:
            meta["metrics_te0"][k] = float(v)
        except Exception:
            meta["metrics_te0"][k] = str(v)

    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(
        f"EXPORT {stem}: wall={wall:.1f}s factor={t_factor:.1f}s "
        f"IL={meta['metrics_te0'].get('IL_avg_dB')} iso={meta['metrics_te0'].get('isolation_avg_dB')}",
        flush=True,
    )
    return out_dir, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--nx", type=int, default=480)
    ap.add_argument("--ny", type=int, default=448)
    ap.add_argument(
        "--jobs",
        default="all",
        help=(
            "Comma list or 'all' (legacy Stage-2 set) or 'stage3': "
            "stage3@4.1/5.75/5.825 + stage2@5.825 for hole-repair movies"
        ),
    )
    args = ap.parse_args()
    FIELD_ROOT.mkdir(parents=True, exist_ok=True)

    catalog = {
        "uniform_4p1": ("uniform", 4.1),
        "stage2_4p1": ("stage2_best", 4.1),
        "stage2_5p65": ("stage2_best", 5.65),
        "stage2_5p75": ("stage2_best", 5.75),
        "stage2_5p825": ("stage2_best", 5.825),
        "stage3_4p1": ("stage3_best", 4.1),
        "stage3_5p75": ("stage3_best", 5.75),
        "stage3_5p825": ("stage3_best", 5.825),
        "stage3_4p3": ("stage3_best", 4.3),
        "stage3_4p775": ("stage3_best", 4.775),
        "stage3_5p8125": ("stage3_best", 5.8125),
        "stage4_4p775": ("stage4_best", 4.775),
        "stage4_5p8125": ("stage4_best", 5.8125),
        "stage4_5p08125": ("stage4_best", 5.08125),
        "stage4_4p725": ("stage4_best", 4.725),
        "stage3_5p08125": ("stage3_best", 5.08125),
        "stage3_4p725": ("stage3_best", 4.725),
    }
    if args.jobs.strip() == "all":
        keys = [
            "uniform_4p1",
            "stage2_4p1",
            "stage2_5p65",
            "stage2_5p75",
        ]
    elif args.jobs.strip() == "stage3":
        keys = [
            "stage3_4p1",
            "stage3_5p75",
            "stage3_5p825",
            "stage3_4p3",
            "stage2_5p825",
        ]
    elif args.jobs.strip() == "stage4":
        keys = [
            "stage4_4p775",
            "stage4_5p8125",
            "stage4_5p08125",
            "stage4_4p725",
            "stage3_5p08125",
            "stage3_4p725",
            "stage3_4p775",
            "stage3_5p8125",
        ]
    else:
        keys = [k.strip() for k in args.jobs.split(",") if k.strip()]

    # Merge into existing INDEX rather than wiping prior exports
    index_path = FIELD_ROOT / "INDEX.json"
    if index_path.exists():
        index = json.loads(index_path.read_text())
        by_key = {e["key"]: e for e in index.get("exports", []) if "key" in e}
    else:
        index = {"exports": [], "field_root": str(FIELD_ROOT)}
        by_key = {}
    for key in keys:
        if key not in catalog:
            raise SystemExit(f"unknown job {key}; choose from {list(catalog)}")
        tag, f_ghz = catalog[key]
        # Skip if export directory already present
        stem = f"{tag}_P1_{f_ghz:.4f}GHz_Bplus".replace(".", "p")
        if (FIELD_ROOT / stem).exists() and (FIELD_ROOT / stem / "meta.json").exists():
            print(f"skip existing {stem}", flush=True)
            meta = json.loads((FIELD_ROOT / stem / "meta.json").read_text())
            by_key[key] = {
                "key": key,
                "path": str(FIELD_ROOT / stem),
                "design_tag": tag,
                "f_GHz": f_ghz,
                "source_port": 1,
                "meta": meta,
            }
            continue
        out_dir, meta = export_one(
            args.grade, tag, f_ghz, source_port=1, nx=args.nx, ny=args.ny
        )
        by_key[key] = {
            "key": key,
            "path": str(out_dir),
            "design_tag": tag,
            "f_GHz": f_ghz,
            "source_port": 1,
            "meta": meta,
        }
    index["exports"] = list(by_key.values())
    index["field_root"] = str(FIELD_ROOT)
    index_path.write_text(json.dumps(index, indent=2) + "\n")
    print("wrote", index_path, flush=True)


if __name__ == "__main__":
    # Match common production launch threading unless overridden
    for k, v in (
        ("OMP_NUM_THREADS", "8"),
        ("OPENBLAS_NUM_THREADS", "8"),
        ("MKL_NUM_THREADS", "8"),
    ):
        os.environ.setdefault(k, v)
    main()
