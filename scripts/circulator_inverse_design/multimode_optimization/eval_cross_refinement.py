#!/usr/bin/env python3
"""Fixed-design crossed-refinement study for opt_uniform_C_warm_M.

Evaluates the SAME s on a sequence of FULL91 grades. No reoptimization.
Records FoMs, power balance, and factorization resource usage.
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "multimode_optimization"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from adjoint_direct_dtn import build_system  # noqa: E402
from direct_dtn_s import build_Q_columns, direct_S_from_factor  # noqa: E402
from full91_geometry import GRADES  # noqa: E402
from multimode_opt_core import (  # noqa: E402
    B_PLUS,
    F_HZ,
    dump_json,
    te0_absorption_mean,
    te0_metrics_from_S,
)

# Frozen bars — do not redefine after seeing results
BARS = {
    "J": 0.05,
    "T_des_avg": 0.04,
    "IL_avg_dB": 0.40,
    "isolation_avg_dB": 1.5,
    "RL_avg_dB": 1.5,
    "TE1_avg": 0.015,
    "T_thr_avg": 0.010,
    "T_other_avg": 0.030,
}

# Default matrix: M-centric + F + F-centric factorable references
DEFAULT_GRADES = ["M", "Mh", "Ma", "Mi", "Mia", "F", "R", "A", "Q", "X2"]
SKIP_GRADES = {"S", "X"}  # SuperLU abort history


def rss_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0**2)


def free_gib() -> float:
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("MemAvailable:"):
            return float(line.split()[1]) / (1024.0**2)
    return float("nan")


def eval_grade(grade: str, s: np.ndarray, f_hz: float, b_tesla: float) -> dict:
    if grade in SKIP_GRADES:
        return {"grade": grade, "skipped": True, "reason": "SuperLU abort history"}
    mesh = ROOT / "outputs" / "validation" / "fem_full_device_validation" / "meshes" / f"full91_{grade}.npz"
    if not mesh.exists():
        return {"grade": grade, "skipped": True, "reason": f"missing mesh {mesh}"}

    data = np.load(mesh)
    n_nodes = int(len(data["points"]))
    n_tris = int(len(data["triangles"]))
    est_lu = 8.0e-6 * n_nodes
    cfg = GRADES.get(grade, {})
    free0 = free_gib()
    if est_lu > 0.55 * free0:
        return {
            "grade": grade,
            "skipped": True,
            "reason": f"est_LU {est_lu:.1f} GiB too close to free {free0:.1f} GiB",
            "nodes": n_nodes,
            "est_LU_GiB": est_lu,
        }

    print(
        f"=== grade {grade} nodes={n_nodes} est_LU={est_lu:.1f}GiB free={free0:.1f}GiB "
        f"cfg={cfg} ===",
        flush=True,
    )
    t0 = time.perf_counter()
    rss0 = rss_gib()
    st = build_system(grade, s, f_hz, b_tesla, measure="projected")
    t_asm = time.perf_counter() - t0
    A = st["A"]
    nnz = int(A.nnz)

    t1 = time.perf_counter()
    lu = splu(A.tocsc())
    t_lu = time.perf_counter() - t1
    rss1 = rss_gib()

    t2 = time.perf_counter()
    cols = build_Q_columns(st["points"], st["ports"], st["modes"], n_prop=2)
    S = direct_S_from_factor(lu, len(st["points"]), cols)
    met = te0_metrics_from_S(S)
    P_mat = te0_absorption_mean(st, lu)
    met["P_material_avg"] = P_mat
    met["rem_avg"] = float(1.0 - met["P_out_avg"] - P_mat)
    t_post = time.perf_counter() - t2

    # Persist S for later field/energy diagnostics if needed
    np.savez_compressed(
        OUT / f"S12_crossref_{grade}_Bp.npz",
        S=S,
        s=s,
        grade=grade,
        f_Hz=f_hz,
        B=b_tesla,
    )

    row = {
        "grade": grade,
        "skipped": False,
        "cfg": cfg,
        "nodes": n_nodes,
        "triangles": n_tris,
        "nnz": nnz,
        "est_LU_GiB": est_lu,
        "free_GiB_before": free0,
        "rss_GiB_before_factor": rss0,
        "rss_GiB_after_factor": rss1,
        "rss_delta_GiB": rss1 - rss0,
        "assemble_s": t_asm,
        "factor_s": t_lu,
        "post_s": t_post,
        "elapsed_s": time.perf_counter() - t0,
        "J": float(met["J"]),
        "T_des_avg": float(met["T_des_avg"]),
        "IL_avg_dB": float(met["IL_avg_dB"]),
        "isolation_avg_dB": float(met["isolation_avg_dB"]),
        "RL_avg_dB": float(met["RL_avg_dB"]),
        "TE1_avg": float(met["TE1_avg"]),
        "T_thr_avg": float(met["T_thr_avg"]),
        "T_other_avg": float(met["T_other_avg"]),
        "R_avg": float(met["R_avg"]),
        "P_out_avg": float(met["P_out_avg"]),
        "P_material_avg": float(P_mat),
        "power_residual": abs(float(met["rem_avg"])),
        "metrics": met,
    }
    print(
        json.dumps(
            {
                "grade": grade,
                "J": row["J"],
                "T_des": row["T_des_avg"],
                "iso": row["isolation_avg_dB"],
                "RL": row["RL_avg_dB"],
                "TE1": row["TE1_avg"],
                "rem": row["power_residual"],
                "factor_s": t_lu,
                "rss_delta_GiB": row["rss_delta_GiB"],
            },
            indent=2,
        ),
        flush=True,
    )
    return row


def classify_vs_ref(rows: dict, ref: str = "F") -> dict:
    """Classify each metric vs frozen bars relative to reference grade."""
    if ref not in rows or rows[ref].get("skipped"):
        return {}
    r0 = rows[ref]
    out = {}
    for g, r in rows.items():
        if r.get("skipped") or g == ref:
            continue
        per = {}
        for k, bar in BARS.items():
            if k == "J":
                d = abs(r["J"] - r0["J"])
            else:
                d = abs(r[k] - r0[k])
            per[k] = {
                "delta": float(d if k == "J" else (r[k] - r0[k])),
                "abs_delta": float(d),
                "bar": bar,
                "within_bar": bool(d <= bar),
            }
        out[g] = {
            "vs": ref,
            "all_within_bar": all(v["within_bar"] for v in per.values()),
            "metrics": per,
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--opt-json",
        default="outputs/circulator_inverse_design/multimode_optimization/opt_uniform_C_warm_M.json",
    )
    ap.add_argument("--grades", default=",".join(DEFAULT_GRADES))
    ap.add_argument("--resume", action="store_true", help="Skip grades already in JSON")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    d = json.loads(Path(args.opt_json).read_text())
    s = np.asarray(d["best"]["s"], float)

    out_path = OUT / "cross_refinement_metrics.json"
    report = {
        "design": args.opt_json,
        "f_Hz": F_HZ,
        "B": B_PLUS,
        "bars": BARS,
        "grades": {},
    }
    if args.resume and out_path.exists():
        report = json.loads(out_path.read_text())
        report.setdefault("grades", {})

    for g in [x.strip() for x in args.grades.split(",") if x.strip()]:
        if args.resume and g in report["grades"] and not report["grades"][g].get("skipped"):
            if "J" in report["grades"][g]:
                print(f"skip existing {g}", flush=True)
                continue
        try:
            report["grades"][g] = eval_grade(g, s, F_HZ, B_PLUS)
        except Exception as e:
            report["grades"][g] = {"grade": g, "skipped": True, "reason": repr(e), "error": True}
            print(f"FAILED {g}: {e}", flush=True)
        dump_json(out_path, report)  # checkpoint after each grade

    report["vs_F"] = classify_vs_ref(report["grades"], "F")
    if "M" in report["grades"] and not report["grades"]["M"].get("skipped"):
        report["vs_M"] = classify_vs_ref(report["grades"], "M")
    dump_json(out_path, report)
    print("wrote", out_path)


if __name__ == "__main__":
    main()
