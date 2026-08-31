#!/usr/bin/env python3
"""
Rotated-horn / Cartesian-grid reciprocity diagnostics (horns_only).

Runs fractional-cell offset sweeps, resolution checks, and wall-model
comparisons via reciprocity_b0_study.py. Writes JSON index + diagnostic figure.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
VAL = ROOT / "scripts" / "validation"
OUT = ROOT / "outputs" / "validation" / "horn_grid"
PY = os.environ.get(
    "PLASMEEP_PYTHON",
    "/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/python",
)
MPIRUN = os.environ.get(
    "PLASMEEP_MPIRUN",
    "/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/mpirun",
)


def _mpi_env() -> dict:
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("FI_PROVIDER", "tcp")
    env.setdefault("MPICH_CH4_NETMOD", "ofi")
    env.setdefault("UCX_TLS", "tcp,self")
    env["PYTHONPATH"] = f"{VAL}:{ROOT / 'scripts'}:" + env.get("PYTHONPATH", "")
    return env


def _cell_tag(v: float) -> str:
    s = f"{v:g}"
    return s.replace("-", "neg").replace(".", "p")


def _run_pair(
    *,
    label: str,
    ports: str,
    res: int,
    run_time: float,
    ranks: int,
    grid_offset: Tuple[float, float] = (0.0, 0.0),
    monitor_offset: Tuple[float, float] = (0.0, 0.0),
    wall_mode: str = "pec",
    formulation: str = "te1_hz_line",
    horn_walls: str = "prism",
    coord_rotation_deg: float = 0.0,
) -> Dict[str, Any]:
    OUT.mkdir(parents=True, exist_ok=True)
    json_out = OUT / f"{label}.json"
    log_out = OUT / f"{label}.log"
    if json_out.is_file():
        with open(json_out, encoding="utf-8") as f:
            data = json.load(f)
        return {
            "label": label,
            "skipped": True,
            "path": str(json_out),
            "P12_dB": data.get("reciprocity", {}).get("max_abs_diff_dB"),
            "ports": data.get("settings", {}).get("ports"),
        }

    ox, oy = grid_offset
    mx, my = monitor_offset
    cmd = [
        MPIRUN,
        "-np",
        str(ranks),
        PY,
        str(VAL / "reciprocity_b0_study.py"),
        "--res",
        str(res),
        "--run-time",
        str(run_time),
        "--ports",
        ports,
        "--port-formulation",
        formulation,
        "--device-mode",
        "horns_only",
        f"--grid-offset-cells={ox},{oy}",
        f"--monitor-offset-cells={mx},{my}",
        "--wall-mode",
        wall_mode,
        "--horn-walls",
        horn_walls,
        "--coord-rotation-deg",
        str(coord_rotation_deg),
        "--label",
        label,
        "--json-out",
        str(json_out),
        "--skip-norm-if-cached",
    ]
    t0 = time.time()
    with open(log_out, "w", encoding="utf-8") as log:
        log.write(f"# {' '.join(cmd)}\n")
        proc = subprocess.run(
            cmd, cwd=str(ROOT), env=_mpi_env(), stdout=log, stderr=subprocess.STDOUT
        )
    row: Dict[str, Any] = {
        "label": label,
        "ports": ports,
        "res": res,
        "run_time": run_time,
        "grid_offset_cells": list(grid_offset),
        "monitor_offset_cells": list(monitor_offset),
        "wall_mode": wall_mode,
        "horn_walls": horn_walls,
        "coord_rotation_deg": coord_rotation_deg,
        "wall_s": time.time() - t0,
        "returncode": proc.returncode,
        "path": str(json_out),
    }
    if json_out.is_file():
        with open(json_out, encoding="utf-8") as f:
            data = json.load(f)
        row["P12_dB"] = data.get("reciprocity", {}).get("max_abs_diff_dB")
        row["incident_rel"] = (data.get("incident_power_mismatch") or {}).get(
            "rel_spread"
        )
    return row


def offset_sweep(
    res: int,
    run_time: float,
    ranks: int,
    steps: List[float],
    ports: str = "0,1",
    horn_walls: str = "prism",
    coord_rotation_deg: float = 0.0,
    label_prefix: str = "offset",
) -> List[Dict[str, Any]]:
    rows = []
    for ox in steps:
        for oy in steps:
            tag = f"{label_prefix}_res{res}_ox{_cell_tag(ox)}_oy{_cell_tag(oy)}_p{ports.replace(',', '')}"
            print(f"\n=== {tag} ===", flush=True)
            row = _run_pair(
                label=tag,
                ports=ports,
                res=res,
                run_time=run_time,
                ranks=ranks,
                grid_offset=(ox, oy),
                horn_walls=horn_walls,
                coord_rotation_deg=coord_rotation_deg,
            )
            rows.append(row)
            print(f"  P12={row.get('P12_dB')} wall={row.get('wall_s', 0):.1f}s", flush=True)
    return rows


CANDIDATES = [
    ("prism_r0", "prism", 0.0),
    ("rotated_blocks_r0", "rotated_blocks", 0.0),
    ("grid_snapped_r0", "grid_snapped_prism", 0.0),
    ("prism_r15", "prism", 15.0),
    ("prism_r30", "prism", 30.0),
    ("rotated_r15", "rotated_blocks", 15.0),
]


def candidate_offset_screen(
    res: int,
    run_time: float,
    ranks: int,
    steps: List[float],
) -> List[Dict[str, Any]]:
    """Screen horn-wall candidates by fractional-cell offset sensitivity."""
    summaries: List[Dict[str, Any]] = []
    for name, walls, rot in CANDIDATES:
        print(f"\n######## candidate {name} ########", flush=True)
        rows = offset_sweep(
            res,
            run_time,
            ranks,
            steps,
            horn_walls=walls,
            coord_rotation_deg=rot,
            label_prefix=f"cand_{name}",
        )
        vals = [r["P12_dB"] for r in rows if r.get("P12_dB") is not None]
        zero = [
            r["P12_dB"]
            for r in rows
            if r.get("grid_offset_cells") == [0.0, 0.0]
        ]
        summary = {
            "candidate": name,
            "horn_walls": walls,
            "coord_rotation_deg": rot,
            "res": res,
            "offset_range_dB": (max(vals) - min(vals)) if vals else None,
            "offset_std_dB": float(np.std(vals)) if vals else None,
            "offset_max_dB": max(vals) if vals else None,
            "offset_min_dB": min(vals) if vals else None,
            "zero_offset_dB": zero[0] if zero else None,
            "runs": rows,
        }
        summaries.append(summary)
        print(
            f"  {name}: zero={summary['zero_offset_dB']:.3f} "
            f"range={summary['offset_range_dB']:.3f} "
            f"std={summary['offset_std_dB']:.3f}",
            flush=True,
        )
    out = OUT / "candidate_screen.json"
    with open(out, "w", encoding="utf-8") as f:
        json.dump({"candidates": summaries}, f, indent=2)
    print(f"Wrote {out}")
    return summaries


def resolution_sweep(
    resolutions: List[int],
    run_time: float,
    ranks: int,
    horn_walls: str = "prism",
    coord_rotation_deg: float = 0.0,
    label_prefix: str = "res",
) -> List[Dict[str, Any]]:
    rows = []
    for res in resolutions:
        for ports, tag_ports in [("0,1", "P1P2"), ("1,2", "P2P3")]:
            label = f"{label_prefix}{res}_{tag_ports}_zero_offset"
            print(f"\n=== {label} ===", flush=True)
            row = _run_pair(
                label=label,
                ports=ports,
                res=res,
                run_time=run_time,
                ranks=ranks,
                horn_walls=horn_walls,
                coord_rotation_deg=coord_rotation_deg,
            )
            rows.append(row)
            print(f"  P12={row.get('P12_dB')}", flush=True)
    return rows


def wall_mode_compare(res: int, run_time: float, ranks: int) -> List[Dict[str, Any]]:
    rows = []
    for wall in ("pec", "high_eps"):
        label = f"wall_{wall}_res{res}_P1P2"
        print(f"\n=== {label} ===", flush=True)
        row = _run_pair(
            label=label,
            ports="0,1",
            res=res,
            run_time=run_time,
            ranks=ranks,
            wall_mode=wall,
        )
        rows.append(row)
    return rows


def monitor_offset_sweep(
    res: int,
    run_time: float,
    ranks: int,
    steps: List[float],
) -> List[Dict[str, Any]]:
    rows = []
    for mx in steps:
        label = f"monoff_res{res}_mx{_cell_tag(mx)}_P1P2"
        print(f"\n=== {label} ===", flush=True)
        row = _run_pair(
            label=label,
            ports="0,1",
            res=res,
            run_time=run_time,
            ranks=ranks,
            monitor_offset=(mx, 0.0),
        )
        rows.append(row)
    return rows


def make_diagnostic_figure(index: Dict[str, Any], out_path: Path) -> Path:
    """2x2 figure: offset heatmap, resolution curves, wall compare, monitor offset."""
    plt.rcParams.update({"figure.dpi": 140, "savefig.dpi": 200, "font.size": 10})
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    # --- offset heatmap (P1P2) ---
    ax = axes[0, 0]
    offset_rows = [
        r
        for r in index.get("runs", [])
        if r.get("label", "").startswith("offset_res")
        and r.get("P12_dB") is not None
        and r.get("grid_offset_cells") is not None
        and r.get("monitor_offset_cells", [0, 0]) == [0.0, 0.0]
    ]
    if offset_rows:
        # pick highest-res offset sweep present
        res_vals = sorted({r.get("res", 32) for r in offset_rows})
        res_use = res_vals[-1]
        subset = [r for r in offset_rows if r.get("res", 32) == res_use]
        ox = sorted({r["grid_offset_cells"][0] for r in subset})
        oy = sorted({r["grid_offset_cells"][1] for r in subset})
        grid = np.full((len(oy), len(ox)), np.nan)
        for r in subset:
            i = oy.index(r["grid_offset_cells"][1])
            j = ox.index(r["grid_offset_cells"][0])
            grid[i, j] = r["P12_dB"]
        im = ax.imshow(
            grid,
            origin="lower",
            aspect="auto",
            extent=[ox[0], ox[-1], oy[0], oy[-1]],
            cmap="magma",
        )
        fig.colorbar(im, ax=ax, label=r"$|P_{12}-P_{21}|$ (dB)")
        ax.set_xlabel("grid offset ox (cells)")
        ax.set_ylabel("grid offset oy (cells)")
        ax.set_title(f"P1↔P2 vs horn grid placement (res={res_use})")
        baseline = [
            r["P12_dB"]
            for r in subset
            if r["grid_offset_cells"] == [0.0, 0.0]
        ]
        if baseline:
            ax.text(
                0.02,
                0.98,
                f"zero-offset: {baseline[0]:.2f} dB",
                transform=ax.transAxes,
                va="top",
                fontsize=9,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
            )
    else:
        ax.text(0.5, 0.5, "no offset sweep data", ha="center", va="center")
        ax.set_axis_off()

    # --- resolution dependence ---
    ax = axes[0, 1]
    res_rows = [r for r in index.get("runs", []) if "zero_offset" in r.get("label", "")]
    p12 = {r["res"]: r["P12_dB"] for r in res_rows if "P1P2" in r["label"]}
    p23 = {r["res"]: r["P12_dB"] for r in res_rows if "P2P3" in r["label"]}
    if p12:
        xs = sorted(p12)
        ax.plot(xs, [p12[r] for r in xs], "o-", color="#E45756", label="P1↔P2")
    if p23:
        xs = sorted(p23)
        ax.plot(xs, [p23[r] for r in xs], "s-", color="#4C78A8", label="P2↔P3")
    ax.axhline(0.2, color="0.4", ls="--", lw=1, label="0.2 dB gate")
    ax.set_xlabel("resolution (px/a)")
    ax.set_ylabel(r"max $|P_{ij}-P_{ji}|$ (dB)")
    ax.set_title("Reciprocity vs resolution (zero offset)")
    ax.legend(loc="best")
    ax.set_ylim(bottom=-0.02)

    # --- wall model ---
    ax = axes[1, 0]
    wall_rows = [
        r
        for r in index.get("runs", [])
        if r.get("label", "").startswith("wall_") and r.get("P12_dB") is not None
    ]
    if wall_rows:
        labels = [r["wall_mode"] for r in wall_rows]
        vals = [r["P12_dB"] for r in wall_rows]
        ax.bar(labels, vals, color=["#54A24B", "#F58518"])
        ax.set_ylabel(r"P1↔P2 $|P_{12}-P_{21}|$ (dB)")
        ax.set_title("PEC vs high-ε metal walls (res32)")
    else:
        ax.text(0.5, 0.5, "no wall-model data", ha="center", va="center")
        ax.set_axis_off()

    # --- monitor offset ---
    ax = axes[1, 1]
    mon_rows = sorted(
        [r for r in index.get("runs", []) if r.get("label", "").startswith("monoff_")],
        key=lambda r: r["monitor_offset_cells"][0],
    )
    if mon_rows:
        xs = [r["monitor_offset_cells"][0] for r in mon_rows]
        ys = [r["P12_dB"] for r in mon_rows]
        ax.plot(xs, ys, "o-", color="#72B7B2")
        ax.set_xlabel("monitor offset along outward (cells)")
        ax.set_ylabel(r"P1↔P2 (dB)")
        ax.set_title("Monitor-only sub-cell shift (walls fixed)")
    else:
        ax.text(0.5, 0.5, "no monitor-offset data", ha="center", va="center")
        ax.set_axis_off()

    fig.suptitle(
        "Is the reciprocity error caused by rotated horns on the Cartesian grid?",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks", type=int, default=32)
    parser.add_argument("--run-time", type=float, default=40.0)
    parser.add_argument(
        "--phase",
        type=str,
        default="all",
        choices=["all", "offset", "resolution", "wall", "monitor", "plot", "candidates"],
    )
    parser.add_argument(
        "--offset-steps",
        type=str,
        default="-0.5,-0.25,0,0.25,0.5",
        help="Fractional-cell values for ox/oy sweeps (use = form if negative, e.g. --offset-steps=-0.5,0,0.5)",
    )
    parser.add_argument(
        "--offset-res",
        type=int,
        nargs="+",
        default=[32, 48],
    )
    parser.add_argument(
        "--resolution-list",
        type=int,
        nargs="+",
        default=[32, 48, 64, 96],
    )
    parser.add_argument("--horn-walls", type=str, default="prism")
    parser.add_argument("--coord-rotation-deg", type=float, default=0.0)
    args = parser.parse_args()

    steps = [float(x) for x in args.offset_steps.split(",") if x.strip()]
    index_path = OUT / "index.json"
    index: Dict[str, Any] = {"runs": []}
    if index_path.is_file():
        with open(index_path, encoding="utf-8") as f:
            index = json.load(f)

    def _append(rows: List[Dict[str, Any]]) -> None:
        existing = {r.get("label") for r in index["runs"]}
        for row in rows:
            if row.get("label") not in existing:
                index["runs"].append(row)
            else:
                for i, old in enumerate(index["runs"]):
                    if old.get("label") == row.get("label"):
                        index["runs"][i] = row
                        break
        index["updated_utc"] = datetime.now(timezone.utc).isoformat()
        with open(index_path, "w", encoding="utf-8") as f:
            json.dump(index, f, indent=2)

    if args.phase in ("all", "offset"):
        for res in args.offset_res:
            _append(offset_sweep(res, args.run_time, args.ranks, steps, ports="0,1"))

    if args.phase in ("all", "resolution"):
        _append(
            resolution_sweep(
                args.resolution_list,
                args.run_time,
                args.ranks,
                horn_walls=args.horn_walls,
                coord_rotation_deg=args.coord_rotation_deg,
                label_prefix=f"res_{args.horn_walls}_r{args.coord_rotation_deg:g}_",
            )
        )

    if args.phase in ("all", "wall"):
        _append(wall_mode_compare(32, args.run_time, args.ranks))

    if args.phase in ("all", "monitor"):
        _append(monitor_offset_sweep(32, args.run_time, args.ranks, steps))

    if args.phase in ("all", "candidates"):
        _append(
            [
                r
                for block in candidate_offset_screen(
                    32, args.run_time, args.ranks, steps
                )
                for r in block["runs"]
            ]
        )

    if args.phase in ("all", "plot"):
        fig_path = OUT / "horn_grid_diagnostic.png"
        make_diagnostic_figure(index, fig_path)
        print(f"Wrote {fig_path}")

    subprocess.run([PY, str(VAL / "validation_registry.py")], cwd=str(ROOT), check=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
