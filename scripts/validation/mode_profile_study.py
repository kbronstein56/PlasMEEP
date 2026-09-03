#!/usr/bin/env python3
"""
Extract numerical Hz mode profiles at monitor planes and compare to te1_hz_line.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import meep as mp
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VAL = os.path.dirname(os.path.abspath(__file__))
for p in (VAL, os.path.join(ROOT, "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from plasmeep.ports.mode_profile import (  # noqa: E402
    ModeLineProfile,
    overlap_metrics,
    te1_cos_amplitude,
)
from plasmeep.ports.numerical_mode import NumericalPortMode  # noqa: E402
from port_formulations import make_te1_hz_sources  # noqa: E402
from sixport_common import (  # noqa: E402
    a,
    build_circulator_device,
    clear_width,
    default_uniform_rho,
    effective_port_dir,
    fs_a,
    horn_for_port,
    monitor_center_for_port,
    physical_resolution_report,
    set_geometry_context,
    source_df,
)
from physical_units import add_resolution_arguments, resolve_simulation_resolution  # noqa: E402

OUT = Path(ROOT) / "outputs" / "validation" / "mode_profiles"


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, complex):
        return {"real": obj.real, "imag": obj.imag}
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def extract_hz_line_profile(
    p_device,
    port_index: int,
    res: int,
    *,
    frequency: float,
    fwidth: float,
    run_time: float,
    n_points: int = 64,
) -> ModeLineProfile:
    """Excite port with te1_hz_line; DFT Hz along monitor tangent."""
    p_device.sources = make_te1_hz_sources(port_index, frequency=frequency, fwidth=fwidth)
    sim = p_device.Get_Sim()

    center = np.asarray(monitor_center_for_port(port_index, res), dtype=float)
    u = np.asarray(effective_port_dir(port_index), dtype=float)
    u = u / np.linalg.norm(u)
    tangent = np.array([-u[1], u[0]])
    span = 0.96 * clear_width
    offsets = np.linspace(-span / 2, span / 2, n_points)
    pts = center[None, :] + offsets[:, None] * tangent[None, :]

    dft_objs = []
    for xy in pts:
        dft_objs.append(
            sim.add_dft_fields(
                [mp.Hz],
                frequency,
                0,
                1,
                center=mp.Vector3(xy[0], xy[1], 0),
                size=mp.Vector3(0, 0, 0),
            )
        )

    sim.run(until_after_sources=run_time)

    hz = np.array(
        [
            complex(np.squeeze(sim.get_dft_array(dft, mp.Hz, 0)))
            for dft in dft_objs
        ],
        dtype=complex,
    )
    return ModeLineProfile(
        port_index=port_index,
        offsets_a=offsets,
        hz_complex=hz,
        span_a=span,
        monitor_xy=(float(center[0]), float(center[1])),
        tangent=(float(tangent[0]), float(tangent[1])),
    )


def plot_profiles(
    profiles: List[ModeLineProfile],
    overlaps: List[Dict[str, Any]],
    out_png: Path,
) -> None:
    labels = {0: "P1 horizontal", 1: "P2 +60°", 5: "P6 −60°"}
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    ax = axes[0]
    for prof in profiles:
        port = prof.port_index
        norm = prof.normalized()
        analytic = te1_cos_amplitude(prof.offsets_a, prof.span_a)
        analytic = analytic / max(np.max(np.abs(analytic)), 1e-30)
        ax.plot(
            prof.offsets_a * a * 1e3,
            np.real(norm),
            label=f"{labels.get(port, f'P{port+1}')} Re(Hz) num",
        )
        ax.plot(
            prof.offsets_a * a * 1e3,
            analytic,
            "--",
            alpha=0.7,
            label=f"{labels.get(port, f'P{port+1}')} te1 cos",
        )
    ax.set_ylabel("normalized amplitude")
    ax.set_title("Numerical Hz profiles vs analytic TE1 cos envelope")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    names = [labels.get(p["port_index"], f"P{p['port_index']+1}") for p in overlaps]
    ov = [p["analytic_te1"]["power_overlap"] for p in overlaps]
    bars = ax.bar(names, ov, color=["#4C78A8", "#F58518", "#54A24B"])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("power overlap |⟨mode|te1⟩|²")
    ax.set_title("Analytic TE1 overlap with numerical guided mode")
    for b, val in zip(bars, ov):
        ax.text(
            b.get_x() + b.get_width() / 2,
            val + 0.02,
            f"{val:.3f}",
            ha="center",
            fontsize=9,
        )
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_xlabel("offset along horn tangent (mm)")

    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_resolution_arguments(parser)
    parser.add_argument("--run-time", type=float, default=40.0)
    parser.add_argument("--ports", type=str, default="0,1,2,3,4,5")
    parser.add_argument("--horn-walls", type=str, default="prism")
    args = parser.parse_args()

    sim_res, points_per_cm, res_report = resolve_simulation_resolution(
        args, a_m=a, default_res=32
    )

    ports = [int(x) for x in args.ports.split(",")]
    set_geometry_context(res=sim_res, horn_walls=args.horn_walls)
    os.makedirs(OUT, exist_ok=True)

    rho = default_uniform_rho()
    B = np.zeros(3)
    _pmm, p_device, _ = build_circulator_device(
        rho, B, res=sim_res, device_mode="horns_only"
    )

    profiles: List[ModeLineProfile] = []
    overlaps: List[Dict[str, Any]] = []
    modes_dir = OUT / f"res{sim_res}"
    modes_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    for port in ports:
        print(f"\n=== mode profile port P{port+1} ===", flush=True)
        prof = extract_hz_line_profile(
            p_device,
            port,
            sim_res,
            frequency=fs_a,
            fwidth=source_df,
            run_time=args.run_time,
        )
        profiles.append(prof)
        analytic = te1_cos_amplitude(prof.offsets_a, prof.span_a)
        metrics = overlap_metrics(prof.hz_complex, analytic)
        overlaps.append(
            {
                "port_index": port,
                "analytic_te1": metrics,
                "peak_hz_abs": float(np.max(np.abs(prof.hz_complex))),
            }
        )
        print(
            f"  power_overlap={metrics['power_overlap']:.4f}  "
            f"|overlap|={metrics['overlap_magnitude']:.4f}  "
            f"phase={metrics['projection_phase_deg']:.2f}°",
            flush=True,
        )

        mode = NumericalPortMode.from_profile(
            prof,
            label=f"P{port+1}_ppc{points_per_cm:g}",
            frequency_a=fs_a,
            res=sim_res,
            outward_dir=effective_port_dir(port),
        )
        mode.save(modes_dir / f"numerical_mode_P{port+1}.json")

    # Cross-port overlaps (numerical modes)
    prof_by_port = {p.port_index: p for p in profiles}
    cross: List[Dict[str, Any]] = []
    for i, pa in enumerate(ports):
        for pb in ports[i + 1 :]:
            fa = prof_by_port[pa].normalized()
            fb = prof_by_port[pb].normalized()
            ov = complex(np.vdot(fa, fb))
            cross.append(
                {
                    "port_a": pa,
                    "port_b": pb,
                    "overlap_magnitude": float(np.abs(ov)),
                    "power_overlap": float(np.abs(ov) ** 2),
                }
            )

    payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "settings": {
            "points_per_cm": points_per_cm,
            "res": sim_res,
            "run_time": args.run_time,
            "ports": ports,
            "horn_walls": args.horn_walls,
            "excitation": "te1_hz_line (for profile extraction only)",
        },
        "physical_resolution": res_report,
        "profiles": [p.as_dict() for p in profiles],
        "analytic_overlaps": overlaps,
        "cross_port_overlaps": cross,
        "timings_s": {"total": time.perf_counter() - t0},
    }

    json_path = OUT / f"mode_profiles_ppc{points_per_cm:g}_res{sim_res}.json"
    with json_path.open("w", encoding="utf-8") as f:
        json.dump(_json_safe(payload), f, indent=2)
        f.write("\n")

    plot_profiles(profiles, overlaps, OUT / f"mode_profiles_ppc{points_per_cm:g}_res{sim_res}.png")
    print(f"\nWrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
