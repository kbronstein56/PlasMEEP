#!/usr/bin/env python3
"""Forward-only dense local audit of a feasibility-band optimum.

Uses 3.125 MHz (default) within the optimized band, reading B / symmetry /
n_var from the opt JSON so C3 and non-production-B probes audit correctly.

Stage-6 NOT STARTED.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
FEAS = OUT / "feasibility"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "jesse_method_campaign"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from audit_bb_stage5_dense import freq_grid, row_from_bundle  # noqa: E402
from multimode_opt_core import B_PLUS, orbits_and_masks  # noqa: E402


def _init_worker(omp_threads: int):
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = str(int(omp_threads))


def _worker_eval(payload: dict) -> dict:
    _init_worker(int(payload["omp_threads"]))
    from eval_forward_only import eval_q_forward

    q = np.asarray(payload["q"], float)
    orbits_payload = payload.get("orbits")
    if orbits_payload is None:
        orbits, _, _ = orbits_and_masks()
    else:
        orbits = [{"rods": list(o["rods"]), "size": int(o.get("size", len(o["rods"])))} for o in orbits_payload]
    f_ghz = float(payload["f_GHz"])
    bun, _ = eval_q_forward(
        payload["grade"],
        q,
        orbits,
        f_hz=f_ghz * 1e9,
        b_tesla=float(payload["b_tesla"]),
    )
    return row_from_bundle(f_ghz, bun)


def _orbits_for_opt(d: dict):
    sym = d.get("symmetry", "C6")
    n_var = int(d.get("n_var", 16))
    if sym == "C6" and n_var == 16:
        return None  # workers rebuild via orbits_and_masks
    from full91_geometry import bulb_centers_device
    import sixport_common as sc
    from threefold_orbits import build_orbits_c3

    centers = bulb_centers_device(sc)
    orbits, _ = build_orbits_c3(centers)
    assert len(orbits) == n_var
    return [{"rods": list(o["rods"] if isinstance(o, dict) else o), "size": len(o["rods"] if isinstance(o, dict) else o)} for o in orbits]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--opt-json", required=True)
    ap.add_argument("--grade", default="C")
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--omp-threads", type=int, default=1)
    ap.add_argument("--df-MHz", type=float, default=3.125)
    ap.add_argument("--b-tesla", type=float, default=None)
    ap.add_argument("--refine-MHz", type=float, default=1.5625, help="local refine around A2<=-5")
    args = ap.parse_args()

    for k, v in (
        ("OMP_NUM_THREADS", str(args.omp_threads)),
        ("OPENBLAS_NUM_THREADS", str(args.omp_threads)),
        ("MKL_NUM_THREADS", str(args.omp_threads)),
        ("PYTHONUNBUFFERED", "1"),
    ):
        os.environ.setdefault(k, v)

    opt_path = Path(args.opt_json)
    d = json.loads(opt_path.read_text())
    q = np.asarray(d["q_best"], float)
    f0, f1 = float(d["f0_GHz"]), float(d["f1_GHz"])
    b_tesla = float(args.b_tesla if args.b_tesla is not None else d.get("b_tesla", B_PLUS))
    orbits_payload = _orbits_for_opt(d)
    df = args.df_MHz * 1e-3
    freqs = freq_grid(f0, f1, df)
    tag = f"audit_{opt_path.stem}_{args.df_MHz:g}MHz"
    resume = FEAS / f"{tag}.json"
    print(
        f"FEASIBILITY AUDIT {opt_path.name} band=[{f0},{f1}] df={args.df_MHz} MHz "
        f"n={len(freqs)} B={b_tesla} T sym={d.get('symmetry','C6')} n_var={len(q)} "
        f"— Stage-6 NOT STARTED",
        flush=True,
    )

    done = {}
    if resume.exists():
        prev = json.loads(resume.read_text())
        for r in prev.get("rows", []):
            done[round(float(r["f_GHz"]), 6)] = r

    missing = [f for f in freqs if round(float(f), 6) not in done]
    rows = [done[round(float(f), 6)] for f in freqs if round(float(f), 6) in done]
    print(f"{tag}: {len(freqs)} requested, {len(rows)} seeded, {len(missing)} to solve", flush=True)

    def _run(missing_freqs, existing_rows):
        if not missing_freqs:
            return existing_rows
        payloads = [
            {
                "grade": args.grade,
                "q": q.tolist(),
                "f_GHz": float(f),
                "omp_threads": int(args.omp_threads),
                "b_tesla": b_tesla,
                "orbits": orbits_payload,
            }
            for f in missing_freqs
        ]
        import multiprocessing as mp

        ctx = mp.get_context("spawn")
        t0 = time.perf_counter()
        with ProcessPoolExecutor(
            max_workers=max(1, int(args.n_workers)),
            mp_context=ctx,
            initializer=_init_worker,
            initargs=(int(args.omp_threads),),
        ) as ex:
            futs = {ex.submit(_worker_eval, p): p["f_GHz"] for p in payloads}
            n_done = 0
            for fut in as_completed(futs):
                row = fut.result()
                existing_rows.append(row)
                n_done += 1
                print(
                    f"{tag} [{n_done}/{len(missing_freqs)}] f={row['f_GHz']:.5f} "
                    f"A2={row['J_A2']:.4f} IL={row['IL_dB']:.2f} "
                    f"iso={row['isolation_dB']:.2f}",
                    flush=True,
                )
                if n_done % 12 == 0 or n_done == len(missing_freqs):
                    _write(existing_rows)
        print(f"{tag}: solved {len(missing_freqs)} in {time.perf_counter()-t0:.1f}s", flush=True)
        return existing_rows

    def _write(rows_in):
        rows_sorted = sorted(rows_in, key=lambda r: r["f_GHz"])
        # dedupe by freq
        byf = {round(float(r["f_GHz"]), 6): r for r in rows_sorted}
        rows_sorted = [byf[k] for k in sorted(byf)]
        resume.write_text(json.dumps({"tag": tag, "rows": rows_sorted}, indent=2) + "\n")
        return rows_sorted

    rows = _run(missing, rows)
    rows = _write(rows)

    # Local refine around deep failures
    deep = [r for r in rows if r["J_A2"] <= -5.0]
    if deep and args.refine_MHz > 0:
        df_r = args.refine_MHz * 1e-3
        refine_freqs = set()
        for r in deep:
            c = float(r["f_GHz"])
            for k in range(-4, 5):
                f = round(c + k * df_r, 6)
                if f0 - 1e-9 <= f <= f1 + 1e-9:
                    refine_freqs.add(f)
        have = {round(float(r["f_GHz"]), 6) for r in rows}
        miss_r = sorted(f for f in refine_freqs if f not in have)
        if miss_r:
            print(f"REFINE {args.refine_MHz} MHz around {len(deep)} A2<=-5 sites: {len(miss_r)} freqs", flush=True)
            rows = _run(miss_r, rows)
            rows = _write(rows)

    worst = min(rows, key=lambda r: r["J_A2"])
    n_le5 = sum(1 for r in rows if r["J_A2"] <= -5)
    n_le7 = sum(1 for r in rows if r["J_A2"] <= -7)
    summary = {
        "opt": str(opt_path),
        "J_BB_best": d.get("J_BB_best"),
        "band": [f0, f1],
        "df_MHz": args.df_MHz,
        "refine_MHz": args.refine_MHz,
        "b_tesla": b_tesla,
        "symmetry": d.get("symmetry", "C6"),
        "n_var": len(q),
        "n_freq": len(rows),
        "worst_A2": {
            "f_GHz": worst["f_GHz"],
            "J_A2": worst["J_A2"],
            "IL_dB": worst["IL_dB"],
            "isolation_dB": worst["isolation_dB"],
            "T_des": worst["T_des"],
        },
        "worst_IL_dB": max(r["IL_dB"] for r in rows),
        "worst_isolation_dB": min(r["isolation_dB"] for r in rows),
        "n_A2_le_m5": n_le5,
        "n_A2_le_m7": n_le7,
        "stage6_status": "NOT STARTED / AWAITING AUTHORIZATION",
    }
    (FEAS / f"{tag}_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        f"TRUE WORST A2={worst['J_A2']:.4f} @{worst['f_GHz']:.5f} "
        f"IL={worst['IL_dB']:.2f} iso={worst['isolation_dB']:.2f} "
        f"n_A2<=-5={n_le5} n_A2<=-7={n_le7}",
        flush=True,
    )
    print("STAGE-6: NOT STARTED / AWAITING AUTHORIZATION", flush=True)


if __name__ == "__main__":
    main()
