"""Full-band dense Stage-4 spectral audit (forward-only, no optimization).

Default: 25 MHz over 4.0-6.0 GHz with process-parallel frequencies
(n_workers=12, omp=1). After the full scan, refine at 12.5 MHz around
suspicious minima (always including 4.775 and 5.75-5.90), then optionally
tighten sharp notches to 6.25 MHz.

Compares quantitatively against Stage-3 dense audit. Does not launch Stage-5.
"""
from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
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

from multimode_opt_core import orbits_and_masks  # noqa: E402
from objective_A2 import objective_A2_from_S  # noqa: E402
from spectral_recon import ensure_summary_fields, summarize  # noqa: E402


def _init_worker(omp_threads: int):
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = str(int(omp_threads))


def _worker_eval(payload: dict) -> dict:
    """Spawn-safe single-frequency forward eval."""
    _init_worker(int(payload["omp_threads"]))
    # Late imports inside spawn worker
    from eval_forward_only import eval_q_forward
    from multimode_opt_core import B_PLUS, orbits_and_masks
    from objective_A2 import objective_A2_from_S

    q = np.asarray(payload["q"], float)
    orbits, _, _ = orbits_and_masks()
    f_ghz = float(payload["f_GHz"])
    bun, _ = eval_q_forward(
        payload["grade"], q, orbits, f_hz=f_ghz * 1e9, b_tesla=B_PLUS
    )
    return row_from_bundle(f_ghz, bun)


def row_from_bundle(f_ghz: float, bun: dict) -> dict:
    m = bun["metrics"]
    J_A2 = float(objective_A2_from_S(bun["S"]))
    per = m.get("per_source") or []
    # Worst competing = reverse path at the source with minimum isolation
    if per:
        worst_iso_src = min(per, key=lambda r: r["isolation_dB"])
        j = int(worst_iso_src["source"])
        competing_port = (j - 1) % 6  # reverse for B_PLUS
        desired_port = (j + 1) % 6
        worst_competing = {
            "source_port": j + 1,
            "desired_port": desired_port + 1,
            "competing_port": competing_port + 1,
            "T_des": float(worst_iso_src["T_des"]),
            "T_rev": float(worst_iso_src["T_rev"]),
            "isolation_dB": float(worst_iso_src["isolation_dB"]),
            "IL_dB": float(worst_iso_src["IL_dB"]),
        }
    else:
        worst_competing = None
    return {
        "f_GHz": float(f_ghz),
        "J_A2": J_A2,
        "J_A": float(bun["J_A"]),
        "J_B": float(bun["J_B"]),
        "T_des": float(m["T_des_avg"]),
        "T_rev": float(m["T_rev_avg"]),
        "IL_dB": float(m["IL_avg_dB"]),
        "IL_worst_dB": float(m.get("IL_worst_dB", m["IL_avg_dB"])),
        "isolation_dB": float(m["isolation_avg_dB"]),
        "isolation_min_dB": float(m.get("isolation_min_dB", m["isolation_avg_dB"])),
        "RL_dB": float(m["RL_avg_dB"]),
        "R": float(m["R_avg"]),
        "TE1": float(m["TE1_avg"]),
        "worst_competing": worst_competing,
        "timings": bun.get("timings"),
    }


def freq_grid(f0: float, f1: float, df: float) -> list[float]:
    n = int(round((f1 - f0) / df)) + 1
    return [round(f0 + i * df, 6) for i in range(n)]


def load_seed_rows(paths: list[Path]) -> dict[float, dict]:
    done: dict[float, dict] = {}
    for p in paths:
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        for r in d.get("rows", []):
            if "J_A2" not in r:
                continue
            key = round(float(r["f_GHz"]), 6)
            # Prefer richer rows if already present
            if key not in done or (
                done[key].get("worst_competing") is None
                and r.get("worst_competing") is not None
            ):
                done[key] = r
    return done


def write_checkpoint(tag: str, grade: str, rows: list[dict], path: Path):
    rows = sorted(rows, key=lambda r: r["f_GHz"])
    summary = summarize(tag, grade, rows)
    summary = ensure_summary_fields(summary, grade=grade)
    # Add A2-specific worst (summarize uses J_A historically)
    worst_a2 = min(rows, key=lambda r: r["J_A2"])
    worst_il = max(rows, key=lambda r: r["IL_dB"])
    worst_iso = min(rows, key=lambda r: r["isolation_dB"])
    out = {
        "tag": tag,
        "grade": grade,
        "rows": rows,
        **summary,
        "worst_A2": {
            "f_GHz": worst_a2["f_GHz"],
            "J_A2": worst_a2["J_A2"],
            "IL_dB": worst_a2["IL_dB"],
            "isolation_dB": worst_a2["isolation_dB"],
            "T_des": worst_a2["T_des"],
        },
        "worst_IL": {
            "f_GHz": worst_il["f_GHz"],
            "IL_dB": worst_il["IL_dB"],
            "J_A2": worst_il["J_A2"],
            "isolation_dB": worst_il["isolation_dB"],
        },
        "worst_isolation": {
            "f_GHz": worst_iso["f_GHz"],
            "isolation_dB": worst_iso["isolation_dB"],
            "J_A2": worst_iso["J_A2"],
            "IL_dB": worst_iso["IL_dB"],
        },
    }
    path.write_text(json.dumps(out, indent=2) + "\n")
    return out


def write_csv(rows: list[dict], path: Path):
    fields = [
        "f_GHz",
        "J_A2",
        "J_A",
        "J_B",
        "T_des",
        "T_rev",
        "IL_dB",
        "IL_worst_dB",
        "isolation_dB",
        "isolation_min_dB",
        "RL_dB",
        "R",
        "TE1",
        "competing_source_port",
        "competing_port",
        "competing_T_rev",
        "competing_iso_dB",
    ]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in sorted(rows, key=lambda x: x["f_GHz"]):
            wc = r.get("worst_competing") or {}
            w.writerow(
                {
                    "f_GHz": r["f_GHz"],
                    "J_A2": r["J_A2"],
                    "J_A": r.get("J_A"),
                    "J_B": r.get("J_B"),
                    "T_des": r["T_des"],
                    "T_rev": r.get("T_rev"),
                    "IL_dB": r["IL_dB"],
                    "IL_worst_dB": r.get("IL_worst_dB"),
                    "isolation_dB": r["isolation_dB"],
                    "isolation_min_dB": r.get("isolation_min_dB"),
                    "RL_dB": r.get("RL_dB"),
                    "R": r.get("R"),
                    "TE1": r.get("TE1"),
                    "competing_source_port": wc.get("source_port"),
                    "competing_port": wc.get("competing_port"),
                    "competing_T_rev": wc.get("T_rev"),
                    "competing_iso_dB": wc.get("isolation_dB"),
                }
            )


def parallel_scan(
    grade: str,
    q: np.ndarray,
    freqs: list[float],
    *,
    tag: str,
    resume_path: Path,
    seed_paths: list[Path],
    n_workers: int,
    omp_threads: int,
):
    done = load_seed_rows(seed_paths + ([resume_path] if resume_path.exists() else []))
    # Keep only requested freqs in working set
    want = {round(float(f), 6) for f in freqs}
    rows = [done[k] for k in sorted(want) if k in done]
    missing = [f for f in freqs if round(float(f), 6) not in done]
    print(
        f"{tag}: {len(freqs)} requested, {len(rows)} seeded, {len(missing)} to solve "
        f"(n_workers={n_workers}, omp={omp_threads})",
        flush=True,
    )
    if not missing:
        out = write_checkpoint(tag, grade, rows, resume_path)
        write_csv(rows, resume_path.with_suffix(".csv"))
        return out

    thread_keys = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
    saved = {k: os.environ.get(k) for k in thread_keys}
    for k in thread_keys:
        os.environ[k] = str(int(omp_threads))
    t0 = time.perf_counter()
    payloads = [
        {
            "grade": grade,
            "q": np.asarray(q, float).tolist(),
            "f_GHz": float(f),
            "omp_threads": int(omp_threads),
        }
        for f in missing
    ]
    ctx = mp.get_context("spawn")
    try:
        with ProcessPoolExecutor(
            max_workers=max(1, int(n_workers)),
            mp_context=ctx,
            initializer=_init_worker,
            initargs=(int(omp_threads),),
        ) as ex:
            futs = {ex.submit(_worker_eval, p): p["f_GHz"] for p in payloads}
            n_done = 0
            for fut in as_completed(futs):
                row = fut.result()
                rows.append(row)
                n_done += 1
                print(
                    f"{tag} [{n_done}/{len(missing)}] f={row['f_GHz']:.3f} "
                    f"A2={row['J_A2']:.4f} IL={row['IL_dB']:.2f} "
                    f"iso={row['isolation_dB']:.2f} Tdes={row['T_des']:.4e}",
                    flush=True,
                )
                write_checkpoint(tag, grade, rows, resume_path)
                write_csv(rows, resume_path.with_suffix(".csv"))
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    wall = time.perf_counter() - t0
    print(f"{tag}: solved {len(missing)} freqs in {wall:.1f}s", flush=True)
    out = write_checkpoint(tag, grade, rows, resume_path)
    write_csv(rows, resume_path.with_suffix(".csv"))
    return out


def local_minima(rows: list[dict], key: str = "J_A2", pad: int = 1) -> list[dict]:
    """Interior points worse than both neighbors (local A2 minima)."""
    rs = sorted(rows, key=lambda r: r["f_GHz"])
    out = []
    for i in range(pad, len(rs) - pad):
        v = rs[i][key]
        if all(v < rs[i + d][key] for d in range(-pad, pad + 1) if d != 0):
            out.append(rs[i])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--omp-threads-per-worker", type=int, default=1)
    ap.add_argument("--f-min", type=float, default=4.0)
    ap.add_argument("--f-max", type=float, default=6.0)
    ap.add_argument("--df-mhz", type=float, default=25.0)
    ap.add_argument(
        "--stage",
        choices=("full25", "refine", "fine625", "all"),
        default="all",
        help="full25 / refine(12.5) / fine625 / all",
    )
    ap.add_argument(
        "--refine-centers",
        default="4.775,5.8125,5.75,5.825",
        help="Comma GHz centers for 12.5 MHz refine (plus auto local minima)",
    )
    ap.add_argument("--refine-half-width-ghz", type=float, default=0.075)
    ap.add_argument(
        "--fine-centers",
        default="5.8125,4.775",
        help="Comma GHz centers for optional 6.25 MHz fine refine of sharp notches",
    )
    ap.add_argument("--fine-half-width-ghz", type=float, default=0.025)
    args = ap.parse_args()

    for k, v in (
        ("OMP_NUM_THREADS", "1"),
        ("OPENBLAS_NUM_THREADS", "1"),
        ("MKL_NUM_THREADS", "1"),
    ):
        os.environ.setdefault(k, v)

    design_path = OUT / "opt_broadband_stage4_A2_C.json"
    d4 = json.loads(design_path.read_text())
    q4 = np.asarray(d4["q_best"], float)
    print(
        f"Stage-4 dense audit design={design_path.name} "
        f"J_BB={d4.get('J_BB_best')} nit={d4.get('accepted_nit_total', d4.get('accepted_nit'))}",
        flush=True,
    )

    seed = [
        OUT / "bb_stage4_dense_25MHz.json",
        OUT / "bb_stage4_dense_12p5MHz_local.json",
    ]

    report = {
        "design": str(design_path),
        "J_BB_best": d4.get("J_BB_best"),
        "accepted_nit": d4.get("accepted_nit_total", d4.get("accepted_nit")),
        "n_workers": args.n_workers,
        "omp_threads_per_worker": args.omp_threads_per_worker,
        "note": "Fresh Stage-4 dense audit; do not launch Stage-5 automatically.",
    }

    a25 = None
    if args.stage in ("full25", "all"):
        df = args.df_mhz / 1000.0
        freqs = freq_grid(args.f_min, args.f_max, df)
        print(
            f"FULL 25 MHz grid: {freqs[0]}–{freqs[-1]} GHz, n={len(freqs)}, df={df}",
            flush=True,
        )
        a25 = parallel_scan(
            args.grade,
            q4,
            freqs,
            tag="bb_stage4_dense_25MHz",
            resume_path=OUT / "bb_stage4_dense_25MHz.json",
            seed_paths=seed,
            n_workers=args.n_workers,
            omp_threads=args.omp_threads_per_worker,
        )
        report["dense_25MHz"] = {
            "n_freq": a25["n_freq"],
            "spacing_GHz": df,
            "worst_A2": a25["worst_A2"],
            "worst_IL": a25["worst_IL"],
            "worst_isolation": a25["worst_isolation"],
            "local_A2_minima": [
                {
                    "f_GHz": r["f_GHz"],
                    "J_A2": r["J_A2"],
                    "IL_dB": r["IL_dB"],
                    "isolation_dB": r["isolation_dB"],
                }
                for r in local_minima(a25["rows"])
            ],
        }

    if args.stage in ("refine", "all"):
        if a25 is None:
            a25 = json.loads((OUT / "bb_stage4_dense_25MHz.json").read_text())
        centers = []
        for tok in args.refine_centers.split(","):
            tok = tok.strip()
            if tok:
                centers.append(float(tok))
        # Auto-add worst A2 and local minima
        centers.append(float(a25["worst_A2"]["f_GHz"]))
        for r in local_minima(a25["rows"]):
            centers.append(float(r["f_GHz"]))
        # Unique within 12.5 MHz
        centers_u = []
        for c in sorted(centers):
            if not centers_u or abs(c - centers_u[-1]) > 0.012:
                centers_u.append(c)
        refine_freqs = set()
        hw = float(args.refine_half_width_ghz)
        for c in centers_u:
            for f in freq_grid(c - hw, c + hw, 0.0125):
                if args.f_min - 1e-9 <= f <= args.f_max + 1e-9:
                    refine_freqs.add(round(f, 6))
        # Always include existing dense points near centers for continuity
        refine_list = sorted(refine_freqs)
        print(
            f"REFINE 12.5 MHz: centers={centers_u} n={len(refine_list)}",
            flush=True,
        )
        a12 = parallel_scan(
            args.grade,
            q4,
            refine_list,
            tag="bb_stage4_dense_12p5MHz_local",
            resume_path=OUT / "bb_stage4_dense_12p5MHz_local.json",
            seed_paths=seed + [OUT / "bb_stage4_dense_25MHz.json"],
            n_workers=args.n_workers,
            omp_threads=args.omp_threads_per_worker,
        )
        report["dense_12p5MHz_local"] = {
            "centers_GHz": centers_u,
            "half_width_GHz": hw,
            "n_freq": a12["n_freq"],
            "worst_A2": a12["worst_A2"],
            "worst_IL": a12["worst_IL"],
            "worst_isolation": a12["worst_isolation"],
        }

    # Optional 6.25 MHz fine refine for sharp notches
    a625 = None
    if args.stage in ("fine625", "all"):
        centers_f = []
        for tok in args.fine_centers.split(","):
            tok = tok.strip()
            if tok:
                centers_f.append(float(tok))
        # Auto-add worst from 12.5 / 25 if very sharp
        src_for_worst = None
        if (OUT / "bb_stage4_dense_12p5MHz_local.json").exists():
            src_for_worst = json.loads((OUT / "bb_stage4_dense_12p5MHz_local.json").read_text())
        elif a25 is not None:
            src_for_worst = a25
        if src_for_worst is not None:
            centers_f.append(float(src_for_worst["worst_A2"]["f_GHz"]))
            for r in local_minima(src_for_worst["rows"]):
                if r["J_A2"] < -5.0:
                    centers_f.append(float(r["f_GHz"]))
        centers_fu = []
        for c in sorted(centers_f):
            if not centers_fu or abs(c - centers_fu[-1]) > 0.006:
                centers_fu.append(c)
        fine_freqs = set()
        hw6 = float(args.fine_half_width_ghz)
        for c in centers_fu:
            for f in freq_grid(c - hw6, c + hw6, 0.00625):
                if args.f_min - 1e-9 <= f <= args.f_max + 1e-9:
                    fine_freqs.add(round(f, 6))
        fine_list = sorted(fine_freqs)
        print(
            f"FINE 6.25 MHz: centers={centers_fu} n={len(fine_list)}",
            flush=True,
        )
        a625 = parallel_scan(
            args.grade,
            q4,
            fine_list,
            tag="bb_stage4_dense_6p25MHz_local",
            resume_path=OUT / "bb_stage4_dense_6p25MHz_local.json",
            seed_paths=seed
            + [
                OUT / "bb_stage4_dense_25MHz.json",
                OUT / "bb_stage4_dense_12p5MHz_local.json",
            ],
            n_workers=args.n_workers,
            omp_threads=args.omp_threads_per_worker,
        )
        report["dense_6p25MHz_local"] = {
            "centers_GHz": centers_fu,
            "half_width_GHz": hw6,
            "n_freq": a625["n_freq"],
            "worst_A2": a625["worst_A2"],
            "worst_IL": a625["worst_IL"],
            "worst_isolation": a625["worst_isolation"],
        }

    # Stage-3 vs Stage-4 quantitative comparison
    s3_25 = OUT / "bb_stage3_dense_25MHz.json"
    s3_12 = OUT / "bb_stage3_dense_12p5MHz_local.json"
    s3_loc58 = OUT / "bb_stage3_localize_5p8_6p25MHz.json"
    s3_loc477 = OUT / "bb_stage3_localize_4p775_6p25MHz.json"
    cmp = {"holes": {}, "notes": []}
    if a25 is None and (OUT / "bb_stage4_dense_25MHz.json").exists():
        a25 = json.loads((OUT / "bb_stage4_dense_25MHz.json").read_text())

    def at(rows, f):
        return min(rows, key=lambda r: abs(r["f_GHz"] - f))

    def pack(r):
        return {
            "f_GHz": r["f_GHz"],
            "J_A2": r.get("J_A2"),
            "IL_dB": r["IL_dB"],
            "isolation_dB": r["isolation_dB"],
            "T_des": r.get("T_des"),
        }

    if a25 is not None and s3_25.exists():
        s3rows = json.loads(s3_25.read_text())["rows"]
        s4rows = a25["rows"]
        for f in (4.1, 4.3, 4.775, 5.75, 5.775, 5.8125, 5.825, 5.85):
            src3 = s3rows
            # Prefer denser Stage-3 localizations where available
            if abs(f - 5.8125) < 1e-9 and s3_loc58.exists():
                src3 = json.loads(s3_loc58.read_text())["rows"]
            elif abs(f - 4.775) < 1e-9 and s3_loc477.exists():
                src3 = json.loads(s3_loc477.read_text())["rows"]
            elif f in (5.75, 5.825, 5.775) and s3_12.exists():
                src3 = json.loads(s3_12.read_text())["rows"]
            # Prefer Stage-4 fine/12.5 when available
            src4 = s4rows
            if a625 is not None and f in (5.8125, 4.775):
                src4 = a625["rows"]
            elif (OUT / "bb_stage4_dense_12p5MHz_local.json").exists() and f in (
                4.775,
                5.75,
                5.775,
                5.8125,
                5.825,
            ):
                src4 = json.loads((OUT / "bb_stage4_dense_12p5MHz_local.json").read_text())[
                    "rows"
                ]
            r3 = at(src3, f)
            r4 = at(src4, f)
            cmp["holes"][f"{f:g}"] = {"stage3": pack(r3), "stage4": pack(r4)}
        cmp["stage3_25_worst_A2"] = json.loads(s3_25.read_text())["worst_A2"]
        cmp["stage4_25_worst_A2"] = a25["worst_A2"]
        if a625 is not None:
            cmp["stage4_6p25_worst_A2"] = a625["worst_A2"]
        # Repair verdicts at localized Stage-3 corpses
        for tag, f, s3path in (
            ("5.8125", 5.8125, s3_loc58),
            ("4.775", 4.775, s3_loc477),
        ):
            if not s3path.exists():
                continue
            r3 = at(json.loads(s3path.read_text())["rows"], f)
            src4 = a625["rows"] if a625 is not None else s4rows
            r4 = at(src4, f)
            repaired = (r4["J_A2"] > -5.0) and (r4["IL_dB"] < 20.0)
            cmp[f"repair_{tag.replace('.', 'p')}"] = {
                "stage3": pack(r3),
                "stage4": pack(r4),
                "repaired": bool(repaired),
            }
    report["stage3_vs_stage4"] = cmp

    # True dense-audit worst across available Stage-4 scans
    true_worst = None
    for key in ("dense_6p25MHz_local", "dense_12p5MHz_local", "dense_25MHz"):
        block = report.get(key)
        if not block:
            continue
        w = block["worst_A2"]
        if true_worst is None or w["J_A2"] < true_worst["J_A2"]:
            true_worst = {**w, "source": key}
    report["true_dense_worst_A2"] = true_worst

    out_report = OUT / "BB_STAGE4_DENSE_AUDIT_SUMMARY.json"
    out_report.write_text(json.dumps(report, indent=2) + "\n")
    print("wrote", out_report, flush=True)
    print(json.dumps({
        "J_BB_best": report.get("J_BB_best"),
        "true_dense_worst_A2": true_worst,
        "dense_25MHz": report.get("dense_25MHz"),
        "stage3_vs_stage4_repairs": {
            k: v for k, v in cmp.items() if k.startswith("repair_")
        },
    }, indent=2), flush=True)

    # Markdown report
    md = OUT / "BB_STAGE4_DENSE_AUDIT.md"
    lines = [
        "# Stage-4 dense broadband audit",
        "",
        f"- Design: `{design_path.name}`",
        f"- Objective-grid J_BB_best: **{report.get('J_BB_best')}**",
        f"- True dense-audit worst A2: **{true_worst}**",
        "",
        "## Stage-3 vs Stage-4 (selected frequencies)",
        "",
    ]
    for k, v in sorted(cmp.get("holes", {}).items()):
        lines.append(f"### {k} GHz")
        lines.append(f"- Stage-3: `{v.get('stage3')}`")
        lines.append(f"- Stage-4: `{v.get('stage4')}`")
        lines.append("")
    for k, v in cmp.items():
        if k.startswith("repair_"):
            lines.append(f"## Repair check `{k}`")
            lines.append(f"- repaired: **{v.get('repaired')}**")
            lines.append(f"- Stage-3: `{v.get('stage3')}`")
            lines.append(f"- Stage-4: `{v.get('stage4')}`")
            lines.append("")
    lines.append("")
    lines.append("Do **not** launch Stage-5 automatically.")
    lines.append("")
    md.write_text("\n".join(lines))
    print("wrote", md, flush=True)


if __name__ == "__main__":
    main()
