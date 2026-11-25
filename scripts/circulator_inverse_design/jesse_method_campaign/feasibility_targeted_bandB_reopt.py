#!/usr/bin/env python3
"""Targeted Band-B middle-ring density reoptimization (NOT Stage-6).

Starts from Band-B C6 Stage-5 q with q[4] -= 0.15. Frees a small subset of
middle-ring C6 orbits identified by the ~5.16 GHz resonance-localization
sensitivity study. Protects the full 4.75–5.25 GHz band via a curated
frequency set (not a single-frequency fit).

Optimizer surrogate: existing A2 softmin hybrid J_BB (τ=0.05, α=0.25).
Selection / acceptance: hard worst-case A2 on the protective grid, then
dense 1.5625 MHz audit. Geometry / feed / B / solver unchanged.

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

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
FEAS = OUT / "feasibility"
TDIR = FEAS / "targeted_bandB_reoptimization"

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

from a2_efficient_core import eval_bb_A2 as eval_bb_A2_core  # noqa: E402
from audit_bb_stage5_dense import row_from_bundle  # noqa: E402
from multimode_opt_core import B_PLUS, S_HI, S_LO, orbits_and_masks  # noqa: E402
from objective_A2 import A2_FREEZE  # noqa: E402

OPT_PATH = FEAS / "B_4p75_5p25_C6_B0.0500_stage5_opt.json"
AUDIT_ORIG = FEAS / "audit_B_4p75_5p25_C6_B0.0500_stage5_opt_1.5625MHz.json"

BAND_F0, BAND_F1 = 4.75, 5.25
B_TESLA = 0.05
TAU, ALPHA = 0.05, 0.25
BOUNDS = (S_LO, S_HI)
# Sensitivity-justified free middle/inner orbits (resonance localization)
FREE_ORBITS = [4, 6, 5, 8, 3]  # prioritize 4; 6/5/8/3 from ±0.15 study
Q4_DELTA = -0.15
# Keep the search *local* — global [S_LO,S_HI] lets L-BFGS slam to corners.
LOCAL_BOX = 0.25
BASELINE_WORST_A2 = -4.964  # dense-audit floor of original Band-B Stage-5


def _env(omp: int = 1):
    for k, v in (
        ("OMP_NUM_THREADS", str(omp)),
        ("OPENBLAS_NUM_THREADS", str(omp)),
        ("MKL_NUM_THREADS", str(omp)),
        ("PYTHONUNBUFFERED", "1"),
        ("PLASMEEP_A2_FAST_KERNEL", "1"),
    ):
        os.environ.setdefault(k, v)


def protective_freqs() -> list[float]:
    """Curated Band-B grid: problem families + goods + edges (not full dense)."""
    freqs = set()
    # coarse band coverage 25 MHz
    f = BAND_F0
    while f <= BAND_F1 + 1e-12:
        freqs.add(round(f, 6))
        f += 0.025
    # ~4.909 family
    for x in (4.900, 4.90625, 4.909375, 4.9125, 4.91875, 4.925):
        freqs.add(x)
    # ~5.13 / 5.16 families + shift neighborhood
    for x in (
        5.120,
        5.1325,
        5.145,
        5.151562,
        5.157812,
        5.164062,
        5.170312,
        5.176562,
        5.182812,
        5.195,
        5.2075,
        5.226562,
    ):
        freqs.add(x)
    # representative goods
    for x in (4.85, 4.967187, 4.992188, 5.0375, 5.059375):
        freqs.add(x)
    return sorted(f for f in freqs if BAND_F0 - 1e-9 <= f <= BAND_F1 + 1e-9)


def atomic_write(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2) + "\n")
    tmp.replace(path)


def make_q0(q_stage5: np.ndarray) -> np.ndarray:
    q = np.asarray(q_stage5, float).copy()
    q[4] = float(np.clip(q[4] + Q4_DELTA, S_LO, S_HI))
    return q


def pack_free(q_full: np.ndarray, free: list[int], x: np.ndarray) -> np.ndarray:
    q = np.asarray(q_full, float).copy()
    for i, v in zip(free, x):
        q[i] = float(np.clip(v, S_LO, S_HI))
    return q


def unpack_free(q_full: np.ndarray, free: list[int]) -> np.ndarray:
    return np.asarray([q_full[i] for i in free], float)


# ---------------------------------------------------------------------------
# Optimization
# ---------------------------------------------------------------------------

def run_optimization(args):
    _env(args.omp_threads)
    TDIR.mkdir(parents=True, exist_ok=True)
    orbits, _, _ = orbits_and_masks()
    opt = json.loads(OPT_PATH.read_text())
    q_stage5 = np.asarray(opt["q_best"], float)
    q_start = make_q0(q_stage5)
    free = list(FREE_ORBITS)
    freqs = protective_freqs()
    freqs_hz = [f * 1e9 for f in freqs]

    meta = {
        "stage6_status": "NOT STARTED",
        "study": "targeted_bandB_middle_ring_reopt",
        "band_GHz": [BAND_F0, BAND_F1],
        "b_tesla": B_TESLA,
        "feed": "original_40mm",
        "geometry_unchanged": True,
        "free_orbits": free,
        "free_rationale": (
            "Resonance-localization ±0.15 study: orbit04_minus best combined "
            "(Δf≈−32 MHz, ΔA2≈+1.46, mild sentinels); orbit06/05/08/03 related "
            "middle/inner controls. Other 11 C6 vars frozen."
        ),
        "q4_delta_start": Q4_DELTA,
        "q_stage5": q_stage5.tolist(),
        "q_start": q_start.tolist(),
        "protective_freqs_GHz": freqs,
        "n_protective_freq": len(freqs),
        "objective_surrogate": {
            "name": "J_BB = A2 softmin hybrid",
            "tau": TAU,
            "alpha": ALPHA,
            "A2": A2_FREEZE,
            "note": (
                "L-BFGS-B maximizes J_BB on the protective frequency set. "
                "J_BB is a smooth surrogate for worst-case A2, not identical to "
                "min_f A2(f). Final candidate selected by hard worst_J on that "
                "set; acceptance requires dense-audit worst A2 ≥ baseline."
            ),
        },
        "baseline_dense_worst_A2": BASELINE_WORST_A2,
        "bounds_global": list(BOUNDS),
        "local_box": LOCAL_BOX,
        "maxiter": args.maxiter,
        "n_workers": args.n_workers,
    }
    atomic_write(TDIR / "opt_meta.json", meta)

    n_eval = 0
    accepted_nit = 0
    history = []
    cache = {"q": None, "rec": None}
    # Track best by hard worst-case A2 (primary scientific criterion)
    best_worst = {"worst_J": -1e300, "J_BB": -1e300, "q": None, "rec": None}
    # Also track best J_BB for diagnostics
    best_jbb = {"J_BB": -1e300, "worst_J": -1e300, "q": None, "rec": None}
    t_wall0 = time.perf_counter()
    ckpt = TDIR / "targeted_opt_ckpt.json"

    def _eval_full(q, track: bool = True):
        nonlocal n_eval, cache
        q = np.asarray(q, float)
        if cache["q"] is not None and np.allclose(q, cache["q"], rtol=0, atol=0):
            return cache["rec"]
        n_eval += 1
        t0 = time.perf_counter()
        rec = eval_bb_A2_core(
            args.grade,
            q,
            orbits,
            freqs_hz,
            tau=TAU,
            alpha=ALPHA,
            n_workers=args.n_workers,
            omp_threads_per_worker=args.omp_threads,
            use_fast_kernel=True,
            b_tesla=B_TESLA,
        )
        wall = time.perf_counter() - t0
        compact = {
            "J_BB": float(rec["J_BB"]),
            "J_mean": float(rec["J_mean"]),
            "worst_J": float(rec["worst_J"]),
            "worst_f_GHz": float(rec["worst_f_GHz"]),
            "grad_norm_full": float(rec["grad_norm"]),
            "grad_q": np.asarray(rec["grad_q"], float),
            "wall_s": float(wall),
            "freqs_GHz": list(rec["freqs_GHz"]),
            "J_per_freq": [float(x) for x in rec["J_per_freq"]],
            "top5_weights": rec.get("top5_weights"),
            "worst_metrics": {
                "T_des_avg": float(rec["worst_metrics"]["T_des_avg"]),
                "IL_avg_dB": float(rec["worst_metrics"]["IL_avg_dB"]),
                "isolation_avg_dB": float(rec["worst_metrics"]["isolation_avg_dB"]),
                "R_avg": float(rec["worst_metrics"]["R_avg"]),
                "TE1_avg": float(rec["worst_metrics"]["TE1_avg"]),
            },
        }
        print(
            f"EVAL #{n_eval} J_BB={compact['J_BB']:.4f} worstA2={compact['worst_J']:.4f} "
            f"@{compact['worst_f_GHz']:.5f} wall={wall:.1f}s "
            f"free={[round(float(q[i]),4) for i in free]}",
            flush=True,
        )
        cache = {"q": q.copy(), "rec": compact}

        if track:
            if compact["worst_J"] > best_worst["worst_J"] + 1e-9 or (
                abs(compact["worst_J"] - best_worst["worst_J"]) < 1e-9
                and compact["J_BB"] > best_worst["J_BB"]
            ):
                best_worst.update(
                    {
                        "worst_J": compact["worst_J"],
                        "J_BB": compact["J_BB"],
                        "q": q.copy(),
                        "rec": {k: v for k, v in compact.items() if k != "grad_q"},
                    }
                )
                print(f"  NEW BEST-by-worstA2: {compact['worst_J']:.4f}", flush=True)
            if compact["J_BB"] > best_jbb["J_BB"]:
                best_jbb.update(
                    {
                        "J_BB": compact["J_BB"],
                        "worst_J": compact["worst_J"],
                        "q": q.copy(),
                        "rec": {k: v for k, v in compact.items() if k != "grad_q"},
                    }
                )
        return compact

    def fun(x):
        q = pack_free(q_start, free, x)
        return -_eval_full(q)["J_BB"]

    def jac(x):
        q = pack_free(q_start, free, x)
        g = _eval_full(q)["grad_q"]
        return -np.asarray([g[i] for i in free], float)

    def callback(xk):
        nonlocal accepted_nit
        accepted_nit += 1
        q = pack_free(q_start, free, xk)
        rec = _eval_full(q)
        history.append(
            {
                "kind": "accepted_iter",
                "iteration": accepted_nit,
                "J_BB": rec["J_BB"],
                "worst_A2": rec["worst_J"],
                "worst_f_GHz": rec["worst_f_GHz"],
                "q_free": {str(i): float(q[i]) for i in free},
                "q": q.tolist(),
            }
        )
        print(
            f"ACCEPTED iter={accepted_nit}/{args.maxiter} "
            f"J_BB={rec['J_BB']:.4f} worstA2={rec['worst_J']:.4f} "
            f"@{rec['worst_f_GHz']:.5f}",
            flush=True,
        )
        atomic_write(
            ckpt,
            {
                "status": "in_progress",
                "accepted_nit": accepted_nit,
                "n_eval": n_eval,
                "elapsed_s": time.perf_counter() - t_wall0,
                "best_by_worst": {
                    "worst_J": best_worst["worst_J"],
                    "J_BB": best_worst["J_BB"],
                    "q": None if best_worst["q"] is None else best_worst["q"].tolist(),
                    "rec": best_worst["rec"],
                },
                "best_by_J_BB": {
                    "J_BB": best_jbb["J_BB"],
                    "worst_J": best_jbb["worst_J"],
                    "q": None if best_jbb["q"] is None else best_jbb["q"].tolist(),
                },
                "history": history,
                "stage6_status": "NOT STARTED",
            },
        )

    print("=" * 72, flush=True)
    print("TARGETED Band-B REOPT — Stage-6 NOT STARTED", flush=True)
    print(
        f"free orbits={free} n_freq_protect={len(freqs)} maxiter={args.maxiter} "
        f"workers={args.n_workers} B={B_TESLA} T",
        flush=True,
    )
    print(
        f"start q[4]={q_start[4]:.4f} (stage5 {q_stage5[4]:.4f}{Q4_DELTA:+.2f}); "
        f"surrogate=J_BB; select=worst_J",
        flush=True,
    )
    print("=" * 72, flush=True)

    # baselines on protective grid (stage5 is reference only — do not enter candidate pool)
    rec_stage5 = _eval_full(q_stage5, track=False)
    history.append(
        {
            "kind": "stage5_baseline_protective",
            "J_BB": rec_stage5["J_BB"],
            "worst_A2": rec_stage5["worst_J"],
            "worst_f_GHz": rec_stage5["worst_f_GHz"],
        }
    )
    rec0 = _eval_full(q_start, track=True)  # initializes best_* from q4−0.15 start
    history.append(
        {
            "kind": "q4_minus_start_protective",
            "J_BB": rec0["J_BB"],
            "worst_A2": rec0["worst_J"],
            "worst_f_GHz": rec0["worst_f_GHz"],
            "q": q_start.tolist(),
        }
    )
    print(
        f"STAGE5 protect: J_BB={rec_stage5['J_BB']:.4f} worst={rec_stage5['worst_J']:.4f} "
        f"@{rec_stage5['worst_f_GHz']:.5f}",
        flush=True,
    )
    print(
        f"START q4-0.15: J_BB={rec0['J_BB']:.4f} worst={rec0['worst_J']:.4f} "
        f"@{rec0['worst_f_GHz']:.5f}",
        flush=True,
    )

    x0 = unpack_free(q_start, free)
    # Local box around the q4−0.15 start (clipped to global density bounds)
    local_bounds = []
    for i, x in zip(free, x0):
        lo = max(S_LO, float(x) - LOCAL_BOX)
        hi = min(S_HI, float(x) + LOCAL_BOX)
        local_bounds.append((lo, hi))
    print(f"Local box ±{LOCAL_BOX} bounds: {list(zip(free, local_bounds))}", flush=True)

    res = minimize(
        fun,
        x0,
        method="L-BFGS-B",
        jac=jac,
        bounds=local_bounds,
        callback=callback,
        options={
            "maxiter": args.maxiter,
            "maxfun": args.maxiter * 8,
            "ftol": 1e-10,
            "maxls": 20,
        },
    )
    elapsed = time.perf_counter() - t_wall0

    # Prefer best-by-worst if it beats start; else keep start
    q_cand = best_worst["q"] if best_worst["q"] is not None else q_start
    rec_cand = best_worst["rec"]
    if rec_cand is None or rec_cand["worst_J"] < rec0["worst_J"] - 1e-6:
        # optimizer failed to improve hard worst — keep q4-0.15 start
        q_cand = q_start
        rec_cand = {k: v for k, v in rec0.items() if k != "grad_q"}
        selection = "fallback_q4_minus_start_no_worst_improvement"
    else:
        selection = "best_by_protective_worst_A2"

    out = {
        "status": "complete",
        "stage6_status": "NOT STARTED",
        "selection": selection,
        "optimizer_success": bool(res.success),
        "optimizer_message": str(res.message),
        "accepted_nit": accepted_nit,
        "n_eval": n_eval,
        "elapsed_s": elapsed,
        "free_orbits": free,
        "q_stage5": q_stage5.tolist(),
        "q_q4_minus": q_start.tolist(),
        "q_best": q_cand.tolist(),
        "delta_q_vs_stage5": (np.asarray(q_cand) - q_stage5).tolist(),
        "delta_q_free_vs_start": {
            str(i): float(q_cand[i] - q_start[i]) for i in free
        },
        "protective_stage5": {
            "J_BB": rec_stage5["J_BB"],
            "worst_A2": rec_stage5["worst_J"],
            "worst_f_GHz": rec_stage5["worst_f_GHz"],
        },
        "protective_q4_minus": {
            "J_BB": rec0["J_BB"],
            "worst_A2": rec0["worst_J"],
            "worst_f_GHz": rec0["worst_f_GHz"],
        },
        "protective_best": rec_cand,
        "best_by_J_BB": {
            "J_BB": best_jbb["J_BB"],
            "worst_J": best_jbb["worst_J"],
            "q": None if best_jbb["q"] is None else best_jbb["q"].tolist(),
        },
        "history": history,
        "meta": meta,
    }
    atomic_write(TDIR / "targeted_opt.json", out)
    print(
        f"OPT DONE selection={selection} worstA2={rec_cand['worst_J']:.4f} "
        f"@{rec_cand['worst_f_GHz']:.5f} n_eval={n_eval} elapsed={elapsed/60:.1f} min",
        flush=True,
    )
    return out


# ---------------------------------------------------------------------------
# Dense audits
# ---------------------------------------------------------------------------

def _init_worker(omp: int):
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = str(int(omp))
    os.environ.setdefault("PLASMEEP_A2_FAST_KERNEL", "1")


def _audit_worker(payload: dict) -> dict:
    _init_worker(int(payload["omp_threads"]))
    from eval_forward_only import eval_q_forward

    q = np.asarray(payload["q"], float)
    orbits, _, _ = orbits_and_masks()
    bun, _ = eval_q_forward(
        payload["grade"],
        q,
        orbits,
        f_hz=float(payload["f_GHz"]) * 1e9,
        b_tesla=float(payload["b_tesla"]),
    )
    row = row_from_bundle(float(payload["f_GHz"]), bun)
    row["design"] = payload["design"]
    return row


def freq_grid(f0, f1, df):
    n = int(round((f1 - f0) / df)) + 1
    return [round(f0 + i * df, 6) for i in range(n)]


def dense_audit(design: str, q: np.ndarray, n_workers: int, omp: int, df_MHz: float = 1.5625):
    freqs = freq_grid(BAND_F0, BAND_F1, df_MHz * 1e-3)
    out_path = TDIR / f"audit_{design}_{df_MHz:g}MHz.json"
    # resume
    done = {}
    if out_path.exists():
        prev = json.loads(out_path.read_text())
        for r in prev.get("rows", []):
            done[round(float(r["f_GHz"]), 6)] = r
    jobs = []
    for f in freqs:
        if round(f, 6) in done:
            continue
        jobs.append(
            {
                "design": design,
                "q": q.tolist(),
                "f_GHz": f,
                "b_tesla": B_TESLA,
                "grade": "C",
                "omp_threads": omp,
            }
        )
    print(
        f"DENSE AUDIT {design}: {len(freqs)} freqs, {len(jobs)} remaining, "
        f"n_workers={n_workers}",
        flush=True,
    )
    t0 = time.time()
    rows = list(done.values())
    if jobs:
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            futs = {ex.submit(_audit_worker, j): j for j in jobs}
            done_n = 0
            for fut in as_completed(futs):
                done_n += 1
                j = futs[fut]
                try:
                    r = fut.result()
                except Exception as e:
                    r = {"design": design, "f_GHz": j["f_GHz"], "error": repr(e)}
                    print(f"FAIL {design} f={j['f_GHz']}: {e}", flush=True)
                rows.append(r)
                if done_n % 20 == 0 or done_n == len(jobs):
                    print(
                        f"  [{done_n}/{len(jobs)}] {design} f={r.get('f_GHz')} "
                        f"A2={r.get('J_A2', 'err')}",
                        flush=True,
                    )
                    atomic_write(
                        out_path,
                        {
                            "design": design,
                            "df_MHz": df_MHz,
                            "rows": sorted(
                                [x for x in rows if "f_GHz" in x],
                                key=lambda x: x["f_GHz"],
                            ),
                            "elapsed_s": time.time() - t0,
                            "partial": True,
                        },
                    )
    rows = sorted([r for r in rows if "J_A2" in r], key=lambda r: r["f_GHz"])
    payload = {
        "design": design,
        "df_MHz": df_MHz,
        "n": len(rows),
        "elapsed_s": time.time() - t0,
        "q": q.tolist(),
        "rows": rows,
        "partial": False,
        "stage6_status": "NOT STARTED",
    }
    atomic_write(out_path, payload)
    return payload


def summarize_audit(design: str, rows: list[dict], q: np.ndarray, q_ref: np.ndarray | None = None):
    ok = [r for r in rows if "J_A2" in r]
    a2 = np.array([r["J_A2"] for r in ok])
    worst = min(ok, key=lambda r: r["J_A2"])
    # channel fields may be nested differently in row_from_bundle
    def get(r, *keys, default=np.nan):
        for k in keys:
            if k in r and r[k] is not None:
                return float(r[k])
        return float(default)

    # Prefer modal fields if present
    # Dense-audit rows use short names (T_des, R, TE1); opt metrics use *_avg.
    des = np.array([get(r, "T_des", "T_des_avg", "TE0_des") for r in ok])
    rev = np.array([get(r, "T_rev", "T_rev_avg") for r in ok])
    R = np.array([get(r, "R", "R_avg", "R_avg_te0") for r in ok])
    te1 = np.array([get(r, "TE1", "TE1_avg", "TE1_total") for r in ok])
    pout = np.array([get(r, "P_out", "P_out_avg") for r in ok])

    wrong = []
    for r in ok:
        if "TE0_wrong" in r:
            wrong.append(float(r["TE0_wrong"]))
        else:
            wrong.append(get(r, "T_rev", "T_rev_avg", default=0.0))
    wrong = np.array(wrong)

    dq = None
    if q_ref is not None:
        dq = (np.asarray(q) - np.asarray(q_ref)).tolist()

    return {
        "design": design,
        "n": len(ok),
        "worst_A2": float(worst["J_A2"]),
        "worst_f_GHz": float(worst["f_GHz"]),
        "mean_A2": float(np.mean(a2)),
        "n_A2_le_m5": int(np.sum(a2 <= -5.0)),
        "n_A2_le_m7": int(np.sum(a2 <= -7.0)),
        "min_T_des": float(np.nanmin(des)),
        "mean_T_des": float(np.nanmean(des)),
        "mean_T_rev": float(np.nanmean(rev)),
        "mean_wrong_TE0_proxy": float(np.nanmean(wrong)),
        "mean_R": float(np.nanmean(R)),
        "mean_TE1": float(np.nanmean(te1)),
        "mean_P_out": float(np.nanmean(pout)),
        "min_P_out": float(np.nanmin(pout)),
        "q": q.tolist(),
        "delta_q_vs_stage5": dq,
        "worst_row": {k: worst[k] for k in worst if k not in ("S",)},
    }


def load_or_copy_orig_audit():
    """Reuse existing 1.5625 MHz Stage-5 audit if present."""
    if AUDIT_ORIG.exists():
        d = json.loads(AUDIT_ORIG.read_text())
        rows = d.get("rows") or d.get("audit_rows") or []
        # write into TDIR for self-contained folder
        opt = json.loads(OPT_PATH.read_text())
        q = np.asarray(opt["q_best"], float)
        payload = {
            "design": "orig_stage5",
            "df_MHz": 1.5625,
            "n": len(rows),
            "rows": rows,
            "q": q.tolist(),
            "source": str(AUDIT_ORIG),
            "copied": True,
        }
        atomic_write(TDIR / "audit_orig_stage5_1.5625MHz.json", payload)
        return payload
    return None


# ---------------------------------------------------------------------------
# Compare / plot / report
# ---------------------------------------------------------------------------

def make_plots(audits: dict, summaries: dict):
    fig, ax = plt.subplots(figsize=(10, 5.2))
    styles = {
        "orig_stage5": ("C0", "-", "Original Stage-5"),
        "q4_minus": ("C1", "--", "q[4]−0.15 (no opt)"),
        "targeted_opt": ("C2", "-", "Targeted local opt"),
    }
    for key, (color, ls, label) in styles.items():
        rows = sorted(
            [r for r in audits[key]["rows"] if "J_A2" in r],
            key=lambda r: r["f_GHz"],
        )
        if not rows:
            continue
        ax.plot(
            [r["f_GHz"] for r in rows],
            [r["J_A2"] for r in rows],
            color=color,
            ls=ls,
            lw=1.4,
            label=label,
        )
        s = summaries[key]
        ax.axvline(s["worst_f_GHz"], color=color, ls=":", lw=0.9, alpha=0.8)
        ax.scatter([s["worst_f_GHz"]], [s["worst_A2"]], color=color, s=40, zorder=5)
    ax.axhline(BASELINE_WORST_A2, color="k", ls=":", lw=0.8, alpha=0.5, label="baseline worst")
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("A2")
    ax.set_title("Band-B dense A2 (1.5625 MHz) — targeted reopt comparison")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(TDIR / "A2_spectra_three_designs.png", dpi=150)
    plt.close(fig)

    # modal overlays if fields exist
    for metrics, ylabel, fname in [
        (("T_des", "T_des_avg"), "Desired TE0", "Tdes_spectra_three_designs.png"),
        (("R", "R_avg"), "Reflection R", "R_spectra_three_designs.png"),
        (("TE1", "TE1_avg"), "TE1", "TE1_spectra_three_designs.png"),
    ]:
        fig, ax = plt.subplots(figsize=(10, 4.5))
        plotted = False
        for key, (color, ls, label) in styles.items():
            rows = sorted(
                [r for r in audits[key]["rows"] if "J_A2" in r],
                key=lambda r: r["f_GHz"],
            )
            ys = []
            xs = []
            for r in rows:
                val = next((r[m] for m in metrics if m in r), None)
                if val is None:
                    continue
                xs.append(r["f_GHz"])
                ys.append(val)
            if not xs:
                continue
            plotted = True
            ax.plot(xs, ys, color=color, ls=ls, lw=1.2, label=label)
        if plotted:
            ax.set_xlabel("Frequency (GHz)")
            ax.set_ylabel(ylabel)
            ax.set_title(f"{ylabel} — three-design comparison")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            fig.savefig(TDIR / fname, dpi=140)
        plt.close(fig)


def write_report(opt_out, summaries, science, elapsed_opt, elapsed_audit):
    lines = []
    lines.append("# Targeted Band-B Middle-Ring Density Reoptimization")
    lines.append("")
    lines.append("**Stage-6 status: NOT STARTED. No 1 GHz / 2 GHz campaign.**")
    lines.append("")
    lines.append(
        "Start: Stage-5 `q` with `q[4] -= 0.15`. Free orbits: "
        f"{FREE_ORBITS}. Surrogate objective: softmin hybrid **J_BB** "
        f"(τ={TAU}, α={ALPHA}) on a protective Band-B frequency set. "
        "Candidate selection by hard worst A2 on that set; acceptance by "
        "dense 1.5625 MHz audit vs baseline."
    )
    lines.append("")
    lines.append("## Dense-audit comparison (1.5625 MHz)")
    lines.append("")
    lines.append(
        "| Design | Worst A2 | Worst f (GHz) | Mean A2 | #≤−5 | #≤−7 | "
        "min Tdes | mean Tdes | mean R | mean TE1 | mean Pout |"
    )
    lines.append(
        "|--------|---------:|--------------:|--------:|-----:|-----:|"
        "---------:|---------:|-------:|---------:|----------:|"
    )
    order = ["orig_stage5", "q4_minus", "targeted_opt"]
    labels = {
        "orig_stage5": "Original Stage-5",
        "q4_minus": "q[4]−0.15 (no opt)",
        "targeted_opt": "Targeted local opt",
    }
    for k in order:
        s = summaries[k]
        lines.append(
            f"| {labels[k]} | {s['worst_A2']:.3f} | {s['worst_f_GHz']:.6f} | "
            f"{s['mean_A2']:.3f} | {s['n_A2_le_m5']} | {s['n_A2_le_m7']} | "
            f"{s['min_T_des']:.3f} | {s['mean_T_des']:.3f} | {s['mean_R']:.3f} | "
            f"{s['mean_TE1']:.3f} | {s['mean_P_out']:.3f} |"
        )
    lines.append("")
    if opt_out:
        lines.append("## Optimization summary")
        lines.append("")
        lines.append(f"- Selection rule: `{opt_out.get('selection')}`")
        lines.append(
            f"- L-BFGS accepted iters: {opt_out.get('accepted_nit')}, "
            f"evals: {opt_out.get('n_eval')}, wall: {elapsed_opt/60:.1f} min"
        )
        lines.append(f"- Optimizer message: {opt_out.get('optimizer_message')}")
        dq = opt_out.get("delta_q_free_vs_start") or {}
        lines.append("- Δq (free) vs q4−0.15 start: " + ", ".join(f"q[{k}]={v:+.4f}" for k, v in dq.items()))
        lines.append(
            f"- Protective-grid worst at best: "
            f"{opt_out.get('protective_best', {}).get('worst_J')} @ "
            f"{opt_out.get('protective_best', {}).get('worst_f_GHz')} GHz"
        )
        lines.append("")
    lines.append("## Scientific questions")
    lines.append("")
    for i in range(1, 7):
        key = f"Q{i}"
        block = science[key]
        lines.append(f"### {key} — {block['question']}")
        lines.append(f"**{block['answer']}**")
        for n in block.get("notes", []):
            lines.append(f"- {n}")
        lines.append("")
    lines.append("Do **not** expand to all 16 variables or launch Stage-6 from this script.")
    lines.append("")
    text = "\n".join(lines) + "\n"
    (TDIR / "TARGETED_BANDB_REOPT_REPORT.md").write_text(text)
    (TDIR / "TARGETED_BANDB_REOPT_EXEC.txt").write_text(
        "\n".join(
            [
                "Targeted Band-B middle-ring reopt. Stage-6 NOT STARTED.",
                f"orig worst={summaries['orig_stage5']['worst_A2']:.3f}",
                f"q4-0.15 worst={summaries['q4_minus']['worst_A2']:.3f}",
                f"targeted worst={summaries['targeted_opt']['worst_A2']:.3f}",
                f"Q1 improve? {science['Q1']['answer']}",
                f"Q6 next: {science['Q6']['answer']}",
            ]
        )
        + "\n"
    )
    return text


def answer_science(summaries, opt_out) -> dict:
    o = summaries["orig_stage5"]
    p = summaries["q4_minus"]
    t = summaries["targeted_opt"]
    improved = t["worst_A2"] > o["worst_A2"] + 0.05
    accepted = t["worst_A2"] >= o["worst_A2"] - 1e-6
    # Did we suppress or shift 5.164?
    # Look at A2 at 5.164 vs new worst frequency
    shift = abs(t["worst_f_GHz"] - 5.164062) > 0.01
    still_near_516 = abs(t["worst_f_GHz"] - 5.164062) < 0.02
    new_limiter = abs(t["worst_f_GHz"] - o["worst_f_GHz"]) > 0.02

    delta = t["worst_A2"] - o["worst_A2"]
    if improved and delta >= 1.0:
        robust = "yes_substantially"
    elif improved and delta >= 0.3:
        robust = "modestly"
    elif accepted and delta > -0.05:
        robust = "not_really_flat"
    else:
        robust = "no"

    if improved:
        density_ctrl = (
            "plasma_density_controllable_locally"
            if delta >= 0.3
            else "weak_density_control_so_far"
        )
    elif p["worst_A2"] > o["worst_A2"]:
        density_ctrl = "perturbation_helps_but_reopt_did_not_improve_further"
    else:
        density_ctrl = "density_knob_insufficient_consider_geometry"

    if improved:
        next_exp = (
            "Repeat dense audit interpretation; if worst is a new family, run a "
            "small localization on the new limiter before any wider opt."
        )
    elif not accepted:
        next_exp = (
            "Do not expand to 16-var Band-B reopt yet. Either (a) probe a "
            "different single-orbit direction (e.g. q[8]+ or q[5]−) with the "
            "same protective protocol, or (b) revisit geometry only if multiple "
            "single-orbit attempts fail to raise the dense floor."
        )
    else:
        next_exp = (
            "Hold geometry fixed; try one more modest free-set tweak "
            "(add/remove one orbit) only if the protective-grid worst clearly "
            "moved but dense audit did not improve."
        )

    return {
        "Q1": {
            "question": "Did targeted reopt improve worst-case Band-B A2?",
            "answer": "yes" if improved else ("no_accepted_tie" if accepted else "no"),
            "notes": [
                f"orig worst={o['worst_A2']:.3f} @{o['worst_f_GHz']:.4f}",
                f"q4−0.15 worst={p['worst_A2']:.3f} @{p['worst_f_GHz']:.4f}",
                f"targeted worst={t['worst_A2']:.3f} @{t['worst_f_GHz']:.4f}",
                f"Δ vs orig = {delta:+.3f}",
            ],
        },
        "Q2": {
            "question": "Suppress original resonance or merely shift it?",
            "answer": (
                "shifted_and_or_reshaped"
                if shift
                else ("still_near_5p164" if still_near_516 else "unclear")
            ),
            "notes": [
                f"new worst frequency {t['worst_f_GHz']:.6f} GHz "
                f"(orig {o['worst_f_GHz']:.6f})",
            ],
        },
        "Q3": {
            "question": "Did another resonance become the new limiting feature?",
            "answer": "yes" if (new_limiter and not still_near_516) else "no_or_same_family",
            "notes": [
                f"worst moved by {abs(t['worst_f_GHz']-o['worst_f_GHz'])*1e3:.1f} MHz",
            ],
        },
        "Q4": {
            "question": "Can 0.5 GHz band become substantially more robust than −4.96?",
            "answer": robust,
            "notes": [f"achieved dense floor {t['worst_A2']:.3f}"],
        },
        "Q5": {
            "question": "Density-controllable vs geometry needed?",
            "answer": density_ctrl,
            "notes": [],
        },
        "Q6": {
            "question": "Smallest justified next experiment?",
            "answer": next_exp,
            "notes": ["Do not start Stage-6 or 1–2 GHz optimization from this result."],
        },
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--omp-threads", type=int, default=1)
    ap.add_argument("--grade", default="C")
    ap.add_argument("--maxiter", type=int, default=12)
    ap.add_argument("--skip-opt", action="store_true")
    ap.add_argument("--skip-audit", action="store_true")
    ap.add_argument("--analyze-only", action="store_true")
    args = ap.parse_args()

    TDIR.mkdir(parents=True, exist_ok=True)
    opt_out = None
    elapsed_opt = 0.0

    if args.analyze_only:
        opt_out = json.loads((TDIR / "targeted_opt.json").read_text())
    elif not args.skip_opt:
        t0 = time.time()
        opt_out = run_optimization(args)
        elapsed_opt = time.time() - t0
    else:
        opt_out = json.loads((TDIR / "targeted_opt.json").read_text())

    opt_s5 = json.loads(OPT_PATH.read_text())
    q_stage5 = np.asarray(opt_s5["q_best"], float)
    q_q4 = make_q0(q_stage5)
    q_targ = np.asarray(opt_out["q_best"], float)

    audits = {}
    elapsed_audit = 0.0
    if not args.skip_audit and not args.analyze_only:
        t0 = time.time()
        # orig: prefer copy
        copied = load_or_copy_orig_audit()
        if copied and copied.get("rows"):
            audits["orig_stage5"] = copied
        else:
            audits["orig_stage5"] = dense_audit(
                "orig_stage5", q_stage5, args.n_workers, args.omp_threads
            )
        audits["q4_minus"] = dense_audit(
            "q4_minus", q_q4, args.n_workers, args.omp_threads
        )
        audits["targeted_opt"] = dense_audit(
            "targeted_opt", q_targ, args.n_workers, args.omp_threads
        )
        elapsed_audit = time.time() - t0
        atomic_write(
            TDIR / "audits_bundle.json",
            {
                "elapsed_s": elapsed_audit,
                "designs": {k: {"n": v["n"], "path": f"audit_{k}_1.5625MHz.json"} for k, v in audits.items()},
            },
        )
    else:
        for k in ("orig_stage5", "q4_minus", "targeted_opt"):
            audits[k] = json.loads((TDIR / f"audit_{k}_1.5625MHz.json").read_text())

    summaries = {
        "orig_stage5": summarize_audit("orig_stage5", audits["orig_stage5"]["rows"], q_stage5),
        "q4_minus": summarize_audit("q4_minus", audits["q4_minus"]["rows"], q_q4, q_stage5),
        "targeted_opt": summarize_audit(
            "targeted_opt", audits["targeted_opt"]["rows"], q_targ, q_stage5
        ),
    }
    # acceptance gate
    if summaries["targeted_opt"]["worst_A2"] < summaries["orig_stage5"]["worst_A2"] - 1e-6:
        summaries["targeted_opt"]["accepted_vs_baseline"] = False
        print(
            "REJECT targeted_opt vs baseline on dense worst A2 "
            f"({summaries['targeted_opt']['worst_A2']:.3f} < "
            f"{summaries['orig_stage5']['worst_A2']:.3f})",
            flush=True,
        )
    else:
        summaries["targeted_opt"]["accepted_vs_baseline"] = True

    science = answer_science(summaries, opt_out)
    make_plots(audits, summaries)
    report = write_report(opt_out, summaries, science, elapsed_opt, elapsed_audit)
    atomic_write(
        TDIR / "comparison_summary.json",
        {
            "summaries": summaries,
            "science": science,
            "elapsed_opt_s": elapsed_opt,
            "elapsed_audit_s": elapsed_audit,
            "stage6_status": "NOT STARTED",
        },
    )
    print(report, flush=True)
    print(
        "\n=== KEY ===\n"
        f"orig worst={summaries['orig_stage5']['worst_A2']:.3f} @ {summaries['orig_stage5']['worst_f_GHz']:.4f}\n"
        f"q4-0.15 worst={summaries['q4_minus']['worst_A2']:.3f} @ {summaries['q4_minus']['worst_f_GHz']:.4f}\n"
        f"targeted worst={summaries['targeted_opt']['worst_A2']:.3f} @ {summaries['targeted_opt']['worst_f_GHz']:.4f}\n"
        f"Q1={science['Q1']['answer']}\n"
        f"Q6={science['Q6']['answer']}\n"
        "Stage-6 NOT STARTED.",
        flush=True,
    )


if __name__ == "__main__":
    main()
