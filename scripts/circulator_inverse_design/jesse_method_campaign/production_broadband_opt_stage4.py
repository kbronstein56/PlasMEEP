#!/usr/bin/env python3
"""Stage-4 4–6 GHz broadband optimization per BB_STAGE4_PROPOSAL.md (AUTHORIZED).

Frozen:
  C6, Objective A2_hybrid, soft-min τ=0.05 α=0.25,
  47-frequency localized grid (12.5 MHz local ladder + 6.25 MHz core at 5.8),
  FULL91-C, L-BFGS-B,
  warm-start from opt_broadband_stage3_A2_C.json q_best,
  process-parallel frequencies (n_workers=12, omp=1).

Does not change FEM physics, A2, mesh, materials, derivatives, ports, or B sign.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

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

from a2_efficient_core import eval_bb_A2 as eval_bb_A2_core  # noqa: E402
from multimode_opt_core import S_HI, S_LO, orbits_and_masks  # noqa: E402
from objective_A2 import A2_FREEZE  # noqa: E402

# FINAL localized Stage-4 set (user-authorized)
F_OBJ_GHZ = [
    4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.55, 4.6, 4.7, 4.75, 4.775, 4.8, 4.825, 4.85, 4.9,
    5.0, 5.1, 5.2, 5.3, 5.35, 5.4, 5.5, 5.5875, 5.6, 5.625, 5.65,
    5.7, 5.7125, 5.725, 5.7375, 5.75, 5.7625, 5.775, 5.7875,
    5.8, 5.80625, 5.8125, 5.81875, 5.825,
    5.8375, 5.85, 5.8625, 5.875, 5.8875, 5.9, 5.95, 6.0,
]
assert len(F_OBJ_GHZ) == 47
assert F_OBJ_GHZ == sorted(F_OBJ_GHZ)
assert 5.8125 in F_OBJ_GHZ and 4.775 in F_OBJ_GHZ
assert 5.80625 in F_OBJ_GHZ and 5.81875 in F_OBJ_GHZ

TAU_FREEZE = 0.05
ALPHA_FREEZE = 0.25
BOUNDS = (S_LO, S_HI)  # [0.05, 1.8]
WARM_START = "opt_broadband_stage3_A2_C.json"

CKPT_PATH = OUT / "bb_stage4_checkpoint.json"
HIST_PATH = OUT / "bb_stage4_history.json"
META_PATH = OUT / "bb_stage4_meta.json"
FINAL_PATH = OUT / "opt_broadband_stage4_A2_C.json"
EARLY_PATH = OUT / "bb_stage4_early_diagnostic.json"
FLAGS_PATH = OUT / "bb_stage4_early_flags.json"


def eval_bb_A2(grade, q, orbits, freqs_hz, tau=TAU_FREEZE, alpha=ALPHA_FREEZE,
               n_workers=1, omp_threads_per_worker=1):
    return eval_bb_A2_core(
        grade, q, orbits, freqs_hz,
        tau=tau, alpha=alpha,
        n_workers=n_workers,
        omp_threads_per_worker=omp_threads_per_worker,
        profile=False,
    )


def _eng_snapshot(met, J_A2):
    return {
        "J_A2": float(J_A2),
        "T_des": float(met["T_des_avg"]),
        "T_rev": float(met["T_rev_avg"]),
        "IL_dB": float(met["IL_avg_dB"]),
        "isolation_dB": float(met["isolation_avg_dB"]),
        "RL_dB": float(met["RL_avg_dB"]),
        "R": float(met["R_avg"]),
        "TE1": float(met["TE1_avg"]),
    }


def _snap_at_f(rec, f_ghz):
    freqs = np.asarray(rec["freqs_GHz"], float)
    i = int(np.argmin(np.abs(freqs - f_ghz)))
    if abs(freqs[i] - f_ghz) > 1e-6:
        return None
    return _eng_snapshot(rec["metrics_per_freq"][i], rec["J_per_freq"][i])


def _record_from_rec(rec, q, iteration, elapsed_s, kind, prev_elapsed=None):
    out = {
        "kind": kind,
        "iteration": int(iteration),
        "elapsed_s": float(elapsed_s),
        "wall_s_since_prev_accepted": (
            None if prev_elapsed is None else float(elapsed_s - prev_elapsed)
        ),
        "q": np.asarray(q, float).tolist(),
        "J_BB": float(rec["J_BB"]),
        "J_mean": float(rec["J_mean"]),
        "worst_A2": float(rec["worst_J"]),
        "worst_f_GHz": float(rec["worst_f_GHz"]),
        "grad_norm": float(rec["grad_norm"]),
        "max_weight": float(rec["max_weight"]),
        "N_eff": float(rec["N_eff"]),
        "top5_weights": rec["top5_weights"],
        "worst": _eng_snapshot(rec["worst_metrics"], rec["worst_J"]),
        "at_4p1": _eng_snapshot(rec["metrics_4p1"], rec["J_4p1"]),
        "at_4p775": _snap_at_f(rec, 4.775),
        "at_5p8125": _snap_at_f(rec, 5.8125),
        "at_5p75": _snap_at_f(rec, 5.75),
        "at_5p825": _snap_at_f(rec, 5.825),
        "J_per_freq": rec["J_per_freq"],
        "weights": rec["weights"],
        "freqs_GHz": rec["freqs_GHz"],
    }
    return out


def _atomic_write(path: Path, obj: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2) + "\n")
    tmp.replace(path)


def load_checkpoint():
    if not CKPT_PATH.exists():
        return None
    return json.loads(CKPT_PATH.read_text())


def evaluate_early_flags(early_rows, tdes0, a2_4p1_start=None):
    """Flag (do not abort) after 5 accepted iterations — preserve 4.1 recovery."""
    flags = []
    if len(early_rows) < 5:
        return {"n": len(early_rows), "flags": flags, "complete": False}
    first = early_rows[0]
    last = early_rows[4]
    tdes5 = last["at_4p1"]["T_des"]
    il5 = last["at_4p1"]["IL_dB"]
    iso5 = last["at_4p1"]["isolation_dB"]
    a2_5 = last["at_4p1"]["J_A2"]
    a2_0 = float(a2_4p1_start) if a2_4p1_start is not None else float(first["at_4p1"]["J_A2"])
    # Regression vs Stage-2 warm-start at 4.1: flag if T_des collapses badly
    if tdes5 < 0.5 * max(tdes0, 1e-12) and tdes5 < 0.05:
        flags.append(
            {
                "code": "4p1_Tdes_regression_risk",
                "detail": f"T_des 4.1 GHz {tdes0:.3e} → {tdes5:.3e} after 5 accepted iters",
            }
        )
    a2_improved = a2_5 > a2_0 + 0.5
    still_dead = (tdes5 < 1e-3) and (il5 > 30.0)
    if a2_improved and still_dead:
        flags.append(
            {
                "code": "A2_vs_engineering_mismatch_risk",
                "detail": (
                    f"A2@4.1 {a2_0:.3f}→{a2_5:.3f} but dead "
                    f"(T_des={tdes5:.3e}, IL={il5:.1f} dB, iso={iso5:.1f} dB)"
                ),
            }
        )
    # Target-notch progress (informational; no auto-abort)
    for tag, key in (("5.8125", "at_5p8125"), ("4.775", "at_4p775")):
        s0 = first.get(key) or {}
        s5 = last.get(key) or {}
        if not s0 or not s5:
            continue
        if s5.get("isolation_dB", 0) < -15 and s5.get("isolation_dB", 0) <= s0.get(
            "isolation_dB", 0
        ) + 1.0:
            flags.append(
                {
                    "code": f"notch_{tag.replace('.', 'p')}_still_poor",
                    "detail": (
                        f"{tag} GHz iso {s0.get('isolation_dB'):.1f}→{s5.get('isolation_dB'):.1f} dB; "
                        f"A2 {s0.get('J_A2'):.2f}→{s5.get('J_A2'):.2f}"
                    ),
                }
            )
    return {
        "n": len(early_rows),
        "complete": True,
        "T_des_4p1_start": float(tdes0),
        "T_des_4p1_iter5": float(tdes5),
        "IL_4p1_iter5": float(il5),
        "iso_4p1_iter5": float(iso5),
        "A2_4p1_start": float(a2_0),
        "A2_4p1_iter5": float(a2_5),
        "flags": flags,
        "abort": False,
        "note": "Flags only — do not auto-abort; do not loosen criteria.",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--maxiter", type=int, default=30)
    ap.add_argument("--tau", type=float, default=TAU_FREEZE)
    ap.add_argument("--alpha", type=float, default=ALPHA_FREEZE)
    ap.add_argument("--resume", action="store_true", default=True)
    ap.add_argument("--no-resume", action="store_false", dest="resume")
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--omp-threads-per-worker", type=int, default=1)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    assert abs(args.tau - TAU_FREEZE) < 1e-15
    assert abs(args.alpha - ALPHA_FREEZE) < 1e-15
    assert abs(BOUNDS[0] - 0.05) < 1e-15 and abs(BOUNDS[1] - 1.8) < 1e-15
    assert len(F_OBJ_GHZ) == 47

    # Parent BLAS/OpenMP = 1 to avoid oversubscription with process workers
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"

    orbits, _, _ = orbits_and_masks()
    freqs = [f * 1e9 for f in F_OBJ_GHZ]

    ckpt = load_checkpoint() if args.resume else None
    history = []
    early_rows = []
    accepted_nit = 0
    best = {"J_BB": -1e300, "q": None, "rec": None, "record": None}
    t0 = time.perf_counter()
    wall0_offset = 0.0
    tdes0_4p1 = None
    a2_4p1_start = None
    flags_written = False
    last_accepted_elapsed = None
    n_eval = 0

    if ckpt is not None:
        q0 = np.asarray(ckpt["q_current"], float)
        accepted_nit = int(ckpt.get("accepted_nit", 0))
        history = list(ckpt.get("history", []))
        early_rows = list(ckpt.get("early_diagnostic", []))
        wall0_offset = float(ckpt.get("elapsed_s", 0.0))
        best["J_BB"] = float(ckpt.get("J_BB_best", -1e300))
        best["q"] = np.asarray(ckpt.get("q_best", q0), float)
        best["rec"] = ckpt.get("rec_best")
        tdes0_4p1 = ckpt.get("T_des_4p1_start")
        a2_4p1_start = ckpt.get("A2_4p1_start")
        flags_written = bool(ckpt.get("early_flags_written", False))
        n_eval = int(ckpt.get("n_eval", 0))
        last_accepted_elapsed = wall0_offset
        print(
            f"RESUME Stage-4 from nit={accepted_nit} J_BB_best={best['J_BB']:.6f} "
            f"elapsed_so_far={wall0_offset:.0f}s",
            flush=True,
        )
    else:
        p = OUT / WARM_START
        if not p.exists():
            raise FileNotFoundError(p)
        d = json.loads(p.read_text())
        q0 = np.asarray(d["q_best"], float)
        print(f"WARM START from {p.name} q_best (Stage-3)", flush=True)

    assert len(q0) == 16, "C6 freeze requires 16 orbit variables"
    remaining = max(0, args.maxiter - accepted_nit)
    if remaining == 0:
        print("Nothing to do: accepted_nit already at maxiter", flush=True)
        return

    cache = {}
    state = {"last_rec": None, "last_q": q0.copy()}

    def elapsed():
        return wall0_offset + (time.perf_counter() - t0)

    def _eval(q):
        nonlocal n_eval
        q = np.asarray(q, float)
        key = q.tobytes()
        if key in cache:
            return cache[key]
        n_eval += 1
        t_ev0 = time.perf_counter()
        print(
            f"EVAL #{n_eval} (accepted_nit={accepted_nit}) "
            f"n_workers={args.n_workers} ...",
            flush=True,
        )
        rec = eval_bb_A2(
            args.grade, q, orbits, freqs,
            tau=args.tau, alpha=args.alpha,
            n_workers=args.n_workers,
            omp_threads_per_worker=args.omp_threads_per_worker,
        )
        print(
            f"EVAL #{n_eval} DONE wall={time.perf_counter()-t_ev0:.1f}s "
            f"J_BB={rec['J_BB']:.6f} worst@{rec['worst_f_GHz']:.3f}",
            flush=True,
        )
        cache.clear()
        cache[key] = rec
        state["last_rec"] = rec
        state["last_q"] = q.copy()
        return rec

    def _compact_rec(rec):
        out = {k: rec[k] for k in rec if k not in ("grad_q", "metrics_per_freq")}
        out["worst_metrics"] = _eng_snapshot(rec["worst_metrics"], rec["worst_J"])
        out["metrics_4p1"] = _eng_snapshot(rec["metrics_4p1"], rec["J_4p1"])
        for f, key in (
            (4.775, "metrics_4p775"),
            (5.8125, "metrics_5p8125"),
            (5.75, "metrics_5p75"),
            (5.825, "metrics_5p825"),
        ):
            snap = _snap_at_f(rec, f)
            if snap:
                out[key] = snap
        return out

    def _save_ckpt(extra=None):
        payload = {
            "spec": "BB_STAGE4_PROPOSAL.md",
            "stage": 4,
            "objective": A2_FREEZE,
            "tau": args.tau,
            "alpha": args.alpha,
            "freqs_GHz": F_OBJ_GHZ,
            "grade": args.grade,
            "accepted_nit": accepted_nit,
            "maxiter": args.maxiter,
            "elapsed_s": elapsed(),
            "n_eval": n_eval,
            "n_workers": args.n_workers,
            "omp_threads_per_worker": args.omp_threads_per_worker,
            "q_current": state["last_q"].tolist(),
            "q_best": (best["q"] if best["q"] is not None else state["last_q"]).tolist(),
            "J_BB_best": best["J_BB"],
            "rec_best": best["rec"],
            "history": history,
            "early_diagnostic": early_rows,
            "T_des_4p1_start": tdes0_4p1,
            "A2_4p1_start": a2_4p1_start,
            "early_flags_written": flags_written,
            "bounds": [BOUNDS[0], BOUNDS[1]],
            "warm_start": WARM_START,
        }
        if extra:
            payload.update(extra)
        _atomic_write(CKPT_PATH, payload)
        _atomic_write(
            HIST_PATH,
            {
                "history": history,
                "best_J_BB": best["J_BB"],
                "accepted_nit": accepted_nit,
                "elapsed_s": elapsed(),
            },
        )
        if best["q"] is not None:
            _atomic_write(
                FINAL_PATH,
                {
                    "meta": meta,
                    "status": payload.get("status", "in_progress"),
                    "elapsed_s": elapsed(),
                    "accepted_nit": accepted_nit,
                    "accepted_nit_total": accepted_nit,
                    "q0_stage4_start": q0_start.tolist(),
                    "q_best": best["q"].tolist(),
                    "J_BB_best": best["J_BB"],
                    "rec_best": best["rec"],
                    "history": history,
                    "early_diagnostic": early_rows,
                },
            )

    meta = {
        "spec": "BB_STAGE4_PROPOSAL.md",
        "stage": 4,
        "grade": args.grade,
        "design_start": WARM_START,
        "n_var": 16,
        "symmetry": "C6",
        "objective": "A2_hybrid",
        "A2_FREEZE": A2_FREEZE,
        "tau": args.tau,
        "alpha": args.alpha,
        "freqs_GHz": F_OBJ_GHZ,
        "n_freq": len(F_OBJ_GHZ),
        "enrichment_GHz": [4.775, 5.80625, 5.8125, 5.81875],
        "maxiter": args.maxiter,
        "bounds": [BOUNDS[0], BOUNDS[1]],
        "warm_start": f"{WARM_START} q_best",
        "n_workers": int(args.n_workers),
        "omp_threads_per_worker": int(args.omp_threads_per_worker),
        "authorization": "user-authorized 2026-11-10",
    }
    _atomic_write(META_PATH, meta)
    q0_start = q0.copy()

    def fun(q):
        nonlocal tdes0_4p1, a2_4p1_start
        rec = _eval(q)
        J = rec["J_BB"]
        if tdes0_4p1 is None:
            tdes0_4p1 = float(rec["metrics_4p1"]["T_des_avg"])
        if a2_4p1_start is None:
            a2_4p1_start = float(rec["J_4p1"])
        record = _record_from_rec(rec, q, accepted_nit, elapsed(), kind="eval")
        if J > best["J_BB"]:
            best["J_BB"] = J
            best["q"] = np.asarray(q, float).copy()
            best["rec"] = _compact_rec(rec)
            best["record"] = record
            s581 = record.get("at_5p8125") or {}
            s477 = record.get("at_4p775") or {}
            print(
                f"NEW BEST J_BB={J:.6f} worstA2={rec['worst_J']:.4f} "
                f"@ {rec['worst_f_GHz']:.4f} GHz "
                f"Tdes@4.1={rec['metrics_4p1']['T_des_avg']:.4e} "
                f"IL@4.1={rec['metrics_4p1']['IL_avg_dB']:.2f} "
                f"A2@5.8125={s581.get('J_A2', float('nan')):.3f} "
                f"iso@5.8125={s581.get('isolation_dB', float('nan')):.1f} "
                f"A2@4.775={s477.get('J_A2', float('nan')):.3f} "
                f"iso@4.775={s477.get('isolation_dB', float('nan')):.1f} "
                f"max_w={rec['max_weight']:.3f}",
                flush=True,
            )
            history.append({**record, "kind": "new_best"})
            _save_ckpt()
        return -J

    def jac(q):
        return -_eval(q)["grad_q"]

    def callback(xk):
        nonlocal accepted_nit, flags_written, last_accepted_elapsed
        accepted_nit += 1
        rec = _eval(xk)
        el = elapsed()
        record = _record_from_rec(
            rec, xk, accepted_nit, el, kind="accepted_iter",
            prev_elapsed=last_accepted_elapsed,
        )
        last_accepted_elapsed = el
        history.append(record)
        s581 = record.get("at_5p8125") or {}
        s477 = record.get("at_4p775") or {}
        dt = record.get("wall_s_since_prev_accepted")
        print(
            f"ACCEPTED iter={accepted_nit}/{args.maxiter} "
            f"J_BB={rec['J_BB']:.6f} worstA2={rec['worst_J']:.4f} "
            f"@ {rec['worst_f_GHz']:.4f} GHz "
            f"Tdes={record['worst'].get('T_des', float('nan')):.4e} "
            f"IL={record['worst'].get('IL_dB', float('nan')):.2f} "
            f"iso={record['worst'].get('isolation_dB', float('nan')):.2f} "
            f"A2@5.8125={s581.get('J_A2', float('nan')):.3f} "
            f"iso@5.8125={s581.get('isolation_dB', float('nan')):.1f} "
            f"A2@4.775={s477.get('J_A2', float('nan')):.3f} "
            f"iso@4.775={s477.get('isolation_dB', float('nan')):.1f} "
            f"max_w={rec['max_weight']:.3f} N_eff={rec['N_eff']:.2f} "
            f"grad={rec['grad_norm']:.3e} "
            f"dt_iter={dt if dt is not None else float('nan'):.0f}s "
            f"t={el:.0f}s",
            flush=True,
        )
        if accepted_nit <= 5:
            early_rows.append(record)
            _atomic_write(EARLY_PATH, {"rows": early_rows, "T_des_4p1_start": tdes0_4p1})
        if accepted_nit >= 5 and not flags_written:
            fl = evaluate_early_flags(
                early_rows,
                tdes0_4p1 if tdes0_4p1 is not None else 0.0,
                a2_4p1_start=a2_4p1_start,
            )
            _atomic_write(FLAGS_PATH, fl)
            flags_written = True
            if fl["flags"]:
                print("EARLY FLAGS (no auto-abort):", json.dumps(fl["flags"], indent=2), flush=True)
            else:
                print("EARLY FLAGS: none", flush=True)
        _save_ckpt()

    print("=" * 72, flush=True)
    print("START Stage-4 BB AUTHORIZED", flush=True)
    print(f"freqs_GHz ({len(F_OBJ_GHZ)}): {F_OBJ_GHZ}", flush=True)
    print(f"warm_start: {WARM_START} q_best", flush=True)
    print(
        f"workers={args.n_workers} omp_per_worker={args.omp_threads_per_worker} "
        f"parent_OMP={os.environ.get('OMP_NUM_THREADS')}",
        flush=True,
    )
    print(
        f"checkpoints: {CKPT_PATH.name}, {HIST_PATH.name}, {FINAL_PATH.name}, {META_PATH.name}",
        flush=True,
    )
    print(f"remaining L-BFGS iters={remaining} (done={accepted_nit})", flush=True)
    print("=" * 72, flush=True)

    rec0 = _eval(q0)
    if tdes0_4p1 is None:
        tdes0_4p1 = float(rec0["metrics_4p1"]["T_des_avg"])
    if a2_4p1_start is None:
        a2_4p1_start = float(rec0["J_4p1"])
    if best["q"] is None or rec0["J_BB"] > best["J_BB"]:
        best["J_BB"] = rec0["J_BB"]
        best["q"] = q0.copy()
        best["rec"] = _compact_rec(rec0)
    history.append(
        _record_from_rec(rec0, q0, accepted_nit, elapsed(), kind="start_eval")
    )
    _save_ckpt({"status": "in_progress"})
    s581 = _snap_at_f(rec0, 5.8125) or {}
    s477 = _snap_at_f(rec0, 4.775) or {}
    print(
        f"START J_BB={rec0['J_BB']:.6f} worstA2={rec0['worst_J']:.4f} "
        f"@ {rec0['worst_f_GHz']:.4f} GHz Tdes@4.1={tdes0_4p1:.4e} "
        f"IL@4.1={rec0['metrics_4p1']['IL_avg_dB']:.2f} "
        f"A2@5.8125={s581.get('J_A2', float('nan')):.3f} "
        f"iso@5.8125={s581.get('isolation_dB', float('nan')):.1f} "
        f"A2@4.775={s477.get('J_A2', float('nan')):.3f} "
        f"iso@4.775={s477.get('isolation_dB', float('nan')):.1f}",
        flush=True,
    )

    res = minimize(
        fun,
        q0,
        method="L-BFGS-B",
        jac=jac,
        bounds=[BOUNDS] * len(q0),
        callback=callback,
        options={
            "maxiter": remaining,
            "maxfun": remaining * 6,
            "ftol": 1e-12,
        },
    )
    wall = elapsed()
    out = {
        "meta": meta,
        "status": "complete",
        "elapsed_s": wall,
        "optimizer_success": bool(res.success),
        "optimizer_message": str(res.message),
        "optimizer_nit_this_launch": int(res.nit),
        "accepted_nit_total": accepted_nit,
        "n_eval": n_eval,
        "q0_stage4_start": q0_start.tolist(),
        "q_final": np.asarray(res.x, float).tolist(),
        "q_best": best["q"].tolist(),
        "J_BB_best": best["J_BB"],
        "rec_best": best["rec"],
        "history": history,
        "early_diagnostic": early_rows,
    }
    _atomic_write(FINAL_PATH, out)
    _save_ckpt({"status": "complete", "optimizer_message": str(res.message)})
    print(
        f"DONE Stage-4 J_BB={best['J_BB']:.6f} accepted_nit={accepted_nit} "
        f"elapsed={wall:.0f}s msg={res.message}",
        flush=True,
    )
    print(
        "POST: run dense Stage-4 audit + Stage-3 vs Stage-4 comparison.",
        flush=True,
    )


if __name__ == "__main__":
    main()
