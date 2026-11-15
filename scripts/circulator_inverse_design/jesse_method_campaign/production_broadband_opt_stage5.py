#!/usr/bin/env python3
"""Stage-5 4–6 GHz broadband optimization (AUTHORIZED).

Frozen scientific package (unchanged from Stages 2–4):
  C6 / 16-var, A2_hybrid, soft-min τ=0.05 α=0.25, bounds [0.05, 1.8],
  FULL91-C, L-BFGS-B maxiter=30,
  deterministic full-band 6.25 MHz grid (n_freq=321),
  warm-start from opt_broadband_stage4_A2_C.json q_best,
  process-parallel frequencies (n_workers=12, omp=1),
  optimized A2 kernel (PLASMEEP_A2_FAST_KERNEL=1 default).

Does not change FEM physics, A2, mesh, materials, derivatives, ports, or B sign.
Does not modify the 321-frequency sampling rule during optimization.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
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
from stage5_frequency_set import stage5_freqs_GHz  # noqa: E402

F_OBJ_GHZ = stage5_freqs_GHz()
assert len(F_OBJ_GHZ) == 321
assert F_OBJ_GHZ[0] == 4.0 and F_OBJ_GHZ[-1] == 6.0
assert abs(max(F_OBJ_GHZ[i + 1] - F_OBJ_GHZ[i] for i in range(320)) - 0.00625) < 1e-12
assert len(set(F_OBJ_GHZ)) == 321

TAU_FREEZE = 0.05
ALPHA_FREEZE = 0.25
BOUNDS = (S_LO, S_HI)  # [0.05, 1.8]
WARM_START = "opt_broadband_stage4_A2_C.json"

# Prior failure / watch frequencies (all on the 6.25 MHz grid)
WATCH_GHZ = [4.1, 4.73125, 4.775, 5.08125, 5.5875, 5.75, 5.8125, 5.825]

CKPT_PATH = OUT / "bb_stage5_checkpoint.json"
HIST_PATH = OUT / "bb_stage5_history.json"
META_PATH = OUT / "bb_stage5_meta.json"
FINAL_PATH = OUT / "opt_broadband_stage5_A2_C.json"
EARLY_PATH = OUT / "bb_stage5_early_diagnostic.json"
FLAGS_PATH = OUT / "bb_stage5_early_flags.json"


def eval_bb_A2(grade, q, orbits, freqs_hz, tau=TAU_FREEZE, alpha=ALPHA_FREEZE,
               n_workers=1, omp_threads_per_worker=1):
    return eval_bb_A2_core(
        grade, q, orbits, freqs_hz,
        tau=tau, alpha=alpha,
        n_workers=n_workers,
        omp_threads_per_worker=omp_threads_per_worker,
        profile=False,
        use_fast_kernel=True,
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


def _watch_snaps(rec):
    out = {}
    for f in WATCH_GHZ:
        key = f"at_{str(f).replace('.', 'p')}"
        out[key] = _snap_at_f(rec, f)
    return out


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
        "J_per_freq": rec["J_per_freq"],
        "weights": rec["weights"],
        "freqs_GHz": rec["freqs_GHz"],
    }
    out.update(_watch_snaps(rec))
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


def git_head() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def print_prelaunch(args, q0):
    gaps = [F_OBJ_GHZ[i + 1] - F_OBJ_GHZ[i] for i in range(len(F_OBJ_GHZ) - 1)]
    print("=" * 72, flush=True)
    print("STAGE-5 PRE-LAUNCH VERIFICATION", flush=True)
    print(f"1. branch/commit: agent/eigenmode-ports @ {git_head()}", flush=True)
    print(f"2. warm-start: {WARM_START} q_best  (n_q={len(q0)})", flush=True)
    print(
        "3. frequency rule: [round(4.0 + i*0.00625, 6) for i in range(321)]",
        flush=True,
    )
    print(f"4. n_freq = {len(F_OBJ_GHZ)}", flush=True)
    print(f"5. first: {F_OBJ_GHZ[:6]}  last: {F_OBJ_GHZ[-6:]}", flush=True)
    print(f"6. max spacing = {max(gaps)*1000:.6f} MHz", flush=True)
    print(
        f"7. A2_hybrid / tau={args.tau} / alpha={args.alpha} / "
        f"bounds=[{BOUNDS[0]}, {BOUNDS[1]}]",
        flush=True,
    )
    print(
        f"8. fast kernel: PLASMEEP_A2_FAST_KERNEL="
        f"{os.environ.get('PLASMEEP_A2_FAST_KERNEL', '1')} "
        f"(eval uses use_fast_kernel=True)",
        flush=True,
    )
    print(f"9. n_workers = {args.n_workers}", flush=True)
    print(
        f"10. OMP={os.environ.get('OMP_NUM_THREADS')} "
        f"OPENBLAS={os.environ.get('OPENBLAS_NUM_THREADS')} "
        f"MKL={os.environ.get('MKL_NUM_THREADS')}",
        flush=True,
    )
    print(f"11. checkpoint: {CKPT_PATH}", flush=True)
    print(f"12. outputs: {OUT}  final={FINAL_PATH.name}", flush=True)
    assert len(F_OBJ_GHZ) == 321
    assert len(set(F_OBJ_GHZ)) == 321
    assert F_OBJ_GHZ[0] == 4.0 and F_OBJ_GHZ[-1] == 6.0
    assert abs(min(gaps) - 0.00625) < 1e-12 and abs(max(gaps) - 0.00625) < 1e-12
    assert args.n_workers == 12
    assert int(args.omp_threads_per_worker) == 1
    assert abs(args.tau - TAU_FREEZE) < 1e-15
    assert abs(args.alpha - ALPHA_FREEZE) < 1e-15
    print("SANITY: frequency list OK (321 unique, 4.0→6.0, Δf=6.25 MHz)", flush=True)
    print("=" * 72, flush=True)


def evaluate_early_flags(early_rows, tdes0, a2_4p1_start=None):
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
    for tag, key in (
        ("5.08125", "at_5p08125"),
        ("4.73125", "at_4p73125"),
        ("5.8125", "at_5p8125"),
    ):
        s0 = first.get(key) or {}
        s5 = last.get(key) or {}
        if not s0 or not s5:
            continue
        if s5.get("J_A2", 0) < -5 and s5.get("J_A2", 0) <= s0.get("J_A2", 0) + 0.5:
            flags.append(
                {
                    "code": f"watch_{tag.replace('.', 'p')}_still_poor",
                    "detail": (
                        f"{tag} GHz A2 {s0.get('J_A2'):.2f}→{s5.get('J_A2'):.2f}; "
                        f"iso {s0.get('isolation_dB'):.1f}→{s5.get('isolation_dB'):.1f} dB"
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

    os.environ.setdefault("PLASMEEP_A2_FAST_KERNEL", "1")
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
            f"RESUME Stage-5 from nit={accepted_nit} J_BB_best={best['J_BB']:.6f} "
            f"elapsed_so_far={wall0_offset:.0f}s n_eval={n_eval}",
            flush=True,
        )
    else:
        p = OUT / WARM_START
        if not p.exists():
            raise FileNotFoundError(p)
        d = json.loads(p.read_text())
        q0 = np.asarray(d["q_best"], float)
        print(f"WARM START from {p.name} q_best (Stage-4 FINAL)", flush=True)

    assert len(q0) == 16, "C6 freeze requires 16 orbit variables"
    print_prelaunch(args, q0)

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
            f"n_workers={args.n_workers} n_freq={len(F_OBJ_GHZ)} ...",
            flush=True,
        )
        rec = eval_bb_A2(
            args.grade, q, orbits, freqs,
            tau=args.tau, alpha=args.alpha,
            n_workers=args.n_workers,
            omp_threads_per_worker=args.omp_threads_per_worker,
        )
        dt = time.perf_counter() - t_ev0
        print(
            f"EVAL #{n_eval} DONE wall={dt:.1f}s "
            f"J_BB={rec['J_BB']:.6f} worstA2={rec['worst_J']:.4f} "
            f"@ {rec['worst_f_GHz']:.5f} GHz "
            f"IL={rec['worst_metrics']['IL_avg_dB']:.2f} "
            f"iso={rec['worst_metrics']['isolation_avg_dB']:.2f}",
            flush=True,
        )
        cache.clear()
        cache[key] = rec
        state["last_rec"] = rec
        state["last_q"] = q.copy()
        return rec

    def _compact_rec(rec):
        # Keep per-freq metrics for post-audit seeding (321 freqs; no grad_q).
        out = {k: rec[k] for k in rec if k != "grad_q"}
        out["worst_metrics"] = _eng_snapshot(rec["worst_metrics"], rec["worst_J"])
        out["metrics_4p1"] = _eng_snapshot(rec["metrics_4p1"], rec["J_4p1"])
        for f in WATCH_GHZ:
            snap = _snap_at_f(rec, f)
            if snap:
                out[f"metrics_{str(f).replace('.', 'p')}"] = snap
        return out

    def _save_ckpt(extra=None):
        payload = {
            "spec": "BB_STAGE5_SAMPLING_STRATEGY.md",
            "stage": 5,
            "objective": A2_FREEZE,
            "tau": args.tau,
            "alpha": args.alpha,
            "freqs_GHz": F_OBJ_GHZ,
            "n_freq": len(F_OBJ_GHZ),
            "grade": args.grade,
            "accepted_nit": accepted_nit,
            "maxiter": args.maxiter,
            "elapsed_s": elapsed(),
            "n_eval": n_eval,
            "n_workers": args.n_workers,
            "omp_threads_per_worker": args.omp_threads_per_worker,
            "fast_kernel": True,
            "git_commit": git_head(),
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
                "n_eval": n_eval,
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
                    "n_eval": n_eval,
                    "q0_stage5_start": q0_start.tolist(),
                    "q_best": best["q"].tolist(),
                    "J_BB_best": best["J_BB"],
                    "rec_best": best["rec"],
                    "history": history,
                    "early_diagnostic": early_rows,
                },
            )

    meta = {
        "spec": "BB_STAGE5_SAMPLING_STRATEGY.md",
        "stage": 5,
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
        "df_MHz": 6.25,
        "maxiter": args.maxiter,
        "bounds": [BOUNDS[0], BOUNDS[1]],
        "warm_start": f"{WARM_START} q_best",
        "n_workers": int(args.n_workers),
        "omp_threads_per_worker": int(args.omp_threads_per_worker),
        "fast_kernel": True,
        "git_commit": git_head(),
        "authorization": "user-authorized 2026-11-15",
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
            s508 = record.get("at_5p08125") or {}
            s473 = record.get("at_4p73125") or {}
            print(
                f"NEW BEST J_BB={J:.6f} worstA2={rec['worst_J']:.4f} "
                f"@ {rec['worst_f_GHz']:.5f} GHz "
                f"IL={record['worst'].get('IL_dB', float('nan')):.2f} "
                f"iso={record['worst'].get('isolation_dB', float('nan')):.2f} "
                f"A2@5.08125={s508.get('J_A2', float('nan')):.3f} "
                f"A2@4.73125={s473.get('J_A2', float('nan')):.3f} "
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
        s508 = record.get("at_5p08125") or {}
        s473 = record.get("at_4p73125") or {}
        s581 = record.get("at_5p8125") or {}
        dt = record.get("wall_s_since_prev_accepted")
        print(
            f"ACCEPTED iter={accepted_nit}/{args.maxiter} "
            f"J_BB={rec['J_BB']:.6f} worstA2={rec['worst_J']:.4f} "
            f"@ {rec['worst_f_GHz']:.5f} GHz "
            f"Tdes={record['worst'].get('T_des', float('nan')):.4e} "
            f"IL={record['worst'].get('IL_dB', float('nan')):.2f} "
            f"iso={record['worst'].get('isolation_dB', float('nan')):.2f} "
            f"A2@5.08125={s508.get('J_A2', float('nan')):.3f} "
            f"A2@4.73125={s473.get('J_A2', float('nan')):.3f} "
            f"A2@5.8125={s581.get('J_A2', float('nan')):.3f} "
            f"max_w={rec['max_weight']:.3f} N_eff={rec['N_eff']:.2f} "
            f"grad={rec['grad_norm']:.3e} "
            f"dt_iter={dt if dt is not None else float('nan'):.0f}s "
            f"n_eval={n_eval} t={el:.0f}s",
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
    print("START Stage-5 BB AUTHORIZED", flush=True)
    print(f"n_freq={len(F_OBJ_GHZ)} df=6.25 MHz  warm_start={WARM_START}", flush=True)
    print(
        f"workers={args.n_workers} omp_per_worker={args.omp_threads_per_worker} "
        f"fast_kernel=True remaining_iters={remaining}",
        flush=True,
    )
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
    print(
        f"START J_BB={rec0['J_BB']:.6f} worstA2={rec0['worst_J']:.4f} "
        f"@ {rec0['worst_f_GHz']:.5f} GHz "
        f"IL={rec0['worst_metrics']['IL_avg_dB']:.2f} "
        f"iso={rec0['worst_metrics']['isolation_avg_dB']:.2f}",
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
        "q0_stage5_start": q0_start.tolist(),
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
        f"DONE Stage-5 J_BB={best['J_BB']:.6f} accepted_nit={accepted_nit} "
        f"n_eval={n_eval} elapsed={wall:.0f}s msg={res.message}",
        flush=True,
    )
    print("POST: run Stage-5 dense audit + 3.125 MHz refinement.", flush=True)


if __name__ == "__main__":
    main()
