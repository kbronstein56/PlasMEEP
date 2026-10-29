#!/usr/bin/env python3
"""Stage-2 4–6 GHz broadband optimization per BROADBAND_PRODUCTION_SPEC_STAGE2_FROZEN.

Frozen:
  C6, Objective A2_hybrid, soft-min τ=0.05 α=0.25,
  28-frequency enriched grid, FULL91-C, L-BFGS-B,
  warm-start from opt_broadband_production_C.json.

Checkpointing:
  Writes after every accepted L-BFGS iteration and on every new best.
  Restartable from latest checkpoint without recomputing completed iters.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.sparse.linalg import splu

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

from adjoint_direct_dtn import build_system  # noqa: E402
from broadband_objective import softmin_bb  # noqa: E402
from diff_fem_core import assemble_dA_rod  # noqa: E402
from direct_dtn_s import build_Q_columns  # noqa: E402
from efficient_adjoint import S_from_U, _pack_C, _solve_multi  # noqa: E402
from multimode_opt_core import (  # noqa: E402
    B_PLUS,
    S_HI,
    S_LO,
    expand_tied,
    orbits_and_masks,
    te0_metrics_from_S,
)
from objective_A2 import (  # noqa: E402
    A2_FREEZE,
    EPS,
    objective_A2_from_S,
    scalar_adjoint_rhs_A2,
)

# Frozen Stage-2 frequency grid (GHz)
F_OBJ_GHZ = [
    4.0,
    4.1,
    4.2,
    4.3,
    4.4,
    4.5,
    4.55,
    4.6,
    4.7,
    4.8,
    4.85,
    4.9,
    5.0,
    5.1,
    5.2,
    5.3,
    5.35,
    5.4,
    5.5,
    5.6,
    5.625,
    5.65,
    5.7,
    5.8,
    5.85,
    5.9,
    5.95,
    6.0,
]
assert len(F_OBJ_GHZ) == 28

TAU_FREEZE = 0.05
ALPHA_FREEZE = 0.25
BOUNDS = (S_LO, S_HI)  # [0.05, 1.8]

CKPT_PATH = OUT / "bb_stage2_checkpoint.json"
HIST_PATH = OUT / "bb_stage2_history.json"
META_PATH = OUT / "bb_stage2_meta.json"
FINAL_PATH = OUT / "opt_broadband_stage2_A2_C.json"
EARLY_PATH = OUT / "bb_stage2_early_diagnostic.json"
FLAGS_PATH = OUT / "bb_stage2_early_flags.json"


def eval_and_grad_A2(grade, q, orbits, f_hz, b_tesla=B_PLUS, eps=EPS):
    """Single-frequency A2 objective + C6 gradient (efficient TE0 adjoint)."""
    s = expand_tied(np.asarray(q, float), orbits)
    st = build_system(grade, s, f_hz, b_tesla, "projected")
    _, _, o_masks = orbits_and_masks(rod_masks=st["rod_masks"])
    points, A = st["points"], st["A"]
    cols = build_Q_columns(points, st["ports"], st["modes"], n_prop=2)
    C, betas = _pack_C(points, cols)
    te0_cols = [2 * j for j in range(6)]
    lu = splu(A.tocsc())
    U = _solve_multi(lu, C[:, te0_cols], trans="N")
    S = S_from_U(C, U, betas, te0_cols)
    J = objective_A2_from_S(S, eps=eps)
    met = te0_metrics_from_S(S)
    R = np.column_stack(
        [scalar_adjoint_rhs_A2(C, betas, S, j, eps=eps) for j in range(6)]
    )
    Lam = _solve_multi(lu, R, trans="T")
    g = np.zeros(len(orbits), float)
    for m, o in enumerate(orbits):
        s_m = float(s[o["rods"][0]])
        dA = assemble_dA_rod(
            points,
            st["tris"],
            o_masks[m],
            st["f_ord"],
            st["fp_ref_ord"],
            st["gamma"],
            st["fc"],
            st["k0"],
            s=s_m,
        )
        acc = 0.0
        for j_loc in range(6):
            dAu = dA @ U[:, j_loc]
            acc += 2.0 * np.real(np.dot(Lam[:, j_loc], dAu))
        g[m] = float(acc / 6.0)
    return {"J": J, "S": S, "metrics": met}, g, s


def eval_bb_A2(grade, q, orbits, freqs_hz, tau=TAU_FREEZE, alpha=ALPHA_FREEZE):
    Js, gs, mets = [], [], []
    for f in freqs_hz:
        bun, g, _ = eval_and_grad_A2(grade, q, orbits, f_hz=float(f))
        Js.append(bun["J"])
        gs.append(g)
        mets.append(bun["metrics"])
    Js = np.asarray(Js, float)
    J_bb, w, J_soft, J_mean, p = softmin_bb(Js, tau=tau, alpha=alpha)
    g_bb = sum(wi * gi for wi, gi in zip(w, gs))
    worst_i = int(np.argmin(Js))
    idx_41 = int(np.argmin(np.abs(np.asarray(freqs_hz) / 1e9 - 4.1)))
    return {
        "J_BB": float(J_bb),
        "J_softmin": float(J_soft),
        "J_mean": float(J_mean),
        "J_per_freq": Js.tolist(),
        "weights": np.asarray(w, float).tolist(),
        "softmax_p": np.asarray(p, float).tolist(),
        "freqs_Hz": [float(f) for f in freqs_hz],
        "freqs_GHz": [float(f) / 1e9 for f in freqs_hz],
        "metrics_per_freq": mets,
        "grad_q": np.asarray(g_bb, float),
        "grad_norm": float(np.linalg.norm(g_bb)),
        "tau": tau,
        "alpha": alpha,
        "worst_J": float(Js[worst_i]),
        "worst_f_Hz": float(freqs_hz[worst_i]),
        "worst_f_GHz": float(freqs_hz[worst_i] / 1e9),
        "worst_metrics": mets[worst_i],
        "metrics_4p1": mets[idx_41],
        "J_4p1": float(Js[idx_41]),
        "max_weight": float(np.max(w)),
        "N_eff": float(1.0 / np.sum(np.asarray(w, float) ** 2)),
        "top5_weights": _top_weights(freqs_hz, w, Js, 5),
    }


def _top_weights(freqs_hz, w, Js, k=5):
    w = np.asarray(w, float)
    order = np.argsort(-w)[:k]
    return [
        {
            "f_GHz": float(freqs_hz[i] / 1e9),
            "weight": float(w[i]),
            "J_A2": float(Js[i]),
        }
        for i in order
    ]


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


def _record_from_rec(rec, q, iteration, elapsed_s, kind):
    return {
        "kind": kind,
        "iteration": int(iteration),
        "elapsed_s": float(elapsed_s),
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
    """Flag (do not abort) after 5 accepted iterations."""
    flags = []
    if len(early_rows) < 5:
        return {"n": len(early_rows), "flags": flags, "complete": False}
    first = early_rows[0]
    last = early_rows[4]  # 5th accepted
    tdes5 = last["at_4p1"]["T_des"]
    il5 = last["at_4p1"]["IL_dB"]
    iso5 = last["at_4p1"]["isolation_dB"]
    a2_5 = last["at_4p1"]["J_A2"]
    a2_0 = float(a2_4p1_start) if a2_4p1_start is not None else float(first["at_4p1"]["J_A2"])
    # Material T_des improvement: ≥2× or absolute gain ≥ 5e-3
    improved = (tdes5 >= 2.0 * max(tdes0, 1e-12)) or ((tdes5 - tdes0) >= 5e-3)
    if not improved:
        flags.append(
            {
                "code": "4p1_Tdes_stalled",
                "detail": f"T_des 4.1 GHz {tdes0:.3e} → {tdes5:.3e} after 5 accepted iters",
            }
        )
    # A2 improved while engineering at 4.1 remains essentially dead
    a2_improved = a2_5 > a2_0 + 0.5
    still_dead = (tdes5 < 1e-3) and (il5 > 30.0)
    if a2_improved and still_dead:
        flags.append(
            {
                "code": "A2_vs_engineering_mismatch_risk",
                "detail": (
                    f"A2@4.1 {a2_0:.3f}→{a2_5:.3f} but still dead "
                    f"(T_des={tdes5:.3e}, IL={il5:.1f} dB, iso={iso5:.1f} dB)"
                ),
            }
        )
    # Oscillation: worst freq flips without worst A2 net improvement
    worst_fs = [r["worst_f_GHz"] for r in early_rows]
    worst_js = [r["worst_A2"] for r in early_rows]
    n_unique = len(set(round(f, 6) for f in worst_fs))
    if n_unique >= 3 and worst_js[-1] <= worst_js[0] + 0.05:
        flags.append(
            {
                "code": "worst_freq_oscillation",
                "detail": (
                    f"worst freq sequence={worst_fs}, "
                    f"worst_A2 {worst_js[0]:.3f}→{worst_js[-1]:.3f}"
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
        "note": "Flags only — do not auto-abort; report before changing tau/alpha/objective.",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--design", default="opt_broadband_production_C")
    ap.add_argument("--maxiter", type=int, default=30)
    ap.add_argument("--tau", type=float, default=TAU_FREEZE)
    ap.add_argument("--alpha", type=float, default=ALPHA_FREEZE)
    ap.add_argument("--resume", action="store_true", default=True)
    ap.add_argument("--no-resume", action="store_false", dest="resume")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    assert abs(args.tau - TAU_FREEZE) < 1e-15
    assert abs(args.alpha - ALPHA_FREEZE) < 1e-15
    assert abs(BOUNDS[0] - 0.05) < 1e-15 and abs(BOUNDS[1] - 1.8) < 1e-15

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
        print(
            f"RESUME Stage-2 from nit={accepted_nit} J_BB_best={best['J_BB']:.6f} "
            f"elapsed_so_far={wall0_offset:.0f}s",
            flush=True,
        )
    else:
        # warm start from Stage-1 production best
        p = OUT / f"{args.design}.json"
        if not p.exists():
            raise FileNotFoundError(p)
        d = json.loads(p.read_text())
        q0 = np.asarray(d["q_best"], float)
        print(f"WARM START from {p.name} q_best", flush=True)

    assert len(q0) == 16, "C6 freeze requires 16 orbit variables"
    remaining = max(0, args.maxiter - accepted_nit)
    if remaining == 0:
        print("Nothing to do: accepted_nit already at maxiter", flush=True)
        return

    cache = {}
    n_eval = 0
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
        print(f"EVAL #{n_eval} (accepted_nit={accepted_nit}) ...", flush=True)
        rec = eval_bb_A2(
            args.grade, q, orbits, freqs, tau=args.tau, alpha=args.alpha
        )
        cache.clear()
        cache[key] = rec
        state["last_rec"] = rec
        state["last_q"] = q.copy()
        return rec

    def _save_ckpt(extra=None):
        payload = {
            "spec": "BROADBAND_PRODUCTION_SPEC_STAGE2_FROZEN.md",
            "objective": A2_FREEZE,
            "tau": args.tau,
            "alpha": args.alpha,
            "freqs_GHz": F_OBJ_GHZ,
            "grade": args.grade,
            "accepted_nit": accepted_nit,
            "maxiter": args.maxiter,
            "elapsed_s": elapsed(),
            "n_eval": n_eval,
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
                    "status": "in_progress",
                    "elapsed_s": elapsed(),
                    "accepted_nit": accepted_nit,
                    "q0_stage2_start": q0_start.tolist(),
                    "q_best": best["q"].tolist(),
                    "J_BB_best": best["J_BB"],
                    "rec_best": best["rec"],
                    "history": history,
                    "early_diagnostic": early_rows,
                },
            )

    meta = {
        "spec": "BROADBAND_PRODUCTION_SPEC_STAGE2_FROZEN.md",
        "grade": args.grade,
        "design_start": args.design,
        "n_var": 16,
        "symmetry": "C6",
        "objective": "A2_hybrid",
        "A2_FREEZE": A2_FREEZE,
        "tau": args.tau,
        "alpha": args.alpha,
        "freqs_GHz": F_OBJ_GHZ,
        "n_freq": len(F_OBJ_GHZ),
        "maxiter": args.maxiter,
        "bounds": [BOUNDS[0], BOUNDS[1]],
        "warm_start": "opt_broadband_production_C.json q_best",
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
            # store rec without huge unused keys duplication
            best["rec"] = {
                k: rec[k]
                for k in rec
                if k
                not in (
                    "grad_q",
                    "metrics_per_freq",
                )
            }
            # keep compact metrics for worst + 4.1
            best["rec"]["worst_metrics"] = _eng_snapshot(
                rec["worst_metrics"], rec["worst_J"]
            )
            best["rec"]["metrics_4p1"] = _eng_snapshot(
                rec["metrics_4p1"], rec["J_4p1"]
            )
            best["record"] = record
            print(
                f"NEW BEST J_BB={J:.6f} worstA2={rec['worst_J']:.4f} "
                f"@ {rec['worst_f_GHz']:.3f} GHz meanA2={rec['J_mean']:.4f} "
                f"Tdes@4.1={rec['metrics_4p1']['T_des_avg']:.4e} "
                f"IL@4.1={rec['metrics_4p1']['IL_avg_dB']:.2f} "
                f"max_w={rec['max_weight']:.3f}",
                flush=True,
            )
            history.append({**record, "kind": "new_best"})
            _save_ckpt()
        return -J

    def jac(q):
        return -_eval(q)["grad_q"]

    def callback(xk):
        nonlocal accepted_nit, flags_written
        accepted_nit += 1
        rec = _eval(xk)
        record = _record_from_rec(
            rec, xk, accepted_nit, elapsed(), kind="accepted_iter"
        )
        history.append(record)
        print(
            f"ACCEPTED iter={accepted_nit}/{args.maxiter} "
            f"J_BB={rec['J_BB']:.6f} worstA2={rec['worst_J']:.4f} "
            f"@ {rec['worst_f_GHz']:.3f} GHz "
            f"Tdes@4.1={rec['metrics_4p1']['T_des_avg']:.4e} "
            f"IL@4.1={rec['metrics_4p1']['IL_avg_dB']:.2f} "
            f"iso@4.1={rec['metrics_4p1']['isolation_avg_dB']:.2f} "
            f"max_w={rec['max_weight']:.3f} N_eff={rec['N_eff']:.2f} "
            f"grad={rec['grad_norm']:.3e} t={elapsed():.0f}s",
            flush=True,
        )
        if accepted_nit <= 5:
            early_rows.append(record)
            _atomic_write(EARLY_PATH, {"rows": early_rows, "T_des_4p1_start": tdes0_4p1})
            print(
                f"EARLY DIAG iter={accepted_nit}: 4.1 "
                f"A2={record['at_4p1']['J_A2']:.4f} "
                f"Tdes={record['at_4p1']['T_des']:.4e} "
                f"IL={record['at_4p1']['IL_dB']:.2f} "
                f"iso={record['at_4p1']['isolation_dB']:.2f} "
                f"Trev={record['at_4p1']['T_rev']:.4e} "
                f"R={record['at_4p1']['R']:.4f} TE1={record['at_4p1']['TE1']:.4f} "
                f"| top5={record['top5_weights']}",
                flush=True,
            )
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
                print("EARLY FLAGS: none — 4.1 progress looks healthy enough to continue", flush=True)
        _save_ckpt()

    print("START Stage-2 BB", meta, flush=True)
    print(f"remaining L-BFGS iters={remaining} (done={accepted_nit})", flush=True)

    # Initial eval / baseline checkpoint before optimize
    rec0 = _eval(q0)
    if tdes0_4p1 is None:
        tdes0_4p1 = float(rec0["metrics_4p1"]["T_des_avg"])
    if a2_4p1_start is None:
        a2_4p1_start = float(rec0["J_4p1"])
    if best["q"] is None or rec0["J_BB"] > best["J_BB"]:
        best["J_BB"] = rec0["J_BB"]
        best["q"] = q0.copy()
        best["rec"] = {
            k: rec0[k]
            for k in rec0
            if k not in ("grad_q", "metrics_per_freq")
        }
        best["rec"]["worst_metrics"] = _eng_snapshot(
            rec0["worst_metrics"], rec0["worst_J"]
        )
        best["rec"]["metrics_4p1"] = _eng_snapshot(rec0["metrics_4p1"], rec0["J_4p1"])
    history.append(
        _record_from_rec(rec0, q0, accepted_nit, elapsed(), kind="start_eval")
    )
    _save_ckpt()
    print(
        f"START J_BB={rec0['J_BB']:.6f} worstA2={rec0['worst_J']:.4f} "
        f"@ {rec0['worst_f_GHz']:.3f} GHz Tdes@4.1={tdes0_4p1:.4e} "
        f"IL@4.1={rec0['metrics_4p1']['IL_avg_dB']:.2f} "
        f"max_w={rec0['max_weight']:.3f} N_eff={rec0['N_eff']:.2f}",
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
        "q0_stage2_start": q0_start.tolist(),
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
        f"DONE Stage-2 J_BB={best['J_BB']:.6f} accepted_nit={accepted_nit} "
        f"elapsed={wall:.0f}s msg={res.message}",
        flush=True,
    )
    print(
        "POST: do NOT auto-start Stage-3. Dense 50/25 MHz audit + Stage-1 vs Stage-2 "
        "comparison required first.",
        flush=True,
    )


if __name__ == "__main__":
    main()
