#!/usr/bin/env python3
"""C6 L-BFGS-B multimode direct-DtN optimization (Phases 3–4).

Uses only the frozen objective/weights/bounds from MULTIMODE_SINGLE_FREQ_OBJECTIVE_FROZEN.md.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from multimode_opt_core import (
    B_PLUS,
    F_HZ,
    OUT,
    S_HI,
    S_LO,
    WEIGHTS,
    dump_json,
    eval_and_grad,
    expand_tied,
    load_start,
    orbits_and_masks,
)


def history_row(it, J, gnorm, q, met, status):
    return {
        "iteration": it,
        "J": float(J),
        "grad_norm": float(gnorm),
        "q": np.asarray(q, float).tolist(),
        "q_min": float(np.min(q)),
        "q_max": float(np.max(q)),
        "T_des_avg": met["T_des_avg"],
        "T_des_min": met["T_des_min"],
        "isolation_avg_dB": met["isolation_avg_dB"],
        "isolation_min_dB": met["isolation_min_dB"],
        "IL_avg_dB": met["IL_avg_dB"],
        "IL_worst_dB": met["IL_worst_dB"],
        "RL_avg_dB": met["RL_avg_dB"],
        "RL_min_dB": met["RL_min_dB"],
        "T_thr_avg": met["T_thr_avg"],
        "T_other_avg": met["T_other_avg"],
        "TE1_avg": met["TE1_avg"],
        "R_avg": met["R_avg"],
        "P_material_avg": met.get("P_material_avg"),
        "P_out_avg": met.get("P_out_avg"),
        "J_std": met.get("J_std"),
        "status": status,
    }


def run_opt(
    grade: str,
    start_label: str,
    start_tag: str,
    maxiter: int,
    b_tesla: float = B_PLUS,
    f_hz: float = F_HZ,
    with_absorption: bool = True,
    tag_suffix: str = "",
):
    OUT.mkdir(parents=True, exist_ok=True)
    orbits, _, _ = orbits_and_masks()
    q0, s0, tag = load_start(start_tag, orbits)
    run_id = f"{start_label}_{grade}{tag_suffix}"
    print(
        f"START run={run_id} orbits={len(orbits)} maxiter={maxiter} "
        f"q0 range=[{q0.min():.4f},{q0.max():.4f}]",
        flush=True,
    )

    history = []
    best = {"J": -1e300, "q": None, "s": None, "metrics": None, "iteration": -1}
    cache = {"q": None, "J": None, "g": None, "met": None, "s": None}
    n_fev = [0]
    n_gev = [0]
    t0 = time.time()

    def _eval(q):
        q = np.asarray(q, float)
        if cache["q"] is not None and np.allclose(q, cache["q"], rtol=0, atol=0):
            return cache["J"], cache["g"], cache["met"], cache["s"]
        bun, gq, s = eval_and_grad(
            grade, q, orbits, f_hz, b_tesla, with_absorption=with_absorption
        )
        cache["q"] = q.copy()
        cache["J"] = float(bun["J"])
        cache["g"] = np.asarray(gq, float)
        cache["met"] = bun["metrics"]
        cache["s"] = s
        return cache["J"], cache["g"], cache["met"], cache["s"]

    def fun(q):
        n_fev[0] += 1
        J, g, met, s = _eval(q)
        gnorm = float(np.linalg.norm(g))
        # maximize J via minimize -J
        if J > best["J"]:
            best["J"] = J
            best["q"] = np.asarray(q, float).copy()
            best["s"] = s.copy()
            best["metrics"] = met
            best["iteration"] = len(history)
            print(
                f"  NEW BEST it≈{len(history)} J={J:.6f} Tdes={met['T_des_avg']:.4f} "
                f"TE1={met['TE1_avg']:.4f} iso={met['isolation_avg_dB']:.2f} "
                f"|g|={gnorm:.3e}",
                flush=True,
            )
        row = history_row(len(history), J, gnorm, q, met, "feval")
        history.append(row)
        # checkpoint every eval
        dump_json(OUT / f"opt_history_{run_id}.json", {"history": history, "best_J": best["J"]})
        return -J

    def jac(q):
        n_gev[0] += 1
        J, g, met, s = _eval(q)
        return -g

    bounds = [(S_LO, S_HI)] * len(orbits)
    # baseline at start
    J0, g0, met0, s_init = _eval(q0)
    history.append(history_row(0, J0, float(np.linalg.norm(g0)), q0, met0, "start"))
    best["J"] = J0
    best["q"] = q0.copy()
    best["s"] = s_init.copy()
    best["metrics"] = met0
    best["iteration"] = 0
    print(
        f"  baseline J={J0:.6f} Tdes={met0['T_des_avg']:.4f} TE1={met0['TE1_avg']:.4f}",
        flush=True,
    )

    res = minimize(
        fun,
        q0,
        method="L-BFGS-B",
        jac=jac,
        bounds=bounds,
        options={"maxiter": maxiter, "ftol": 1e-9, "gtol": 1e-6, "disp": True},
    )
    elapsed = time.time() - t0
    q_final = np.asarray(res.x, float)
    Jf, gf, metf, sf = _eval(q_final)
    final_row = history_row(len(history), Jf, float(np.linalg.norm(gf)), q_final, metf, str(res.message))
    history.append(final_row)

    out = {
        "run_id": run_id,
        "start_label": start_label,
        "start_tag": tag,
        "grade": grade,
        "B": b_tesla,
        "f_Hz": f_hz,
        "weights": WEIGHTS,
        "bounds": [S_LO, S_HI],
        "maxiter": maxiter,
        "n_fev": n_fev[0],
        "n_gev": n_gev[0],
        "elapsed_s": elapsed,
        "optimizer_success": bool(res.success),
        "optimizer_message": str(res.message),
        "optimizer_nit": int(res.nit),
        "q0": q0.tolist(),
        "J0": float(J0),
        "metrics0": met0,
        "best": {
            "J": float(best["J"]),
            "q": best["q"].tolist(),
            "s": best["s"].tolist(),
            "metrics": best["metrics"],
            "iteration": best["iteration"],
            "q_min": float(np.min(best["q"])),
            "q_max": float(np.max(best["q"])),
        },
        "final": {
            "J": float(Jf),
            "q": q_final.tolist(),
            "s": sf.tolist(),
            "metrics": metf,
            "q_min": float(np.min(q_final)),
            "q_max": float(np.max(q_final)),
        },
        "history": history,
        "S_at_best_saved": False,
    }
    # save S at best
    bun_best, _, _ = eval_and_grad(
        grade, best["q"], orbits, f_hz, b_tesla, with_absorption=True
    )
    np.savez_compressed(
        OUT / f"S12_opt_{run_id}_Bp.npz",
        S=bun_best["S"],
        q=best["q"],
        s=best["s"],
    )
    out["best"]["metrics"] = bun_best["metrics"]
    out["best"]["J"] = float(bun_best["J"])
    out["S_at_best_saved"] = True
    dump_json(OUT / f"opt_{run_id}.json", out)
    print(
        f"DONE {run_id}: best J={out['best']['J']:.6f} "
        f"Tdes={out['best']['metrics']['T_des_avg']:.4f} "
        f"TE1={out['best']['metrics']['TE1_avg']:.4f} "
        f"elapsed={elapsed/3600:.2f} h",
        flush=True,
    )
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument(
        "--start",
        default="pref_peak",
        choices=("uniform", "pref_peak", "pref_broadband", "custom"),
    )
    ap.add_argument("--maxiter", type=int, default=25)
    ap.add_argument("--B", type=float, default=B_PLUS)
    ap.add_argument("--no-absorption", action="store_true")
    ap.add_argument("--tag-suffix", default="")
    ap.add_argument(
        "--q-json",
        default=None,
        help="Optional JSON with best.q or q array (for --start custom)",
    )
    ap.add_argument("--start-label", default=None)
    args = ap.parse_args()
    tag_map = {
        "uniform": "uniform",
        "pref_peak": "opt_single_uniform_M_Bp0.0500",
        "pref_broadband": "opt_broadband_C_Bp0.0500",
    }
    orbits, _, _ = orbits_and_masks()
    if args.start == "custom":
        if not args.q_json:
            raise SystemExit("--q-json required with --start custom")
        d = json.loads(Path(args.q_json).read_text())
        if "best" in d and "q" in d["best"]:
            q0 = np.asarray(d["best"]["q"], float)
        elif "q" in d:
            q0 = np.asarray(d["q"], float)
        else:
            raise SystemExit("q-json missing q")
        label = args.start_label or "custom"
        # monkey-patch via run_opt internals: write temporary start through expand
        from multimode_opt_core import expand_tied

        s0 = expand_tied(q0, orbits)
        # stash as a design file load_s can read
        tmp = OUT / "_tmp_custom_start.json"
        dump_json(tmp, {"best": {"s": s0.tolist(), "q": q0.tolist()}})
        # load_s expects DES path basename without .json — write into DES
        from multimode_opt_core import DES

        des_tmp = DES / "_tmp_custom_start.json"
        dump_json(des_tmp, {"best": {"s": s0.tolist(), "q": q0.tolist()}})
        run_opt(
            grade=args.grade,
            start_label=label,
            start_tag="_tmp_custom_start",
            maxiter=args.maxiter,
            b_tesla=args.B,
            with_absorption=not args.no_absorption,
            tag_suffix=args.tag_suffix,
        )
        return
    run_opt(
        grade=args.grade,
        start_label=args.start_label or args.start,
        start_tag=tag_map[args.start],
        maxiter=args.maxiter,
        b_tesla=args.B,
        with_absorption=not args.no_absorption,
        tag_suffix=args.tag_suffix,
    )


if __name__ == "__main__":
    main()
