#!/usr/bin/env python3
"""Feasibility study: band-limited A2 optimization (NOT Stage-6).

Same frozen physics package as Stage-5 (C6/16, A2_hybrid, τ/α, bounds, FULL91-C,
optimized kernel) but restricted deterministic 6.25 MHz grids on narrower bands.

Does not launch Stage-6. Does not change FEM / materials / ports.
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
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign" / "feasibility"
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
from multimode_opt_core import B_PLUS, S_HI, S_LO, orbits_and_masks  # noqa: E402
from objective_A2 import A2_FREEZE  # noqa: E402

TAU = 0.05
ALPHA = 0.25
BOUNDS = (S_LO, S_HI)
DF = 0.00625  # 6.25 MHz

BANDS = {
    "A_4p9_5p1": (4.9, 5.1),
    "B_4p75_5p25": (4.75, 5.25),
    "C_4p5_5p5": (4.5, 5.5),
    "D_4p25_5p75": (4.25, 5.75),
}


def band_freqs(f0: float, f1: float, df: float = DF) -> list[float]:
    n = int(round((f1 - f0) / df))
    freqs = [round(f0 + i * df, 6) for i in range(n + 1)]
    assert abs(freqs[0] - f0) < 1e-12 and abs(freqs[-1] - f1) < 1e-12
    assert abs(max(freqs[i + 1] - freqs[i] for i in range(len(freqs) - 1)) - df) < 1e-12
    return freqs


def _atomic_write(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2) + "\n")
    tmp.replace(path)


def load_q(start: str) -> tuple[np.ndarray, dict]:
    if start == "uniform":
        orbits, _, _ = orbits_and_masks()
        return np.ones(len(orbits), float), {"start": "uniform_q=1"}
    if start == "stage5":
        p = ROOT / "outputs/circulator_inverse_design/jesse_method_campaign/opt_broadband_stage5_A2_C.json"
        d = json.loads(p.read_text())
        return np.asarray(d["q_best"], float), {
            "start": "stage5_q_best",
            "J_BB_stage5": d.get("J_BB_best"),
            "source": str(p),
        }
    if start.startswith("file:"):
        p = Path(start[5:])
        d = json.loads(p.read_text())
        q = np.asarray(d.get("q_best", d.get("q")), float)
        return q, {"start": f"file:{p.name}", "source": str(p)}
    raise SystemExit(f"unknown --start {start}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", required=True, choices=list(BANDS) + ["custom"])
    ap.add_argument("--f0", type=float, default=None)
    ap.add_argument("--f1", type=float, default=None)
    ap.add_argument("--start", default="stage5", help="stage5|uniform|file:PATH")
    ap.add_argument("--tag", default=None, help="output tag suffix")
    ap.add_argument("--grade", default="C")
    ap.add_argument("--maxiter", type=int, default=20)
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--omp-threads", type=int, default=1)
    ap.add_argument("--b-tesla", type=float, default=None)
    ap.add_argument("--symmetry", choices=["C6", "C3"], default="C6")
    args = ap.parse_args()

    for k, v in (
        ("OMP_NUM_THREADS", str(args.omp_threads)),
        ("OPENBLAS_NUM_THREADS", str(args.omp_threads)),
        ("MKL_NUM_THREADS", str(args.omp_threads)),
        ("PYTHONUNBUFFERED", "1"),
        ("PLASMEEP_A2_FAST_KERNEL", "1"),
    ):
        os.environ.setdefault(k, v)

    if args.band == "custom":
        if args.f0 is None or args.f1 is None:
            raise SystemExit("custom band requires --f0 --f1")
        f0, f1 = float(args.f0), float(args.f1)
        band_name = f"custom_{f0:g}_{f1:g}"
    else:
        f0, f1 = BANDS[args.band]
        band_name = args.band

    b_tesla = float(B_PLUS if args.b_tesla is None else args.b_tesla)
    freqs = band_freqs(f0, f1)
    freqs_hz = [f * 1e9 for f in freqs]

    q0, start_meta = load_q(args.start)
    if args.symmetry == "C6":
        orbits, _, _ = orbits_and_masks()
        assert len(q0) == len(orbits)
    else:
        from full91_geometry import bulb_centers_device
        import sixport_common as sc
        from threefold_orbits import build_orbits_c3, map_c6_q_to_c3

        centers = bulb_centers_device(sc)
        orbits, _ = build_orbits_c3(centers)
        if len(q0) == 16:
            orbits_c6, _, _ = orbits_and_masks()
            q0 = map_c6_q_to_c3(q0, orbits_c6, orbits)
            start_meta = {**start_meta, "mapped": "C6→C3 orbit average"}
        assert len(q0) == len(orbits) == 31

    tag = args.tag or f"{band_name}_{args.symmetry}_B{b_tesla:.4f}_{args.start.replace(':','_')}"
    stem = OUT / tag
    ckpt = stem.with_name(stem.name + "_ckpt.json")
    final = stem.with_name(stem.name + "_opt.json")
    hist_path = stem.with_name(stem.name + "_history.json")

    print("=" * 72, flush=True)
    print("FEASIBILITY BAND OPT — Stage-6 NOT STARTED", flush=True)
    print(
        f"band={band_name} [{f0},{f1}] GHz n_freq={len(freqs)} df=6.25 MHz "
        f"symmetry={args.symmetry} n_var={len(q0)} B={b_tesla} T "
        f"start={start_meta.get('start')} maxiter={args.maxiter}",
        flush=True,
    )
    print(
        f"workers={args.n_workers} omp={args.omp_threads} A2={A2_FREEZE} "
        f"tau={TAU} alpha={ALPHA} bounds={BOUNDS}",
        flush=True,
    )
    print("=" * 72, flush=True)

    t_wall0 = time.perf_counter()
    n_eval = 0
    accepted_nit = 0
    history = []
    cache = {"q": None, "rec": None}
    best = {"J_BB": -1e300, "q": None, "rec": None}

    if ckpt.exists():
        d = json.loads(ckpt.read_text())
        if d.get("status") == "complete":
            print(f"already complete: {ckpt}", flush=True)
            return
        accepted_nit = int(d.get("accepted_nit", 0))
        best["J_BB"] = float(d.get("J_BB_best", best["J_BB"]))
        if d.get("q_best") is not None:
            best["q"] = np.asarray(d["q_best"], float)
            q0 = best["q"].copy()
        history = d.get("history") or []
        n_eval = int(d.get("n_eval", 0))
        print(f"RESUME nit={accepted_nit} J_BB={best['J_BB']:.6f} n_eval={n_eval}", flush=True)

    remaining = max(0, args.maxiter - accepted_nit)
    if remaining == 0:
        print("nothing to do", flush=True)
        return

    def _eval(q):
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
            b_tesla=b_tesla,
        )
        wall = time.perf_counter() - t0
        # compact metrics (drop huge arrays in checkpoint)
        mets = rec["metrics_per_freq"]
        compact = {
            "J_BB": float(rec["J_BB"]),
            "J_mean": float(rec["J_mean"]),
            "worst_J": float(rec["worst_J"]),
            "worst_f_GHz": float(rec["worst_f_GHz"]),
            "grad_norm": float(rec["grad_norm"]),
            "max_weight": float(rec["max_weight"]),
            "N_eff": float(rec["N_eff"]),
            "wall_s": float(wall),
            "worst_metrics": {
                "T_des": float(mets[int(np.argmin(np.abs(np.asarray(rec['freqs_GHz']) - rec['worst_f_GHz'])))]["T_des_avg"]),
                "IL_dB": float(rec["worst_metrics"]["IL_avg_dB"]),
                "isolation_dB": float(rec["worst_metrics"]["isolation_avg_dB"]),
            },
            "grad_q": np.asarray(rec["grad_q"], float),
            "freqs_GHz": list(rec["freqs_GHz"]),
            "J_per_freq": [float(x) for x in rec["J_per_freq"]],
            "metrics_per_freq": [
                {
                    "T_des_avg": float(m["T_des_avg"]),
                    "IL_avg_dB": float(m["IL_avg_dB"]),
                    "isolation_avg_dB": float(m["isolation_avg_dB"]),
                }
                for m in mets
            ],
        }
        print(
            f"EVAL #{n_eval} J_BB={compact['J_BB']:.6f} worstA2={compact['worst_J']:.4f} "
            f"@ {compact['worst_f_GHz']:.5f} IL={compact['worst_metrics']['IL_dB']:.2f} "
            f"iso={compact['worst_metrics']['isolation_dB']:.2f} "
            f"grad={compact['grad_norm']:.3e} wall={wall:.1f}s",
            flush=True,
        )
        cache = {"q": q.copy(), "rec": compact}
        return compact

    def fun(q):
        rec = _eval(q)
        J = rec["J_BB"]
        if J > best["J_BB"]:
            best["J_BB"] = J
            best["q"] = np.asarray(q, float).copy()
            best["rec"] = {k: v for k, v in rec.items() if k != "grad_q"}
            print(f"NEW BEST J_BB={J:.6f}", flush=True)
        return -J

    def jac(q):
        return -_eval(q)["grad_q"]

    def _save(status="in_progress", extra=None):
        payload = {
            "status": status,
            "study": "feasibility_bandwidth_ladder",
            "stage6_status": "NOT STARTED / AWAITING AUTHORIZATION",
            "band": band_name,
            "f0_GHz": f0,
            "f1_GHz": f1,
            "n_freq": len(freqs),
            "df_MHz": DF * 1e3,
            "symmetry": args.symmetry,
            "n_var": len(q0),
            "b_tesla": b_tesla,
            "start": start_meta,
            "A2": A2_FREEZE,
            "tau": TAU,
            "alpha": ALPHA,
            "bounds": list(BOUNDS),
            "maxiter": args.maxiter,
            "accepted_nit": accepted_nit,
            "n_eval": n_eval,
            "J_BB_best": best["J_BB"],
            "q_best": None if best["q"] is None else best["q"].tolist(),
            "rec_best": best["rec"],
            "history": history,
            "elapsed_s": time.perf_counter() - t_wall0,
            "n_workers": args.n_workers,
        }
        if extra:
            payload.update(extra)
        _atomic_write(ckpt, payload)
        _atomic_write(hist_path, {"history": history})

    def callback(xk):
        nonlocal accepted_nit
        accepted_nit += 1
        rec = _eval(xk)
        record = {
            "kind": "accepted_iter",
            "iteration": accepted_nit,
            "J_BB": rec["J_BB"],
            "worst_A2": rec["worst_J"],
            "worst_f_GHz": rec["worst_f_GHz"],
            "grad_norm": rec["grad_norm"],
            "IL_dB": rec["worst_metrics"]["IL_dB"],
            "isolation_dB": rec["worst_metrics"]["isolation_dB"],
            "T_des": rec["worst_metrics"]["T_des"],
            "wall_s_eval": rec["wall_s"],
            "q": np.asarray(xk, float).tolist(),
            "n_lo": int(np.sum(np.asarray(xk) <= S_LO + 0.02)),
            "n_hi": int(np.sum(np.asarray(xk) >= S_HI - 0.02)),
        }
        history.append(record)
        print(
            f"ACCEPTED iter={accepted_nit}/{args.maxiter} J_BB={rec['J_BB']:.6f} "
            f"worstA2={rec['worst_J']:.4f} @{rec['worst_f_GHz']:.5f} "
            f"IL={rec['worst_metrics']['IL_dB']:.2f} "
            f"iso={rec['worst_metrics']['isolation_dB']:.2f} "
            f"grad={rec['grad_norm']:.3e} bounds(lo,hi)=({record['n_lo']},{record['n_hi']})",
            flush=True,
        )
        _save()

    rec0 = _eval(q0)
    if best["q"] is None or rec0["J_BB"] > best["J_BB"]:
        best["J_BB"] = rec0["J_BB"]
        best["q"] = q0.copy()
        best["rec"] = {k: v for k, v in rec0.items() if k != "grad_q"}
    history.append(
        {
            "kind": "start_eval",
            "iteration": accepted_nit,
            "J_BB": rec0["J_BB"],
            "worst_A2": rec0["worst_J"],
            "worst_f_GHz": rec0["worst_f_GHz"],
            "grad_norm": rec0["grad_norm"],
            "q": q0.tolist(),
        }
    )
    _save()
    print(
        f"START J_BB={rec0['J_BB']:.6f} worstA2={rec0['worst_J']:.4f} "
        f"@{rec0['worst_f_GHz']:.5f}",
        flush=True,
    )

    res = minimize(
        fun,
        q0,
        method="L-BFGS-B",
        jac=jac,
        bounds=[BOUNDS] * len(q0),
        callback=callback,
        options={"maxiter": remaining, "maxfun": remaining * 6, "ftol": 1e-12},
    )
    wall = time.perf_counter() - t_wall0
    out = {
        "status": "complete",
        "study": "feasibility_bandwidth_ladder",
        "stage6_status": "NOT STARTED / AWAITING AUTHORIZATION",
        "band": band_name,
        "f0_GHz": f0,
        "f1_GHz": f1,
        "n_freq": len(freqs),
        "df_MHz": DF * 1e3,
        "symmetry": args.symmetry,
        "n_var": len(q0),
        "b_tesla": b_tesla,
        "start": start_meta,
        "A2": A2_FREEZE,
        "tau": TAU,
        "alpha": ALPHA,
        "bounds": list(BOUNDS),
        "maxiter": args.maxiter,
        "accepted_nit": accepted_nit,
        "n_eval": n_eval,
        "elapsed_s": wall,
        "optimizer_success": bool(res.success),
        "optimizer_message": str(res.message),
        "J_BB_best": best["J_BB"],
        "q_best": best["q"].tolist(),
        "rec_best": best["rec"],
        "history": history,
        "n_workers": args.n_workers,
        "q_near_lo": int(np.sum(best["q"] <= S_LO + 0.02)),
        "q_near_hi": int(np.sum(best["q"] >= S_HI - 0.02)),
    }
    _atomic_write(final, out)
    _save("complete", {"optimizer_message": str(res.message)})
    print(
        f"DONE FEASIBILITY {tag} J_BB={best['J_BB']:.6f} nit={accepted_nit} "
        f"n_eval={n_eval} elapsed={wall:.0f}s msg={res.message}",
        flush=True,
    )
    print("STAGE-6: NOT STARTED / AWAITING AUTHORIZATION", flush=True)


if __name__ == "__main__":
    main()
