#!/usr/bin/env python3
"""Forward-only modal autopsy of dense-audit A2 notches.

No optimization. Stage-6 NOT STARTED.

Diagnoses whether notches are predominantly:
  A) desired TE0 -> TE1 conversion
  B) reflection resonance
  C) wrong-port TE0 routing
  D) mixed
  E) inconsistent power accounting / numerics
"""
from __future__ import annotations

import argparse
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
FEAS = OUT / "feasibility"
AUTO = FEAS / "notch_autopsy"
sys_paths = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "jesse_method_campaign"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sys

sys.path[:0] = sys_paths

from eval_forward_only import eval_q_forward  # noqa: E402
from matplotlib.tri import LinearTriInterpolator, Triangulation  # noqa: E402
from mode_objectives import _port_powers  # noqa: E402
from multimode_opt_core import orbits_and_masks  # noqa: E402
from objective_A2 import objective_A2_from_S  # noqa: E402

DESIGNS = {
    "A_uniform": {
        "label": "Band A 0.2 GHz B=0.05 uniform-start (clean control)",
        "opt": FEAS / "A_4p9_5p1_C6_B0.0500_uniform_opt.json",
        "audit": FEAS / "audit_A_4p9_5p1_C6_B0.0500_uniform_opt_3.125MHz.json",
        "b_tesla": 0.05,
        "band": (4.9, 5.1),
    },
    "B_stage5": {
        "label": "Band B 0.5 GHz B=0.05 C6 Stage-5 (borderline)",
        "opt": FEAS / "B_4p75_5p25_C6_B0.0500_stage5_opt.json",
        "audit": FEAS / "audit_B_4p75_5p25_C6_B0.0500_stage5_opt_1.5625MHz.json",
        "b_tesla": 0.05,
        "band": (4.75, 5.25),
    },
    "C_B010_Bs5": {
        "label": "Band C 1.0 GHz B=0.10 from Band-B start (best wider-band)",
        "opt": FEAS / "C_4p5_5p5_C6_B0.1000_Bs5_opt.json",
        "audit": FEAS / "audit_C_4p5_5p5_C6_B0.1000_Bs5_opt_3.125MHz.json",
        "b_tesla": 0.10,
        "band": (4.5, 5.5),
    },
}


def _init_worker(omp: int = 1):
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = str(int(omp))
    os.environ.setdefault("PLASMEEP_A2_FAST_KERNEL", "1")


def _avg_port_powers(S: np.ndarray) -> dict:
    keys = [
        "des0",
        "des1",
        "rev0",
        "rev1",
        "thr0",
        "thr1",
        "other0",
        "other1",
        "R0",
        "R1",
    ]
    acc = {k: 0.0 for k in keys}
    per = []
    for j in range(6):
        p = _port_powers(S, j)
        per.append({"source": j, **{k: float(p[k]) for k in keys}})
        for k in keys:
            acc[k] += float(p[k])
    for k in keys:
        acc[k] /= 6.0
    te1_by_port = {
        "des1": acc["des1"],
        "rev1": acc["rev1"],
        "thr1": acc["thr1"],
        "other1": acc["other1"],
        "R1": acc["R1"],
    }
    te0_wrong = acc["rev0"] + acc["thr0"] + acc["other0"]
    te1_total = sum(te1_by_port.values())
    refl = acc["R0"] + acc["R1"]
    pout = sum(acc.values())
    return {
        "avg": acc,
        "te1_by_role": te1_by_port,
        "TE0_des": acc["des0"],
        "TE0_rev": acc["rev0"],
        "TE0_thr": acc["thr0"],
        "TE0_other": acc["other0"],
        "TE0_wrong": te0_wrong,
        "TE1_total": te1_total,
        "TE1_at_des": acc["des1"],
        "R0": acc["R0"],
        "R1": acc["R1"],
        "R_tot": refl,
        "P_out": pout,
        "missing": 1.0 - pout,
        "per_source": per,
    }


def classify(diag: dict) -> dict:
    """Heuristic classification from average modal power fractions."""
    d = diag
    des0 = d["TE0_des"]
    te1 = d["TE1_total"]
    rtot = d["R_tot"]
    wrong = d["TE0_wrong"]
    miss = abs(d["missing"])
    # Dominant deficit relative to a healthy Tdes~0.3–0.5 circulator
    # Use absolute fractions of incident (unit) power.
    scores = {
        "A_TE1_conversion": te1,
        "B_reflection": rtot,
        "C_wrong_TE0_port": wrong,
    }
    # If desired TE0 is OK but isolation bad, wrong-port still matters
    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    top, top_v = ranked[0]
    second_v = ranked[1][1]
    if miss > 0.15:
        label = "E_power_budget_inconsistent"
        reason = f"|1-P_out|={miss:.3f} exceeds 0.15"
    elif top_v < 0.12 and des0 > 0.15:
        label = "ok_or_mild"
        reason = f"no large loss channel; des0={des0:.3f}"
    elif second_v >= 0.55 * top_v and top_v >= 0.15:
        label = "D_mixed"
        reason = f"top={top}:{top_v:.3f}, second={ranked[1][0]}:{second_v:.3f}"
    else:
        label = top
        reason = f"dominant {top}={top_v:.3f} (des0={des0:.3f})"
    return {
        "class": label,
        "reason": reason,
        "scores": scores,
        "des0": des0,
        "iso_te0_proxy_dB": float(
            10.0 * np.log10(max(des0, 1e-30) / max(d["TE0_rev"], 1e-30))
        ),
    }


def diagnose_bundle(f_ghz: float, bun: dict) -> dict:
    S = bun["S"]
    met = bun["metrics"]
    powers = _avg_port_powers(S)
    # worst source by A2-relevant isolation / IL
    per = met["per_source"]
    worst_iso = min(per, key=lambda r: r["isolation_dB"])
    worst_il = max(per, key=lambda r: r["IL_dB"])
    cls = classify(powers)
    return {
        "f_GHz": float(f_ghz),
        "J_A2": float(objective_A2_from_S(S)),
        "J_A": float(bun["J_A"]),
        "J_B": float(bun["J_B"]),
        "IL_avg_dB": float(met["IL_avg_dB"]),
        "isolation_avg_dB": float(met["isolation_avg_dB"]),
        "RL_avg_dB": float(met["RL_avg_dB"]),
        "T_des_avg": float(met["T_des_avg"]),
        "T_rev_avg": float(met["T_rev_avg"]),
        "R_avg_te0": float(met["R_avg"]),
        "TE1_avg": float(met["TE1_avg"]),
        "P_out_avg_metrics": float(met["P_out_avg"]),
        "powers": {k: v for k, v in powers.items() if k != "per_source"},
        "powers_per_source": powers["per_source"],
        "worst_iso_source": {
            "source": int(worst_iso["source"]),
            "T_des": float(worst_iso["T_des"]),
            "T_rev": float(worst_iso["T_rev"]),
            "R": float(worst_iso["R"]),
            "TE1_total": float(worst_iso["TE1_total"]),
            "P_out": float(worst_iso["P_out"]),
            "IL_dB": float(worst_iso["IL_dB"]),
            "isolation_dB": float(worst_iso["isolation_dB"]),
        },
        "worst_il_source": {
            "source": int(worst_il["source"]),
            "T_des": float(worst_il["T_des"]),
            "T_rev": float(worst_il["T_rev"]),
            "R": float(worst_il["R"]),
            "TE1_total": float(worst_il["TE1_total"]),
            "P_out": float(worst_il["P_out"]),
            "IL_dB": float(worst_il["IL_dB"]),
            "isolation_dB": float(worst_il["isolation_dB"]),
        },
        "classification": cls,
        "timings": bun.get("timings"),
    }


def _worker(payload: dict) -> dict:
    _init_worker(int(payload.get("omp_threads", 1)))
    from multimode_opt_core import orbits_and_masks as _oam
    from eval_forward_only import eval_q_forward as _eqf

    q = np.asarray(payload["q"], float)
    orbits, _, _ = _oam()
    bun, _ = _eqf(
        payload["grade"],
        q,
        orbits,
        f_hz=float(payload["f_GHz"]) * 1e9,
        b_tesla=float(payload["b_tesla"]),
    )
    # Drop heavy objects before classify (re-import helpers in worker)
    from feasibility_notch_autopsy import diagnose_bundle  # type: ignore

    out = diagnose_bundle(float(payload["f_GHz"]), bun)
    out["design"] = payload["design"]
    out["role"] = payload["role"]
    out["b_tesla"] = float(payload["b_tesla"])
    out["note"] = payload.get("note", "")
    return out


def dither(f0: float, deltas_mhz=( -2.0, -1.0, 0.0, 1.0, 2.0)) -> list[float]:
    return [round(f0 + d * 1e-3, 6) for d in deltas_mhz]


def pick_sites() -> dict:
    """Choose notch + control frequencies from existing dense audits."""
    sites = {}

    # --- A ---
    rows = json.loads(DESIGNS["A_uniform"]["audit"].read_text())["rows"]
    rows = sorted(rows, key=lambda r: r["f_GHz"])
    soft = min(rows, key=lambda r: r["J_A2"])  # 5.090625
    best = max(rows, key=lambda r: r["J_A2"])
    # nearest milder neighbors as controls
    sites["A_uniform"] = {
        "notches": [
            {
                "f_GHz": soft["f_GHz"],
                "role": "soft_min",
                "audit_A2": soft["J_A2"],
                "dither": True,
            }
        ],
        "controls": [
            {"f_GHz": best["f_GHz"], "role": "best_in_band", "audit_A2": best["J_A2"]},
            {"f_GHz": 5.0, "role": "band_center", "audit_A2": None},
            {"f_GHz": 5.05625, "role": "near_soft_side", "audit_A2": None},
        ],
    }

    # --- B ---
    rows = json.loads(DESIGNS["B_stage5"]["audit"].read_text())["rows"]
    rows = sorted(rows, key=lambda r: r["f_GHz"])
    a = np.array([r["J_A2"] for r in rows])
    local = []
    for i in range(1, len(a) - 1):
        if a[i] <= a[i - 1] and a[i] <= a[i + 1]:
            local.append(rows[i])
    local = sorted(local, key=lambda r: r["J_A2"])
    best = max(rows, key=lambda r: r["J_A2"])
    # top scientifically interesting: two deepest A2 mins + a high-R min
    notches = []
    for r in local[:2]:
        notches.append(
            {"f_GHz": r["f_GHz"], "role": "A2_local_min", "audit_A2": r["J_A2"], "dither": True}
        )
    r_peak = max(local, key=lambda r: r["R"])
    if abs(r_peak["f_GHz"] - notches[0]["f_GHz"]) > 0.01:
        notches.append(
            {
                "f_GHz": r_peak["f_GHz"],
                "role": "high_R_local_min",
                "audit_A2": r_peak["J_A2"],
                "dither": False,
            }
        )
    # side controls near deepest notch
    f0 = notches[0]["f_GHz"]
    left = [r for r in rows if r["f_GHz"] < f0 - 0.02]
    right = [r for r in rows if r["f_GHz"] > f0 + 0.02]
    controls = [{"f_GHz": best["f_GHz"], "role": "best_in_band", "audit_A2": best["J_A2"]}]
    if left:
        bl = max(left, key=lambda r: r["J_A2"])
        controls.append({"f_GHz": bl["f_GHz"], "role": "left_of_notch", "audit_A2": bl["J_A2"]})
    if right:
        br = max(right, key=lambda r: r["J_A2"])
        controls.append({"f_GHz": br["f_GHz"], "role": "right_of_notch", "audit_A2": br["J_A2"]})
    sites["B_stage5"] = {"notches": notches, "controls": controls}

    # --- C ---
    rows = json.loads(DESIGNS["C_B010_Bs5"]["audit"].read_text())["rows"]
    rows = sorted(rows, key=lambda r: r["f_GHz"])
    a = np.array([r["J_A2"] for r in rows])
    local = []
    for i in range(1, len(a) - 1):
        if a[i] <= a[i - 1] and a[i] <= a[i + 1]:
            local.append(rows[i])
    local = sorted(local, key=lambda r: r["J_A2"])
    best = max(rows, key=lambda r: r["J_A2"])
    deep = local[0]
    notches = [
        {
            "f_GHz": deep["f_GHz"],
            "role": "A2_le_m5_notch",
            "audit_A2": deep["J_A2"],
            "dither": True,
        }
    ]
    # secondary: high-R min and next A2 min if distinct
    for r, role in [
        (local[1], "A2_secondary_min"),
        (max(local[:5], key=lambda x: x["R"]), "high_R_local_min"),
    ]:
        if all(abs(r["f_GHz"] - n["f_GHz"]) > 0.015 for n in notches):
            notches.append(
                {"f_GHz": r["f_GHz"], "role": role, "audit_A2": r["J_A2"], "dither": False}
            )
    f0 = deep["f_GHz"]
    left = [r for r in rows if 4.55 <= r["f_GHz"] < f0 - 0.015]
    right = [r for r in rows if f0 + 0.015 < r["f_GHz"] <= 4.70]
    controls = [{"f_GHz": best["f_GHz"], "role": "best_in_band", "audit_A2": best["J_A2"]}]
    if left:
        bl = max(left, key=lambda r: r["J_A2"])
        controls.append({"f_GHz": bl["f_GHz"], "role": "left_of_notch", "audit_A2": bl["J_A2"]})
    if right:
        br = max(right, key=lambda r: r["J_A2"])
        controls.append({"f_GHz": br["f_GHz"], "role": "right_of_notch", "audit_A2": br["J_A2"]})
    sites["C_B010_Bs5"] = {"notches": notches, "controls": controls}
    return sites


def build_jobs(sites: dict) -> list[dict]:
    jobs = []
    qs = {}
    for key, meta in DESIGNS.items():
        d = json.loads(meta["opt"].read_text())
        qs[key] = np.asarray(d["q_best"], float).tolist()

    def add(design, f, role, b=None, note="", dither_parent=None):
        jobs.append(
            {
                "design": design,
                "q": qs[design],
                "f_GHz": float(f),
                "b_tesla": float(DESIGNS[design]["b_tesla"] if b is None else b),
                "grade": "C",
                "role": role,
                "note": note,
                "omp_threads": 1,
                "dither_parent": dither_parent,
            }
        )

    for design, pack in sites.items():
        for n in pack["notches"]:
            if n.get("dither"):
                for f in dither(n["f_GHz"]):
                    role = "notch" if abs(f - n["f_GHz"]) < 1e-9 else "dither"
                    add(
                        design,
                        f,
                        f"{role}:{n['role']}",
                        note=f"audit_A2={n['audit_A2']}",
                        dither_parent=n["f_GHz"],
                    )
            else:
                add(design, n["f_GHz"], f"notch:{n['role']}", note=f"audit_A2={n['audit_A2']}")
        for c in pack["controls"]:
            add(design, c["f_GHz"], f"control:{c['role']}", note=f"audit_A2={c['audit_A2']}")

    # Cross-band / B-swap forward probes (fixed q, no opt)
    # Does Band-A design host the Band-B/C notch frequencies?
    for f, note in [
        (4.609375, "C_notch_freq"),
        (4.909375, "B_deepest_min_freq"),
        (5.164062, "B_second_min_freq"),
    ]:
        add("A_uniform", f, "crossband_probe", note=note)
    # Band-B design at C notch
    add("B_stage5", 4.609375, "crossband_probe", note="C_notch_freq")
    # Fixed-q B swap on C design around its notch
    for f in dither(4.609375):
        add(
            "C_B010_Bs5",
            f,
            "Bswap_0p05:" + ("notch" if abs(f - 4.609375) < 1e-9 else "dither"),
            b=0.05,
            note="same_q_at_B=0.05",
            dither_parent=4.609375,
        )
    # dedupe exact (design,f,b,role-prefix) keeping first
    seen = set()
    uniq = []
    for j in jobs:
        key = (j["design"], round(j["f_GHz"], 6), round(j["b_tesla"], 4), j["role"].split(":")[0])
        # allow both notch and dither etc — use full role
        key2 = (j["design"], round(j["f_GHz"], 6), round(j["b_tesla"], 4), j["role"], j.get("note"))
        if key2 in seen:
            continue
        seen.add(key2)
        uniq.append(j)
    return uniq


def run_jobs(jobs: list[dict], n_workers: int = 12) -> list[dict]:
    AUTO.mkdir(parents=True, exist_ok=True)
    results = []
    t0 = time.time()
    print(f"NOTCH AUTOPSY: {len(jobs)} forward evals, n_workers={n_workers}", flush=True)
    # Import diagnose in workers via module path
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = {ex.submit(_worker_standalone, j): j for j in jobs}
        done = 0
        for fut in as_completed(futs):
            done += 1
            j = futs[fut]
            try:
                r = fut.result()
            except Exception as e:
                r = {
                    "design": j["design"],
                    "f_GHz": j["f_GHz"],
                    "b_tesla": j["b_tesla"],
                    "role": j["role"],
                    "error": repr(e),
                }
                print(f"FAIL {j['design']} f={j['f_GHz']} {e}", flush=True)
            results.append(r)
            if done % 5 == 0 or done == len(jobs):
                print(
                    f"  [{done}/{len(jobs)}] last={r.get('design')} f={r.get('f_GHz')} "
                    f"A2={r.get('J_A2', 'err')} class={((r.get('classification') or {}).get('class'))}",
                    flush=True,
                )
                (AUTO / "notch_autopsy_partial.json").write_text(
                    json.dumps({"rows": results, "elapsed_s": time.time() - t0}, indent=2) + "\n"
                )
    results = sorted(
        results,
        key=lambda r: (r.get("design") or "", r.get("b_tesla") or 0, r.get("f_GHz") or 0, r.get("role") or ""),
    )
    return results


def _worker_standalone(payload: dict) -> dict:
    """Worker that does not import this module (avoids pickling issues)."""
    _init_worker(int(payload.get("omp_threads", 1)))
    import numpy as np
    from multimode_opt_core import orbits_and_masks, te0_metrics_from_S
    from eval_forward_only import eval_q_forward
    from mode_objectives import _port_powers, objective_A_from_S, objective_B_from_S
    from objective_A2 import objective_A2_from_S

    q = np.asarray(payload["q"], float)
    orbits, _, _ = orbits_and_masks()
    bun, _ = eval_q_forward(
        payload["grade"],
        q,
        orbits,
        f_hz=float(payload["f_GHz"]) * 1e9,
        b_tesla=float(payload["b_tesla"]),
    )
    S = bun["S"]
    met = bun["metrics"]
    keys = ["des0", "des1", "rev0", "rev1", "thr0", "thr1", "other0", "other1", "R0", "R1"]
    acc = {k: 0.0 for k in keys}
    per_src = []
    for j in range(6):
        p = _port_powers(S, j)
        per_src.append({"source": j, **{k: float(p[k]) for k in keys}})
        for k in keys:
            acc[k] += float(p[k])
    for k in keys:
        acc[k] /= 6.0
    te1_total = acc["des1"] + acc["rev1"] + acc["thr1"] + acc["other1"] + acc["R1"]
    wrong = acc["rev0"] + acc["thr0"] + acc["other0"]
    rtot = acc["R0"] + acc["R1"]
    pout = sum(acc.values())
    miss = 1.0 - pout
    scores = {
        "A_TE1_conversion": te1_total,
        "B_reflection": rtot,
        "C_wrong_TE0_port": wrong,
    }
    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    top, top_v = ranked[0]
    second_v = ranked[1][1]
    des0 = acc["des0"]
    if miss > 0.15:
        label = "E_power_budget_inconsistent"
        reason = f"|1-P_out|={miss:.3f} exceeds 0.15"
    elif top_v < 0.12 and des0 > 0.15:
        label = "ok_or_mild"
        reason = f"no large loss channel; des0={des0:.3f}"
    elif second_v >= 0.55 * top_v and top_v >= 0.15:
        label = "D_mixed"
        reason = f"top={top}:{top_v:.3f}, second={ranked[1][0]}:{second_v:.3f}"
    else:
        label = top
        reason = f"dominant {top}={top_v:.3f} (des0={des0:.3f})"
    per = met["per_source"]
    worst_iso = min(per, key=lambda r: r["isolation_dB"])
    worst_il = max(per, key=lambda r: r["IL_dB"])
    return {
        "design": payload["design"],
        "f_GHz": float(payload["f_GHz"]),
        "b_tesla": float(payload["b_tesla"]),
        "role": payload["role"],
        "note": payload.get("note", ""),
        "dither_parent": payload.get("dither_parent"),
        "J_A2": float(objective_A2_from_S(S)),
        "J_A": float(objective_A_from_S(S)),
        "J_B": float(objective_B_from_S(S)),
        "IL_avg_dB": float(met["IL_avg_dB"]),
        "isolation_avg_dB": float(met["isolation_avg_dB"]),
        "RL_avg_dB": float(met["RL_avg_dB"]),
        "T_des_avg": float(met["T_des_avg"]),
        "T_rev_avg": float(met["T_rev_avg"]),
        "R_avg_te0": float(met["R_avg"]),
        "TE1_avg": float(met["TE1_avg"]),
        "P_out_avg_metrics": float(met["P_out_avg"]),
        "powers": {
            "avg": acc,
            "TE0_des": des0,
            "TE0_rev": acc["rev0"],
            "TE0_thr": acc["thr0"],
            "TE0_other": acc["other0"],
            "TE0_wrong": wrong,
            "TE1_total": te1_total,
            "TE1_at_des": acc["des1"],
            "TE1_by_role": {
                "des1": acc["des1"],
                "rev1": acc["rev1"],
                "thr1": acc["thr1"],
                "other1": acc["other1"],
                "R1": acc["R1"],
            },
            "R0": acc["R0"],
            "R1": acc["R1"],
            "R_tot": rtot,
            "P_out": pout,
            "missing": miss,
        },
        "powers_per_source": per_src,
        "worst_iso_source": {
            "source": int(worst_iso["source"]),
            "T_des": float(worst_iso["T_des"]),
            "T_rev": float(worst_iso["T_rev"]),
            "R": float(worst_iso["R"]),
            "TE1_total": float(worst_iso["TE1_total"]),
            "P_out": float(worst_iso["P_out"]),
            "IL_dB": float(worst_iso["IL_dB"]),
            "isolation_dB": float(worst_iso["isolation_dB"]),
        },
        "classification": {
            "class": label,
            "reason": reason,
            "scores": scores,
            "des0": des0,
            "iso_te0_proxy_dB": float(
                10.0 * np.log10(max(des0, 1e-30) / max(acc["rev0"], 1e-30))
            ),
        },
        "timings": bun.get("timings"),
    }


def export_field_maps(design_key: str, freqs: list[tuple[float, str]], b_tesla: float):
    """Save comparable |Hz| maps for Band-C notch vs good frequency."""
    from adjoint_direct_dtn import build_system
    from direct_dtn_s import build_Q_columns
    from efficient_adjoint import S_from_U, _pack_C, _solve_multi
    from multimode_opt_core import expand_tied, orbits_and_masks, te0_metrics_from_S
    from scipy.sparse.linalg import splu

    _init_worker(1)
    d = json.loads(DESIGNS[design_key]["opt"].read_text())
    q = np.asarray(d["q_best"], float)
    orbits, _, _ = orbits_and_masks()
    s = expand_tied(q, orbits)
    out_dirs = []
    # Determine shared color scale after both solves
    grids = []
    for f_ghz, tag in freqs:
        f_hz = f_ghz * 1e9
        st = build_system("C", s, f_hz, b_tesla, "projected")
        points, A = st["points"], st["A"]
        cols = build_Q_columns(points, st["ports"], st["modes"], n_prop=2)
        C, betas = _pack_C(points, cols)
        te0_cols = [2 * j for j in range(6)]
        lu = splu(A.tocsc())
        U = _solve_multi(lu, C[:, te0_cols], trans="N")
        S = S_from_U(C, U, betas, te0_cols)
        met = te0_metrics_from_S(S)
        # Use worst-IL source for visualization of failure mode
        per = met["per_source"]
        j_src = int(max(per, key=lambda r: r["IL_dB"])["source"])
        Hz = np.asarray(U[:, j_src], complex)
        # interpolate |Hz|
        tri = Triangulation(points[:, 0], points[:, 1], st["tris"])
        xmin, xmax = float(points[:, 0].min()), float(points[:, 0].max())
        ymin, ymax = float(points[:, 1].min()), float(points[:, 1].max())
        nx, ny = 420, 400
        x = np.linspace(xmin, xmax, nx)
        y = np.linspace(ymin, ymax, ny)
        X, Y = np.meshgrid(x, y)
        abs_i = LinearTriInterpolator(tri, np.abs(Hz))
        Mag = np.ma.filled(abs_i(X, Y), np.nan).astype(np.float32)
        grids.append(
            {
                "f_GHz": f_ghz,
                "tag": tag,
                "x": x,
                "y": y,
                "Mag": Mag,
                "source": j_src + 1,
                "met": {
                    "IL_avg_dB": float(met["IL_avg_dB"]),
                    "isolation_avg_dB": float(met["isolation_avg_dB"]),
                    "T_des_avg": float(met["T_des_avg"]),
                    "T_rev_avg": float(met["T_rev_avg"]),
                    "R_avg": float(met["R_avg"]),
                    "TE1_avg": float(met["TE1_avg"]),
                    "J_A2": float(objective_A2_from_S(S)),
                },
            }
        )
        print(
            f"FIELD {tag} f={f_ghz} src=P{j_src+1} IL={met['IL_avg_dB']:.2f} "
            f"TE1={met['TE1_avg']:.3f} R={met['R_avg']:.3f}",
            flush=True,
        )

    vmax = float(
        np.nanpercentile(
            np.concatenate([g["Mag"][~np.isnan(g["Mag"])].ravel() for g in grids]), 99.0
        )
    )
    fig, axes = plt.subplots(1, len(grids), figsize=(5.2 * len(grids), 4.6), constrained_layout=True)
    if len(grids) == 1:
        axes = [axes]
    for ax, g in zip(axes, grids):
        im = ax.imshow(
            g["Mag"],
            origin="lower",
            extent=[g["x"][0], g["x"][-1], g["y"][0], g["y"][-1]],
            cmap="magma",
            vmin=0.0,
            vmax=vmax,
            aspect="equal",
        )
        ax.set_title(
            f"{g['tag']}\n{g['f_GHz']:.6f} GHz | P{g['source']} |A2|={g['met']['J_A2']:.2f}\n"
            f"Tdes={g['met']['T_des_avg']:.3f} Trev={g['met']['T_rev_avg']:.3f} "
            f"R={g['met']['R_avg']:.3f} TE1={g['met']['TE1_avg']:.3f}",
            fontsize=9,
        )
        ax.set_xlabel("x [cm]")
        ax.set_ylabel("y [cm]")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="|Hz| (arb.)")
    fig.suptitle(
        f"Band-C B=0.10 Band-B-start — |Hz| notch vs good (identical scale, vmax={vmax:.3g})",
        fontsize=11,
    )
    png = AUTO / "C_B010_Bs5_field_notch_vs_good.png"
    fig.savefig(png, dpi=150)
    plt.close(fig)
    # also save npz compact
    np.savez_compressed(
        AUTO / "C_B010_Bs5_field_notch_vs_good.npz",
        **{f"{g['tag']}_Mag": g["Mag"] for g in grids},
        x=grids[0]["x"],
        y=grids[0]["y"],
        freqs=np.array([g["f_GHz"] for g in grids]),
        tags=np.array([g["tag"] for g in grids]),
        vmax=vmax,
    )
    meta = {
        "design": design_key,
        "b_tesla": b_tesla,
        "vmax_absHz": vmax,
        "frames": [
            {"tag": g["tag"], "f_GHz": g["f_GHz"], "source_port": g["source"], "metrics": g["met"]}
            for g in grids
        ],
        "png": str(png),
        "note": "Identical magma scale; source = worst-IL TE0 drive for that frequency.",
    }
    (AUTO / "C_B010_Bs5_field_notch_vs_good_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def write_report(sites, results):
    by = {}
    for r in results:
        if "error" in r:
            continue
        by.setdefault(r["design"], []).append(r)

    lines = [
        "# Feasibility notch autopsy (forward-only)",
        "",
        "**Stage-6: NOT STARTED / AWAITING AUTHORIZATION**",
        "**No new optimization.**",
        "",
        "Modal powers are averages over the six TE0 source directions, extracted from the",
        "production 12×12 eigenmode-port S matrix (`|S_{mode,port; TE0 source}|^2`).",
        "`P_out` is the sum of all outgoing TE0+TE1 modal powers (power-budget closure vs 1).",
        "",
        "Classification heuristic uses absolute power fractions:",
        "- **A_TE1_conversion** — large total TE1",
        "- **B_reflection** — large R0+R1",
        "- **C_wrong_TE0_port** — large rev0+thr0+other0",
        "- **D_mixed** — two channels both large",
        "- **E_power_budget_inconsistent** — |1−P_out| > 0.15",
        "",
        "## Numerical modal-power table",
        "",
        "| design | B | f [GHz] | role | A2 | des0 | rev0 | wrong0 | R0 | R1 | TE1 | P_out | miss | class |",
        "|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for r in results:
        if "error" in r:
            lines.append(
                f"| {r.get('design')} | {r.get('b_tesla')} | {r.get('f_GHz')} | {r.get('role')} | ERR | | | | | | | | | {r.get('error')} |"
            )
            continue
        p = r["powers"]
        cls = r["classification"]["class"]
        lines.append(
            f"| {r['design']} | {r['b_tesla']:.3f} | {r['f_GHz']:.6f} | {r['role']} | "
            f"{r['J_A2']:.3f} | {p['TE0_des']:.3f} | {p['TE0_rev']:.3f} | {p['TE0_wrong']:.3f} | "
            f"{p['R0']:.3f} | {p['R1']:.3f} | {p['TE1_total']:.3f} | {p['P_out']:.3f} | "
            f"{p['missing']:.3f} | {cls} |"
        )

    # Dither smoothness
    lines += ["", "## Dither smoothness (±1, ±2 MHz)", ""]
    for design, rows in by.items():
        dithers = [r for r in rows if r.get("dither_parent") is not None]
        parents = sorted({r["dither_parent"] for r in dithers if r.get("dither_parent") is not None})
        for parent in parents:
            seq = sorted(
                [r for r in dithers if abs(r.get("dither_parent", -1) - parent) < 1e-9 and abs(r["b_tesla"] - DESIGNS[design]["b_tesla"]) < 1e-9],
                key=lambda r: r["f_GHz"],
            )
            # also include Bswap separately
            if not seq:
                continue
            a2 = [r["J_A2"] for r in seq]
            te1 = [r["powers"]["TE1_total"] for r in seq]
            wrong = [r["powers"]["TE0_wrong"] for r in seq]
            rtot = [r["powers"]["R_tot"] for r in seq]
            smooth = (
                max(a2) - min(a2) < 8.0
                and all(np.isfinite(a2))
                and (np.max(np.abs(np.diff(te1))) < 0.35 if len(te1) > 1 else True)
            )
            lines.append(
                f"- **{design}** parent {parent:.6f} GHz: A2={['%.3f'%x for x in a2]}, "
                f"TE1={['%.3f'%x for x in te1]}, wrong0={['%.3f'%x for x in wrong]}, "
                f"Rtot={['%.3f'%x for x in rtot]}; "
                f"{'smooth/continuous' if smooth else 'abrupt — inspect'} redistribution"
            )
        # Bswap dither for C
        if design == "C_B010_Bs5":
            seq = sorted(
                [r for r in rows if str(r["role"]).startswith("Bswap_0p05") and r.get("dither_parent")],
                key=lambda r: r["f_GHz"],
            )
            if seq:
                lines.append(
                    f"- **{design} B=0.05 swap** parent 4.609375: A2={['%.3f'%r['J_A2'] for r in seq]}, "
                    f"TE1={['%.3f'%r['powers']['TE1_total'] for r in seq]}, "
                    f"wrong0={['%.3f'%r['powers']['TE0_wrong'] for r in seq]}, "
                    f"Rtot={['%.3f'%r['powers']['R_tot'] for r in seq]}"
                )

    # Interpretation scaffold filled from numbers
    def row_at(design, f, b=None, role_substr=None):
        cands = by.get(design, [])
        out = []
        for r in cands:
            if abs(r["f_GHz"] - f) > 1e-6:
                continue
            if b is not None and abs(r["b_tesla"] - b) > 1e-9:
                continue
            if role_substr and role_substr not in r["role"]:
                continue
            out.append(r)
        return out[0] if out else None

    a_soft = row_at("A_uniform", sites["A_uniform"]["notches"][0]["f_GHz"], role_substr="notch")
    a_best = None
    for r in by.get("A_uniform", []):
        if r["role"].startswith("control:best"):
            a_best = r
            break
    b_deep = row_at("B_stage5", sites["B_stage5"]["notches"][0]["f_GHz"], role_substr="notch")
    c_notch = row_at("C_B010_Bs5", 4.609375, b=0.10, role_substr="notch")
    c_good = None
    for r in by.get("C_B010_Bs5", []):
        if r["role"].startswith("control:best") and abs(r["b_tesla"] - 0.10) < 1e-9:
            c_good = r
            break
    a_at_c = row_at("A_uniform", 4.609375, role_substr="crossband")
    a_at_b = row_at("A_uniform", 4.909375, role_substr="crossband")
    c_swap = row_at("C_B010_Bs5", 4.609375, b=0.05, role_substr="Bswap")

    lines += [
        "",
        "## Physical interpretation",
        "",
        "### 1. What are the narrow failures?",
    ]
    if c_notch:
        lines.append(
            f"- Band-C deep notch @ {c_notch['f_GHz']:.6f} GHz classifies as "
            f"**{c_notch['classification']['class']}** "
            f"(des0={c_notch['powers']['TE0_des']:.3f}, wrong0={c_notch['powers']['TE0_wrong']:.3f}, "
            f"Rtot={c_notch['powers']['R_tot']:.3f}, TE1={c_notch['powers']['TE1_total']:.3f}, "
            f"P_out={c_notch['powers']['P_out']:.3f}). "
            f"{c_notch['classification']['reason']}."
        )
    if b_deep:
        lines.append(
            f"- Band-B deepest min @ {b_deep['f_GHz']:.6f} GHz → "
            f"**{b_deep['classification']['class']}** "
            f"(des0={b_deep['powers']['TE0_des']:.3f}, wrong0={b_deep['powers']['TE0_wrong']:.3f}, "
            f"Rtot={b_deep['powers']['R_tot']:.3f}, TE1={b_deep['powers']['TE1_total']:.3f})."
        )
    if a_soft:
        lines.append(
            f"- Band-A softest in-band min @ {a_soft['f_GHz']:.6f} GHz → "
            f"**{a_soft['classification']['class']}** "
            f"(des0={a_soft['powers']['TE0_des']:.3f}, Rtot={a_soft['powers']['R_tot']:.3f}, "
            f"TE1={a_soft['powers']['TE1_total']:.3f}) — not an A2≤−5 corpse."
        )

    lines += ["", "### 2. Does Band A avoid the resonance family, or place it outside the band?"]
    if a_at_c and a_at_b:
        lines.append(
            f"- Same Band-A uniform q evaluated at Band-C notch 4.609375 GHz: "
            f"A2={a_at_c['J_A2']:.3f}, class={a_at_c['classification']['class']}, "
            f"TE1={a_at_c['powers']['TE1_total']:.3f}, wrong0={a_at_c['powers']['TE0_wrong']:.3f}, "
            f"Rtot={a_at_c['powers']['R_tot']:.3f}."
        )
        lines.append(
            f"- Same q at Band-B deepest min 4.909375 GHz: "
            f"A2={a_at_b['J_A2']:.3f}, class={a_at_b['classification']['class']}, "
            f"TE1={a_at_b['powers']['TE1_total']:.3f}, wrong0={a_at_b['powers']['TE0_wrong']:.3f}, "
            f"Rtot={a_at_b['powers']['R_tot']:.3f}."
        )
        # interpret
        if a_at_c["J_A2"] < -5 or a_at_b["J_A2"] < -5:
            lines.append(
                "- **Evidence:** the A-uniform density pattern still hosts deep failures at out-of-band "
                "notch frequencies → Band A success is partly **band placement / avoidance**, not total "
                "elimination of the resonance family for that q."
            )
        elif a_at_c["J_A2"] < a_soft["J_A2"] - 1 or a_at_b["J_A2"] < a_soft["J_A2"] - 1:
            lines.append(
                "- **Evidence:** out-of-band probes are worse than the in-band soft min but not catastrophic "
                "A2≤−5 → partial avoidance; resonances may be weakened or detuned relative to B/C optima."
            )
        else:
            lines.append(
                "- **Evidence:** out-of-band probes remain mild → this q appears to avoid the deep failure "
                "family, not merely push it outside 4.9–5.1 GHz."
            )

    lines += ["", "### 3. Does Band B put the first such resonance at/near the usable-band edge?"]
    if b_deep:
        f0, f1 = DESIGNS["B_stage5"]["band"]
        edge_dist = min(abs(b_deep["f_GHz"] - f0), abs(b_deep["f_GHz"] - f1))
        lines.append(
            f"- Deepest Band-B min at {b_deep['f_GHz']:.6f} GHz is {edge_dist*1e3:.1f} MHz from the "
            f"nearer band edge ({f0}–{f1}); A2={b_deep['J_A2']:.3f} (soft floor, not A2≤−5 swarm)."
        )
        lines.append(
            "- Interpretation depends on class above: if TE1/wrong-port dominated mid-band, the "
            "soft floor is an in-band modal tradeoff; if near-edge, usable BW may be edge-clipped."
        )

    lines += ["", "### 4. Does increasing B move/split resonances or change modal character?"]
    if c_notch and c_swap:
        lines.append(
            f"- Same Band-C q at notch freq: B=0.10 → A2={c_notch['J_A2']:.3f}, "
            f"class={c_notch['classification']['class']}, TE1={c_notch['powers']['TE1_total']:.3f}, "
            f"wrong0={c_notch['powers']['TE0_wrong']:.3f}, Rtot={c_notch['powers']['R_tot']:.3f}."
        )
        lines.append(
            f"- Same q, B=0.05 swap: A2={c_swap['J_A2']:.3f}, class={c_swap['classification']['class']}, "
            f"TE1={c_swap['powers']['TE1_total']:.3f}, wrong0={c_swap['powers']['TE0_wrong']:.3f}, "
            f"Rtot={c_swap['powers']['R_tot']:.3f}."
        )
        lines.append(
            "- Compare also the B=0.05 dither curve vs B=0.10 dither curve in the dither section: "
            "a shifted A2 minimum indicates resonance motion; a class change indicates qualitative "
            "modal redistribution."
        )

    lines += [
        "",
        "### 5. Smallest justified next physical change",
        "",
        "Pending the numerical classes above, priority order for a *small* next test (still no Stage-6):",
        "1. If notches are **A (TE1)** dominated → feed/horn TE1 cutoff or mode-filter geometry probe.",
        "2. If **C (wrong TE0 port)** dominated → circulation routing / gyrotropy (justified higher B) or port-coupling geometry.",
        "3. If **B (reflection)** dominated → matching / density patterning near ports; not denser Stage-6 sampling.",
        "4. If **D mixed** → do not expect a single knob; combine B with a targeted geometry change.",
        "5. Basin/initialization remains relevant for narrow bands (A Stage-5 vs uniform) but is secondary once a physical channel is identified.",
        "",
        "## Artifacts",
        "",
        f"- `{AUTO/'notch_autopsy_results.json'}`",
        f"- `{AUTO/'notch_autopsy_table.csv'}`",
        f"- `{AUTO/'C_B010_Bs5_field_notch_vs_good.png'}`",
        "",
        "Stage-6 remains: **NOT STARTED / AWAITING AUTHORIZATION.**",
        "",
    ]
    path = AUTO / "NOTCH_AUTOPSY_REPORT.md"
    path.write_text("\n".join(lines) + "\n")
    return path


def write_csv(results):
    import csv

    path = AUTO / "notch_autopsy_table.csv"
    fields = [
        "design",
        "b_tesla",
        "f_GHz",
        "role",
        "J_A2",
        "TE0_des",
        "TE0_rev",
        "TE0_wrong",
        "R0",
        "R1",
        "TE1_total",
        "TE1_at_des",
        "P_out",
        "missing",
        "class",
        "IL_avg_dB",
        "isolation_avg_dB",
        "note",
    ]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in results:
            if "error" in r:
                w.writerow(
                    {
                        "design": r.get("design"),
                        "b_tesla": r.get("b_tesla"),
                        "f_GHz": r.get("f_GHz"),
                        "role": r.get("role"),
                        "note": r.get("error"),
                    }
                )
                continue
            p = r["powers"]
            w.writerow(
                {
                    "design": r["design"],
                    "b_tesla": r["b_tesla"],
                    "f_GHz": r["f_GHz"],
                    "role": r["role"],
                    "J_A2": r["J_A2"],
                    "TE0_des": p["TE0_des"],
                    "TE0_rev": p["TE0_rev"],
                    "TE0_wrong": p["TE0_wrong"],
                    "R0": p["R0"],
                    "R1": p["R1"],
                    "TE1_total": p["TE1_total"],
                    "TE1_at_des": p["TE1_at_des"],
                    "P_out": p["P_out"],
                    "missing": p["missing"],
                    "class": r["classification"]["class"],
                    "IL_avg_dB": r["IL_avg_dB"],
                    "isolation_avg_dB": r["isolation_avg_dB"],
                    "note": r.get("note"),
                }
            )
    return path


def plot_dither(results):
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7), constrained_layout=True)
    series = [
        ("B_stage5", 0.05, "B dither @ deepest min"),
        ("C_B010_Bs5", 0.10, "C dither @ A2≤−5 notch"),
        ("C_B010_Bs5", 0.05, "C same-q B=0.05 swap dither"),
        ("A_uniform", 0.05, "A dither @ soft min"),
    ]
    for ax, (design, b, title) in zip(axes.ravel(), series):
        rows = [
            r
            for r in results
            if r.get("design") == design
            and abs(r.get("b_tesla", -1) - b) < 1e-9
            and r.get("dither_parent") is not None
            and "error" not in r
        ]
        # group by parent — take the parent with most points / deepest
        if not rows:
            ax.set_title(title + " (missing)")
            continue
        parent = sorted(rows, key=lambda r: r["J_A2"])[0]["dither_parent"]
        seq = sorted([r for r in rows if abs(r["dither_parent"] - parent) < 1e-9], key=lambda r: r["f_GHz"])
        f = [r["f_GHz"] for r in seq]
        ax.plot(f, [r["J_A2"] for r in seq], "o-", label="A2")
        ax.set_ylabel("A2")
        ax2 = ax.twinx()
        ax2.plot(f, [r["powers"]["TE1_total"] for r in seq], "s--", color="#b85c38", label="TE1")
        ax2.plot(f, [r["powers"]["TE0_wrong"] for r in seq], "^:", color="#2a6f4e", label="wrong0")
        ax2.plot(f, [r["powers"]["R_tot"] for r in seq], "x-.", color="#6b3fa0", label="Rtot")
        ax.set_title(f"{title}\nparent={parent:.6f}")
        ax.grid(True, alpha=0.3)
        ax.set_xlabel("f [GHz]")
        ax2.set_ylabel("power fraction")
    fig.suptitle("Fixed-q frequency dither — modal redistribution")
    fig.savefig(AUTO / "notch_autopsy_dither.png", dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--skip-fields", action="store_true")
    ap.add_argument("--fields-only", action="store_true")
    args = ap.parse_args()
    AUTO.mkdir(parents=True, exist_ok=True)

    sites = pick_sites()
    (AUTO / "notch_autopsy_sites.json").write_text(json.dumps(sites, indent=2) + "\n")
    print("SITES:", json.dumps(sites, indent=2), flush=True)

    results = []
    if not args.fields_only:
        jobs = build_jobs(sites)
        (AUTO / "notch_autopsy_jobs.json").write_text(json.dumps(jobs, indent=2) + "\n")
        t0 = time.time()
        results = run_jobs(jobs, n_workers=args.n_workers)
        elapsed = time.time() - t0
        payload = {
            "stage6_status": "NOT STARTED / AWAITING AUTHORIZATION",
            "n_jobs": len(jobs),
            "n_workers": args.n_workers,
            "elapsed_s": elapsed,
            "sites": sites,
            "rows": results,
        }
        (AUTO / "notch_autopsy_results.json").write_text(json.dumps(payload, indent=2) + "\n")
        write_csv(results)
        plot_dither(results)
        report = write_report(sites, results)
        print(f"wrote {report} elapsed={elapsed/60:.1f} min", flush=True)
        # print compact table to terminal
        print("\n=== MODAL POWER TABLE ===", flush=True)
        print(
            f"{'design':14} {'B':5} {'f':10} {'role':28} {'A2':8} {'des0':6} {'wrong':6} {'Rtot':6} {'TE1':6} {'Pout':6} {'class'}",
            flush=True,
        )
        for r in results:
            if "error" in r:
                continue
            p = r["powers"]
            print(
                f"{r['design']:14} {r['b_tesla']:5.3f} {r['f_GHz']:10.6f} {r['role'][:28]:28} "
                f"{r['J_A2']:8.3f} {p['TE0_des']:6.3f} {p['TE0_wrong']:6.3f} {p['R_tot']:6.3f} "
                f"{p['TE1_total']:6.3f} {p['P_out']:6.3f} {r['classification']['class']}",
                flush=True,
            )

    if not args.skip_fields:
        # deepest C notch vs best-in-band control
        c_sites = sites["C_B010_Bs5"]
        notch_f = c_sites["notches"][0]["f_GHz"]
        good_f = c_sites["controls"][0]["f_GHz"]
        print(f"FIELD MAPS notch={notch_f} good={good_f}", flush=True)
        export_field_maps(
            "C_B010_Bs5",
            [(notch_f, "notch"), (good_f, "good")],
            b_tesla=0.10,
        )

    # Finalize interpretation paragraph with field note
    if results and not args.fields_only:
        # rewrite report once more after fields exist (paths already referenced)
        write_report(sites, results)

    print("STAGE-6: NOT STARTED / AWAITING AUTHORIZATION", flush=True)


if __name__ == "__main__":
    main()
