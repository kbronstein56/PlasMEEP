#!/usr/bin/env python3
"""Offline Objective A alignment audit + A2 candidate ranking on saved spectra.

No FEM. Uses bb_final_audit_*.json and opt_broadband_production_C.json metrics.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"

# Production Objective A weights (frozen)
W_DES, W_REV, W_THR, W_LEAK, W_REFL, W_TE1 = 1.0, 2.0, 1.5, 1.0, 0.25, 1.0
EPS = 1e-6


def load_merged_spectrum():
    """Merge 50 MHz full-band + 25 MHz local midpoints (prefer denser where present)."""
    a50 = json.loads((OUT / "bb_final_audit_50MHz.json").read_text())
    a25 = json.loads((OUT / "bb_final_audit_25MHz_local.json").read_text())
    by_f = {}
    for r in a50["rows"]:
        by_f[round(r["f_GHz"], 6)] = dict(r)
    for r in a25["rows"]:
        by_f[round(r["f_GHz"], 6)] = dict(r)  # overwrite with 25 MHz where available
    rows = [by_f[k] for k in sorted(by_f)]
    # Need T_rev etc — from 50/25 audit we have T_des, TE1, isolation, IL, RL, rev0
    # For seeded 100 MHz nodes may lack des_tot; reconstruct engineering fields
    for r in rows:
        Tdes = float(r.get("T_des") or r.get("des0") or 0.0)
        Trev = float(r.get("rev0") if r.get("rev0") is not None else r.get("T_rev", 0.0) or 0.0)
        # If T_rev missing, infer from isolation: iso_dB = 10 log10(Tdes/Trev)
        if (r.get("rev0") is None and "T_rev" not in r) or Trev == 0.0:
            iso = r.get("isolation_dB")
            if iso is not None and Tdes > 0:
                Trev = Tdes / (10 ** (iso / 10.0))
            elif iso is not None and Tdes <= 0:
                # Tdes tiny: use RL and power rough — leave as nan-safe eps estimate
                Trev = max(Tdes, EPS) / (10 ** (iso / 10.0)) if np.isfinite(iso) else EPS
        R = None
        if r.get("RL_dB") is not None:
            R = 10 ** (-float(r["RL_dB"]) / 10.0)
        r["T_des"] = Tdes
        r["T_rev"] = float(Trev)
        r["R"] = float(R) if R is not None else float("nan")
        r["TE1"] = float(r.get("TE1") or 0.0)
        r["IL_dB"] = float(r["IL_dB"])
        r["isolation_dB"] = float(r["isolation_dB"])
        r["RL_dB"] = float(r["RL_dB"])
        # J_A from file if present
        if r.get("J_A") is not None:
            r["J_A"] = float(r["J_A"])
    return rows


def load_obj_grid_full_metrics():
    """Full per-source-averaged metrics from production rec_best (has T_thr, other)."""
    opt = json.loads((OUT / "opt_broadband_production_C.json").read_text())
    rec = opt["rec_best"]
    rows = []
    for fHz, J, m in zip(rec["freqs_Hz"], rec["J_per_freq"], rec["metrics_per_freq"]):
        rows.append(
            {
                "f_GHz": fHz / 1e9,
                "J_A": float(J),
                "T_des": m["T_des_avg"],
                "T_rev": m["T_rev_avg"],
                "T_thr": m["T_thr_avg"],
                "T_other": m["T_other_avg"],
                "R": m["R_avg"],
                "TE1": m["TE1_avg"],
                "IL_dB": m["IL_avg_dB"],
                "isolation_dB": m["isolation_avg_dB"],
                "RL_dB": m["RL_avg_dB"],
                "P_out": m["P_out_avg"],
            }
        )
    return rows, opt


def decompose_JA(r):
    """J_A = W_des Tdes - W_rev Trev - W_thr Tthr - W_leak other - W_refl R - W_te1 TE1."""
    Tdes = r["T_des"]
    Trev = r["T_rev"]
    Tthr = r.get("T_thr", float("nan"))
    other = r.get("T_other", float("nan"))
    R = r["R"]
    te1 = r["TE1"]
    terms = {
        "reward_des": W_DES * Tdes,
        "pen_rev": -W_REV * Trev,
        "pen_thr": -W_THR * Tthr if np.isfinite(Tthr) else float("nan"),
        "pen_other": -W_LEAK * other if np.isfinite(other) else float("nan"),
        "pen_refl": -W_REFL * R if np.isfinite(R) else float("nan"),
        "pen_te1": -W_TE1 * te1,
    }
    if all(np.isfinite(v) for v in terms.values()):
        terms["J_recon"] = sum(terms.values())
    else:
        # partial: des, rev, refl, te1 only
        terms["J_partial"] = (
            terms["reward_des"]
            + terms["pen_rev"]
            + (terms["pen_refl"] if np.isfinite(terms["pen_refl"]) else 0.0)
            + terms["pen_te1"]
        )
    return terms


# ---------------------------------------------------------------------------
# Candidate objectives (maximize). All differentiable in (Tdes, Trev, ...).
# ---------------------------------------------------------------------------

def J_A(r):
    return float(r["J_A"])


def J_A2_log_ratio(r, eps=EPS):
    """log(Tdes+eps) + log((Tdes+eps)/(Trev+eps)) - penalties on R, TE1, thr, other.

    Physical meaning:
      log(Tdes): rewards useful transmission (IL-like, unbounded below as T→0)
      log(Tdes/Trev): rewards isolation in dB/10*ln(10) units
    Scale: coefficients chosen so a 40 dB IL failure (Tdes~1e-4) is ~ -9.2 from
    log term alone, dominating typical A-scale O(1) scores.
    """
    Td, Tr = r["T_des"] + eps, r["T_rev"] + eps
    R = r["R"] if np.isfinite(r["R"]) else 0.0
    te1 = r["TE1"]
    thr = r.get("T_thr", 0.0) or 0.0
    oth = r.get("T_other", 0.0) or 0.0
    if not np.isfinite(thr):
        thr = 0.0
    if not np.isfinite(oth):
        oth = 0.0
    return float(
        np.log(Td)
        + np.log(Td / Tr)
        - 0.25 * R
        - 1.0 * te1
        - 1.5 * thr
        - 1.0 * oth
    )


def J_A2_log_des_iso(r, eps=EPS):
    """α log(Tdes+eps) + β * 0.1*ln(10) * iso_nat where iso = ln((Td)/(Tr)).

    Equivalent to weighted log-des + isolation in nepers.
    Uses α=1, β=1 → same as J_A2_log_ratio without absolute R/TE1 if absent.
    Keep same leakage penalties as A for TE1/thr/other/R.
    """
    return J_A2_log_ratio(r, eps=eps)


def J_A2_softplus_floor(r, eps=EPS, t_floor=0.1, iso_floor_lin=10.0):
    """Soft constraint form (thresholded) — for comparison only.

    t_floor=0.1 ↔ IL≈10 dB (illustrative, not frozen).
    iso_floor_lin=10 ↔ isolation 10 dB (Tdes/Trev ≥ 10).
    softplus(x)=ln(1+e^x); we maximize -softplus(floor - metric).
    """
    Td, Tr = r["T_des"] + eps, r["T_rev"] + eps
    ratio = Td / Tr
    # maximize: -softplus(t_floor - Tdes) - softplus(log(iso_floor) - log(ratio))
    sp = lambda x: np.log1p(np.exp(np.clip(x, -40, 40)))
    return float(
        -sp(t_floor - r["T_des"])
        - sp(np.log(iso_floor_lin) - np.log(ratio))
        - 0.25 * (r["R"] if np.isfinite(r["R"]) else 0.0)
        - r["TE1"]
    )


def J_A2_hybrid(r, eps=EPS):
    """Preferred candidate: absolute-power A terms replaced by log-des + log-ratio.

    J = log(Tdes+ε) + log((Tdes+ε)/(Trev+ε))
        - W_thr Tthr - W_leak Tother - W_refl R - W_te1 TE1

    Keeps linear penalties on thr/other/R/TE1 (same weights as A) so those
    channels remain discouraged; the pathological Tdes≈Trev≈0 case is hit by
    the log(Tdes) term (→ −∞ as Tdes→0), not escapable by shrinking reverse.
    """
    Td, Tr = r["T_des"] + eps, r["T_rev"] + eps
    thr = r.get("T_thr", 0.0)
    oth = r.get("T_other", 0.0)
    if not np.isfinite(thr):
        thr = 0.0
    if not np.isfinite(oth):
        oth = 0.0
    R = r["R"] if np.isfinite(r["R"]) else 0.0
    return float(
        np.log(Td)
        + np.log(Td / Tr)
        - W_THR * thr
        - W_LEAK * oth
        - W_REFL * R
        - W_TE1 * r["TE1"]
    )


CANDIDATES = {
    "A_current": J_A,
    "A2_log_ratio": J_A2_log_ratio,
    "A2_hybrid": J_A2_hybrid,
    "A2_softplus_floor": J_A2_softplus_floor,
}


def spearman(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3:
        return float("nan")
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def pearson(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def rank_ascending(scores):
    """Rank 1 = worst (lowest score) for maximization objectives."""
    order = np.argsort(scores)  # ascending: worst first
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(scores) + 1)
    return ranks


def analyze():
    OUT.mkdir(parents=True, exist_ok=True)
    obj_rows, opt = load_obj_grid_full_metrics()
    dense = load_merged_spectrum()

    # --- Phase 1: decompose 4.1 and worst-J freqs on obj grid ---
    r41 = next(r for r in obj_rows if abs(r["f_GHz"] - 4.1) < 1e-9)
    # find worst J_A on dense
    dense_ja = [r for r in dense if r.get("J_A") is not None]
    worst_dense = min(dense_ja, key=lambda r: r["J_A"])
    r5625 = next((r for r in dense_ja if abs(r["f_GHz"] - 5.625) < 1e-9), None)

    decomp_41 = decompose_JA(r41)
    # For 5.625 we may lack T_thr/other — use partial + note
    decomp_5625 = decompose_JA(r5625) if r5625 else None

    # Pathological test: synthetic tiny both
    synth = {
        "T_des": 1e-4,
        "T_rev": 1e-3,
        "T_thr": 0.0,
        "T_other": 0.0,
        "R": 0.9,
        "TE1": 0.05,
        "J_A": W_DES * 1e-4 - W_REV * 1e-3 - W_REFL * 0.9 - W_TE1 * 0.05,
    }

    # Rank disagreements on obj grid: IL vs J_A
    ils = [r["IL_dB"] for r in obj_rows]
    isos = [r["isolation_dB"] for r in obj_rows]
    jas = [r["J_A"] for r in obj_rows]
    # Higher IL = worse; lower J_A = worse. Rank IL descending vs J ascending
    rank_il = rank_ascending([-il for il in ils])  # worst IL → rank 1
    rank_ja = rank_ascending(jas)
    disagreements = []
    for i, r in enumerate(obj_rows):
        disagreements.append(
            {
                "f_GHz": r["f_GHz"],
                "J_A": r["J_A"],
                "IL_dB": r["IL_dB"],
                "isolation_dB": r["isolation_dB"],
                "T_des": r["T_des"],
                "T_rev": r["T_rev"],
                "R": r["R"],
                "TE1": r["TE1"],
                "rank_J_A_worst": int(rank_ja[i]),
                "rank_IL_worst": int(rank_il[i]),
                "rank_gap_IL_vs_J": int(rank_il[i] - rank_ja[i]),
                "decomp": decompose_JA(r),
            }
        )
    disagreements.sort(key=lambda d: -abs(d["rank_gap_IL_vs_J"]))

    # --- Phase 3: offline candidate ranking on dense (need T_thr for hybrid — fill 0 if missing) ---
    cand_reports = {}
    for name, fn in CANDIDATES.items():
        scores = []
        for r in dense_ja:
            rr = dict(r)
            if "T_thr" not in rr or rr.get("T_thr") is None:
                rr["T_thr"] = 0.0
            if "T_other" not in rr or rr.get("T_other") is None:
                rr["T_other"] = 0.0
            if name == "A_current":
                scores.append(fn(rr))
            else:
                scores.append(fn(rr))
        scores = np.asarray(scores, float)
        ranks = rank_ascending(scores)
        # map f -> index
        freqs = [r["f_GHz"] for r in dense_ja]
        idx = {round(f, 6): i for i, f in enumerate(freqs)}
        worst10 = [
            {
                "f_GHz": freqs[i],
                "score": float(scores[i]),
                "IL_dB": dense_ja[i]["IL_dB"],
                "isolation_dB": dense_ja[i]["isolation_dB"],
                "T_des": dense_ja[i]["T_des"],
                "T_rev": dense_ja[i]["T_rev"],
            }
            for i in np.argsort(scores)[:10]
        ]

        def rank_of(f):
            i = idx.get(round(f, 6))
            return None if i is None else int(ranks[i])

        # pathology check: 4.1 must be among worst 5 for A2
        r41_rank = rank_of(4.1)
        r5625_rank = rank_of(5.625)
        # correlation with IL (higher IL ↔ lower score → negative pearson with IL)
        ils_d = [r["IL_dB"] for r in dense_ja]
        isos_d = [r["isolation_dB"] for r in dense_ja]
        cand_reports[name] = {
            "worst_10": worst10,
            "rank_4p1_GHz": r41_rank,
            "rank_5p625_GHz": r5625_rank,
            "pearson_score_vs_neg_IL": pearson(scores, [-x for x in ils_d]),
            "pearson_score_vs_isolation": pearson(scores, isos_d),
            "spearman_score_vs_neg_IL": spearman(scores, [-x for x in ils_d]),
            "spearman_score_vs_isolation": spearman(scores, isos_d),
            "n_freq": len(scores),
            "score_at_4p1": float(scores[idx[4.1]]) if 4.1 in idx else None,
            "score_at_5p625": float(scores[idx[5.625]]) if 5.625 in idx else None,
        }

    # Synthetic pathology: Tdes=Trev→0
    synth_scores = {
        name: (fn(synth) if name != "A_current" else synth["J_A"])
        for name, fn in CANDIDATES.items()
    }
    # Good reference point
    good = {
        "T_des": 0.5,
        "T_rev": 0.01,
        "T_thr": 0.05,
        "T_other": 0.02,
        "R": 0.05,
        "TE1": 0.02,
        "J_A": W_DES * 0.5 - W_REV * 0.01 - W_THR * 0.05 - W_LEAK * 0.02 - W_REFL * 0.05 - W_TE1 * 0.02,
    }
    good_scores = {
        name: (fn(good) if name != "A_current" else good["J_A"])
        for name, fn in CANDIDATES.items()
    }

    # Accept/reject
    decisions = {}
    for name, rep in cand_reports.items():
        # Must: 4.1 in worst 5 (rank <= 5), and synth << good
        ok_41 = rep["rank_4p1_GHz"] is not None and rep["rank_4p1_GHz"] <= 5
        ok_path = synth_scores[name] < good_scores[name] - 1.0
        # Correlation with -IL should be positive (higher score ↔ lower IL)
        ok_corr = rep["spearman_score_vs_neg_IL"] > 0.3
        decisions[name] = {
            "accept": bool(ok_41 and ok_path and ok_corr),
            "4p1_in_worst5": ok_41,
            "pathology_synth_worse_than_good": ok_path,
            "corr_with_neg_IL_ok": ok_corr,
            "synth_score": synth_scores[name],
            "good_score": good_scores[name],
        }

    report = {
        "phase1_why_4p1_not_worst_JA": {
            "row_4p1": r41,
            "decomp_4p1": decomp_41,
            "explanation": (
                "J_A = 1·Tdes − 2·Trev − 1.5·Tthr − 1·Tother − 0.25·R − 1·TE1. "
                "At 4.1 GHz Tdes≈3e-4 and Trev is also tiny (isolation −19 dB means "
                "Trev ≈ Tdes/10^(−1.88) ≈ 20×Tdes still tiny absolute). Absolute-power "
                "reverse/through penalties therefore vanish. Remaining score is dominated by "
                "−0.25·R − TE1 ≈ −0.25·(~0.7) − 0.22 ≈ −0.4, not by the 40 dB IL failure. "
                "Both desired and reverse can collapse toward zero without a severe J_A hit."
            ),
            "worst_dense_JA": {"f_GHz": worst_dense["f_GHz"], "J_A": worst_dense["J_A"]},
            "row_5p625": r5625,
            "decomp_5p625": decomp_5625,
            "synthetic_tiny_both": {"point": synth, "J_A": synth["J_A"]},
        },
        "rank_disagreements_obj_grid": disagreements[:10],
        "correlations_obj_grid": {
            "pearson_JA_vs_Tdes": pearson(jas, [r["T_des"] for r in obj_rows]),
            "pearson_JA_vs_neg_IL": pearson(jas, [-r["IL_dB"] for r in obj_rows]),
            "pearson_JA_vs_isolation": pearson(jas, isos),
            "spearman_JA_vs_neg_IL": spearman(jas, [-r["IL_dB"] for r in obj_rows]),
            "spearman_JA_vs_isolation": spearman(jas, isos),
        },
        "candidate_ranking": cand_reports,
        "candidate_decisions": decisions,
        "recommended": "A2_hybrid",
        "eps": EPS,
    }
    (OUT / "BB_OBJECTIVE_ALIGNMENT_OFFLINE.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    rep = analyze()
    print("recommended", rep["recommended"])
    print("decisions", json.dumps(rep["candidate_decisions"], indent=2))
    print("4.1 decomp", json.dumps(rep["phase1_why_4p1_not_worst_JA"]["decomp_4p1"], indent=2))
    print("rank 4.1 / 5.625:")
    for k, v in rep["candidate_ranking"].items():
        print(f"  {k}: rank4.1={v['rank_4p1_GHz']} rank5.625={v['rank_5p625_GHz']} "
              f"spearman_negIL={v['spearman_score_vs_neg_IL']:.3f}")
