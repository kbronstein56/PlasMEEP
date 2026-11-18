#!/usr/bin/env python3
"""Consolidate feasibility-study results into comparison table, notch
diagnostics, tradeoff plots, and CIRCULATOR_FEASIBILITY_STUDY.md.

No new optimizations. Stage-6 NOT STARTED.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
FEAS = OUT / "feasibility"
BOUNDS = (0.05, 1.8)
LO_EPS, HI_EPS = 0.05001, 1.7999


def _fmt(x, nd=4):
    if x is None:
        return "—"
    if isinstance(x, float):
        if math.isnan(x):
            return "—"
        return f"{x:.{nd}f}"
    return str(x)


def _start_label(start) -> str:
    if start is None:
        return "—"
    if isinstance(start, str):
        return start
    if isinstance(start, dict):
        s = str(start.get("start") or "")
        src = str(start.get("source") or "")
        # File warm-starts first (filenames may contain "stage5").
        if s.startswith("file:") or src.endswith("_opt.json"):
            name = Path(src or s.replace("file:", "")).name
            return f"file:{name}"
        if "uniform" in s.lower():
            return "uniform q=1"
        if s in ("stage5", "stage5_q_best") or s.lower() == "stage5_q_best":
            return "Stage-5 q_best"
        if "stage5_q_best" in s.lower():
            return "Stage-5 q_best"
        return (s or src)[:60] or "—"
    return str(start)[:60]


def _bound_counts(q) -> tuple[int, int, int]:
    q = np.asarray(q, float)
    if q.size == 0:
        return (None, None, None)  # type: ignore
    n_lo = int(np.sum(q <= LO_EPS))
    n_hi = int(np.sum(q >= HI_EPS))
    n_int = int(q.size - n_lo - n_hi)
    return n_lo, n_hi, n_int


def _notch_regions(freqs, a2, thr=-5.0):
    """Contiguous regions with A2 <= thr; return list of dicts."""
    freqs = np.asarray(freqs, float)
    a2 = np.asarray(a2, float)
    order = np.argsort(freqs)
    f = freqs[order]
    y = a2[order]
    regions = []
    i = 0
    n = len(f)
    while i < n:
        if y[i] > thr:
            i += 1
            continue
        j = i
        while j + 1 < n and y[j + 1] <= thr:
            j += 1
        sl = slice(i, j + 1)
        k = i + int(np.argmin(y[sl]))
        # width: span of contiguous failing samples; add half-bin edges
        if j > i:
            df = float(np.median(np.diff(f[i : j + 1]))) if j > i else 0.0
        else:
            df = float(np.median(np.diff(f))) if n > 1 else 0.0
        width = float(f[j] - f[i]) + df
        regions.append(
            {
                "f_lo_GHz": float(f[i]),
                "f_hi_GHz": float(f[j]),
                "width_GHz": width,
                "width_MHz": width * 1e3,
                "n_pts": int(j - i + 1),
                "min_A2": float(y[k]),
                "f_min_GHz": float(f[k]),
            }
        )
        i = j + 1
    return regions


def _load_audit_rows(opt_stem: str):
    """Prefer densest SUMMARY+rows for an opt stem."""
    sums = sorted(FEAS.glob(f"audit_{opt_stem}_*_SUMMARY.json"))
    if not sums:
        return None, None
    # prefer finer spacing in name if multiple
    def _df_key(p: Path):
        name = p.name
        for tok in ("1.5625MHz", "3.125MHz", "6.25MHz"):
            if tok in name:
                return {"1.5625MHz": 0, "3.125MHz": 1, "6.25MHz": 2}[tok]
        return 9

    sums = sorted(sums, key=_df_key)
    summary = json.loads(sums[0].read_text())
    rows_path = Path(str(sums[0]).replace("_SUMMARY.json", ".json"))
    rows = None
    if rows_path.exists():
        payload = json.loads(rows_path.read_text())
        rows = payload.get("rows") or payload
    return summary, rows


def _stage5_band_e():
    opt = json.loads((OUT / "opt_broadband_stage5_A2_C.json").read_text())
    q = np.asarray(opt.get("q_best") or [], float)
    n_lo, n_hi, n_int = _bound_counts(q)
    aud = json.loads((OUT / "BB_STAGE5_DENSE_AUDIT_SUMMARY.json").read_text())
    tw = aud.get("true_worst") or {}
    # Prefer local 3.125 MHz for spectrum counts if available
    local = aud.get("dense_3p125MHz_local") or {}
    w = local.get("worst_A2") or tw
    # Load rows for notch + spectrum
    rows = None
    for cand in [
        OUT / "bb_stage5_dense_3p125MHz_local.json",
        OUT / "spectral_bb_stage5_dense_3p125MHz_local.json",
    ]:
        if cand.exists():
            payload = json.loads(cand.read_text())
            rows = payload.get("rows") or payload
            break
    n5 = n7 = None
    if rows:
        n5 = sum(1 for r in rows if r["J_A2"] <= -5)
        n7 = sum(1 for r in rows if r["J_A2"] <= -7)
        # true worst may be finer than 3.125 grid
        worst_row = min(rows, key=lambda r: r["J_A2"])
        if tw.get("J_A2") is not None and float(tw["J_A2"]) < float(worst_row["J_A2"]):
            w = {
                "J_A2": tw.get("J_A2") or tw.get("A2"),
                "f_GHz": tw.get("f_GHz"),
                "IL_dB": tw.get("IL_dB"),
                "isolation_dB": tw.get("isolation_dB"),
            }
        else:
            w = worst_row
    return {
        "tag": "E_4p0_6p0_C6_B0.0500_stage5",
        "band": "E_4p0_6p0",
        "f0_GHz": 4.0,
        "f1_GHz": 6.0,
        "bw_GHz": 2.0,
        "b_tesla": 0.05,
        "symmetry": "C6",
        "n_var": 16,
        "init": "Stage-5 campaign (production)",
        "accepted_nit": opt.get("accepted_nit"),
        "n_eval": opt.get("n_eval"),
        "J_BB": opt.get("J_BB_best"),
        "audit_df_MHz": 3.125,
        "audit_note": "3.125 MHz local + finer true_worst refine",
        "audit_worst_A2": w.get("J_A2"),
        "audit_worst_f": w.get("f_GHz"),
        "audit_IL": w.get("IL_dB"),
        "audit_iso": w.get("isolation_dB"),
        "audit_n_le5": n5 if n5 is not None else local.get("n_A2_le_m5"),
        "audit_n_le7": n7 if n7 is not None else local.get("n_A2_le_m7"),
        "n_lo": n_lo,
        "n_hi": n_hi,
        "n_int": n_int,
        "elapsed_h": (opt.get("elapsed_s") or 0) / 3600.0,
        "optimizer_message": opt.get("optimizer_message"),
        "status": "complete",
        "rows": rows,
        "opt_worst_A2": (opt.get("rec_best") or {}).get("worst_J"),
    }


def load_all_cases():
    cases = []
    for opt_path in sorted(FEAS.glob("*_opt.json")):
        d = json.loads(opt_path.read_text())
        if d.get("status") != "complete":
            continue
        q = np.asarray(d.get("q_best") or [], float)
        n_lo, n_hi, n_int = _bound_counts(q)
        # also trust stored near-bound flags if present
        if d.get("q_near_lo") is not None:
            n_lo = int(d["q_near_lo"])
        if d.get("q_near_hi") is not None:
            n_hi = int(d["q_near_hi"])
            n_int = int(d.get("n_var") or q.size) - n_lo - n_hi
        summary, rows = _load_audit_rows(opt_path.stem)
        w = (summary or {}).get("worst_A2") or {}
        # parse df from summary or filename
        df = (summary or {}).get("df_MHz")
        if df is None and summary is not None:
            for tok in ("1.5625", "3.125", "6.25"):
                if tok in str(summary.get("opt") or ""):
                    df = float(tok)
                    break
        case = {
            "tag": opt_path.stem.replace("_opt", ""),
            "band": d.get("band"),
            "f0_GHz": d.get("f0_GHz"),
            "f1_GHz": d.get("f1_GHz"),
            "bw_GHz": float(d["f1_GHz"]) - float(d["f0_GHz"]) if d.get("f0_GHz") is not None else None,
            "b_tesla": d.get("b_tesla"),
            "symmetry": d.get("symmetry"),
            "n_var": d.get("n_var"),
            "init": _start_label(d.get("start")),
            "accepted_nit": d.get("accepted_nit"),
            "n_eval": d.get("n_eval"),
            "J_BB": d.get("J_BB_best"),
            "opt_worst_A2": (d.get("rec_best") or {}).get("worst_J"),
            "audit_df_MHz": df or (summary or {}).get("df_MHz"),
            "audit_worst_A2": w.get("J_A2"),
            "audit_worst_f": w.get("f_GHz"),
            "audit_IL": w.get("IL_dB"),
            "audit_iso": w.get("isolation_dB"),
            "audit_n_le5": (summary or {}).get("n_A2_le_m5"),
            "audit_n_le7": (summary or {}).get("n_A2_le_m7"),
            "n_lo": n_lo,
            "n_hi": n_hi,
            "n_int": n_int,
            "elapsed_h": (d.get("elapsed_s") or 0) / 3600.0,
            "optimizer_message": d.get("optimizer_message"),
            "status": d.get("status"),
            "rows": rows,
        }
        if rows:
            freqs = [r["f_GHz"] for r in rows]
            a2 = [r["J_A2"] for r in rows]
            case["notches_le5"] = _notch_regions(freqs, a2, -5.0)
            case["notches_le7"] = _notch_regions(freqs, a2, -7.0)
            # opt-grid vs audit: compare opt worst to audit worst
            ow = case.get("opt_worst_A2")
            aw = case.get("audit_worst_A2")
            if ow is not None and aw is not None:
                case["hidden_gap"] = float(aw) - float(ow)  # negative => audit worse
                case["hidden_failure"] = bool(aw < ow - 0.5)  # audit meaningfully worse
        cases.append(case)

    cases.append(_stage5_band_e())
    # notches for E if rows present
    for c in cases:
        if c["tag"].startswith("E_") and c.get("rows") and "notches_le5" not in c:
            freqs = [r["f_GHz"] for r in c["rows"]]
            a2 = [r["J_A2"] for r in c["rows"]]
            c["notches_le5"] = _notch_regions(freqs, a2, -5.0)
            c["notches_le7"] = _notch_regions(freqs, a2, -7.0)
            ow = c.get("opt_worst_A2")
            aw = c.get("audit_worst_A2")
            if ow is not None and aw is not None:
                c["hidden_gap"] = float(aw) - float(ow)
                c["hidden_failure"] = bool(aw < ow - 0.5)
    return cases


def write_table(cases):
    keys = [
        ("tag", "case/tag"),
        ("band_str", "band / BW"),
        ("b_tesla", "B [T]"),
        ("sym_str", "sym / n_var"),
        ("init", "init"),
        ("nit_str", "nit / n_eval"),
        ("J_BB", "J_BB"),
        ("audit_df_MHz", "audit Δf [MHz]"),
        ("audit_worst_A2", "audit worst A2"),
        ("audit_worst_f", "f_worst [GHz]"),
        ("audit_IL", "IL [dB]"),
        ("audit_iso", "iso [dB]"),
        ("audit_n_le5", "n(A2≤−5)"),
        ("audit_n_le7", "n(A2≤−7)"),
        ("bound_str", "q lo/hi/int"),
    ]
    enriched = []
    for c in cases:
        e = dict(c)
        e["band_str"] = f"{c.get('f0_GHz')}–{c.get('f1_GHz')} ({_fmt(c.get('bw_GHz'),3)} GHz)"
        e["sym_str"] = f"{c.get('symmetry')}/{c.get('n_var')}"
        e["nit_str"] = f"{c.get('accepted_nit')}/{c.get('n_eval')}"
        e["bound_str"] = f"{c.get('n_lo')}/{c.get('n_hi')}/{c.get('n_int')}"
        enriched.append(e)

    # JSON (drop rows to keep compact)
    compact = []
    for e in enriched:
        row = {k: e.get(k) for k, _ in keys}
        row.update(
            {
                "bw_GHz": e.get("bw_GHz"),
                "f0_GHz": e.get("f0_GHz"),
                "f1_GHz": e.get("f1_GHz"),
                "symmetry": e.get("symmetry"),
                "n_var": e.get("n_var"),
                "accepted_nit": e.get("accepted_nit"),
                "n_eval": e.get("n_eval"),
                "n_lo": e.get("n_lo"),
                "n_hi": e.get("n_hi"),
                "n_int": e.get("n_int"),
                "opt_worst_A2": e.get("opt_worst_A2"),
                "hidden_gap": e.get("hidden_gap"),
                "hidden_failure": e.get("hidden_failure"),
                "elapsed_h": e.get("elapsed_h"),
                "notches_le5": e.get("notches_le5"),
                "notches_le7": e.get("notches_le7"),
                "optimizer_message": e.get("optimizer_message"),
            }
        )
        compact.append(row)
    (FEAS / "FEASIBILITY_COMPARISON_TABLE.json").write_text(json.dumps(compact, indent=2) + "\n")

    # Markdown
    hdr = [lab for _, lab in keys]
    lines = [
        "# Feasibility comparison table (all completed cases)",
        "",
        "Stage-6: **NOT STARTED / AWAITING AUTHORIZATION**",
        "",
        "Dense-audit metrics are authoritative for notch severity; optimizer J_BB is sampled on the design grid only.",
        "",
        "| " + " | ".join(hdr) + " |",
        "| " + " | ".join(["---"] * len(hdr)) + " |",
    ]
    for e in sorted(enriched, key=lambda x: (x.get("bw_GHz") or 99, str(x.get("tag")))):
        vals = []
        for k, _ in keys:
            v = e.get(k)
            if isinstance(v, float):
                vals.append(_fmt(v, 4 if abs(v) < 100 else 2))
            else:
                vals.append("—" if v is None else str(v))
        lines.append("| " + " | ".join(vals) + " |")
    (FEAS / "FEASIBILITY_COMPARISON_TABLE.md").write_text("\n".join(lines) + "\n")
    return compact


def write_notch_report(cases):
    lines = [
        "# Notch character (dense audit, A2 ≤ −5 regions)",
        "",
        "Stage-6: NOT STARTED",
        "",
    ]
    for c in sorted(cases, key=lambda x: (x.get("bw_GHz") or 99, str(x.get("tag")))):
        notches = c.get("notches_le5") or []
        lines.append(f"## {c['tag']}")
        lines.append(
            f"- band {c.get('f0_GHz')}–{c.get('f1_GHz')} GHz | B={c.get('b_tesla')} | "
            f"{c.get('symmetry')}/{c.get('n_var')} | init={c.get('init')}"
        )
        lines.append(
            f"- audit worst A2={_fmt(c.get('audit_worst_A2'))} @ {_fmt(c.get('audit_worst_f'))} GHz; "
            f"IL={_fmt(c.get('audit_IL'),2)} dB; iso={_fmt(c.get('audit_iso'),2)} dB"
        )
        lines.append(
            f"- hidden vs opt-grid: gap={_fmt(c.get('hidden_gap'))} "
            f"(audit−opt; negative ⇒ denser audit worse); "
            f"hidden_failure={c.get('hidden_failure')}"
        )
        if not notches:
            lines.append("- **no** contiguous A2≤−5 regions on the dense audit grid")
        else:
            lines.append(f"- **{len(notches)}** A2≤−5 region(s):")
            for n in notches:
                lines.append(
                    f"  - {_fmt(n['f_lo_GHz'],4)}–{_fmt(n['f_hi_GHz'],4)} GHz "
                    f"(~{_fmt(n['width_MHz'],1)} MHz, {n['n_pts']} pts); "
                    f"min A2={_fmt(n['min_A2'])} @ {_fmt(n['f_min_GHz'],4)} GHz"
                )
        lines.append("")
    (FEAS / "FEASIBILITY_NOTCH_REPORT.md").write_text("\n".join(lines) + "\n")


def make_plots(cases):
    # 1) worst A2 vs bandwidth for C6 B=0.05 comparable cases
    # Prefer Stage-5-started ladder; include uniform A and mark separately; exclude Auni collapse as non-comparable warm-start failure for the ladder floor
    ladder_tags_preferred = {
        "A_4p9_5p1_C6_B0.0500_uniform",
        "A_4p9_5p1_C6_B0.0500_stage5",
        "B_4p75_5p25_C6_B0.0500_stage5",
        "C_4p5_5p5_C6_B0.0500_stage5",
        "D_4p25_5p75_C6_B0.0500_Bs5",
        "E_4p0_6p0_C6_B0.0500_stage5",
    }
    ladder = [c for c in cases if c["tag"] in ladder_tags_preferred]
    ladder = sorted(ladder, key=lambda c: c.get("bw_GHz") or 99)

    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    for c in ladder:
        ax.scatter([c["bw_GHz"]], [c["audit_worst_A2"]], s=70, zorder=3)
        ax.annotate(
            c["tag"].split("_C6")[0],
            (c["bw_GHz"], c["audit_worst_A2"]),
            textcoords="offset points",
            xytext=(6, 4),
            fontsize=7,
        )
    if len(ladder) >= 2:
        xs = [c["bw_GHz"] for c in ladder]
        ys = [c["audit_worst_A2"] for c in ladder]
        ax.plot(xs, ys, "-", color="#1f4e79", alpha=0.5, zorder=1)
    ax.axhline(-5, color="#b85c38", ls="--", lw=1, label="A2=−5")
    ax.set_xlabel("Required bandwidth [GHz]")
    ax.set_ylabel("Dense-audit worst A2")
    ax.set_title("C6, B=0.05 T — worst dense A2 vs bandwidth")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.savefig(FEAS / "feasibility_worstA2_vs_bandwidth.png", dpi=150)
    plt.close(fig)

    # 2) A2(f) overlays for A/B/C/D/E C6 B=0.05 primary designs
    overlay_order = [
        ("A_4p9_5p1_C6_B0.0500_uniform", "A 0.2 GHz (uniform)"),
        ("B_4p75_5p25_C6_B0.0500_stage5", "B 0.5 GHz (Stage-5 start)"),
        ("C_4p5_5p5_C6_B0.0500_stage5", "C 1.0 GHz (Stage-5 start)"),
        ("D_4p25_5p75_C6_B0.0500_Bs5", "D 1.5 GHz (B warm start)"),
        ("E_4p0_6p0_C6_B0.0500_stage5", "E 2.0 GHz (Stage-5)"),
    ]
    by_tag = {c["tag"]: c for c in cases}
    fig, ax = plt.subplots(figsize=(10, 5.2), constrained_layout=True)
    colors = ["#2a6f4e", "#1f4e79", "#b85c38", "#6b3fa0", "#7a1f1f"]
    for (tag, lab), col in zip(overlay_order, colors):
        c = by_tag.get(tag)
        if not c or not c.get("rows"):
            continue
        rows = sorted(c["rows"], key=lambda r: r["f_GHz"])
        f = [r["f_GHz"] for r in rows]
        y = [r["J_A2"] for r in rows]
        ax.plot(f, y, lw=1.2, color=col, label=lab, alpha=0.9)
    ax.axhline(-5, color="k", ls="--", lw=0.8, alpha=0.6)
    ax.axhline(-7, color="k", ls=":", lw=0.8, alpha=0.5)
    ax.set_xlim(4.0, 6.0)
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel("A2")
    ax.set_title("Dense-audit A2(f) overlays — C6 B=0.05 T ladder")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")
    fig.savefig(FEAS / "feasibility_A2_overlays_ABCDE.png", dpi=150)
    plt.close(fig)

    # 3) B-field comparison spectra on Band C
    fig, ax = plt.subplots(figsize=(9, 4.8), constrained_layout=True)
    for tag, lab, col in [
        ("C_4p5_5p5_C6_B0.0500_stage5", "B=0.05 T (Stage-5 start)", "#1f4e79"),
        ("C_4p5_5p5_C6_B0.1000_stage5", "B=0.10 T (Stage-5 start)", "#b85c38"),
        ("C_4p5_5p5_C6_B0.1000_Bs5", "B=0.10 T (Band-B warm start)", "#2a6f4e"),
    ]:
        c = by_tag.get(tag)
        if not c or not c.get("rows"):
            continue
        rows = sorted(c["rows"], key=lambda r: r["f_GHz"])
        ax.plot([r["f_GHz"] for r in rows], [r["J_A2"] for r in rows], lw=1.2, color=col, label=lab)
    ax.axhline(-5, color="k", ls="--", lw=0.8, alpha=0.6)
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel("A2")
    ax.set_title("Band C — B-field dense-audit A2(f)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.savefig(FEAS / "feasibility_A2_overlays_Bfield_C.png", dpi=150)
    plt.close(fig)

    # 4) DOF C3 vs C6 on Band B
    fig, ax = plt.subplots(figsize=(9, 4.8), constrained_layout=True)
    for tag, lab, col in [
        ("B_4p75_5p25_C6_B0.0500_stage5", "C6/16 Stage-5 start", "#1f4e79"),
        ("B_4p75_5p25_C3_B0.0500_stage5", "C3/31 Stage-5 start", "#b85c38"),
        ("B_4p75_5p25_C3_B0.0500_Bs5", "C3/31 from C6 Band-B", "#2a6f4e"),
    ]:
        c = by_tag.get(tag)
        if not c or not c.get("rows"):
            continue
        rows = sorted(c["rows"], key=lambda r: r["f_GHz"])
        ax.plot([r["f_GHz"] for r in rows], [r["J_A2"] for r in rows], lw=1.2, color=col, label=lab)
    ax.axhline(-5, color="k", ls="--", lw=0.8, alpha=0.6)
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel("A2")
    ax.set_title("Band B — C6 vs C3 dense-audit A2(f)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.savefig(FEAS / "feasibility_A2_overlays_DOF_B.png", dpi=150)
    plt.close(fig)

    # Keep summarize-style tradeoff panels
    fig, axes = plt.subplots(2, 2, figsize=(9, 7), constrained_layout=True)
    # best-per-bandwidth among C6/0.05 excluding Auni
    best = {}
    for c in cases:
        if c.get("symmetry") != "C6":
            continue
        if abs(float(c.get("b_tesla") or 0) - 0.05) > 1e-9:
            continue
        if "Auni" in c["tag"]:
            continue
        bw = c.get("bw_GHz")
        if bw is None or c.get("audit_worst_A2") is None:
            continue
        prev = best.get(bw)
        if prev is None or c["audit_worst_A2"] > prev["audit_worst_A2"]:
            best[bw] = c
    bw_rows = [best[k] for k in sorted(best)]
    xs = [r["bw_GHz"] for r in bw_rows]
    for ax, key, title in [
        (axes[0, 0], "audit_worst_A2", "Worst A2 (audit)"),
        (axes[0, 1], "audit_IL", "IL at worst A2 [dB]"),
        (axes[1, 0], "audit_iso", "Isolation at worst A2 [dB]"),
        (axes[1, 1], "J_BB", "Optimized J_BB"),
    ]:
        ys = [r.get(key) for r in bw_rows]
        ax.plot(xs, ys, "o-", color="#1f4e79")
        ax.set_xlabel("Required bandwidth [GHz]")
        ax.set_ylabel(title)
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
    fig.suptitle("Best-per-bandwidth C6 B=0.05 T (Auni excluded)")
    fig.savefig(FEAS / "feasibility_tradeoff_bandwidth.png", dpi=150)
    plt.close(fig)


def write_study_md(cases, compact):
    by_tag = {c["tag"]: c for c in cases}

    def g(tag, key, nd=4):
        c = by_tag.get(tag) or {}
        return _fmt(c.get(key), nd)

    # Bandwidth ladder narrative numbers
    a_u = by_tag.get("A_4p9_5p1_C6_B0.0500_uniform")
    a_s = by_tag.get("A_4p9_5p1_C6_B0.0500_stage5")
    b_s = by_tag.get("B_4p75_5p25_C6_B0.0500_stage5")
    b_a = by_tag.get("B_4p75_5p25_C6_B0.0500_Auni")
    c05 = by_tag.get("C_4p5_5p5_C6_B0.0500_stage5")
    c10s = by_tag.get("C_4p5_5p5_C6_B0.1000_stage5")
    c10b = by_tag.get("C_4p5_5p5_C6_B0.1000_Bs5")
    d = by_tag.get("D_4p25_5p75_C6_B0.0500_Bs5")
    e = by_tag.get("E_4p0_6p0_C6_B0.0500_stage5")
    c3s = by_tag.get("B_4p75_5p25_C3_B0.0500_stage5")
    c3b = by_tag.get("B_4p75_5p25_C3_B0.0500_Bs5")

    # Classify bandwidth threshold
    # robust: audit A2 > -5 and n_le5==0
    # transitional: around -5
    # notch-limited: << -5 with many notches
    def floor_label(c):
        if not c or c.get("audit_worst_A2") is None:
            return "missing"
        a2 = c["audit_worst_A2"]
        n5 = c.get("audit_n_le5") or 0
        if a2 > -1.0 and n5 == 0:
            return "robustly feasible (tested)"
        if a2 > -5.2 and n5 <= 2:
            return "borderline / soft floor near A2≈−5"
        return "notch-limited (dense audit)"

    # B-field conclusion
    b_concl = (
        "Doubling B from 0.05→0.10 T on Band C **improves** the dense-audit floor for the better start "
        f"(0.05: A2={g('C_4p5_5p5_C6_B0.0500_stage5','audit_worst_A2')}, n≤−5={g('C_4p5_5p5_C6_B0.0500_stage5','audit_n_le5',0)}; "
        f"0.10 Band-B warm start: A2={g('C_4p5_5p5_C6_B0.1000_Bs5','audit_worst_A2')}, n≤−5={g('C_4p5_5p5_C6_B0.1000_Bs5','audit_n_le5',0)}), "
        "but **does not eliminate** A2≤−5 structure — notches remain, with changed depth/count/location. "
        "The Stage-5-start 0.10 T run remains deeply notched "
        f"(A2={g('C_4p5_5p5_C6_B0.1000_stage5','audit_worst_A2')}). "
        "Conclusion: increased magnetization helps the Pareto tradeoff but is **not sufficient alone** for clean 1 GHz."
    )

    # DOF conclusion
    dof_concl = (
        f"On Band B, C6/16 Stage-5 start reaches audit A2={g('B_4p75_5p25_C6_B0.0500_stage5','audit_worst_A2')} "
        f"(n≤−5={g('B_4p75_5p25_C6_B0.0500_stage5','audit_n_le5',0)}). "
        f"C3/31 from Stage-5: A2={g('B_4p75_5p25_C3_B0.0500_stage5','audit_worst_A2')} "
        f"(n≤−5={g('B_4p75_5p25_C3_B0.0500_stage5','audit_n_le5',0)}, bounds lo/hi="
        f"{g('B_4p75_5p25_C3_B0.0500_stage5','n_lo',0)}/{g('B_4p75_5p25_C3_B0.0500_stage5','n_hi',0)}). "
        f"C3 from C6 Band-B: A2={g('B_4p75_5p25_C3_B0.0500_Bs5','audit_worst_A2')}. "
        "Extra DOFs do **not** materially beat the C6 soft floor; one C3 start saturates many bounds without better audit performance. "
        "Evidence argues **against** C6/16 as the primary 0.5–2 GHz limiter."
    )

    # Next diagnostics (cheap)
    next_diag = [
        "1. **Narrow-band physics audit of a Band-A success vs Band-C failure:** at the Band-A uniform q_best and Band-C 0.05 q_best, extract port S-parameters + field energy / modal indicators on a fine local grid around each A2≤−5 notch (forward-only). Goal: decide whether failures are transmission zeros / trapped resonances of the fixed 91-element+horn geometry.",
        "2. **Single-frequency continuous minimax probe (no new multi-day opt):** on Band C at B=0.05, take the densest audit minima and run a few dozen forward evaluations with ±1–2 MHz dither around each notch under fixed q — confirm notch width/depth are physical resonances, not FEM/noise artifacts.",
        "3. **Only if (1–2) show movable zeros:** a *small* justified-B ladder on Band B (already nearly clean at 0.05) at one additional documented B between 0.05 and 0.10, forward+short reopt, to map whether the soft A2≈−5 floor lifts with gyrotropy before any Stage-6.",
        "4. **Do not** launch Stage-6 on 4–6 GHz with the current formulation: denser sampling already failed through 6.25 MHz / 321 freqs; feasibility ladder shows collapse by ~1 GHz required BW.",
    ]

    lines = []
    lines += [
        "# Circulator feasibility study",
        "",
        "**Status:** COMPLETE (queued ladder + B/DOF probes finished; no further expensive opts)",
        "**Stage-6:** **NOT STARTED / AWAITING AUTHORIZATION**",
        "**Physics package (frozen):** FULL91-C, magnetized plasma, C6/16 unless noted, density bounds [0.05, 1.8],",
        "A2_hybrid τ=0.05 α=0.25, validated fast kernel, n_workers=12, OMP/OpenBLAS/MKL=1.",
        "",
        "Companion artifacts:",
        "- `CIRCULATOR_FEASIBILITY_BASELINE.md`",
        "- `CIRCULATOR_PHYSICAL_LIMITS_CONTEXT.md`",
        "- `feasibility/FEASIBILITY_COMPARISON_TABLE.md` (+ `.json`)",
        "- `feasibility/FEASIBILITY_NOTCH_REPORT.md`",
        "- `feasibility/feasibility_worstA2_vs_bandwidth.png`",
        "- `feasibility/feasibility_A2_overlays_ABCDE.png`",
        "- `feasibility/feasibility_A2_overlays_Bfield_C.png`",
        "- `feasibility/feasibility_A2_overlays_DOF_B.png`",
        "",
        "---",
        "",
        "## Executive answer",
        "",
        "**What limits broadband performance?** Evidence from the completed ladder points primarily to a",
        "**required-bandwidth / resonant-geometry limit** of the present 91-element + C6 + B=0.05 T package,",
        "with a secondary, real but incomplete assist from magnetic bias. C6→C3 DOF expansion and density-bound",
        "saturation are **not** supported as the dominant cause. Denser objective sampling alone is already falsified",
        "by Stages 2–5.",
        "",
        f"- **~0.2 GHz (Band A):** convincingly achievable — uniform start audit A2≈{g('A_4p9_5p1_C6_B0.0500_uniform','audit_worst_A2')}, IL≈{g('A_4p9_5p1_C6_B0.0500_uniform','audit_IL',1)} dB, iso≈{g('A_4p9_5p1_C6_B0.0500_uniform','audit_iso',1)} dB, zero A2≤−5.",
        f"- **~0.5 GHz (Band B):** transitional soft floor — best C6 audit A2≈{g('B_4p75_5p25_C6_B0.0500_stage5','audit_worst_A2')} with n(A2≤−5)={g('B_4p75_5p25_C6_B0.0500_stage5','audit_n_le5',0)}; IL/iso still poor at the worst site.",
        f"- **≥1.0 GHz (Bands C/D/E):** not demonstrated achievable — deep narrow notches return (C audit A2≈{g('C_4p5_5p5_C6_B0.0500_stage5','audit_worst_A2')}; E true-worst A2≈{g('E_4p0_6p0_C6_B0.0500_stage5','audit_worst_A2')}"
        + (f"; D audit A2≈{g('D_4p25_5p75_C6_B0.0500_Bs5','audit_worst_A2')}" if d else "")
        + ").",
        "",
        "**Stage-6 on the current 4–6 GHz formulation is not justified by this study.** Next effort should change",
        "the scientific ask (narrower band and/or higher justified B / geometry), not add another dense full-band L-BFGS stage.",
        "",
        "---",
        "",
        "## Master comparison table",
        "",
        "See `feasibility/FEASIBILITY_COMPARISON_TABLE.md` for the full per-case table (band, B, symmetry/n_var,",
        "init, nit/n_eval, J_BB, audit spacing, worst A2/f/IL/iso, n≤−5/≤−7, bound counts).",
        "",
        "### Snapshot (dense-audit authoritative)",
        "",
        "| case | BW [GHz] | B | DOF | init | J_BB | audit A2 | f_worst | IL | iso | n≤−5 | n≤−7 | lo/hi |",
        "|---|---:|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]

    for c in sorted(cases, key=lambda x: (x.get("bw_GHz") or 99, x["tag"])):
        lines.append(
            f"| {c['tag']} | {_fmt(c.get('bw_GHz'),3)} | {_fmt(c.get('b_tesla'),3)} | "
            f"{c.get('symmetry')}/{c.get('n_var')} | {c.get('init')} | {_fmt(c.get('J_BB'))} | "
            f"{_fmt(c.get('audit_worst_A2'))} | {_fmt(c.get('audit_worst_f'),4)} | "
            f"{_fmt(c.get('audit_IL'),1)} | {_fmt(c.get('audit_iso'),1)} | "
            f"{_fmt(c.get('audit_n_le5'),0)} | {_fmt(c.get('audit_n_le7'),0)} | "
            f"{c.get('n_lo')}/{c.get('n_hi')} |"
        )

    lines += [
        "",
        "---",
        "",
        "## 1. Bandwidth scaling (A→B→C→D→E)",
        "",
        f"| Band | classification | notes |",
        f"|---|---|---|",
        f"| A 0.2 GHz uniform | {floor_label(a_u)} | best positive control |",
        f"| A 0.2 GHz Stage-5 start | {floor_label(a_s)} | bad basin on a solvable band |",
        f"| B 0.5 GHz Stage-5 start | {floor_label(b_s)} | soft floor; 1.5625 MHz refine agreed |",
        f"| B 0.5 GHz A-uniform start | {floor_label(b_a)} | same BW, bad basin → collapse |",
        f"| C 1.0 GHz @0.05 | {floor_label(c05)} | many narrow deep notches |",
    ]
    if d:
        lines.append(f"| D 1.5 GHz @0.05 | {floor_label(d)} | Band-B warm start |")
    lines.append(f"| E 2.0 GHz Stage-5 | {floor_label(e)} | production full-band |")
    lines += [
        "",
        "**Threshold:** performance changes from *robustly feasible* at **~0.2 GHz** to a *soft ~A2≈−5 floor*",
        "at **~0.5 GHz**, then to *deep multi-notch failure* by **~1.0 GHz** and remains failed at 1.5–2.0 GHz.",
        "The A-uniform→B collapse at fixed physics shows bandwidth pressure is not only an optimizer quirk.",
        "",
        "Plot: `feasibility/feasibility_worstA2_vs_bandwidth.png`, overlays `feasibility_A2_overlays_ABCDE.png`.",
        "",
        "---",
        "",
        "## 2. B-field sensitivity (Band C, dense audits)",
        "",
        b_concl,
        "",
        "Plot: `feasibility/feasibility_A2_overlays_Bfield_C.png`.",
        "",
        "---",
        "",
        "## 3. Design DOF / symmetry (Band B, C3 vs C6)",
        "",
        dof_concl,
        "",
        "Plot: `feasibility/feasibility_A2_overlays_DOF_B.png`.",
        "",
        "---",
        "",
        "## 4. Notch character",
        "",
        "Full region list: `feasibility/FEASIBILITY_NOTCH_REPORT.md`.",
        "",
        "Pattern across bandwidth:",
        "- Narrow-band success (A uniform) shows **no** A2≤−5 regions.",
        "- Borderline Band B shows at most shallow/touching structure near A2≈−5 without a swarm of deep notches.",
        "- Wider bands (C/E and D if failed) add **more** deep, still **narrow** notches distributed across the band",
        "  rather than one broad stopband — consistent with migrating/newly admitted resonant zeros.",
        "",
        "---",
        "",
        "## 5. Optimization-grid artifact vs design limitation",
        "",
    ]

    for c in sorted(cases, key=lambda x: x["tag"]):
        if c.get("hidden_failure") is None:
            continue
        flag = "HIDDEN FAIL" if c.get("hidden_failure") else "grid≈audit"
        lines.append(
            f"- `{c['tag']}`: {flag}; opt_worst={_fmt(c.get('opt_worst_A2'))}, "
            f"audit_worst={_fmt(c.get('audit_worst_A2'))}, gap={_fmt(c.get('hidden_gap'))}"
        )

    lines += [
        "",
        "Interpretation: Stage-5-style full-band designs remain vulnerable to inter-node notches (H1 partial),",
        "but the ladder also contains cases where dense audits **agree** with a poor floor (C, E) and a case",
        "where a different init solves a narrow band the Stage-5 basin cannot (A). Together this means sampling",
        "artifacts are real **and** insufficient to explain the bandwidth collapse.",
        "",
        "---",
        "",
        "## 6. Hypotheses (H1–H6)",
        "",
        "| H | Verdict | Confidence |",
        "|---|---|---|",
        "| H1 sampling/opt only | **Weak as sole explanation** — denser grids failed; A works with same optimizer; Auni→B collapses | medium-high |",
        "| H2 limited DOF/C6 | **Not supported as primary** — C3 ≉ better on Band B | medium |",
        "| H3 limited nonreciprocity | **Partial** — B=0.10 helps but leaves notches on 1 GHz | medium |",
        "| H4 density bounds | **Unlikely primary** — best A/B solutions largely interior; saturated C3 did not win | medium-high |",
        "| H5 geometry/modal | **Leading** for ≥1 GHz ask under current passive platform | medium |",
        "| H6 combination | **Likely** — bandwidth×geometry primary, B secondary, basin tertiary | medium-high |",
        "",
        "---",
        "",
        "## 7. Smallest next diagnostics (no Stage-6)",
        "",
    ]
    lines += next_diag
    lines += [
        "",
        "---",
        "",
        "## 8. Final recommendations",
        "",
        "1. **Convincing BW of current package:** ~**200 MHz** centered near 5 GHz; ~**500 MHz** only at a weak soft floor.",
        "2. **Degradation with BW:** yes — abrupt between 0.2 and 1.0 GHz.",
        "3. **4–6 GHz:** not demonstrated achievable in the tested design space.",
        "4. **Most important constraint:** required **bandwidth** interacting with **fixed geometry / modal structure**; B helps secondarily; DOF/bounds/optimizer secondary.",
        "5. **One change most likely to help:** **narrow the target band** (or retarget usable sub-bands) **or** raise B within documented hardware limits **and** reconsider geometry — not denser Stage-6 sampling.",
        "6. **Next major effort:** **E (narrow target)** and/or **C (justified higher B)** and/or **D (geometry)** — **not A (Stage-6 on current 4–6 GHz)**.",
        "",
        "Stage-6 remains: **NOT STARTED / AWAITING AUTHORIZATION.**",
        "",
    ]
    path = OUT / "CIRCULATOR_FEASIBILITY_STUDY.md"
    path.write_text("\n".join(lines) + "\n")
    # also copy a pointer into feasibility/
    (FEAS / "CIRCULATOR_FEASIBILITY_STUDY.md").write_text(
        "Canonical copy: `../CIRCULATOR_FEASIBILITY_STUDY.md` (updated by feasibility_consolidate.py).\n"
    )
    return path


def print_conclusions(cases):
    by = {c["tag"]: c for c in cases}
    print("=" * 72)
    print("FEASIBILITY STUDY — KEY NUMERICAL CONCLUSIONS")
    print("Stage-6: NOT STARTED / AWAITING AUTHORIZATION")
    print("=" * 72)
    order = [
        "A_4p9_5p1_C6_B0.0500_uniform",
        "A_4p9_5p1_C6_B0.0500_stage5",
        "B_4p75_5p25_C6_B0.0500_stage5",
        "B_4p75_5p25_C6_B0.0500_Auni",
        "B_4p75_5p25_C3_B0.0500_stage5",
        "B_4p75_5p25_C3_B0.0500_Bs5",
        "C_4p5_5p5_C6_B0.0500_stage5",
        "C_4p5_5p5_C6_B0.1000_stage5",
        "C_4p5_5p5_C6_B0.1000_Bs5",
        "D_4p25_5p75_C6_B0.0500_Bs5",
        "E_4p0_6p0_C6_B0.0500_stage5",
    ]
    print(
        f"{'tag':48} {'BW':>5} {'B':>5} {'DOF':>6} {'J_BB':>8} {'audA2':>8} {'n-5':>4} {'n-7':>4}"
    )
    for tag in order:
        c = by.get(tag)
        if not c:
            print(f"{tag:48} MISSING")
            continue
        print(
            f"{tag:48} {_fmt(c.get('bw_GHz'),2):>5} {_fmt(c.get('b_tesla'),3):>5} "
            f"{str(c.get('n_var')):>6} {_fmt(c.get('J_BB')):>8} {_fmt(c.get('audit_worst_A2')):>8} "
            f"{_fmt(c.get('audit_n_le5'),0):>4} {_fmt(c.get('audit_n_le7'),0):>4}"
        )
    print("-" * 72)
    print("BANDWIDTH THRESHOLD: ~0.2 GHz robust; ~0.5 GHz soft floor; ≥1.0 GHz notch-limited.")
    print("B-FIELD: 0.10 T helps Band C but does not clear notches.")
    print("DOF: C3 ≉ better than C6 on Band B → not primarily DOF-limited.")
    print("BOUNDS: best narrow-band solutions interior → widening bounds unlikely primary fix.")
    print("STAGE-6: NOT JUSTIFIED on current 4–6 GHz ask; change band/B/geometry first.")
    print("=" * 72)


def main():
    FEAS.mkdir(parents=True, exist_ok=True)
    cases = load_all_cases()
    # require D if opt exists complete; if still missing, still write partial but warn
    tags = {c["tag"] for c in cases}
    if "D_4p25_5p75_C6_B0.0500_Bs5" not in tags:
        print("WARNING: Band D opt/audit not yet in completed set")
    compact = write_table(cases)
    write_notch_report(cases)
    make_plots(cases)
    # also refresh ladder table via existing summarizer logic (inline call)
    path = write_study_md(cases, compact)
    print_conclusions(cases)
    print(f"wrote {path}")
    print(f"wrote {FEAS/'FEASIBILITY_COMPARISON_TABLE.md'}")
    print(f"wrote {FEAS/'FEASIBILITY_NOTCH_REPORT.md'}")
    print("figures: feasibility_worstA2_vs_bandwidth.png, feasibility_A2_overlays_*.png")


if __name__ == "__main__":
    main()
