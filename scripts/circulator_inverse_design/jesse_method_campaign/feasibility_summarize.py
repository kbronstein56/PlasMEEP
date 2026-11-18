#!/usr/bin/env python3
"""Aggregate feasibility ladder results into compact tables + tradeoff plots.

Stage-6 NOT STARTED.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
FEAS = OUT / "feasibility"


def load_case(opt_path: Path) -> dict:
    d = json.loads(opt_path.read_text())
    hist = d.get("history") or []
    last = next((h for h in reversed(hist) if h.get("kind") == "accepted_iter"), None)
    rec = d.get("rec_best") or {}
    wm = rec.get("worst_metrics") or {}
    q = np.asarray(d.get("q_best") or [], float)
    row = {
        "tag": opt_path.stem.replace("_opt", ""),
        "band": d.get("band"),
        "f0": d.get("f0_GHz"),
        "f1": d.get("f1_GHz"),
        "bw_GHz": float(d["f1_GHz"]) - float(d["f0_GHz"]) if d.get("f0_GHz") is not None else None,
        "n_freq": d.get("n_freq"),
        "symmetry": d.get("symmetry"),
        "n_var": d.get("n_var"),
        "b_tesla": d.get("b_tesla"),
        "J_BB": d.get("J_BB_best"),
        "opt_worst_A2": rec.get("worst_J"),
        "opt_worst_f": rec.get("worst_f_GHz"),
        "opt_IL": wm.get("IL_dB") if wm else (last or {}).get("IL_dB"),
        "opt_iso": wm.get("isolation_dB") if wm else (last or {}).get("isolation_dB"),
        "grad_norm": rec.get("grad_norm") or (last or {}).get("grad_norm"),
        "n_lo": int(np.sum(q <= 0.07)) if q.size else None,
        "n_hi": int(np.sum(q >= 1.78)) if q.size else None,
        "n_eval": d.get("n_eval"),
        "accepted_nit": d.get("accepted_nit"),
        "elapsed_h": (d.get("elapsed_s") or 0) / 3600.0,
        "optimizer_message": d.get("optimizer_message"),
        "status": d.get("status"),
    }
    # attach audit summary if present
    audits = sorted(FEAS.glob(f"audit_{opt_path.stem}_*_SUMMARY.json"))
    if audits:
        a = json.loads(audits[-1].read_text())
        w = a.get("worst_A2") or {}
        row.update(
            {
                "audit_worst_A2": w.get("J_A2"),
                "audit_worst_f": w.get("f_GHz"),
                "audit_IL": w.get("IL_dB"),
                "audit_iso": w.get("isolation_dB"),
                "audit_n_le5": a.get("n_A2_le_m5"),
                "audit_n_le7": a.get("n_A2_le_m7"),
                "audit_worst_IL": a.get("worst_IL_dB"),
                "audit_worst_iso": a.get("worst_isolation_dB"),
            }
        )
    return row


def main():
    FEAS.mkdir(parents=True, exist_ok=True)
    opts = sorted(FEAS.glob("*_opt.json"))
    rows = [load_case(p) for p in opts if json.loads(p.read_text()).get("status") == "complete"]

    # Stage-5 / Band E reference from campaign outputs
    s5 = OUT / "opt_broadband_stage5_A2_C.json"
    if s5.exists():
        d = json.loads(s5.read_text())
        e = {
            "tag": "E_4p0_6p0_C6_B0.0500_stage5",
            "band": "E_4p0_6p0",
            "f0": 4.0,
            "f1": 6.0,
            "bw_GHz": 2.0,
            "n_freq": 321,
            "symmetry": "C6",
            "n_var": 16,
            "b_tesla": 0.05,
            "J_BB": d.get("J_BB_best"),
            "n_eval": d.get("n_eval"),
            "accepted_nit": d.get("accepted_nit"),
            "elapsed_h": (d.get("elapsed_s") or 0) / 3600.0,
            "status": "complete",
        }
        # prefer dense audit summary if present
        for cand in [
            OUT / "BB_STAGE5_DENSE_AUDIT.json",
            OUT / "spectral_bb_stage5_dense_3p125MHz_local.json",
            OUT / "spectral_bb_stage5_dense_6p25MHz.json",
        ]:
            if not cand.exists():
                continue
            a = json.loads(cand.read_text())
            rows_a = a.get("rows") or []
            if not rows_a and "worst_A2" in a:
                w = a["worst_A2"]
                e.update(
                    {
                        "audit_worst_A2": w.get("J_A2") or w.get("A2"),
                        "audit_worst_f": w.get("f_GHz"),
                        "audit_IL": w.get("IL_dB"),
                        "audit_iso": w.get("isolation_dB"),
                    }
                )
                break
            if rows_a:
                worst = min(rows_a, key=lambda r: r["J_A2"])
                e.update(
                    {
                        "audit_worst_A2": worst["J_A2"],
                        "audit_worst_f": worst["f_GHz"],
                        "audit_IL": worst["IL_dB"],
                        "audit_iso": worst["isolation_dB"],
                        "audit_n_le5": sum(1 for r in rows_a if r["J_A2"] <= -5),
                        "audit_n_le7": sum(1 for r in rows_a if r["J_A2"] <= -7),
                        "audit_source": cand.name,
                    }
                )
                break
        rows.append(e)

    (FEAS / "FEASIBILITY_LADDER_TABLE.json").write_text(json.dumps(rows, indent=2) + "\n")

    # Markdown table
    keys = [
        "tag",
        "bw_GHz",
        "n_freq",
        "n_var",
        "b_tesla",
        "J_BB",
        "audit_worst_A2",
        "audit_IL",
        "audit_iso",
        "audit_n_le5",
        "audit_n_le7",
        "elapsed_h",
        "n_lo",
        "n_hi",
    ]
    lines = [
        "# Feasibility ladder table",
        "",
        "Stage-6: NOT STARTED / AWAITING AUTHORIZATION",
        "",
        "| " + " | ".join(keys) + " |",
        "| " + " | ".join(["---"] * len(keys)) + " |",
    ]
    for r in sorted(rows, key=lambda x: (x.get("bw_GHz") or 99, x.get("tag") or "")):
        lines.append(
            "| "
            + " | ".join(
                ("—" if r.get(k) is None else (f"{r[k]:.4f}" if isinstance(r.get(k), float) else str(r.get(k))))
                for k in keys
            )
            + " |"
        )
    (FEAS / "FEASIBILITY_LADDER_TABLE.md").write_text("\n".join(lines) + "\n")

    # Bandwidth tradeoff plots (C6 / B=0.05 only)
    bw_rows = [
        r
        for r in rows
        if r.get("symmetry", "C6") == "C6"
        and abs(float(r.get("b_tesla") or 0.05) - 0.05) < 1e-9
        and r.get("bw_GHz") is not None
        and "uniform" not in (r.get("tag") or "")
    ]
    bw_rows = sorted(bw_rows, key=lambda r: r["bw_GHz"])
    if bw_rows:
        fig, axes = plt.subplots(2, 2, figsize=(9, 7), constrained_layout=True)
        xs = [r["bw_GHz"] for r in bw_rows]
        series = [
            (axes[0, 0], "audit_worst_A2", "Worst A2 (audit)"),
            (axes[0, 1], "audit_IL", "Worst-site IL [dB]"),
            (axes[1, 0], "audit_iso", "Worst-site isolation [dB]"),
            (axes[1, 1], "J_BB", "Optimized J_BB"),
        ]
        for ax, key, title in series:
            ys = [r.get(key) for r in bw_rows]
            ax.plot(xs, ys, "o-", color="#1f4e79")
            ax.set_xlabel("Required bandwidth [GHz]")
            ax.set_ylabel(title)
            ax.set_title(title)
            ax.grid(True, alpha=0.3)
        fig.suptitle("Feasibility ladder — C6, B=0.05 T (Stage-6 NOT STARTED)")
        fig.savefig(FEAS / "feasibility_tradeoff_bandwidth.png", dpi=140)
        plt.close(fig)

    # B sensitivity
    b_rows = [r for r in rows if r.get("band") in ("C_4p5_5p5",) or (r.get("tag") or "").startswith("C_4p5_5p5")]
    if len(b_rows) >= 1:
        fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
        for r in sorted(b_rows, key=lambda x: float(x.get("b_tesla") or 0)):
            ax.scatter(
                [r.get("b_tesla")],
                [r.get("audit_worst_A2", r.get("opt_worst_A2"))],
                s=80,
                label=r["tag"],
            )
        ax.set_xlabel("B [T]")
        ax.set_ylabel("Worst A2")
        ax.set_title("B-field probe (band C)")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
        fig.savefig(FEAS / "feasibility_tradeoff_B.png", dpi=140)
        plt.close(fig)

    # DOF
    dof = [r for r in rows if "B_4p75_5p25" in (r.get("tag") or "")]
    if dof:
        fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
        for r in dof:
            ax.scatter([r.get("n_var")], [r.get("audit_worst_A2", r.get("opt_worst_A2"))], s=80, label=r["tag"])
        ax.set_xlabel("Independent density DOFs")
        ax.set_ylabel("Worst A2")
        ax.set_title("Design-DOF probe (band B)")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
        fig.savefig(FEAS / "feasibility_tradeoff_DOF.png", dpi=140)
        plt.close(fig)

    print(f"wrote {FEAS/'FEASIBILITY_LADDER_TABLE.md'} n_rows={len(rows)}")
    print("STAGE-6: NOT STARTED / AWAITING AUTHORIZATION")


if __name__ == "__main__":
    main()
