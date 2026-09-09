#!/usr/bin/env python3
"""
Finalize overnight S-param audit into MASTER_REPORT §20 results tables.

Reads outputs/validation/sparam_audit/*.json case files and rewrites the
results subsection. Safe to re-run as cases complete.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "outputs" / "validation" / "sparam_audit"
REPORT = ROOT / "outputs" / "validation" / "MASTER_REPORT.md"


def _load_cases() -> List[Dict[str, Any]]:
    cases = []
    for p in sorted(AUDIT.glob("horns_dual_*.json")) + sorted(AUDIT.glob("full_dual_*.json")):
        try:
            cases.append(json.loads(p.read_text(encoding="utf-8")))
        except Exception:
            continue
    return cases


def _row(data: Dict[str, Any]) -> str:
    label = data.get("label", "?")
    res = data.get("results", {})
    pairs = res.get("pair_summaries") or []
    pair = pairs[0] if pairs else {}
    disc = (res.get("discrete_control") or {}).get("amp_err_dB")
    flux = pair.get("flux_diff_dB")
    modal = pair.get("modal_power_diff_dB")
    healthy = (res.get("flux_subtraction_health") or {}).get("healthy")
    skip = data.get("settings", {}).get("skip_flux_subtraction")
    mon = data.get("settings", {}).get("monitor_outward_cells", 0)
    mode = data.get("settings", {}).get("device_mode")
    wall = (data.get("timings_s") or {}).get("total")
    return (
        f"| `{label}` | {mode} | skip_sub={skip} mon={mon} | "
        f"{_fmt(flux)} | {_fmt(modal)} | {_fmt(disc)} | {healthy} | {_fmt(wall, '.0f')}s |"
    )


def _fmt(v: Any, fmt: str = ".4g") -> str:
    if v is None:
        return "—"
    try:
        return format(float(v), fmt)
    except Exception:
        return str(v)


def build_section(cases: List[Dict[str, Any]]) -> str:
    lines = [
        "### Overnight results (auto-filled)",
        "",
        "| Case | Device | Settings | Flux Δ dB | Modal Δ dB | Discrete dB | Flux-sub healthy | Wall |",
        "|---|---|---|---:|---:|---:|---|---:|",
    ]
    if not cases:
        lines.append("| *(no case JSON yet)* | | | | | | | |")
    else:
        for c in cases:
            lines.append(_row(c))

    # Decision bullets
    horns = [c for c in cases if c.get("settings", {}).get("device_mode") == "horns_only"]
    full = [c for c in cases if c.get("settings", {}).get("device_mode") == "full"]

    def first_pair_flux(cs, tag):
        for c in cs:
            if tag in c.get("label", ""):
                pairs = c.get("results", {}).get("pair_summaries") or []
                if pairs:
                    return pairs[0]
        return None

    lines += ["", "### Decision-tree status", ""]
    hp = first_pair_flux(horns, "P1P2")
    fp = first_pair_flux(full, "P1P2")
    if hp:
        lines.append(
            f"- Horns-only P1↔P2: flux Δ={_fmt(hp.get('flux_diff_dB'))} dB, "
            f"modal Δ={_fmt(hp.get('modal_power_diff_dB'))} dB"
        )
    else:
        lines.append("- Horns-only P1↔P2: pending")
    if fp:
        lines.append(
            f"- Full PMM P1↔P2: flux Δ={_fmt(fp.get('flux_diff_dB'))} dB, "
            f"modal Δ={_fmt(fp.get('modal_power_diff_dB'))} dB"
        )
        # Compare amplification
        if hp and hp.get("flux_diff_dB") is not None and fp.get("flux_diff_dB") is not None:
            h = abs(float(hp["flux_diff_dB"]))
            f = abs(float(fp["flux_diff_dB"]))
            if h < 0.2 and f > 1.0:
                lines.append("- **Interpretation:** horns nearly OK; full PMM amplifies a small port error.")
            elif h > 1.0:
                lines.append("- **Interpretation:** port method already biased with horns only — fix ports before more PMM.")
            modal_f = fp.get("modal_power_diff_dB")
            flux_f = fp.get("flux_diff_dB")
            if modal_f is not None and flux_f is not None:
                if abs(float(modal_f)) < 0.2 and abs(float(flux_f)) > 1.0:
                    lines.append("- **Interpretation:** direct modal overlap recovers reciprocity → flux integration path is the culprit (not subtraction, per Phase 1).")
                elif abs(float(modal_f)) > 1.0 and abs(float(flux_f)) > 1.0:
                    lines.append("- **Interpretation:** both flux and modal fail → launch/grid/mode definition, not flux subtraction.")
    else:
        lines.append("- Full PMM P1↔P2: pending")

    # nosub check
    sub = next((c for c in full if c.get("label", "").endswith("_sub") and "P1P2" in c.get("label", "")), None)
    nosub = next((c for c in full if c.get("label", "").endswith("_nosub")), None)
    if sub and nosub:
        ps = (sub.get("results", {}).get("pair_summaries") or [{}])[0]
        pn = (nosub.get("results", {}).get("pair_summaries") or [{}])[0]
        lines.append(
            f"- Full P1↔P2 sub vs nosub flux Δ: {_fmt(ps.get('flux_diff_dB'))} vs {_fmt(pn.get('flux_diff_dB'))} "
            "(expect nearly identical transmission)"
        )

    lines.append("")
    return "\n".join(lines)


def main() -> None:
    cases = _load_cases()
    section = build_section(cases)
    text = REPORT.read_text(encoding="utf-8")
    marker = "### Overnight results (auto-filled)"
    end_marker = "### Phase 7 (modes)"
    if marker in text:
        pre = text.split(marker)[0]
        if end_marker in text:
            post = end_marker + text.split(end_marker, 1)[1]
        else:
            post = ""
        new = pre + section + "\n" + post
    else:
        # Insert before Phase 7 or at end of §20
        if end_marker in text:
            pre, post = text.split(end_marker, 1)
            new = pre + section + "\n" + end_marker + post
        else:
            new = text.rstrip() + "\n\n" + section + "\n"
    REPORT.write_text(new, encoding="utf-8")
    print(f"Updated {REPORT} with {len(cases)} cases")


if __name__ == "__main__":
    main()
