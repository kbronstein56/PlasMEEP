#!/usr/bin/env python3
"""Exact differences between sentinel_repro.json and the stored validation records."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"
NEW = json.loads((OUT / "sentinel_repro.json").read_text())


def diff(name, new, old):
    return {"name": name, "new": new, "stored": old, "delta": new - old}


def main():
    rows = []
    if "mms" in NEW:
        level = NEW["mms"]
        rows += [
            diff("mms_L2", level["L2"], 0.00038972516086799877),
            diff("mms_H1", level["H1"], 0.06603468609492105),
            diff("mms_max_nodal", level["max_nodal"], 0.00021657680783101593),
        ]
    if "planar" in NEW:
        g = NEW["planar"]["pec_guide_h0025"]
        rows += [
            diff("pec_L2", g["L2"], 0.0004988688378469143),
            diff("pec_rel_L2", g["rel_L2"], 0.0014262934809217242),
        ]
        s = NEW["planar"]["plasma_slab_0.40_h002"]
        rows += [
            diff("slab_T_db", s["T_db"], -0.0015444837813689536),
            diff("slab_T_phase_deg", s["T_phase_deg"], 0.0253770245186856),
            diff("slab_R_db", s["R_db"], -0.0010216396477051711),
        ]
    if "oblique" in NEW:
        stored = {0.0: 0.00013435979587412896, 25.0: 0.00013693452336939016, 60.0: 0.00021776852108119522}
        for rec in NEW["oblique"]:
            rows.append(diff(f"oblique_R_{rec['angle_deg']}", rec["R_abs"], stored[rec["angle_deg"]]))
    if "recip" in NEW:
        rows.append(diff("recip_rel", NEW["recip"]["rel"], 0.0))
        rows.append(diff("recip_abs_diff", NEW["recip"]["abs_diff"], 0.0))
    if "power" in NEW:
        rows.append(diff("power_R_plus_T", NEW["power"]["R_plus_T"], 0.9996778442034607))
    if "guide" in NEW:
        g = NEW["guide"]
        rows += [
            diff("guide_hz_l2", g["hz_l2"], 7.784834586015421e-05),
            diff("guide_grad_l2", g["grad_l2"], 0.0008357086854731789),
            diff("guide_EH_max", g["element_EH_max"], 0.0008966596176084519),
            diff("guide_power_rel", g["power_rel"], 0.00011155132664959487),
        ]
    if "coated" in NEW:
        fwd = NEW["coated"]["probes"]["forward"]
        near = NEW["coated"]["probes"]["near_quartz"]
        rows += [
            diff("coated_forward_db", fwd["db"], -0.011769872989058922),
            diff("coated_forward_phase", fwd["phase_deg"], -0.4611246863224272),
            diff("coated_near_db", near["db"], 0.047414214391662046),
            diff("coated_near_phase", near["phase_deg"], -0.04603344227972082),
            diff("coated_ring_l2", NEW["coated"]["ring_l2"], 0.004238897679460683),
            diff("coated_dofs", float(NEW["coated"]["dofs"]), 5940235.0),
        ]
    scatter = json.loads((OUT / "scatter_partial.json").read_text())

    def stored_probe(case, h, h_edge, probe, key):
        for rec in scatter:
            if rec["case"] == case and abs(rec["h"] - h) < 1e-9 and abs(rec["h_edge"] - h_edge) < 1e-9:
                return rec["probes"][probe][key]
        raise KeyError(case)

    if "bare" in NEW:
        fwd = NEW["bare"]["probes"]["forward"]
        rows += [
            diff("bare_forward_db", fwd["db"], stored_probe("bare_prod", 0.04, 0.012, "forward", "db")),
            diff("bare_forward_phase", fwd["phase_deg"], stored_probe("bare_prod", 0.04, 0.012, "forward", "phase_deg")),
        ]
    if "three" in NEW:
        fwd = NEW["three"]["probes"]["forward"]
        rows += [
            diff("three_forward_db", fwd["db"], stored_probe("coated_3", 0.04, 0.012, "forward", "db")),
            diff("three_forward_phase", fwd["phase_deg"], stored_probe("coated_3", 0.04, 0.012, "forward", "phase_deg")),
        ]
    if "seven" in NEW:
        fwd = NEW["seven"]["probes"]["forward"]
        near = NEW["seven"]["probes"]["near_quartz"]
        rows += [
            diff("seven_forward_db", fwd["db"], -0.008965716930644118),
            diff("seven_forward_phase", fwd["phase_deg"], -0.24347867232748555),
            diff("seven_near_db", near["db"], -0.020395264772488195),
            diff("seven_near_phase", near["phase_deg"], -0.3352348109653991),
            diff("seven_dofs", float(NEW["seven"]["dofs"]), 618177.0),
        ]
    (OUT / "sentinel_diffs.json").write_text(json.dumps(rows, indent=2) + "\n")
    for row in rows:
        print(f"{row['name']}: delta={row['delta']:.3e} new={row['new']} stored={row['stored']}")


if __name__ == "__main__":
    main()
