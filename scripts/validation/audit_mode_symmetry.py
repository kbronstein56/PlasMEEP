#!/usr/bin/env python3
"""Audit P2/P3 numerical mode symmetry (launch weights, tangents, phases)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VAL = ROOT / "scripts" / "validation"
sys.path.insert(0, str(VAL))

from plasmeep.ports.mode_registry import audit_mode_cache, get_numerical_mode  # noqa: E402
from plasmeep.ports.numerical_launch import port_tangent  # noqa: E402
from sixport_common import effective_port_dir, fs_a, set_geometry_context  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "port_gates" / "mode_symmetry_audit.json"


def main() -> int:
    set_geometry_context(res=32, horn_walls="prism")
    ports = [1, 2]  # P2, P3
    report = {"res": 32, "ports": {}}
    for p in ports:
        audit = audit_mode_cache(p, res=32, frequency_a=fs_a)
        mode = get_numerical_mode(p, res=32, frequency_a=fs_a)
        u = effective_port_dir(p)
        t = port_tangent(u)
        amps = mode.source_amplitudes()
        report["ports"][f"P{p+1}"] = {
            "cache": audit.as_dict(),
            "outward_dir": u.tolist(),
            "tangent": t.tolist(),
            "mode_tangent": list(mode.tangent),
            "tangent_dot": float(mode.tangent[0] * t[0] + mode.tangent[1] * t[1]),
            "n_samples": len(mode.offsets_a),
            "offsets_match": mode.offsets_a.tolist()[:3],
            "weight_rms": float((abs(amps) ** 2).sum() ** 0.5),
            "phase_span_deg": float(
                (max(abs(__import__("numpy").angle(amps))) - min(abs(__import__("numpy").angle(amps))))
                * 180 / 3.14159
            ),
        }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
