#!/usr/bin/env python3
"""Regression: compare_grids must tolerate checkpoint audits without worst_* keys.

Does not run FEM. Reproduces the KeyError: 'worst_isolation_dB' failure mode.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "scripts" / "circulator_inverse_design" / "jesse_method_campaign")]

from spectral_recon import compare_grids, ensure_summary_fields  # noqa: E402


def _fake_rows():
    return [
        {"f_GHz": 5.50, "isolation_dB": 1.0, "IL_dB": 5.0, "RL_dB": 8.0, "J_A": -0.2, "des_tot": 0.1, "rev_tot": 0.05},
        {"f_GHz": 5.55, "isolation_dB": -2.0, "IL_dB": 12.0, "RL_dB": 3.0, "J_A": -0.5, "des_tot": 0.05, "rev_tot": 0.08},
        {"f_GHz": 5.60, "isolation_dB": 0.5, "IL_dB": 6.0, "RL_dB": 7.0, "J_A": -0.3, "des_tot": 0.08, "rev_tot": 0.04},
    ]


def test_ensure_summary_fields_from_rows_only():
    blob = {"tag": "checkpoint", "grade": "C", "rows": _fake_rows()}
    # intentionally missing worst_isolation_dB
    assert "worst_isolation_dB" not in blob
    fixed = ensure_summary_fields(blob)
    assert "worst_isolation_dB" in fixed
    assert fixed["worst_isolation_dB"]["f_GHz"] == 5.55
    assert fixed["worst_IL_dB"]["f_GHz"] == 5.55


def test_compare_grids_accepts_rows_only_checkpoint():
    coarse = {
        "tag": "50",
        "grade": "C",
        "rows": _fake_rows(),
        "worst_isolation_dB": {"value": -2.0, "f_GHz": 5.55},
        "worst_IL_dB": {"value": 12.0, "f_GHz": 5.55},
        "worst_RL_dB": {"value": 3.0, "f_GHz": 5.55},
    }
    # fine is a mid-scan checkpoint WITHOUT worst_* — previously KeyError
    fine = {"tag": "25", "grade": "C", "rows": _fake_rows()}
    out = compare_grids(coarse, fine)
    assert "d_worst_iso_dB" in out
    assert out["d_worst_iso_dB"] == 0.0


def test_compare_grids_partial_summary_missing_RL():
    """Matches on-disk bb_final_audit_25MHz_local.json schema at failure time:
    had worst_isolation_dB / worst_IL_dB but no worst_RL_dB.
    """
    coarse = {
        "tag": "50",
        "grade": "C",
        "rows": _fake_rows(),
        "worst_isolation_dB": {"value": -2.0, "f_GHz": 5.55},
        "worst_IL_dB": {"value": 12.0, "f_GHz": 5.55},
        "worst_RL_dB": {"value": 3.0, "f_GHz": 5.55},
    }
    fine = {
        "tag": "25",
        "grade": "C",
        "rows": _fake_rows(),
        "worst_isolation_dB": {"value": -2.0, "f_GHz": 5.55},
        "worst_IL_dB": {"value": 12.0, "f_GHz": 5.55},
        # intentionally omit worst_RL_dB
    }
    out = compare_grids(coarse, fine)
    assert out["d_worst_RL_dB"] == 0.0


def test_compare_grids_empty_rows_raises():
    try:
        compare_grids({"rows": []}, {"rows": []})
    except KeyError as e:
        assert "rows" in str(e) or "worst" in str(e).lower() or "missing" in str(e).lower()
    else:
        raise AssertionError("expected KeyError for empty audit blobs")


if __name__ == "__main__":
    test_ensure_summary_fields_from_rows_only()
    test_compare_grids_accepts_rows_only_checkpoint()
    test_compare_grids_partial_summary_missing_RL()
    test_compare_grids_empty_rows_raises()
    print("test_audit_summary_schema: PASS")
