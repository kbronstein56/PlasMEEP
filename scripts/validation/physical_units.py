"""
Physical resolution and normalization-length helpers for validation harnesses.

User-facing resolution is ``points_per_cm`` (physical grid density). Internal Meep
``res`` is derived from the normalization length ``a_m``:

  dx_mm = 10 / points_per_cm
  res   = round(a_m / dx_m)
"""

from __future__ import annotations

import argparse
from typing import Any, Dict, Optional, Tuple

# Meep length unit = 20 mm lattice center-to-center pitch (Rodriguez et al. 2026).
A_M_DEFAULT = 0.020


def dx_mm_from_points_per_cm(points_per_cm: float) -> float:
    return 10.0 / float(points_per_cm)


def dx_m_from_points_per_cm(points_per_cm: float) -> float:
    return dx_mm_from_points_per_cm(points_per_cm) * 1e-3


def meep_resolution_from_points_per_cm(
    points_per_cm: float, *, a_m: float = A_M_DEFAULT
) -> int:
    dx_m = dx_m_from_points_per_cm(points_per_cm)
    return max(1, int(round(float(a_m) / dx_m)))


def points_per_cm_from_meep_resolution(res: int, *, a_m: float = A_M_DEFAULT) -> float:
    dx_m = float(a_m) / float(res)
    return 10.0 / (dx_m * 1000.0)


def physical_resolution_report(
    *,
    points_per_cm: Optional[float] = None,
    res: Optional[int] = None,
    a_m: float = A_M_DEFAULT,
    fs_Hz: float = 3.85e9,
) -> Dict[str, float]:
    """Return the canonical resolution block for JSON logs and reports."""
    if points_per_cm is None and res is None:
        raise ValueError("provide points_per_cm or res")
    if points_per_cm is None:
        res_i = int(res)
        ppc = points_per_cm_from_meep_resolution(res_i, a_m=a_m)
    else:
        ppc = float(points_per_cm)
        res_i = meep_resolution_from_points_per_cm(ppc, a_m=a_m)
    dx_mm = dx_mm_from_points_per_cm(ppc)
    dx_a = 1.0 / float(res_i)
    a_cm = float(a_m) * 100.0
    c_mps = 2.99792458e8
    lambda0_mm = c_mps / fs_Hz * 1000.0
    return {
        "points_per_cm": ppc,
        "dx_mm": dx_mm,
        "meep_resolution": float(res_i),
        "a_m": float(a_m),
        "a_cm": a_cm,
        "grid_spacing_a": dx_a,
        "grid_spacing_mm": dx_mm,
        "pixels_per_lattice_20mm": 20.0 / dx_mm,
        "pixels_per_free_space_lambda0": lambda0_mm / dx_mm,
        "lambda0_mm": lambda0_mm,
        # Legacy alias retained for older plotting scripts.
        "pixels_per_cm": ppc,
    }


def add_resolution_arguments(parser: argparse.ArgumentParser) -> None:
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--points-per-cm",
        type=float,
        dest="points_per_cm",
        help="Physical grid density (primary). dx_mm = 10/points_per_cm.",
    )
    group.add_argument(
        "--res",
        type=int,
        help="Internal Meep resolution (legacy; prefer --points-per-cm).",
    )


def resolve_simulation_resolution(
    args: argparse.Namespace,
    *,
    a_m: float = A_M_DEFAULT,
    default_points_per_cm: Optional[float] = None,
    default_res: Optional[int] = None,
) -> Tuple[int, float, Dict[str, float]]:
    """
    Return (meep_res, points_per_cm, physical_resolution_report dict).
    """
    if getattr(args, "points_per_cm", None) is not None:
        ppc = float(args.points_per_cm)
        res = meep_resolution_from_points_per_cm(ppc, a_m=a_m)
    elif getattr(args, "res", None) is not None:
        res = int(args.res)
        ppc = points_per_cm_from_meep_resolution(res, a_m=a_m)
    elif default_points_per_cm is not None:
        ppc = float(default_points_per_cm)
        res = meep_resolution_from_points_per_cm(ppc, a_m=a_m)
    elif default_res is not None:
        res = int(default_res)
        ppc = points_per_cm_from_meep_resolution(res, a_m=a_m)
    else:
        raise ValueError("no resolution specified")
    report = physical_resolution_report(points_per_cm=ppc, a_m=a_m)
    return res, ppc, report
