"""
Load cached NumericalPortMode profiles for six-port validation.

Each device port uses its own reference profile extracted at that port's
orientation.  Symmetry classes (+60°/−60° pairs) do *not* share one JSON
file: P2 (60°) and P3 (120°) are adjacent hex faces with different tangents.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np

from plasmeep.ports.numerical_mode import NumericalPortMode, cache_key, load_cached_mode

logger = logging.getLogger(__name__)


def default_mode_cache_root() -> Path:
    env = os.environ.get("PLASMEEP_MODE_CACHE")
    if env:
        return Path(env)
    return (
        Path(__file__).resolve().parents[2]
        / "outputs"
        / "validation"
        / "mode_profiles"
    )


def canonical_profile_port(port_index: int) -> int:
    """Device port index equals cached reference port index (no sharing)."""
    if not 0 <= port_index <= 5:
        raise ValueError(f"port_index must be 0..5, got {port_index}")
    return port_index


def mode_cache_dir(
    res: int,
    *,
    base_dir: Optional[Path | str] = None,
    horn_walls: str = "prism",
    grid_offset_cells: Tuple[float, float] = (0.0, 0.0),
    coord_rotation_deg: float = 0.0,
    frequency_a: Optional[float] = None,
) -> Path:
    """
    Directory holding cached numerical modes for one geometry/resolution context.

    Legacy layouts store files directly under ``res{N}/``; newer keyed layouts
  use ``res{N}/<cache_key>/``.  ``resolve_numerical_mode`` checks both.
    """
    root = Path(base_dir or default_mode_cache_root())
    base = root / f"res{res}"
    if frequency_a is None:
        return base
    key = cache_key(
        res=res,
        frequency_a=frequency_a,
        horn_walls=horn_walls,
        grid_offset_cells=grid_offset_cells,
        coord_rotation_deg=coord_rotation_deg,
    )
    keyed = base / key
    return keyed if keyed.is_dir() else base


@dataclass(frozen=True)
class ModeCacheAudit:
    """Full identity of a numerical port-mode cache entry."""

    device_port_index: int
    profile_port_index: int
    cache_path: str
    res: int
    frequency_a: float
    horn_walls: str
    grid_offset_cells: Tuple[float, float]
    coord_rotation_deg: float
    mode_label: str
    mode_res: int
    mode_frequency_a: float
    outward_dir: Tuple[float, float]
    tangent: Tuple[float, float]
    n_samples: int
    res_match: bool
    frequency_match: bool

    def as_dict(self) -> Dict:
        return {
            "device_port": f"P{self.device_port_index + 1}",
            "profile_port": f"P{self.profile_port_index + 1}",
            "cache_path": self.cache_path,
            "cache_key": cache_key(
                res=self.res,
                frequency_a=self.frequency_a,
                horn_walls=self.horn_walls,
                grid_offset_cells=self.grid_offset_cells,
                coord_rotation_deg=self.coord_rotation_deg,
            ),
            "requested": {
                "res": self.res,
                "frequency_a": self.frequency_a,
                "horn_walls": self.horn_walls,
                "grid_offset_cells": list(self.grid_offset_cells),
                "coord_rotation_deg": self.coord_rotation_deg,
            },
            "stored_mode": {
                "label": self.mode_label,
                "res": self.mode_res,
                "frequency_a": self.mode_frequency_a,
                "outward_dir": list(self.outward_dir),
                "tangent": list(self.tangent),
                "n_samples": self.n_samples,
            },
            "res_match": self.res_match,
            "frequency_match": self.frequency_match,
        }

    def log(self) -> None:
        d = self.as_dict()
        logger.info(
            "mode cache P%d -> %s (res_match=%s freq_match=%s)",
            self.device_port_index + 1,
            self.cache_path,
            self.res_match,
            self.frequency_match,
        )
        print(
            f"  [mode-cache] P{self.device_port_index + 1} "
            f"file={self.cache_path} "
            f"label={self.mode_label} stored_res={self.mode_res} "
            f"req_res={self.res} res_match={self.res_match}",
            flush=True,
        )


def _candidate_paths(
    profile_port: int,
    *,
    res: int,
    frequency_a: float,
    cache_dir: Optional[Path | str],
    horn_walls: str,
    grid_offset_cells: Tuple[float, float],
    coord_rotation_deg: float,
) -> list[Path]:
    root = Path(cache_dir or default_mode_cache_root())
    fname = f"numerical_mode_P{profile_port + 1}.json"
    keyed = (
        root
        / f"res{res}"
        / cache_key(
            res=res,
            frequency_a=frequency_a,
            horn_walls=horn_walls,
            grid_offset_cells=grid_offset_cells,
            coord_rotation_deg=coord_rotation_deg,
        )
        / fname
    )
    legacy = root / f"res{res}" / fname
    paths = []
    if keyed.is_file():
        paths.append(keyed)
    if legacy.is_file() and legacy not in paths:
        paths.append(legacy)
    return paths


def resolve_numerical_mode(
    port_index: int,
    *,
    res: int,
    frequency_a: float,
    cache_dir: Optional[Path | str] = None,
    horn_walls: str = "prism",
    grid_offset_cells: Tuple[float, float] = (0.0, 0.0),
    coord_rotation_deg: float = 0.0,
) -> NumericalPortMode:
    """Return cached numerical mode for a device port."""
    profile_port = canonical_profile_port(port_index)
    paths = _candidate_paths(
        profile_port,
        res=res,
        frequency_a=frequency_a,
        cache_dir=cache_dir,
        horn_walls=horn_walls,
        grid_offset_cells=grid_offset_cells,
        coord_rotation_deg=coord_rotation_deg,
    )
    if not paths:
        legacy = (
            Path(cache_dir or default_mode_cache_root())
            / f"res{res}"
            / f"numerical_mode_P{profile_port + 1}.json"
        )
        raise FileNotFoundError(
            f"Missing numerical mode cache for P{port_index + 1} "
            f"(profile P{profile_port + 1}): searched {legacy}"
        )
    mode = load_cached_mode(paths[0])
    if mode is None:
        raise FileNotFoundError(f"Failed to load mode cache: {paths[0]}")
    return mode


def audit_mode_cache(
    port_index: int,
    *,
    res: int,
    frequency_a: float,
    cache_dir: Optional[Path | str] = None,
    horn_walls: str = "prism",
    grid_offset_cells: Tuple[float, float] = (0.0, 0.0),
    coord_rotation_deg: float = 0.0,
) -> ModeCacheAudit:
    profile_port = canonical_profile_port(port_index)
    paths = _candidate_paths(
        profile_port,
        res=res,
        frequency_a=frequency_a,
        cache_dir=cache_dir,
        horn_walls=horn_walls,
        grid_offset_cells=grid_offset_cells,
        coord_rotation_deg=coord_rotation_deg,
    )
    if not paths:
        raise FileNotFoundError(
            f"No mode cache for P{port_index + 1} at res{res}"
        )
    mode = load_cached_mode(paths[0])
    if mode is None:
        raise FileNotFoundError(paths[0])
    return ModeCacheAudit(
        device_port_index=port_index,
        profile_port_index=profile_port,
        cache_path=str(paths[0].resolve()),
        res=res,
        frequency_a=frequency_a,
        horn_walls=horn_walls,
        grid_offset_cells=grid_offset_cells,
        coord_rotation_deg=coord_rotation_deg,
        mode_label=mode.label,
        mode_res=mode.res,
        mode_frequency_a=mode.frequency_a,
        outward_dir=mode.outward_dir,
        tangent=mode.tangent,
        n_samples=len(mode.offsets_a),
        res_match=(mode.res == res),
        frequency_match=(abs(mode.frequency_a - frequency_a) < 1e-9),
    )


def validate_mode_alignment(
    mode: NumericalPortMode,
    port_outward: np.ndarray,
    port_tangent: np.ndarray,
    *,
    port_index: int,
    atol: float = 1e-3,
) -> None:
    """Raise if cached mode orientation does not match the device port."""
    u = np.asarray(port_outward, dtype=float).reshape(2)
    u = u / max(np.linalg.norm(u), 1e-30)
    t = np.asarray(port_tangent, dtype=float).reshape(2)
    t = t / max(np.linalg.norm(t), 1e-30)
    mo = np.asarray(mode.outward_dir, dtype=float)
    mt = np.asarray(mode.tangent, dtype=float)
    outward_dot = float(np.dot(u, mo))
    tangent_dot = float(np.dot(t, mt))
    if outward_dot < 1.0 - atol or tangent_dot < 1.0 - atol:
        raise ValueError(
            f"P{port_index + 1}: cached mode orientation mismatch "
            f"(outward_dot={outward_dot:.4f}, tangent_dot={tangent_dot:.4f}). "
            f"Expected profile from this port, not a symmetry-shared reference."
        )


_in_memory: Dict[str, NumericalPortMode] = {}


def get_numerical_mode(
    port_index: int,
    *,
    res: int,
    frequency_a: float,
    cache_dir: Optional[Path | str] = None,
    horn_walls: str = "prism",
    grid_offset_cells: Tuple[float, float] = (0.0, 0.0),
    coord_rotation_deg: float = 0.0,
    reload: bool = False,
    validate_alignment: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    log_audit: bool = False,
) -> NumericalPortMode:
    key = (
        f"{port_index}:{res}:{frequency_a}:{cache_dir}:"
        f"{horn_walls}:{grid_offset_cells}:{coord_rotation_deg}"
    )
    if not reload and key in _in_memory:
        mode = _in_memory[key]
    else:
        mode = resolve_numerical_mode(
            port_index,
            res=res,
            frequency_a=frequency_a,
            cache_dir=cache_dir,
            horn_walls=horn_walls,
            grid_offset_cells=grid_offset_cells,
            coord_rotation_deg=coord_rotation_deg,
        )
        audit = audit_mode_cache(
            port_index,
            res=res,
            frequency_a=frequency_a,
            cache_dir=cache_dir,
            horn_walls=horn_walls,
            grid_offset_cells=grid_offset_cells,
            coord_rotation_deg=coord_rotation_deg,
        )
        if not audit.res_match:
            raise ValueError(
                f"P{port_index + 1}: mode file res={audit.mode_res} != "
                f"simulation res={res} ({audit.cache_path})"
            )
        if not audit.frequency_match:
            raise ValueError(
                f"P{port_index + 1}: mode frequency mismatch "
                f"({audit.mode_frequency_a} vs {frequency_a})"
            )
        if log_audit:
            audit.log()
        _in_memory[key] = mode

    if validate_alignment is not None:
        validate_mode_alignment(
            mode,
            validate_alignment[0],
            validate_alignment[1],
            port_index=port_index,
        )
    return mode
