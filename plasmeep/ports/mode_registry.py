"""
Load cached NumericalPortMode profiles for six-port validation.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Optional

from plasmeep.ports.numerical_mode import NumericalPortMode, load_cached_mode


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
    """Map device port index to the cached reference profile port."""
    if port_index == 0:
        return 0
    if port_index in (1, 2):
        return 1
    if port_index in (4, 5):
        return 5
    if port_index == 3:
        return 0
    return port_index


def mode_cache_dir(
    res: int,
    *,
    base_dir: Optional[Path | str] = None,
) -> Path:
    return Path(base_dir or default_mode_cache_root()) / f"res{res}"


def resolve_numerical_mode(
    port_index: int,
    *,
    res: int,
    frequency_a: float,
    cache_dir: Optional[Path | str] = None,
) -> NumericalPortMode:
    """Return cached numerical mode for a device port."""
    cache_dir = mode_cache_dir(res, base_dir=cache_dir)
    profile_port = canonical_profile_port(port_index)
    path = cache_dir / f"numerical_mode_P{profile_port + 1}.json"
    mode = load_cached_mode(path)
    if mode is None:
        raise FileNotFoundError(
            f"Missing numerical mode cache: {path} (device port P{port_index + 1})"
        )
    return mode


_in_memory: Dict[str, NumericalPortMode] = {}


def get_numerical_mode(
    port_index: int,
    *,
    res: int,
    frequency_a: float,
    cache_dir: Optional[Path | str] = None,
    reload: bool = False,
) -> NumericalPortMode:
    key = f"{port_index}:{res}:{frequency_a}:{cache_dir}"
    if not reload and key in _in_memory:
        return _in_memory[key]
    mode = resolve_numerical_mode(
        port_index, res=res, frequency_a=frequency_a, cache_dir=cache_dir
    )
    _in_memory[key] = mode
    return mode
