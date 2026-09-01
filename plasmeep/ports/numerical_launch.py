"""
Complex numerical-port-mode Meep source construction.
"""

from __future__ import annotations

from typing import List, Optional

import meep as mp
import numpy as np

from plasmeep.ports.numerical_mode import NumericalPortMode


def make_numerical_hz_sources(
    mode: NumericalPortMode,
    *,
    center_xy: np.ndarray,
    tangent_xy: np.ndarray,
    frequency: float,
    fwidth: float,
    target_power: float = 1.0,
) -> List[mp.Source]:
    """
  Launch a cached numerical Hz port mode.

  Normalization convention
  ------------------------
  Let φ_k be the L2-normalized discrete mode samples on offsets s_k.
  Source amplitudes are a_k = sqrt(P0) * conj(φ_k) / ||φ||_2 so that
  ∑_k |a_k|^2 = P0 (default P0=1 before incident normalization scaling).

  Complex phase of φ is preserved on each Meep source.
  """
    center_xy = np.asarray(center_xy, dtype=float).reshape(2)
    tangent = np.asarray(tangent_xy, dtype=float).reshape(2)
    tangent = tangent / max(np.linalg.norm(tangent), 1e-30)
    weights = mode.source_amplitudes(target_power=target_power)
    offsets = np.asarray(mode.offsets_a, dtype=float)
    if len(weights) != len(offsets):
        raise ValueError("mode offsets and weights length mismatch")

    sources: List[mp.Source] = []
    for s, amp in zip(offsets, weights):
        xy = center_xy + s * tangent
        sources.append(
            mp.Source(
                src=mp.GaussianSource(frequency=frequency, fwidth=fwidth),
                component=mp.Hz,
                center=mp.Vector3(xy[0], xy[1], 0),
                amplitude=complex(amp),
            )
        )
    return sources


def port_tangent(outward_dir: np.ndarray) -> np.ndarray:
    u = np.asarray(outward_dir, dtype=float).reshape(2)
    u = u / max(np.linalg.norm(u), 1e-30)
    return np.array([-u[1], u[0]])
