"""
Analytic TE1-like profiles and overlap metrics for port-mode validation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Tuple

import numpy as np


def te1_cos_weights(offsets: np.ndarray, span: float) -> np.ndarray:
    """Normalized cos(π s/L) weights matching port_formulations.make_te1_hz_sources."""
    offsets = np.asarray(offsets, dtype=float)
    span = float(span)
    weights = np.cos(np.pi * offsets / span)
    weights = np.clip(weights, 0.0, None)
    wsum = float(np.sum(weights))
    if wsum <= 0:
        return np.ones_like(offsets) / max(len(offsets), 1)
    return weights / wsum


def te1_cos_amplitude(offsets: np.ndarray, span: float) -> np.ndarray:
    """Unnormalized cos envelope (peak = 1 at center)."""
    offsets = np.asarray(offsets, dtype=float)
    span = float(span)
    return np.cos(np.pi * offsets / span).clip(min=0.0)


@dataclass(frozen=True)
class ModeLineProfile:
    """Complex Hz sampled along horn tangent at the monitor plane."""

    port_index: int
    offsets_a: np.ndarray
    hz_complex: np.ndarray
    span_a: float
    monitor_xy: Tuple[float, float]
    tangent: Tuple[float, float]

    def normalized(self) -> np.ndarray:
        """Unit L2 norm on sampled line (discrete)."""
        f = np.asarray(self.hz_complex, dtype=complex)
        nrm = float(np.sqrt(np.sum(np.abs(f) ** 2)))
        if nrm <= 0:
            return f
        return f / nrm

    def as_dict(self) -> Dict[str, Any]:
        return {
            "port_index": self.port_index,
            "offsets_a": self.offsets_a.tolist(),
            "hz_complex": [
                {"real": z.real, "imag": z.imag} for z in self.hz_complex
            ],
            "span_a": self.span_a,
            "monitor_xy": list(self.monitor_xy),
            "tangent": list(self.tangent),
        }


def complex_overlap(a: np.ndarray, b: np.ndarray) -> complex:
    """Discrete ⟨a|b⟩ with both vectors L2-normalized."""
    a = np.asarray(a, dtype=complex).reshape(-1)
    b = np.asarray(b, dtype=complex).reshape(-1)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch {a.shape} vs {b.shape}")
    na = float(np.sqrt(np.sum(np.abs(a) ** 2)))
    nb = float(np.sqrt(np.sum(np.abs(b) ** 2)))
    if na <= 0 or nb <= 0:
        return 0.0 + 0.0j
    return complex(np.vdot(a / na, b / nb))


def overlap_metrics(
    numerical: np.ndarray,
    analytic_real: np.ndarray,
) -> Dict[str, float]:
    """
    Overlap between numerical Hz profile and a real analytic envelope.

    Returns |⟨n|a⟩|, power overlap, and phase of projection.
    """
    n = np.asarray(numerical, dtype=complex).reshape(-1)
    a = np.asarray(analytic_real, dtype=float).reshape(-1)
    ov = complex_overlap(n, a.astype(complex))
    power_frac = float(np.abs(ov) ** 2)
    return {
        "overlap_magnitude": float(np.abs(ov)),
        "power_overlap": power_frac,
        "projection_phase_deg": float(np.degrees(np.angle(ov))),
        "orthogonal_power_frac": float(max(0.0, 1.0 - power_frac)),
    }
