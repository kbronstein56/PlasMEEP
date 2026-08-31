"""
Port-specific numerical mode references from isolated/reference-guide runs.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from plasmeep.ports.mode_profile import ModeLineProfile, complex_overlap, te1_cos_amplitude


@dataclass(frozen=True)
class NumericalPortMode:
    """Stored complex Hz profile for one horn orientation."""

    port_index: int
    label: str
    offsets_a: np.ndarray
    hz_complex: np.ndarray
    span_a: float
    frequency_a: float
    res: int
    monitor_xy: Tuple[float, float]
    tangent: Tuple[float, float]
    outward_dir: Tuple[float, float]

    @classmethod
    def from_profile(
        cls,
        profile: ModeLineProfile,
        *,
        label: str,
        frequency_a: float,
        res: int,
        outward_dir: Sequence[float],
    ) -> "NumericalPortMode":
        return cls(
            port_index=profile.port_index,
            label=label,
            offsets_a=np.asarray(profile.offsets_a, dtype=float),
            hz_complex=np.asarray(profile.hz_complex, dtype=complex),
            span_a=profile.span_a,
            frequency_a=frequency_a,
            res=res,
            monitor_xy=profile.monitor_xy,
            tangent=profile.tangent,
            outward_dir=tuple(float(x) for x in outward_dir),
        )

    def normalized_field(self) -> np.ndarray:
        f = np.asarray(self.hz_complex, dtype=complex)
        nrm = float(np.sqrt(np.sum(np.abs(f) ** 2)))
        return f if nrm <= 0 else f / nrm

    def overlap_with_analytic_te1(self) -> Dict[str, float]:
        from plasmeep.ports.mode_profile import overlap_metrics

        analytic = te1_cos_amplitude(self.offsets_a, self.span_a)
        return overlap_metrics(self.hz_complex, analytic)

    def overlap_with(self, other: "NumericalPortMode") -> complex:
        return complex_overlap(self.normalized_field(), other.normalized_field())

    def source_amplitudes(self, target_power: float = 1.0) -> np.ndarray:
        """Discrete source weights proportional to conj(mode), unit launch power."""
        mode = self.normalized_field()
        weights = np.conj(mode)
        wsum = float(np.sum(np.abs(weights) ** 2))
        if wsum <= 0:
            return np.ones_like(weights) / max(len(weights), 1)
        scale = np.sqrt(target_power / wsum)
        return weights * scale

    def receiver_coefficient(self, field: np.ndarray) -> complex:
        """Mode overlap coefficient for a sampled Hz line field."""
        return complex_overlap(self.normalized_field(), field)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "port_index": self.port_index,
            "label": self.label,
            "offsets_a": self.offsets_a.tolist(),
            "hz_complex": [
                {"real": z.real, "imag": z.imag} for z in self.hz_complex
            ],
            "span_a": self.span_a,
            "frequency_a": self.frequency_a,
            "res": self.res,
            "monitor_xy": list(self.monitor_xy),
            "tangent": list(self.tangent),
            "outward_dir": list(self.outward_dir),
            "analytic_te1_overlap": self.overlap_with_analytic_te1(),
        }

    def save(self, path: Path | str) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as f:
            json.dump(self.as_dict(), f, indent=2)
            f.write("\n")

    @classmethod
    def load(cls, path: Path | str) -> "NumericalPortMode":
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        hz = np.array(
            [complex(z["real"], z["imag"]) for z in data["hz_complex"]],
            dtype=complex,
        )
        return cls(
            port_index=int(data["port_index"]),
            label=str(data["label"]),
            offsets_a=np.asarray(data["offsets_a"], dtype=float),
            hz_complex=hz,
            span_a=float(data["span_a"]),
            frequency_a=float(data["frequency_a"]),
            res=int(data["res"]),
            monitor_xy=tuple(data["monitor_xy"]),
            tangent=tuple(data["tangent"]),
            outward_dir=tuple(data["outward_dir"]),
        )


def symmetry_equivalent_ports(port_index: int) -> Tuple[int, ...]:
    """Ports sharing magnitude of feed angle (+60° / −60° reuse one profile)."""
    if port_index in (1, 2):
        return (1, 2)
    if port_index in (4, 5):
        return (4, 5)
    return (port_index,)


def pick_reference_port(ports: Iterable[int]) -> int:
    """Choose canonical port for a symmetry class."""
    return min(ports)
