"""
Full electromagnetic line-mode overlap for 2D Hz / Ex / Ey ports.

Implements a Meep-style modal bilinear form on a 1D port cross-section:

  <psi_m, psi> = ∫ [E_t,m* H_z + E_t H_z,m*] dl

with forward/backward amplitudes for a reciprocal uniform guide:

  alpha+ = 0.5 * (I1 + I2)
  alpha- = 0.5 * (I1 - I2)

where
  I1 = ∫ E_t,m* H_z dl
  I2 = ∫ E_t H_z,m* dl

The reference mode is normalized to unit modal power

  P = Re ∫ E_t* H_z dl = 1

so that |alpha±|^2 are mode powers in the same convention as Meep's
get_eigenmode_coefficients documentation.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Sequence

import meep as mp
import numpy as np


def trapezoid_weights(offsets: np.ndarray) -> np.ndarray:
    s = np.asarray(offsets, dtype=float).reshape(-1)
    if s.size == 0:
        return np.array([], dtype=float)
    if s.size == 1:
        return np.array([1.0], dtype=float)
    ds = np.diff(s)
    w = np.zeros_like(s)
    w[0] = 0.5 * ds[0]
    w[-1] = 0.5 * ds[-1]
    if s.size > 2:
        w[1:-1] = 0.5 * (ds[:-1] + ds[1:])
    return w


@dataclass(frozen=True)
class EMLineMode:
    port_index: int
    label: str
    offsets_a: np.ndarray
    ex_complex: np.ndarray
    ey_complex: np.ndarray
    hz_complex: np.ndarray
    center_xy: tuple[float, float]
    tangent_xy: tuple[float, float]
    outward_dir: tuple[float, float]
    mode_power: float

    def tangent(self) -> np.ndarray:
        return np.asarray(self.tangent_xy, dtype=float)

    def outward(self) -> np.ndarray:
        return np.asarray(self.outward_dir, dtype=float)

    def et_complex(self) -> np.ndarray:
        t = self.tangent()
        return t[0] * self.ex_complex + t[1] * self.ey_complex

    def quadrature_weights(self) -> np.ndarray:
        return trapezoid_weights(self.offsets_a)

    def power_from_fields(self) -> float:
        et = self.et_complex()
        hz = np.asarray(self.hz_complex, dtype=complex)
        w = self.quadrature_weights()
        return float(np.real(np.sum(et * np.conj(hz) * w)))

    def assert_finite(self) -> None:
        for name, arr in (
            ("Ex", self.ex_complex),
            ("Ey", self.ey_complex),
            ("Hz", self.hz_complex),
        ):
            if not np.all(np.isfinite(arr)):
                raise ValueError(f"{self.label}: non-finite {name} samples")
        if not np.isfinite(self.mode_power):
            raise ValueError(f"{self.label}: non-finite mode_power")
        if self.mode_power <= 0:
            raise ValueError(f"{self.label}: non-positive mode_power={self.mode_power}")
        if self.mode_power > 1e12:
            raise ValueError(f"{self.label}: absurdly large mode_power={self.mode_power}")

    def normalized_to_unit_power(self) -> "EMLineMode":
        self.assert_finite()
        scale = np.sqrt(self.mode_power)
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError(f"{self.label}: invalid scale={scale}")
        out = EMLineMode(
            port_index=self.port_index,
            label=self.label,
            offsets_a=np.asarray(self.offsets_a, dtype=float),
            ex_complex=np.asarray(self.ex_complex, dtype=complex) / scale,
            ey_complex=np.asarray(self.ey_complex, dtype=complex) / scale,
            hz_complex=np.asarray(self.hz_complex, dtype=complex) / scale,
            center_xy=self.center_xy,
            tangent_xy=self.tangent_xy,
            outward_dir=self.outward_dir,
            mode_power=float(self.mode_power / (scale * scale)),
        )
        # Recompute after normalization to reduce drift.
        return out.with_recomputed_power()

    def with_recomputed_power(self) -> "EMLineMode":
        p = self.power_from_fields()
        return EMLineMode(
            port_index=self.port_index,
            label=self.label,
            offsets_a=np.asarray(self.offsets_a, dtype=float),
            ex_complex=np.asarray(self.ex_complex, dtype=complex),
            ey_complex=np.asarray(self.ey_complex, dtype=complex),
            hz_complex=np.asarray(self.hz_complex, dtype=complex),
            center_xy=self.center_xy,
            tangent_xy=self.tangent_xy,
            outward_dir=self.outward_dir,
            mode_power=float(p),
        )

    def overlap_amplitudes(
        self,
        *,
        ex_field: Sequence[complex],
        ey_field: Sequence[complex],
        hz_field: Sequence[complex],
    ) -> Dict[str, complex | float]:
        ex = np.asarray(ex_field, dtype=complex).reshape(-1)
        ey = np.asarray(ey_field, dtype=complex).reshape(-1)
        hz = np.asarray(hz_field, dtype=complex).reshape(-1)
        if ex.shape != self.ex_complex.shape or ey.shape != self.ey_complex.shape or hz.shape != self.hz_complex.shape:
            raise ValueError("field/mode sample shape mismatch")
        t = self.tangent()
        et_mode = self.et_complex()
        et_field = t[0] * ex + t[1] * ey
        w = self.quadrature_weights()
        i1 = np.sum(np.conj(et_mode) * hz * w)
        i2 = np.sum(et_field * np.conj(self.hz_complex) * w)
        alpha_plus = 0.5 * (i1 + i2)
        alpha_minus = 0.5 * (i1 - i2)
        return {
            "alpha_plus": complex(alpha_plus),
            "alpha_minus": complex(alpha_minus),
            "power_plus": float(np.abs(alpha_plus) ** 2),
            "power_minus": float(np.abs(alpha_minus) ** 2),
            "i1": complex(i1),
            "i2": complex(i2),
        }

    def as_dict(self) -> Dict[str, Any]:
        def _carr(arr: np.ndarray) -> List[Dict[str, float]]:
            return [{"real": float(z.real), "imag": float(z.imag)} for z in arr]

        return {
            "port_index": self.port_index,
            "label": self.label,
            "offsets_a": self.offsets_a.tolist(),
            "center_xy": list(self.center_xy),
            "tangent_xy": list(self.tangent_xy),
            "outward_dir": list(self.outward_dir),
            "mode_power": float(self.mode_power),
            "ex_complex": _carr(np.asarray(self.ex_complex, dtype=complex)),
            "ey_complex": _carr(np.asarray(self.ey_complex, dtype=complex)),
            "hz_complex": _carr(np.asarray(self.hz_complex, dtype=complex)),
        }

    def save(self, path: Path | str) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as f:
            json.dump(self.as_dict(), f, indent=2)
            f.write("\n")

    @classmethod
    def load(cls, path: Path | str) -> "EMLineMode":
        data = json.loads(Path(path).read_text(encoding="utf-8"))

        def _carr(items: Sequence[Dict[str, float]]) -> np.ndarray:
            return np.array(
                [complex(z["real"], z["imag"]) for z in items], dtype=complex
            )

        return cls(
            port_index=int(data["port_index"]),
            label=str(data["label"]),
            offsets_a=np.asarray(data["offsets_a"], dtype=float),
            ex_complex=_carr(data["ex_complex"]),
            ey_complex=_carr(data["ey_complex"]),
            hz_complex=_carr(data["hz_complex"]),
            center_xy=tuple(float(x) for x in data["center_xy"]),
            tangent_xy=tuple(float(x) for x in data["tangent_xy"]),
            outward_dir=tuple(float(x) for x in data["outward_dir"]),
            mode_power=float(data["mode_power"]),
        )


def add_em_line_monitor(
    sim: mp.Simulation,
    *,
    center_xy: np.ndarray,
    tangent_xy: np.ndarray,
    offsets_a: Sequence[float],
    frequency: float,
) -> Dict[str, Any]:
    center = np.asarray(center_xy, dtype=float).reshape(2)
    tangent = np.asarray(tangent_xy, dtype=float).reshape(2)
    tangent = tangent / max(np.linalg.norm(tangent), 1e-30)
    offsets = np.asarray(offsets_a, dtype=float)
    dft_objs = []
    for s in offsets:
        xy = center + float(s) * tangent
        dft_objs.append(
            sim.add_dft_fields(
                [mp.Ex, mp.Ey, mp.Hz],
                frequency,
                0,
                1,
                center=mp.Vector3(float(xy[0]), float(xy[1]), 0),
                size=mp.Vector3(),
            )
        )
    return {
        "dft_objs": dft_objs,
        "center_xy": center,
        "tangent_xy": tangent,
        "offsets_a": offsets,
        "frequency": frequency,
    }


def sample_em_line_fields(sim: mp.Simulation, mon_info: Dict[str, Any]) -> Dict[str, np.ndarray]:
    ex = []
    ey = []
    hz = []
    for dft in mon_info["dft_objs"]:
        ex.append(complex(np.squeeze(sim.get_dft_array(dft, mp.Ex, 0))))
        ey.append(complex(np.squeeze(sim.get_dft_array(dft, mp.Ey, 0))))
        hz.append(complex(np.squeeze(sim.get_dft_array(dft, mp.Hz, 0))))
    return {
        "ex": np.asarray(ex, dtype=complex),
        "ey": np.asarray(ey, dtype=complex),
        "hz": np.asarray(hz, dtype=complex),
    }


def build_em_mode_from_samples(
    *,
    port_index: int,
    label: str,
    center_xy: np.ndarray,
    tangent_xy: np.ndarray,
    outward_dir: np.ndarray,
    offsets_a: Sequence[float],
    ex: Sequence[complex],
    ey: Sequence[complex],
    hz: Sequence[complex],
) -> EMLineMode:
    mode = EMLineMode(
        port_index=int(port_index),
        label=str(label),
        offsets_a=np.asarray(offsets_a, dtype=float),
        ex_complex=np.asarray(ex, dtype=complex),
        ey_complex=np.asarray(ey, dtype=complex),
        hz_complex=np.asarray(hz, dtype=complex),
        center_xy=(float(center_xy[0]), float(center_xy[1])),
        tangent_xy=(float(tangent_xy[0]), float(tangent_xy[1])),
        outward_dir=(float(outward_dir[0]), float(outward_dir[1])),
        mode_power=0.0,
    ).with_recomputed_power()
    return mode
