"""
Orientation-aware six-port source/monitor formulations for B=0 validation.

Isolated from production circulator physics. Registry entries swap only how
ports are excited and measured; geometry, rho, and B stay fixed.

Formulations
------------
baseline_hz_line
    Notebook path: Hz point line on true tangent + axis-aligned FluxRegion.

aligned_hz_line
    Hz points on the axis-aligned chord matching FluxRegion geometry.

te1_hz_line
    Hz on true tangent with TE1 cos amplitude + axis-aligned FluxRegion.

te1_guide_normal
    TE1 launch + guide-normal Poynting flux: point FluxRegions along the true
    feed tangent with weights (n_x, n_y) so Φ ≈ ∫ S·n̂ ds (same for all angles).

te1_guide_normal_te1w
    Same as te1_guide_normal but quadrature weights include a TE1 cos envelope
    (mode-weighted aperture flux; still not full modal decomposition).

eigenmode
    EigenModeSource candidate (failed on PEC horns; kept for reference only).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import meep as mp
import numpy as np

import sixport_common as sc


FORMULATION_NAMES = (
    "baseline_hz_line",
    "aligned_hz_line",
    "te1_hz_line",
    "te1_ez_line",
    "te1_guide_normal",
    "te1_guide_normal_te1w",
    "te1_guide_normal_dense",
    "te1_axis_at_source",
    "te1_dft_sdotn",
    "baseline_guide_normal",
    "eigenmode",
)

FluxSpec = Tuple[List[mp.FluxRegion], float]  # regions, outward_sign


def _unit(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("zero-length direction")
    return v / n


EIGENMODE_WALL_EPS = 1.0e6


def uses_finite_metal_walls(formulation: str) -> bool:
    return formulation == "eigenmode"


def axis_aligned_port_line(
    center_xy: np.ndarray,
    outward_dir: np.ndarray,
    span_factor: float = 0.98,
) -> Dict[str, Any]:
    """Same axis-aligned line construction as sixport_common.make_flux_region."""
    center_xy = np.asarray(center_xy, dtype=float)
    u = _unit(outward_dir)
    n = np.array([-u[1], u[0]])

    if abs(u[0]) >= abs(u[1]):
        line_dir = np.array([0.0, 1.0])
        alignment = abs(np.dot(line_dir, n))
        span = span_factor * sc.clear_width / alignment
        size = mp.Vector3(0, span, 0)
        direction = mp.X
        normal = np.array([1.0, 0.0])
    else:
        line_dir = np.array([1.0, 0.0])
        alignment = abs(np.dot(line_dir, n))
        span = span_factor * sc.clear_width / alignment
        size = mp.Vector3(span, 0, 0)
        direction = mp.Y
        normal = np.array([0.0, 1.0])

    outward_sign = float(np.sign(np.dot(u, normal)))
    return {
        "center": center_xy,
        "size": size,
        "direction": direction,
        "normal": normal,
        "line_dir": line_dir,
        "span": float(span),
        "outward_sign": outward_sign,
        "outward_dir": u,
    }


def _hz_points_along_line(
    center: np.ndarray,
    line_dir: np.ndarray,
    span: float,
    n_points: int,
    frequency: float,
    fwidth: float,
) -> List[mp.Source]:
    line_dir = _unit(line_dir)
    offsets = np.linspace(-span / 2, +span / 2, n_points)
    sources: List[mp.Source] = []
    for s in offsets:
        xy = center + s * line_dir
        sources.append(
            mp.Source(
                src=mp.GaussianSource(frequency=frequency, fwidth=fwidth),
                component=mp.Hz,
                center=mp.Vector3(xy[0], xy[1], 0),
                amplitude=1.0 / n_points,
            )
        )
    return sources


def make_baseline_sources(
    port_index: int,
    frequency: Optional[float] = None,
    fwidth: Optional[float] = None,
    n_points: int = 31,
) -> List[mp.Source]:
    return sc.make_port_source(
        port_index, frequency=frequency, fwidth=fwidth, n_points=n_points
    )


def make_aligned_hz_sources(
    port_index: int,
    frequency: Optional[float] = None,
    fwidth: Optional[float] = None,
    n_points: int = 31,
) -> List[mp.Source]:
    if frequency is None:
        frequency = sc.fs_a
    if fwidth is None:
        fwidth = sc.source_df
    center = np.asarray(sc.horn_for_port(port_index, sc.current_res())["source_center"], dtype=float)
    geom = axis_aligned_port_line(center, sc.effective_port_dir(port_index), span_factor=0.96)
    return _hz_points_along_line(
        center, geom["line_dir"], geom["span"], n_points, frequency, fwidth
    )


def make_te1_hz_sources(
    port_index: int,
    frequency: Optional[float] = None,
    fwidth: Optional[float] = None,
    n_points: int = 31,
) -> List[mp.Source]:
    if frequency is None:
        frequency = sc.fs_a
    if fwidth is None:
        fwidth = sc.source_df

    center = np.asarray(sc.horn_for_port(port_index, sc.current_res())["source_center"], dtype=float)
    u = _unit(sc.effective_port_dir(port_index))
    tangent = np.array([-u[1], u[0]])
    source_span = 0.96 * sc.clear_width
    offsets = np.linspace(-source_span / 2, +source_span / 2, n_points)
    weights = np.cos(np.pi * offsets / source_span)
    weights = np.clip(weights, 0.0, None)
    wsum = float(np.sum(weights))
    if wsum <= 0:
        weights = np.ones_like(weights)
        wsum = float(n_points)
    weights = weights / wsum

    sources: List[mp.Source] = []
    for s, amp in zip(offsets, weights):
        xy = center + s * tangent
        sources.append(
            mp.Source(
                src=mp.GaussianSource(frequency=frequency, fwidth=fwidth),
                component=mp.Hz,
                center=mp.Vector3(xy[0], xy[1], 0),
                amplitude=float(amp),
            )
        )
    return sources


def make_te1_ez_sources(
    port_index: int,
    frequency: Optional[float] = None,
    fwidth: Optional[float] = None,
    n_points: int = 31,
) -> List[mp.Source]:
    """Paper/Ceviche-equivalent: Ez polarized along discharge axis (2D: mp.Ez)."""
    if frequency is None:
        frequency = sc.fs_a
    if fwidth is None:
        fwidth = sc.source_df

    center = np.asarray(sc.horn_for_port(port_index, sc.current_res())["source_center"], dtype=float)
    u = _unit(sc.effective_port_dir(port_index))
    tangent = np.array([-u[1], u[0]])
    source_span = 0.96 * sc.clear_width
    offsets = np.linspace(-source_span / 2, +source_span / 2, n_points)
    weights = np.cos(np.pi * offsets / source_span)
    weights = np.clip(weights, 0.0, None)
    wsum = float(np.sum(weights))
    if wsum <= 0:
        weights = np.ones_like(weights)
        wsum = float(n_points)
    weights = weights / wsum

    sources: List[mp.Source] = []
    for s, amp in zip(offsets, weights):
        xy = center + s * tangent
        sources.append(
            mp.Source(
                src=mp.GaussianSource(frequency=frequency, fwidth=fwidth),
                component=mp.Ez,
                center=mp.Vector3(xy[0], xy[1], 0),
                amplitude=float(amp),
            )
        )
    return sources


def make_eigenmode_sources(
    port_index: int,
    frequency: Optional[float] = None,
    fwidth: Optional[float] = None,
) -> List[mp.Source]:
    if frequency is None:
        frequency = sc.fs_a
    if fwidth is None:
        fwidth = sc.source_df

    center = np.asarray(sc.horn_for_port(port_index, sc.current_res())["source_center"], dtype=float)
    u = _unit(sc.effective_port_dir(port_index))
    inward = -u
    geom = axis_aligned_port_line(center, u, span_factor=0.90)
    size = geom["size"]
    kpoint = mp.Vector3(inward[0], inward[1], 0)

    if abs(u[0]) >= 0.95:
        kwargs = dict(direction=mp.X, eig_kpoint=kpoint)
    elif abs(u[1]) >= 0.95:
        kwargs = dict(direction=mp.Y, eig_kpoint=kpoint)
    else:
        kwargs = dict(direction=mp.NO_DIRECTION, eig_kpoint=kpoint)

    src = mp.EigenModeSource(
        src=mp.GaussianSource(frequency=frequency, fwidth=fwidth),
        center=mp.Vector3(center[0], center[1], 0),
        size=size,
        eig_band=1,
        eig_parity=mp.ODD_Z,
        eig_match_freq=True,
        **kwargs,
    )
    return [src]


def make_axis_aligned_flux_spec(
    center_xy: np.ndarray,
    outward_dir: np.ndarray,
) -> FluxSpec:
    region, sign = sc.make_flux_region(center_xy, outward_dir)
    return [region], float(sign)


def make_guide_normal_flux_spec(
    center_xy: np.ndarray,
    outward_dir: np.ndarray,
    n_points: int = 31,
    span_factor: float = 0.96,
    weight_profile: str = "uniform",
) -> FluxSpec:
    """
    Approximate ∫ S·n̂ ds along the true feed cross-section.

    Meep FluxRegions are axis-aligned, so we place point monitors along the
    *guide tangent* and weight X/Y Poynting components by the outward normal
    components n̂ = (n_x, n_y). Positive flux = power in the +n̂ direction.
    """
    center_xy = np.asarray(center_xy, dtype=float)
    u = _unit(outward_dir)
    tangent = np.array([-u[1], u[0]])
    span = span_factor * sc.clear_width
    offsets = np.linspace(-span / 2, +span / 2, n_points)
    if n_points == 1:
        ds_weights = np.array([span])
    else:
        ds = span / (n_points - 1)
        ds_weights = np.full(n_points, ds)
        ds_weights[0] *= 0.5
        ds_weights[-1] *= 0.5

    if weight_profile == "te1":
        env = np.cos(np.pi * offsets / span)
        env = np.clip(env, 0.0, None)
        # Renormalize so total quadrature weight still equals physical span
        # (preserves units of power ≈ flux through the aperture).
        raw = ds_weights * env
        raw_sum = float(np.sum(raw))
        if raw_sum > 0:
            ds_weights = raw * (span / raw_sum)
        # else keep uniform ds_weights

    regions: List[mp.FluxRegion] = []
    for s, w_ds in zip(offsets, ds_weights):
        xy = center_xy + s * tangent
        c = mp.Vector3(xy[0], xy[1], 0)
        if abs(u[0]) > 1e-14:
            regions.append(
                mp.FluxRegion(
                    center=c,
                    size=mp.Vector3(),
                    direction=mp.X,
                    weight=float(u[0] * w_ds),
                )
            )
        if abs(u[1]) > 1e-14:
            regions.append(
                mp.FluxRegion(
                    center=c,
                    size=mp.Vector3(),
                    direction=mp.Y,
                    weight=float(u[1] * w_ds),
                )
            )

    if not regions:
        raise ValueError("guide-normal flux produced no regions")
    # Flux already oriented along outward n̂.
    return regions, 1.0


def port_measure_center(formulation: str, port_index: int) -> np.ndarray:
    """Geometric center used for flux/DFT monitors."""
    res = sc.current_res()
    if formulation == "te1_axis_at_source":
        return np.asarray(
            sc.horn_for_port(port_index, res)["source_center"], dtype=float
        )
    return sc.monitor_center_for_port(port_index, res)


def make_flux_region_for_formulation(
    formulation: str,
    center_xy: np.ndarray,
    outward_dir: np.ndarray,
) -> FluxSpec:
    if formulation in (
        "te1_guide_normal",
        "baseline_guide_normal",
    ):
        return make_guide_normal_flux_spec(
            center_xy, outward_dir, n_points=31, weight_profile="uniform"
        )
    if formulation == "te1_guide_normal_dense":
        return make_guide_normal_flux_spec(
            center_xy, outward_dir, n_points=61, weight_profile="uniform"
        )
    if formulation in ("te1_guide_normal_te1w",):
        return make_guide_normal_flux_spec(
            center_xy, outward_dir, weight_profile="te1"
        )
    return make_axis_aligned_flux_spec(center_xy, outward_dir)


def add_flux_monitor(sim, regions: Sequence[mp.FluxRegion]):
    """Register one or more FluxRegions as a single Meep flux object."""
    return sim.add_flux(sc.fs_a, 0, 1, *list(regions))


def add_dft_sdotn_monitor(sim, center_xy: np.ndarray, outward_dir: np.ndarray):
    """
    Axis-aligned DFT volume covering the rotated aperture chord.

    Power is extracted later by interpolating Ex,Ey,Hz along the true tangent
    and integrating (1/2) Re(E × H*) · n̂.
    """
    center_xy = np.asarray(center_xy, dtype=float)
    u = _unit(outward_dir)
    tangent = np.array([-u[1], u[0]])
    span = 0.96 * sc.clear_width
    # Bounding box of the rotated chord, padded by ~2 pixels worth of a-units
    # (caller should use a reasonable resolution; pad in absolute a-units).
    # Pad enough that axis-aligned ports (thin in the guide-normal direction)
    # still have several Yee cells for bilinear interpolation.
    pad = 0.20
    corners = np.array(
        [
            center_xy + 0.5 * span * tangent,
            center_xy - 0.5 * span * tangent,
        ]
    )
    xmin, ymin = corners.min(axis=0) - pad
    xmax, ymax = corners.max(axis=0) + pad
    # Enforce a minimum normal thickness so horizontal/vertical ports are not
    # a single-cell DFT slab.
    min_thick = 0.40
    if (xmax - xmin) < min_thick:
        mid = 0.5 * (xmin + xmax)
        xmin, xmax = mid - 0.5 * min_thick, mid + 0.5 * min_thick
    if (ymax - ymin) < min_thick:
        mid = 0.5 * (ymin + ymax)
        ymin, ymax = mid - 0.5 * min_thick, mid + 0.5 * min_thick
    vol = mp.Volume(
        center=mp.Vector3(0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0),
        size=mp.Vector3(xmax - xmin, ymax - ymin, 0),
    )
    dft = sim.add_dft_fields([mp.Ex, mp.Ey, mp.Hz], [sc.fs_a], where=vol)
    return {
        "dft": dft,
        "volume": vol,
        "center": center_xy,
        "outward": u,
        "tangent": tangent,
        "span": span,
    }


def _bilinear(arr: np.ndarray, x: np.ndarray, y: np.ndarray, xq: float, yq: float) -> complex:
    """Bilinear interpolate complex array on a rectilinear x,y grid (arr[ix, iy])."""
    if arr.ndim != 2:
        raise ValueError(f"expected 2d DFT array, got shape {arr.shape}")
    nx, ny = arr.shape
    if nx < 1 or ny < 1:
        return 0.0 + 0.0j
    if nx == 1 and ny == 1:
        return complex(arr[0, 0])
    if nx == 1:
        # Linear in y only
        if yq <= y[0]:
            return complex(arr[0, 0])
        if yq >= y[-1]:
            return complex(arr[0, -1])
        iy0 = int(np.searchsorted(y, yq) - 1)
        iy1 = iy0 + 1
        ty = 0.0 if y[iy1] == y[iy0] else (yq - y[iy0]) / (y[iy1] - y[iy0])
        return (1 - ty) * complex(arr[0, iy0]) + ty * complex(arr[0, iy1])
    if ny == 1:
        if xq <= x[0]:
            return complex(arr[0, 0])
        if xq >= x[-1]:
            return complex(arr[-1, 0])
        ix0 = int(np.searchsorted(x, xq) - 1)
        ix1 = ix0 + 1
        tx = 0.0 if x[ix1] == x[ix0] else (xq - x[ix0]) / (x[ix1] - x[ix0])
        return (1 - tx) * complex(arr[ix0, 0]) + tx * complex(arr[ix1, 0])
    if xq <= x[0]:
        ix0 = 0
    elif xq >= x[-1]:
        ix0 = nx - 2
    else:
        ix0 = int(np.searchsorted(x, xq) - 1)
    if yq <= y[0]:
        iy0 = 0
    elif yq >= y[-1]:
        iy0 = ny - 2
    else:
        iy0 = int(np.searchsorted(y, yq) - 1)
    ix1 = ix0 + 1
    iy1 = iy0 + 1
    tx = 0.0 if x[ix1] == x[ix0] else (xq - x[ix0]) / (x[ix1] - x[ix0])
    ty = 0.0 if y[iy1] == y[iy0] else (yq - y[iy0]) / (y[iy1] - y[iy0])
    a00 = arr[ix0, iy0]
    a10 = arr[ix1, iy0]
    a01 = arr[ix0, iy1]
    a11 = arr[ix1, iy1]
    return (1 - tx) * (1 - ty) * a00 + tx * (1 - ty) * a10 + (1 - tx) * ty * a01 + tx * ty * a11


def _dft_axes(ex: np.ndarray, mon_info: Dict[str, Any], sim) -> tuple:
    """Return 1d x,y axes matching ex.shape = (nx, ny)."""
    try:
        meta = sim.get_array_metadata(dft_cell=mon_info["dft"])
        x_coords = np.asarray(meta[0], dtype=float)
        y_coords = np.asarray(meta[1], dtype=float)
    except Exception:  # noqa: BLE001
        x_coords = np.array([])
        y_coords = np.array([])

    nx, ny = ex.shape
    # Prefer unique sorted axes when metadata is a mesh.
    if x_coords.size and y_coords.size:
        xu = np.unique(x_coords.ravel())
        yu = np.unique(y_coords.ravel())
        if len(xu) == nx and len(yu) == ny:
            return xu, yu
        if x_coords.ndim == 1 and len(x_coords) == nx and len(y_coords) == ny:
            return x_coords, y_coords

    vol = mon_info["volume"]
    c = vol.center
    s = vol.size
    x_1d = np.linspace(c.x - 0.5 * s.x, c.x + 0.5 * s.x, nx) if nx > 1 else np.array([c.x])
    y_1d = np.linspace(c.y - 0.5 * s.y, c.y + 0.5 * s.y, ny) if ny > 1 else np.array([c.y])
    return x_1d, y_1d


def extract_dft_sdotn_power(sim, mon_info: Dict[str, Any], n_points: int = 61) -> float:
    """Return signed outward power ≈ ∫ (1/2) Re(E×H*)·n̂ ds along the chord."""
    dft = mon_info["dft"]
    center = mon_info["center"]
    u = mon_info["outward"]
    tangent = mon_info["tangent"]
    span = mon_info["span"]

    ex = np.asarray(sim.get_dft_array(dft, mp.Ex, 0))
    ey = np.asarray(sim.get_dft_array(dft, mp.Ey, 0))
    hz = np.asarray(sim.get_dft_array(dft, mp.Hz, 0))
    if ex.ndim == 1:
        # Degenerate 1d DFT slab — reshape using volume aspect
        vol = mon_info["volume"]
        # Prefer interpreting as variation along the longer volume axis
        if vol.size.x >= vol.size.y:
            ex = ex.reshape((-1, 1))
            ey = ey.reshape((-1, 1))
            hz = hz.reshape((-1, 1))
        else:
            ex = ex.reshape((1, -1))
            ey = ey.reshape((1, -1))
            hz = hz.reshape((1, -1))
    x_1d, y_1d = _dft_axes(ex, mon_info, sim)

    offsets = np.linspace(-span / 2, span / 2, n_points)
    if n_points == 1:
        w = np.array([span])
    else:
        ds = span / (n_points - 1)
        w = np.full(n_points, ds)
        w[0] *= 0.5
        w[-1] *= 0.5

    total = 0.0
    for s, ws in zip(offsets, w):
        xy = center + s * tangent
        Ex = _bilinear(ex, x_1d, y_1d, float(xy[0]), float(xy[1]))
        Ey = _bilinear(ey, x_1d, y_1d, float(xy[0]), float(xy[1]))
        Hz = _bilinear(hz, x_1d, y_1d, float(xy[0]), float(xy[1]))
        Sx = 0.5 * np.real(Ey * np.conj(Hz))
        Sy = -0.5 * np.real(Ex * np.conj(Hz))
        total += (u[0] * Sx + u[1] * Sy) * ws
    return float(total)


def extract_flux_powers(sim, monitors, signs) -> np.ndarray:
    fluxes = np.zeros(len(monitors))
    for k, mon in enumerate(monitors):
        fluxes[k] = signs[k] * mp.get_fluxes(mon)[0]
    return fluxes


def eigenmode_kpoint(port_index: int) -> mp.Vector3:
    u = _unit(sc.effective_port_dir(port_index))
    inward = -u
    return mp.Vector3(inward[0], inward[1], 0)


def extract_eigenmode_powers(
    sim,
    monitors,
    port_indices: Sequence[int],
    signs: Sequence[float],
) -> np.ndarray:
    powers = np.zeros(len(monitors))
    for k, (mon, port_index, sign) in enumerate(zip(monitors, port_indices, signs)):
        kpt = eigenmode_kpoint(port_index)
        res = sim.get_eigenmode_coefficients(
            mon,
            [1],
            eig_parity=mp.ODD_Z,
            direction=mp.NO_DIRECTION,
            kpoint_func=lambda f, n, _k=kpt: _k,
        )
        a_plus = res.alpha[0, 0, 0]
        a_minus = res.alpha[0, 0, 1]
        p_plus = float(np.abs(a_plus) ** 2)
        p_minus = float(np.abs(a_minus) ** 2)
        if sign >= 0:
            powers[k] = p_plus - p_minus
        else:
            powers[k] = p_minus - p_plus
    return powers


@dataclass(frozen=True)
class Formulation:
    name: str
    measurement: str  # "flux" | "dft_sdotn" | "eigenmode"
    make_sources: Callable[..., List[mp.Source]]


_REGISTRY: Dict[str, Formulation] = {
    "baseline_hz_line": Formulation(
        name="baseline_hz_line",
        measurement="flux",
        make_sources=make_baseline_sources,
    ),
    "aligned_hz_line": Formulation(
        name="aligned_hz_line",
        measurement="flux",
        make_sources=make_aligned_hz_sources,
    ),
    "te1_hz_line": Formulation(
        name="te1_hz_line",
        measurement="flux",
        make_sources=make_te1_hz_sources,
    ),
    "te1_ez_line": Formulation(
        name="te1_ez_line",
        measurement="flux",
        make_sources=make_te1_ez_sources,
    ),
    "te1_guide_normal": Formulation(
        name="te1_guide_normal",
        measurement="flux",
        make_sources=make_te1_hz_sources,
    ),
    "te1_guide_normal_te1w": Formulation(
        name="te1_guide_normal_te1w",
        measurement="flux",
        make_sources=make_te1_hz_sources,
    ),
    "te1_guide_normal_dense": Formulation(
        name="te1_guide_normal_dense",
        measurement="flux",
        make_sources=make_te1_hz_sources,
    ),
    "te1_axis_at_source": Formulation(
        name="te1_axis_at_source",
        measurement="flux",
        make_sources=make_te1_hz_sources,
    ),
    "te1_dft_sdotn": Formulation(
        name="te1_dft_sdotn",
        measurement="dft_sdotn",
        make_sources=make_te1_hz_sources,
    ),
    "baseline_guide_normal": Formulation(
        name="baseline_guide_normal",
        measurement="flux",
        make_sources=make_baseline_sources,
    ),
    "eigenmode": Formulation(
        name="eigenmode",
        measurement="flux",
        make_sources=make_eigenmode_sources,
    ),
}


def get_formulation(name: str) -> Formulation:
    if name not in _REGISTRY:
        raise KeyError(
            f"Unknown port formulation {name!r}; choose from {sorted(_REGISTRY)}"
        )
    return _REGISTRY[name]


def list_formulations() -> List[str]:
    return list(FORMULATION_NAMES)
