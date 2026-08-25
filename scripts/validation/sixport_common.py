"""
Minimal six-port PMM circulator forward model for B=0 reciprocity tests.

Extracted from scripts/Sketchbook_PMMCirculator.ipynb (standardized geometry /
source / simulate / normalize path). Physics is kept faithful to the notebook;
this module does not attempt to "fix" reciprocity residuals.
"""

from __future__ import annotations

import os
import pickle
import sys
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import meep as mp
import numpy as np
from scipy.spatial import ConvexHull

# Import PMMI from scripts/
_SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
_ROOT = os.path.abspath(os.path.join(_SCRIPTS_DIR, ".."))
for _p in (_ROOT, _SCRIPTS_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from PMMCirculatorInverse import PMMI  # noqa: E402

# ---------------------------------------------------------------------------
# Physical / numerical constants (notebook cells 4, 22, 23, 27, 28)
# ---------------------------------------------------------------------------

a = 0.028  # 1 Meep length unit = 2.8 cm
dpml = 1.0

fs_Hz = 3.85e9
fp_Hz = 8.00e9
gamma_Hz = 1.00e6

# Nondimensional frequencies via PMMI (same c/a scaling as notebook)
_freq_helper = PMMI(a=a, res=32, nx=10, ny=10, dpml=1)
fs_a = _freq_helper.f_a(fs_Hz / 1e9)
fp_a = _freq_helper.f_a(fp_Hz / 1e9)
gamma_a = _freq_helper.f_a(gamma_Hz / 1e9)

source_df = 0.10 * fs_a

# Experimental bulb / array dimensions (cell 22)
d_exp = 0.020 / a  # 20 mm center-to-center
r_plasma = 0.005 / a  # ~5 mm plasma radius
r_bulb_inner = 0.0065 / a  # 13 mm ID
r_bulb_outer = 0.0075 / a  # 15 mm OD

# Horn dimensions from PMMInverse TM orientation (cell 23)
wall_thickness = 0.004 / a  # 4 mm
width_open = 0.104 / a  # 104 mm aperture
width_base = 0.048 / a  # 48 mm throat
horn_depth = 0.089 / a  # 89 mm flare depth

# Numerical straight-feed length (cell 27)
feed_length_m = 0.060
feed_length = feed_length_m / a

clear_width = width_base - 2 * wall_thickness  # cell 28

# In-memory normalization cache: (res, run_time) -> dict
_NORM_CACHE: Dict[Tuple[int, float], Dict[str, Any]] = {}


def default_uniform_rho() -> np.ndarray:
    """91-vector with every element at fp corresponding to 8 GHz."""
    return np.ones(91, dtype=float) * fp_a


# ---------------------------------------------------------------------------
# Hexagon faces / outward port directions (cells 21–22), then full horns (27)
# ---------------------------------------------------------------------------

def _compute_port_geometry_from_locs(locs: np.ndarray):
    """Outward port normals ordered CCW starting near +x."""
    locs = np.asarray(locs, dtype=float)
    array_center = np.mean(locs, axis=0)

    hull = ConvexHull(locs)
    corners = locs[hull.vertices]
    assert len(corners) == 6

    face_centers = []
    port_dirs = []
    for i in range(6):
        p1 = corners[i]
        p2 = corners[(i + 1) % 6]
        face_center = (p1 + p2) / 2
        tangent = p2 - p1
        tangent = tangent / np.linalg.norm(tangent)
        normal = np.array([tangent[1], -tangent[0]])
        if np.dot(normal, face_center - array_center) < 0:
            normal = -normal
        face_centers.append(face_center)
        port_dirs.append(normal)

    face_centers = np.asarray(face_centers, dtype=float)
    port_dirs = np.asarray(port_dirs, dtype=float)

    angles = np.mod(np.arctan2(port_dirs[:, 1], port_dirs[:, 0]), 2 * np.pi)
    order = np.argsort(angles)
    face_centers = face_centers[order]
    port_dirs = port_dirs[order]
    return array_center, face_centers, port_dirs


def get_full_horn(xy_open_cen, outward_dir, feed_length_local):
    """Complete horn wall polygons + source/monitor centers (cell 27)."""
    outward_dir = np.asarray(outward_dir, dtype=float)
    outward_dir = outward_dir / np.linalg.norm(outward_dir)
    horn_dir = -outward_dir
    horn_dir_orth = np.array([horn_dir[1], -horn_dir[0]])

    left_open = np.array(
        [
            xy_open_cen + horn_dir_orth * width_open / 2,
            xy_open_cen + horn_dir_orth * (width_open / 2 - wall_thickness),
            xy_open_cen
            + horn_dir_orth * (width_base / 2 - wall_thickness)
            - horn_dir * horn_depth,
            xy_open_cen + horn_dir_orth * width_base / 2 - horn_dir * horn_depth,
        ]
    )
    left_feed = np.array(
        [
            xy_open_cen
            + horn_dir_orth * (width_base / 2 - wall_thickness)
            - horn_dir * horn_depth,
            xy_open_cen + horn_dir_orth * width_base / 2 - horn_dir * horn_depth,
            xy_open_cen
            + horn_dir_orth * width_base / 2
            - horn_dir * (horn_depth + feed_length_local),
            xy_open_cen
            + horn_dir_orth * (width_base / 2 - wall_thickness)
            - horn_dir * (horn_depth + feed_length_local),
        ]
    )
    right_open = np.array(
        [
            xy_open_cen - horn_dir_orth * width_open / 2,
            xy_open_cen - horn_dir_orth * (width_open / 2 - wall_thickness),
            xy_open_cen
            - horn_dir_orth * (width_base / 2 - wall_thickness)
            - horn_dir * horn_depth,
            xy_open_cen - horn_dir_orth * width_base / 2 - horn_dir * horn_depth,
        ]
    )
    right_feed = np.array(
        [
            xy_open_cen
            - horn_dir_orth * (width_base / 2 - wall_thickness)
            - horn_dir * horn_depth,
            xy_open_cen - horn_dir_orth * width_base / 2 - horn_dir * horn_depth,
            xy_open_cen
            - horn_dir_orth * width_base / 2
            - horn_dir * (horn_depth + feed_length_local),
            xy_open_cen
            - horn_dir_orth * (width_base / 2 - wall_thickness)
            - horn_dir * (horn_depth + feed_length_local),
        ]
    )

    throat_center = xy_open_cen + outward_dir * horn_depth
    monitor_center = throat_center + outward_dir * (0.30 * feed_length_local)
    source_center = throat_center + outward_dir * (0.70 * feed_length_local)

    return {
        "left_flare": left_open,
        "right_flare": right_open,
        "left_feed": left_feed,
        "right_feed": right_feed,
        "throat_center": throat_center,
        "monitor_center": monitor_center,
        "source_center": source_center,
    }


def _init_canonical_geometry():
    """Build port_dirs, full_horns, and six-port domain sizes once."""
    # Temporary array placement (cell 22) — only locations matter for faces.
    pmm_tmp = PMMI(
        a=a,
        res=32,
        nx=14,
        ny=12,
        dpml=dpml,
        B=np.array([0.0, 0.0, 0.0]),
    )
    pmm_tmp.Rod_Array_Hexagon_train(
        xy_cen=np.array([14 / 2, 12 / 2]),
        side_dim=6,
        r=r_plasma,
        d=d_exp,
        bulbs=False,
        uniform=True,
    )
    locs = np.asarray(pmm_tmp.train_elem_locs, dtype=float)
    assert len(locs) == 91

    array_center, face_centers, port_dirs = _compute_port_geometry_from_locs(locs)

    # Horns in array-centered coordinates (cell 27)
    face_offsets = face_centers - array_center
    horn_open_offsets = face_offsets + r_bulb_outer * port_dirs

    full_horns = [
        get_full_horn(horn_open_offsets[i], port_dirs[i], feed_length)
        for i in range(6)
    ]

    all_wall_points = np.vstack(
        [
            horn[name]
            for horn in full_horns
            for name in ("left_flare", "right_flare", "left_feed", "right_feed")
        ]
    )
    max_x = np.max(np.abs(all_wall_points[:, 0]))
    max_y = np.max(np.abs(all_wall_points[:, 1]))

    dpml_ports = 2 * dpml
    air_buffer = 0.5
    nx_ports = int(np.ceil(2 * (max_x + dpml_ports + air_buffer)))
    ny_ports = int(np.ceil(2 * (max_y + dpml_ports + air_buffer)))

    port_monitor_centers = np.array(
        [horn["monitor_center"] for horn in full_horns], dtype=float
    )

    return {
        "array_center_ref": array_center,
        "face_centers_ref": face_centers,
        "port_dirs": port_dirs,
        "full_horns": full_horns,
        "dpml_ports": dpml_ports,
        "nx_ports": nx_ports,
        "ny_ports": ny_ports,
        "port_monitor_centers": port_monitor_centers,
    }


_GEO = _init_canonical_geometry()
port_dirs = _GEO["port_dirs"]
full_horns = _GEO["full_horns"]
dpml_ports = _GEO["dpml_ports"]
nx_ports = _GEO["nx_ports"]
ny_ports = _GEO["ny_ports"]
port_monitor_centers = _GEO["port_monitor_centers"]


# ---------------------------------------------------------------------------
# Geometry helpers (cells 28, 35, 36)
# ---------------------------------------------------------------------------

def make_flux_region(center_xy, outward_dir):
    """Axis-aligned flux plane; positive sign = outward from PMM (cell 28)."""
    u = np.asarray(outward_dir, dtype=float)
    u = u / np.linalg.norm(u)
    n = np.array([-u[1], u[0]])

    if abs(u[0]) >= abs(u[1]):
        line_dir = np.array([0.0, 1.0])
        alignment = abs(np.dot(line_dir, n))
        span = 0.98 * clear_width / alignment
        region = mp.FluxRegion(
            center=mp.Vector3(center_xy[0], center_xy[1], 0),
            size=mp.Vector3(0, span, 0),
            direction=mp.X,
        )
        normal = np.array([1.0, 0.0])
    else:
        line_dir = np.array([1.0, 0.0])
        alignment = abs(np.dot(line_dir, n))
        span = 0.98 * clear_width / alignment
        region = mp.FluxRegion(
            center=mp.Vector3(center_xy[0], center_xy[1], 0),
            size=mp.Vector3(span, 0, 0),
            direction=mp.Y,
        )
        normal = np.array([0.0, 1.0])

    outward_sign = np.sign(np.dot(u, normal))
    return region, outward_sign


def make_port_source(
    port_index,
    frequency=None,
    fwidth=None,
    n_points=31,
):
    """Hz point-source line across any port feed (cell 35; all six ports)."""
    if frequency is None:
        frequency = fs_a
    if fwidth is None:
        fwidth = source_df

    center = np.asarray(full_horns[port_index]["source_center"], dtype=float)
    u = np.asarray(port_dirs[port_index], dtype=float)
    u = u / np.linalg.norm(u)
    tangent = np.array([-u[1], u[0]])
    source_span = 0.96 * clear_width
    offsets = np.linspace(-source_span / 2, +source_span / 2, n_points)

    sources = []
    for s in offsets:
        xy = center + s * tangent
        sources.append(
            mp.Source(
                src=mp.GaussianSource(frequency=frequency, fwidth=fwidth),
                component=mp.Hz,
                center=mp.Vector3(xy[0], xy[1], 0),
                amplitude=1.0 / n_points,
            )
        )
    return sources


def rotated_wall(centerline, axis, normal, normal_offset, length, thickness):
    """PEC wall rectangle for straight-feed reference guides (cell 36)."""
    wall_center = centerline + normal_offset * normal
    p1 = wall_center - axis * length / 2 - normal * thickness / 2
    p2 = wall_center + axis * length / 2 - normal * thickness / 2
    p3 = wall_center + axis * length / 2 + normal * thickness / 2
    p4 = wall_center - axis * length / 2 + normal * thickness / 2
    vertices = np.zeros((4, 3))
    vertices[:, 0:2] = np.array([p1, p2, p3, p4])
    return vertices


# ---------------------------------------------------------------------------
# Device builder (cell 34) — resolution is an argument
# ---------------------------------------------------------------------------

def build_circulator_device(rho, B, res=64):
    """Build six-port PlasMEEP geometry; no sources/monitors/run."""
    rho = np.asarray(rho, dtype=float).flatten()
    B = np.asarray(B, dtype=float).flatten()
    if len(rho) != 91:
        raise ValueError(f"Expected 91 rho values, got {len(rho)}.")
    if len(B) != 3:
        raise ValueError(f"B must contain [Bx, By, Bz], got {B}.")

    pmm_device = PMMI(
        a=a,
        res=res,
        nx=nx_ports,
        ny=ny_ports,
        dpml=dpml_ports,
        B=B,
    )
    array_center_device = np.array([nx_ports / 2, ny_ports / 2])
    pmm_device.Rod_Array_Hexagon_train(
        xy_cen=array_center_device,
        side_dim=6,
        r=r_plasma,
        d=d_exp,
        bulbs=False,
        uniform=True,
    )
    assert len(pmm_device.train_elems) == 91

    wp_values, _elem_locations = pmm_device.Scale_Rho_wp(
        rho, w_src=fs_a, wp_max=0, gamma=gamma_a
    )
    assert len(wp_values) == 91

    P_device = pmm_device.Build_Sim()
    locs_device = np.asarray(pmm_device.train_elem_locs, dtype=float)
    for i in range(91):
        center_meep = np.array(
            [
                locs_device[i, 0] - nx_ports / 2,
                locs_device[i, 1] - ny_ports / 2,
                0.0,
            ]
        )
        P_device.Add_Bulb(
            r_bulb=(r_bulb_inner, r_bulb_outer),
            center=center_meep,
            wp=wp_values[i],
            gamma=gamma_a,
            axis=np.array([0, 0, 1]),
            profile=0,
        )

    for horn in full_horns:
        for name in ("left_flare", "right_flare", "left_feed", "right_feed"):
            wall = horn[name]
            vertices = np.zeros((4, 3))
            vertices[:, 0] = wall[:, 0]
            vertices[:, 1] = wall[:, 1]
            P_device.Add_Prism(
                vertices=vertices,
                axis=np.array([0, 0, 1]),
                PEC=True,
            )

    return pmm_device, P_device, wp_values


# ---------------------------------------------------------------------------
# Circulator objective (cell 49) — used by simulate_circulator when 6x6
# ---------------------------------------------------------------------------

def circulator_objective(
    power_matrix,
    direction="CCW",
    w_reverse=1.0,
    w_reflection=0.5,
    w_leakage=0.25,
):
    power_matrix = np.asarray(power_matrix, dtype=float)
    if power_matrix.shape != (6, 6):
        raise ValueError("power_matrix must have shape (6, 6)")

    direction = direction.upper()
    if direction not in ("CCW", "CW"):
        raise ValueError("direction must be 'CCW' or 'CW'")

    desired_power = np.zeros(6)
    reverse_power = np.zeros(6)
    reflection_power = np.zeros(6)
    leakage_power = np.zeros(6)

    for input_port in range(6):
        if direction == "CCW":
            desired_output = (input_port + 1) % 6
            reverse_output = (input_port - 1) % 6
        else:
            desired_output = (input_port - 1) % 6
            reverse_output = (input_port + 1) % 6

        desired_power[input_port] = power_matrix[desired_output, input_port]
        reverse_power[input_port] = power_matrix[reverse_output, input_port]
        reflection_power[input_port] = power_matrix[input_port, input_port]
        leakage_outputs = [
            op
            for op in range(6)
            if op not in (desired_output, reverse_output, input_port)
        ]
        leakage_power[input_port] = np.sum(power_matrix[leakage_outputs, input_port])

    desired_mean = np.mean(desired_power)
    reverse_mean = np.mean(reverse_power)
    reflection_mean = np.mean(reflection_power)
    leakage_mean = np.mean(leakage_power)
    objective = (
        desired_mean
        - w_reverse * reverse_mean
        - w_reflection * reflection_mean
        - w_leakage * leakage_mean
    )

    isolation_dB = np.full(6, np.nan)
    for i in range(6):
        if desired_power[i] > 0 and reverse_power[i] > 0:
            isolation_dB[i] = 10 * np.log10(desired_power[i] / reverse_power[i])

    return {
        "objective": objective,
        "desired_power": desired_power,
        "reverse_power": reverse_power,
        "reflection_power": reflection_power,
        "leakage_power": leakage_power,
        "desired_mean": desired_mean,
        "reverse_mean": reverse_mean,
        "reflection_mean": reflection_mean,
        "leakage_mean": leakage_mean,
        "isolation_dB": isolation_dB,
        "mean_isolation_dB": np.nanmean(isolation_dB),
    }


# ---------------------------------------------------------------------------
# Normalization (cell 52) — make_port_source for every port
# ---------------------------------------------------------------------------

def normalize_port(source_port, res=64, run_time=80, verbose=True):
    """Straight-feed reference run for one port; returns (power, flux_data)."""
    if verbose:
        print()
        print("=" * 54)
        print(f"NORMALIZING PORT P{source_port + 1}")
        print("=" * 54)

    u = np.asarray(port_dirs[source_port], dtype=float)
    u = u / np.linalg.norm(u)
    n = np.array([-u[1], u[0]])

    pmm_ref = PMMI(
        a=a,
        res=res,
        nx=nx_ports,
        ny=ny_ports,
        dpml=dpml_ports,
        B=np.array([0.0, 0.0, 0.0]),
    )
    P_ref = pmm_ref.Build_Sim()

    reference_center = np.asarray(
        full_horns[source_port]["source_center"], dtype=float
    )
    guide_length = 2.0 * np.hypot(nx_ports, ny_ports)
    wall_center_offset = clear_width / 2 + wall_thickness / 2

    wall_A = rotated_wall(
        centerline=reference_center,
        axis=u,
        normal=n,
        normal_offset=+wall_center_offset,
        length=guide_length,
        thickness=wall_thickness,
    )
    wall_B = rotated_wall(
        centerline=reference_center,
        axis=u,
        normal=n,
        normal_offset=-wall_center_offset,
        length=guide_length,
        thickness=wall_thickness,
    )
    P_ref.Add_Prism(vertices=wall_A, axis=np.array([0, 0, 1]), PEC=True)
    P_ref.Add_Prism(vertices=wall_B, axis=np.array([0, 0, 1]), PEC=True)

    P_ref.sources = make_port_source(source_port)
    ref_sim = P_ref.Get_Sim()

    region, sign = make_flux_region(
        port_monitor_centers[source_port], port_dirs[source_port]
    )
    monitor = ref_sim.add_flux(fs_a, 0, 1, region)

    if verbose:
        print("Running reference...")
    ref_sim.run(until_after_sources=run_time)

    raw_flux = mp.get_fluxes(monitor)[0]
    signed_flux = sign * raw_flux
    incident_power = abs(signed_flux)
    incident_data = ref_sim.get_flux_data(monitor)

    if verbose:
        print("signed incident flux =", signed_flux)
        print("incident power =", incident_power)
        if signed_flux < 0:
            print("Direction check: PASS")
        else:
            print("WARNING: expected inward incident flux.")

    return incident_power, incident_data


def ensure_normalizations(
    res,
    run_time=80,
    force=False,
    ports: Optional[Sequence[int]] = None,
    cache_path: Optional[str] = None,
    verbose=True,
) -> Dict[str, Any]:
    """
    Return incident powers and flux-subtraction data for the requested ports.

    Caches in memory by (res, run_time). Optional pickle at cache_path for
    cross-process reuse (--skip-norm-if-cached).
    """
    key = (int(res), float(run_time))
    port_list = list(range(6) if ports is None else ports)

    if cache_path and (not force) and os.path.isfile(cache_path):
        with open(cache_path, "rb") as f:
            loaded = pickle.load(f)
        if (
            loaded.get("res") == int(res)
            and loaded.get("run_time") == float(run_time)
            and all(p in loaded["incident_power_by_port"] for p in port_list)
        ):
            _NORM_CACHE[key] = loaded
            if verbose:
                print(f"Loaded normalizations from {cache_path}")
            return loaded

    if (not force) and key in _NORM_CACHE:
        cached = _NORM_CACHE[key]
        if all(p in cached["incident_power_by_port"] for p in port_list):
            return cached

    # Start from existing cache entry if present (partial fill)
    if key in _NORM_CACHE and not force:
        incident_power_by_port = dict(_NORM_CACHE[key]["incident_power_by_port"])
        incident_flux_data_by_port = dict(
            _NORM_CACHE[key]["incident_flux_data_by_port"]
        )
    else:
        incident_power_by_port = {}
        incident_flux_data_by_port = {}

    for source_port in port_list:
        if (not force) and source_port in incident_power_by_port:
            continue
        power_i, data_i = normalize_port(
            source_port, res=res, run_time=run_time, verbose=verbose
        )
        incident_power_by_port[source_port] = power_i
        incident_flux_data_by_port[source_port] = data_i

    result = {
        "res": int(res),
        "run_time": float(run_time),
        "incident_power_by_port": incident_power_by_port,
        "incident_flux_data_by_port": incident_flux_data_by_port,
    }
    _NORM_CACHE[key] = result

    if cache_path:
        os.makedirs(os.path.dirname(os.path.abspath(cache_path)) or ".", exist_ok=True)
        with open(cache_path, "wb") as f:
            pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)
        if verbose:
            print(f"Wrote normalizations to {cache_path}")

    return result


# ---------------------------------------------------------------------------
# Forward model (cell 50) + reciprocity metrics (cell 40 style)
# ---------------------------------------------------------------------------

def reciprocity_metrics(power_matrix_dB, port_ids=None) -> Dict[str, Any]:
    """Mean/max |Pij - Pji| in dB and worst unique port pairs."""
    M = np.asarray(power_matrix_dB, dtype=float)
    n = M.shape[0]
    if port_ids is None:
        port_ids = list(range(n))
    else:
        port_ids = list(port_ids)

    err = np.abs(M - M.T)
    np.fill_diagonal(err, np.nan)
    finite = err[np.isfinite(err)]
    mean_abs = float(np.mean(finite)) if finite.size else float("nan")
    max_abs = float(np.max(finite)) if finite.size else float("nan")

    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            d = err[i, j]
            pairs.append(
                {
                    "i": int(port_ids[i]),
                    "j": int(port_ids[j]),
                    "label": f"P{port_ids[i] + 1}<->P{port_ids[j] + 1}",
                    "abs_diff_dB": None if not np.isfinite(d) else float(d),
                    "T_ji_dB": float(M[i, j]) if np.isfinite(M[i, j]) else None,
                    "T_ij_dB": float(M[j, i]) if np.isfinite(M[j, i]) else None,
                }
            )
    pairs_sorted = sorted(
        [p for p in pairs if p["abs_diff_dB"] is not None],
        key=lambda p: p["abs_diff_dB"],
        reverse=True,
    )
    return {
        "mean_abs_diff_dB": mean_abs,
        "max_abs_diff_dB": max_abs,
        "worst_pairs": pairs_sorted[:5],
        "all_pairs": pairs,
        "error_matrix_dB": err,
    }


def simulate_circulator(
    rho,
    B,
    res=64,
    run_time=80,
    direction="CCW",
    verbose=True,
    incident_cache=None,
    ports: Optional[Sequence[int]] = None,
):
    """
    Full six-port (or port-subset) forward model.

    Matrix convention: row = OUTPUT, column = INPUT.
    Uses make_port_source for every excited port (cell 52 standardization).
    """
    rho = np.asarray(rho, dtype=float).flatten()
    B = np.asarray(B, dtype=float).flatten()
    if len(rho) != 91:
        raise ValueError(f"Expected 91 rho values, got {len(rho)}.")
    if len(B) != 3:
        raise ValueError("B must be [Bx, By, Bz].")

    port_list = list(range(6) if ports is None else ports)
    n_ports = len(port_list)
    # Map global port index -> local matrix index
    local_index = {p: k for k, p in enumerate(port_list)}

    if incident_cache is None:
        incident_cache = ensure_normalizations(
            res, run_time=run_time, force=False, ports=port_list, verbose=verbose
        )
    incident_power_by_port = incident_cache["incident_power_by_port"]
    incident_flux_data_by_port = incident_cache["incident_flux_data_by_port"]
    for p in port_list:
        if p not in incident_power_by_port or p not in incident_flux_data_by_port:
            raise KeyError(
                f"incident_cache missing port {p}; call ensure_normalizations first."
            )

    if verbose:
        print("=" * 54)
        print("SIX-PORT CIRCULATOR FORWARD MODEL")
        print("=" * 54)
        print("B =", B, "T")
        print("res =", res, "  run_time =", run_time)
        print("rho elements =", len(rho))
        print("ports =", [p + 1 for p in port_list])

    _pmm_device, P_device, wp_values = build_circulator_device(rho, B, res=res)

    power_matrix = np.full((n_ports, n_ports), np.nan)

    for source_port in port_list:
        if verbose:
            print()
            print("-" * 54)
            print(f"Input P{source_port + 1}")
            print("-" * 54)

        P_device.sources = make_port_source(source_port)
        sim_i = P_device.Get_Sim()

        monitors_i = []
        signs_i = []
        for output_port in port_list:
            region, sign = make_flux_region(
                port_monitor_centers[output_port],
                port_dirs[output_port],
            )
            monitor = sim_i.add_flux(fs_a, 0, 1, region)
            monitors_i.append(monitor)
            signs_i.append(sign)

        # Incident-field subtraction at SOURCE port (reflection on diagonal)
        src_local = local_index[source_port]
        sim_i.load_minus_flux_data(
            monitors_i[src_local],
            incident_flux_data_by_port[source_port],
        )

        if verbose:
            print("Running Meep...")
        sim_i.run(until_after_sources=run_time)

        flux_i = np.zeros(n_ports)
        for k, output_port in enumerate(port_list):
            raw_flux = mp.get_fluxes(monitors_i[k])[0]
            flux_i[k] = signs_i[k] * raw_flux

        incident_i = incident_power_by_port[source_port]
        column_i = flux_i / incident_i
        power_matrix[:, src_local] = column_i

        if verbose:
            print(f"P{source_port + 1} column:")
            for k, output_port in enumerate(port_list):
                print(f"  -> P{output_port + 1}: {column_i[k]:.6e}")

    power_matrix_dB = np.full((n_ports, n_ports), np.nan)
    positive = power_matrix > 0
    power_matrix_dB[positive] = 10 * np.log10(power_matrix[positive])

    negative_entries = np.argwhere(power_matrix < 0)
    if len(negative_entries) > 0:
        print()
        print("WARNING:")
        print(len(negative_entries), "negative normalized power entries were found.")
        print(
            "Do not send this result to an optimizer until those entries are understood."
        )
        objective_result = None
    elif n_ports == 6 and port_list == list(range(6)):
        objective_result = circulator_objective(power_matrix, direction=direction)
    else:
        objective_result = None

    if verbose:
        print()
        print("=" * 54)
        print("SIMULATION COMPLETE")
        print("=" * 54)
        print("Power matrix [dB]:")
        print(np.round(power_matrix_dB, 2))
        if objective_result is not None:
            print()
            print(direction, "circulator objective =", objective_result["objective"])
            print(
                "Mean adjacent isolation =",
                objective_result["mean_isolation_dB"],
                "dB",
            )

    return {
        "rho": rho.copy(),
        "B": B.copy(),
        "wp": np.asarray(wp_values).copy(),
        "ports": list(port_list),
        "power_matrix": power_matrix,
        "power_matrix_dB": power_matrix_dB,
        "objective": objective_result,
        "negative_entries": negative_entries,
        "res": int(res),
        "run_time": float(run_time),
    }


def geometry_summary() -> str:
    lines = [
        f"a = {a} m",
        f"nx_ports = {nx_ports}, ny_ports = {ny_ports}, dpml_ports = {dpml_ports}",
        f"feed_length = {feed_length_m * 1e3:.1f} mm",
        f"clear_width = {clear_width * a * 1e3:.3f} mm",
        f"fs = {fs_Hz / 1e9:.2f} GHz -> {fs_a:.6f} c/a",
        f"fp = {fp_Hz / 1e9:.2f} GHz -> {fp_a:.6f} c/a",
        f"n_elements (assert 91) = 91",
    ]
    for i in range(6):
        u = port_dirs[i]
        lines.append(
            f"  P{i + 1}: dir=[{u[0]:.4f}, {u[1]:.4f}]  "
            f"monitor={np.round(port_monitor_centers[i], 3)}"
        )
    return "\n".join(lines)


if __name__ == "__main__":
    print(geometry_summary())
    rho0 = default_uniform_rho()
    assert len(rho0) == 91, f"expected 91 elements, got {len(rho0)}"
    print("default_uniform_rho length assert: 91 OK")
