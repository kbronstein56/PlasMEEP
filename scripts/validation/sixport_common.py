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
    from plasmeep.ports.horn import build_horn_geometry, horn_to_legacy_dict

    horn = build_horn_geometry(
        xy_open_cen,
        outward_dir,
        a_m=a,
        feed_length_m=feed_length_local * a,
    )
    return horn_to_legacy_dict(horn)


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


# Active geometry context (grid registration relative to Yee mesh)
_GEOM_CTX: Dict[str, Any] = {
    "grid_offset_cells": (0.0, 0.0),
    "monitor_offset_cells": (0.0, 0.0),
    "res": None,
    "horn_walls": "prism",
    "coord_rotation_deg": 0.0,
    "flare_steps": 4,
}


def set_geometry_context(
    *,
    grid_offset_cells: Optional[Tuple[float, float]] = None,
    monitor_offset_cells: Optional[Tuple[float, float]] = None,
    res: Optional[int] = None,
    horn_walls: Optional[str] = None,
    coord_rotation_deg: Optional[float] = None,
    flare_steps: Optional[int] = None,
) -> None:
    """Set geometry registration / horn-wall representation context."""
    if grid_offset_cells is not None:
        _GEOM_CTX["grid_offset_cells"] = (
            float(grid_offset_cells[0]),
            float(grid_offset_cells[1]),
        )
    if monitor_offset_cells is not None:
        _GEOM_CTX["monitor_offset_cells"] = (
            float(monitor_offset_cells[0]),
            float(monitor_offset_cells[1]),
        )
    if res is not None:
        _GEOM_CTX["res"] = int(res)
    if horn_walls is not None:
        _GEOM_CTX["horn_walls"] = str(horn_walls)
    if coord_rotation_deg is not None:
        _GEOM_CTX["coord_rotation_deg"] = float(coord_rotation_deg)
    if flare_steps is not None:
        _GEOM_CTX["flare_steps"] = int(flare_steps)


def get_horn_walls() -> str:
    return str(_GEOM_CTX["horn_walls"])


def get_coord_rotation_deg() -> float:
    return float(_GEOM_CTX["coord_rotation_deg"])


def get_flare_steps() -> int:
    return int(_GEOM_CTX["flare_steps"])


def effective_port_dir(port_index: int) -> np.ndarray:
    """Outward port normal, including optional simulation-frame rotation."""
    from plasmeep.ports.horn import rotate_dir2

    u = np.asarray(port_dirs[port_index], dtype=float)
    angle = np.deg2rad(get_coord_rotation_deg())
    if abs(angle) < 1e-15:
        return u
    return rotate_dir2(u, angle)


def current_res() -> int:
    r = _GEOM_CTX.get("res")
    if r is None:
        raise RuntimeError("geometry res not set; call set_geometry_context(res=...)")
    return int(r)


def get_grid_offset_cells() -> Tuple[float, float]:
    return tuple(_GEOM_CTX["grid_offset_cells"])


def get_monitor_offset_cells() -> Tuple[float, float]:
    return tuple(_GEOM_CTX["monitor_offset_cells"])


def grid_offset_a(res: int) -> np.ndarray:
    from plasmeep.ports.horn import offset_cells_to_a

    return offset_cells_to_a(get_grid_offset_cells(), res)


def horn_for_port(port_index: int, res: int) -> Dict[str, np.ndarray]:
    from plasmeep.ports.horn import (
        rotate_legacy_horn_about_origin,
        snap_legacy_horn_to_grid,
        translate_legacy_horn,
    )

    horn = full_horns[port_index]
    angle = np.deg2rad(get_coord_rotation_deg())
    if abs(angle) > 1e-15:
        horn = rotate_legacy_horn_about_origin(horn, angle)
    delta = grid_offset_a(res)
    if not np.allclose(delta, 0):
        horn = translate_legacy_horn(horn, delta)
    if get_horn_walls() == "grid_snapped_prism":
        horn = snap_legacy_horn_to_grid(horn, res)
    return horn


def horns_for_device(res: int) -> List[Dict[str, np.ndarray]]:
    return [horn_for_port(i, res) for i in range(6)]


def monitor_center_for_port(port_index: int, res: int) -> np.ndarray:
    """Monitor plane center with optional monitor-only sub-cell offset."""
    from plasmeep.ports.horn import offset_cells_to_a

    center = np.asarray(horn_for_port(port_index, res)["monitor_center"], dtype=float)
    mo = get_monitor_offset_cells()
    if np.allclose(mo, 0):
        return center
    u = effective_port_dir(port_index)
    u = u / np.linalg.norm(u)
    tangent = np.array([-u[1], u[0]])
    delta_a = offset_cells_to_a(mo, res)
    # monitor_offset_cells[0] along outward, [1] along tangent
    return center + mo[0] * u / res + mo[1] * tangent / res


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

    center = np.asarray(horn_for_port(port_index, current_res())["source_center"], dtype=float)
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


def _mount_horn_walls(
    p_device,
    horn: Dict[str, np.ndarray],
    *,
    res: int,
    wall_pec: bool,
    wall_eps: float,
) -> None:
    """Add PEC/metal horn walls using the active horn-wall representation."""
    from plasmeep.ports.horn import legacy_horn_oriented_blocks

    rep = get_horn_walls()
    medium = (
        mp.perfect_electric_conductor
        if wall_pec
        else p_device.Get_Med(wall_eps, PEC=False)
    )

    if rep == "rotated_blocks":
        for block in legacy_horn_oriented_blocks(
            horn, n_flare_steps=get_flare_steps()
        ):
            p_device.geometry.append(
                mp.Block(
                    mp.Vector3(block.length, block.thickness, mp.inf),
                    center=mp.Vector3(
                        block.center[0], block.center[1], block.center[2]
                    ),
                    e1=mp.Vector3(block.e1[0], block.e1[1], block.e1[2]),
                    e2=mp.Vector3(block.e2[0], block.e2[1], block.e2[2]),
                    material=medium,
                )
            )
        return

    for name in ("left_flare", "right_flare", "left_feed", "right_feed"):
        wall = horn[name]
        vertices = np.zeros((4, 3))
        vertices[:, 0] = wall[:, 0]
        vertices[:, 1] = wall[:, 1]
        if wall_pec:
            p_device.Add_Prism(
                vertices=vertices,
                axis=np.array([0, 0, 1]),
                PEC=True,
            )
        else:
            p_device.Add_Prism(
                vertices=vertices,
                axis=np.array([0, 0, 1]),
                eps=wall_eps,
                PEC=False,
            )


# ---------------------------------------------------------------------------
# Device builder (cell 34) — resolution is an argument
# ---------------------------------------------------------------------------

def build_circulator_device(
    rho,
    B,
    res=64,
    wall_pec: bool = True,
    wall_eps: float = 1.0e6,
    device_mode: str = "full",
):
    """Build six-port PlasMEEP geometry; no sources/monitors/run.

    device_mode:
      full       — 91-element PMM array + horns (default)
      horns_only — PEC horns in vacuum (no plasma array / bulbs)
    """
    if device_mode not in ("full", "horns_only"):
        raise ValueError(f"device_mode must be 'full' or 'horns_only', got {device_mode!r}")

    rho = np.asarray(rho, dtype=float).flatten()
    B = np.asarray(B, dtype=float).flatten()
    if device_mode == "full" and len(rho) != 91:
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

    wp_values = []
    if device_mode == "full":
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

    if device_mode == "full":
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

    for horn in horns_for_device(res):
        _mount_horn_walls(
            P_device,
            horn,
            res=res,
            wall_pec=wall_pec,
            wall_eps=wall_eps,
        )

    return pmm_device, P_device, np.asarray(wp_values, dtype=float)


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

def normalize_port(
    source_port,
    res=64,
    run_time=80,
    verbose=True,
    formulation: str = "baseline_hz_line",
):
    """Straight-feed reference run for one port; returns (power, flux_data)."""
    from port_formulations import (
        EIGENMODE_WALL_EPS,
        add_dft_sdotn_monitor,
        add_flux_monitor,
        extract_dft_sdotn_power,
        extract_eigenmode_powers,
        get_formulation,
        make_flux_region_for_formulation,
        port_measure_center,
        uses_finite_metal_walls,
    )

    set_geometry_context(res=res)
    form = get_formulation(formulation)
    finite_metal = uses_finite_metal_walls(form.name)
    if verbose:
        print()
        print("=" * 54)
        print(f"NORMALIZING PORT P{source_port + 1}  [{form.name}]")
        print("=" * 54)

    u = np.asarray(effective_port_dir(source_port), dtype=float)
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
        horn_for_port(source_port, res)["source_center"], dtype=float
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
    if finite_metal:
        P_ref.Add_Prism(
            vertices=wall_A, axis=np.array([0, 0, 1]), eps=EIGENMODE_WALL_EPS
        )
        P_ref.Add_Prism(
            vertices=wall_B, axis=np.array([0, 0, 1]), eps=EIGENMODE_WALL_EPS
        )
    else:
        P_ref.Add_Prism(vertices=wall_A, axis=np.array([0, 0, 1]), PEC=True)
        P_ref.Add_Prism(vertices=wall_B, axis=np.array([0, 0, 1]), PEC=True)

    P_ref.sources = form.make_sources(source_port)
    ref_sim = P_ref.Get_Sim()

    measure_xy = port_measure_center(form.name, source_port)

    if form.measurement == "dft_sdotn":
        mon_info = add_dft_sdotn_monitor(ref_sim, measure_xy, effective_port_dir(source_port))
        if verbose:
            print("Running reference (DFT S·n)...")
        ref_sim.run(until_after_sources=run_time)
        signed_flux = extract_dft_sdotn_power(ref_sim, mon_info)
        # Reference run: power into the guide toward -outward (into device for
        # device ports). For a straight reference, expect negative outward flux.
        incident_power = abs(signed_flux)
        incident_data = None
    else:
        regions, sign = make_flux_region_for_formulation(
            form.name, measure_xy, effective_port_dir(source_port)
        )
        monitor = add_flux_monitor(ref_sim, regions)
        if verbose:
            print("Running reference...")
        ref_sim.run(until_after_sources=run_time)
        if form.measurement == "eigenmode":
            powers = extract_eigenmode_powers(
                ref_sim, [monitor], [source_port], [sign]
            )
            signed_flux = float(powers[0])
            incident_power = abs(signed_flux)
            incident_data = ref_sim.get_flux_data(monitor)
        else:
            raw_flux = mp.get_fluxes(monitor)[0]
            signed_flux = sign * raw_flux
            incident_power = abs(signed_flux)
            incident_data = ref_sim.get_flux_data(monitor)

    if verbose:
        print("signed incident flux/power =", signed_flux)
        print("incident power =", incident_power)
        if signed_flux < 0:
            print("Direction check: PASS")
        else:
            print("WARNING: expected inward incident flux/power.")

    return incident_power, incident_data


def ensure_normalizations(
    res,
    run_time=80,
    force=False,
    ports: Optional[Sequence[int]] = None,
    cache_path: Optional[str] = None,
    verbose=True,
    formulation: str = "baseline_hz_line",
) -> Dict[str, Any]:
    """
    Return incident powers and flux-subtraction data for the requested ports.

    Caches in memory by (res, run_time, formulation). Optional pickle at
    cache_path for cross-process reuse (--skip-norm-if-cached).
    """
    key = (
        int(res),
        float(run_time),
        str(formulation),
        get_grid_offset_cells(),
        get_monitor_offset_cells(),
        get_horn_walls(),
        get_coord_rotation_deg(),
    )
    port_list = list(range(6) if ports is None else ports)

    if cache_path and (not force) and os.path.isfile(cache_path):
        with open(cache_path, "rb") as f:
            loaded = pickle.load(f)
        if (
            loaded.get("res") == int(res)
            and loaded.get("run_time") == float(run_time)
            and loaded.get("formulation", "baseline_hz_line") == str(formulation)
            and loaded.get("grid_offset_cells", (0.0, 0.0)) == get_grid_offset_cells()
            and loaded.get("monitor_offset_cells", (0.0, 0.0))
            == get_monitor_offset_cells()
            and loaded.get("horn_walls", "prism") == get_horn_walls()
            and float(loaded.get("coord_rotation_deg", 0.0))
            == get_coord_rotation_deg()
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
            source_port,
            res=res,
            run_time=run_time,
            verbose=verbose,
            formulation=formulation,
        )
        incident_power_by_port[source_port] = power_i
        incident_flux_data_by_port[source_port] = data_i

    result = {
        "res": int(res),
        "run_time": float(run_time),
        "formulation": str(formulation),
        "grid_offset_cells": get_grid_offset_cells(),
        "monitor_offset_cells": get_monitor_offset_cells(),
        "horn_walls": get_horn_walls(),
        "coord_rotation_deg": get_coord_rotation_deg(),
        "incident_power_by_port": incident_power_by_port,
        "incident_flux_data_by_port": incident_flux_data_by_port,
    }
    _NORM_CACHE[key] = result

    if cache_path:
        os.makedirs(os.path.dirname(os.path.abspath(cache_path)) or ".", exist_ok=True)
        # Rank-0 only: avoid MPI pickle races.
        try:
            from mpi4py import MPI

            rank = int(MPI.COMM_WORLD.Get_rank())
        except Exception:
            rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", "0"))
        if rank == 0:
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
    formulation: str = "baseline_hz_line",
    device_mode: str = "full",
    wall_pec: bool = True,
):
    """
    Full six-port (or port-subset) forward model.

    Matrix convention: row = OUTPUT, column = INPUT.
    Port excitation/measurement selected by `formulation`.
    """
    from port_formulations import (
        EIGENMODE_WALL_EPS,
        add_dft_sdotn_monitor,
        add_flux_monitor,
        extract_dft_sdotn_power,
        extract_eigenmode_powers,
        extract_flux_powers,
        get_formulation,
        make_flux_region_for_formulation,
        port_measure_center,
        uses_finite_metal_walls,
    )

    form = get_formulation(formulation)
    set_geometry_context(res=res)
    rho = np.asarray(rho, dtype=float).flatten()
    B = np.asarray(B, dtype=float).flatten()
    if device_mode == "full" and len(rho) != 91:
        raise ValueError(f"Expected 91 rho values, got {len(rho)}.")
    if len(B) != 3:
        raise ValueError("B must be [Bx, By, Bz].")

    port_list = list(range(6) if ports is None else ports)
    n_ports = len(port_list)
    # Map global port index -> local matrix index
    local_index = {p: k for k, p in enumerate(port_list)}

    if incident_cache is None:
        incident_cache = ensure_normalizations(
            res,
            run_time=run_time,
            force=False,
            ports=port_list,
            verbose=verbose,
            formulation=formulation,
        )
    if incident_cache.get("formulation", "baseline_hz_line") != form.name:
        raise ValueError(
            f"incident_cache formulation "
            f"{incident_cache.get('formulation')!r} != {form.name!r}"
        )
    incident_power_by_port = incident_cache["incident_power_by_port"]
    incident_flux_data_by_port = incident_cache["incident_flux_data_by_port"]
    for p in port_list:
        if p not in incident_power_by_port:
            raise KeyError(
                f"incident_cache missing port {p}; call ensure_normalizations first."
            )
        # DFT formulations do not provide flux-subtraction data.
        if form.measurement != "dft_sdotn" and p not in incident_flux_data_by_port:
            raise KeyError(
                f"incident_cache missing flux data for port {p}."
            )

    if verbose:
        print("=" * 54)
        print("SIX-PORT CIRCULATOR FORWARD MODEL")
        print("=" * 54)
        print("B =", B, "T")
        print("res =", res, "  run_time =", run_time)
        print("formulation =", form.name)
        print("device_mode =", device_mode)
        print("rho elements =", len(rho))
        print("ports =", [p + 1 for p in port_list])

    finite_metal = uses_finite_metal_walls(form.name)
    _pmm_device, P_device, wp_values = build_circulator_device(
        rho,
        B,
        res=res,
        wall_pec=wall_pec and (not finite_metal),
        wall_eps=EIGENMODE_WALL_EPS,
        device_mode=device_mode,
    )

    power_matrix = np.full((n_ports, n_ports), np.nan)

    for source_port in port_list:
        if verbose:
            print()
            print("-" * 54)
            print(f"Input P{source_port + 1}")
            print("-" * 54)

        P_device.sources = form.make_sources(source_port)
        sim_i = P_device.Get_Sim()

        monitors_i = []
        signs_i = []
        dft_infos = []
        for output_port in port_list:
            measure_xy = port_measure_center(form.name, output_port)
            if form.measurement == "dft_sdotn":
                dft_infos.append(
                    add_dft_sdotn_monitor(
                        sim_i, measure_xy, effective_port_dir(output_port)
                    )
                )
                monitors_i.append(None)
                signs_i.append(1.0)
            else:
                regions, sign = make_flux_region_for_formulation(
                    form.name,
                    measure_xy,
                    effective_port_dir(output_port),
                )
                monitor = add_flux_monitor(sim_i, regions)
                monitors_i.append(monitor)
                signs_i.append(sign)
                dft_infos.append(None)

        # Incident-field subtraction at SOURCE port (reflection on diagonal)
        src_local = local_index[source_port]
        if form.measurement != "dft_sdotn":
            sim_i.load_minus_flux_data(
                monitors_i[src_local],
                incident_flux_data_by_port[source_port],
            )

        if verbose:
            print("Running Meep...")
        sim_i.run(until_after_sources=run_time)

        if form.measurement == "dft_sdotn":
            flux_i = np.array(
                [extract_dft_sdotn_power(sim_i, info) for info in dft_infos]
            )
        elif form.measurement == "eigenmode":
            flux_i = extract_eigenmode_powers(
                sim_i, monitors_i, port_list, signs_i
            )
        else:
            flux_i = extract_flux_powers(sim_i, monitors_i, signs_i)

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
        "formulation": form.name,
        "device_mode": device_mode,
    }


def geometry_dimensions_table() -> Dict[str, Any]:
    """Important lengths in both a-units and SI for audit reports."""
    def mm(x_a: float) -> float:
        return float(x_a * a * 1000)

    return {
        "a_m": float(a),
        "a_cm": float(a * 100),
        "lattice_pitch_a": float(d_exp),
        "lattice_pitch_mm": 20.0,
        "quartz_od_mm": mm(r_bulb_outer * 2),
        "quartz_id_mm": mm(r_bulb_inner * 2),
        "wall_thickness_mm": mm(wall_thickness),
        "horn_aperture_mm": mm(width_open),
        "horn_throat_mm": mm(width_base),
        "horn_depth_mm": mm(horn_depth),
        "feed_length_mm": mm(feed_length),
        "clear_width_mm": mm(clear_width),
        "plasma_radius_mm": mm(r_plasma),
        "domain_nx_ny": [int(nx_ports), int(ny_ports)],
        "domain_extent_mm": [mm(nx_ports), mm(ny_ports)],
    }


def physical_resolution_report(res: int) -> Dict[str, float]:
    """Map Meep resolution to physical grid metrics (a is the Meep length unit)."""
    a_cm = float(a) * 100.0
    a_mm = float(a) * 1000.0
    dx_a = 1.0 / float(res)
    dx_mm = dx_a * a_mm
    pixels_per_cm = float(res) / a_cm
    lattice_mm = 20.0
    pixels_per_lattice = lattice_mm / dx_mm
    c_mps = 2.99792458e8
    lambda0_mm = c_mps / fs_Hz * 1000.0
    pixels_per_lambda0 = lambda0_mm / dx_mm
    return {
        "a_m": float(a),
        "a_cm": a_cm,
        "meep_resolution": float(res),
        "grid_spacing_a": dx_a,
        "grid_spacing_mm": dx_mm,
        "pixels_per_cm": pixels_per_cm,
        "pixels_per_lattice_20mm": pixels_per_lattice,
        "pixels_per_free_space_lambda0": pixels_per_lambda0,
        "lambda0_mm": lambda0_mm,
    }


def geometry_summary() -> str:
    lines = [
        f"a = {a} m = {a * 100:.2f} cm  (Meep length unit; NOT the 20 mm lattice pitch)",
        f"lattice pitch = 20 mm = {0.020 / a:.6f} a-units",
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
