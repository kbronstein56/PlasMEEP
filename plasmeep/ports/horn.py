"""
PEC horn flare + straight feed geometry (2D, extruded along z).

Dimensions match the circulator validation harness (PMMInverse TM horn sizes).
Polarization is selected at source/monitor construction time (Ez or Hz), not here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Literal, Sequence, Tuple

import numpy as np

HornWallRepresentation = Literal[
    "prism",  # mp.Prism extrusions (baseline, staircased)
    "rotated_blocks",  # oriented mp.Block walls (less oblique bias)
    "grid_snapped_prism",  # prism vertices snapped to Yee nodes at res
]

DEFAULT_FLARE_STEPS = 4


@dataclass(frozen=True)
class HornWallSet:
  left_flare: np.ndarray
  right_flare: np.ndarray
  left_feed: np.ndarray
  right_feed: np.ndarray


@dataclass(frozen=True)
class HornGeometry:
  """One oriented horn in Meep xy coordinates (a-units)."""

  outward_dir: np.ndarray
  throat_center: np.ndarray
  monitor_center: np.ndarray
  source_center: np.ndarray
  walls: HornWallSet

  @property
  def tangent(self) -> np.ndarray:
    u = self.outward_dir / np.linalg.norm(self.outward_dir)
    return np.array([-u[1], u[0]])


def build_horn_geometry(
  open_center_xy: np.ndarray,
  outward_dir: np.ndarray,
  *,
  a_m: float,
  wall_thickness_m: float = 0.004,
  width_open_m: float = 0.104,
  width_base_m: float = 0.048,
  horn_depth_m: float = 0.089,
  feed_length_m: float = 0.060,
  monitor_frac: float = 0.30,
  source_frac: float = 0.70,
) -> HornGeometry:
  """
  Build horn polygons and probe centers.

  All lengths are supplied in metres and converted to Meep a-units via `a_m`.
  """
  a_m = float(a_m)
  wall_thickness = wall_thickness_m / a_m
  width_open = width_open_m / a_m
  width_base = width_base_m / a_m
  horn_depth = horn_depth_m / a_m
  feed_length = feed_length_m / a_m

  xy_open_cen = np.asarray(open_center_xy, dtype=float)
  outward_dir = np.asarray(outward_dir, dtype=float)
  outward_dir = outward_dir / np.linalg.norm(outward_dir)
  horn_dir = -outward_dir
  horn_dir_orth = np.array([horn_dir[1], -horn_dir[0]])

  def _quad(
    p0: np.ndarray,
    p1: np.ndarray,
    p2: np.ndarray,
    p3: np.ndarray,
  ) -> np.ndarray:
    return np.array([p0, p1, p2, p3])

  left_open = _quad(
    xy_open_cen + horn_dir_orth * width_open / 2,
    xy_open_cen + horn_dir_orth * (width_open / 2 - wall_thickness),
    xy_open_cen
    + horn_dir_orth * (width_base / 2 - wall_thickness)
    - horn_dir * horn_depth,
    xy_open_cen + horn_dir_orth * width_base / 2 - horn_dir * horn_depth,
  )
  left_feed = _quad(
    xy_open_cen
    + horn_dir_orth * (width_base / 2 - wall_thickness)
    - horn_dir * horn_depth,
    xy_open_cen + horn_dir_orth * width_base / 2 - horn_dir * horn_depth,
    xy_open_cen
    + horn_dir_orth * width_base / 2
    - horn_dir * (horn_depth + feed_length),
    xy_open_cen
    + horn_dir_orth * (width_base / 2 - wall_thickness)
    - horn_dir * (horn_depth + feed_length),
  )
  right_open = _quad(
    xy_open_cen - horn_dir_orth * width_open / 2,
    xy_open_cen - horn_dir_orth * (width_open / 2 - wall_thickness),
    xy_open_cen
    - horn_dir_orth * (width_base / 2 - wall_thickness)
    - horn_dir * horn_depth,
    xy_open_cen - horn_dir_orth * width_base / 2 - horn_dir * horn_depth,
  )
  right_feed = _quad(
    xy_open_cen
    - horn_dir_orth * (width_base / 2 - wall_thickness)
    - horn_dir * horn_depth,
    xy_open_cen - horn_dir_orth * width_base / 2 - horn_dir * horn_depth,
    xy_open_cen
    - horn_dir_orth * width_base / 2
    - horn_dir * (horn_depth + feed_length),
    xy_open_cen
    - horn_dir_orth * (width_base / 2 - wall_thickness)
    - horn_dir * (horn_depth + feed_length),
  )

  throat_center = xy_open_cen + outward_dir * horn_depth
  monitor_center = throat_center + outward_dir * (monitor_frac * feed_length)
  source_center = throat_center + outward_dir * (source_frac * feed_length)

  return HornGeometry(
    outward_dir=outward_dir,
    throat_center=throat_center,
    monitor_center=monitor_center,
    source_center=source_center,
    walls=HornWallSet(
      left_flare=left_open,
      right_flare=right_open,
      left_feed=left_feed,
      right_feed=right_feed,
    ),
  )


def horn_to_legacy_dict(horn: HornGeometry) -> Dict[str, np.ndarray]:
  """Match keys used by sixport_common.full_horns entries."""
  w = horn.walls
  return {
    "left_flare": w.left_flare,
    "right_flare": w.right_flare,
    "left_feed": w.left_feed,
    "right_feed": w.right_feed,
    "throat_center": horn.throat_center,
    "monitor_center": horn.monitor_center,
    "source_center": horn.source_center,
  }


def offset_cells_to_a(offset_cells: Tuple[float, float], res: int) -> np.ndarray:
  """Convert fractional Yee-cell offsets to Meep a-units."""
  ox, oy = offset_cells
  r = float(res)
  return np.array([ox / r, oy / r], dtype=float)


def translate_horn_geometry(horn: HornGeometry, delta_xy: np.ndarray) -> HornGeometry:
  """Shift horn polygons and probe centers by delta_xy (a-units)."""
  d = np.asarray(delta_xy, dtype=float).reshape(2)

  def _shift(poly: np.ndarray) -> np.ndarray:
    out = np.asarray(poly, dtype=float).copy()
    out[:, 0:2] += d
    return out

  w = horn.walls
  return HornGeometry(
    outward_dir=np.asarray(horn.outward_dir, dtype=float),
    throat_center=horn.throat_center + d,
    monitor_center=horn.monitor_center + d,
    source_center=horn.source_center + d,
    walls=HornWallSet(
      left_flare=_shift(w.left_flare),
      right_flare=_shift(w.right_flare),
      left_feed=_shift(w.left_feed),
      right_feed=_shift(w.right_feed),
    ),
  )


def translate_legacy_horn(
  horn_dict: Dict[str, np.ndarray], delta_xy: np.ndarray
) -> Dict[str, np.ndarray]:
  """Shift legacy sixport horn dict keys by delta_xy (a-units)."""
  d = np.asarray(delta_xy, dtype=float).reshape(2)
  out: Dict[str, np.ndarray] = {}
  for key, val in horn_dict.items():
    arr = np.asarray(val, dtype=float).copy()
    if arr.ndim == 2 and arr.shape[1] >= 2:
      arr[:, 0:2] += d
    elif arr.ndim == 1 and arr.size >= 2:
      arr = arr.copy()
      arr[0:2] += d
    out[key] = arr
  return out


@dataclass(frozen=True)
class OrientedWallBlock:
  """Thin PEC wall as an oriented Meep Block (2D extruded along z)."""

  center: np.ndarray  # (3,)
  length: float
  thickness: float
  e1: np.ndarray  # unit vector along wall length
  e2: np.ndarray  # unit vector along wall thickness


def _rotation_matrix(angle_rad: float) -> np.ndarray:
  c, s = np.cos(angle_rad), np.sin(angle_rad)
  return np.array([[c, -s], [s, c]], dtype=float)


def rotate_xy(points: np.ndarray, angle_rad: float) -> np.ndarray:
  """Rotate 2D points (..., 2) about the origin."""
  pts = np.asarray(points, dtype=float)
  r = _rotation_matrix(angle_rad)
  flat = pts.reshape(-1, 2)
  out = flat @ r.T
  return out.reshape(pts.shape)


def rotate_dir2(dir_xy: np.ndarray, angle_rad: float) -> np.ndarray:
  v = np.asarray(dir_xy, dtype=float).reshape(2)
  return _rotation_matrix(angle_rad) @ v


def quad_to_oriented_block(vertices: np.ndarray) -> OrientedWallBlock:
  """
  Convert a wall quadrilateral (4,2) to center/length/thickness/e1/e2.

  Vertex order matches legacy horn walls:
    outer-open, inner-open, inner-throat, outer-throat.
  """
  v = np.asarray(vertices, dtype=float)[:, :2]
  mid_open = 0.5 * (v[0] + v[1])
  mid_throat = 0.5 * (v[2] + v[3])
  center_xy = 0.5 * (mid_open + mid_throat)
  e1 = mid_throat - mid_open
  length = float(np.linalg.norm(e1))
  if length <= 0:
    raise ValueError("degenerate wall length")
  e1 = e1 / length
  e2 = v[1] - v[0]
  thickness = float(np.linalg.norm(e2))
  if thickness <= 0:
    raise ValueError("degenerate wall thickness")
  e2 = e2 / thickness
  return OrientedWallBlock(
    center=np.array([center_xy[0], center_xy[1], 0.0]),
    length=length,
    thickness=thickness,
    e1=np.array([e1[0], e1[1], 0.0]),
    e2=np.array([e2[0], e2[1], 0.0]),
  )


def subdivide_wall_quad(vertices: np.ndarray, n_steps: int) -> List[np.ndarray]:
  """Split a wall quad into n_steps sub-quads along the horn depth direction."""
  if n_steps <= 1:
    return [np.asarray(vertices, dtype=float)]
  v = np.asarray(vertices, dtype=float)
  quads: List[np.ndarray] = []
  for k in range(n_steps):
    t0 = k / n_steps
    t1 = (k + 1) / n_steps
    p0 = (1 - t0) * v[0] + t0 * v[3]
    p1 = (1 - t0) * v[1] + t0 * v[2]
    p2 = (1 - t1) * v[1] + t1 * v[2]
    p3 = (1 - t1) * v[0] + t1 * v[3]
    quads.append(np.array([p0, p1, p2, p3]))
  return quads


def horn_oriented_wall_blocks(
  horn: HornGeometry,
  *,
  n_flare_steps: int = DEFAULT_FLARE_STEPS,
) -> List[OrientedWallBlock]:
  """PEC walls as oriented blocks; flare trapezoids are subdivided."""
  blocks: List[OrientedWallBlock] = []
  w = horn.walls
  for quad in (w.left_feed, w.right_feed):
    blocks.append(quad_to_oriented_block(quad))
  for quad in (w.left_flare, w.right_flare):
    for sub in subdivide_wall_quad(quad, n_flare_steps):
      blocks.append(quad_to_oriented_block(sub))
  return blocks


def snap_legacy_horn_to_grid(horn_dict: Dict[str, np.ndarray], res: int) -> Dict[str, np.ndarray]:
  """Snap wall vertices and probe centers to the Yee grid at resolution res."""
  r = float(res)

  def _snap(arr: np.ndarray) -> np.ndarray:
    out = np.asarray(arr, dtype=float).copy()
    if out.ndim == 2 and out.shape[1] >= 2:
      out[:, 0] = np.round(out[:, 0] * r) / r
      out[:, 1] = np.round(out[:, 1] * r) / r
    elif out.ndim == 1 and out.size >= 2:
      out[0] = round(out[0] * r) / r
      out[1] = round(out[1] * r) / r
    return out

  return {key: _snap(val) for key, val in horn_dict.items()}


def legacy_horn_oriented_blocks(
  horn_dict: Dict[str, np.ndarray],
  *,
  n_flare_steps: int = DEFAULT_FLARE_STEPS,
) -> List[OrientedWallBlock]:
  """Oriented blocks from a legacy sixport horn dict."""
  blocks: List[OrientedWallBlock] = []
  for name in ("left_feed", "right_feed"):
    blocks.append(quad_to_oriented_block(horn_dict[name]))
  for name in ("left_flare", "right_flare"):
    for sub in subdivide_wall_quad(horn_dict[name], n_flare_steps):
      blocks.append(quad_to_oriented_block(sub))
  return blocks
  """Snap wall vertices and probe centers to the Yee grid at resolution res."""
  r = float(res)

  def _snap(arr: np.ndarray) -> np.ndarray:
    out = np.asarray(arr, dtype=float).copy()
    if out.ndim == 2 and out.shape[1] >= 2:
      out[:, 0] = np.round(out[:, 0] * r) / r
      out[:, 1] = np.round(out[:, 1] * r) / r
    elif out.ndim == 1 and out.size >= 2:
      out[0] = round(out[0] * r) / r
      out[1] = round(out[1] * r) / r
    return out

  return {key: _snap(val) for key, val in horn_dict.items()}


def rotate_legacy_horn_about_origin(
  horn_dict: Dict[str, np.ndarray], angle_rad: float
) -> Dict[str, np.ndarray]:
  """Rotate horn dict about the simulation origin."""
  out: Dict[str, np.ndarray] = {}
  for key, val in horn_dict.items():
    arr = np.asarray(val, dtype=float).copy()
    if arr.ndim == 2 and arr.shape[1] >= 2:
      arr[:, :2] = rotate_xy(arr[:, :2], angle_rad)
    elif arr.ndim == 1 and arr.size >= 2:
      arr[:2] = rotate_xy(arr[:2], angle_rad)
    out[key] = arr
  return out
