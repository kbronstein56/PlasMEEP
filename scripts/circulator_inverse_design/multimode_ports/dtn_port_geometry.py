"""Locate outer straight-feed port segments on the FULL91 mesh."""
from __future__ import annotations

import numpy as np

from feed_constants import clear_width


def port_outer_center(horn: dict, outward: np.ndarray) -> np.ndarray:
    """Center of outer feed mouth (beyond legacy source at 0.7·feed)."""
    throat = np.asarray(horn["throat_center"], float)
    u = np.asarray(outward, float)
    u = u / np.linalg.norm(u)
    # horn legacy uses throat + outward * frac * feed; outer end at full feed from throat
    # throat_center = open + outward*horn_depth; outer = throat + outward*feed_length
    src = np.asarray(horn["source_center"], float)
    mon = np.asarray(horn["monitor_center"], float)
    # source at 0.7 feed from throat, monitor 0.3 → feed_length = (src - throat)/u / 0.7
    feed_len = float(np.dot(src - throat, u) / 0.7)
    return throat + u * feed_len


def mesh_port_nodes(
    points: np.ndarray,
    center: np.ndarray,
    outward: np.ndarray,
    tangent: np.ndarray,
    *,
    W: float = clear_width,
    plane_tol: float = 0.06,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Nodes on the clear aperture at the outer port plane (mesh coordinates)."""
    u = outward / np.linalg.norm(outward)
    t = tangent / np.linalg.norm(tangent)
    rel = points - center[None, :]
    along = rel @ u
    cross = rel @ t
    mask = (np.abs(along) <= plane_tol) & (np.abs(cross) <= 0.5 * W + 1e-9)
    idx = np.flatnonzero(mask)
    order = idx[np.argsort(cross[idx])]
    y = cross[order]
    if len(y) == 0:
        return order, y, np.array([])
    if len(y) == 1:
        w = np.array([W])
    else:
        ds = np.diff(y)
        w = np.zeros(len(y))
        w[0] = 0.5 * ds[0]
        w[-1] = 0.5 * ds[-1]
        w[1:-1] = 0.5 * (ds[:-1] + ds[1:])
    return order, y, w
