#!/usr/bin/env python3
"""True 1-D P1 boundary edges for DtN ports after exterior excision."""
from __future__ import annotations

from typing import Sequence

import numpy as np
from scipy import sparse

from feed_constants import clear_width
from modes import normalize_on_line, transverse_mode


def physical_triangle_mask(rho_xx, rho_xy, rho_yx, rho_yy) -> np.ndarray:
    """Triangles that remain in the FEM (not ρ→0 excised/PEC)."""
    return ~(
        (np.abs(rho_xx) < 1e-10)
        & (np.abs(rho_xy) < 1e-10)
        & (np.abs(rho_yx) < 1e-10)
        & (np.abs(rho_yy) < 1e-10)
    )


def boundary_edges_from_physical(tris: np.ndarray, physical: np.ndarray) -> list[tuple[int, int]]:
    """Undirected edges belonging to exactly one physical triangle."""
    edge_count: dict[tuple[int, int], int] = {}
    for t, ok in zip(tris, physical):
        if not ok:
            continue
        nodes = [int(t[0]), int(t[1]), int(t[2])]
        for a, b in ((0, 1), (1, 2), (2, 0)):
            e = (nodes[a], nodes[b])
            key = (min(e), max(e))
            edge_count[key] = edge_count.get(key, 0) + 1
    return [e for e, c in edge_count.items() if c == 1]


def select_port_boundary_edges(
    points: np.ndarray,
    edges: Sequence[tuple[int, int]],
    center: np.ndarray,
    outward: np.ndarray,
    tangent: np.ndarray,
    *,
    W: float = clear_width,
    plane_tol: float = 0.04,
) -> list[tuple[int, int]]:
    """Boundary edges whose midpoint lies on the clear-aperture port plane."""
    u = outward / np.linalg.norm(outward)
    t = tangent / np.linalg.norm(tangent)
    c = np.asarray(center, float)
    out = []
    for i, j in edges:
        mid = 0.5 * (points[i] + points[j])
        along = float((mid - c) @ u)
        cross = float((mid - c) @ t)
        if abs(along) <= plane_tol and abs(cross) <= 0.5 * W + 1e-8:
            evec = points[j] - points[i]
            elen = float(np.linalg.norm(evec))
            if elen < 1e-14:
                continue
            # Keep staircase edges on the cut; only drop edges almost parallel to
            # the guide axis (feed sidewalls), which are not aperture-crossing.
            if abs(float(evec @ u)) / elen > 0.92:
                continue
            out.append((i, j))
    return out


def order_port_polyline(
    points: np.ndarray,
    edges: Sequence[tuple[int, int]],
    center: np.ndarray,
    tangent: np.ndarray,
) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """Order edges into a polyline sorted by tangent coordinate; orient consistently."""
    if not edges:
        return np.array([], dtype=int), []
    t = tangent / np.linalg.norm(tangent)
    # Build adjacency
    adj: dict[int, list[int]] = {}
    for i, j in edges:
        adj.setdefault(i, []).append(j)
        adj.setdefault(j, []).append(i)
    # Endpoints: degree 1; if loop, start at min cross
    ends = [n for n, nb in adj.items() if len(nb) == 1]
    if len(ends) >= 1:
        # Choose the end with smaller cross as start
        start = min(ends, key=lambda n: float((points[n] - center) @ t))
    else:
        start = min(adj.keys(), key=lambda n: float((points[n] - center) @ t))
    ordered_nodes = [start]
    ordered_edges: list[tuple[int, int]] = []
    prev = None
    cur = start
    seen_e = set()
    while True:
        nbs = [n for n in adj[cur] if n != prev]
        if not nbs:
            break
        # Prefer unused edge
        nxt = None
        for n in nbs:
            key = (min(cur, n), max(cur, n))
            if key not in seen_e:
                nxt = n
                break
        if nxt is None:
            break
        key = (min(cur, nxt), max(cur, nxt))
        seen_e.add(key)
        ordered_edges.append((cur, nxt))
        ordered_nodes.append(nxt)
        prev, cur = cur, nxt
        if len(ordered_edges) > len(edges) + 2:
            break
    return np.asarray(ordered_nodes, dtype=int), ordered_edges


def p1_edge_modal_load(
    points: np.ndarray,
    edges: Sequence[tuple[int, int]],
    phi_at_nodes: dict[int, float],
) -> tuple[np.ndarray, np.ndarray, float]:
    """Assemble c_i = Σ_edges ∫_e φ N_i ds for P1 hats on boundary edges.

    Returns (node_indices, c_values, total_length).
    """
    contrib: dict[int, complex] = {}
    total_L = 0.0
    for i, j in edges:
        pi, pj = points[i], points[j]
        L = float(np.linalg.norm(pj - pi))
        if L < 1e-18:
            continue
        total_L += L
        phi_i = float(phi_at_nodes[i])
        phi_j = float(phi_at_nodes[j])
        # ∫_e φ N_i ds with φ linear, N_i linear: (L/6)*(2 φ_i + φ_j)
        contrib[i] = contrib.get(i, 0.0) + (L / 6.0) * (2.0 * phi_i + phi_j)
        contrib[j] = contrib.get(j, 0.0) + (L / 6.0) * (2.0 * phi_j + phi_i)
    if not contrib:
        return np.array([], dtype=int), np.array([], dtype=np.complex128), 0.0
    idx = np.array(sorted(contrib.keys()), dtype=int)
    c = np.array([contrib[i] for i in idx], dtype=np.complex128)
    return idx, c, total_L


def build_port_from_boundary(
    points: np.ndarray,
    tris: np.ndarray,
    rho,
    center: np.ndarray,
    outward: np.ndarray,
    tangent: np.ndarray,
    *,
    W: float = clear_width,
    plane_tol: float = 0.04,
    port: int = 0,
) -> dict:
    """Construct a DtN port with true boundary edges and P1 edge loads."""
    physical = physical_triangle_mask(*rho)
    all_bedges = boundary_edges_from_physical(tris, physical)
    pedges = select_port_boundary_edges(
        points, all_bedges, center, outward, tangent, W=W, plane_tol=plane_tol
    )
    # Use ALL selected edges for P1 loads (staircase cut may be multi-component).
    # Polyline ordering is only diagnostic.
    nodes_ord, _ordered = order_port_polyline(points, pedges, center, tangent)
    node_set = sorted({n for e in pedges for n in e})
    nodes = np.asarray(node_set, dtype=int)
    t = tangent / np.linalg.norm(tangent)
    c0 = np.asarray(center, float)
    y = np.array([float((points[i] - c0) @ t) for i in nodes]) if len(nodes) else np.array([])
    total_L = float(sum(np.linalg.norm(points[j] - points[i]) for i, j in pedges))
    # Projected aperture span along tangent (should be ≈ W)
    if len(nodes):
        crosses = np.array([float((points[i] - c0) @ t) for i in nodes])
        y_span = float(crosses.max() - crosses.min())
    else:
        y_span = 0.0
    return {
        "port": port,
        "outward": np.asarray(outward, float) / np.linalg.norm(outward),
        "tangent": t,
        "center": c0,
        "nodes": nodes,
        "y": y,
        "edges": list(pedges),
        "n_nodes": int(len(nodes)),
        "n_edges": int(len(pedges)),
        "total_length": total_L,
        "y_span": y_span,
        "polyline_nodes": nodes_ord,
    }


def modal_c_vector(
    points: np.ndarray,
    port: dict,
    m: int,
    W: float = clear_width,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (global_node_indices, c) for mode m on port boundary edges."""
    if port["n_edges"] == 0:
        return np.array([], dtype=int), np.array([], dtype=np.complex128)
    # Unnormalized φ then L2-normalize using edge quadrature
    edges = port["edges"]
    # First pass: compute ∫ φ_raw² for normalization
    raw_at = {}
    for i, j in edges:
        for n in (i, j):
            if n not in raw_at:
                y = float((points[n] - port["center"]) @ port["tangent"])
                raw_at[n] = float(transverse_mode(np.array([y]), m, W)[0])
    # ∫ φ² ds with linear φ: on edge (L/6)*(φi²+φiφj+φj²)*2? 
    # Exact for linear: ∫_0^L (φi(1-s)+φj s)² ds = L/3 (φi²+φiφj+φj²)
    nrm2 = 0.0
    for i, j in edges:
        L = float(np.linalg.norm(points[j] - points[i]))
        a, b = raw_at[i], raw_at[j]
        nrm2 += (L / 3.0) * (a * a + a * b + b * b)
    nrm = float(np.sqrt(max(nrm2, 1e-30)))
    phi_at = {n: raw_at[n] / nrm for n in raw_at}
    idx, c, _L = p1_edge_modal_load(points, edges, phi_at)
    return idx, c


def add_modal_robin_from_c(
    A: sparse.spmatrix,
    idx: np.ndarray,
    c: np.ndarray,
    beta: complex,
) -> sparse.spmatrix:
    """Add frozen (−i β) c c^T Robin block (propagating) / (−α) for α>0 evanescent.

    Propagating: beta_eff = −i β  (matches validated straight-guide sign).
    Evanescent β=iα: first-order outgoing ∂n H = −α H → add (−α) after replacing
    natural ∫(∂n)v; with α>0 use beta_eff = −α = −Im(β) when β=iα... 
    Existing code used +iα for evanescent which equals −(−iα)= −β_eff_prop pattern.
    Keep identical frozen map as add_modal_robin.
    """
    beta = complex(beta)
    if abs(beta.imag) > 1e-12 and beta.imag > 0:
        beta_eff = 1j * beta.imag  # same as dtn_assembly.add_modal_robin
    else:
        beta_eff = -1j * float(beta.real)
    block = beta_eff * np.outer(c, c)
    A = A.tolil()
    for i, p in enumerate(idx):
        for j, q in enumerate(idx):
            A[p, q] = A[p, q] + block[i, j]
    return A.tocsr()


def add_incident_rhs_from_c(
    b: np.ndarray,
    idx: np.ndarray,
    c: np.ndarray,
    beta: complex,
    a_field: complex,
) -> np.ndarray:
    """RHS += 2 i β a_field c  (frozen, matches straight-guide)."""
    coeff = 2j * complex(beta) * complex(a_field)
    for i, p in enumerate(idx):
        b[int(p)] += coeff * c[i]
    return b
