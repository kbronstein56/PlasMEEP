#!/usr/bin/env python3
"""Single plasma inclusion: FEM fixed epsilon versus Meep Drude.

The disk (or an equal-area square) sits in air. A point Hz source is far enough
upstream that the wave across the inclusion is nearly planar. Vacuum and
inclusion runs share the source, so transmission and scattering ratios do not
depend on the source amplitude.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from fem_validated_solver import (  # noqa: E402
    LAST_OPERATOR,
    ElementSampler,
    assemble_anisotropic,
    eps_tensor_at_bias,
    meep_sigma_max,
    rho_from_eps,
    solve_system,
)

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "plasma_interface"
LX, LY, DPML = 18.0, 14.0, 2.0
# Source and inclusion in origin-centered coordinates. Mesh adds (LX/2, LY/2).
SRC_XY = np.array([-6.0, 0.0])
R_DISK = 4.6 * float(sc.r_bulb_inner) / 6.5


def origin_to_mesh(xy: np.ndarray) -> np.ndarray:
    return np.asarray(xy, dtype=float) + np.array([LX / 2.0, LY / 2.0])


def mesh_to_origin(xy: np.ndarray) -> np.ndarray:
    return np.asarray(xy, dtype=float) - np.array([LX / 2.0, LY / 2.0])


def plasma_eps() -> complex:
    return complex(eps_tensor_at_bias(0.0)[0])


def square_half() -> float:
    """Half-side of the equal-area square."""
    return 0.5 * R_DISK * np.sqrt(np.pi)


def _tri():
    vend = ROOT / "scripts" / "validation" / "fem_meep_validation" / "_vendor"
    if str(vend) not in sys.path:
        sys.path.insert(0, str(vend))
    import triangle as tr

    return tr


def build_mesh(h: float, h_edge: float, kind: str):
    """Constrained mesh of the PML box with one inclusion boundary."""
    tr = _tri()
    verts: list[list[float]] = []
    segs: list[list[int]] = []
    index: dict[tuple[float, float], int] = {}

    def add_vertex(p) -> int:
        key = (round(float(p[0]), 8), round(float(p[1]), 8))
        if key in index:
            return index[key]
        index[key] = len(verts)
        verts.append([key[0], key[1]])
        return index[key]

    def add_chain(pts, closed: bool) -> None:
        ids = [add_vertex(p) for p in pts]
        pairs = list(zip(ids, ids[1:]))
        if closed:
            pairs.append((ids[-1], ids[0]))
        for a, b in pairs:
            if a != b:
                segs.append([a, b])

    def linspace_edge(p0, p1, spacing):
        dist = float(np.linalg.norm(np.asarray(p1) - np.asarray(p0)))
        n = max(2, int(np.ceil(dist / spacing)))
        t = np.linspace(0.0, 1.0, n)
        return (1 - t)[:, None] * np.asarray(p0, float) + t[:, None] * np.asarray(p1, float)

    corners = [(0.0, 0.0), (LX, 0.0), (LX, LY), (0.0, LY)]
    for i in range(4):
        add_chain(linspace_edge(corners[i], corners[(i + 1) % 4], h), closed=False)
    center = origin_to_mesh(np.array([0.0, 0.0]))
    def add_circle(radius: float) -> None:
        n = max(32, int(np.ceil(2 * np.pi * radius / h_edge)))
        ang = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
        ring = center + radius * np.column_stack([np.cos(ang), np.sin(ang)])
        add_chain(ring, closed=True)

    if kind == "disk":
        add_circle(R_DISK)
    elif kind == "bulb":
        # Production bulb: plasma core, vacuum gap, 1 mm quartz shell.
        add_circle(R_DISK)
        add_circle(float(sc.r_bulb_inner))
        add_circle(float(sc.r_bulb_outer))
    elif kind == "square":
        s = square_half()
        c = center
        poly = [
            c + np.array([-s, -s]),
            c + np.array([s, -s]),
            c + np.array([s, s]),
            c + np.array([-s, s]),
        ]
        loop = []
        for i in range(4):
            edge = linspace_edge(poly[i], poly[(i + 1) % 4], h_edge)
            loop.append(edge[:-1])
        add_chain(np.vstack(loop), closed=True)
    else:
        raise ValueError(kind)
    amax = 0.45 * h * h
    # Triangle does not parse scientific notation after the area flag.
    mesh = tr.triangulate(
        {"vertices": np.asarray(verts, float), "segments": np.asarray(segs, np.int32)},
        f"pq20a{amax:.8f}",
    )
    pts = np.asarray(mesh["vertices"], float)
    tris = np.asarray(mesh["triangles"], int)
    a = pts[tris[:, 1]] - pts[tris[:, 0]]
    b = pts[tris[:, 2]] - pts[tris[:, 0]]
    area = 0.5 * np.abs(a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]).sum()
    if abs(area - LX * LY) / (LX * LY) > 0.02:
        raise RuntimeError(f"mesh area {area} != domain {LX * LY}")
    if len(pts) < 0.2 * (LX / h) * (LY / h):
        raise RuntimeError(f"mesh has only {len(pts)} nodes; the area constraint did not refine")
    return pts, tris


def inclusion_mask(points, tris, kind: str) -> np.ndarray:
    cents = points[tris].mean(1)
    origin = mesh_to_origin(cents)
    r = np.sqrt(np.sum(origin**2, axis=1))
    if kind == "disk":
        return r <= R_DISK
    if kind == "bulb":
        # Mask is unused for multi-material assignment; plasma core only.
        return r <= R_DISK
    s = square_half()
    return (np.abs(origin[:, 0]) <= s) & (np.abs(origin[:, 1]) <= s)


def rho_on(points, tris, kind: str, filled: bool):
    T = len(tris)
    rho = [
        np.ones(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.ones(T, dtype=np.complex128),
    ]
    if kind == "bulb":
        origin = mesh_to_origin(points[tris].mean(1))
        rad = np.sqrt(np.sum(origin**2, axis=1))
        plasma = rad <= R_DISK
        quartz = (rad > float(sc.r_bulb_inner)) & (rad <= float(sc.r_bulb_outer))
        if filled:
            rxx, rxy, ryx, ryy = rho_from_eps(*eps_tensor_at_bias(0.0))
            rho[0][plasma] = rxx
            rho[1][plasma] = rxy
            rho[2][plasma] = ryx
            rho[3][plasma] = ryy
            rho[0][quartz] = 1.0 / 3.8
            rho[3][quartz] = 1.0 / 3.8
        return rho, plasma
    mask = inclusion_mask(points, tris, kind)
    if filled:
        rxx, rxy, ryx, ryy = rho_from_eps(*eps_tensor_at_bias(0.0))
        rho[0][mask] = rxx
        rho[1][mask] = rxy
        rho[2][mask] = ryx
        rho[3][mask] = ryy
    return rho, mask


def diameter_points(n: int = 801) -> np.ndarray:
    x = np.linspace(-4.0, 4.0, n)
    return np.column_stack([x, np.zeros_like(x)])


def contour_points(half: float = 1.0, n_side: int = 81):
    """Closed square, origin-centered, traversed counterclockwise. Returns xy and outward normals."""
    s = np.linspace(-half, half, n_side)
    # bottom y=-half, outward normal -y, then right, top, left
    parts = []
    normals = []
    parts.append(np.column_stack([s[:-1], np.full(n_side - 1, -half)]))
    normals.append(np.tile(np.array([0.0, -1.0]), (n_side - 1, 1)))
    parts.append(np.column_stack([np.full(n_side - 1, half), s[:-1]]))
    normals.append(np.tile(np.array([1.0, 0.0]), (n_side - 1, 1)))
    parts.append(np.column_stack([s[:-1][::-1], np.full(n_side - 1, half)]))
    normals.append(np.tile(np.array([0.0, 1.0]), (n_side - 1, 1)))
    parts.append(np.column_stack([np.full(n_side - 1, -half), s[:-1][::-1]]))
    normals.append(np.tile(np.array([-1.0, 0.0]), (n_side - 1, 1)))
    return np.vstack(parts), np.vstack(normals)


def poynting_out(hz, ex, ey, normals, dl: float) -> float:
    """Outward power. Sx=1/2 Re(Ey conj(Hz)), Sy=-1/2 Re(Ex conj(Hz))."""
    sx = 0.5 * np.real(ey * np.conj(hz))
    sy = -0.5 * np.real(ex * np.conj(hz))
    return float(np.sum((sx * normals[:, 0] + sy * normals[:, 1]) * dl))


def solve_pair(points, tris, kind: str):
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = LX, LY, DPML
    k0 = 2 * np.pi * float(sc.fs_a)
    pec = np.zeros(len(points), dtype=bool)
    sampler = ElementSampler(points, tris)
    src = origin_to_mesh(SRC_XY)
    b = np.zeros(len(points), dtype=np.complex128)
    t = int(sampler.locate(src[None, :])[0])
    if t < 0:
        raise RuntimeError("source missed the mesh")
    w = sampler._bary(t, src)
    for k in range(3):
        b[int(tris[t, k])] += complex(w[k])
    out = {}
    for filled, tag in ((False, "vac"), (True, "disk")):
        rho, mask = rho_on(points, tris, kind, filled)
        A = assemble_anisotropic(points, tris, *rho, k0, pec)
        x, _lu, t_fac, t_sol, resid = solve_system(A, b, pec)
        if tag == "disk":
            op = dict(LAST_OPERATOR)
            if op.get("pml_profile") != "meep_quadratic_R1e-15" or op["sigma_max"] < 20:
                raise RuntimeError(f"PML fell back: {op}")
            if abs(op["sigma_max"] - meep_sigma_max(DPML)) > 1e-9:
                raise RuntimeError("sigma_max mismatch")
            out["operator"] = op
            out["factor_s"] = t_fac
            out["solve_s"] = t_sol
            out["residual"] = resid
            out["n_inclusion"] = int(mask.sum())
        xy = diameter_points()
        hz, ex, ey = sampler.fields(x, *rho, k0, origin_to_mesh(xy))
        out[f"Hz_{tag}"] = hz
        out[f"Ex_{tag}"] = ex
        out[f"Ey_{tag}"] = ey
        cxy, normals = contour_points()
        chz, cex, cey = sampler.fields(x, *rho, k0, origin_to_mesh(cxy))
        dl = 2.0 / 80.0
        out[f"P_out_{tag}"] = poynting_out(chz, cex, cey, normals, dl)
    out["x"] = diameter_points()[:, 0]
    if sampler.misses:
        raise RuntimeError(f"field samples missed the mesh: {sampler.misses}")
    return out


def fem_main(kind: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    levels = [("FEM-C", 0.030, 0.012), ("FEM-F", 0.018, 0.007)]
    rows = []
    for name, h, he in levels:
        print(f"=== FEM {kind} {name} ===", flush=True)
        t0 = time.perf_counter()
        points, tris = build_mesh(h, he, kind)
        t_mesh = time.perf_counter() - t0
        sol = solve_pair(points, tris, kind)
        rec = {
            "kind": kind,
            "level": name,
            "h": h,
            "h_edge": he,
            "n_nodes": int(len(points)),
            "n_tris": int(len(tris)),
            "mesh_s": t_mesh,
            "r_disk": R_DISK,
            "square_half": float(square_half()),
            "eps": [plasma_eps().real, plasma_eps().imag],
            "n_inclusion": sol["n_inclusion"],
            "factor_s": sol["factor_s"],
            "solve_s": sol["solve_s"],
            "residual": sol["residual"],
            "P_out_vac": sol["P_out_vac"],
            "P_out_disk": sol["P_out_disk"],
            "operator_pml": sol["operator"]["pml_profile"],
            "sigma_max": sol["operator"]["sigma_max"],
        }
        ratios = field_ratios(sol["x"], sol["Hz_vac"], sol["Hz_disk"])
        rec.update(ratios)
        print(name, "DOFs", rec["n_nodes"], "T", rec["T_forward"], "R", rec["R_back"], flush=True)
        rows.append(rec)
        np.savez_compressed(
            OUT / f"fem_{kind}_{name}.npz",
            x=sol["x"],
            Hz_vac=sol["Hz_vac"],
            Hz_disk=sol["Hz_disk"],
            Ex_disk=sol["Ex_disk"],
            Ey_disk=sol["Ey_disk"],
            Ex_vac=sol["Ex_vac"],
            Ey_vac=sol["Ey_vac"],
        )
        (OUT / f"fem_{kind}.json").write_text(json.dumps(rows, indent=2) + "\n")
    if len(rows) == 2:
        for key in ("T_forward", "R_back", "inside_over_vac", "just_out_over_vac"):
            a = complex(*rows[0][key])
            b = complex(*rows[1][key])
            if abs(a) > 0:
                print(f"{key} fine/coarse dB", float(20 * np.log10(abs(b) / abs(a))), flush=True)


def value_at(x, z, x0) -> complex:
    return complex(np.interp(x0, x, z.real) + 1j * np.interp(x0, x, z.imag))


def field_ratios(x, hz_vac, hz_disk) -> dict:
    def amp(z, x0):
        return value_at(x, z, x0)

    vac_f = amp(hz_vac, 3.0)
    disk_f = amp(hz_disk, 3.0)
    vac_b = amp(hz_vac, -3.0)
    disk_b = amp(hz_disk, -3.0)
    vac_c = amp(hz_vac, 0.0)
    disk_c = amp(hz_disk, 0.0)
    vac_o = amp(hz_vac, 0.30)
    disk_o = amp(hz_disk, 0.30)
    inside = (np.abs(x) < 0.12)
    return {
        "T_forward": [ (disk_f / vac_f).real, (disk_f / vac_f).imag ],
        "R_back": [ ((disk_b - vac_b) / vac_b).real, ((disk_b - vac_b) / vac_b).imag ],
        "inside_over_vac": [ (disk_c / vac_c).real, (disk_c / vac_c).imag ],
        "just_out_over_vac": [ (disk_o / vac_o).real, (disk_o / vac_o).imag ],
        "mean_abs_Hz_inside": float(np.mean(np.abs(hz_disk[inside]))),
        "mean_abs_Hz_vac_inside": float(np.mean(np.abs(hz_vac[inside]))),
    }


def mie_coefficients(k0: float, eps: complex, radius: float, nmax: int = 12):
    """Hz cylinder coefficients for e^{-iωt}.

    Outside: i^n [J_n(k0 r) + a_n H_n(k0 r)] e^{inθ}
    Continuity of Hz and (1/ε) ∂Hz/∂r. ε_out = 1.
    """
    from scipy.special import hankel1, h1vp, jv, jvp

    a = np.zeros(2 * nmax + 1, dtype=np.complex128)
    ns = np.arange(-nmax, nmax + 1)
    kp = k0 * np.sqrt(eps)
    for i, n in enumerate(ns):
        n = int(n)
        j = jv(n, k0 * radius)
        jp = jvp(n, k0 * radius, 1)
        h = hankel1(n, k0 * radius)
        hp = h1vp(n, k0 * radius, 1)
        jin = jv(n, kp * radius)
        jnp = jvp(n, kp * radius, 1)
        right = (kp / eps) * (jnp / jin)
        a[i] = (right * j - k0 * jp) / (k0 * hp - right * h)
    return ns, a


def mie_self_check() -> dict:
    k0 = 2 * np.pi * float(sc.fs_a)
    eps = plasma_eps()
    ns, a_air = mie_coefficients(k0, 1.0 + 0j, R_DISK)
    ns, a = mie_coefficients(k0, eps, R_DISK)
    return {
        "air_max_abs_a": float(np.max(np.abs(a_air))),
        "plasma_max_abs_a": float(np.max(np.abs(a))),
        "n_terms": int(len(ns)),
        "eps": [eps.real, eps.imag],
        "ka": float(k0 * R_DISK),
    }


def meep_main(kind: str) -> None:
    import meep as mp
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.Get_rank()
    ppc = float(os.environ.get("DISK_PPC", "25"))
    res = int(round(ppc * float(sc.a) * 100.0))  # a=0.02 m → res = 2*ppc
    # physical_units: res = a_m / dx_m, dx_m = 0.01/ppc, a_m=0.02 → res = 2*ppc
    averaging = os.environ.get("DISK_AVG", "1") != "0"
    until = float(os.environ.get("DISK_UNTIL", "20"))
    fs = float(sc.fs_a)
    fp = float(sc.fp_a)
    gamma = float(sc.gamma_a)
    cell = mp.Vector3(LX, LY)
    plasma = mp.Medium(
        epsilon=1.0,
        E_susceptibilities=[mp.DrudeSusceptibility(frequency=fp, gamma=gamma, sigma=1.0)],
    )
    if kind == "disk":
        geom = [mp.Cylinder(radius=R_DISK, center=mp.Vector3(0, 0), material=plasma)]
    elif kind == "bulb":
        quartz = mp.Medium(epsilon=3.8)
        vacuum = mp.Medium(epsilon=1.0)
        geom = [
            mp.Cylinder(radius=float(sc.r_bulb_outer), center=mp.Vector3(0, 0), material=quartz),
            mp.Cylinder(radius=float(sc.r_bulb_inner), center=mp.Vector3(0, 0), material=vacuum),
            mp.Cylinder(radius=R_DISK, center=mp.Vector3(0, 0), material=plasma),
        ]
    elif kind == "square":
        s = float(square_half())
        geom = [
            mp.Block(
                size=mp.Vector3(2 * s, 2 * s, mp.inf),
                center=mp.Vector3(0, 0),
                material=plasma,
            )
        ]
    else:
        raise ValueError(kind)
    # Confirm the Drude frequency before time stepping.
    sus = plasma.E_susceptibilities[0]
    if abs(float(sus.frequency) - fp) > 1e-9:
        raise RuntimeError(f"Drude frequency {sus.frequency} != fp {fp}")
    eps_m = complex(plasma.epsilon(fs)[0, 0])
    if abs(eps_m - plasma_eps()) > 1e-8:
        raise RuntimeError(f"Meep Drude epsilon {eps_m} != FEM {plasma_eps()}")

    def one(material_geometry):
        sim = mp.Simulation(
            cell_size=cell,
            resolution=res,
            boundary_layers=[mp.PML(DPML)],
            geometry=material_geometry,
            sources=[
                mp.Source(
                    mp.GaussianSource(fs, fwidth=0.20 * fs),
                    component=mp.Hz,
                    center=mp.Vector3(SRC_XY[0], SRC_XY[1]),
                )
            ],
            default_material=mp.Medium(epsilon=1, mu=1),
            eps_averaging=averaging,
        )
        vol = mp.Volume(center=mp.Vector3(0, 0), size=mp.Vector3(8.2, 4.2))
        dft = sim.add_dft_fields([mp.Hz, mp.Ex, mp.Ey], fs, 0, 1, where=vol)
        t0 = time.perf_counter()

        def _probe(sim_):
            pts = ((-0.23, 0.0), (0.0, 0.0), (0.23, 0.0), (3.0, 0.0))
            vals = []
            for px, py in pts:
                hz = sim_.get_field_point(mp.Hz, mp.Vector3(px, py))
                vals.append(f"{abs(hz):.4e}")
            if rank == 0:
                print(
                    f"PROBE t={sim_.meep_time():.2f} |Hz|=" + ",".join(vals),
                    flush=True,
                )

        step = [mp.at_every(10, _probe)] if os.environ.get("DISK_PROBE") == "1" else []
        sim.run(*step, until_after_sources=until)
        wall = time.perf_counter() - t0
        fields = {}
        for comp, name in ((mp.Hz, "Hz"), (mp.Ex, "Ex"), (mp.Ey, "Ey")):
            raw = np.array(sim.get_dft_array(dft, comp, 0))
            fields[name] = raw
        xs, ys, *_ = sim.get_array_metadata(dft_cell=dft)
        return fields, np.asarray(xs, float), np.asarray(ys, float), wall

    if rank == 0:
        print(f"MEEP {kind} ppc={ppc} res={res} avg={averaging} eps={eps_m}", flush=True)
    vac, xs, ys, wall_vac = one([])
    disk, xs2, ys2, wall_disk = one(geom)
    if rank != 0:
        return
    tag = f"meep_{kind}_ppc{ppc:g}" + ("" if averaging else "_noavg")
    if abs(until - 20.0) > 1e-9:
        tag += f"_u{until:g}"
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        OUT / f"{tag}.npz",
        xs=xs,
        ys=ys,
        Hz_vac=vac["Hz"],
        Hz_disk=disk["Hz"],
        Ex_disk=disk["Ex"],
        Ey_disk=disk["Ey"],
        Ex_vac=vac["Ex"],
        Ey_vac=vac["Ey"],
    )
    # Sample the diameter from the DFT grid for a quick ratio.
    x_line = np.linspace(-4.0, 4.0, 801)
    hz_v = interp_grid(xs, ys, vac["Hz"], x_line, np.zeros_like(x_line))
    hz_d = interp_grid(xs, ys, disk["Hz"], x_line, np.zeros_like(x_line))
    ratios = field_ratios(x_line, hz_v, hz_d)
    rec = {
        "kind": kind,
        "ppc": ppc,
        "res": res,
        "eps_averaging": averaging,
        "until_after_sources": until,
        "eps": [eps_m.real, eps_m.imag],
        "wall_vac_s": wall_vac,
        "wall_disk_s": wall_disk,
        "dx_mm": 10.0 / ppc,
        **ratios,
    }
    (OUT / f"{tag}.json").write_text(json.dumps(rec, indent=2) + "\n")
    print("RATIOS", json.dumps(ratios), flush=True)
    print(f"{tag} EXIT", flush=True)


def interp_grid(xs, ys, field, xq, yq):
    field = np.asarray(field)
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    if field.shape == (len(ys), len(xs)):
        field = field.T
    elif field.shape != (len(xs), len(ys)):
        raise RuntimeError(f"DFT shape {field.shape} vs nx {len(xs)} ny {len(ys)}")
    ix = np.interp(xq, xs, np.arange(len(xs)))
    iy = np.interp(yq, ys, np.arange(len(ys)))
    i0 = np.clip(np.floor(ix).astype(int), 0, len(xs) - 2)
    j0 = np.clip(np.floor(iy).astype(int), 0, len(ys) - 2)
    tx = np.clip(ix - i0, 0, 1)
    ty = np.clip(iy - j0, 0, 1)
    return (
        (1 - tx) * (1 - ty) * field[i0, j0]
        + tx * (1 - ty) * field[i0 + 1, j0]
        + (1 - tx) * ty * field[i0, j0 + 1]
        + tx * ty * field[i0 + 1, j0 + 1]
    )


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) > 1 else "fem"
    kind = sys.argv[2] if len(sys.argv) > 2 else "disk"
    if mode == "fem":
        check = mie_self_check()
        print("MIE", json.dumps(check), flush=True)
        if check["air_max_abs_a"] > 1e-8:
            raise RuntimeError(f"Mie air coefficients are not zero: {check}")
        (OUT).mkdir(parents=True, exist_ok=True)
        (OUT / "mie_check.json").write_text(json.dumps(check, indent=2) + "\n")
        fem_main(kind)
        return 0
    if mode == "meep":
        meep_main(kind)
        return 0
    raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    raise SystemExit(main())
