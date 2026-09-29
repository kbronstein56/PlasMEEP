#!/usr/bin/env python3
"""Compare the single-inclusion FEM and Meep fields and draw the diameter plots."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "plasma_interface"


def cplx(pair):
    return complex(pair[0], pair[1])


def db_amp(a, b):
    """Power-equivalent dB between two complex field ratios."""
    if abs(b) == 0:
        return None
    return float(20 * np.log10(abs(a) / abs(b)))


def load_fem(kind, level):
    rows = json.loads((OUT / f"fem_{kind}.json").read_text())
    rec = next(r for r in rows if r["level"] == level)
    data = np.load(OUT / f"fem_{kind}_{level}.npz")
    return rec, data


def load_meep(kind, ppc, noavg=False):
    tag = f"meep_{kind}_ppc{ppc:g}" + ("_noavg" if noavg else "")
    rec = json.loads((OUT / f"{tag}.json").read_text())
    data = np.load(OUT / f"{tag}.npz")
    return rec, data


def diameter_from_meep(data):
    import sys

    sys.path[:0] = [str(ROOT / "scripts" / "validation" / "fem_meep_validation")]
    from plasma_interface import interp_grid

    x = np.linspace(-4.0, 4.0, 801)
    y = np.zeros_like(x)
    return {
        "x": x,
        "Hz_vac": interp_grid(data["xs"], data["ys"], data["Hz_vac"], x, y),
        "Hz_disk": interp_grid(data["xs"], data["ys"], data["Hz_disk"], x, y),
        "Ex_disk": interp_grid(data["xs"], data["ys"], data["Ex_disk"], x, y),
        "Ey_disk": interp_grid(data["xs"], data["ys"], data["Ey_disk"], x, y),
    }


def scale_to_fem(x, fem_vac, meep_vac):
    window = (x > -4.0) & (x < -1.5)
    a = meep_vac[window]
    b = fem_vac[window]
    scale = np.vdot(a, b) / np.vdot(a, a)
    rel = float(np.linalg.norm(scale * a - b) / np.linalg.norm(b))
    return complex(scale), rel


def mie_ratio(x):
    """Point-source Mie ratio Hz(with disk)/Hz(vacuum) on the x axis.

    The source sits at x=-6. Graf's addition theorem expands that outgoing
    Hankel wave in regular Bessel waves about the disk, and each partial
    wave uses the same cylinder coefficient as the plane-wave series.
    """
    import sys

    sys.path[:0] = [
        str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
        str(ROOT / "scripts" / "validation"),
        str(ROOT / "scripts"),
    ]
    import sixport_common as sc
    from plasma_interface import R_DISK, mie_coefficients, plasma_eps
    from scipy.special import hankel1, jv

    k0 = 2 * np.pi * float(sc.fs_a)
    eps = plasma_eps()
    ns, coeff = mie_coefficients(k0, eps, R_DISK, nmax=40)
    kp = k0 * np.sqrt(complex(eps))
    rs = 6.0
    cn = np.array([hankel1(int(n), k0 * rs) * ((-1) ** int(n)) for n in ns])
    ratio = np.zeros(len(x), dtype=np.complex128)
    for i, xv in enumerate(x):
        r = abs(float(xv))
        sign = 1.0 if xv >= 0.0 else -1.0
        vac = hankel1(0, k0 * abs(float(xv) + rs))
        tot = 0j
        for n, a, c in zip(ns, coeff, cn):
            phase = sign ** int(n)
            if r >= R_DISK:
                tot += phase * c * (jv(int(n), k0 * r) + a * hankel1(int(n), k0 * r))
            else:
                j_b = jv(int(n), k0 * R_DISK)
                h_b = hankel1(int(n), k0 * R_DISK)
                jin_b = jv(int(n), kp * R_DISK)
                b = c * (j_b + a * h_b) / jin_b
                tot += phase * b * jv(int(n), kp * r)
        ratio[i] = tot / vac
    return ratio


def compare(kind: str) -> None:
    fem_rows = json.loads((OUT / f"fem_{kind}.json").read_text())
    fine = fem_rows[-1]
    print(f"\n{kind} FEM {fine['level']} DOFs {fine['n_nodes']}")
    keys = ("T_forward", "R_back", "inside_over_vac", "just_out_over_vac")
    if len(fem_rows) >= 2:
        for key in keys:
            a = cplx(fem_rows[-2][key])
            b = cplx(fem_rows[-1][key])
            print(f"  FEM {key} fine vs coarse {db_amp(b, a):+.4f} dB")
    meep_files = sorted(OUT.glob(f"meep_{kind}_ppc*.json"))
    for path in meep_files:
        rec = json.loads(path.read_text())
        print(f"  Meep {path.name}")
        for key in keys:
            dB = db_amp(cplx(rec[key]), cplx(fine[key]))
            rel = abs(cplx(rec[key]) - cplx(fine[key])) / max(abs(cplx(fine[key])), 1e-30)
            print(f"    {key:22s} Meep-FEM {dB:+.3f} dB   |Δ|/|FEM| {rel:.4f}")


def plot_disk() -> None:
    rec, fem = load_fem("disk", "FEM-F")
    x = fem["x"]
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.2), sharex=True)
    colors = {25: "C0", 35: "C1", 50: "C2"}
    mie = mie_ratio(x)
    fem_ratio = fem["Hz_disk"] / fem["Hz_vac"]
    axes[0, 0].plot(x, np.abs(fem_ratio), "k", lw=2, label="FEM-F")
    axes[0, 0].plot(x, np.abs(mie), "k--", lw=1, label="Mie")
    for ppc, color in colors.items():
        path = OUT / f"meep_disk_ppc{ppc}.json"
        if not path.exists():
            continue
        data = np.load(OUT / f"meep_disk_ppc{ppc}.npz")
        m = diameter_from_meep(data)
        ratio = m["Hz_disk"] / m["Hz_vac"]
        axes[0, 0].plot(m["x"], np.abs(ratio), color=color, label=f"Meep {ppc}")
        scale, rel = scale_to_fem(x, fem["Hz_vac"], m["Hz_vac"])
        err = fem["Hz_disk"] - scale * m["Hz_disk"]
        axes[1, 0].plot(x, np.abs(err) / np.maximum(np.abs(fem["Hz_disk"]), 1e-12), color=color, label=f"{ppc} rel {rel:.3f}")
        axes[0, 1].plot(x, np.abs(fem["Hz_disk"]), "k", lw=2)
        axes[0, 1].plot(m["x"], np.abs(scale * m["Hz_disk"]), color=color)
    axes[0, 0].axvline(0.23, color="0.6", lw=0.6)
    axes[0, 0].axvline(-0.23, color="0.6", lw=0.6)
    axes[0, 0].set_ylabel("|Hz disk / Hz vacuum|")
    axes[0, 0].set_title("Field ratio along the diameter")
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].set_xlim(-1.2, 1.2)
    axes[1, 0].set_xlim(-1.2, 1.2)
    axes[1, 0].set_ylabel("|FEM − scaled Meep| / |FEM|")
    axes[1, 0].set_xlabel("x (a)")
    axes[1, 0].legend(fontsize=8)
    axes[0, 1].set_xlim(-1.2, 1.2)
    axes[0, 1].set_ylabel("|Hz|")
    axes[0, 1].set_title("Total field, one complex scale from vacuum")
    # zoom of the ratio error near the interface
    axes[1, 1].plot(x, np.abs(fem_ratio - mie), "k", label="|FEM ratio − Mie|")
    axes[1, 1].axvline(0.23, color="0.6", lw=0.6)
    axes[1, 1].axvline(-0.23, color="0.6", lw=0.6)
    axes[1, 1].set_xlim(-0.6, 0.6)
    axes[1, 1].set_xlabel("x (a)")
    axes[1, 1].set_ylabel("|ratio error|")
    axes[1, 1].legend(fontsize=8)
    axes[1, 1].set_title("FEM versus the point-source Mie ratio")
    for ax in axes.ravel():
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "disk_diameter.png", dpi=140)
    print("wrote", OUT / "disk_diameter.png")
    # Where does FEM-Mie error jump?
    err = np.abs(fem_ratio - mie)
    for x0 in (-0.5, -0.23, 0.0, 0.23, 0.5, 2.0):
        i = int(np.argmin(np.abs(x - x0)))
        print(f"  |FEM-Mie| at x={x[i]:+.3f}: {err[i]:.4f}")


def plot_convergence() -> None:
    fem = json.loads((OUT / "fem_disk.json").read_text())[-1]
    keys = ("T_forward", "R_back", "inside_over_vac", "just_out_over_vac")
    rows = []
    for path in sorted(OUT.glob("meep_disk_ppc*.json")):
        rec = json.loads(path.read_text())
        if rec.get("until_after_sources", 20) not in (20, 20.0) and "u40" not in path.name and "u100" not in path.name:
            continue
        rows.append((path.name, rec))
    if not rows:
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    markers = {"T_forward": "o", "R_back": "s", "inside_over_vac": "^", "just_out_over_vac": "D"}
    for key in keys:
        xs, ys = [], []
        for name, rec in rows:
            if "noavg" in name or "u100" in name or "u40" in name:
                continue
            xs.append(rec["dx_mm"])
            ys.append(abs(cplx(rec[key]) - cplx(fem[key])))
        if xs:
            ax.plot(xs, ys, marker=markers[key], label=key)
    ax.set_xlabel("dx (mm)")
    ax.set_ylabel("|Meep − FEM-F| of the complex ratio")
    ax.set_title("Bare plasma disk, until_after_sources = 20")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "disk_convergence.png", dpi=140)
    print("wrote", OUT / "disk_convergence.png")


def main():
    for kind in ("disk", "bulb", "square"):
        if (OUT / f"fem_{kind}.json").exists():
            compare(kind)
    if (OUT / "fem_disk_FEM-F.npz").exists():
        plot_disk()
        plot_convergence()


if __name__ == "__main__":
    main()
