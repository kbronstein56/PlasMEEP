#!/usr/bin/env python
"""
Small two-port dielectric-waveguide calibration for Add_Port conventions.

Uses existing Add_Horn / Add_Port only. No plasma, no scaffold, no
five-port assembly. One frequency. Short FDTD run.

Cases:
  --aligned   axis-aligned guide (default)
  --oblique   72-degree guide (five-port angular spacing)

Prints modal coefficients, flux, and stored-normalization S-parameters.
"""

import argparse
import math
import sys

import meep as mp
import numpy as np

from plasmeep.lib import (
    PORT_AWAY_FROM_APERTURE,
    PORT_TOWARD_APERTURE,
    Plasmeep,
)
from plasmeep.normalization import IncidentNormalization


A = 0.01
EPS = 4.0
FEED_W = 0.50
FEED_L = 8.0
FCEN = 0.25
FWIDTH = 0.20 * FCEN
EIG_BAND = 1
RES = 16
DPML = 2.0
SOURCE_OFFSET = 1.0


def _v3(arr):
    arr = np.asarray(arr, dtype=float).reshape(-1)
    z = float(arr[2]) if arr.size > 2 else 0.0
    return mp.Vector3(float(arr[0]), float(arr[1]), z)


def _xyz(vec):
    return np.array([vec.x, vec.y, vec.z], dtype=float)


def build_two_port(theta, dpml=DPML, res=RES, nx=None, ny=None):
    """
    Two facing straight feeds forming one uniform dielectric guide.

    Input horn at - (FEED_L/2) * u, axis +u (toward the joint / device).
    Output horn at + (FEED_L/2) * u, axis -u (toward the joint / device).
    The guide occupies [-FEED_L, +FEED_L] along u and extends into PML
    when the cell is sized from the feed ends plus dpml.
    """
    u = np.array([math.cos(theta), math.sin(theta), 0.0])
    n = np.array([-math.sin(theta), math.cos(theta), 0.0])
    half = 0.5 * FEED_L
    # Along-axis: cell edge at the feed ends so dielectric enters PML.
    # Transverse: extra air so the snapped port line and feed walls
    # are not sitting in the side PML (important at 72 deg).
    along = FEED_L * np.abs(u)
    trans = 0.5 * FEED_W * np.abs(n) + dpml + 1.5
    if nx is None:
        nx = 2.0 * max(along[0], trans[0])
    if ny is None:
        ny = 2.0 * max(along[1], trans[1])

    model = Plasmeep(A, res, dpml, nx, ny, verbose=False)
    horn_in = model.Add_Horn(-half * u[:2], theta, FEED_W, FEED_L, EPS)
    horn_out = model.Add_Horn(
        half * u[:2], theta + math.pi, FEED_W, FEED_L, EPS
    )
    common = dict(
        frequency=FCEN, fwidth=FWIDTH, eig_band=EIG_BAND, eig_parity=mp.ODD_Z,
        source_offset=SOURCE_OFFSET,
    )
    port_in = model.Add_Port(horn_in, PORT_TOWARD_APERTURE, **common)
    port_out = model.Add_Port(horn_out, PORT_TOWARD_APERTURE, **common)
    model.Add_Port_Source(port_in)
    return model, horn_in, horn_out, port_in, port_out


def _mode_region(port):
    return mp.ModeRegion(
        center=_v3(port["port_center"]), size=_v3(port["volume_size"])
    )


def _kfunc(kvec):
    k = mp.Vector3(float(kvec[0]), float(kvec[1]), float(kvec[2]))
    return lambda freq, mode: k


def measure_alphas(sim, dft_mon, port, eig_parity):
    """Both overlap indices for toward and away kpoint_func."""
    out = {}
    for name, kvec in (
        ("toward", port["kpoint_toward_aperture"]),
        ("away", port["kpoint_away_from_aperture"]),
    ):
        res = sim.get_eigenmode_coefficients(
            dft_mon,
            [port["eig_band"]],
            direction=mp.NO_DIRECTION,
            kpoint_func=_kfunc(kvec),
            eig_parity=eig_parity,
        )
        out[name] = {
            "idx0": complex(res.alpha[0, 0, 0]),
            "idx1": complex(res.alpha[0, 0, 1]),
            "kpoint": _xyz(res.kpoints[0]),
        }
    return out


def run_case(label, theta, dpml=DPML, extra_pml=0.0):
    dpml = float(dpml) + float(extra_pml)
    model, horn_in, horn_out, port_in, port_out = build_two_port(
        theta, dpml=dpml
    )
    sim = model.Get_Sim()
    mon_in = sim.add_mode_monitor(FCEN, 0, 1, _mode_region(port_in))
    mon_out = sim.add_mode_monitor(FCEN, 0, 1, _mode_region(port_out))
    flux_in = sim.add_flux(FCEN, 0, 1, _mode_region(port_in))
    flux_out = sim.add_flux(FCEN, 0, 1, _mode_region(port_out))

    mid = 0.5 * (horn_in["feed_front"] + horn_out["feed_front"])
    # Downstream of the source along +input axis; independent launch check.
    down = horn_in["port_center"] + 1.25 * horn_in["axis"]
    up = horn_in["port_center"] - 1.25 * horn_in["axis"]
    flux_down = sim.add_flux(
        FCEN, 0, 1,
        mp.FluxRegion(center=_v3(down), size=_v3(port_in["volume_size"])),
    )
    flux_up = sim.add_flux(
        FCEN, 0, 1,
        mp.FluxRegion(center=_v3(up), size=_v3(port_in["volume_size"])),
    )
    flux_mid = sim.add_flux(
        FCEN, 0, 1,
        mp.FluxRegion(center=_v3(mid), size=_v3(port_in["volume_size"])),
    )

    decay_pt = _v3(horn_out["port_center"])
    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(20, mp.Ez, decay_pt, 1e-6)
    )

    parity = port_in["eig_parity"]
    a_in = measure_alphas(sim, mon_in, port_in, parity)
    a_out = measure_alphas(sim, mon_out, port_out, parity)

    # Add_Port_Monitor path (same kpoint_func / overlap-index 0).
    adj_in_plus = model.Add_Port_Monitor(port_in, sim, PORT_TOWARD_APERTURE)
    adj_in_minus = model.Add_Port_Monitor(port_in, sim, PORT_AWAY_FROM_APERTURE)
    adj_out_minus = model.Add_Port_Monitor(port_out, sim, PORT_AWAY_FROM_APERTURE)
    for mon in (adj_in_plus, adj_in_minus, adj_out_minus):
        mon.register_monitors(np.array([FCEN]))
    # DFT monitors just created are empty; reuse the already-filled raw alphas
    # for physics. Cross-check only that kpoint_func matches.
    k_plus = _xyz(adj_in_plus.kpoint_func(FCEN, EIG_BAND))
    k_minus = _xyz(adj_in_minus.kpoint_func(FCEN, EIG_BAND))

    c_plus = a_in["toward"]["idx0"]
    c_refl = a_in["away"]["idx0"]
    c_tran = a_out["away"]["idx0"]
    c_out_toward = a_out["toward"]["idx0"]

    norm = IncidentNormalization.from_port(
        port_in, np.array([c_plus], dtype=complex), source_id="in"
    )
    S11 = norm.normalize(np.array([c_refl], dtype=complex), source_id="in")[0]
    S21 = norm.normalize(np.array([c_tran], dtype=complex), source_id="in")[0]

    f_in = mp.get_fluxes(flux_in)[0]
    f_out = mp.get_fluxes(flux_out)[0]
    f_down = mp.get_fluxes(flux_down)[0]
    f_up = mp.get_fluxes(flux_up)[0]
    f_mid = mp.get_fluxes(flux_mid)[0]

    result = {
        "label": label,
        "theta": theta,
        "dpml": dpml,
        "cell": (float(model.cell.x), float(model.cell.y)),
        "axis_in": horn_in["axis"].copy(),
        "c_plus": c_plus,
        "c_refl": c_refl,
        "c_tran": c_tran,
        "c_out_toward": c_out_toward,
        "alphas_in": a_in,
        "alphas_out": a_out,
        "S11": S11,
        "S21": S21,
        "flux_in": f_in,
        "flux_out": f_out,
        "flux_down": f_down,
        "flux_up": f_up,
        "flux_mid": f_mid,
        "kpoint_monitor_plus": k_plus,
        "kpoint_monitor_minus": k_minus,
        "launch_k_dot_axis": float(
            a_in["toward"]["kpoint"] @ horn_in["axis"]
        ),
    }
    return result


def _fmt_c(z):
    return "{: .6e} {:+.6e}j  |z|={:.6e}".format(z.real, z.imag, abs(z))


def report(result):
    print("=" * 72)
    print("{}  theta={:.3f} rad ({:.1f} deg)  dpml={}  cell={}".format(
        result["label"], result["theta"], math.degrees(result["theta"]),
        result["dpml"], result["cell"],
    ))
    print("input axis u =", result["axis_in"])
    print()
    print("Driven port, toward_aperture idx0  (intended c_i^+):")
    print(" ", _fmt_c(result["c_plus"]))
    print("Driven port, away_from_aperture idx0 (intended reflection):")
    print(" ", _fmt_c(result["c_refl"]))
    print("Receive port, away_from_aperture idx0 (intended c_j^-):")
    print(" ", _fmt_c(result["c_tran"]))
    print("Receive port, toward_aperture idx0 (should be small):")
    print(" ", _fmt_c(result["c_out_toward"]))
    print()
    print("Raw overlap pairs at input (toward k):")
    print("  idx0", _fmt_c(result["alphas_in"]["toward"]["idx0"]))
    print("  idx1", _fmt_c(result["alphas_in"]["toward"]["idx1"]))
    print("  MPB k · u = {:.4f}".format(result["launch_k_dot_axis"]))
    print("Raw overlap pairs at input (away k):")
    print("  idx0", _fmt_c(result["alphas_in"]["away"]["idx0"]))
    print("  idx1", _fmt_c(result["alphas_in"]["away"]["idx1"]))
    print()
    print("Flux (MEEP +cartesian normal of each region):")
    print("  input port plane {: .6e}".format(result["flux_in"]))
    print("  output port plane {: .6e}".format(result["flux_out"]))
    print("  downstream of src {: .6e}".format(result["flux_down"]))
    print("  upstream of src   {: .6e}".format(result["flux_up"]))
    print("  mid-guide         {: .6e}".format(result["flux_mid"]))
    print()
    print("|c_plus|^2        {:.6e}".format(abs(result["c_plus"]) ** 2))
    print("|c_tran|^2        {:.6e}".format(abs(result["c_tran"]) ** 2))
    print("|c_refl|^2        {:.6e}".format(abs(result["c_refl"]) ** 2))
    print("|S11|={:.4f}  |S21|={:.4f}  |S11|^2+|S21|^2={:.4f}".format(
        abs(result["S11"]), abs(result["S21"]),
        abs(result["S11"]) ** 2 + abs(result["S21"]) ** 2,
    ))
    print("S11", _fmt_c(result["S11"]))
    print("S21", _fmt_c(result["S21"]))
    print("monitor k_plus · u  {:.4f}".format(
        result["kpoint_monitor_plus"] @ result["axis_in"]
    ))
    print("monitor k_minus · u {:.4f}".format(
        result["kpoint_monitor_minus"] @ result["axis_in"]
    ))
    print("=" * 72)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--aligned", action="store_true", help="run the axis-aligned case"
    )
    parser.add_argument(
        "--oblique", action="store_true", help="run the 72-degree case"
    )
    parser.add_argument(
        "--extra-pml", type=float, default=0.0,
        help="add this thickness to the default PML (oblique reflection check)"
    )
    args = parser.parse_args(argv)
    if not args.aligned and not args.oblique:
        args.aligned = True
        args.oblique = True

    mp.verbosity(1)
    results = []
    if args.aligned:
        results.append(report(run_case("aligned", 0.0, extra_pml=args.extra_pml)))
    if args.oblique:
        results.append(report(run_case(
            "oblique-72", 2.0 * math.pi / 5.0, extra_pml=args.extra_pml
        )))
    return results


if __name__ == "__main__":
    main()
