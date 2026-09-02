"""
The module composed in this file is meant to serve as a set of tools to make
simulating the EM response of plasma-based optical devices a bit easier.
Jesse A Rodriguez, 11/11/2021
"""

import glob
import h5py
import numpy as np
import matplotlib as mpl
from matplotlib import animation
import matplotlib.pylab as plt
import meep as mp
from meep import mpb
from meep import Harminv, after_sources
import os
from PIL import Image
import shutil
import sys
import scipy.special as sp
from scipy.optimize import brentq
import csv

###############################################################################
## Utility Functions and Globals
###############################################################################
c = 299792458
e = 1.60217662*10**(-19)
epso = 8.8541878128*10**(-12)
me = 9.1093837015*10**(-31)
J01 = 2.40482555768577
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',\
          '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

def WP(n):
    """
    Function for calculating plasma frequency given density

    n: electron density in m^(-3)
    """
    return (n*e**2/me/epso)**(1/2)

def n_e(wp):
    """
    Function for calculating electron density given plasma frequency

    wp: plasma frequency in rad/s
    """
    return wp**2*me*epso/e**2


def _as_float_array(x):
    """Return x as a float ndarray, preserving scalar-ness for callers."""
    return np.asarray(x, dtype=float)


def _is_scalar(x):
    return np.isscalar(x) or isinstance(x, np.ndarray) and x.ndim == 0


def _validate_beta(beta):
    """
    beta is the dimensionless wall argument in n_e(r) = n0 J0(beta r / R_p).

    Physical first-lobe range: (0, J01].
    beta -> 0+ is the flat-top (collisional) limit.
    beta = J01 is the Schottky solution (first J0 zero at the wall).
    """
    beta_arr = _as_float_array(beta)
    if np.any(~np.isfinite(beta_arr)):
        raise ValueError("beta must be finite.")
    if np.any(beta_arr <= 0.0) or np.any(beta_arr > J01):
        raise ValueError(
            "beta must lie in the first positive J0 lobe: 0 < beta <= J01 "
            "(J01 = {:.12f}). beta = 0 is the excluded flat-top limit; "
            "beta > J01 leaves the first lobe and can make n_e(r) negative."
            .format(J01)
        )
    return beta_arr


def _validate_h(h):
    """
    h is the edge-to-center density ratio, not a Bessel argument.

    h = n_e(R_p) / n0 = J0(beta), so h is in [0, 1).
    h = 0 is Schottky (beta = J01).
    h -> 1- is the flat-top limit (beta -> 0+).
    The Stage A/B prior samples h uniformly in about [0.05, 0.8].
    """
    h_arr = _as_float_array(h)
    if np.any(~np.isfinite(h_arr)):
        raise ValueError("h must be finite.")
    if np.any(h_arr < 0.0) or np.any(h_arr >= 1.0):
        raise ValueError(
            "h = J0(beta) = n_e(R_p)/n0 must satisfy 0 <= h < 1 on the first "
            "J0 lobe. h = 0 is Schottky; h -> 1 is the excluded flat-top limit. "
            "Do not pass beta as h."
        )
    return h_arr


def h_from_beta(beta):
    """
    Convert wall argument beta to edge-to-center ratio h.

    h = J0(beta) = n_e(R_p) / n0.

    beta: dimensionless, 0 < beta <= J01. Not a density ratio.
    Returns h in [0, 1). J0(J01) is flushed to 0 to absorb float error.
    """
    scalar = _is_scalar(beta)
    beta_arr = _validate_beta(beta)
    # J0 is nonnegative on the first lobe; clip the signed residual at J01.
    # PlasMEEP's J01 is the library Schottky wall argument (not scipy's jn_zeros).
    h = np.clip(sp.j0(beta_arr), 0.0, 1.0)
    h = np.where(np.isclose(beta_arr, J01, rtol=0.0, atol=1e-12), 0.0, h)
    if scalar:
        return float(np.asarray(h).reshape(()))
    return h


def beta_from_h(h):
    """
    Invert h = J0(beta) on the first positive J0 lobe.

    Returns the unique beta in (0, J01] such that J0(beta) = h.
    Uses scipy.optimize.brentq. h = 0 maps exactly to J01.

    h: edge-to-center ratio n_e(R_p)/n0, 0 <= h < 1. Not a Bessel argument.
    """
    scalar = _is_scalar(h)
    h_arr = _validate_h(h)
    out = np.empty(h_arr.shape, dtype=float)

    def _root(hv):
        if hv == 0.0:
            return J01
        # brentq on [0, J01]; xtol tightened because J0'(0) = 0 makes
        # the flat-top end mildly ill-conditioned.
        return brentq(lambda b: sp.j0(b) - hv, 0.0, J01, xtol=1e-14)

    for idx in np.ndindex(h_arr.shape):
        out[idx] = _root(h_arr[idx])
    if scalar:
        return float(out.reshape(()))
    return out


def Bessel_profile(r, n0, R_p, beta=None, h=None):
    """
    Radial Bessel / diffusion-family electron-density profile.

    n_e(r) = n0 * J0(beta * r / R_p)

    n0 is the on-axis density (same units as the return value, typically
    m^{-3}). This helper does not convert n_e to plasma frequency and does
    not nondimensionalize; use WP(n) and Nondimensionalize_Freq separately.

    Specify exactly one of:
        beta: dimensionless wall argument, 0 < beta <= J01.
              beta -> 0+ : flat-top. beta = J01 : Schottky (n_e(R_p) = 0).
        h:    edge-to-center ratio n_e(R_p)/n0 = J0(beta), 0 <= h < 1.
              The tomography prior samples h, not beta, uniformly in
              about [0.05, 0.8], then converts with beta_from_h.

    r and R_p must share length units. The physical column is 0 <= r <= R_p;
    coordinates outside that interval are rejected.
    """
    if (beta is None) == (h is None):
        raise ValueError("Specify exactly one of beta or h, not both or neither.")
    if h is not None:
        if not _is_scalar(h):
            raise ValueError("Bessel_profile expects a scalar h.")
        beta = beta_from_h(h)
    else:
        if not _is_scalar(beta):
            raise ValueError("Bessel_profile expects a scalar beta.")
        _validate_beta(beta)

    try:
        n0 = float(n0)
        R_p = float(R_p)
    except (TypeError, ValueError) as exc:
        raise ValueError("n0 and R_p must be real scalars.") from exc
    if not np.isfinite(n0) or n0 < 0.0:
        raise ValueError("n0 must be a finite on-axis density, n0 >= 0.")
    if not np.isfinite(R_p) or R_p <= 0.0:
        raise ValueError("R_p must be a finite plasma radius, R_p > 0.")

    scalar_r = _is_scalar(r)
    r_arr = _as_float_array(r)
    if np.any(~np.isfinite(r_arr)):
        raise ValueError("Radial coordinates r must be finite.")
    if np.any(r_arr < 0.0) or np.any(r_arr > R_p):
        raise ValueError(
            "Radial coordinates must satisfy 0 <= r <= R_p. "
            "Bessel_profile is the interior column n_e(r); it is not defined "
            "outside the plasma radius."
        )

    beta = float(np.asarray(beta, dtype=float).reshape(()))
    # Interior uses J0(beta r / R_p). The wall is set from h = J0(beta)
    # so n(R_p) = h n0 exactly, including Schottky (h = 0 at beta = J01),
    # where scipy.j0(J01) is a ~1e-12 residual times n0.
    n = n0 * sp.j0(beta * r_arr / R_p)
    n = np.maximum(n, 0.0)
    n = np.where(r_arr >= R_p, n0 * h_from_beta(beta), n)
    if scalar_r:
        return float(np.asarray(n).reshape(()))
    return n

def get_trilat_shift(loc,\
              a_basis = np.array([[np.sqrt(3)/2,0.5,0],[np.sqrt(3)/2,-0.5,0]])):
    """
    Returns an array that tells you where to put the source copy in a triangular
    lattice supercell band diagram calculation given the location of the
    original source. Supports sqrt(3)-by-1 (x-by-y) supercells with the default
    basis vector order given here.
    """
    if loc[0] < 0 and loc[1] < 0:
        return a_basis[0,:]
    elif loc[0] < 0 and loc[1] > 0:
        return a_basis[1,:]
    elif loc[0] > 0 and loc[1] > 0:
        return -1*a_basis[0,:]
    elif loc[0] > 0 and loc[1] < 0:
        return -1*a_basis[1,:]
    else:
        raise RuntimeError("Try running again; one of the sources was placed in \
                           an inconvenient position.")
        
def filter_freq_lists(freqs1, freqs2, amps1, amps2, phi_correct, tol = 10**-2):
    """
    Takes two frequency lists from different harminv points and uses them to
    filter for the correct frequencies at a single k_point
    """
    if phi_correct < 0:
        phi_cor = phi_correct + 2*np.pi
    else:
        phi_cor = phi_correct
        
    freqs = []
    i_f = 0
    for f in freqs1:
        searching = True
        keep_f = True
        i = 0
        # First, check for a match
        while searching:
            if np.abs((f.real-freqs2[i].real)/f.real) < tol:
                searching = False
                i_match = i
            i += 1
            if i == len(freqs2):
                searching = False
                keep_f = False
                print("Did not find matching frequency for freq:")
                print(f)
                
        # Next, make sure the phase shift is correct
        if keep_f:
            phi = np.angle(amps2[i_match]/amps1[i_f])
            if phi < 0:
                phi = phi + 2*np.pi
            if np.abs(phi_cor-phi) < tol*2*np.pi:
                freqs.append(f)
            else:
                print("Phase shift did not match for freq:")
                print(f)

        i_f += 1
    
    return freqs

class Logger(object):
    """
    Object for piping stdout outputs to both file and terminal
    """
    def __init__(self, savefile):
        sys.stdout.flush()
        self.terminal = sys.stdout
        self.log = open(savefile, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def close(self):
        self.log.close()

    def flush(self):
        self.terminal.flush()

###############################################################################
## Plasma-focused meep class 
###############################################################################
class Plasmeep:
    def __init__(self, a, res, dpml, nx, ny, nz = 0, B = np.array([0,0,0]),\
                 units = "T", verbose = True):
        """
        Initialize plasmeep object

        Args:
            a: dimensionalized unit length (e.g. 0.01 m) for translating to and
               from non-dimensionalized units
            res: pixels per a unit
            dpml: thickness of PML boundary in a units
            nx: number of a units in x-direction for active sim region
            ny: number of a units in y-direction for active sim region
            nz: number of a units in z-direction for active sim region
            B: externally applied b-field, vector with cyclotron frequency in a
               units as magnitude
            verbose: bool
        """
        self.a = a
        self.res = res
        self.cell = mp.Vector3(nx, ny, nz)
        self.pml = [mp.PML(thickness=dpml)]
        if units == "T":
            B = e*B/me/2/np.pi/c*a
        elif units == "Hz":
            B = B*a/c
        elif units == "rad/s":
            B = B/2/np.pi/c*a
        elif units == "a":
            B = B
        else:
            raise RuntimeError("Unrecognized units.")
        self.b = mp.Vector3(B[0],B[1],B[2])
        if B[0] != 0 or B[1] != 0 or B[2] != 0:
            self.is_magnetized = True
        else:
            self.is_magnetized = False
        self.geometry = [] # Objects in simulation domain
        self.sources = []
        self.verbose = verbose

    ###########################################################################
    ## Simulation Domain Stuff
    ###########################################################################
    def Get_Med(self, eps, wp = 0, gamma = 0, PEC = False, PMC = False):
        """
        Return a meep medium object with the correct susceptibility.

        Args:
            eps: 0th order relative permittivity
            wp: plasma frequency in a units
            gamma: collision frequency in a units
            PEC/PMC: bools, is perfect electric/magnetic conductor?
        """
        if PEC:
            return mp.perfect_electric_conductor
        if PMC:
            return mp.perfect_magnetic_conductor

        if wp > 0:
            if self.is_magnetized:
                gyro_susc = [mp.GyrotropicDrudeSusceptibility(frequency=wp,\
                             gamma=gamma, sigma=1, bias=self.b)]
                return mp.Medium(epsilon=eps, mu=1,\
                                 E_susceptibilities=gyro_susc)
            else:
                drude = [mp.DrudeSusceptibility(frequency=wp, gamma=gamma,\
                         sigma=1)]
                return mp.Medium(epsilon=eps, mu=1,\
                                 E_susceptibilities=drude)
        else:
            return mp.Medium(epsilon=eps, mu=1)
    
    
    def Add_Rod(self, r, center, eps = 1, wp = 0, gamma = 0,\
                axis = [0,0,1], height = mp.inf):
        """
        Add a single rod with radius r at center to geometry.

        Args:
            r: radius of the rod in a units
            center: x,y,z coords of the rod center in a units
            eps: 0th order relative permittivity
            wp: plasma frequency in a units
            gamma: collision frequency in a units
            axis: 1-hot vector selecting cylinder axis
            height: cylinder height, for 3D designs.
        """
        medium = self.Get_Med(eps, wp, gamma)
        self.geometry.append(mp.Cylinder(r, material=medium,\
            center=mp.Vector3(center[0],center[1],center[2]),\
            axis=mp.Vector3(axis[0],axis[1],axis[2]), height=height))


    def density_to_meep_wp(self, n):
        """
        Convert electron density to a MEEP Drude frequency (a units).

        n: electron density in m^{-3}.
        WP(n) returns the angular plasma frequency omega_p in rad/s.
        Get_Med / DrudeSusceptibility.frequency is a cyclic frequency
        in units of c/a, same convention as Nondimensionalize_Freq (Hz -> a).

        Therefore:
            f_p [Hz]  = WP(n) / (2*pi)
            wp_meep   = Nondimensionalize_Freq(f_p) = f_p * a / c

        Do not pass WP(n) directly to Nondimensionalize_Freq; that skips 2*pi
        and is inconsistent with Bands_MagPlasma (fp in Hz) and Get_Med.
        """
        n_arr = np.asarray(n, dtype=float)
        if np.any(~np.isfinite(n_arr)) or np.any(n_arr < 0.0):
            raise ValueError("Electron density n must be finite and n >= 0.")
        f_hz = WP(n_arr) / (2.0 * np.pi)
        wp = self.Nondimensionalize_Freq(f_hz)
        if np.isscalar(n) or getattr(n, "ndim", 1) == 0:
            return float(np.asarray(wp).reshape(()))
        return wp


    def Add_Rod_radial_shells(self, R_p, center, n_e_profile, n_shells = 12,\
                              eps = 1, gamma = 0, axis = [0,0,1],\
                              height = mp.inf):
        """
        Add a plasma cylinder of radius R_p as concentric Drude shells.

        This is a reusable geometry helper: n_e_profile may be any callable
        n_e(r) -> density in m^{-3}. Bessel_profile is one valid choice.

        Args:
            R_p: plasma radius in a units (same length convention as Add_Rod)
            center: x,y,z coords of the cylinder center in a units
            n_e_profile: callable. Argument r is the radial coordinate in a
                units (same as R_p). Return value is electron density in
                m^{-3}. Sampled at each shell midpoint.
            n_shells: number of equal-width rings covering 0 <= r <= R_p.
                Default 12, in the guide's ~10-15 recommendation.
            eps: 0th-order relative permittivity of each shell (usually 1)
            gamma: collision frequency in a units, same as Add_Rod / Get_Med.
                Convert from Hz with Nondimensionalize_Freq at the call site,
                as Bands_MagPlasma does. Do not pass rad/s.
            axis: 1-hot vector selecting cylinder axis
            height: cylinder height, for 3D designs.

        Geometry / precedence:
            Equal-width shells i = 0..n_shells-1 occupy
            [i, i+1] * (R_p / n_shells), sampled at the midpoint.
            Cylinders are appended largest-first, matching Add_Bulb's J0
            shells. MEEP later objects overwrite earlier ones, so the
            last (smallest) cylinder is the core.

        Density -> MEEP frequency (does not change Get_Med):
            wp = density_to_meep_wp(n) = Nondimensionalize_Freq(WP(n)/(2*pi))

        Returns:
            list of dicts, inner shell first: r_inner, r_outer, r_mid, n_e, wp
        """
        try:
            R_p = float(R_p)
        except (TypeError, ValueError) as exc:
            raise ValueError("R_p must be a real scalar.") from exc
        if isinstance(n_shells, bool) or not isinstance(n_shells, (int, np.integer)):
            raise ValueError("n_shells must be an integer >= 1.")
        n_shells = int(n_shells)
        if not np.isfinite(R_p) or R_p <= 0.0:
            raise ValueError("R_p must be a finite plasma radius, R_p > 0.")
        if n_shells < 1:
            raise ValueError("n_shells must be an integer >= 1.")
        if not callable(n_e_profile):
            raise ValueError("n_e_profile must be a callable n_e(r) in m^{-3}.")

        dr = R_p / n_shells
        shells = []
        for i in range(n_shells):
            r_inner = i * dr
            r_outer = (i + 1) * dr
            r_mid = 0.5 * (r_inner + r_outer)
            n_mid = n_e_profile(r_mid)
            try:
                n_mid = float(np.asarray(n_mid, dtype=float).reshape(()))
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "n_e_profile(r) must return a scalar electron density "
                    "in m^{-3}."
                ) from exc
            if not np.isfinite(n_mid) or n_mid < 0.0:
                raise ValueError(
                    "n_e_profile(r) must return a finite density n >= 0; "
                    "got n_e({0}) = {1}.".format(r_mid, n_mid)
                )
            wp = self.density_to_meep_wp(n_mid)
            shells.append({
                "r_inner": r_inner,
                "r_outer": r_outer,
                "r_mid": r_mid,
                "n_e": n_mid,
                "wp": wp,
            })

        # Largest cylinder first so later, smaller cylinders own the core.
        # Same overlap convention as Add_Bulb profile == 6 / 'J0'.
        for shell in reversed(shells):
            self.Add_Rod(shell["r_outer"], center, eps, shell["wp"], gamma,\
                         axis, height)
        return shells


    def Add_Horn(self, center, theta, feed_width, feed_length, eps,\
                 flare_length = 0, aperture_width = None, height = mp.inf):
        """
        Add a 2D dielectric waveguide feed with an optional flare.

        This is a reusable geometry helper. It places one feed at an
        arbitrary location and orientation. Five 72-degree tomography
        horns are the caller's responsibility, not this function's.

        Material: dielectric epsilon (tomography feeds use eps ~ 4).
        eps is required. There is no default: PlasMEEP's usual eps=1
        would be invisible against vacuum, and PMM-Design's eps=-1000
        PEC walls are the wrong object for this project.

        All lengths are in a units.

        Args:
            center: x,y[,z] of the straight-feed centroid in a units.
                The recommended future port plane passes through this
                point. z defaults to 0.
            theta: feed-axis angle in radians, measured from +x toward
                +y. The unit axis u = (cos theta, sin theta) points from
                the back of the feed toward the aperture (launch
                direction, toward the device when the horn faces inward).
            feed_width: full width of the uniform waveguide, a units.
            feed_length: length of the uniform waveguide along u, a units.
            eps: relative permittivity of the feed and flare. Must be
                finite and > 0. Do not pass PMM-Design's -1000.
            flare_length: length of the optional flare along u. 0 means
                no flare (straight feed only).
            aperture_width: full width at the flare mouth. Required if
                flare_length > 0, and must exceed feed_width so the
                aperture expands toward the mouth.
            height: extrusion along z (default mp.inf; 2D TM).

        Geometry:
            Feed occupies center +/- (feed_length/2) * u and is a
            rectangle of width feed_width. If flare_length > 0, a
            trapezoid is attached on the aperture side (the +u end).
            Representation: mp.Prism polygons (MEEP-native, robust at
            arbitrary theta). Feed first, then flare.

        Returns a metadata dict for a future Add_Port:
            center, theta, axis, normal, feed_width, feed_length,
            flare_length, aperture_width, eps, feed_back, feed_front,
            aperture_center, port_center, port_p1, port_p2.
            port_* is the mid-feed cross-section, inside the uniform
            dielectric, perpendicular to the axis.
        """
        center = np.asarray(center, dtype=float).reshape(-1)
        if center.size == 2:
            center = np.array([center[0], center[1], 0.0], dtype=float)
        elif center.size == 3:
            center = center.astype(float)
        else:
            raise ValueError("center must be an (x, y) or (x, y, z) position.")
        if not np.all(np.isfinite(center)):
            raise ValueError("center must be finite.")

        try:
            theta = float(theta)
            feed_width = float(feed_width)
            feed_length = float(feed_length)
            flare_length = float(flare_length)
            eps = float(eps)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "theta, feed_width, feed_length, flare_length, and eps "
                "must be real scalars."
            ) from exc
        if not np.isfinite(theta):
            raise ValueError("theta must be a finite angle in radians.")
        if not np.isfinite(feed_width) or feed_width <= 0.0:
            raise ValueError("feed_width must be finite and > 0.")
        if not np.isfinite(feed_length) or feed_length <= 0.0:
            raise ValueError("feed_length must be finite and > 0.")
        if not np.isfinite(flare_length) or flare_length < 0.0:
            raise ValueError("flare_length must be finite and >= 0.")
        if not np.isfinite(eps) or eps <= 0.0:
            raise ValueError(
                "eps must be a finite dielectric constant > 0. "
                "Require it explicitly; do not use PEC-like negative eps."
            )

        if flare_length == 0.0:
            if aperture_width is not None and float(aperture_width) != feed_width:
                raise ValueError(
                    "aperture_width must be omitted or equal feed_width "
                    "when flare_length is 0."
                )
            aperture_width = feed_width
        else:
            if aperture_width is None:
                raise ValueError(
                    "aperture_width is required when flare_length > 0."
                )
            try:
                aperture_width = float(aperture_width)
            except (TypeError, ValueError) as exc:
                raise ValueError("aperture_width must be a real scalar.")\
                    from exc
            if not np.isfinite(aperture_width) or aperture_width <= feed_width:
                raise ValueError(
                    "aperture_width must be finite and > feed_width so the "
                    "flare expands toward the mouth."
                )

        axis = np.array([np.cos(theta), np.sin(theta), 0.0], dtype=float)
        # Transverse unit vector: axis rotated +90 deg CCW in the xy plane.
        normal = np.array([-np.sin(theta), np.cos(theta), 0.0], dtype=float)

        feed_back = center - 0.5 * feed_length * axis
        feed_front = center + 0.5 * feed_length * axis
        aperture_center = feed_front + flare_length * axis

        def _quad(along_a, half_a, along_b, half_b):
            # CCW when viewed from +z: start at -normal on the back edge.
            return np.vstack([
                along_a - half_a * normal,
                along_b - half_b * normal,
                along_b + half_b * normal,
                along_a + half_a * normal,
            ])

        feed_verts = _quad(feed_back, 0.5 * feed_width,
                           feed_front, 0.5 * feed_width)
        self.Add_Prism(feed_verts, axis=[0, 0, 1], height=height, eps=eps)

        if flare_length > 0.0:
            flare_verts = _quad(feed_front, 0.5 * feed_width,
                                aperture_center, 0.5 * aperture_width)
            self.Add_Prism(flare_verts, axis=[0, 0, 1], height=height, eps=eps)
        else:
            flare_verts = None

        port_p1 = center - 0.5 * feed_width * normal
        port_p2 = center + 0.5 * feed_width * normal
        return {
            "center": center.copy(),
            "theta": theta,
            "axis": axis,
            "normal": normal,
            "feed_width": feed_width,
            "feed_length": feed_length,
            "flare_length": flare_length,
            "aperture_width": aperture_width,
            "eps": eps,
            "feed_back": feed_back,
            "feed_front": feed_front,
            "aperture_center": aperture_center,
            "port_center": center.copy(),
            "port_p1": port_p1,
            "port_p2": port_p2,
            "feed_vertices": feed_verts,
            "flare_vertices": flare_verts,
        }


    def Add_Bulb(self, r_bulb, center, wp = 0, gamma = 0,\
                 axis = [0,0,1], height = mp.inf, profile = 0):
        """
        Add a single plasma bulb to the geometry list

        Args:
            r_bulb: (inner, outer) radius of the rod in a units
            center: x,y,z coords of the rod center in a units
            eps: 0th order relative permittivity
            wp: plasma frequency in a units
            gamma: collision frequency in a units
            axis: 1-hot vector selecting cylinder axis
            height: cylinder height, for 3D designs.
            profile: int or string selecting density profile
        """
        quartz = self.Get_Med(3.8)
        if profile == 0:
            plasma = self.Get_Med(1, wp, gamma)
            vac = self.Get_Med(1)
            self.geometry.append(mp.Cylinder(r_bulb[1], material=quartz,\
                center=mp.Vector3(center[0],center[1],center[2]),\
                axis=mp.Vector3(axis[0],axis[1],axis[2]), height=height))
            self.geometry.append(mp.Cylinder(r_bulb[0], material=vac,\
                center=mp.Vector3(center[0],center[1],center[2]),\
                axis=mp.Vector3(axis[0],axis[1],axis[2]), height=height))
            self.geometry.append(mp.Cylinder(4.6*r_bulb[0]/6.5, material=plasma,\
                center=mp.Vector3(center[0],center[1],center[2]),\
                axis=mp.Vector3(axis[0],axis[1],axis[2]), height=height))

        elif profile == 6:
            self.geometry.append(mp.Cylinder(r_bulb[1], material=quartz,\
                center=mp.Vector3(center[0],center[1],center[2]),\
                axis=mp.Vector3(axis[0],axis[1],axis[2]), height=height))
            shells = int(self.res*r_bulb[0])
            for i in range(shells):
                r = r_bulb[0]*(1-i/shells)
                fp = wp*(4*(4.6/6.5)**2*(1-(r/r_bulb[0])**6)/3)**0.5
                shell = self.Get_Med(1, fp, gamma)
                self.geometry.append(mp.Cylinder(r, material=shell,\
                    center=mp.Vector3(center[0],center[1],center[2]),\
                    axis=mp.Vector3(axis[0],axis[1],axis[2]), height=height))

        elif profile == 'J0':
            self.geometry.append(mp.Cylinder(r_bulb[1], material=quartz,\
                center=mp.Vector3(center[0],center[1],center[2]),\
                axis=mp.Vector3(axis[0],axis[1],axis[2]), height=height))
            shells = int(self.res*r_bulb[0])
            for i in range(shells):
                r = r_bulb[0]*(1-i/shells)
                fp = wp*((J01**2/2/1.24846)*sp.jv(0,J01*r/r_bulb[0]))**0.5
                shell = self.Get_Med(1, fp, gamma)
                self.geometry.append(mp.Cylinder(r, material=shell,\
                    center=mp.Vector3(center[0],center[1],center[2]),\
                    axis=mp.Vector3(axis[0],axis[1],axis[2]), height=height))

        else:
            raise RuntimeError("Density profile not implemented yet.")


    def Add_Block(self, low_left, size, eps = 1, wp = 0, gamma = 0, PEC = False,\
                  PMC = False):
        """
        Add a single rod with radius r and rel. permittivity eps to epsr.

        Args:
            center: x,y,z coords of the bottom left corner in a units
            size: z, y, and z size of the block in a units
            eps: relative permittivity of the block
            wp: plasma frequency in a units
            gamma: collision frequency in a units
        """
        medium = self.Get_Med(eps, wp, gamma, PEC, PMC)
        center = (low_left[0]+size[0]/2, low_left[1]+size[1]/2, low_left[2]+size[2]/2)
        self.geometry.append(mp.Block(mp.Vector3(size[0],size[1],size[2]),\
                center=mp.Vector3(center[0],center[1],center[2]),\
                material=medium))


    def Add_Prism(self, vertices, axis = [0,0,1], height = mp.inf, eps = 1,\
                  wp = 0, gamma = 0, PEC = False, PMC = False):
        """
        Add a single rod with radius r and rel. permittivity eps to the static
        elems array.

        Args:
            vertices: x,y,z coords of the vertices of the prism in a units
            height: height of prism in axis direction
            axis: one-hot vector, direction of height dimension
            eps: relative permittivity of the block
            wp: plasma frequency in a units
            gamma: collision frequency in a units
        """
        medium = self.Get_Med(eps, wp, gamma, PEC, PMC)
        Vertices = []
        for i in range(vertices.shape[0]):
            Vertices.append(mp.Vector3(vertices[i,0], vertices[i,1],\
                                       vertices[i,2]))
        Axis = mp.Vector3(axis[0],axis[1],axis[2])
        self.geometry.append(mp.Prism(vertices=Vertices, height=height,\
                axis=Axis,material=medium))


    def Rod_Array(self, r, xyz_start, a_l, rod_eps, axis, height = mp.inf):
        """
        Add a 2D cubic rod array to the domain.

        Args:
            r: radius of rods in a units
            xyz_start: coords of the bottom left of the array in a units
            axis: axis of the rods e.g. [0,0,1] for z
            rod_eps: np array of size (nrods_x, nrods_y, 3) giving the 
                     permittivity parameters of each of the rods. e.g. 
                     rod_eps[i,j,:] = [eps, fp, gamma]
            a_l: lattice spacing
            height: height of the rods in a units
        """
        for k in range(3):
            if axis[k] == 1 and k == 0:
                pass
            elif axis[k] == 1 and k == 1:
                pass
            elif  axis[k] == 1 and k == 2:
                for i in range(rod_eps.shape[0]):
                    for j in range(rod_eps.shape[1]):
                        x = xyz_start[0] + i*a_l
                        y = xyz_start[1] + j*a_l
                        self.Add_Rod(r, [x,y,xyz_start[2]], rod_eps[i,j,0],\
                                rod_eps[i,j,1], rod_eps[i,j,2], axis, height)
            else:
                pass


    def Rod_Array_Hexagon(self, xyz_cen, side_dim, r = 0.5, a = 1, eps = 1,\
                          a_basis = np.array([[0,1,0],[np.sqrt(3)/2,1./2,0]]),\
                          bulbs = False, r_bulb = (0, 0), wp = 0, gamma = 0,\
                          axis = np.array([0,0,1]), height = mp.inf,\
                          profile = 0):
        """
        Add a hexagonal triangular array of rods or bulbs to the domain

        Args:
            xyz_cen: np.array, center of hexagon
            side_dim: number of rods along one side of hexagon
            r: radius of the rod in a units
            a: array spacing in a units
            eps: 0th order relative permittivity of rods
            a_basis: np.array, contains basis vectors that determine orientation
                     of hexagon
            bulbs: bool, true if using bulbs
            r_bulb: (inner, outer) radius of the rod in a units
            wp: plasma frequency in a units
            gamma: collision frequency in a units
            axis: 1-hot vector selecting cylinder axis
            height: cylinder height, for 3D designs
            profile: int or string selecting density profile
        """
        for i in range(side_dim):
            for j in range(side_dim+i):
                b1_loc = xyz_cen-(side_dim-1-i)*a*a_basis[0,:]-(i-j)*a*a_basis[1,:]
                b2_loc = xyz_cen+(side_dim-1-i)*a*a_basis[0,:]+(i-j)*a*a_basis[1,:]
                if bulbs:
                    self.Add_Bulb(r_bulb, b1_loc, wp = wp, gamma = gamma,\
                                          axis = axis, height = height,\
                                          profile = profile)
                    if i < side_dim - 1:
                        self.Add_Bulb(r_bulb, b2_loc, wp = wp, gamma = gamma,\
                                          axis = axis, height = height,\
                                          profile = profile)
                else:
                    self.Add_Rod(r, b1_loc, eps = eps, wp = wp, gamma = gamma,\
                                 axis = axis, height = height)
                    if i < side_dim - 1:
                        self.Add_Rod(r, b2_loc, eps = eps, wp = wp, gamma = gamma,\
                                     axis = axis, height = height)
            
        return


    ###########################################################################
    ## Simulation
    ###########################################################################
    def Add_Cont_Source(self, fs, Pol, center, size, E = True):
        """
        Add a source to the domain.

        Args:
            fs: Source frequency
            Pol: np array, unit vector pointing in the polarization direction,
                 also controls amplitude
            center: location of the center of the source
            size: size of source (can choose (0,0,0) for point-dipole)
            E: bool determining whether source is E or H.
        """
        for i in range(3):
            if np.abs(Pol[i]) > 0:
                if E:
                    if i == 0:
                        comp = mp.Ex
                    elif i == 1:
                        comp = mp.Ey
                    else:
                        comp = mp.Ez
                else:
                    if i == 0:
                        comp = mp.Hx
                    elif i == 1:
                        comp = mp.Hy
                    else:
                        comp = mp.Hz

                self.sources.append(mp.Source(\
                mp.ContinuousSource(frequency=fs, is_integrated=True),\
                component=comp,\
                center=mp.Vector3(center[0],center[1],center[2]),\
                size=mp.Vector3(size[0],size[1],size[2]), amplitude=Pol[i]))


    def Add_Gauss_Source(self, fs, df, Pol, center, size, E = True):
        """
        Add a source to the domain.

        Args:
            fs: Source frequency
            Pol: np array, unit vector pointing in the polarization direction,
                 also controls amplitude
            center: location of the center of the source
            size: size of source (can choose (0,0,0) for point-dipole)
        """
        for i in range(3):
            if np.abs(Pol[i]) > 0:
                if E:
                    if i == 0:
                        comp = mp.Ex
                    elif i == 1:
                        comp = mp.Ey
                    else:
                        comp = mp.Ez
                else:
                    if i == 0:
                        comp = mp.Hx
                    elif i == 1:
                        comp = mp.Hy
                    else:
                        comp = mp.Hz

                self.sources.append(mp.Source(\
                mp.GaussianSource(frequency=fs, fwidth=df),\
                component=comp,\
                center=mp.Vector3(center[0],center[1],center[2]),\
                size=mp.Vector3(size[0],size[1],size[2]), amplitude=Pol[i]))

                
    def Wipe_sources(self):
        """
        Deletes source list
        """
        self.sources = []
                

    def Get_Sim(self):
        return mp.Simulation(cell_size=self.cell, boundary_layers=self.pml,\
                geometry=self.geometry, sources=self.sources,\
                resolution=self.res, default_material=mp.Medium(epsilon=1,mu=1))


    def Bands2D(self, k_points, fs=1, Pol=np.array([0,0,1]),\
                size = mp.Vector3(0,0,0), E = True,\
                freq_range = np.array([0,1.5]), savefile = 'bands.csv',\
                lattice = 'square',\
                a_basis = np.array([[1,0,0],[0,1,0],[0,0,1]]), plot = True,\
                MPB = False, num_bands = 10, TETM = 'both', tol = 10**-2):
        """
        Create a 2-D band diagram given the domain as a single cell of a square
        or triangular lattice photonic crystal

        Args:
            k_points: number of points to interpolate between Gamma-X, X-M, and
                      M-Gamma. (k_points*3+4 total points will be simulated)
            fs: Source frequency, best set to 1, but could set to a small value
                if one is interested in resolving bands close to f=0, for ex.
            Pol: np array, unit vector pointing in the polarization direction,
                 also controls amplitude
            size: size of source (can choose (0,0,0) for point-dipole)
            E: bool determining whether source is E or H.
            freq_range: array, [min, max] of frequency range for band diagram
            savefile: Filename for bands data, must be .csv
        """
        if MPB:
            return self.Bands2D_MPB(k_points, savefile, lattice, plot,\
                                    num_bands, TETM)
        if lattice == 'triangular':
            return self.Bands2DTri(k_points, a_basis, fs, Pol, size, E,\
                                   freq_range, savefile, plot, tol = tol)
        if lattice != 'square':
            raise RuntimeError("Band diagrams with lattices that are not \
                                square or triangular must be produced using MPB \
                                (and therefore cannot have dispersive elements).")

        f_scales = np.array([[0.5,1.5], [0.2,0.2], [0.5,0.75], [1.5,3]])
        for i in range(4):
            r = np.clip(0.5*np.random.randn(3), -0.45, 0.45)
            loc = [self.cell[0]*r[0], self.cell[1]*r[1], self.cell[2]*r[2]]
            self.Add_Gauss_Source(fs*f_scales[i,0], fs*f_scales[i,1], Pol, loc,\
                                  size, E)

        if self.cell[2] != 0:
            raise RuntimeError('Can only generate band diagrams for 2D X-Y \
                                systems right now')

        sim = self.Get_Sim()
        if os.path.exists(savefile):
            os.remove(savefile)

        freqs = sim.run_k_points(300, mp.interpolate(k_points,\
                         [mp.Vector3(0,0,0), mp.Vector3(0.5,0,0),\
                         mp.Vector3(0.5,0.5,0), mp.Vector3(0,0,0)]))

        for i in range(len(freqs)):
            for j in range(len(freqs[i])):
                freqs[i][j] = freqs[i][j].real
                
        with open(savefile,'w') as sfile:
            writer = csv.writer(sfile)
            for i in range(len(freqs)):
                writer.writerow(freqs[i])

        if plot:
            self.BandPlot(k_points, freq_range, savefile)

        return

    
    def Bands2DTri(self, k_points,\
                a_basis = np.array([[np.sqrt(3)/2,0.5,0],[np.sqrt(3)/2,-0.5,0]]),\
                fs=1, Pol=np.array([0,0,1]), size = mp.Vector3(0,0,0), E = True,\
                freq_range = np.array([0,1.5]), savefile = 'bands.csv',\
                plot = True, tol = 10**-2):
        """
        Create a 2-D band diagram given the domain as a single cell of a
        triangular lattice photonic crystal

        Args:
            k_points: number of points to interpolate between Gamma-X, X-M, and
                      M-Gamma. (k_points*3+4 total points will be simulated)
            fs: Source frequency, best set to 1, but could set to a small value
                if one is interested in resolving bands close to f=0, for ex.
            Pol: np array, unit vector pointing in the polarization direction,
                 also controls amplitude
            size: size of source (can choose (0,0,0) for point-dipole)
            E: bool determining whether source is E or H.
            freq_range: array, [min, max] of frequency range for band diagram
            savefile: Filename for bands data, must be .csv
        """
        if a_basis[0,0] != np.sqrt(3)/2:
            raise RuntimeError('For now, the triangular lattice supercell must \
                                be sqrt(3)-by-1 in size (thats x-by-y), and the \
                                basis vecs must be in the order [/,\]')
        
        f_scales = np.array([[0.5,1.2], [0.2,0.2], [0.5,0.75], [1.5,3]])
        src_params = []
        shifts = []
        for i in range(1):
            r = np.clip(0.5*np.random.randn(3), -0.45, 0.45)
            loc = [self.cell[0]*r[0], self.cell[1]*r[1], self.cell[2]*r[2]]
            shifts.append(get_trilat_shift(loc, a_basis))
            self.Add_Gauss_Source(fs*f_scales[i,0], fs*f_scales[i,1], Pol, loc,\
                                  size, E)
            src_params.append([fs*f_scales[i,0],fs*f_scales[i,1],Pol,loc,size,E])

        if self.cell[2] != 0:
            raise RuntimeError('Can only generate band diagrams for 2D X-Y \
                                systems right now')

        sim = self.Get_Sim()
        if os.path.exists(savefile):
            os.remove(savefile)

        k_points_list = mp.interpolate(k_points, [mp.Vector3(0,0,0),\
                                       mp.Vector3(1./np.sqrt(3),1./3,0),\
                                       mp.Vector3(1./np.sqrt(3),0,0),\
                                       mp.Vector3(0,0,0)])
        freqs = self.supercell_k_points(sim, 300, k_points_list, src_params,\
                                   shifts, tol = tol)

        for i in range(len(freqs)):
            for j in range(len(freqs[i])):
                freqs[i][j] = freqs[i][j].real
                
        with open(savefile,'w') as sfile:
            writer = csv.writer(sfile)
            for i in range(len(freqs)):
                writer.writerow(freqs[i])

        if plot:
            self.BandPlot(k_points, freq_range, savefile, lattice = 'triangular')

        return
    

    def Bands2D_MPB(self, k_points, savefile = 'bands.csv', lattice = 'square',\
                    plot = True, num_bands = 10, TETM = 'both', dims = False,\
                    dimlabel = 'GHz'):
        """
        Create a 2-D band diagram given the domain as a single cell of a specified
        lattice photonic crystal using MPB (meep photonic bands)

        Args:
            k_points: number of points to interpolate between Gamma-X, X-M, and
                      M-Gamma. (k_points*3+4 total points will be simulated)
            fs: Source frequency, best set to 1, but could set to a small value
                if one is interested in resolving bands close to f=0, for ex.
            Pol: np array, unit vector pointing in the polarization direction,
                 also controls amplitude
            size: size of source (can choose (0,0,0) for point-dipole)
            E: bool determining whether source is E or H.
            freq_range: array, [min, max] of frequency range for band diagram
            savefile: Filename for bands data, must be .csv
            lattice: str, type of lattice, e.g. 'square' or 'triangular'
            plot: bool, determines whether or not the bands are plotted.
            num_bands: number of bands you would like to solve for.
            TETM: str, determines whether just TE, TM bands are plotted or both.
        """
        if lattice == 'square':
            geometry_lattice = mp.Lattice(size=mp.Vector3(1,1,0),\
                                basis1=mp.Vector3(1,0,0),\
                                basis2=mp.Vector3(0,1,0))
            k_pts = [mp.Vector3(), mp.Vector3(0.5,0,0), mp.Vector3(0.5,0.5,0),\
                     mp.Vector3()]
            k_pts = mp.interpolate(k_points, k_pts)
            
        if lattice == 'triangular':
            geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1),\
                       basis1=mp.Vector3(np.sqrt(3)/2, 0.5),\
                       basis2=mp.Vector3(np.sqrt(3)/2, -0.5))
            k_pts = [mp.Vector3(), mp.Vector3(-1./3, 1./3, 0),\
                     mp.Vector3(-0.5, 0, 0),\
                     mp.Vector3()]
            k_pts = mp.interpolate(k_points, k_pts)
            
        modes = mpb.ModeSolver(geometry=self.geometry,\
                geometry_lattice=geometry_lattice, k_points=k_pts,\
                resolution=self.res, num_bands=num_bands)

        if TETM == 'TE':
            modes.run_te()
            te_freqs = modes.all_freqs
            te_gaps = modes.gap_list
            return self.BandPlot_MPB([te_freqs, te_freqs], [te_gaps, te_gaps],\
                                     savefile, TETM, dims, dimlabel,\
                                     lattice = lattice)

        elif TETM == 'TM':
            modes.run_tm()
            tm_freqs = modes.all_freqs
            tm_gaps = modes.gap_list
            return self.BandPlot_MPB([tm_freqs, tm_freqs], [tm_gaps, tm_gaps],\
                                     savefile, TETM, dims, dimlabel,\
                                     lattice = lattice)

        elif TETM == 'both':
            modes.run_te()
            te_freqs = modes.all_freqs
            te_gaps = modes.gap_list
            modes.run_tm()
            tm_freqs = modes.all_freqs
            tm_gaps = modes.gap_list
            return self.BandPlot_MPB([tm_freqs, te_freqs], [tm_freqs, te_gaps],\
                                     savefile, TETM, dims, dimlabel,\
                                     lattice = lattice)
        
        
    def supercell_k_points(self, sim, T, k_points_list, main_src_params,\
                       lattice_shifts, tol = 10**-2):
        """
        Alternate version of the run_k_points function that works for supercell
        simulations, needed for band diagrams for PCs with non-rectangular unit
        cells
        
        Args:
            sim: meep sim object for lattice supercell
            T: time for which to run sources (300 is a good choice)
            k_points_list: list of mp.Vector3's for each k-point you want to
                           run, use mp.interpolate()
            main_src_params: list of lists countaining arguments for original
                             source(s)
            lattice_shifts: list of lists containing the lattice vector shifts
                            needed to create the appropriate copies of the
                            original source(s)
        """
        all_freqs = []
        k_index = 0
        for k in k_points_list:
            k_index += 1
            sim.restart_fields()
            self.Wipe_sources()
            for i in range(len(main_src_params)):
                fs = main_src_params[i][0]
                df = main_src_params[i][1]
                Pol = main_src_params[i][2]
                loc = main_src_params[i][3]
                size = main_src_params[i][4]
                E = main_src_params[i][5]
                self.Add_Gauss_Source(fs, df, Pol, loc, size, E)
                x_l = lattice_shifts[i]
                x = mp.Vector3(x_l[0], x_l[1], x_l[2])
                phi = 2*np.pi*x.dot(k)
                loc_shift = [loc[0]+x_l[0],loc[1]+x_l[1],loc[2]+x_l[2]]
                self.Add_Gauss_Source(fs, df, Pol*np.exp(1j*phi),\
                                      loc_shift, size, E)
            sim.change_sources(self.sources)
            h = self.supercell_run_k_point(sim, T, k, lattice_shifts)
            freqs1 = [complex(m.freq, m.decay) for m in h[0].modes]
            amps1 = [m.amp for m in h[0].modes]
            freqs2 = [complex(m.freq, m.decay) for m in h[1].modes]
            amps2 = [m.amp for m in h[1].modes]
            
            x_l = lattice_shifts[0]
            x = mp.Vector3(x_l[0], x_l[1], x_l[2])
            phi = 2*np.pi*x.dot(k)
            freqs = filter_freq_lists(freqs1, freqs2, amps1, amps2,\
                                      phi_correct = phi, tol = tol)

            if self.verbose:
                print("freqs:, {}, {}, {}, {}, "\
                      .format(k_index, k.x, k.y, k.z), end='')
                print(', '.join([str(f.real) for f in freqs1]))
                print("freqs-im:, {}, {}, {}, {}, "\
                      .format(k_index, k.x, k.y, k.z), end='')
                print(', '.join([str(f.imag) for f in freqs1]))
                print("freqs_shift:, {}, {}, {}, {}, "\
                      .format(k_index, k.x, k.y, k.z), end='')
                print(', '.join([str(f.real) for f in freqs2]))
                print("freqs_shift-im:, {}, {}, {}, {}, "\
                      .format(k_index, k.x, k.y, k.z), end='')
                print(', '.join([str(f.imag) for f in freqs2]))

            all_freqs.append(freqs)
        
        return all_freqs

    
    def supercell_run_k_point(self, sim, t = None, k = None, lattice_shifts = []):
        """
        Lower level function called by `supercell_run_k_points` that runs a simulation
        for a single *k* point `k_point` and returns a `Harminv` instance for each of
        the periodic copies of the lead source.
        """
        components = [s.component for s in sim.sources]
        pts = [s.center for s in sim.sources]

        src_freqs_min = min(
            s.src.frequency - 1 / s.src.width / 2
            if isinstance(s.src, mp.GaussianSource)
            else mp.inf
            for s in sim.sources
        )
        fmin = max(0, src_freqs_min)

        fmax = max(
            s.src.frequency + 1 / s.src.width / 2
            if isinstance(s.src, mp.GaussianSource)
            else 0
            for s in sim.sources
        )

        if not components or fmin > fmax:
            raise ValueError("Running with k_points requires a 'GaussianSource' \
                              source")

        sim.change_k_point(k)
        
        h = []
        sim.restart_fields()

        h.append(Harminv(components[0], pts[0], 0.5 * (fmin + fmax), fmax - fmin))
        sim.run(after_sources(h[0]), until_after_sources=t)
        
        sim.restart_fields()
 
        x_l = lattice_shifts[0]
        x = mp.Vector3(x_l[0], x_l[1], x_l[2])
        phase_shift = np.exp(1j*2*np.pi*x.dot(k))
        h.append(Harminv(components[0], pts[0]+x, 0.5 * (fmin + fmax), fmax - fmin))
        sim.run(after_sources(h[1]), until_after_sources=t)

        return h
    

    ############################################################################
    ## Visualization
    ############################################################################
    def BandPlot_MPB(self, bands, gaps, savefile, TETM = 'both', dims = False,\
                 dimlabel = 'GHz', lattice = 'square'):
        """
        Plots band diagram given output from MPB band solver

        Args:
            bands: list containing te, tm bands.
            gaps: list containing te, tm bandgaps.
            savefile: savefile, .csv format
            TETM: str, which polarizations you would like to plot
            dims: bool, dimensionalized frequency if true
            dimlabel: str, label for frequency units
        """
        fig, ax = plt.subplots()
        if dims:
            bands[0] = self.Dimensionalize_Freq(bands[0])/10**9
            bands[1] = self.Dimensionalize_Freq(bands[1])/10**9
            ax.set_ylabel('$f$ ('+dimlabel+')', fontsize = 16)
        else:
            ax.set_ylabel('$\omega/(2\pi c/a)$', fontsize = 16)
        x = range(len(bands[0]))
        # Plot bands
        for xz, tmz, tez in zip(x, bands[0], bands[1]):
            ax.scatter([xz]*len(tmz), tmz, color='blue')
            ax.scatter([xz]*len(tez), tez, color='red', facecolors='none')
        if TETM == 'both' or TETM == 'TM':
            ax.plot(bands[0], color='blue')
        if TETM == 'both' or TETM == 'TE':
            ax.plot(bands[1], color='red')
        ax.set_ylim([0, 1])
        ax.set_xlim([x[0], x[-1]])

        # Plot gaps
        if TETM == 'both' or TETM == 'TM':
            for gap in gaps[0]:
                if gap[0] > 1:
                    ax.fill_between(x, gap[1], gap[2], color='blue', alpha=0.2)
        if TETM == 'both' or TETM == 'TE':
            for gap in gaps[1]:
                if gap[0] > 1:
                    ax.fill_between(x, gap[1], gap[2], color='red', alpha=0.2)

        # Plot labels
        if TETM == 'both' or TETM == 'TM':
            ax.text(12, 0.04, 'TM bands', color='blue', size=15)
        if TETM == 'both' or TETM == 'TE':
            ax.text(13.05, 0.235, 'TE bands', color='red', size=15)

        points_in_between = (len(bands[0]) - 4) / 3
        tick_locs = [i*points_in_between+i for i in range(4)]
        if lattice == 'square':
            tick_labs = ['$\Gamma$', 'X', 'M', '$\Gamma$']
        elif lattice == 'triangular':
            tick_labs = ['$\Gamma$', 'K', 'M', '$\Gamma$']
        ax.set_xticks(tick_locs)
        ax.set_xticklabels(tick_labs, size=16)
        ax.grid(True)

        
        if dims:
            plt.savefig(savefile[:-4]+'_'+dimlabel+'.pdf', dpi=1500)
        else:
            plt.savefig(savefile[:-4]+'.pdf', dpi=1500)

        return fig, ax


    def BandPlot(self, k_points, freq_range, savefile = 'noname.csv',\
                 freq_list = [], dims = False, dimlabel = 'GHz',\
                 lattice = 'square'):
        """
        Plots band diagram given .csv output from Bands2D. See that function for
        arg descriptions.
        """
        # TODO: Find a more elegant way of checking if multiple band diagrams are
        # being plotted.
        if len(freq_list) == k_points*3+4:
            mult = False
            freqs = np.array([])
            indices = np.array([])
        else:
            if freq_list != []:
                mult = True
            else:
                mult = False
            freqs = []
            indices = []
        idx = 0
        no_bands = 0
        
        if freq_list != []:
            if mult:
                j = 0
                for f_list in freq_list:
                    idx = 0
                    freqs.append(np.array([]))
                    indices.append(np.array([]))
                    for i in range(len(f_list)):
                        num = len(f_list[i])
                        if num == 0:
                            no_bands += 1
                        else:
                            freqs[j] = np.append(freqs[j],\
                                                 np.array(f_list[i]).astype(float))
                            indices[j] = np.append(indices[j], idx*np.ones(num))
                        idx += 1
                    j += 1
            else:
                for i in range(len(freq_list)):
                    num = len(freq_list[i])
                    if num == 0:
                        no_bands += 1
                    else:
                        freqs = np.append(freqs,\
                                          np.array(freq_list[i]).astype(float))
                        indices = np.append(indices, idx*np.ones(num))
                    idx += 1        
        else:
            with open(savefile) as csvfile:
                readCSV = csv.reader(csvfile, delimiter=',')
                for row in readCSV:
                    num = len(row)
                    if row == [' ']:
                        no_bands += 1
                    else:
                        freqs = np.append(freqs, np.array(row).astype(float))
                        indices = np.append(indices, idx*np.ones(num))
                    idx += 1

        fig, ax = plt.subplots(1,1,figsize=(8,6))

        if lattice == 'square':
            xlabels = ['']*(3*(k_points+1)+1)
            xlabels[0] = '$\Gamma$'
            xlabels[k_points+1] = '$X$'
            xlabels[2*(k_points+1)] = '$M$'
            xlabels[3*(k_points+1)] = '$\Gamma$'
        elif lattice == 'triangular':
            xlabels = ['']*(3*(k_points+1)+1)
            xlabels[0] = '$\Gamma$'
            xlabels[k_points+1] = '$K$'
            xlabels[2*(k_points+1)] = '$M$'
            xlabels[3*(k_points+1)] = '$\Gamma$'

        if dims:                       
            ax.set_ylabel('$f$ ('+dimlabel+')', fontsize = 18)
        else:
            ax.set_ylabel('$\omega/(2\pi c/a)$', fontsize = 18)
        ax.tick_params(labelsize = 16)
        ax.set_xticks(range(3*(k_points+1)+1), minor=False)
        ax.set_xticklabels(xlabels, minor=False)
        if not mult:
            if dims:
                freqs = self.Dimensionalize_Freq(freqs)/10**9
            ax.scatter(indices,freqs)
        else:
            for k in range(len(freqs)):
                if dims:
                    freqs = self.Dimensionalize_Freq(freqs)/10**9
                ax.scatter(indices[k], freqs[k], color = colors[k])
                
        ax.set_xlim([-(1/50)*(3*(k_points+1)),(51/50)*(3*(k_points+1))])
        ax.set_ylim(freq_range)

        fig.tight_layout()
        if dims:
            plt.savefig(savefile[:-4]+'_'+dimlabel+'.pdf', dpi=1500)
        else:
            plt.savefig(savefile[:-4]+'.pdf', dpi=1500)

        return fig, ax

    
    def Viz_Domain(self, output_dir, freq = 1):
        Sim = self.Get_Sim()
        Sim.use_output_directory(output_dir)
        Sim.init_sim()
        mp.output_epsilon(Sim, frequency=freq)
        
        for filename in os.listdir(output_dir):
            if filename.endswith('.h5') and 'eps' in filename:
                f = h5py.File(output_dir+'/'+filename, 'r')
                eps = f['eps.r'][...]
                plt.imshow(np.flipud(eps.T), cmap='magma')
                plt.savefig(output_dir+'/domain.png')
                os.remove(output_dir+'/'+filename)
    
    def Sim_And_Plot(self, at_every, field, until, output_dir,\
                     component = 'hz', out_prefix = '', remove = True,\
                     all_frames = False, intensity = 1.0):
        Sim = self.Get_Sim()
        
        if field == 'E':
            meep_output = mp.output_efield
        elif field == 'H':
            meep_output = mp.output_hfield
        elif field == 'S':
            meep_output = mp.output_poynting
            raise RuntimeError("Haven't worked out plotting for \
                                Poynting vector yet.")

        data_dir = output_dir+out_prefix+'-out/'
        if os.path.isdir(data_dir):
            shutil.rmtree(data_dir)
            os.mkdir(data_dir)
        else:
            os.mkdir(data_dir)
            
        Sim.use_output_directory(data_dir)
        Sim.run(mp.at_every(at_every , meep_output), until = until)
        
        frames = []
        for filename in os.listdir(data_dir):
            if filename.endswith('.h5'):
                frame = int(filename[filename.find('-')+1:filename.find('.')])
                frames.append(frame)
        frames = np.sort(np.array(frames))
        
        for filename in os.listdir(data_dir):
            if filename.endswith('.h5'):
                frame = int(filename[filename.find('-')+1:filename.find('.')])
                if all_frames:
                    f = h5py.File(data_dir+'/'+filename, 'r')
                    field = f[component][...]
                    print('Plotting from '+filename)
                    plt.imshow(intensity*np.flipud(field.T), cmap='RdBu',\
                               norm=mpl.colors.Normalize(vmin=-1, vmax=1))
                    plt.axis('off')
                    plt.savefig(data_dir+'/'+str(frame)+'.png', dpi = 250)
                else:
                    if frame == frames[len(frames)-1]:
                        f = h5py.File(data_dir+'/'+filename, 'r')
                        field = f[component][...]
                        print('Plotting from '+filename)
                        plt.imshow(intensity*np.flipud(field.T), cmap='RdBu',\
                                   norm=mpl.colors.Normalize(vmin=-1, vmax=1))
                        plt.axis('off')
                        plt.savefig(data_dir+'/'+str(frame)+'.png', dpi = 250)
            
            if remove and filename.endswith('.h5'):
                os.remove(data_dir+'/'+filename)
                
                
    def Sim_And_Plot_Gif(self, at_every, field, until, output_dir,\
                     component = 'hz', out_prefix = '', remove = False,\
                     filetype = 'gif', intensity = 1.0):
        Sim = self.Get_Sim()
        
        if field == 'E':
            meep_output = mp.output_efield
        elif field == 'H':
            meep_output = mp.output_hfield
        elif field == 'S':
            meep_output = mp.output_poynting
            raise RuntimeError("Haven't worked out plotting for \
                                Poynting vector yet.")

        data_dir = output_dir+out_prefix+'-out/'
        img_dir = output_dir+out_prefix+'-out/frames'
        if os.path.isdir(data_dir):
            shutil.rmtree(data_dir)
            os.mkdir(data_dir)
        else:
            os.mkdir(data_dir)
        if os.path.isdir(img_dir):
            shutil.rmtree(img_dir)
            os.mkdir(img_dir)
        else:
            os.mkdir(img_dir)
            
        Sim.use_output_directory(data_dir)
        Sim.run(mp.at_every(at_every , meep_output), until = until)
        
        frames = []
        for filename in os.listdir(data_dir):
            if filename.endswith('.h5'):
                frame = int(filename[filename.find('-')+1:filename.find('.')])
                frames.append(frame)
                f = h5py.File(data_dir+'/'+filename, 'r')
                field = f[component][...]
                print('Plotting from '+filename)
                plt.imshow(intensity*np.flipud(field.T), cmap='RdBu',\
                           norm=mpl.colors.Normalize(vmin=-1, vmax=1))
                plt.axis('off')
                plt.savefig(img_dir+'/'+str(frame)+'.png', dpi = 250)
            
            if remove and filename.endswith('.h5'):
                os.remove(data_dir+'/'+filename)
               
        # Create new figure for GIF
        fig, ax = plt.subplots()

        # Adjust figure so GIF does not have extra whitespace
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
        ax.axis('off')
        
        ims = []
        frames = np.sort(np.array(frames))
        for frame in frames:
            im = ax.imshow(plt.imread(img_dir+'/'+str(frame)+'.png'), animated = True)
            ims.append([im])

        ani = animation.ArtistAnimation(fig, ims)
        if filetype != 'gif':
            FFwriter = animation.FFMpegWriter(fps=6)
            ani.save(data_dir+'/fields.'+filetype, writer=FFwriter)
        else:
            ani.save(data_dir+'/fields.'+filetype)
            
                
    def Plot_Results(self, at_every, field, until, output_dir,\
                     component = 'hz', out_prefix = '', remove = False):
        data_dir = output_dir+out_prefix+'-out/'
        
        for filename in os.listdir(data_dir):
            if filename.endswith('.h5'):
                f = h5py.File(data_dir+'/'+filename, 'r')
                field = f[component][...]
                print('Plotting from '+filename)
                plt.imshow(field, cmap='RdBu',\
                           norm=mpl.colors.Normalize(vmin=-1, vmax=1))
                plt.savefig(data_dir+'/'+filename[:-3]+'.pdf', dpi = 1500)
            
            if remove:
                os.remove(data_dir+'/'+filename)

    ############################################################################
    ## Dimensionalization
    ############################################################################
    def Nondimensionalize_Freq(self, f):
        """
        Gives frequency in a units given f in Hz
        """
        f0 = c/self.a

        return f/f0

    
    def Dimensionalize_Freq(self, f):
        """
        Gives frequency in Hz given f in a units
        """
        f0 = c/self.a

        return f*f0


    def Nondimensionalize_B(self, B):
        """
        Gives cyclotron frequency in a units given B in T
        """
        fc = e*B/me/2/np.pi
        f0 = c/self.a

        return fc/f0


    def Dimensionalize_B(self, wc):
        fc = wc*c/self.a
        
        return fc*me*2*np.pi/e


