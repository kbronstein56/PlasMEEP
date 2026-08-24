"""
PMMCirculatorInverse.py

PlasMEEP / Meep forward model for the in-silico PMM circulator.

The first goal is to reproduce the PMM forward model with B = 0,
then turn on the Helmholtz-coil magnetic field B and calculate nonreciprocal
S-parameters.

IMPORTANT:
    - Ceviche is NOT used here.  Meep is the electromagnetic solver.
    - PlasMEEP supplies the plasma material/geometry layer on top of Meep.
    - rho still parameterizes the plasma state; rho -> plasma frequency is kept
      in the same style as the existing PMM code.
    - B is a physical magnetic-field input to the PlasMEEP simulation.
    - The port/S-parameter methods are left as explicit TODOs until the original
      PMMInverse.py implementation is ported exactly.
"""

import os
import numpy as np

from plasmeep.lib import Plasmeep as pm


###############################################################################
## Utility functions and globals
###############################################################################

c = 299792458
e = 1.60217662e-19
epso = 8.8541878128e-12
muo = 4*np.pi*1e-7
me = 9.1093837015e-31


###############################################################################
## In-silico inverse design class
###############################################################################

class PMMI:
    """
    In-silico PMM circulator model using PlasMEEP / Meep.

    The naming is intentionally kept close to PMMInverse.py / PMMInSitu.py so
    it is easy to compare the old Ceviche simulation, the new Meep simulation,
    and the eventual in-situ implementation.
    """

    def __init__(self, a=0.01, res=50, dpml=1, nx=10, ny=8,
                 B=np.array([0.0, 0.0, 0.0])):
        """
        Args:
            a: dimensionalized unit length in meters
            res: Meep simulation resolution (pixels per a)
            dpml: PML thickness in a units
            nx: simulation x-size in a units
            ny: simulation y-size in a units
            B: applied magnetic-field vector. For the circulator this will
               normally be [0, 0, Bz].
        """

        self.a = a
        self.res = res
        self.dpml = dpml
        self.nx = nx
        self.ny = ny
        self.B = np.asarray(B, dtype=float)

        self.sim = None

        # Keep track of the current PMM configuration.
        self.rho = None
        self.wp = None
        self.gamma = None
        self.bulb_centers = None

        # Added to stay compatible with PMMInverse.py-style trainable arrays.
        self.train_elems = []
        self.train_elem_locs = []


    ###########################################################################
    ## Frequency / plasma-parameter helpers
    ###########################################################################

    def f_GHz(self, f):
        """
        Returns dimensionalized frequency in GHz.

        Args:
            f: frequency in a units
        """
        return f*c/self.a/1e9


    def f_a(self, f):
        """
        Returns nondimensionalized frequency in a units.

        Args:
            f: frequency in GHz
        """
        return f*1e9/c*self.a


    def Scale_Rho_wp(self, rho, w_src=None, wp_max=0, gamma=0):
        """
        Uses the same arctan-barrier style as the existing PMM code to map rho
        to NONDIMENSIONALIZED plasma frequency.

        Args:
            rho: parameters being optimized
            w_src: source frequency in a units. Kept for compatibility with
                   PMMInverse.py; it is not needed to calculate wp itself.
            wp_max: approximate maximum nondimensionalized plasma frequency
            gamma: collision frequency in a units. Kept for compatibility with
                   PMMInverse.py; it is not needed to calculate wp itself.

        Returns:
            wp: nondimensionalized plasma frequency for each trainable element
            elem_locations: x,y locations of the trainable elements

        Notes:
            In the old PMMInverse.py, rho was mapped to wp and then immediately
            mapped to scalar epsilon for Ceviche. Here we stop at wp because
            PlasMEEP/Meep uses wp, gamma, and B to build the plasma material.
        """
        rho = np.asarray(rho, dtype=float).flatten()

        if wp_max > 0:
            wp = (wp_max/1.5)*np.arctan(rho/(wp_max/7.5))
        else:
            wp = rho

        wp = np.abs(wp)

        elem_locations = np.asarray(self.train_elem_locs, dtype=float)

        return wp, elem_locations


    def Scale_Rho_fp(self, rho, wp_max):
        """
        Uses an arctan barrier to map optimal parameters to dimensionalized
        plasma frequency in GHz.

        Args:
            rho: parameters being optimized
            wp_max: approximate maximum nondimensionalized plasma frequency
        """
        fp, _ = self.Scale_Rho_wp(rho, wp_max=wp_max)
        fp_dim = fp*c/self.a/1e9
        return fp_dim


    def Scale_Rho_ne(self, rho, wp_max):
        """
        Maps rho to electron density in m^-3 through the plasma frequency.

        Args:
            rho: parameters being optimized
            wp_max: approximate maximum nondimensionalized plasma frequency
        """
        wp, _ = self.Scale_Rho_wp(rho, wp_max=wp_max)

        # PlasMEEP/Meep frequency normalization:
        # nondimensional f -> Hz = f*c/a.
        # Convert ordinary plasma frequency fp to angular frequency omega_p.
        fp_Hz = wp*c/self.a
        omega_p = 2*np.pi*fp_Hz

        ne = omega_p**2 * me * epso / e**2
        return ne


    ###########################################################################
    ## Trainable array geometry
    ###########################################################################

    def Add_Rod_train(self, r, center):
        """
        Add a trainable plasma rod location.

        This keeps the same general role/name as PMMInverse.py, but for now
        stores geometry information instead of building a Ceviche pixel mask.
        """
        center = np.asarray(center, dtype=float)

        self.train_elem_locs.append([center[0], center[1]])
        self.train_elems.append({
            "r": r,
            "center": center
        })

        return


    def Rod_Array_Hexagon_train(self, xy_cen, side_dim, r, d,
                           a_basis=np.array([[0,1],
                                             [np.sqrt(3)/2,1./2]]),
                           bulbs=False, r_bulb=(0, 0), eps_bulb=3.8,
                           uniform=True):
        """
        Add a hexagonal triangular array of trainable plasma elements.

        This keeps the same placement logic and function name as PMMInverse.py.

        Args:
            xy_cen: np.array, center of hexagon
            side_dim: number of rods along one side of hexagon
            r: radius of the rod in a units
            d: array spacing in a units
            a_basis: basis vectors that determine orientation of hexagon
            bulbs: kept for compatibility; physical bulb geometry is added
                   separately with Add_Bulb_Array in this early version
            r_bulb: kept for compatibility
            eps_bulb: kept for compatibility
            uniform: kept for compatibility
        """
        xy_cen = np.asarray(xy_cen, dtype=float)
        a_basis = np.asarray(a_basis, dtype=float)

        for i in range(side_dim):
            for j in range(side_dim+i):
                b1_loc = xy_cen-(side_dim-1-i)*d*a_basis[0,:]-\
                        (i-j)*d*a_basis[1,:]

                b2_loc = xy_cen+(side_dim-1-i)*d*a_basis[0,:]+\
                        (i-j)*d*a_basis[1,:]

                self.Add_Rod_train(r, b1_loc)

                if i < side_dim - 1:
                    self.Add_Rod_train(r, b2_loc)

        return


    ###########################################################################
    ## Build the PlasMEEP simulation
    ###########################################################################

    def Build_Sim(self, B=None):
        """
        Creates a fresh PlasMEEP simulation.

        A fresh simulation is useful because both rho and B can change between
        inverse-design evaluations.

        Args:
            B: optional new magnetic-field vector. If None, use self.B.
        """
        if B is not None:
            self.B = np.asarray(B, dtype=float)

        self.sim = pm(
            a=self.a,
            res=self.res,
            dpml=self.dpml,
            nx=self.nx,
            ny=self.ny,
            B=self.B
        )

        return self.sim


    def Add_Bulb_Array(self, rho, wp_max, bulb_centers, r_bulb,
                       gamma=0.0, axis=np.array([0, 0, 1]), profile=0):
        """
        Adds the full plasma-bulb array to the current PlasMEEP simulation.

        Args:
            rho: flattened array of PMM design parameters, one per bulb
            wp_max: approximate maximum nondimensionalized plasma frequency
            bulb_centers: iterable of [x, y, z] bulb-center coordinates in a units
            r_bulb: bulb-radius tuple expected by PlasMEEP Add_Bulb()
            gamma: collision frequency in PlasMEEP/Meep units
            axis: bulb cylinder axis
            profile: PlasMEEP bulb-profile option
        """
        if self.sim is None:
            self.Build_Sim()

        rho = np.asarray(rho, dtype=float).flatten()
        bulb_centers = np.asarray(bulb_centers, dtype=float)

        if rho.shape[0] != bulb_centers.shape[0]:
            raise ValueError(
                "rho must contain one value for every bulb center: "
                f"{rho.shape[0]} rho values vs. {bulb_centers.shape[0]} centers."
            )

        wp, _ = self.Scale_Rho_wp(rho, wp_max=wp_max)

        self.rho = np.copy(rho)
        self.wp = np.copy(wp)
        self.gamma = gamma
        self.bulb_centers = np.copy(bulb_centers)

        for i in range(rho.shape[0]):
            self.sim.Add_Bulb(
                r_bulb=r_bulb,
                center=bulb_centers[i],
                wp=wp[i],
                gamma=gamma,
                axis=axis,
                profile=profile
            )

        return


    ###########################################################################
    ## Sources / probes
    ###########################################################################

    def Add_Source(self, freq, pol, src_loc, src_size):
        """
        Adds a continuous-wave source.

        This preserves the familiar Add_Source-style interface for initial
        B=0 validation and field visualization.

        IMPORTANT:
            The final circulator S-parameter model will probably replace this
            with the exact source/port treatment used in PMMInverse.py, ported
            to Meep.
        """
        if self.sim is None:
            self.Build_Sim()

        self.sim.Add_Cont_Source(freq, pol, src_loc, src_size)
        return


    def Add_Probe(self, *args, **kwargs):
        """
        TODO: Port PMMInverse.py Add_Probe() to a Meep monitor/eigenmode monitor.

        This is deliberately NOT guessed here because the exact old probe
        geometry and S-parameter normalization need to remain comparable to
        PMMInverse.py.
        """
        raise NotImplementedError(
            "Add_Probe() needs the original PMMInverse.py implementation "
            "before it is ported to Meep."
        )


    def Get_Trans_Denom(self, *args, **kwargs):
        """
        TODO: Port PMMInverse.py transmission-normalization method.

        The normalization philosophy should remain the same as the existing
        PMM simulation, while the underlying field solver changes to Meep.
        """
        raise NotImplementedError(
            "Get_Trans_Denom() needs to be ported from PMMInverse.py."
        )


    def Calc_SParams_opt(self, *args, **kwargs):
        """
        TODO: Port PMMInverse.py S-parameter extraction to Meep.

        This will ultimately return the complex S-parameters required by the
        circulator objective, e.g. S21, S12, S32, S23, S13, S31.
        """
        raise NotImplementedError(
            "Calc_SParams_opt() needs to be ported from PMMInverse.py."
        )


    ###########################################################################
    ## Visualization / I/O
    ###########################################################################

    def Viz_Domain(self, output_dir, freq=0):
        """
        Wrapper around PlasMEEP Viz_Domain().
        """
        if self.sim is None:
            raise RuntimeError("Build the simulation before calling Viz_Domain().")

        os.makedirs(output_dir, exist_ok=True)

        # PlasMEEP's tutorial uses Viz_Domain(output_dir) without freq.
        # Pass freq only when it is explicitly nonzero.
        if freq == 0:
            return self.sim.Viz_Domain(output_dir)
        else:
            return self.sim.Viz_Domain(output_dir, freq=freq)


    def Save_Params(self, rho, savepath):
        """
        Wrapper for np.savetxt, matching the existing PMM code.
        """
        return np.savetxt(savepath, rho, delimiter=',')


    def Read_Params(self, readpath, iteration=0):
        """
        Reads optimization parameters, matching the existing PMM code style.

        Args:
            readpath: path to CSV file
            iteration: 0 -> full array, positive integer -> chosen row,
                       'last' -> final row
        """
        if iteration > 0:
            rho = np.loadtxt(readpath, delimiter=',')
            return rho[iteration-1, :]
        elif iteration == 'last':
            rho = np.loadtxt(readpath, delimiter=',')
            return rho[rho.shape[0]-1, :]
        else:
            return np.loadtxt(readpath, delimiter=',')


###############################################################################
## Minimal smoke test
###############################################################################

if __name__ == "__main__":
    # This only checks that the class can create a PlasMEEP simulation.
    # It does NOT yet run the full circulator or calculate S-parameters.

    pmm = PMMI(
        a=0.01,
        res=50,
        dpml=1,
        nx=10,
        ny=8,
        B=np.array([0.0, 0.0, 0.0])
    )

    pmm.Build_Sim()

    print("PMMCirculatorInverse initialized successfully.")
    print("B =", pmm.B)
