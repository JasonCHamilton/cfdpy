"""
Flow Variables.

    The file contains the clases used to store the simulation data.
"""

import numpy as np


# Class to cotain 2D grid data
class Fluid2D(object):
    """Class used for 2D structed simulation variables."""

    def __init__(self, nx, ny, ngc):
        """Initialize Grid2D object."""
        # ngc is number of ghost cells, nxc and nyc are number of cells
        # including ghost nodes
        nxc = nx + 2*ngc - 1
        nyc = ny + 2*ngc - 1

        self.nxc = nxc
        self.nyc = nyc
        self.ngc = ngc

        # Prandtl number with defualt values set to 1
        self.Pr = 1.0

        # Specific Gas Constant with default set to Air
        self.Ru = 8314.4598  # (g*m^2)/(s^2*K*mol)

        # Q is the solution vector for the conservative values
        # Q[0] = rho
        # Q[1] = rho*u
        # Q[2] = rho*v
        # Q[3] = rho*Et
        self.Q = np.zeros((4, nxc, nyc))

        # F is the flux vector for the conservative values in the i direction
        self.F = np.zeros((4, nxc, nyc))

        # E is the flux vector for the conservative values in the j direction
        self.E = np.zeros((4, nxc, nyc))

        # S is the volumetric source term for the conservative values
        # at the cell centers
        self.S = np.zeros((4, nxc, nyc))

        # u is the velocity vector at cell centers
        # u[0] = u  (x component)
        # u[1] = v  (y component)
        self.u = np.zeros((2, nxc, nyc))

        # E is the total Energy per unit mass
        self.Et = np.zeros((nxc, nyc))

        # ein is the internal energy per unit mass
        self.ein = np.zeros((nxc, nyc))

        # P is the pressure at cell centers
        self.P = np.zeros((nxc, nyc))

        # T is the temperature at cell centers
        self.T = np.zeros((nxc, nyc))

        # mw is the molecular weight at cell centers
        self.mw = np.zeros((nxc, nyc))

        # R is the specific gas constant at cell centers
        self.R = np.zeros((nxc, nyc))

        # mu is the viscosity at cell centers
        self.mu = np.zeros((nxc, nyc))

        # cp is the specific heat at constant pressure at cell centers
        self.cp = np.zeros((nxc, nyc))

        # cv is the specific heat at constant volume at cell centers
        self.cv = np.zeros((nxc, nyc))

        # gamma is the specific heat ratio (cp/cv) at cell centers
        self.gamma = np.zeros((nxc, nyc))

        # kf is the thermal conductivity at cell centers
        self.kf = np.zeros((nxc, nyc))

        # tauxx is the normal stress in the x direction
        self.tauxx = np.zeros((nxc, nyc))

        # tauyy is the normal stress in the y direction
        self.tauyy = np.zeros((nxc, nyc))

        # tauxy=tauyx is the normal stress in the xy plane
        self.tauxy = np.zeros((nxc, nyc))

        # qx is the heat transfer in the x direction
        self.qx = np.zeros((nxc, nyc))

        # qy is the heat transfer in the y direction
        self.qy = np.zeros((nxc, nyc))

        # dictionary of variables names
        self.varnames = [
                         "Density[kg/m^3]",
                         "Velocity[m/s]",
                         "TotalEnergy[m^2/s^2]",
                         "InternalEnergy[m^2/s^2]",
                         "Pressure[Pa]",
                         "Temperature[K]",
                         "MolecularWeight[g/mol]",
                         "SpecificGasConstant[J/(kg*K)]",
                         "Viscosity[kg/(m*s)]",
                         "Cp[J/(kg*K)]",
                         "Cv[J/(kg*K)]",
                         "SpecificHeatRatio",
                         "ThermalConductivity[W/(m*K)]",
                         "Tau_xx[Pa]",
                         "Tau_yy[Pa]",
                         "Tau_xy[Pa]",
                         "q_x[W/m^2]",
                         "q_y[W/m^2]"
                        ]
        self.vartype = [
                        "SCALARS",
                        "VECTORS",
                        "SCALARS",
                        "SCALARS",
                        "SCALARS",
                        "SCALARS",
                        "SCALARS",
                        "SCALARS",
                        "SCALARS",
                        "SCALARS",
                        "SCALARS",
                        "SCALARS",
                        "SCALARS",
                        "SCALARS",
                        "SCALARS",
                        "SCALARS",
                        "SCALARS",
                        "SCALARS"
                       ]
        self.datatype = [
                         "float",
                         "float",
                         "float",
                         "float",
                         "float",
                         "float",
                         "float",
                         "float",
                         "float",
                         "float",
                         "float",
                         "float",
                         "float",
                         "float",
                         "float",
                         "float",
                         "float",
                         "float"
                        ]


# Class to cotain 2D grid data
class BC2D(object):
    """Class used for 2D Boundary Conditions."""

    def __init__(self):
        """Initialize Boundary Conditions"""
        # BC types:
        # 0 - subsonic velocity inlet
        # 1 - supersonic velocity intlet
        # 2 - subsonic mass flow inlet
        # 3 - supersonic mass flow intlet
        # 4 - subsonic outlet
        # 5 - supersonic outlet
        # 6 - slip wall
        # 7 - no slip wall
        # all boundaries deflaut to no slip walls
        # imin boundary
        self.imin_type = 7
        self.imin_u = 1.0
        self.imin_v = 0.0
        self.imin_T = 300.0
        self.imin_P = 101325.0
        self.imin_massflow = 1.0
        # imax boundary
        self.imax_type = 7
        self.imax_u = 1.0
        self.imax_v = 0.0
        self.imax_T = 300.0
        self.imax_P = 101325.0
        self.imax_massflow = 1.0
        # jmin boundary
        self.jmin_type = 7
        self.jmin_u = 1.0
        self.jmin_v = 0.0
        self.jmin_T = 300.0
        self.jmin_P = 101325.0
        self.jmin_massflow = 1.0
        # imax boundary
        self.jmax_type = 7
        self.jmax_u = 1.0
        self.jmax_v = 0.0
        self.jmax_T = 300.0
        self.jmax_P = 101325.0
        self.jmax_massflow = 1.0
