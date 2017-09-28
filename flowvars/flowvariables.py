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
        self.varnames = {
                         "Density[kg/m^3]": self.Q[0],
                         "Velocity[m/s]": self.u,
                         "TotalEnergy[m^2/s^2]": self.Et,
                         "InternalEnergy[m^2/s^2]": self.ein,
                         "Pressure[Pa]": self.P,
                         "Temperature[K]": self.T,
                         "MolecularWeight[g/mol]": self.mw,
                         "SpecificGasConstant[J/(kg*K)]": self.R,
                         "Viscosity[kg/(m*s)]": self.mu,
                         "Cp[J/(kg*K)]": self.cp,
                         "Cv[J/(kg*K)]": self.cv,
                         "SpecificHeatRatio": self.gamma,
                         "ThermalConductivity[W/(m*K)]": self.kf,
                         "Tau_xx[Pa]": self.tauxx,
                         "Tau_yy[Pa]": self.tauyy,
                         "Tau_xy[Pa]": self.tauxy,
                         "q_x[W/m^2]": self.qx,
                         "q_y[W/m^2]": self.qy,
                         }
        self.vartype = {
                         "Density[kg/m^3]": "SCALARS",
                         "Velocity[m/s]": "VECTORS",
                         "TotalEnergy[m^2/s^2]": "SCALARS",
                         "InternalEnergy[m^2/s^2]": "SCALARS",
                         "Pressure[Pa]": "SCALARS",
                         "Temperature[K]": "SCALARS",
                         "MolecularWeight[g/mol]": "SCALARS",
                         "SpecificGasConstant[J/(kg*K)]": "SCALARS",
                         "Viscosity[kg/(m*s)]": "SCALARS",
                         "Cp[J/(kg*K)]": "SCALARS",
                         "Cv[J/(kg*K)]": "SCALARS",
                         "SpecificHeatRatio": "SCALARS",
                         "ThermalConductivity[W/(m*K)]": "SCALARS",
                         "Tau_xx[Pa]": "SCALARS",
                         "Tau_yy[Pa]": "SCALARS",
                         "Tau_xy[Pa]": "SCALARS",
                         "q_x[W/m^2]": "SCALARS",
                         "q_y[W/m^2]": "SCALARS",
                         }
        self.datatype = {
                         "Density[kg/m^3]": "float",
                         "Velocity[m/s]": "float",
                         "TotalEnergy[m^2/s^2]": "float",
                         "InternalEnergy[m^2/s^2]": "float",
                         "Pressure[Pa]": "float",
                         "Temperature[K]": "float",
                         "MolecularWeight[g/mol]": "float",
                         "SpecificGasConstant[J/(kg*K)]": "float",
                         "Viscosity[kg/(m*s)]": "float",
                         "Cp[J/(kg*K)]": "float",
                         "Cv[J/(kg*K)]": "float",
                         "SpecificHeatRatio": "float",
                         "ThermalConductivity[W/(m*K)]": "float",
                         "Tau_xx[Pa]": "float",
                         "Tau_yy[Pa]": "float",
                         "Tau_xy[Pa]": "float",
                         "q_x[W/m^2]": "float",
                         "q_y[W/m^2]": "float",
                         }
