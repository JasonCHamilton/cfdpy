"""
Initial Conditions.

    The file contains the routines used to compute initial conditions using
    user inputs.
"""


from ..flowvars.flowvariables import Fluid2D

from ..thermodynamics.cpg import (computedensity, computeinternalenergy,
                                  computetotalenergy)


def initializefluid(nx, ny, ngc, Tref, Pref, mwref, gamma, Pr_ref, Uref, Vref):
    """Create and initialize solution values from user inputs."""
    fluid = Fluid2D(nx, ny, ngc)

    # Prandtl number with defualt values set to 1
    fluid.Pr = Pr_ref

    # mw is the molecular weight at cell centers
    fluid.mw[:, :] = mwref

    # gamma is the specific heat ratio (cp/cv) at cell centers
    fluid.gamma[:, :] = gamma

    # R is the specific gas constant at cell centers
    fluid.R = fluid.Ru/fluid.mu  # (m^2)/(s^2*K)

    # P is the pressure at cell centers
    fluid.P[:, :] = Pref

    # T is the temperature at cell centers
    fluid.T[:, :] = Tref

    # u is the velocity vector at cell centers
    # u[0] = u  (x component)
    # u[1] = v  (y component)
    fluid.u[0, :, :] = Uref
    fluid.u[1, :, :] = Vref

    # Q[0] is density at cell ceners
    fluid.Q[0] = computedensity(fluid.P, fluid.T, fluid.R)
    # Q[1] is x momentum (rhou) at cell ceners
    fluid.Q[1] = fluid.Q[0]*fluid.u[0]
    # Q[2] is y momentum (rhou) at cell ceners
    fluid.Q[2] = fluid.Q[0]*fluid.u[1]

    # ein is the internal energy per unit mass
    fluid.ein = computeinternalenergy(fluid.P, fluid.gamma, fluid.Q[0])

    # Et is the total Energy per unit mass
    fluid.Et = computetotalenergy(fluid.ein, fluid.u[0], fluid.u[1])

    # Q[3] is total energy per at cell ceners
    fluid.Q[3] = fluid.Q[0]*fluid.u[1]

    # F is the flux vector for the conservative values in the i direction
    fluid.F = np.zeros((4, nxc, nyc))

    # E is the flux vector for the conservative values in the j direction
    fluid.E = np.zeros((4, nxc, nyc))

    # S is the volumetric source term for the conservative values
    # at the cell centers
    fluid.S = np.zeros((4, nxc, nyc))

    # cp is the specific heat at constant pressure at cell centers
    fluid.cp = np.zeros((nxc, nyc))

    # cv is the specific heat at constant volume at cell centers
    fluid.cv = np.zeros((nxc, nyc))

    # kf is the thermal conductivity at cell centers
    fluid.kf = np.zeros((nxc, nyc))

    # mu is the viscosity at cell centers
    fluid.mu[:, :]

    return fluid
