"""
Initial Conditions.

    The file contains the routines used to compute initial conditions using
    user inputs.
"""


from flowvars.flowvariables import Fluid2D

from thermodynamics.cpg import (computedensity, computeinternalenergy,
                                computetotalenergy, computecp, computecv)

from transport.sutherland import (computeviscosity,
                                  computethermalconductivity)

from finitevolume.maccormackflux import (computeviscousstress, computehearflux)


def initializefluid(nx, ny, ngc, T0, P0, U0, V0, mw, gamma, T_ref, Pr_ref,
                    S_ref, mu_ref, Cmu, dzetadx, detadx, dzetady, detady):
    """Create and initialize solution values from user inputs."""
    fluid = Fluid2D(nx, ny, ngc)

    # Prandtl number with defualt values set to 1
    fluid.Pr = Pr_ref

    # mw is the molecular weight at cell centers
    fluid.mw[:, :] = mw

    # gamma is the specific heat ratio (cp/cv) at cell centers
    fluid.gamma[:, :] = gamma

    # R is the specific gas constant at cell centers
    fluid.R = fluid.Ru/fluid.mw  # (m^2)/(s^2*K)

    # P is the pressure at cell centers
    fluid.P[:, :] = P0

    # T is the temperature at cell centers
    fluid.T[:, :] = T0

    # u is the velocity vector at cell centers
    # u[0] = u  (x component)
    # u[1] = v  (y component)
    fluid.u[0, :, :] = U0
    fluid.u[1, :, :] = V0

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
    fluid.Q[3] = fluid.Q[0]*fluid.Et

    # cp is the specific heat at constant pressure at cell centers
    fluid.cp = computecp(fluid.gamma, fluid.R)

    # cv is the specific heat at constant volume at cell centers
    fluid.cv = computecv(fluid.gamma, fluid.R)

    # mu is the viscosity at cell centers
    fluid.mu = computeviscosity(fluid.T, Cmu, S_ref)

    # kf is the thermal conductivity at cell centers
    fluid.kf = computethermalconductivity(fluid.Pr, fluid.cp, fluid.mu)

    # F is the flux vector for the conservative values in the i direction
    fluid.F[:, :, :] = 0.0  # Fluxes are initilized with zero

    # E is the flux vector for the conservative values in the j direction
    fluid.E[:, :, :] = 0.0  # Fluzes are initilized with zero

    # S is the volumetric source term for the conservative values
    # at the cell centers
    fluid.S[:, :] = 0.0  # Source terms are initilized with zero

    # tauxx is the normal stress in the x directio
    # tauyy is the normal stress in the y direction
    # tauxy=tauyx is the normal stress in the xy plane
    fluid.tauxx, fluid.tauyy, fluid.tauxy = computeviscousstress(fluid.mu,
                                                                 fluid.u[0],
                                                                 fluid.u[1],
                                                                 dzetadx,
                                                                 detadx,
                                                                 dzetady,
                                                                 detady)

    # qx is the heat transfer in the x direction
    # qy is the heat transfer in the y direction
    fluid.qx, fluid.qy = computehearflux(fluid.T, fluid.kf, dzetadx,
                                         detadx, dzetady, detady)

    return fluid
