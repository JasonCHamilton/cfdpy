"""
MacCormack Flux.

    The file contains routines to compute the viscous and inviscid fluxes
    using MacCormack's central differencing predictor corrector methods.
"""


import numpy as np


def computefluxvectors(fluid, grid):
    """Compute flux vectors at cell centers."""
    # Compute inviscid flux vector in the x and y direction
    # last argument is a direction flag, 0 = x, 1 = y
    EI = computeinviscidflux(fluid, 0)
    FI = computeinviscidflux(fluid, 1)

    # Compute viscous flux vector in the x and y direction
    # last argument is a direction flag, 0 = x, 1 = y
    EV = computeviscousflux(fluid, grid, 0)
    FV = computeviscousflux(fluid, grid, 1)

    # Combine the inviscid and visous flux vectors
    fluid.E = EI + EV
    fluid.F = FI + FV

    return fluid


def computeinviscidflux(fluid, direction):
    """Compute inviscid flux vector at cell centers."""
    # Initialize flux vector
    flux = np.zeros(np.shape(fluid.E))
    # check which direction to compute flux at cell centers
    if direction == 0:
        # Compute mass flux
        flux[0] = fluid.Q[1]  # rhou
        # Compute x momentum flux in the x direction
        flux[1] = fluid.Q[1]*fluid.u[0] + fluid.P
        # Compute y momentum flux in the x direction
        flux[2] = fluid.Q[2]*fluid.u[0]
        # Compute Energy flux in the x direction
        flux[3] = fluid.Q[3]*fluid.u[0] + fluid.P*fluid.u[0]
    elif direction == 1:
        # Compute mass flux
        flux[0] = fluid.Q[2]  # rhou
        # Compute x momentum flux in the y direction
        flux[1] = fluid.Q[1]*fluid.u[1]
        # Compute y momentum flux in the y direction
        flux[2] = fluid.Q[2]*fluid.u[1] + fluid.P
        # Compute Energy flux in the y direction
        flux[3] = fluid.Q[3]*fluid.u[1] + fluid.P*fluid.u[1]

    return flux


def computeviscousflux(fluid, grid, direction):
    """Compute inviscid flux vector at cell centers."""
    # Initialize flux vector
    flux = np.zeros(np.shape(fluid.E))
    # check which direction to compute flux at cell centers
    if direction == 0:
        # Compute mass flux
        # flux[0] = 0.0  # no viscous terms in the continuity equation
        # Compute x momentum flux in the x direction
        flux[1] = -fluid.tauxx
        # Compute y momentum flux in the x direction
        flux[2] = -fluid.tauxy
        # Compute Energy flux in the x direction
        flux[3] = -fluid.u[0]*fluid.tauxx - fluid.u[1]*fluid.tauxy + fluid.qx
    elif direction == 1:
        # Compute mass flux
        # flux[0] =  0.0
        # Compute x momentum flux in the y direction
        flux[1] = -fluid.tauxy
        # Compute y momentum flux in the y direction
        flux[2] = -fluid.tauyy
        # Compute Energy flux in the y direction
        flux[3] = -fluid.u[0]*fluid.tauxy - fluid.u[1]*fluid.tauyy + fluid.qy

    return flux


def computeviscousstress(mu, u, v, dzetadx, detadx, dzetady, detady):
    """Compute viscous stress at cell centers."""
    # Compute velocity gradients
    dudzeta, dudeta = computecellgradients(u)
    dvdzeta, dvdeta = computecellgradients(v)

    # Compute shear stress in the xy plan
    tauxy = mu*(dzetady*dudzeta + detady*dudeta
                + dzetadx*dvdzeta + detadx*dudeta)
    tauyy = mu*(2.0/3.0)*(2.0*(dzetady*dvdzeta + detady*dvdeta)
                          - (dzetadx*dudzeta + detadx*dudeta))
    tauxx = mu*(2.0/3.0)*(2.0*(dzetadx*dudzeta + detadx*dudeta)
                          - (dzetady*dvdzeta + detady*dvdeta))

    return tauxy, tauxx, tauyy


def computehearflux(T, kf, dzetadx, detadx, dzetady, detady):
    """Compute hear flux at cell centers."""
    # Compute temperature gradients
    dTdzeta, dTdeta = computecellgradients(T)

    # Compute heat flux
    qx = -kf*(dzetadx*dTdzeta + detadx*dTdeta)
    qy = -kf*(dzetady*dTdzeta + detady*dTdeta)

    return qx, qy


def computecellgradients(phi):
    """Compute gradient of cell centered variable w.r.t. grid coordinates."""
    # Initialize gradient Variables
    dphidzeta = np.zeros(np.shape(phi))
    dphideta = np.zeros(np.shape(phi))
    # Compute gradient for interier cell using 2nd order central differencing
    dphidzeta[1:-1, :] = 0.5*(phi[2:, :]-phi[:-2, :])
    dphideta[:, 1:-1] = 0.5*(phi[:, 2:]-phi[:, :-2])
    # Compute gradient for imin and jmin cells
    # using 2nd order forward differencing
    dphidzeta[0, :] = 0.5*(-3.0*phi[0, :] + 4.0*phi[1, :] - phi[2, :])
    dphideta[:, 0] = 0.5*(-3.0*phi[:, 0] + 4.0*phi[:, 1] - phi[:, 2])

    # Compute gradient for imin and jmin cells
    # using 2nd order forward differencing
    dphidzeta[-1, :] = 0.5*(3.0*phi[-1, :] - 4.0*phi[-2, :] + phi[-3, :])
    dphideta[:, -1] = 0.5*(3.0*phi[:, -1] - 4.0*phi[:, -2] + phi[:, -3])

    return dphidzeta, dphideta
