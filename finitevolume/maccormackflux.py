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
    computeshearstress(fluid, grid)
    computehearflux(fluid, grid)
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


def computeshearstress(mu, u, u, dzetadx, detadx, dzetady, detady):
    """Compute shear stress at cell centers."""
    # Compute velocity gradients
    dudzeta = 0.5*(u[2:, 1:-1]-u[:-2, 1:-1])
    dudeta = 0.5*(u[1:-1, 2:]-u[1:-1, :-2])
    dydzeta = 0.5*(v[2:, 1:-1]-v[:-2, 1:-1])
    dudeta = 0.5*(v[1:-1, 2:]-v[1:-1, :-2])
    # Compute shear stress in the xy plan

    return tauxy, tauxx, tauyy


def computehearflux(fluid, grid):
    """Compute hear flux at cell centers."""
    return
