"""
Boundary Conditions.

    The file contains routines to apply boundary conditions to ghost cells.
"""


from thermodynamics.cpg import (computedensity, computeinternalenergy,
                                computetotalenergy, computecp, computecv)

from transport.sutherland import (computeviscosity,
                                  computethermalconductivity)

from finitevolume.maccormackflux import (computeviscousstress, computehearflux)


def applyboundaryconditions(fluid, sim_bc, Cmu, S_ref, dzetadx, detadx,
                            dzetady, detady):
    """Apply boundary conditions to ghost cells."""
    # Apply imin bc

    return fluid


def constantvelocityinlet(fluid, u, v, T, P, Cmu, S_ref, dzetadx, detadx,
                          dzetady, detady, ig, jg):
    """Apply constant velocity inlet boundary conditions to ghost cells."""
    fluid.u[0, ig, jg] = u
    fluid.v[1, ig, jg] = v
    fluid.T[ig, jg] = T
    fluid.P[ig, jg] = P
    fluid.Q[0, ig, jg] = computedensity(fluid.P[ig, jg], fluid.T[ig, jg],
                                        fluid.R[ig, jg])
    fluid.Q[1, ig, jg] = fluid.Q[0, ig, jg]*fluid.u[0, ig, jg]
    fluid.Q[2, ig, jg] = fluid.Q[0, ig, jg]*fluid.v[0, ig, jg]
    fluid.ein[ig, jg] = computeinternalenergy(fluid.P[ig, jg],
                                              fluid.gamma[ig, jg],
                                              fluid.Q[0, ig, jg])
    fluid.Et[ig, jg] = computetotalenergy(fluid.ein[ig, jg],
                                          fluid.u[0, ig, jg],
                                          fluid.u[1, ig, jg])
    fluid.Q[3, ig, jg] = fluid.Q[0, ig, jg]*fluid.Et[ig, jg]
    fluid.cp[ig, jg] = computecp(fluid.gamma[ig, jg], fluid.R[ig, jg])
    fluid.cv[ig, jg] = computecv(fluid.gamma[ig, jg], fluid.R[ig, jg])
    fluid.mu[ig, jg] = computeviscosity(fluid.T[ig, jg], Cmu, S_ref)
    fluid.kf[ig, jg] = computethermalconductivity(fluid.Pr, fluid.cp[ig, jg],
                                                  fluid.mu[ig, jg])
    (fluid.tauxx[ig, jg], fluid.tauyy[ig, jg], fluid.tauxy[ig, jg]) = \
        computeviscousstress(fluid.mu[ig, jg], fluid.u[0, ig, jg],
                             fluid.u[1, ig, jg], dzetadx[ig, jg],
                             detadx[ig, jg], dzetady[ig, jg], detady[ig, jg])
    fluid.qx[ig, jg], fluid.qy[ig, jg] = computehearflux(fluid.T[ig, jg],
                                                         fluid.kf[ig, jg],
                                                         dzetadx[ig, jg],
                                                         detadx[ig, jg],
                                                         dzetady[ig, jg],
                                                         detady[ig, jg])
    return fluid


def noslipwall(fluid, u, v, T, P, Cmu, S_ref, dzetadx, detadx, dzetady, detady,
               ig, jg, ic, jc):
    """No slip adiabatic wall boundary conditions to ghost cells."""
    fluid.u[0, ig, jg] = fluid.u[0, ic, jc]
    fluid.v[1, ig, jg] = fluid.v[1, ic, jc]
    fluid.T[ig, jg] = fluid.T[ic, jc]
    fluid.P[ig, jg] = fluid.P[ic, jc]
    fluid.Q[0, ig, jg] = fluid.Q[0, ic, jc]
    fluid.Q[1, ig, jg] = fluid.Q[1, ic, jc]
    fluid.Q[2, ig, jg] = fluid.Q[2, ic, jc]
    fluid.ein[ig, jg] = fluid.ein[ic, jc]
    fluid.Et[ig, jg] = fluid.Et[ic, jc]
    fluid.Q[3, ig, jg] = fluid.Q[3, ic, jc]
    fluid.cp[ig, jg] = fluid.cp[ic, jc]
    fluid.cv[ig, jg] = fluid.cv[ic, jc]
    fluid.mu[ig, jg] = fluid.mu[ic, jc]
    fluid.kf[ig, jg] = fluid.kf[ic, jc]
    fluid.tauxx[ig, jg] = fluid.tauxx[ic, jc]
    fluid.tauyy[ig, jg] = fluid.tauyy[ic, jc]
    fluid.tauxy[ig, jg] = fluid.tauxy[ic, jc]
    fluid.qx[ig, jg] = fluid.qx[ic, jc]
    fluid.qy[ig, jg] = fluid.qy[ic, jc]

    return fluid
