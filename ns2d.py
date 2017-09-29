"""
2D Compressible Navier-Stokes Finite Volume Solver.

    This code is a DNS solver for the 2 dimentional Navier-Stokes equations
    using the finite volume method. The solver uses single block structured
    grids.

    Temporal Schemes
    - 2nd Order MacCormack Predictor-Corrector

    Spacial Schemes
    - 2nd Order Central MacCormack

    Thermodynamics Models
    - Calorical Perfect Gas (CPG)

    Transport Models
    - Sutherland's Law

    Boundar Conditions
    - Inlets
        + Dirchlet: prescribe constant values for velocity, temperature,
                    and pressure; density is computed from thermodynamics
    - Outlets
        + Neumann: gradients at outlet face set to zero.
                   ghost cell set to interier cells
    - Walls
        + No slip adiabatic: velocity at face set to zero,
                             temperature gradients at face set to zero.
        + Slip adiabatic: velocity gradients at face set to zero,
                          temperature gradients at face set to zero.
    Development Notes:
    - Add Species Equation
    - Add TPG to Thermodynamics Models
    - Add 3rd Order MUSCL to Spacial Schemes
    - Add chemical reactions/kinetics
    - Add RK4 to Temporal Schemes
    - Add Peng-Robinson to Thermodynamics Models
    - Add Chungs's (1988) high density model
    - Add input file support (most likely xml format)

    This is code is being developed as a learning tool. Will likely be
    rewitten as sepereate program when transitioned into 3D. Once 3D single
    block had been validated/verified then multiblock functionality will be
    added. The goal is grain a better understanding of core principles in CFD
    and scientific computing, as well as improving programming efficiency.

    The end product should provide a platform for rapid model development and
    testing. Once new models/features are tested using python solvers they can
    be implimented into fortran solvers.
"""


if __name__ == "__main__":

    import numpy as np
    import grids.gridgenerator as gridgen
    import fileio.outputvtk as outputvtk
    import simsetup.initialconditions as ics
    from flowvars.flowvariables import BC2D

    # Inputs Section #
    # Grid Inputs
    nx = 4        # number of nodes in the x direction
    ny = 4        # number of nodes in the y direction
    ngc = 1       # number of ghost cells
    Lx = 0.00168  # Length in the x direction in physical units
    Ly = 0.00168  # Length in the y direction in physical units
    datatype = 'ASCII'
    # Initial Condition Inputs
    T0 = 300.0       # temperature in K
    P0 = 101325.0    # pressure in Pa
    U0 = 0.0         # x velocity in m/s
    V0 = 0.0         # y velocity in m/s
    # Boundary Conditions
    sim_bc = BC2D()  # initilize simulation boundary conditions
    # imin boundary condition
    sim_bc.imin_type = 7
    # imax boundary condition
    sim_bc.imax_type = 7
    # jmin boundary condition
    sim_bc.imin_type = 7
    # jmax boundary condition
    sim_bc.imax_type = 1
    sim_bc.jmax_u = 1.0
    sim_bc.jmax_v = 0.0
    sim_bc.jmax_T = 300.0
    sim_bc.jmax_P = 101325.0
    # Thermodynamic properies
    mw = 28.938        # molecular weight in g/mol
    gamma = 1.4        # heat capacity ratio
    # Transport properies
    T_ref = 273.15     # Reference temperature in K
    Pr_ref = 0.72       # Prandtl number
    S_ref = 110.4      # Sutherland temperature in K
    mu_ref = 1.716e-5  # Reference viscosity
    Cmu = mu_ref*(T_ref+S_ref)/(T_ref*np.sqrt(T_ref))  # Sutherland Coefficient

    grid = gridgen.generategrid_2d_uniform(Lx, Ly, nx, ny, ngc)

    outputvtk.writegridtovtk(datatype, grid.ngx, grid.ngy, grid.x, grid.y)

    fluid = ics.initializefluid(nx, ny, ngc, T0, P0, U0, V0, mw, gamma,
                                T_ref, Pr_ref, S_ref, mu_ref, Cmu,
                                grid.dzetadx_cell, grid.detadx_cell,
                                grid.dzetady_cell, grid.detady_cell)

    outputvtk.writesolutiontovtk(0, datatype, grid.ngx, grid.ngy,
                                 grid.x, grid.y, fluid)
