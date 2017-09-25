''' 2D Compressible Navier-Stokes Finite Volume Solver

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
'''


# import numpy as np
import grids.gridgenerator as gridgen
import fileio.outputvtk as outputvtk


# Inputs Section #
# Grid Inputs
nx = 4    # number of nodes in the x direction
ny = 4    # number of nodes in the y direction
ngc = 1   # number of ghost cells
Lx = 3.0  # Length in the x direction in physical units
Ly = 3.0  # Length in the y direction in physical units
datatype = 'ASCII'

grid = gridgen.generategrid_2d_uniform(Lx, Ly, nx, ny, ngc)

outputvtk.writegridtovtk(datatype, grid.ngx, grid.ngy, grid.x, grid.y)
