''' Output VTK files
    The file contains routines for outputting grid and solution data
    in the VTK file format(s).
'''


import numpy as np


def writegridtovtk(datatype, nx, ny, x, y):
    filename = '../output/grid.vtk'
    header = '# vtk DataFile Version 2.0'
    title = '2D structured grid'
    datatype = 'ASCII'  # For now we will for output to be in ASCII
    n = nx * ny
    dummy = np.zeros((nx, ny))

    # open file
    f = open(filename, 'w')

    f.write(header+'\n')
    f.write(title+'\n')
    f.write(datatype+'\n')
    f.write('DATASET STRUCTURED_GRID'+'\n')
    f.write('DIMENSIONS %d %d 1\n' % (nx, ny))
    f.write('POINTS %d float\n' % n)
    for j in xrange(ny):
        for i in xrange(nx):
            f.write('%f %f %f\n' % (x[i, j], y[i, j], dummy[i, j]))
    f.close()


def writesolutiontovtk(solnum, datatype, nx, ny, x, y, fluid):
    filename = '../output/solution_%05d.vtk' % solnum
    header = '# vtk DataFile Version 2.0'
    title = '2D structured grid solution'
    datatype = 'ASCII'  # For now we will for output to be in ASCII
    n = nx * ny
    nc = (nx-1)*(ny-1)
    dummy = np.zeros((nx, ny))

    # open file
    f = open(filename, 'w')

    f.write(header+'\n')
    f.write(title+'\n')
    f.write(datatype+'\n')
    f.write('DATASET STRUCTURED_GRID'+'\n')
    f.write('DIMENSIONS %d %d 1\n' % (nx, ny))
    f.write('POINTS %d float\n' % n)
    for j in xrange(ny):
        for i in xrange(nx):
            f.write('%f %f %f\n' % (x[i, j], y[i, j], dummy[i, j]))
    f.write('CELL_DATA %d\n' % nc)
    f.write('SCALARS Density[kg/m^3] float\n')
    f.write('LOOKUP_TABLE default\n')
    for j in xrange(ny-1):
        for i in xrange(nx-1):
            f.write('%f\n' % fluid.Q[0, i, j])
    f.write('SCALARS Pressure[Pa] float\n')
    f.write('LOOKUP_TABLE default\n')
    for j in xrange(ny-1):
        for i in xrange(nx-1):
            f.write('%f\n' % fluid.P[i, j])
    f.write('SCALARS Temperature[K] float\n')
    f.write('LOOKUP_TABLE default\n')
    for j in xrange(ny-1):
        for i in xrange(nx-1):
            f.write('%f\n' % fluid.T[i, j])
    f.write('SCALARS Viscosity[kg/(m*s)] float\n')
    f.write('LOOKUP_TABLE default\n')
    for j in xrange(ny-1):
        for i in xrange(nx-1):
            f.write('%f\n' % fluid.mu[i, j])
    f.write('SCALARS ThermalConductivity[W/(m*K)] float\n')
    f.write('LOOKUP_TABLE default\n')
    for j in xrange(ny-1):
        for i in xrange(nx-1):
            f.write('%f\n' % fluid.kf[i, j])
    f.write('SCALARS TotalEnergy[m^2/s^2] float\n')
    f.write('LOOKUP_TABLE default\n')
    for j in xrange(ny-1):
        for i in xrange(nx-1):
            f.write('%f\n' % fluid.Et[i, j])
    f.write('VECTORS Velocity[m/s] float\n')
    for j in xrange(ny-1):
        for i in xrange(nx-1):
            f.write('%f %f %f\n' % (fluid.u[0, i, j], fluid.u[1, i, j], 0.0))
    f.close()
