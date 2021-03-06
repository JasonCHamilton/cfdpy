"""
Output VTK files.

    The file contains routines for outputting grid and solution data
    in the VTK file format(s).
"""


import numpy as np


def writegridtovtk(datatype, nx, ny, x, y):
    """Output grid data to vtk structured grid legacy format."""
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
    """Output fluid and grid data to vtk structured grid legacy format."""
    outvars = [
                fluid.Q[0],
                fluid.u,
                fluid.Et,
                fluid.ein,
                fluid.P,
                fluid.T,
                fluid.mw,
                fluid.R,
                fluid.mu,
                fluid.cp,
                fluid.cv,
                fluid.gamma,
                fluid.kf,
                fluid.tauxx,
                fluid.tauyy,
                fluid.tauxy,
                fluid.qx,
                fluid.qy
              ]

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
    for vin in xrange(len(outvars)):
        if fluid.vartype[vin] == 'SCALARS':
            writescalar(f, outvars[vin], fluid.varnames[vin],
                        fluid.datatype[vin], ny-1, nx-1)
        elif fluid.vartype[vin] == 'VECTORS':
            writevector(f, outvars[vin], fluid.varnames[vin],
                        fluid.datatype[vin], ny-1, nx-1)
    f.close()


def writescalar(f, var, varname, datatype, ny, nx):
    """Write variable data to file as integer or float."""
    f.write('SCALARS %s %s\n' % (varname, datatype))
    f.write('LOOKUP_TABLE default\n')
    for j in xrange(ny):
        for i in xrange(nx):
            if datatype == 'float':
                f.write('%f\n' % var[i, j])
            elif datatype == 'integer':
                f.write('%d\n' % var[i, j])
    return None


def writevector(f, var, varname, datatype, ny, nx):
    """Write vector data to file as integer or float."""
    f.write('VECTORS %s %s\n' % (varname, datatype))
    for j in xrange(ny):
        for i in xrange(nx):
            if datatype == 'float':
                f.write('%f %f %f\n' % (var[0, i, j], var[1, i, j], 0.0))
            elif datatype == 'integer':
                f.write('%d %d %d\n' % (var[0, i, j], var[1, i, j], 0))
    return None
