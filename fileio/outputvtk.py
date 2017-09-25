''' Output VTK files
    The file contains routines for outputting grid and solution data
    in the VTK file format(s).
'''


import numpy as np


def writegridtovtk(datatype, nx, ny, x, y):
    filename = 'grid.vtk'
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
