"""
Grid Generators.

    The file contains the routines used to generate structured grids.
    Currently only supportes 2D structed grids.
"""

import numpy as np


# Class to cotain 2D grid data
class Grid2D(object):
    """Class used for 2D structed grids."""

    def __init__(self, nx, ny, ngc):
        """Initialize Grid2D object."""
        # ngc is number of ghost cells, ngx and ngy are number of nodes
        # including ghost nodes
        ngx = nx + 2*ngc
        ngy = ny + 2*ngc

        self.nx = nx
        self.ny = ny
        self.ngc = ngc

        self.ngx = ngx
        self.ngy = ngy

        # x and y locations of nodes including ghost nodes
        self.x = np.empty((ngx, ngy))
        self.y = np.empty((ngx, ngy))

        # x and y locations of cell centers including ghost cells
        self.xc = np.empty((ngx-1, ngy-1))
        self.yc = np.empty((ngx-1, ngy-1))

        # x and y locations of cell faces centers along j node lines
        self.xiface = np.empty((ngx, ngy-1))
        self.yiface = np.empty((ngx, ngy-1))

        # x and y locations of cell faces centers along i node lines
        self.xjface = np.empty((ngx-1, ngy))
        self.yjface = np.empty((ngx-1, ngy))

        # cell Jacobian (cell Volume/Area)
        self.JacobianCell = np.empty((ngx-1, ngy-1))

        # inverse cell Jacobian
        self.InvJacobianCell = np.empty((ngx-1, ngy-1))

        # cell grid metrics
        self.dxdzeta_cell = np.empty((ngx-1, ngy-1))
        self.dxdeta_cell = np.empty((ngx-1, ngy-1))
        self.dydzeta_cell = np.empty((ngx-1, ngy-1))
        self.dydeta_cell = np.empty((ngx-1, ngy-1))

        # cell inverse grid metrics
        self.dzetadx_cell = np.empty((ngx-1, ngy-1))
        self.detadx_cell = np.empty((ngx-1, ngy-1))
        self.dzetady_cell = np.empty((ngx-1, ngy-1))
        self.detady_cell = np.empty((ngx-1, ngy-1))

        # grid metrics computed on cell faces normal to the i direction
        self.dxdzeta_iface = np.empty((ngx, ngy-1))
        self.dxdeta_iface = np.empty((ngx, ngy-1))
        self.dydzeta_iface = np.empty((ngx, ngy-1))
        self.dydeta_iface = np.empty((ngx, ngy-1))

        # grid metrics computed on cell faces normal to the j direction
        self.dxdzeta_jface = np.empty((ngx-1, ngy))
        self.dxdeta_jface = np.empty((ngx-1, ngy))
        self.dydzeta_jface = np.empty((ngx-1, ngy))
        self.dydeta_jface = np.empty((ngx-1, ngy))


def generategrid_2d_uniform(Lx, Ly, nx, ny, ngc):
    """Compute uniform 2D grid from inputs."""
    # nx and ny are number of nodes in the zeta and eta computation coordinates
    # ngc is the number of ghost cell to be added to the grid
    # Lx and Ly are the x and y lenghts of the grid

    # Initialize grid with Grid2D class
    grid = Grid2D(nx, nx, ngc)

    # Compute the uniform steps in x and y directions
    dx = Lx/(nx-1)
    dy = Ly/(ny-1)

    for i in xrange(nx):
        for j in xrange(ny):
            grid.x[i+ngc, j+ngc] = i*dx
            grid.y[i+ngc, j+ngc] = j*dy

    grid.x, grid.y = compute2dghostnodes(grid.x, grid.y, grid.ngc)

    grid.xc, grid.yc = compute2dcellcenters(grid.x, grid.y)

    grid.xiface, grid.yiface, grid.xjface, grid.yjface = compute2dfacecenters(
                                                         grid.x, grid.y)

    grid.dxdzeta_cell, grid.dxdeta_cell, grid.dydzeta_cell, grid.dydeta_cell\
        = compute2dcellgridmetrics(grid.xiface, grid.yiface,
                                   grid.xjface, grid.yjface)

    grid.dxdzeta_iface, grid.dxdeta_iface, grid.dydzeta_iface, \
        grid.dydeta_iface, grid.dxdzeta_jface, grid.dxdeta_jface, \
        grid.dydzeta_jface, grid.dydeta_jface = compute2dfacegridmetrics(
            grid.x, grid.y, grid.xiface, grid.yiface, grid.xjface, grid.yjface)

    grid.JacobianCell = compute2djacobian(grid.dxdzeta_cell, grid.dxdeta_cell,
                                          grid.dydzeta_cell, grid.dydeta_cell)

    grid.InvJacobianCell = compute2dinversejacobian(grid.JacobianCell)

    grid.dzetadx_cell, grid.dzetady_cell, grid.detadx_cell, grid.detady_cell\
        = compute2dcellinvmetrics(grid.InvJacobianCell, grid.dxdzeta_cell,
                                  grid.dxdeta_cell, grid.dydzeta_cell,
                                  grid.dydeta_cell)

    return grid


def compute2dcellinvmetrics(invjacobian, dxdzeta, dxdeta, dydzeta, dydeta):
    """Compute inverse grid metrics at cell center for 2D grid."""
    # compute inverse grid metrics using jacobian
    # and grid metrics computed at the cell centers
    dzetadx = invjacobian*dydeta
    dzetady = -1.0*invjacobian*dxdeta
    detadx = -1.0*invjacobian*dydzeta
    detady = invjacobian*dxdzeta
    return dzetadx, dzetady, detadx, detady


def compute2dinversejacobian(jacobian):
    """Compute inverse jacobian from jacobian."""
    # stole 1/J for later uses
    invjacobian = 1.0/jacobian
    return invjacobian


def compute2djacobian(dxdzeta_cell, dxdeta_cell, dydzeta_cell, dydeta_cell):
    """Compute jacobian for 2D grid."""
    # inputs are the grid metrics computed at the cell centers using the
    # cell face centers

    # compute jacobian at all cell centers
    jacobian = dxdzeta_cell*dydeta_cell - dxdeta_cell*dydzeta_cell

    return jacobian


def compute2dfacegridmetrics(x, y, xiface, yiface, xjface, yjface):
    """Compute grid metrics at cell faces for 2D grid."""
    # x and y are coordinates for the grid nodes
    nx = len(x[:, 0])  # number of points in the i direction
    ny = len(x[0, :])  # number of points in the j direction

    # grid metrics i normal cell faces
    dxdzeta_iface = np.empty((nx, ny-1))
    dxdeta_iface = np.empty((nx, ny-1))
    dydzeta_iface = np.empty((nx, ny-1))
    dydeta_iface = np.empty((nx, ny-1))

    # grid metrics j normal cell faces
    dxdzeta_jface = np.empty((nx-1, ny))
    dxdeta_jface = np.empty((nx-1, ny))
    dydzeta_jface = np.empty((nx-1, ny))
    dydeta_jface = np.empty((nx-1, ny))

    # compute derivatives on the ifaces
    for j in xrange(ny-1):
        for i in xrange(nx):
            # second order accurate derivative w.r.t. zeta
            if i == 0:
                # second order finite forward difference
                dxdzeta_iface[i, j] = 0.5*(-1.0*xiface[i+2, j]
                                           + 4.0*xiface[i+1, j]
                                           - 3.0*xiface[i, j])
                dydzeta_iface[i, j] = 0.5*(-1.0*yiface[i+2, j]
                                           + 4.0*yiface[i+1, j]
                                           - 3.0*yiface[i, j])
            elif i == nx-1:
                # second order finite backward difference
                dxdzeta_iface[i, j] = 0.5*(3.0*xiface[i, j]
                                           - 4.0*xiface[i-1, j]
                                           + xiface[i-2, j])
                dydzeta_iface[i, j] = 0.5*(3.0*yiface[i, j]
                                           - 4.0*yiface[i-1, j]
                                           + yiface[i-2, j])
            else:
                # second order central difference
                dxdzeta_iface[i, j] = 0.5*(xiface[i+1, j] - xiface[i-1, j])
                dydzeta_iface[i, j] = 0.5*(yiface[i+1, j] - yiface[i-1, j])
            # second order accurate derivative w.r.t. eta
            if j == 0:
                # second order finite forward difference
                dxdeta_iface[i, j] = 0.5*(-1.0*xiface[i, j+2]
                                          + 4.0*xiface[i, j+1]
                                          - 3.0*xiface[i, j])
                dydeta_iface[i, j] = 0.5*(-1.0*yiface[i, j+2]
                                          + 4.0*yiface[i, j+1]
                                          - 3.0*yiface[i, j])
            elif j == ny-2:
                # second order finite backward difference
                dxdeta_iface[i, j] = 0.5*(3.0*xiface[i, j]
                                          - 4.0*xiface[i, j-1]
                                          + xiface[i, j-2])
                dydeta_iface[i, j] = 0.5*(3.0*yiface[i, j]
                                          - 4.0*yiface[i, j-1]
                                          + yiface[i, j-2])
            else:
                # second order central difference
                dxdeta_iface[i, j] = 0.5*(xiface[i, j+1] - xiface[i, j-1])
                dydeta_iface[i, j] = 0.5*(yiface[i, j+1] - yiface[i, j-1])

    # compute derivatives on the jfaces
    for j in xrange(ny):
        for i in xrange(nx-1):
            # second order accurate derivative w.r.t. zeta
            if i == 0:
                # second order finite forward difference
                dxdzeta_jface[i, j] = 0.5*(-1.0*xjface[i+2, j]
                                           + 4.0*xjface[i+1, j]
                                           - 3.0*xjface[i, j])
                dydzeta_jface[i, j] = 0.5*(-1.0*yjface[i+2, j]
                                           + 4.0*yjface[i+1, j]
                                           - 3.0*yjface[i, j])
            elif i == nx-2:
                # second order finite backward difference
                dxdzeta_jface[i, j] = 0.5*(3.0*xjface[i, j]
                                           - 4.0*xjface[i-1, j]
                                           + xjface[i-2, j])
                dydzeta_jface[i, j] = 0.5*(3.0*yjface[i, j]
                                           - 4.0*yjface[i-1, j]
                                           + yjface[i-2, j])
            else:
                # second order central difference
                dxdzeta_jface[i, j] = 0.5*(xjface[i+1, j] - xjface[i-1, j])
                dydzeta_jface[i, j] = 0.5*(yjface[i+1, j] - yjface[i-1, j])
            # second order accurate derivative w.r.t. eta
            if j == 0:
                # second order finite forward difference
                dxdeta_jface[i, j] = 0.5*(-1.0*xjface[i, j+2]
                                          + 4.0*xjface[i, j+1]
                                          - 3.0*xjface[i, j])
                dydeta_jface[i, j] = 0.5*(-1.0*yjface[i, j+2]
                                          + 4.0*yjface[i, j+1]
                                          - 3.0*yjface[i, j])
            elif j == ny-1:
                # second order finite backward difference
                dxdeta_jface[i, j] = 0.5*(3.0*xjface[i, j]
                                          - 4.0*xjface[i, j-1]
                                          + xjface[i, j-2])
                dydeta_jface[i, j] = 0.5*(3.0*yjface[i, j]
                                          - 4.0*yjface[i, j-1]
                                          + yjface[i, j-2])
            else:
                # second order central difference
                dxdeta_jface[i, j] = 0.5*(xjface[i, j+1] - xjface[i, j-1])
                dydeta_jface[i, j] = 0.5*(yjface[i, j+1] - yjface[i, j-1])
    return dxdzeta_iface, dxdeta_iface, dydzeta_iface, dydeta_iface, \
        dxdzeta_jface, dxdeta_jface, dydzeta_jface, dydeta_jface


def compute2dcellgridmetrics(xiface, yiface, xjface, yjface):
    """Compute grid metrics at cell center for 2D grid."""
    # xiface and yiface are coordinates of the cell faces along j node lines
    # xjface and yjface are coordinates of the cell faces along i node lines
    nx = len(xjface[:, 0])  # number of points in the i direction
    ny = len(xiface[0, :])  # number of points in the j direction

    # x direction grid metrics computed at the cell centers using face values
    dxdzeta_cell = np.empty((nx, ny))
    dxdeta_cell = np.empty((nx, ny))
    # y direction grid metrics computed at the cell centers using face values
    dydzeta_cell = np.empty((nx, ny))
    dydeta_cell = np.empty((nx, ny))
    for j in xrange(ny):
        for i in xrange(nx):
            # all dzeta and deta are equal to 1 by definition
            dxdzeta_cell[i, j] = xiface[i+1, j] - xiface[i, j]
            dydzeta_cell[i, j] = yiface[i+1, j] - yiface[i, j]
            dxdeta_cell[i, j] = xjface[i, j+1] - xjface[i, j]
            dydeta_cell[i, j] = yjface[i, j+1] - yjface[i, j]

    return dxdzeta_cell, dxdeta_cell, dydzeta_cell, dydeta_cell


def compute2dfacecenters(x, y):
    """Compute cell face centers for 2D grid."""
    # nx and ny are number of nodes in the x and y directions
    # including ghost nodes
    nx = len(x[:, 0])
    ny = len(x[0, :])
    # xiface and yiface are the x and y coordinates of the face centers
    # along constant j node lines
    xiface = np.empty((nx, ny-1))
    yiface = np.empty((nx, ny-1))
    # xjface and yjface are the x and y coordinates of the face centers
    # along constant i node lines
    xjface = np.empty((nx-1, ny))
    yjface = np.empty((nx-1, ny))

    for j in xrange(ny-1):
        for i in xrange(ny):
            xiface[i, j] = 0.5*(x[i, j] + x[i, j+1])
            yiface[i, j] = 0.5*(y[i, j] + y[i, j+1])
    for j in xrange(ny):
        for i in range(nx-1):
            xjface[i, j] = 0.5*(x[i, j] + x[i+1, j])
            yjface[i, j] = 0.5*(y[i, j] + y[i+1, j])

    return xiface, yiface, xjface, yjface


def compute2dghostnodes(x, y, ngc):
    """Compute ghost node locations for 2D grid."""
    # x and y are the node coordinates in the x and y directions
    # x and y contain indices for the ghost cell
    # ngc is the number of ghost cell layers

    # loop will work from inner most ghost cells to outer most ghost cells
    for offset in xrange(ngc):
        buf = ngc - offset
        # compute outer interior cell normals normals
        boundary = ['imin', 'imax', 'jmin', 'jmax']
        # loop through all block bounaries
        for loc in boundary:
            if loc == 'imin':
                # get x and y values along imin bourder
                xb = x[buf, buf:-buf]
                yb = y[buf, buf:-buf]
                # set norm signs so that normal is facing outward
                sx = 1.0
                sy = 1.0
                # compute the dx and dy between bourder nodes and
                # first interior nodes
                dx = x[buf, buf:-buf] - x[buf+1, buf:-buf]
                dy = y[buf, buf:-buf] - y[buf+1, buf:-buf]
            elif loc == 'imax':
                # get x and y values along imax bourder
                xb = x[-(buf+1), buf:-buf]
                yb = y[-(buf+1), buf:-buf]
                # set norm signs so that normal is facing outward
                sx = 1.0
                sy = 1.0
                # compute the dx and dy between bourder nodes and
                # first interior nodes
                dx = x[-(buf+1), buf:-buf] - x[-(buf+2), buf:-buf]
                dy = y[-(buf+1), buf:-buf] - y[-(buf+2), buf:-buf]
            elif loc == 'jmin':
                # get x and y values along jmin bourder
                xb = x[buf:-buf, buf]
                yb = y[buf:-buf, buf]
                # set norm signs so that normal is facing outward
                sx = 1.0
                sy = 1.0
                # compute the dx and dy between bourder nodes and
                # first interior nodes
                dx = x[buf:-buf, buf] - x[buf:-buf, buf+1]
                dy = y[buf:-buf, buf] - y[buf:-buf, buf+1]
            elif loc == 'jmax':
                # get x and y values along jmax bourder
                xb = x[buf:-buf, -(buf+1)]
                yb = y[buf:-buf, -(buf+1)]
                # set norm signs so that normal is facing outward
                sx = 1.0
                sy = 1.0
                # compute the dx and dy between bourder nodes and
                # first interior nodes
                dx = x[buf:-buf, -(buf+1)] - x[buf:-buf, -(buf+2)]
                dy = y[buf:-buf, -(buf+1)] - y[buf:-buf, -(buf+2)]
            # compute face normals using boundary node values
            facenormx, facenormy = compute2dfacenormals(xb, yb, sx, sy)
            # compute node normals at boundary using face normals
            nodenormx, nodenormy = compute2dnodenormals(facenormx, facenormy)
            # compute ghost node x and y coordinates
            if loc == 'imin':
                x[buf-1, buf:-buf] = x[buf+1, buf:-buf] + 2.0*dx*nodenormx
                y[buf-1, buf:-buf] = y[buf+1, buf:-buf] + 2.0*dy*nodenormy
            if loc == 'imax':
                x[-buf, buf:-buf] = x[-(buf+2), buf:-buf] + 2.0*dx*nodenormx
                y[-buf, buf:-buf] = y[-(buf+2), buf:-buf] + 2.0*dy*nodenormy
            if loc == 'jmin':
                x[buf:-buf, buf-1] = x[buf:-buf, buf+1] + 2.0*dx*nodenormx
                y[buf:-buf, buf-1] = y[buf:-buf, buf+1] + 2.0*dy*nodenormy
            if loc == 'jmax':
                x[buf:-buf, -buf] = x[buf:-buf, -(buf+2)] + 2.0*dx*nodenormx
                y[buf:-buf, -buf] = y[buf:-buf, -(buf+2)] + 2.0*dy*nodenormy

        # compute corner ghost node x and y coordinates
        # imin,jmin corner
        x[buf-1, buf-1] = x[buf-1, buf]
        y[buf-1, buf-1] = y[buf, buf-1]
        # imin,jmax corner
        x[buf-1, -buf] = x[buf-1, -(buf+1)]
        y[buf-1, -buf] = y[buf, -buf]
        # imax,jmin corner
        x[-buf, buf-1] = x[-buf, buf]
        y[-buf, buf-1] = y[-(buf+1), buf-1]
        # imax,jmax corner
        x[-buf, -buf] = x[-buf, -(buf+1)]
        y[-buf, -buf] = y[-(buf+1), -buf]
    return x, y


def compute2dfacenormals(x, y, sx, sy):
    """Compute cell face normal vectors for 2D grid."""
    # initialize facenormal x and y components
    normx = np.ones(len(x)-1)
    normy = np.ones(len(x)-1)
    # step though boundary nodes to compute face normal components
    for i in xrange(len(x)-1):
        # x and y tangent components at cell face
        Tx = x[i+1]-x[i]
        Ty = y[i+1]-y[i]
        # magnitude of face tangent vector
        Tmag = np.sqrt(Tx**2.0 + Ty**2.0)
        # compute normal component with signs adjust so normal point outward
        normx[i] = sx*Ty/Tmag
        normy[i] = sy*Tx/Tmag
    return normx, normy


def compute2dnodenormals(fnx, fny):
    """Compute node normal vectors for 2D grid."""
    # fnx and fny are face normal x and y components
    # initialize node normal components
    normx = np.ones(len(fnx)+1)
    normy = np.ones(len(fny)+1)
    # set edge node normals equal to inner face normals
    normx[0] = fnx[0]
    normy[0] = fny[0]
    normx[-1] = fnx[-1]
    normy[-1] = fny[-1]
    # set node normals equal to the mean of neighboring cell normals
    for i in range(1, len(normx)-1):
        normx[i] = 0.5*(fnx[i-1] + fnx[i])
        normy[i] = 0.5*(fny[i-1] + fny[i])
    return normx, normy


def compute2dcellcenters(x, y):
    """Compute cell centers for 2D grid."""
    # x and y are the coordinates at the nodes
    # nx and ny and the number of nodes in the x and y directions
    nx = len(x[:, 0])
    ny = len(y[0, :])

    # xc and yc are the x and y values of the cell centers
    xc = np.empty((nx-1, ny-1))
    yc = np.empty((nx-1, ny-1))

    # To compute cell center of quadralateral we will compute the centriod
    # of triangle ABC and triangle BCD and then use the area weighted average
    # to compute the centroid of the quadralateral
    # C o----------o D
    #   |          |
    #   |    *Pc   |
    #   |          |
    # A o----------o B

    for j in xrange(ny-1):
        for i in xrange(nx-1):
            # Point coordinates
            ax = x[i, j]
            ay = y[i, j]
            bx = x[i+1, j]
            by = y[i+1, j]
            cx = x[i, j+1]
            cy = y[i, j+1]
            dx = x[i+1, j+1]
            dy = y[i+1, j+1]
            # Area of triangle ABC
            T1 = 0.5*abs(ax*(by-cy) + bx*(cy-ay) + cx*(ay-by))
            # Area of triangle ABC
            T2 = 0.5*abs(bx*(cy-dy) + cx*(dy-by) + dx*(by-cy))
            # Centroid coordinates of triangle ABC
            P1x = (ax+bx+cx)/3.0
            P1y = (ay+by+cy)/3.0
            # Centroid coordinates of triangle BCD
            P2x = (bx+cx+dx)/3.0
            P2y = (by+cy+dy)/3.0
            # Centroid of quadralateral
            xc[i, j] = (T1*P1x + T2*P2x)/(T1+T2)
            yc[i, j] = (T1*P1y + T2*P2y)/(T1+T2)
    return xc, yc
