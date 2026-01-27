"""
Geometry-related functions for ISSM

This module contains functions to compute geometric properties on ISSM model meshes.
"""

import numpy as np
from pyissm import model


def slope(md, field = None):
    """
    Compute the slope of a given field over the mesh elements.

    Parameters
    ----------
    md : object
        Model object containing mesh and geometry information.
    field : array-like, optional
        Field values at mesh nodes. If not provided, uses the surface elevation
        from md.geometry.surface. Default is None.

    Returns
    -------
    sx : ndarray
        Slope component in the x-direction at element centers.
    sy : ndarray
        Slope component in the y-direction at element centers.
    s : ndarray
        Magnitude of the slope at element centers.

    Raises
    ------
    NotImplementedError
        If the mesh dimension is 3D, as 3D meshes are not yet supported.

    Notes
    -----
    The slope is computed using nodal basis functions N(x, y) = alpha*x + beta*y + gamma.
    For 2D meshes, elements and coordinates are taken directly from the mesh.
    For 3D meshes, 2D elements and coordinates are used instead.

    Examples
    --------
    >>> sx, sy, s = slope(md)
    >>> sx, sy, s = slope(md, field=md.geometry.bed)
    """

    # Condtionally assign mesh variables based on dimension
    if md.mesh.dimension() == 2:
        elements = md.mesh.elements
        x = md.mesh.x
        y = md.mesh.y
    else:
        elements = md.mesh.elements2d
        x = md.mesh.x2d
        y = md.mesh.y2d

    # If field is not provided, use surface
    if field is None:
        field = md.geometry.surface

    # Compute nodal functions coefficients N(x, y) = alpha x + beta y + gamma
    alpha, beta, _ = model.mesh.get_nodal_functions_coeff(elements, x, y)

    # Compute slope components at element centers
    ones = np.ones((3, 1))
    sx = np.dot(field[elements - 1] * alpha, ones).ravel()
    sy = np.dot(field[elements - 1] * beta, ones).ravel()
    s = np.sqrt(sx**2 + sy**2)

    # Project to 3D if necessary
    if md.mesh.dimension() == 3:
        raise NotImplementedError("pyissm.tools.geometry.slope: 3D meshes not supported yet.")
    
    return sx, sy, s

def nowicki_profile(x):
    """
    Create a theoretical ice profile at the transition zone
    based on Sophie Nowicki's thesis.

    Parameters
    ----------
    x : array_like
        Along-flow coordinate

    Returns
    -------
    b : ndarray
        Ice base
    h : ndarray
        Ice thickness
    sea : float
        Sea level
    """
    x = np.asarray(x)
    n = len(x)
    mid = n // 2

    # Physical constants
    delta = 0.1 # (rho_w / rho_i) - 1
    hg = 1.0 # ice thickness at grounding line
    lam = 0.1 # deviatoric stress / water pressure
    beta = 5.0 # friction coefficient
    ms = 0.005 # surface accumulation rate
    q = 0.801 # ice mass flux

    sea = hg / (1 + delta)

    # Allocate arrays
    b = np.zeros(n)
    h = np.zeros(n)

    # Upstream of grounding line
    for i in range(mid):
        coeffs = [1,
                  4 * lam * beta,
                  0,
                  0,
                  6 * lam * ms * x[i]**2 + 12 * lam * q * x[i] - hg**4 - 4 * lam * beta * hg**3]

        roots = np.roots(coeffs)
        s = roots[(roots.imag == 0) & (roots.real > 0)].real[0]

        h[i] = s
        b[i] = 0.0

    # Downstream of grounding line
    xd = x[mid:]
    h[mid:] = (xd / (4 * (delta + 1) * q) + hg**-2) ** -0.5
    b[mid:] = sea - h[mid:] / (1 + delta)

    return b, h, sea

