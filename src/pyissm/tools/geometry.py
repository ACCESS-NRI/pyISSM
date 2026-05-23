"""
Geometry-related functions for ISSM

This module contains functions to compute geometric properties on ISSM model meshes.
"""

import warnings
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


def effectivepressure(md, head=None):
    """
    Calculate the effective basal pressure N from md.geometry and the
    effective pressure coupling rule in md.friction.

    Parameters
    ----------
    md : object
        ISSM model object.
    head : array-like, optional
        Hydraulic head at mesh vertices (m). Only used when
        ``md.friction.coupling == 4``.

    Returns
    -------
    N : ndarray
        Effective pressure at the base (Pa).

    Raises
    ------
    ValueError
        If ``coupling == 4`` and ``md.hydrology`` is not a
        ``hydrologyprescribe`` instance and no ``head`` is provided.
    ValueError
        If an unsupported coupling value is given.

    Notes
    -----
    Coupling modes:

    * 0 — Uniform sheet; negative water pressure permitted (default).
    * 1 — Effective pressure equals overburden pressure.
    * 2 — Uniform sheet; water pressure clamped to >= 0.
    * 3 — Prescribed effective pressure from ``md.friction.effective_pressure``.
    * 4 — Dynamically coupled to the hydrology model (``md.hydrology.head``
          or a supplied ``head`` array).

    See also
    --------
    basalstress
    """
    from pyissm.model.classes.hydrology import shreve as hydrologyprescribe  # closest prescribe analogue

    if hasattr(md.friction, 'coupling'):
        coupling = md.friction.coupling
    else:
        warnings.warn('md.friction.coupling not found; defaulting to coupling=0.')
        coupling = 0

    sealevel = 0.0  # sea level reference elevation (m)

    p_ice   = md.constants.g * md.materials.rho_ice   * md.geometry.thickness
    p_water = md.constants.g * md.materials.rho_water * (sealevel - md.geometry.base)

    if coupling == 0:
        N = p_ice - p_water
    elif coupling == 1:
        N = p_ice
    elif coupling == 2:
        p_water = np.maximum(p_water, 0.0)
        N = p_ice - p_water
    elif coupling == 3:
        N = np.maximum(md.friction.effective_pressure, 0.0)
    elif coupling == 4:
        if head is not None:
            head = np.asarray(head)
            if head.size != md.mesh.numberofvertices:
                raise ValueError(
                    f'head size ({head.size}) does not match numberofvertices '
                    f'({md.mesh.numberofvertices}).'
                )
            p_water = md.constants.g * md.materials.rho_freshwater * (head - md.geometry.base)
        else:
            if hasattr(md.hydrology, 'head'):
                p_water = md.constants.g * md.materials.rho_freshwater * (md.hydrology.head - md.geometry.base)
            else:
                raise ValueError(
                    f'coupling=4 is not supported for {type(md.hydrology).__name__}. '
                    'Provide a head array or use a hydrology class with a head field.'
                )
        N = p_ice - p_water
    else:
        raise ValueError(f'Unsupported coupling value: {coupling}')

    return N

