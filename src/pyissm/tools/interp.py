"""
Interpolation functions for ISSM

This module contains various interpolation functions that can be used in conjunction with ISSM models.
"""

import numpy as np
from scipy.sparse import csc_matrix
from pyissm import model


def averaging(md, data, iterations, layer = 0):
    """
    Smooth input data over a mesh using area/volume-weighted nodal averaging.

    This function performs an iterative smoothing of data defined either on mesh
    elements or on mesh vertices. For element-defined data, it first distributes
    element values to nodes using element areas/volumes as weights. For vertex-
    defined data, it uses the provided nodal values as the starting point. Each
    iteration computes element-averaged values from current nodal values and then
    recomputes nodal values as the area/volume-weighted average of adjacent
    element averages. The algorithm supports 2-D meshes, full 3-D meshes (using
    element volumes), and extraction of a single 2-D layer from a 3-D mesh.
    
    Parameters
    ----------
    md : ISSM model object
    data : array_like
        1-D array-like of scalar values defined either on elements or on
        vertices. Its length must be equal to either:
        - md.mesh.numberofelements  (data on elements), or
        - md.mesh.numberofvertices  (data on vertices).
        If `layer` is provided for a 3-D mesh, the vertex-based length should be
        md.mesh.numberofvertices2d and element-based length md.mesh.numberofelements2d.
    iterations : int
        Number of smoothing iterations to perform. Must be a non-negative
        integer. A value of 0 will still perform the initial element->node
        distribution if `data` is provided on elements.
    layer : int, optional
        Layer index to extract when working with a 3-D mesh. Default is 0,
        meaning operate on the full mesh. If non-zero, it must satisfy
        1 <= layer <= md.mesh.numberoflayers and the function will operate on
        the corresponding 2-D layer (using mesh.elements2d, x2d, y2d, etc.).

    Returns
    -------
    average : numpy.ndarray, shape (n, 1)
        Dense column vector containing the smoothed nodal values. The length n is
        md.mesh.numberofvertices when `layer == 0`, or md.mesh.numberofvertices2d
        when a 2-D layer is extracted.

    Raises
    ------
    Exception
        If the length of `data` does not match the expected number of elements or
        vertices for the mesh/layer combination.
    ValueError
        If a non-zero `layer` is provided for a 3-D mesh but is outside the valid
        range (not in 1..md.mesh.numberoflayers).

    Notes
    -----
    - Weighting: element areas (2-D) or volumes (3-D) are used to weight the
      contribution of each element to its corner nodes when distributing element
      values to nodes and when recomputing nodal values from element averages.
    - Indexing: mesh connectivity arrays are expected to be 1-based in the model
      and are converted to 0-based indices internally.
    - Output is returned as a dense array (not a sparse matrix) for compatibility
      with downstream code (e.g. C routines) that expect full arrays.

    Examples
    --------
    >>> # Smooth element data on a 2-D mesh for 5 iterations
    >>> avg = averaging(md, element_data, iterations=5)
    >>> # Smooth nodal data on layer 2 of a 3-D mesh
    >>> avg2d = averaging(md, vertex_data_2d, iterations=3, layer=2)
    """
    ## NOTE: Code taken directly from $ISSM_DIR/src/m/interp/averaging.py

    # Error checks
    if len(data) not in (md.mesh.numberofelements, md.mesh.numberofvertices):
        raise Exception('pyissm.tools.interp.averaging: Data must be defined on elements or vertices.')
    
    if (md.mesh.dimension() == 3) & (layer != 0):
        if (layer <= 0) or (layer > md.mesh.numberoflayers):
            raise ValueError('pyissm.tools.interp.averaging: Layer should be between 1 and md.mesh.numberoflayers for 3D meshes.')
    else:
        layer = 0

    # Initialization
    if layer == 0:
        weights = np.zeros((md.mesh.numberofvertices, ))
        data = np.asarray(data).flatten()
    else:
        weights = np.zeros((md.mesh.numberofvertices2d, ))
        # Extract the specified layer from data
        data = data[(layer - 1) * md.mesh.numberofvertices2d + 1:layer * md.mesh.numberofvertices2d - 1, :]

    # Define mesh variables
    if layer == 0:
        index = md.mesh.elements
        numberofnodes = md.mesh.numberofvertices
        numberofelements = md.mesh.numberofelements
    else:
        index = md.mesh.elements2d
        numberofnodes = md.mesh.numberofvertices2d
        numberofelements = md.mesh.numberofelements2d

    # Build some variables
    if (md.mesh.dimension() == 3) & (layer == 0):
        ## 3D mesh, get volumes
        rep = 6
        areas = model.mesh.get_element_areas_volumes(index, md.mesh.x, md.mesh.y, md.mesh.z)
    elif md.mesh.dimension() == 2:
        ## 2D mesh, get areas
        rep = 3
        areas = model.mesh.get_element_areas_volumes(index, md.mesh.x, md.mesh.y)
    else:
        ## 3D mesh with 2D layer extracted, get areas
        rep = 3
        areas = model.mesh.get_element_areas_volumes(index, md.mesh.x2d, md.mesh.y2d)

    index = index - 1  # Python indexes from zero
    line = index.T.flatten()
    areas = np.vstack(areas).reshape(-1, )
    summation = 1. / rep * np.ones((rep,1) )
    linesize = rep * numberofelements

    # Update weights that hold the volume of all the element holding the node i
    weights = csc_matrix((np.tile(areas, (1, rep)).reshape(-1,), (line, np.zeros(linesize, ))), shape=(numberofnodes, 1))

    # Initialization
    if len(data) == numberofelements:
        average_node = csc_matrix((np.tile(np.multiply(areas,data), (1, rep)).reshape(-1, ), (line, np.zeros(linesize, ))), shape=(numberofnodes, 1))
        average_node = np.divide(average_node,weights)
        average_node = csc_matrix(average_node)
    else:
        average_node = csc_matrix(data.reshape(-1, 1))

    # Loop over iteration
    for i in np.arange(1, iterations + 1):
        average_el = np.asarray(average_node.todense()[index].reshape(numberofelements, rep)*summation).reshape(-1, )
        average_node = csc_matrix((np.tile(np.multiply(areas,average_el.reshape(-1)), (1, rep)).reshape(-1, ), (line, np.zeros(linesize, ))), shape=(numberofnodes, 1))
        average_node = np.divide(average_node,weights)
        average_node = csc_matrix(average_node)

    # Return output as a full matrix (C code does not like sparse matrices)
    average = np.expand_dims(np.asarray(average_node.todense()).reshape(-1, ),axis=1)

    return average