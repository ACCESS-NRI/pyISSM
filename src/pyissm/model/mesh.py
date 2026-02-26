"""

Functions for building and interacting with an ISSM model mesh.

"""
import numpy as np
import os
import collections
import matplotlib.tri as tri
from scipy.interpolate import griddata
import warnings
from pyissm import model, tools

def get_mesh(mesh_x,
             mesh_y,
             mesh_elements):
    """
    Create a triangular mesh from an unstructured model object.

    Extracts node coordinates and element connectivity from the model and
    constructs a `matplotlib.tri.Triangulation` object for downstream
    operations.

    Parameters
    ----------
    mesh_x : 1d array
        x-coordinates of the mesh nodes.
    mesh_y : 1d array
        y-coordinates of the mesh nodes.
    mesh_elements : 2d array
        element connectivity.

    Returns
    -------
    mesh : matplotlib.tri.Triangulation
        Triangular mesh object representing the model domain.

    Notes
    -----
    - If necessary, the function adjusts the element indexing from 1-based to 0-based indexing
     required by `Triangulation`.

    See Also
    --------
    matplotlib.tri.Triangulation : Triangular mesh representation.
    make_gridded_domain_mask : Uses this mesh to determine in-domain points.
    grid_model_field : Interpolates data onto a regular grid using the mesh structure.
    """

    ## Check mesh_elements uses 0-based indexing & is defined correctly
    ## -------------------------------------
    if mesh_elements.min() == 1:
        mesh_elements = mesh_elements - 1

    ## Create triangulation feature
    mesh = tri.Triangulation(mesh_x, mesh_y, mesh_elements)

    return mesh

def process_mesh(md):
    """
    Process ISSM model mesh

    This function processes key elements of an ISSM model mesh and is
    used in several core pyISSM functions to ensure consistency

    Parameters
    ----------
    md : ISSM Model object
        ISSM Model object from which the mesh should be extracted/processed.

    Returns
    -------
    mesh : matplotlib.tri.Triangulation
        Triangular mesh object representing the model domai
    mesh_x : 1d array
        x-coordinates of the mesh nodes.
    mesh_y : 1d array
        y-coordinates of the mesh nodes.
    mesh_elements : 2d array
        element connectivity.
    is3d : bool
        'True' if elements2d exists (and the model is 3D), 'False' otherwise

    Example
    -------
    mesh, mesh_x, mesh_y, mesh_elements, is3d = process_mesh(md)
    """

    ## Set default values
    is3d = False

    ## Process a 3D model
    if tools.general.has_nested_attr(md, 'mesh', 'elements2d'):

        # Create mesh object
        mesh = get_mesh(md.mesh.x2d, md.mesh.y2d, md.mesh.elements2d)

        # Extract X/Y Coordinates for 2D mesh
        mesh_x = md.mesh.x2d
        mesh_y = md.mesh.y2d

        # Extract and adjust 2D element numbering to be 0-indexed
        mesh_elements = md.mesh.elements2d - 1

        # Return is3d and display warning for 3D mesh
        is3d = True
        warnings.warn('process_mesh: 3D model found. Processing as 2D mesh.')

    ## Process a 2D model
    else:

        # Create mesh object
        mesh = get_mesh(md.mesh.x, md.mesh.y, md.mesh.elements)

        # Extract X/Y Coordinates
        mesh_x = md.mesh.x
        mesh_y = md.mesh.y

        # Extract and adjust element numbering to be 0-indexed
        mesh_elements = md.mesh.elements - 1

    return mesh, mesh_x, mesh_y, mesh_elements, is3d

def find_node_types(md,
                    ice_levelset,
                    ocean_levelset):
    """
    Identify node types (ice, ice-front, ocean, floating ice, grounded ice) from level set data.

    This function processes level set fields for ice and ocean and classifies mesh nodes into
    several categories based on their sign.

    For 3D meshes, only the surface layer (where `vertexonsurface == 1`) is processed.

    Parameters
    ----------
    md : ISSM Model object
        ISSM Model object containing the mesh and geometry information. Must have attributes:
        'md.mesh.*' used by process_mesh().
    ice_levelset : ndarray
        1D array of values from the ice level set:
            Ice < 0
            Ice front = 0
            No Ice > 0
    ocean_levelset : ndarray
        1D array of values from the ocean level set:
            Ocean < 0
            No ocean > 0

    Returns
    -------
    dict of str -> ndarray
        Dictionary with boolean arrays (same length as number of surface nodes), indicating:

        - 'ice_nodes' : Nodes with ice (ice_levelset < 0)
        - 'ice_front_nodes' : Nodes on the ice front (ice_levelset == 0)
        - 'ocean_nodes' : Nodes with ocean (ocean_levelset < 0)
        - 'floating_ice_nodes' : Nodes with floating ice (ice_levelset < 0 & ocean_levelset < 0)
        - 'grounded_ice_nodes' : Nodes with grounded ice (ice_levelset < 0 & ocean_levelset >= 0)

    Warnings
    --------
    If a 3D mesh is detected, only the surface layer is used. A warning is issued.

    Example
    --------
    model_node_types = find_node_types(md, md.mask.ice_levelset, md.mask.ocean_levelset)
    model_node_types = find_node_types(md, md.results.TransientSolution.MaskIceLevelset[34], md.results.TransientSolution.MaskOceanLevelset[34])
    """

    ## Process model mesh
    mesh, mesh_x, mesh_y, mesh_elements, is3d = process_mesh(md)

    ## Identify ice/ocean nodes from surface layer of 3D model
    if is3d:
        ice_nodes = ice_levelset[md.mesh.vertexonsurface == 1] < 0
        ice_front_nodes = ice_levelset[md.mesh.vertexonsurface == 1] == 0
        ocean_nodes = ocean_levelset[md.mesh.vertexonsurface == 1] < 0

        warnings.warn('find_node_types: 3D model found. Processing surface layer only.')

    ## Identify ice/ocean nodes
    else:
        ice_nodes = ice_levelset < 0
        ice_front_nodes = ice_levelset == 0
        ocean_nodes = ocean_levelset < 0

    ## Identify floating and grounded ice nodes
    floating_ice_nodes = ice_nodes & ocean_nodes
    grounded_ice_nodes = ice_nodes  & ~ocean_nodes

    ## Compile dictionary of node types to return
    output_dict = {
        'ice_nodes': ice_nodes,
        'ice_front_nodes': ice_front_nodes,
        'ocean_nodes': ocean_nodes,
        'floating_ice_nodes': floating_ice_nodes,
        'grounded_ice_nodes': grounded_ice_nodes
    }

    return output_dict


def find_element_types(md,
                       ice_levelset,
                       ocean_levelset):
    """
    Identify element types (ice, ice-front, ocean, floating ice, grounded ice, grounding line) from level set data.

    This function processes level set fields for ice and ocean and classifies mesh elements into
    several categories based on their sign.

    For 3D meshes, only the surface layer (where `vertexonsurface == 1`) is processed (see find_node_types()).

    Parameters
    ----------
    md : ISSM Model object
        ISSM Model object containing the mesh and geometry information. Must have attributes:
        'md.mesh.*' used by process_mesh().
    ice_levelset : ndarray
        1D array of values from the ice level set:
            Ice < 0
            Ice front = 0
            No Ice > 0
    ocean_levelset : ndarray
        1D array of values from the ocean level set:
            Ocean < 0
            No ocean > 0

    Returns
    -------
    dict of str -> ndarray
        Dictionary with boolean arrays (same length as number of surface nodes), indicating:

        - 'ice_elements' : Elements with ice
        - 'ice_front_elements' : Elements on the ice front
        - 'ocean_elements' : Elements with ocean
        - 'floating_ice_elements' : Elements with floating ice
        - 'grounded_ice_elements' : Elements with grounded ice
        - 'grounding_line_elements' : Elements on the grounding line

    Warnings
    --------
    If a 3D mesh is detected, only the surface layer is used. A warning is issued.

    Example
    --------
    model_element_types = find_element_types(md, md.mask.ice_levelset, md.mask.ocean_levelset)
    model_element_types = find_element_types(md, md.results.TransientSolution.MaskIceLevelset[34], md.results.TransientSolution.MaskOceanLevelset[34])
    """

    ## Process model mesh
    mesh, mesh_x, mesh_y, mesh_elements, is3d = process_mesh(md)

    ## Identify ice/ocean nodes model
    ## NOTE: This accounts for 3D models internally.
    node_types = find_node_types(md,
                                 ice_levelset,
                                 ocean_levelset)

    ## Isolate individual node types
    ice_nodes = node_types['ice_nodes']
    ice_front_nodes = node_types['ice_front_nodes']
    ocean_nodes = node_types['ocean_nodes']
    floating_ice_nodes = node_types['floating_ice_nodes']
    grounded_ice_nodes = node_types['grounded_ice_nodes']

    ## Identify "no-ice" nodes
    no_ice_nodes = ~ice_nodes

    ## Identify required element types
    ice_elements = np.sum(ice_nodes[mesh_elements], axis = 1)
    ocean_elements = np.sum(ocean_nodes[mesh_elements], axis = 1)
    no_ice_elements = np.sum(no_ice_nodes[mesh_elements], axis = 1)
    zero_ice_elements = np.sum(ice_front_nodes[mesh_elements], axis = 1)
    floating_ice_elements = np.sum(floating_ice_nodes[mesh_elements], axis=1)
    grounded_ice_elements = np.sum(grounded_ice_nodes[mesh_elements], axis=1)

    ## Identify custom elements types
    ice_front_elements = (ice_elements.astype(bool) & no_ice_elements.astype(bool)) & ~((ice_elements == 2) & zero_ice_elements.astype(bool))
    grounding_line_elements = (ocean_elements != np.max(ocean_elements)) & (ocean_elements != 0)

    ## Compile dictionary of element types to return
    output_dict = {
        'ice_elements': ice_elements,
        'ocean_elements': ocean_elements,
        'floating_ice_elements': floating_ice_elements,
        'grounded_ice_elements': grounded_ice_elements,
        'ice_front_elements': ice_front_elements,
        'grounding_line_elements': grounding_line_elements
    }
    return output_dict

def make_gridded_domain_mask(mesh_x,
                             mesh_y,
                             mesh_elements,
                             grid_x,
                             grid_y):
    """
    Generate a binary domain mask on a regular grid based on an unstructured mesh.

    This function identifies which points in a regular 2D grid fall within an
    unstructured triangular mesh. Points outside the mesh are marked as `False`,
    and those inside are `True`.

    Parameters
    ----------
    mesh_x : ndarray
        1D array of x-coordinates for mesh nodes.
    mesh_y : ndarray
        1D array of y-coordinates for mesh nodes.
    mesh_elements : ndarray
        2D array of element connectivity.
    grid_x : ndarray
        2D array of x-coordinates from `np.meshgrid` defining the regular grid.
    grid_y : ndarray
        2D array of y-coordinates from `np.meshgrid` defining the regular grid.

    Returns
    -------
    domain_mask : ndarray of bool
        Boolean mask array of the same shape as `grid_x` and `grid_y`, where
        `True` indicates that the grid point lies inside the mesh domain,
        and `False` indicates it lies outside.

    Notes
    -----
    - Internally uses `matplotlib.tri.Triangulation` and its `get_trifinder()` method
      to determine point inclusion.

    See Also
    --------
    grid_data : Interpolates data onto a regular grid, using this mask by default.
    """

    ## Check mesh_elements uses 0-based indexing & is defined correctly
    ## -------------------------------------
    if mesh_elements.min() == 1:
        mesh_elements = mesh_elements - 1

    ## Make mesh (triangulation)
    mesh = get_mesh(mesh_x, mesh_y, mesh_elements)

    ## Get XY points of regular grid
    grid_points = np.column_stack((grid_x.ravel(), grid_y.ravel()))

    ## Get ID of each triangle in the mesh. ID values < 0 are outside the mesh
    mask = mesh.get_trifinder()

    ## Identify points that have an ID >= 0 (i.e. are inside the mesh)
    mask = mask(grid_points[:, 0], grid_points[:, 1]) >= 0

    ## Reshape to match XY points
    domain_mask = mask.reshape(grid_x.shape)

    return domain_mask

def grid_model_field(md,
                     model_field,
                     grid_x,
                     grid_y,
                     method = 'linear',
                     domain_mask = None,
                     fill_value = np.nan):
    """
    Interpolate unstructured model data onto a regular 2D grid.

    This function handles both time-varying and static fields defined on either
    mesh nodes or elements. It optionally applies a domain mask to exclude
    values outside the desired region (e.g. outside the mesh or ice-covered areas).

    Parameters
    ----------
    md : object
        A model object containing the unstructured mesh with attributes:
        - `md.mesh.x` (1D array): x-coordinates of mesh nodes.
        - `md.mesh.y` (1D array): y-coordinates of mesh nodes.
        - `md.mesh.elements` (2D array): triangular elements (1-based indexing).
    model_field : ndarray
        The field to be interpolated. Should be either:
        - (npoints,) for static data
        - (nt, npoints) for time-varying data
        Where `npoints` must match the number of mesh nodes or elements.
    grid_x : ndarray
        2D array of x-coordinates from `np.meshgrid` defining the regular output grid.
    grid_y : ndarray
        2D array of y-coordinates from `np.meshgrid` defining the regular output grid.
    method : str, optional
        Interpolation method to use. Options are:
        - 'linear' (default)
        - 'nearest'
        - 'cubic'
    domain_mask : ndarray of bool, optional
        Optional binary mask array (same shape as `grid_x`/`grid_y`) indicating valid
        interpolation region. If not provided, a mask will be automatically generated
        based model mesh. Values where `domain_mask == False` will be set to fill_value.

    fill_value : float, optional
        Value to be used to fill masked area. Default value is np.nan

    Returns
    -------
    gridded_model_field : ndarray
        Interpolated data on the regular grid. Shape is:
        - (ny, nx) for static fields
        - (nt, ny, nx) for time-varying fields
        Invalid/masked regions are set to `np.nan`.

    Raises
    ------
    ValueError
        If the shape of `model_field` does not match number of mesh nodes or elements.
        If a custom `domain_mask` is provided and its shape does not match `grid_x`.

    Notes
    -----
    - Element-based fields are interpolated using element centroids.
    - Time-varying fields are interpolated one time step at a time.

    See Also
    --------
    make_gridded_domain_mask : Generates a domain mask on a regular grid.
    """

    ## Get model_field information
    ## -------------------------------------
    # Does it vary in time, or is it static?
    if model_field.ndim == 1:
        # If static, add a false time dimension (1)
        model_field = model_field[np.newaxis, :]

    # Get model_field dimensions
    nt, npoints = model_field.shape

    # Is it defined on vertices or elements?
    if npoints == md.mesh.numberofvertices:
        # Convert elements to 0-based indexing
        mesh_elements = md.mesh.elements - 1

        # model_field is on vertices. Take native coordinates for interpolation
        mesh_x = md.mesh.x
        mesh_y = md.mesh.y

    elif npoints == md.mesh.numberofelements:
        # Convert elements to 0-based indexing
        mesh_elements = md.mesh.elements - 1

        # model_field is on elements. Take element centroid for interpolation
        mesh_x = np.mean(md.mesh.x[mesh_elements], axis = 1)
        mesh_y = np.mean(md.mesh.y[mesh_elements], axis = 1)
    else:
        raise ValueError('grid_model_field: model_field must be defined on vertices or elements')

    ## Initialise output container
    ## -------------------------------------
    ngrid_x, ngrid_y = grid_x.shape
    gridded_model_field = np.full((nt, ngrid_x, ngrid_y), np.nan)

    ## Interpolate data, one time-step at a time
    ## -------------------------------------
    for t in range(nt):
        gridded_model_field[t] = griddata(
            points = np.column_stack((mesh_x, mesh_y)),
            values = model_field[t],
            xi = (grid_x, grid_y),
            method = method
        )

    ## Mask the gridded data
    ## -------------------------------------
    # By default griddata() returns data for the convex hull of the x/y points
    # Here, we mask data to the mesh (default) or a provided domain_mask
    if domain_mask is None:
        # By default, mask data to mesh (use native md.mesh coordinates, not element centroids)
        domain_mask = make_gridded_domain_mask(md.mesh.x, md.mesh.y, mesh_elements, grid_x, grid_y)
    elif domain_mask.shape != grid_x.shape:
        # If a custom domain_mask is supplied, it must be defined on grid_x / grid_y
        raise ValueError('grid_model_field: domain_mask should be defined on grid_x / grid_y.')
    elif domain_mask.dtype != bool:
        # If a custom domain_mask is supplied, it must be boolean
        raise TypeError('grid_model_field: domain_mask should be boolean')

    # Apply mask
    gridded_model_field[:, ~domain_mask] = fill_value

    ## Squeeze output to remove time dimension if it's static
    ## -------------------------------------
    gridded_model_field = gridded_model_field.squeeze()

    return gridded_model_field

def get_element_areas_volumes(index,
                              x,
                              y,
                              z = np.array([])):
    """
    Computes areas of triangular elements or volumes of pentahedrons.

    Parameters
    ----------
    index : ndarray
        Element connectivity array. For 2D meshes, should have 3 columns.
        For 3D meshes, should have 6 columns.
    x : ndarray
        1D array of x-coordinates of mesh nodes.
    y : ndarray
        1D array of y-coordinates of mesh nodes.
    z : ndarray, optional
        1D array of z-coordinates of mesh nodes. If provided, volumes are computed.
        Default is empty array (areas computed).

    Returns
    -------
    areas : ndarray
        1D array of element areas (2D) or volumes (3D).

    Raises
    ------
    TypeError
        If x, y, and z arrays don't have the same length.
        If index contains values above the number of nodes.
        If index doesn't have the correct number of columns for the mesh type.

    Examples
    --------
    Compute areas of triangular elements:

    >>> areas = get_element_areas(md.mesh.elements, md.mesh.x, md.mesh.y)

    Compute volumes of pentahedral elements:

    >>> volumes = get_element_areas(md.mesh.elements, md.mesh.x, md.mesh.y, md.mesh.z)
    """

    ## Convert to 0-based indexing
    index = index - 1

    ## Get number of elements and number of nodes
    num_elements = np.shape(index)[0]
    num_nodes = np.shape(x)[0]

    ## Check dimensions of inputs
    if (np.shape(y)[0] != num_nodes) or (z.size > 0 and np.shape(z)[0] != num_nodes):
        raise TypeError('get_element_areas: x, y and z do not have the same length.')
    if np.max(index) > num_nodes:
        raise TypeError('get_element_areas: index should not have values above {}.'.format(num_nodes))
    if z.size == 0 and np.shape(index)[1] != 3:
        raise TypeError('get_element_areas: index should have 3 columns for 2D meshes.')
    if z.size > 0 and np.shape(index)[1] != 6:
        raise TypeError('get_element_areas: index should have 6 columns for 3D meshes.')

    ## Initialise x/y points
    areas = np.zeros(num_elements)
    x1 = x[index[:, 0]]
    x2 = x[index[:, 1]]
    x3 = x[index[:, 2]]
    y1 = y[index[:, 0]]
    y2 = y[index[:, 1]]
    y3 = y[index[:, 2]]

    ## Compute areas of each element (surface of the triangle)
    if z.size == 0:
        output = (0.5 * ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)))

    ## Compute volumes of each element (surface of the triangle * thickness)
    else:
        thickness = np.mean(z[index[:, 3:6]]) - np.mean(z[index[:, 0:3]])
        output = (0.5 * ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1))) * thickness

    return output

def get_nodal_functions_coeff(index, x, y):
    """
    Compute the coefficients alpha, beta and gamma of 2D triangular elements.
    For each triangular element, the nodal functions are defined as:
        N_i(x, y) = alpha_i * x + beta_i * y + gamma_i
    
    Parameters
    ----------
    index : numpy.ndarray
        Element connectivity array of shape (num_elements, 3). Each row contains
        the indices of the three nodes that form a triangular element. Indices
        are 1-based.
    x : numpy.ndarray
        X-coordinates of mesh nodes. Will be reshaped to column vector.
    y : numpy.ndarray
        Y-coordinates of mesh nodes. Will be reshaped to column vector.
    
    Returns
    -------
    alpha : numpy.ndarray
        Alpha coefficients of shape (num_elements, 3). Each row contains the
        alpha coefficients for the three nodal functions of an element.
    beta : numpy.ndarray
        Beta coefficients of shape (num_elements, 3). Each row contains the
        beta coefficients for the three nodal functions of an element.
    gamma : numpy.ndarray
        Gamma coefficients of shape (num_elements, 3). Each row contains the
        gamma coefficients for the three nodal functions of an element.
    
    Raises
    ------
    TypeError
        If x and y arrays have different lengths.
    TypeError
        If any index value exceeds the number of nodes.
    TypeError
        If index array does not have exactly 3 columns (non-triangular elements).
    
    Notes
    -----
    This function is specifically designed for 2D triangular finite element meshes.
    The nodal functions form a complete linear basis over each triangular element.
    """
    
    # Convert to columns
    x = x.reshape(-1)
    y = y.reshape(-1)

    # Get number of elements and number of nodes
    num_elements = np.size(index, axis = 0)
    num_nodes = np.size(x)

    # Check dimensions of inputs
    if np.size(y) != num_nodes:
        raise TypeError("pyissm.model.mesh.get_nodal_functions_coeff: x and y do not have the same length.")
    if np.max(index) > num_nodes:
        raise TypeError(f"pyissm.model.mesh.get_nodal_functions_coeff: index should not have values above {num_nodes}.")
    if np.size(index, axis=1) != 3:
        raise TypeError("pyissm.model.mesh.get_nodal_functions_coeff: only 2d meshes supported. index should have 3 columns.")

    # Initialize output
    alpha = np.zeros((num_elements, 3))
    beta = np.zeros((num_elements, 3))
    gamma = np.zeros((num_elements, 3))

    # Compute nodal functions coefficients N(x, y) = alpha x + beta y + gamma
    x1 = x[index[:, 0] - 1]
    x2 = x[index[:, 1] - 1]
    x3 = x[index[:, 2] - 1]
    y1 = y[index[:, 0] - 1]
    y2 = y[index[:, 1] - 1]
    y3 = y[index[:, 2] - 1]
    invdet = 1. / (x1 * (y2 - y3) - x2 * (y1 - y3) + x3 * (y1 - y2))

    # Get alpha, beta, and gamma
    alpha = np.vstack(((invdet * (y2 - y3)).reshape(-1, ), (invdet * (y3 - y1)).reshape(-1, ), (invdet * (y1 - y2)).reshape(-1, ))).T
    beta = np.vstack(((invdet * (x3 - x2)).reshape(-1, ), (invdet * (x1 - x3)).reshape(-1, ), (invdet * (x2 - x1)).reshape(-1, ))).T
    gamma = np.vstack(((invdet * (x2 * y3 - x3 * y2)).reshape(-1, ), (invdet * (y1 * x3 - y3 * x1)).reshape(-1, ), (invdet * (x1 * y2 - x2 * y1)).reshape(-1, ))).T

    return alpha, beta, gamma

def compute_hessian(index,
                    x,
                    y,
                    field,
                    type):
    """
    Compute the Hessian matrix from a field.

    Computes the Hessian matrix of a given field and returns the three components 
    Hxx, Hxy, Hyy for each element or each node.

    Parameters
    ----------
    index : ndarray
        Element connectivity matrix defining the triangular mesh elements.
        Shape: (num_elements, 3) with 1-based indexing.
    x : ndarray
        X-coordinates of the mesh nodes. Shape: (num_nodes,).
    y : ndarray
        Y-coordinates of the mesh nodes. Shape: (num_nodes,).
    field : ndarray
        Field values defined either on nodes or elements.
        Shape: (num_nodes,) or (num_elements,).
    type : str
        Type of output desired. Must be either 'node' or 'element'.

    Returns
    -------
    ndarray
        Hessian matrix components. Shape depends on `type`:
        
        - If type is 'element': (num_elements, 3) with columns [Hxx, Hxy, Hyy] 
          for each element.
        - If type is 'node': (num_nodes, 3) with columns [Hxx, Hxy, Hyy] 
          for each node.

    Raises
    ------
    TypeError
        If the field is not defined on nodes or elements, or if type is not
        'node' or 'element'.

    Examples
    --------
    >>> hessian = compute_hessian(md.mesh.elements, md.mesh.x, md.mesh.y, 
    ...                          md.inversion.vel_obs, 'node')
    >>> hessian = compute_hessian(md.mesh.elements, md.mesh.x, md.mesh.y,
    ...                          md.thermal.temperature, 'element')

    Notes
    -----
    The Hessian computation uses finite element nodal functions and area-weighted
    averaging for nodal values. For element-based fields, values are first
    interpolated to nodes before computing gradients and Hessian components.
    
    The Hessian matrix H has components:
    H = [[Hxx, Hxy], [Hxy, Hyy]]
    
    This function returns the three unique components as [Hxx, Hxy, Hyy].
    """

    num_nodes = np.size(x)
    num_elements = np.size(index, axis=0)

    # Error checks
    if np.size(field) not in [num_nodes, num_elements]:
        raise TypeError("Field should be defined on nodes or elements.")
    if type.lower() not in ['node', 'element']:
        raise TypeError("Type must be one of 'node' or 'element'.")

    # Flatten element connectivity for nodal accumulation
    line = index.reshape(-1, order='F')
    line_size = 3 * num_elements

    # Get areas and nodal function coefficients
    alpha, beta, gamma = get_nodal_functions_coeff(index, x, y)
    areas = get_element_areas_volumes(index, x, y)

    # Compute weights: sum of areas around each node
    weights = np.zeros(num_nodes)
    np.add.at(weights, line - 1, np.tile(areas, 3))

    # If field is element-based, interpolate to nodes
    if field.size == num_elements:
        node_field = np.zeros(num_nodes)
        np.add.at(node_field, line - 1, np.tile(areas * field, 3))
        field = node_field / weights

    # Compute gradient for each element
    grad_elx = np.sum(field[index - 1] * alpha, axis=1)
    grad_ely = np.sum(field[index - 1] * beta, axis=1)

    # Compute gradient for each node (average of surrounding elements)
    gradx = np.zeros(num_nodes)
    grady = np.zeros(num_nodes)
    np.add.at(gradx, line - 1, np.tile(areas * grad_elx, 3))
    np.add.at(grady, line - 1, np.tile(areas * grad_ely, 3))

    gradx /= weights
    grady /= weights

    # Compute Hessian for each element
    hessian = np.vstack((
        np.sum(gradx[index - 1] * alpha, axis=1),
        np.sum(grady[index - 1] * alpha, axis=1),
        np.sum(grady[index - 1] * beta, axis=1)
    )).T

    # If type is 'node', average Hessian over surrounding elements
    if type.lower() == 'node':
        hessian_node = np.zeros((num_nodes, 3))
        np.add.at(hessian_node[:, 0], line - 1, np.tile(areas * hessian[:, 0], 3))
        np.add.at(hessian_node[:, 1], line - 1, np.tile(areas * hessian[:, 1], 3))
        np.add.at(hessian_node[:, 2], line - 1, np.tile(areas * hessian[:, 2], 3))
        hessian = hessian_node / weights[:, None]

    return hessian

def compute_metric(hessian,
                   scale,
                   epsilon,
                   hmin,
                   hmax,
                   pos):
    """
    Calculates anisotropic metric tensors used for adaptive mesh 
    generation based on Hessian matrices. The metric tensor controls element 
    size and orientation in the mesh by analyzing eigenvalues and eigenvectors
    of the Hessian matrix.

    Parameters
    ----------
    hessian : numpy.ndarray
        Array of shape (n, 3) containing Hessian matrix components for each node.
        Columns represent [H11, H12, H22] where H is the 2x2 Hessian matrix.
    scale : float
        Scaling factor for the metric computation.
    epsilon : float
        Tolerance parameter used in the metric scaling calculation.
    hmin : float
        Minimum allowed element size in the mesh.
    hmax : float
        Maximum allowed element size in the mesh.
    pos : numpy.ndarray
        Array of indices corresponding to water elements or special boundary 
        conditions that require uniform metric treatment.

    Returns
    -------
    numpy.ndarray
        Array of shape (n, 3) containing the computed metric tensor components
        [M11, M12, M22] for each node, where M is the 2x2 symmetric metric tensor.

    Raises
    ------
    RuntimeError
        If NaN values persist in the metric tensor after all correction attempts.

    Notes
    -----
    The function performs the following key operations:
    1. Computes eigenvalues and eigenvectors of the Hessian matrix
    2. Applies size constraints using hmin and hmax parameters
    3. Handles special cases (zero eigenvalues, water elements)
    4. Uses numpy.linalg.eig as a fallback for numerical issues
    5. Ensures the resulting metric is free of NaN values
    The metric tensor M is used in adaptive mesh generation where element
    sizes are controlled by the relationship: h^T * M * h = 1, where h
    represents the edge vector in the mesh.

    Examples
    --------
    >>> hessian = compute_hessian(md.mesh.elements, md.mesh.x, md.mesh.y, 
    ...                          md.inversion.vel_obs, 'node')
    >>> metric = compute_metric(hessian, 1.0, 0.01, 0.1, 10.0, np.array([]))
    """

    # Find the eigen values of each line of H = [hessian(i, 1) hessian(i, 2); hessian(i, 2) hessian(i, 3)]
    a = hessian[:, 0]
    b = hessian[:, 1]
    d = hessian[:, 2]
    lambda1 = 0.5 * ((a + d) + np.sqrt(4. * b**2 + (a - d)**2))
    lambda2 = 0.5 * ((a + d) - np.sqrt(4. * b**2 + (a - d)**2))
    pos1 = np.nonzero(lambda1 == 0.)[0]
    pos2 = np.nonzero(lambda2 == 0.)[0]
    pos3 = np.nonzero(np.logical_and(b == 0., lambda1 == lambda2))[0]

    # Modify eigen values to control the shape of the elements
    lambda1 = np.minimum(np.maximum(np.abs(lambda1) * scale / epsilon, 1. / hmax**2), 1. / hmin**2)
    lambda2 = np.minimum(np.maximum(np.abs(lambda2) * scale / epsilon, 1. / hmax**2), 1. / hmin**2)

    # Compute eigen vectors
    norm1 = np.sqrt(8. * b**2 + 2. * (d - a)**2 + 2. * (d - a) * np.sqrt((a - d)**2 + 4. * b**2))
    v1x = 2. * b / norm1
    v1y = ((d - a) + np.sqrt((a - d)**2 + 4. * b**2)) / norm1
    norm2 = np.sqrt(8. * b**2 + 2. * (d - a)**2 - 2. * (d - a) * np.sqrt((a - d)**2 + 4. * b**2))
    v2x = 2. * b / norm2
    v2y = ((d - a) - np.sqrt((a - d)**2 + 4. * b**2)) / norm2

    v1x[pos3] = 1.
    v1y[pos3] = 0.
    v2x[pos3] = 0.
    v2y[pos3] = 1.

    # Compute new metric (for each node M = V * Lambda * V^-1)
    metric = np.vstack((((v1x * v2y - v1y * v2x)**(-1) * (lambda1 * v2y * v1x - lambda2 * v1y * v2x)).reshape(-1, ),
                        ((v1x * v2y - v1y * v2x)**(-1) * (lambda1 * v1y * v2y - lambda2 * v1y * v2y)).reshape(-1, ),
                        ((v1x * v2y - v1y * v2x)**(-1) * (-lambda1 * v2x * v1y + lambda2 * v1x * v2y)).reshape(-1, ))).T

    # Corrections for 0 eigen values
    metric[pos1, :] = np.tile(np.array([[1. / hmax**2, 0., 1. / hmax**2]]), (np.size(pos1), 1))
    metric[pos2, :] = np.tile(np.array([[1. / hmax**2, 0., 1. / hmax**2]]), (np.size(pos2), 1))

    # Handle water elements
    metric[pos, :] = np.tile(np.array([[1. / hmax**2, 0., 1. / hmax**2]]), (np.size(pos), 1))

    # Handle NaNs if any (use Numpy eig in a loop)
    pos = np.nonzero(np.isnan(metric))[0]
    if np.size(pos):
        print((f"pyissm.model.mesh.compute_metric: {np.size(pos)} NaNs found in the metric. Use Numpy routine to fix them."))
        for posi in pos:
            H = np.array([[hessian[posi, 0], hessian[posi, 1]], [hessian[posi, 1], hessian[posi, 2]]])
            [v, u] = np.linalg.eig(H)
            v = np.diag(v)
            lambda1 = v[0, 0]
            lambda2 = v[1, 1]
            v[0, 0] = np.minimum(np.maximum(np.abs(lambda1) * scale / epsilon, 1. / hmax**2), 1. / hmin**2)
            v[1, 1] = np.minimum(np.maximum(np.abs(lambda2) * scale / epsilon, 1. / hmax**2), 1. / hmin**2)

            metricTria = np.dot(np.dot(u, v), np.linalg.inv(u))
            metric[posi, :] = np.array([metricTria[0, 0], metricTria[0, 1], metricTria[1, 1]])

    if np.any(np.isnan(metric)):
        raise RuntimeError("pyissm.model.mesh.compute_metric: NaN in the metric despite our efforts...")
    
    return metric

def elements_from_edge(elements, A, B):
    """
    Find elements connected to one edge defined by nodes A and B.

    Parameters
    ----------
    elements : array_like
        Array of element connectivity information where each row represents 
        an element and columns represent the nodes that define the element.
    A : int
        First node ID defining the edge.
    B : int
        Second node ID defining the edge.

    Returns
    -------
    ndarray
        1D array of element IDs (1-based indexing) that contain the edge 
        defined by nodes A and B.

    Examples
    --------
    >>> edgeelements = elements_from_edge(md.mesh.elements, node1, node2)
    """

    # Broadcast A vs B to all 3 combinations of node pairs in each triangle
    mask = ((elements == A)[:, :, None] & (elements == B)[:, None, :]).any(axis=(1, 2))

    # Convert to 1-based indexing
    edgeelements = np.nonzero(mask)[0] + 1

    return edgeelements

def export_gmsh():
    """
    Export the model mesh to a Gmsh .msh file.
    
    Raises
    ------
    NotImplementedError
        Function is not yet implemented.
    """
    raise NotImplementedError("pyissm.model.mesh.export_gmsh: This functionality is not yet implemented. Please contact ACCESS-NRI for support.")

def find_segments():
    """
    Build segments model field

    Raises
    ------
    NotImplementedError
        Function is not yet implemented.
    """

    raise NotImplementedError("pyissm.model.mesh.find_segments: This functionality is not yet implemented. Please contact ACCESS-NRI for support.")

def flag_elements(md, region = 'all', inside = True):
    """
    Flag elements based on their location within the model domain.

    This function allows users to flag elements in the mesh based on whether they
    are inside or outside a specified domain. The region can be the entire mesh,
    no elements, a region specified by a provided *.exp file, or defined by boolean arrays.

    Parameters
    ----------
    md : ISSM Model object
        ISSM Model object containing the mesh
    region : str or ndarray, optional
        Region specification. Options are:
        
        - 'all' (default): Flag all elements in the mesh.
        - '': Flag no elements.
        - Path to a *.exp file: Flag elements inside or outside the polygon defined in the file.
        - ndarray: Boolean array of size (numberofelements,) or (numberofvertices,).
          If vertices array, elements are flagged when all vertices are flagged.
    inside : bool, optional
        If `region` is a polygon file or array, this parameter specifies whether to flag
        elements inside (`True`, default) or outside (`False`) the region.
        Default is True.

    Returns
    -------
    ndarray of bool
        Boolean array of length `md.mesh.numberofelements` where `True` indicates
        flagged elements and `False` indicates unflagged elements.

    Raises
    ------
    RuntimeError
        If python wrappers are not installed and a *.exp file is provided.
    ValueError
        If region array does not match number of elements or vertices.
    TypeError
        If region is neither a string nor an array.

    Examples
    --------
    Flag all elements in the mesh:

    >>> flags = flag_elements(md)
    >>> flags = flag_elements(md, region='all')

    Flag no elements:

    >>> flags = flag_elements(md, region='')

    Flag elements inside a polygon:

    >>> flags = flag_elements(md, region='path/to/polygon.exp')

    Flag elements outside a polygon:

    >>> flags = flag_elements(md, region='path/to/polygon.exp', inside=False)
    """

    # If region is None, flag no elements
    if region is None:
        flag = np.zeros(md.mesh.numberofelements, dtype=bool)    

    # If region is a string, check if it's 'all', or a file path
    elif isinstance(region, str):
        ## If 'all', flag all elements
        if region.lower() == 'all':
            flag = np.ones(md.mesh.numberofelements, dtype=bool)

        ## If a file path, load polygon and flag elements inside or outside
        elif region.endswith('.exp'):
            if tools.wrappers.check_wrappers_installed():
                flag = tools.wrappers.ContourToMesh(md.mesh.elements[:, 0:3], md.mesh.x, md.mesh.y, region, 'element', 1).astype(bool)
            else:
                raise RuntimeError('pyissm.model.mesh.flag_elements: Python wrappers not installed. Cannot flag elements from *.exp file.')

        ## If inside is False, invert the flag to get outside elements
        if not inside:
            flag = np.logical_not(flag)

    # If region is an array, it must the same size as number of elements or vertices
    elif isinstance(region, np.ndarray):
        ## If region is an array of elements, use it directly
        if region.shape[0] == md.mesh.numberofelements:
            flag = region.astype(bool)
        
        ## If region is an array of vertices, flag elements where all vertices are in the region
        elif region.shape[0] == md.mesh.numberofvertices:
            flag = (np.sum(region[md.mesh.elements -1] > 0, axis=1) == md.mesh.elements.shape[1]).astype(bool)
        else:
            raise ValueError("Flag list for region must be of same size as number of elements or vertices.")
        
        ## If inside is False, invert the flag to get outside elements
        if not inside:
            flag = np.logical_not(flag)
    
    # If region is neither a string nor an array, raise an error
    else:
        raise TypeError("Region must be None, a string ('all' or path to *.exp file), or a boolean array.")

    return flag

def triangle(md,
             domain_name,
             resolution,
             rift_name = None):
    
    """
    Create a triangular mesh for an ISSM model using Triangle mesh generator.

    This function generates a triangular mesh based on a domain outline and optional
    rift constraints. It uses the Triangle mesh generator to create high-quality
    Delaunay triangulations with specified resolution constraints.

    Parameters
    ----------
    md : object
        ISSM model object containing the mesh structure to be populated.
    domain_name : str
        Path to the file containing the domain outline geometry.
    resolution : float
        Target mesh resolution in meters. This represents the characteristic
        edge length for mesh elements. The actual mesh area constraint is
        calculated as resolution squared.
    rift_name : str, optional
        Path to the file containing rift constraint geometry. If provided,
        these constraints will be incorporated into the mesh generation.
        Default is None.

    Returns
    -------
    md : object
        The input ISSM model object with updated mesh properties including:
        - mesh.x, mesh.y: Node coordinates
        - mesh.elements: Element connectivity matrix
        - mesh.segments: Boundary segment definitions
        - mesh.segmentmarkers: Boundary segment markers
        - mesh.numberofvertices: Total number of mesh vertices
        - mesh.numberofelements: Total number of mesh elements
        - mesh.vertexonboundary: Boolean array indicating boundary vertices
        - mesh.vertexconnectivity: Vertex-to-vertex connectivity
        - mesh.elementconnectivity: Element-to-element connectivity

    Raises
    ------
    IOError
        If the domain outline file or rift file (when specified) does not exist.
    RuntimeError
        If the Triangle Python wrappers are not installed and mesh creation fails.

    Warnings
    --------
    - Issues a warning if the existing mesh is not empty and will be overwritten.
    - Issues a warning if orphaned nodes are found and removed from the mesh.

    Notes
    -----
    This function requires the Triangle Python wrappers to be installed for
    mesh generation. If wrappers are not available, the function raises a RuntimeError.
    The function automatically handles orphaned nodes (nodes not belonging to
    any element) by removing them and updating the connectivity accordingly.

    Examples
    --------
    >>> import pyissm
    >>> md = pyissm.Model()
    >>> md = pyissm.model.mesh.triangle(md, 'domain.exp', 1000.0)
    >>> md = pyissm.model.mesh.triangle(md, 'domain.exp', 500.0, 'rifts.exp')
    """

    # Error checks
    ## Check if md mesh is empty
    if md.mesh.numberofelements:
        raise RuntimeError('md.mesh is not empty. Use md.mesh = pyissm.model.classes.mesh.mesh2d() to reset the mesh.')

    ## Check if file(s) exist
    if not os.path.exists(domain_name):
        raise IOError(f"Domain outline file {domain_name} does not exist.")
    
    if rift_name is not None and not os.path.exists(rift_name):
        raise IOError(f"Rift file {rift_name} does not exist.")
    
    ## Check if wrappers are installed
    if not tools.wrappers.check_wrappers_installed():
        raise RuntimeError('pyissm.model.mesh.triangle: Python wrappers not installed. Cannot create mesh.')


    # Calculate characteristic area. Resolution is node-oriented. A 1000 m resolution = 1000 * 1000 area
    area = resolution ** 2

    # Make the mesh using Triangle_python
    ## NOTE: Check for wrappers already done above
    elements, x, y, segments, segmentmarkers = tools.wrappers.Triangle(domain_name, rift_name, area)

    # Check that all created nodes belong to at least one element
    ## NOTE: Orphan node removal taken from $ISSM_DIR/src/m/mesh/triangle.m
    uniqueelements = np.sort(np.unique(elements))
    orphans = np.nonzero((~np.isin(range(1, len(x)), uniqueelements)).astype(int))[0]
    if len(orphans):
        warnings.warn(f'pyissm.model.mesh.triangle: {len(orphans)} orphaned nodes found. These nodes do not belong to any element.\n'
                      'Removing orphaned nodes from the mesh.')
        
        for i in range(0, len(orphans)):
            print('WARNING: removing orphans')
            # Get rid of the orphan node i
            # Update x and y
            x = np.concatenate((x[0:(orphans[i] - i)], x[(orphans[i] - i + 1):]))
            y = np.concatenate((y[0:(orphans[i] - i)], y[(orphans[i] - i + 1):]))
            # Update elements
            pos = np.nonzero((elements > (orphans[i] - i)).flatten(order = 'F'))[0]
            elementstmp = elements.flatten(order = 'F')
            elementstmp[pos] -= 1
            elements = elementstmp.reshape(np.shape(elements), order = 'F')
            # Update segments
            pos1 = np.nonzero(segments[:,0] > (orphans[i] - i))[0]
            pos2 = np.nonzero(segments[:,1] > (orphans[i] - i))[0]
            segments[pos1, 0] -= 1
            segments[pos2, 1] -= 1
    
    # Assign to md structure
    md.mesh = model.classes.mesh.mesh2d()
    md.mesh.x = x
    md.mesh.y = y
    md.mesh.elements = elements.astype(int)
    md.mesh.segments = segments.astype(int)
    md.mesh.segmentmarkers = segmentmarkers.astype(int)
    md.mesh.numberofvertices = len(x)
    md.mesh.numberofelements = len(elements)
    md.mesh.vertexonboundary = np.zeros(md.mesh.numberofvertices, int)
    md.mesh.vertexonboundary[md.mesh.segments[:,0:2] - 1] = 1

    ## Build connectivity arrays
    ### NOTE: Check for wrappers already done above
    md.mesh.vertexconnectivity = tools.wrappers.NodeConnectivity(md.mesh.elements, md.mesh.numberofvertices)
    md.mesh.elementconnectivity = tools.wrappers.ElementConnectivity(md.mesh.elements, md.mesh.vertexconnectivity)

    return md

def square_mesh(md,
                Lx,
                Ly,
                nx,
                ny):
    """
    Create a structured triangular mesh on a rectangular domain.
    
    This function generates a structured triangular mesh on a rectangular domain
    and populates the mesh fields of an ISSM model object. The mesh consists of
    triangular elements arranged in a regular grid pattern.

    Parameters
    ----------
    md : object
        ISSM model object whose mesh fields will be populated.
    Lx : float
        Length of the domain in the x-direction.
    Ly : float
        Length of the domain in the y-direction.
    nx : int
        Number of nodes in the x-direction.
    ny : int
        Number of nodes in the y-direction.
    
    Raises
    ------
    RuntimeError
        If Python wrappers are not installed.
    
    Warnings
    --------
    UserWarning
        If the model mesh is not empty and will be overwritten.
    
    Notes
    -----
    The function creates a structured triangular mesh by:
    1. Generating node coordinates on a regular grid
    2. Creating triangular elements by splitting each grid cell into two triangles
    3. Defining boundary segments around the domain perimeter
    4. Computing connectivity arrays for mesh topology
    The resulting mesh will have (nx-1)*(ny-1)*2 triangular elements and nx*ny nodes.

    Examples
    --------
    >>> import pyissm
    >>> md = pyissm.Model()
    >>> md = pyissm.model.mesh.square_mesh(md, Lx=100.0, Ly=50.0, nx=11, ny=6)
    """  

    # Error checks
    ## Check if md mesh is empty
    if md.mesh.numberofelements:
        warnings.warn('md.mesh is not empty. Overwriting existing mesh.')

    ## Check if wrappers are installed
    if not tools.wrappers.check_wrappers_installed():
        raise RuntimeError('pyissm.model.mesh.square_mesh: Python wrappers not installed. Cannot create mesh.')

    # Get number of elements and nodes
    n_elements = (nx - 1) * (ny - 1) * 2
    n_nodes = nx * ny

    # Initialize arrays
    index = np.zeros((n_elements, 3), int)
    x = np.zeros((n_nodes))
    y = np.zeros((n_nodes))

    # Create coordinates
    for n in range(0, nx):
        for m in range(0, ny):
            x[n * ny + m] = float(n)
            y[n * ny + m] = float(m)

    # Create index
    for n in range(0, nx - 1):
        for m in range(0, ny - 1):
            A = n * ny + (m + 1)
            B = A + 1
            C = (n + 1) * ny + (m + 1)
            D = C + 1
            index[n * (ny - 1) * 2 + 2 * m, :] = [A, C, B]
            index[n * (ny - 1) * 2 + 2 * (m + 1) - 1, :] = [B, C, D]

    # Scale x and y
    x = x / np.max(x) * Lx
    y = y / np.max(y) * Ly

    # Create segments
    segments = np.zeros((2 * (nx - 1) + 2 * (ny - 1), 3), int)
    segments[0:ny - 1, :] = np.vstack((np.arange(2, ny + 1), np.arange(1, ny), (2 * np.arange(1, ny) - 1))).T
    segments[ny - 1:2 * (ny - 1), :] = np.vstack((np.arange(ny * (nx - 1) + 1, nx * ny), np.arange(ny * (nx - 1) + 2, nx * ny + 1), 2 * np.arange((ny - 1) * (nx - 2) + 1, (nx - 1) * (ny - 1) + 1))).T
    segments[2 * (ny - 1):2 * (ny - 1) + (nx - 1), :] = np.vstack((np.arange(2 * ny, ny * nx + 1, ny), np.arange(ny, ny * (nx - 1) + 1, ny), np.arange(2 * (ny - 1), 2 * (nx - 1) * (ny - 1) + 1, 2 * (ny - 1)))).T
    segments[2 * (ny - 1) + (nx - 1):2 * (nx - 1) + 2 * (ny - 1), :] = np.vstack((np.arange(1, (nx - 2) * ny + 2, ny), np.arange(ny + 1, ny * (nx - 1) + 2, ny), np.arange(1, 2 * (nx - 2) * (ny - 1) + 2, 2 * (ny - 1)))).T

    # Populate md structure
    md.mesh = model.classes.mesh.mesh2d()
    md.mesh.x = x
    md.mesh.y = y
    md.mesh.elements = index.astype(int)
    md.mesh.segments = segments.astype(int)
    md.mesh.numberofvertices = len(x)
    md.mesh.numberofelements = len(index)
    md.mesh.vertexonboundary = np.zeros(md.mesh.numberofvertices, int)
    md.mesh.vertexonboundary[md.mesh.segments[:,0:2] - 1] = 1

    ## Build connectivity arrays
    ### NOTE: Check for wrappers already done above
    md.mesh.vertexconnectivity = tools.wrappers.NodeConnectivity(md.mesh.elements, md.mesh.numberofvertices)
    md.mesh.elementconnectivity = tools.wrappers.ElementConnectivity(md.mesh.elements, md.mesh.vertexconnectivity)

    return md

def round_mesh(md,
               radius,
               resolution,
               exp_output_name = None,
               keep_exp = False):
    
    """
    Create a structured triangular mesh on a circular domain.
    
    This function generates a triangular mesh on a circular domain by first
    creating a domain outline file and then using the Triangle mesh generator.
    The mesh is created with approximately uniform resolution around the circle
    perimeter and moves the closest node to the origin.

    Parameters
    ----------
    md : object
        ISSM model object whose mesh fields will be populated.
    radius : float
        Radius of the circular domain in meters.
    resolution : float
        Target mesh resolution in meters. This represents the characteristic
        edge length for mesh elements around the circle perimeter.
    exp_output_name : str, optional
        Path for the output domain outline file (.exp format). If None,
        defaults to 'round_mesh.exp'. Default is None.
    keep_exp : bool, optional
        Whether to keep the temporary domain outline file after mesh creation.
        If False, the file is automatically deleted. Default is False.

    Returns
    -------
    md : object
        The input ISSM model object with updated mesh properties including:
        - mesh.x, mesh.y: Node coordinates
        - mesh.elements: Element connectivity matrix
        - mesh.segments: Boundary segment definitions
        - mesh.segmentmarkers: Boundary segment markers
        - mesh.numberofvertices: Total number of mesh vertices
        - mesh.numberofelements: Total number of mesh elements
        - mesh.vertexonboundary: Boolean array indicating boundary vertices
        - mesh.vertexconnectivity: Vertex-to-vertex connectivity
        - mesh.elementconnectivity: Element-to-element connectivity

    Raises
    ------
    IOError
        If the specified exp_output_name file already exists.
    RuntimeError
        If the Triangle Python wrappers are not installed.

    Warnings
    --------
    UserWarning
        If the model mesh is not empty and will be overwritten.

    Notes
    -----
    The function creates a circular mesh by:
    1. Generating points uniformly distributed around the circle perimeter
    2. Writing these points to a domain outline file (.exp format)
    3. Using the Triangle mesh generator to create the triangular mesh
    4. Moving the closest node to the origin (0,0) for convenience
    5. Optionally removing the temporary domain outline file

    The number of points on the circle perimeter is calculated based on the
    target resolution to ensure approximately uniform spacing.

    Examples
    --------
    >>> import pyissm
    >>> md = pyissm.Model()
    >>> md = pyissm.model.mesh.round_mesh(md, radius=5000.0, resolution=500.0)
    >>> md = pyissm.model.mesh.round_mesh(md, radius=1000.0, resolution=100.0, 
    ...                                   exp_output_name='circle.exp', keep_exp=True)
    """

    # Internal helper function
    def _round_sig_fig(x, n):
        nonzeros = np.where(x != 0)
        digits = np.ceil(np.log10(np.abs(x[nonzeros])))
        x[nonzeros] = x[nonzeros] / 10.**digits
        x[nonzeros] = np.round(x[nonzeros], decimals=n)
        x[nonzeros] = x[nonzeros] * 10.**digits
        return x
    
    ## ----------------------------------------------------------

    # Error checks
    ## NOTE: Existing mesh & wrapper installation checks handled in triangle()

    ## Check if file(s) exist. Do not overwrite existing files.
    if exp_output_name is not None and os.path.exists(exp_output_name):
        raise IOError(f"exp_output_name file {exp_output_name} already exists.")
    
    # Create the domain outline file
    if exp_output_name is None:
        exp_output_name = 'round_mesh.exp'

    # Construct the mesh
    ## Get number of points on the circle
    pointsonedge = int(np.floor((2. * np.pi * radius) / resolution) + 1)  # +1 to close the outline

    ## Calculate the Cartesian coordinates of the points
    theta = np.linspace(0., 2. * np.pi, pointsonedge)
    x_list = _round_sig_fig(radius * np.cos(theta), 12)
    y_list = _round_sig_fig(radius * np.sin(theta), 12)

    ## Create the contour dictionary
    contour = collections.OrderedDict()
    contour['x'] = x_list
    contour['y'] = y_list
    contour['density'] = 1.

    ## Write the contour to an exp file
    tools.exp.exp_write(contour, exp_output_name)

    ## Create the mesh using triangle
    md = triangle(md, exp_output_name, resolution)

    ## Move closest node to (0,0)
    min_pos = np.argmin(np.sqrt(md.mesh.x**2 + md.mesh.y**2))
    md.mesh.x[min_pos] = 0.
    md.mesh.y[min_pos] = 0.

    # Remove temporary exp file if not required
    if not keep_exp:
        os.remove(exp_output_name)

    return md

def _bamg_geom(**kwargs):
    """
    Initialize a BAMG geometry dictionary with default empty numpy arrays.

    This function initializes a geometry dictionary with default empty numpy arrays
    for various BAMG (Bidimensional Anisotropic Mesh Generator) geometry components
    and allows users to override any of these defaults through keyword arguments.

    NOTE: Intended for internal used within BAMG meshing functions here only.

    Parameters
    ----------
    **kwargs : dict, optional
        User-specified geometry options to override defaults. Valid keys include:
        
        - 'Vertices' : array_like, shape (n, 3)
            Vertex coordinates [x, y, marker]
        - 'Edges' : array_like, shape (n, 3) 
            Edge definitions [vertex1, vertex2, marker]
        - 'TangentAtEdges' : array_like, shape (n, 4)
            Tangent vectors at edges [edge_id, tx, ty, marker]
        - 'Corners' : array_like, shape (n, 1)
            Corner vertex indices
        - 'RequiredVertices' : array_like, shape (n, 1)
            Indices of vertices that must be preserved
        - 'RequiredEdges' : array_like, shape (n, 1)
            Indices of edges that must be preserved
        - 'CrackedEdges' : array_like, shape (n, 0)
            Cracked edge definitions
        - 'SubDomains' : array_like, shape (n, 4)
            Subdomain specifications [x, y, marker, area_constraint]

    Returns
    -------
    dict
        BAMG geometry dictionary with default empty arrays updated by user options.
        
    Examples
    --------
    >>> # Create default geometry
    >>> geom = _bamg_geom()
    >>> 
    >>> # Create geometry with custom vertices
    >>> vertices = np.array([[0, 0, 1], [1, 0, 1], [0, 1, 1]])
    >>> geom = _bamg_geom(Vertices=vertices)
    """

    geom = {
        'Vertices': np.empty((0, 3)),
        'Edges': np.empty((0, 3)),
        'TangentAtEdges': np.empty((0, 4)),
        'Corners': np.empty((0, 1)),
        'RequiredVertices': np.empty((0, 1)),
        'RequiredEdges': np.empty((0, 1)),
        'CrackedEdges': np.empty((0, 0)),
        'SubDomains': np.empty((0, 4)),
    }

    # Update defaults with user-specified options
    geom.update(kwargs)

    return geom

def _bamg_mesh(**kwargs):
    """
    Initialize a BAMG mesh dictionary with default empty numpy arrays.

    This function initializes a mesh dictionary with default empty numpy arrays
    for various BAMG (Bidimensional Anisotropic Mesh Generator) mesh components
    and allows users to override any of these defaults through keyword arguments.

    NOTE: Intended for internal use within BAMG meshing functions here only.

    Parameters
    ----------
    **kwargs : dict, optional
        User-specified mesh options to override defaults. Valid keys include:
        
        - 'Vertices' : array_like, shape (n, 3)
            Vertex coordinates [x, y, marker]
        - 'Edges' : array_like, shape (n, 3) 
            Edge definitions [vertex1, vertex2, marker]
        - 'Triangles' : array_like, shape (n, 0)
            Triangle element definitions
        - 'IssmEdges' : array_like, shape (n, 0)
            ISSM-specific edge data
        - 'IssmSegments' : array_like, shape (n, 0)
            ISSM-specific boundary segment data
        - 'VerticesOnGeomVertex' : array_like, shape (n, 0)
            Vertices located on geometry vertices
        - 'VerticesOnGeomEdge' : array_like, shape (n, 0)
            Vertices located on geometry edges
        - 'EdgesOnGeomEdge' : array_like, shape (n, 0)
            Edges located on geometry edges
        - 'SubDomains' : array_like, shape (n, 4)
            Subdomain specifications [x, y, marker, area_constraint]
        - 'SubDomainsFromGeom' : array_like, shape (n, 0)
            Subdomains derived from geometry
        - 'ElementConnectivity' : array_like, shape (n, 0)
            Element-to-element connectivity matrix
        - 'NodalConnectivity' : array_like, shape (n, 0)
            Node-to-node connectivity matrix
        - 'NodalElementConnectivity' : array_like, shape (n, 0)
            Node-to-element connectivity matrix
        - 'CrackedVertices' : array_like, shape (n, 0)
            Vertices involved in crack features
        - 'CrackedEdges' : array_like, shape (n, 0)
            Edges involved in crack features

    Returns
    -------
    dict
        BAMG mesh dictionary with default empty arrays updated by user options.
        
    Examples
    --------
    >>> # Create default mesh dictionary
    >>> mesh = _bamg_mesh()
    >>> 
    >>> # Create mesh dictionary with custom vertices
    >>> vertices = np.array([[0, 0, 1], [1, 0, 1], [0, 1, 1]])
    >>> mesh = _bamg_mesh(Vertices=vertices)
    """

    mesh = {
        'Vertices': np.empty((0, 3)),
        'Edges': np.empty((0, 3)),
        'Triangles': np.empty((0, 0)),
        'IssmEdges': np.empty((0, 0)),
        'IssmSegments': np.empty((0, 0)),
        'VerticesOnGeomVertex': np.empty((0, 0)),
        'VerticesOnGeomEdge': np.empty((0, 0)),
        'EdgesOnGeomEdge': np.empty((0, 0)),
        'SubDomains': np.empty((0, 4)),
        'SubDomainsFromGeom': np.empty((0, 0)),
        'ElementConnectivity': np.empty((0, 0)),
        'NodalConnectivity': np.empty((0, 0)),
        'NodalElementConnectivity': np.empty((0, 0)),
        'CrackedVertices': np.empty((0, 0)),
        'CrackedEdges': np.empty((0, 0))
    }

    # Update defaults with user-specified options
    mesh.update(kwargs)

    return mesh

def bamg(md, **kwargs):
    """
    Create a triangular mesh using the BAMG (Bidimensional Anisotropic Mesh Generator) algorithm.

    This function generates high-quality anisotropic triangular meshes for complex 
    geometries using the BAMG mesh generator. It supports various mesh constraints 
    including domain boundaries, holes, subdomains, rifts, and anisotropic metrics.

    Parameters
    ----------
    md : object
        ISSM model object whose mesh fields will be populated with the generated mesh.

    kwargs : dict, optional
        Additional keyword arguments to customize mesh generation. Supported options include:
        - anisomax : float, maximum anisotropy ratio allowed in the mesh (default= 1e30)
        - coeff : float, global mesh size coefficient multiplier (default= 1.0)
        - Crack : int, enable crack processing (0=disabled, 1=enabled) (default= 0)
        - cutoff : float, cutoff value for metric interpolation (default= 1e-5)
        - domain : str or list, path to domain outline file or list of domain contours (default= None)
        - err : float, target interpolation error for mesh adaptation (default= 0.01)
        - errg : float, target geometric error for mesh adaptation (default= 0.1)
        - field : array_like, field values for metric computation (default= empty array)
        - gradation : float, mesh size gradation parameter (default= 1.5)
        - Hessiantype : int, type of Hessian computation (0=P1, 1=P2) (default= 0)
        - hmax : float, maximum allowed edge length (default= 1e100)
        - hmin : float, minimum allowed edge length (default= 1e-100)
        - hmaxVertices : array_like, maximum edge lengths at specific vertices (default= empty array)
        - hminVertices : array_like, minimum edge lengths at specific vertices (default= empty array)
        - holes : str or list, path to holes file or list of hole contours (default= None)
        - hVertices : array_like, target edge lengths at specific vertices (default= empty array)
        - KeepVertices : int, keep vertices from previous mesh (0=no, 1=yes) (default= 1)
        - Markers : array_like, edge markers for boundary identification (default= None)
        - maxnbv : float, maximum number of vertices allowed (default= 1.0e6)
        - maxsubdiv : float, maximum number of edge subdivisions (default= 10.0)
        - metric : array_like, anisotropic metric tensor field (default= empty array)
        - Metrictype : int, type of metric (0=isotropic, 1=anisotropic) (default= 0)
        - nbjacobi : int, number of Jacobi smoothing iterations (default= 1)
        - nbsmooth : int, number of mesh smoothing iterations (default= 3)
        - NoBoundaryRefinement : int, disable boundary refinement for domain edges (0=allow, 1=disable) (default= 0)
        - NoBoundaryRefinementAllBoundaries : int, disable boundary refinement for all edges (0=allow, 1=disable) (default= 0)
        - omega : float, relaxation parameter for smoothing (default= 1.8)
        - power : float, power for metric computation (default= 1.0)
        - RequiredVertices : array_like, coordinates of vertices that must be included in the mesh (default= None)
        - rifts : str, path to rifts file for fracture modeling (default= None)
        - splitcorners : int, split corners in mesh generation (0=no, 1=yes) (default= 1)
        - subdomains : str or list, path to subdomains file or list of subdomain contours (default= None)
        - tol : float, tolerance for geometric operations (default= None)
        - toltip : float, tolerance for rift tip processing (default= None)
        - tracks : str or array_like, path to tracks file or track coordinates (default= None)
        - verbose : int, verbosity level (0=quiet, 1=verbose) (default= 1)
        - vertical : int, create 2D vertical mesh (0=standard 2D, 1=vertical) (default= 0)
        - 3dsurface : int, create 3D surface mesh (0=standard 2D, 1=3D surface) (default= 0)

    Returns
    -------
    md : object
        The input ISSM model object with updated mesh properties including:
        - mesh.x, mesh.y: Node coordinates
        - mesh.elements: Element connectivity matrix
        - mesh.edges: Edge connectivity matrix
        - mesh.segments: Boundary segment definitions
        - mesh.segmentmarkers: Boundary segment markers
        - mesh.numberofvertices: Total number of mesh vertices
        - mesh.numberofelements: Total number of mesh elements
        - mesh.numberofedges: Total number of mesh edges
        - mesh.vertexonboundary: Boolean array indicating boundary vertices
        - mesh.elementconnectivity: Element-to-element connectivity
        - private.bamg: BAMG-specific mesh and geometry data

    Raises
    ------
    IOError
        If specified input files (domain, holes, subdomains, rifts) do not exist.
    RuntimeError
        If mesh generation fails or if incompatible options are specified.
    TypeError
        If input arguments are of incorrect type.

    Notes
    -----
    This function is a comprehensive interface to the BAMG mesh generator, supporting
    complex geometries with multiple constraints. The mesh can be adapted based on
    metric fields for anisotropic meshing. Special handling is provided for rifts
    and fractures in ice sheet modeling applications.

    Examples
    --------
    >>> import pyissm
    >>> md = pyissm.Model()
    >>> # Basic domain meshing
    >>> md = pyissm.model.mesh.bamg(md, domain='outline.exp', hmax=1000.0)
    >>> # Anisotropic meshing with metric
    >>> md = pyissm.model.mesh.bamg(md, domain='outline.exp', 
    ...                             metric=metric_field, err=0.005)
    >>> # Mesh with holes and subdomains
    >>> md = pyissm.model.mesh.bamg(md, domain='outline.exp', 
    ...                             holes='holes.exp', subdomains='regions.exp')
    """

    # Define default options
    defaults = {
        "anisomax": 1e30,
        "coeff": 1.0,
        "Crack": 0,
        "cutoff": 1e-5,
        "domain": None,
        "err": 0.01,
        "errg": 0.1,
        "field": np.empty((0, 1)),
        "gradation": 1.5,
        "Hessiantype": 0,
        "hmax": 1e100,
        "hmin": 1e-100,
        "hmaxVertices": np.empty((0, 1)),
        "hminVertices": np.empty((0, 1)),
        "holes": None,
        "hVertices": np.empty((0, 1)),
        "KeepVertices": 1,
        "Markers": None,
        "maxnbv": 1.0e6,
        "maxsubdiv": 10.0,
        "metric": np.empty((0, 1)),
        "Metrictype": 0,
        "nbjacobi": 1,
        "nbsmooth": 3,
        "NoBoundaryRefinement": 0,
        "NoBoundaryRefinementAllBoundaries": 0,
        "omega": 1.8,
        "power": 1.0,
        "RequiredVertices": None,
        "rifts": None,
        "splitcorners": 1,
        "subdomains": None,
        "tol": None,
        "toltip": None,
        "tracks": None,
        "verbose": 1,
        "vertical": 0,
        "3dsurface": 0
    }

    # Define helper functions
    def _load_spatial_components(component):
        """
        Load spatial components from file or directly from input.
        """

        if component is not None:

            ## Check the file exists if a filename is provided
            if isinstance(component, str):
                if not os.path.exists(component):
                    raise IOError(f"BAMG spatial component file {component} does not exist.")
                return tools.exp.exp_read(component)
            
            ## If a list or dict is provided, it must be a list of dictionaries
            elif isinstance(component, list):
                if len(component):
                    if all(isinstance(c, (dict, collections.OrderedDict)) for c in component):
                        return component
                    else:
                        raise Exception("pyissm.model.mesh.bamg: if a spatial component is a list, its elements must be of type dict or OrderedDict")
                else:
                    return component
            
            ## Single contour provided as a dict, return list with one element
            elif isinstance(component, (dict, collections.OrderedDict)):
                return [component]
            else:
                raise Exception("pyissm.model.mesh.bamg: spatial components must be a filename (str), a list of contours (list of dicts), or a single contour (dict).")
        return []
    
    def _process_spatial_component(component, domain, ref_counter, count, is_hole = False, is_subdomain = False):
        """
        Process a spatial component (domain, holes, subdomains) for BAMG.
        """

        ## For each contour in the component
        for i in range(len(component)):

            ### Check the contour is closed
            if component[i]['x'][0] != component[i]['x'][-1] or component[i]['y'][0] != component[i]['y'][-1]:
                raise Exception("pyissm.model.mesh.bamg: each contour must be closed (first and last points must be identical).")

            ### Check all holes and subdomain contours are INSIDE the principle domain
            ## NOTE: This check differs from existing bamg.m and bamg.py which only checks the domain contours (not holes)
            if is_hole or is_subdomain:
                for c in component:
                    flags = tools.wrappers.ContourToNodes(c['x'], c['y'], [domain[0]], 0)[0]
                    if np.any(np.logical_not(flags)):
                        raise Exception("pyissm.model.mesh.bamg: all contours in 'holes' and 'subdomains' must be inside the first contour of 'domain'.")

            ### Check the orientation of the contour
            nods = component[i]['nods'] - 1
            test = np.sum((component[i]['x'][1:nods + 1] - component[i]['x'][0:nods]) * (component[i]['y'][1:nods + 1] + component[i]['y'][0:nods]))

            if (is_hole and test < 0) or (is_subdomain and test > 0):
                ## TODO: Why is subdomain orientation opposite to holes?
                print("At least one contour in 'holes' or 'subdomains' was not correctly oriented and has been re-oriented")
                component[i]['x'] = np.flipud(component[i]['x'])
                component[i]['y'] = np.flipud(component[i]['y'])

            ## If processing domain, add info to bamg_geom
            if not is_hole and not is_subdomain:
                ### Flag how many edges we have so far
                edge_length = len(bamg_geom['Edges'])

                ### Add all points to bamg_geom
                bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'], np.vstack((component[i]['x'][0:nods], component[i]['y'][0:nods], np.ones((nods)))).T))
                bamg_geom['Edges'] = np.vstack((bamg_geom['Edges'], np.vstack((np.arange(count + 1, count + nods + 1), np.hstack((np.arange(count + 2, count + nods + 1), count + 1)), 1. * np.ones((nods)))).T))

                new_edge_length = len(bamg_geom['Edges'])
                edges_required = np.asarray(range((edge_length + 1), (new_edge_length + 1)))  # NOTE: Upper bound of range is non-inclusive (compare to src/m/mesh/bamg.m)
                if i > 0:
                    bamg_geom['SubDomains'] = np.vstack((bamg_geom['SubDomains'], [2, count + 1, 1, -ref_counter]))
                    ref_counter = ref_counter + 1
                else:
                    bamg_geom['SubDomains'] = np.vstack((bamg_geom['SubDomains'], [2, count + 1, 1, 0]))

            ## If processing holes, add info to bamg_geom
            ## TODO: Why negative ref_counter for holes?
            if is_hole:
                bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'], np.vstack((holes[i]['x'][0:nods], holes[i]['y'][0:nods], np.ones((nods)))).T))
                bamg_geom['Edges'] = np.vstack((bamg_geom['Edges'], np.vstack((np.arange(count + 1, count + nods + 1), np.hstack((np.arange(count + 2, count + nods + 1), count + 1)), 1. * np.ones((nods)))).T))
                bamg_geom['SubDomains'] = np.vstack((bamg_geom['SubDomains'], [2, count + 1, 1, -ref_counter]))

                ref_counter = ref_counter + 1
            
            ## If processing subdomains, add info to bamg_geom
            ## TODO: Why positive ref_counter for subdomains?
            if is_subdomain:
                bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'], np.vstack((subdomains[i]['x'][0:nods], subdomains[i]['y'][0:nods], np.ones((nods)))).T))
                bamg_geom['Edges'] = np.vstack((bamg_geom['Edges'], np.vstack((np.arange(count + 1, count + nods + 1), np.hstack((np.arange(count + 2, count + nods + 1), count + 1)), 1. * np.ones((nods)))).T))
                bamg_geom['SubDomains'] = np.vstack((bamg_geom['SubDomains'], [2, count + 1, 1, ref_counter]))

                ref_counter = ref_counter + 1

            # Update counter
            count += nods

        if not (is_hole or is_subdomain):
            return component, ref_counter, count, bamg_geom, edges_required
        else:
            return component, ref_counter, count, bamg_geom
    
    def _seg_intersect(seg1, seg2):
        """
        Check if 2 segments intersect.
        NOTE: Taken directly from $ISSM_DIR/src/m/geometry/SegIntersect.py
        """

        bval = 1

        xA = seg1[0, 0]
        yA = seg1[0, 1]
        xB = seg1[1, 0]
        yB = seg1[1, 1]
        xC = seg2[0, 0]
        yC = seg2[0, 1]
        xD = seg2[1, 0]
        yD = seg2[1, 1]

        O2A = np.array([xA, yA]) - np.array([xD / 2. + xC / 2., yD / 2. + yC / 2.])
        O2B = np.array([xB, yB]) - np.array([xD / 2. + xC / 2., yD / 2. + yC / 2.])
        O1C = np.array([xC, yC]) - np.array([xA / 2. + xB / 2., yB / 2. + yA / 2.])
        O1D = np.array([xD, yD]) - np.array([xA / 2. + xB / 2., yB / 2. + yA / 2.])

        n1 = np.array([yA - yB, xB - xA])  #normal vector to segA
        n2 = np.array([yC - yD, xD - xC])  #normal vector to segB

        test1 = np.dot(n2, O2A)
        test2 = np.dot(n2, O2B)

        if test1 * test2 > 0:
            bval = 0
            return bval

        test3 = np.dot(n1, O1C)
        test4 = np.dot(n1, O1D)

        if test3 * test4 > 0:
            bval = 0
            return bval

        #if colinear
        if test1 * test2 == 0 and test3 * test4 == 0 and np.linalg.det(np.hstack((n1.reshape((-1, )), n2.reshape(-1, )))) == 0:

            #projection on the axis O1O2
            O2O1 = np.array([xA / 2. + xB / 2., yB / 2. + yA / 2.]) - np.array([xD / 2. + xC / 2., yD / 2. + yC / 2.])
            O1A = np.dot(O2O1, (O2A - O2O1))
            O1B = np.dot(O2O1, (O2B - O2O1))
            O1C = np.dot(O2O1, O1C)
            O1D = np.dot(O2O1, O1D)

        #test if one point is included in the other segment (-> bval = 1)
            if (O1C - O1A) * (O1D - O1A) < 0:
                bval = 1
                return bval
            if (O1C - O1B) * (O1D - O1B) < 0:
                bval = 1
                return bval
            if (O1A - O1C) * (O1B - O1C) < 0:
                bval = 1
                return bval
            if (O1A - O1D) * (O1B - O1D) < 0:
                bval = 1
                return bval

        #test if the 2 segments have the same middle (-> bval = 1)
            if O2O1 == 0:
                bval = 1
                return bval

        #else
            bval = 0
            return bval

        return bval

    def _process_rifts(rift_file):
        """
        Process rifts for BAMG.
        """
        
        # Error checks
        if not isinstance(rift_file, (str)):
            raise TypeError("pyissm.model.mesh.bamg: rifts must be a filename (str).")
        if not os.path.exists(rift_file):
            raise IOError(f"pyissm.model.mesh.bamg: rift file {rift_file} does not exist.")

        # Read rift file
        rift = tools.exp.exp_read(rift_file)

        # Process each rift
        for i in range(len(rift)):

            ## Check whether all points of the rift are inside the domain
            flags = tools.wrappers.ContourToNodes(rift[i]['x'], rift[i]['y'], [domain[0]], 0)[0]
            if np.all(np.logical_not(flags)):
                raise RuntimeError("pyissm.model.mesh.bamg: one rift has all its points outside of the domain outline")
            
            ## Check if rift tip is outside of the domain
            elif np.any(np.logical_not(flags)):
                # We have LOTS of work to do
                print('Rift tip outside of or on the domain has been detected and is being processed...')

                ## Check that only one point is outside (for now)
                if np.sum(np.logical_not(flags).astype(int)) != 1:
                    raise RuntimeError("pyissm.model.mesh.bamg: only one point outside of the domain is supported at this time")

                ## Move tip outside to the first position
                if not flags[0]:
                    # OK, first point is outside (do nothing),
                    pass
                elif not flags[-1]:
                    rift[i]['x'] = np.flipud(rift[i]['x'])
                    rift[i]['y'] = np.flipud(rift[i]['y'])
                else:
                    raise RuntimeError('pyissm.model.mesh.bamg: only a rift tip can be outside of the domain')

                # Get coordinate of intersection point
                x1 = rift[i]['x'][0]
                y1 = rift[i]['y'][0]
                x2 = rift[i]['x'][1]
                y2 = rift[i]['y'][1]
                for j in range(0, np.size(domain[0]['x']) - 1):
                    if _seg_intersect(np.array([[x1, y1], [x2, y2]]), np.array([[domain[0]['x'][j], domain[0]['y'][j]], [domain[0]['x'][j + 1], domain[0]['y'][j + 1]]])):

                        # Get position of the two nodes of the edge in domain
                        i1 = j
                        i2 = j + 1

                        # Rift is crossing edge [i1, i2] of the domain
                        # Get coordinate of intersection point (http://mathworld.wolfram.com/Line-LineIntersection.html)
                        x3 = domain[0]['x'][i1]
                        y3 = domain[0]['y'][i1]
                        x4 = domain[0]['x'][i2]
                        y4 = domain[0]['y'][i2]
                        x = np.linalg.det(np.array([[np.linalg.det(np.array([[x1, y1], [x2, y2]])), x1 - x2], [np.linalg.det(np.array([[x3, y3], [x4, y4]])), x3 - x4]])) / np.linalg.det(np.array([[x1 - x2, y1 - y2], [x3 - x4, y3 - y4]]))
                        y = np.linalg.det(np.array([[np.linalg.det(np.array([[x1, y1], [x2, y2]])), y1 - y2], [np.linalg.det(np.array([[x3, y3], [x4, y4]])), y3 - y4]])) / np.linalg.det(np.array([[x1 - x2, y1 - y2], [x3 - x4, y3 - y4]]))

                        segdis = np.sqrt((x4 - x3)**2 + (y4 - y3)**2)
                        tipdis = np.array([np.sqrt((x - x3)**2 + (y - y3)**2), np.sqrt((x - x4)**2 + (y - y4)**2)])

                        if np.min(tipdis) / segdis < options['toltip']:
                            print('moving tip-domain intersection point')

                            # Get position of the closer point
                            if tipdis[0] > tipdis[1]:
                                pos = i2
                            else:
                                pos = i1

                            # This point is only in Vertices (number pos).
                            # OK, now we can add our own rift
                            nods = rift[i]['nods'] - 1
                            bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'], np.hstack((rift[i]['x'][1:].reshape(-1, ), rift[i]['y'][1:].reshape(-1, ), np.ones((nods, 1))))))
                            bamg_geom['Edges'] = np.vstack((
                                bamg_geom['Edges'],
                                np.array([[pos, count + 1, (1 + i)]]),
                                np.hstack((np.arange(count + 1, count + nods).reshape(-1, ), np.arange(count + 2, count + nods + 1).reshape(-1, ), (1 + i) * np.ones((nods - 1, 1))))
                            ))
                            count += nods
                            break
                        else:
                            # Add intersection point to Vertices
                            bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'],
                                np.array([[x, y, 1]])
                            ))
                            count += 1

                            # Decompose the crossing edge into 2 subedges
                            pos = np.nonzero(np.logical_and(bamg_geom['Edges'][:, 0] == i1, bamg_geom['Edges'][:, 1] == i2))[0]
                            if not pos:
                                raise RuntimeError('pyissm.model.mesh.bamg: a problem occurred...')
                            bamg_geom['Edges'] = np.vstack((
                                bamg_geom['Edges'][0:pos - 1, :],
                                np.array([[
                                    bamg_geom['Edges'][pos, 0],
                                    count,
                                    bamg_geom['Edges'][pos, 2]
                                ]]),
                                np.array([[
                                    count,
                                    bamg_geom['Edges'][pos, 1],
                                    bamg_geom['Edges'][pos, 2]
                                ]]),
                                bamg_geom['Edges'][pos + 1:, :]
                            ))

                            # OK, now we can add our own rift
                            nods = rift[i]['nods'] - 1
                            bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'],
                                np.hstack((
                                    rift[i]['x'][1:].reshape(-1, ),
                                    rift[i]['y'][1:].reshape(-1, ),
                                    np.ones((nods, 1))
                                ))
                            ))
                            bamg_geom['Edges'] = np.vstack((
                                bamg_geom['Edges'],
                                np.array([[count, count + 1, 2]]),
                                np.hstack((
                                    np.arange(count + 1, count + nods).reshape(-1, ),
                                    np.arange(count + 2, count + nods + 1).reshape(-1, ),
                                    (1 + i) * np.ones((nods - 1, 1))
                                ))
                            ))
                            count += nods
                            break
            else:
                nods = rift[i]['nods'] - 1
                bamg_geom['Vertices'] = np.vstack((
                    bamg_geom['Vertices'],
                    np.hstack((
                        rift[i]['x'][:],
                        rift[i]['y'][:],
                        np.ones((nods + 1, 1))
                    ))
                ))
                bamg_geom['Edges'] = np.vstack((
                    bamg_geom['Edges'],
                    np.hstack((
                        np.arange(count + 1, count + nods).reshape(-1, ),
                        np.arange(count + 2, count + nods + 1).reshape(-1, ),
                        i * np.ones((nods, 1))
                    ))
                ))
                count += (nods + 1)

        return count, bamg_geom
    
    def _process_tracks(track):
        """
        Process tracks for BAMG.
        """
                
        # Read tracks
        if all(isinstance(track, str)):
            track = tools.exp.expread(track)
            track = np.hstack((track.x.reshape(-1, ), track.y.reshape(-1, )))
        else:
            track = float(track)

        if np.size(track, axis=1) == 2:
            track = np.hstack((track, 3. * np.ones((np.size(track, axis=0), 1))))

        # Only keep those inside
        flags = tools.wrappers.ContourToNodes(track[:, 0], track[:, 1], [domain[0]], 0)[0]
        track = track[np.nonzero(flags), :]

        # Add all points to bamg_geometry
        nods = np.size(track, axis=0)
        bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'], track))
        bamg_geom['Edges'] = np.vstack((
            bamg_geom['Edges'],
            np.hstack((
                np.arange(count + 1, count + nods).reshape(-1, ),
                np.arange(count + 2, count + nods + 1).reshape(-1, ),
                3. * np.ones((nods - 1, 1))
            ))
        ))

        # Update counter
        count += nods

        return count, bamg_geom
    
    def _process_required_vertices(required_vertices):
        """
        Process required vertices for BAMG.
        """

        if np.size(required_vertices, axis=1) == 2:
            required_vertices = np.hstack((required_vertices, 4. * np.ones((np.size(required_vertices, axis=0), 1))))

        # Only keep those inside
        flags = tools.wrappers.ContourToNodes(required_vertices[:, 0], required_vertices[:, 1], [domain[0]], 0)[0]
        required_vertices = required_vertices[np.nonzero(flags)[0], :]

        # Add all points to bamg_geom
        nods = np.size(required_vertices, axis=0)
        bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'], required_vertices))

        # Update counter
        count += nods

        return count, bamg_geom
    
    # ---------------------------------------------------------------
    
    # Update defaults with user-specified options
    options = collections.OrderedDict(defaults)
    options.update(kwargs)

    # Initialize geometry and mesh structures
    bamg_geom = _bamg_geom()
    bamg_mesh = _bamg_mesh()

    # Initialise counters
    subdomain_ref = 1
    hole_ref = 1

    # Build BAMG mesh from domain, holes, subdomains
    if options['domain'] is not None:
        domain_file = options['domain']
        domain = _load_spatial_components(domain_file)

        holes = []
        if options['holes'] is not None:
            hole_file = options['holes']
            holes = _load_spatial_components(hole_file)

        subdomains = []
        if options['subdomains'] is not None:
            subdomain_file = options['subdomains']
            subdomains = _load_spatial_components(subdomain_file)

        # Build geometry
        count = 0

        ## Process domain
        domain, subdomain_ref, count, bamg_geom, edges_required = _process_spatial_component(domain, domain, subdomain_ref, count, is_hole = False, is_subdomain = False)

        ## Process holes
        holes, hole_ref, count, bamg_geom = _process_spatial_component(holes, domain, hole_ref, count, is_hole = True, is_subdomain = False)

        ## Process subdomains
        subdomains, subdomain_ref, count, bamg_geom = _process_spatial_component(subdomains, domain, subdomain_ref, count, is_hole = False, is_subdomain = True)

        # Process vertical options
        if options['vertical'] == 1:
            if np.size(options['Markers']) != np.size(bamg_geom['Edges'], 0):
                edges_size = np.size(bamg_geom['Edges'], 0)
                raise RuntimeError(f'for 2d vertical mesh, \'Markers\' option is required, and should be of size {edges_size}')
        
        if np.size(options['Markers']) == np.size(bamg_geom['Edges'], 0):
            bamg_geom['Edges'][:, 2] = options['Markers']

        # Process rifts
        if options['rifts'] is not None:
            count, bamg_geom = _process_rifts(options['rifts'])

        # Process tracks
        if options['tracks'] is not None:
            count, bamg_geom = _process_tracks(options['tracks'])

        # Proces required vertices
        if options['RequiredVertices'] is not None:
            count, bamg_geom = _process_required_vertices(options['RequiredVertices'])

        # Process RequiredEdges
        if options['NoBoundaryRefinement'] == 1:
            bamg_geom['RequiredEdges'] = edges_required
        elif options['NoBoundaryRefinementAllBoundaries'] == 1:
                bamg_geom['RequiredEdges'] = np.arange(1, bamg_geom['Edges'].shape[0]).T
    
    # If a geometry is already provided, use it
    elif isinstance(md.private.bamg, dict) and 'geometry' in md.private.bamg:
        bamg_geom = _bamg_geom(**md.private.bamg['geometry'])
    else:
        # Do nothing...
        pass

    # If domain is not specified, check for existing mesh
    if md.mesh.numberofvertices and md.mesh.element_type() == 'Tria':
        ## If there is an existing BAMG mesh, use it
        if isinstance(md.private.bamg, dict) and 'mesh' in md.private.bamg:
            bamg_mesh = _bamg_mesh(**md.private.bamg['mesh'])
        else:
            ## If there is an existing non-BAMG mesh, convert it
            bamg_mesh['Vertices'] = np.vstack((
                md.mesh.x,
                md.mesh.y,
                np.ones((md.mesh.numberofvertices))
            )).T
            bamg_mesh['Triangles'] = np.hstack((md.mesh.elements, np.ones((md.mesh.numberofelements, 1))))

        ## If there are rifts in the model, raise an error (not supported yet)
        if isinstance(md.rifts.riftstruct, dict):
            raise TypeError('pyissm.model.mesh.bamg: rifts not supported yet. Do meshprocessrift after bamg.')
    
    # Call the BAMG mesher
    bamg_mesh_out, bamg_geom_out = tools.wrappers.BamgMesher(bamg_mesh, bamg_geom, options)

    # Populate md structure
    if options['vertical'] == 1:
        ## Create 2D vertical mesh
        md.mesh = model.classes.mesh.mesh2dvertical()
        md.mesh.x = bamg_mesh_out['Vertices'][:, 0].copy()
        md.mesh.y = bamg_mesh_out['Vertices'][:, 1].copy()
        md.mesh.elements = bamg_mesh_out['Triangles'][:, 0:3].astype(int)
        md.mesh.edges = bamg_mesh_out['IssmEdges'].astype(int)
        md.mesh.segments = bamg_mesh_out['IssmSegments'][:, 0:3].astype(int)
        md.mesh.segmentmarkers = bamg_mesh_out['IssmSegments'][:, 3].astype(int)

        md.mesh.numberofelements = np.size(md.mesh.elements, axis=0)
        md.mesh.numberofvertices = np.size(md.mesh.x)
        md.mesh.numberofedges = np.size(md.mesh.edges, axis=0)
        md.mesh.vertexonboundary = np.zeros(md.mesh.numberofvertices, int)
        md.mesh.vertexonboundary[md.mesh.segments[:, 0:2] - 1] = 1

    elif options['3dsurface'] == 1:
        md.mesh = model.classes.mesh.mesh3dsurface()
        md.mesh.x = bamg_mesh_out['Vertices'][:, 0].copy()
        md.mesh.y = bamg_mesh_out['Vertices'][:, 1].copy()
        md.mesh.z = md.mesh.x
        md.mesh.z[:] = 0
        md.mesh.elements = bamg_mesh_out['Triangles'][:, 0:3].astype(int)
        md.mesh.edges = bamg_mesh_out['IssmEdges'].astype(int)
        md.mesh.segments = bamg_mesh_out['IssmSegments'][:, 0:3].astype(int)
        md.mesh.segmentmarkers = bamg_mesh_out['IssmSegments'][:, 3].astype(int)

        # Fill in rest of fields
        md.mesh.numberofelements = np.size(md.mesh.elements, axis=0)
        md.mesh.numberofvertices = np.size(md.mesh.x)
        md.mesh.numberofedges = np.size(md.mesh.edges, axis=0)
        md.mesh.vertexonboundary = np.zeros(md.mesh.numberofvertices, int)
        md.mesh.vertexonboundary[md.mesh.segments[:, 0:2] - 1] = 1

    else:
        md.mesh = model.classes.mesh.mesh2d()
        md.mesh.x = bamg_mesh_out['Vertices'][:, 0].copy()
        md.mesh.y = bamg_mesh_out['Vertices'][:, 1].copy()
        md.mesh.elements = bamg_mesh_out['Triangles'][:, 0:3].astype(int)
        md.mesh.edges = bamg_mesh_out['IssmEdges'].astype(int)
        md.mesh.segments = bamg_mesh_out['IssmSegments'][:, 0:3].astype(int)
        md.mesh.segmentmarkers = bamg_mesh_out['IssmSegments'][:, 3].astype(int)

        # Fill in rest of fields
        md.mesh.numberofelements = np.size(md.mesh.elements, axis=0)
        md.mesh.numberofvertices = np.size(md.mesh.x)
        md.mesh.numberofedges = np.size(md.mesh.edges, axis=0)
        md.mesh.vertexonboundary = np.zeros(md.mesh.numberofvertices, int)
        md.mesh.vertexonboundary[md.mesh.segments[:, 0:2] - 1] = 1
    

    # BAMG private fields
    md.private.bamg = collections.OrderedDict()
    md.private.bamg['mesh'] = _bamg_mesh(**bamg_mesh_out)
    md.private.bamg['geometry'] = _bamg_geom(**bamg_geom_out)
    md.mesh.elementconnectivity = md.private.bamg['mesh']['ElementConnectivity']
    md.mesh.elementconnectivity[np.nonzero(np.isnan(md.mesh.elementconnectivity))] = 0
    md.mesh.elementconnectivity = md.mesh.elementconnectivity.astype(int)

    # Check for orphan vertices
    nodes = np.arange(1, md.mesh.numberofvertices + 1)
    if np.any(~np.isin(nodes, md.mesh.elements)):
        raise RuntimeError('Output mesh has orphans. Check your Domain and/or RequiredVertices')

    return md

def bamg_flowband(md,
                 x,
                 surf,
                 base,
                 **kwargs):
    
    """
    Create a flowband mesh using BAMG (Bidimensional Anisotropic Mesh Generator).

    This function generates a triangular mesh for a flowband (vertical 2D slice) using 
    the BAMG mesh generator. The flowband is defined by surface and base profiles 
    along a specified coordinate path, creating a vertical mesh suitable for ice flow 
    modeling in the vertical plane.

    Parameters
    ----------
    md : object
        ISSM model object whose mesh fields will be populated with the generated flowband mesh.
    x : array_like
        1D array of coordinates along the flowband path. These represent the horizontal
        positions where the surface and base elevations are defined.
    surf : array_like
        1D array of surface elevations corresponding to each x coordinate. Must have
        the same length as x.
    base : array_like
        1D array of base (bed) elevations corresponding to each x coordinate. Must have
        the same length as x.
    **kwargs : dict, optional
        Additional keyword arguments passed to the bamg function. See bamg() documentation
        for supported options.

    Returns
    -------
    md : object
        A new ISSM model object with a vertical 2D mesh populated, including:
        - mesh.x, mesh.y: Node coordinates in the flowband coordinate system
        - mesh.elements: Element connectivity matrix for triangular elements
        - mesh.edges: Edge connectivity matrix
        - mesh.segments: Boundary segment definitions with markers
        - mesh.segmentmarkers: Boundary segment markers (1=base, 2=right, 3=surface, 4=left)
        - mesh.numberofvertices: Total number of mesh vertices
        - mesh.numberofelements: Total number of mesh elements
        - mesh.numberofedges: Total number of mesh edges
        - mesh.vertexonboundary: Boolean array indicating boundary vertices
        - mesh.vertexonbase: Boolean array indicating vertices on the base boundary
        - mesh.vertexonsurface: Boolean array indicating vertices on the surface boundary
        - mesh.elementconnectivity: Element-to-element connectivity

    Raises
    ------
    ValueError
        If x, surf, and base arrays do not have the same length.
    RuntimeError
        If BAMG mesh generation fails or if incompatible meshing options are specified.

    Notes
    -----
    This function creates a vertical 2D mesh by:
    1. Constructing a closed domain from the surface and base profiles
    2. Assigning boundary markers: 1=base, 2=right side, 3=surface, 4=left side
    3. Calling the BAMG mesh generator with vertical=1 option (to convert to 2D vertical mesh)
    4. Post-processing to identify vertices on base and surface boundaries

    The resulting mesh is suitable for flowband modeling where ice flow is assumed
    to be primarily in the vertical plane defined by the x-coordinate path.

    Examples
    --------
    >>> md = pyissm.Model()
    >>> x = np.arange(1, 3001, 100).T
    >>> h = np.linspace(1000, 300, np.size(x)).T
    >>> b = -917. / 1023. * h
    >>> md = pyissm.model.mesh.bamg_flowband(md, x = x, surf = b + h, base = b, hmax = 80.)
    """

    # Create domain structure
    domain = collections.OrderedDict()
    domain['x'] = np.concatenate((x, np.flipud(x), [x[0]]))
    domain['y'] = np.concatenate((base, np.flipud(surf), [base[0]]))
    domain['nods'] = np.size(domain['x'])

    # Create markers (base, right side, top surface, left side)
    m = np.ones((np.size(domain['x']) - 1, ))
    m[np.size(x) - 1] = 2
    m[np.size(x):2 * np.size(x) - 1] = 3
    m[2 * np.size(x) - 1] = 4

    # Call bamg
    md = model.Model()
    md = bamg(md, domain = [domain], Markers = m, vertical = 1, **kwargs)

    # Deal with vertices on bed
    md.mesh.vertexonbase = np.zeros((md.mesh.numberofvertices, ))
    md.mesh.vertexonbase[np.where(md.mesh.vertex_flags(1))] = 1
    md.mesh.vertexonsurface = np.zeros((md.mesh.numberofvertices, ))
    md.mesh.vertexonsurface[np.where(md.mesh.vertex_flags(3))] = 1

    return md

def mesh_convert(md, **kwargs):
    """
    Convert mesh to BAMG format for advanced mesh operations.

    This function converts an existing mesh to BAMG (Bidimensional Anisotropic Mesh Generator)
    format, enabling access to BAMG's advanced mesh manipulation capabilities. The conversion
    creates internal BAMG data structures while preserving the original mesh geometry and
    connectivity.

    Parameters
    ----------
    md : object
        ISSM model object containing the mesh to be converted. The mesh should have valid
        elements, coordinates, and connectivity information.
    **kwargs : dict, optional
        Additional keyword arguments to customize the conversion:
        - index : array_like, optional
            Element connectivity matrix. Defaults to md.mesh.elements.
        - x : array_like, optional
            X-coordinates of mesh vertices. Defaults to md.mesh.x.
        - y : array_like, optional
            Y-coordinates of mesh vertices. Defaults to md.mesh.y.

    Returns
    -------
    md : object
        The input ISSM model object with updated mesh properties and BAMG data structures:
        - mesh.x, mesh.y: Node coordinates
        - mesh.elements: Element connectivity matrix
        - mesh.edges: Edge connectivity matrix
        - mesh.segments: Boundary segment definitions
        - mesh.segmentmarkers: Boundary segment markers
        - mesh.numberofvertices: Total number of mesh vertices
        - mesh.numberofelements: Total number of mesh elements
        - mesh.numberofedges: Total number of mesh edges
        - mesh.vertexonboundary: Boolean array indicating boundary vertices
        - mesh.elementconnectivity: Element-to-element connectivity
        - private.bamg: BAMG-specific mesh and geometry data structures

    Notes
    -----
    This function is primarily used to prepare existing meshes for advanced BAMG operations
    such as mesh adaptation, refinement, or anisotropic meshing. The conversion creates
    internal BAMG data structures that enable seamless integration with other BAMG-based
    mesh operations.

    The function preserves all original mesh properties while adding BAMG-specific data
    structures to the model's private fields. This allows subsequent BAMG operations
    to work efficiently without data conversion overhead.

    Examples
    --------
    >>> md = pyissm.Model()
    >>> md = pyissm.model.mesh.triangle(md, 'domain.exp', 1000.0)
    >>> md = pyissm.model.mesh.meshconvert(md)
    >>> md = pyissm.model.mesh.meshconvert(md, x=custom_x, y=custom_y)
    """

    # Default arguments
    options = {
        'index': md.mesh.elements,
        'x': md.mesh.x,
        'y': md.mesh.y,
    }

    # Update options with any user-specified arguments
    options.update(kwargs)

    # Call the BAMG mesh converter
    bamg_mesh_out, bamg_geom_out = tools.wrappers.BamgConvertMesh(options['index'],
                                                                  options['x'],
                                                                  options['y'])
    
    # Populate md structure
    md.private.bamg = collections.OrderedDict()
    md.private.bamg['mesh'] = _bamg_mesh(**bamg_mesh_out)
    md.private.bamg['geometry'] = _bamg_geom(**bamg_geom_out)
    md.mesh = model.classes.mesh.mesh2d()
    md.mesh.x = bamg_mesh_out['Vertices'][:, 0].copy()
    md.mesh.y = bamg_mesh_out['Vertices'][:, 1].copy()
    md.mesh.elements = bamg_mesh_out['Triangles'][:, 0:3].astype(int)
    md.mesh.edges = bamg_mesh_out['IssmEdges'].astype(int)
    md.mesh.segments = bamg_mesh_out['IssmSegments'][:, 0:3].astype(int)
    md.mesh.segmentmarkers = bamg_mesh_out['IssmSegments'][:, 3].astype(int)

    md.mesh.numberofelements = np.size(md.mesh.elements, axis=0)
    md.mesh.numberofvertices = np.size(md.mesh.x)
    md.mesh.numberofedges = np.size(md.mesh.edges, axis=0)
    md.mesh.vertexonboundary = np.zeros(md.mesh.numberofvertices, int)
    md.mesh.vertexonboundary[md.mesh.segments[:, 0:1] - 1] = 1
    md.mesh.elementconnectivity = md.private.bamg['mesh']['ElementConnectivity']
    md.mesh.elementconnectivity[np.where(np.isnan(md.mesh.elementconnectivity))[0]] = 0

    return md

def model_intersect_3d():
    """
    Intersect a 3D model with a plane defined by points (xs, ys, zs).

    Raises
    ------
    NotImplementedError
        Function is not yet implemented.
    """

    raise NotImplementedError('pyissm.model.mesh.model_intersect_3d:  This functionality is not yet implemented. Please contact ACCESS-NRI for support.')

def model_merge_3d():
    """
    Merge two 3D models into a single model.

    Raises
    ------
    NotImplementedError
        Function is not yet implemented.
    """

    raise NotImplementedError('pyissm.model.mesh.model_merge_3d:  This functionality is not yet implemented. Please contact ACCESS-NRI for support.')

def twod_to_3d():
    """
    Convert 2D mesh to 3D surface mesh.

    Raises
    ------
    NotImplementedError
        Function is not yet implemented.
    """

    raise NotImplementedError('pyissm.model.mesh.twod_to_3d:  This functionality is not yet implemented. Please contact ACCESS-NRI for support.')

def project_3d(md,
               vector,
               type = 'node',
               layer = 0,
               padding = 0,
               degree = 0.0):
    """
    Vertically project a vector from 2D mesh into a 3D mesh.

    This function extrudes a 2D vector (defined on nodes or elements) into
    a 3D mesh by replicating it across layers. The vector can be extruded
    uniformly across all layers or applied to a specific layer with optional
    polynomial interpolation from bottom to top.

    Parameters
    ----------
    md : ISSM Model object
        ISSM Model object containing a 3D mesh. Must have valid 3D mesh
        attributes (e.g., md.mesh.numberofvertices, md.mesh.numberoflayers).
    vector : ndarray
        2D vector to be projected into 3D. Can be:
        - (md.mesh.numberofvertices2d,) for node-based vector
        - (md.mesh.numberofvertices2d + 1,) for node vector with extra element
        - (md.mesh.numberofelements2d,) for element-based vector
        - (md.mesh.numberofelements2d + 1,) for element vector with extra element
        - 2D array with shape (n, m) for multiple fields
    type : str, optional
        Type of vector projection. Options are:
        - 'node': Project node-based vector (default)
        - 'element': Project element-based vector
        - 'poly': Project with polynomial interpolation from bottom to top
    layer : int, optional
        Layer number (1-indexed) where vector should keep its values. If 0 (default),
        all layers adopt the value of the 2D vector. If specified, only the
        designated layer is populated; other layers are set to padding value.
    padding : float, optional
        Value used to fill 3D layers not being directly projected.
        Default is 0.
    degree : float, optional
        Degree of polynomial used for extrusion when type='poly'.
        Controls the vertical variation of the vector from bottom to top.
        Default is 0.0 (uniform extrusion).

    Returns
    -------
    projected_vector : ndarray
        3D vector with shape:
        - (md.mesh.numberofvertices,) or (md.mesh.numberofvertices + 1,) for node type
        - (md.mesh.numberofelements,) or (md.mesh.numberofelements + 1,) for element type
        - (n, m) preserved for input with multiple fields
        Data type is preserved from input vector.

    Raises
    ------
    TypeError
        If md is not provided or does not contain a 3D mesh.
    TypeError
        If vector size does not match expected dimensions for the specified type.
    TypeError
        If type is not one of 'node', 'element', or 'poly'.

    Examples
    --------
    Project a 2D node vector uniformly across all layers:

    >>> extruded_vector = project_3d(md, vector=vector2d, type='node')

    Project a 2D vector to a specific layer:

    >>> extruded_vector = project_3d(md, vector=vector2d, type='node', layer=1)

    Project with polynomial interpolation:

    >>> extruded_vector = project_3d(md, vector=vector2d, type='poly', degree=2.0)

    Notes
    -----
    The extra element in vector arrays (e.g., md.mesh.numberofvertices2d + 1)
    is preserved in the output as the last element.
    """

    ## NOTE: This function is taken directly from $ISSM_DIR/src/m/extrusion/project3d.py with only minor modifications for pyISSM integration.

    # Error checks
    if not md:
        raise TypeError("pyissm.model.mesh.project_3d: md must be provided")
    if md.mesh.domain_type().lower() != '3d':
        raise TypeError("pyissm.model.mesh.project_3d: md must contain a 3D mesh")

    #Handle special case where vector2d is single element (differs from representation in MATLAB)
    if isinstance(vector, (bool, int, float)):
        projected_vector = vector

    if np.size(vector) == 1:
        projected_vector = vector

    elif type.lower() == 'node':
        #Initialize 3d vector
        if np.ndim(vector) == 1:
            if vector.shape[0] == md.mesh.numberofvertices2d:
                projected_vector = (padding * np.ones((md.mesh.numberofvertices))).astype(vector.dtype)
            elif vector.shape[0] == md.mesh.numberofvertices2d + 1:
                projected_vector = (padding * np.ones((md.mesh.numberofvertices + 1))).astype(vector.dtype)
                projected_vector[-1] = vector[-1]
                vector = vector[:-1]
            else:
                raise TypeError("pyissm.model.mesh.project_3d: Vector must be the length of md.mesh.numberofvertices2d or md.mesh.numberofvertices2d + 1")
            #Fill in
            if layer == 0:
                for i in range(md.mesh.numberoflayers):
                    projected_vector[(i * md.mesh.numberofvertices2d):((i + 1) * md.mesh.numberofvertices2d)] = vector
            else:
                projected_vector[((layer - 1) * md.mesh.numberofvertices2d):(layer * md.mesh.numberofvertices2d)] = vector
        else:
            if vector.shape[0] == md.mesh.numberofvertices2d:
                projected_vector = (padding * np.ones((md.mesh.numberofvertices, np.size(vector, axis=1)))).astype(vector.dtype)
            elif vector.shape[0] == md.mesh.numberofvertices2d + 1:
                projected_vector = (padding * np.ones((md.mesh.numberofvertices + 1, np.size(vector, axis=1)))).astype(vector.dtype)
                projected_vector[-1, :] = vector[-1, :]
                vector = vector[:-1, :]
            else:
                raise TypeError("pyissm.model.mesh.project_3d: Vector must be the length of md.mesh.numberofvertices2d or md.mesh.numberofvertices2d + 1")
            #Fill in
            if layer == 0:
                for i in range(md.mesh.numberoflayers):
                    projected_vector[(i * md.mesh.numberofvertices2d):((i + 1) * md.mesh.numberofvertices2d), :] = vector
            else:
                projected_vector[((layer - 1) * md.mesh.numberofvertices2d):(layer * md.mesh.numberofvertices2d), :] = vector

    elif type.lower() == 'element':
        #Initialize 3d vector
        if np.ndim(vector) == 1:
            if vector.shape[0] == md.mesh.numberofelements2d:
                projected_vector = (padding * np.ones((md.mesh.numberofelements))).astype(vector.dtype)
            elif vector.shape[0] == md.mesh.numberofelements2d + 1:
                projected_vector = (padding * np.ones((md.mesh.numberofelements + 1))).astype(vector.dtype)
                projected_vector[-1] = vector[-1]
                vector = vector[:-1]
            else:
                raise TypeError("pyissm.model.mesh.project_3d: Vector must be the length of md.mesh.numberofelements2d or md.mesh.numberofelements2d + 1")
            #Fill in
            if layer == 0:
                for i in range(md.mesh.numberoflayers - 1):
                    projected_vector[(i * md.mesh.numberofelements2d):((i + 1) * md.mesh.numberofelements2d)] = vector
            else:
                projected_vector[((layer - 1) * md.mesh.numberofelements2d):(layer * md.mesh.numberofelements2d)] = vector
        else:
            if vector.shape[0] == md.mesh.numberofelements2d:
                projected_vector = (padding * np.ones((md.mesh.numberofelements, np.size(vector, axis=1)))).astype(vector.dtype)
            elif vector.shape[0] == md.mesh.numberofelements2d + 1:
                projected_vector = (padding * np.ones((md.mesh.numberofelements + 1, np.size(vector, axis=1)))).astype(vector.dtype)
                projected_vector[-1, :] = vector[-1, :]
                vector = vector[:-1, :]
            else:
                raise TypeError("pyissm.model.mesh.project_3d: Vector must be the length of md.mesh.numberofelements2d or md.mesh.numberofelements2d + 1")
            #Fill in
            if layer == 0:
                for i in range(md.mesh.numberoflayers - 1):
                    projected_vector[(i * md.mesh.numberofelements2d):((i + 1) * md.mesh.numberofelements2d), :] = vector
            else:
                projected_vector[((layer - 1) * md.mesh.numberofelements2d):(layer * md.mesh.numberofelements2d), :] = vector
    elif type.lower() == 'poly':
        #Initialize 3d vector
        if np.ndim(vector) == 1:
            if vector.shape[0] == md.mesh.numberofvertices2d:
                projected_vector = (padding * np.ones((md.mesh.numberofvertices))).astype(vector.dtype)
            elif vector.shape[0] == md.mesh.numberofvertices2d + 1:
                projected_vector = (padding * np.ones((md.mesh.numberofvertices + 1))).astype(vector.dtype)
                projected_vector[-1] = vector[-1]
                vector = vector[:-1]
            else:
                raise TypeError("pyissm.model.mesh.project_3d: Vector must be the length of md.mesh.numberofvertices2d or md.mesh.numberofvertices2d + 1")
            #Fill in
            if layer == 0:
                for i in range(md.mesh.numberoflayers - 1):
                    projected_vector[(i * md.mesh.numberofvertices2d):((i + 1) * md.mesh.numberofvertices2d)] = vector * (1.0 - (1.0 - i / (md.mesh.numberoflayers - 1.0))**degree)
            else:
                projected_vector[((layer - 1) * md.mesh.numberofvertices2d):(layer * md.mesh.numberofvertices2d)] = vector * (1.0 - (1.0 - layer / (md.mesh.numberoflayers - 1.0))**degree)
        else:
            if vector.shape[0] == md.mesh.numberofvertices2d:
                projected_vector = (padding * np.ones((md.mesh.numberofvertices, np.size(vector, axis=1)))).astype(vector.dtype)
            elif vector.shape[0] == md.mesh.numberofvertices2d + 1:
                projected_vector = (padding * np.ones((md.mesh.numberofvertices + 1, np.size(vector, axis=1)))).astype(vector.dtype)
                projected_vector[-1, :] = vector[-1, :]
                vector = vector[:-1, :]
            else:
                raise TypeError("pyissm.model.mesh.project_3d: Vector must be the length of md.mesh.numberofvertices2d or md.mesh.numberofvertices2d + 1")
            #Fill in
            if layer == 0:
                for i in range(md.mesh.numberoflayers - 1):
                    projected_vector[(i * md.mesh.numberofvertices2d):((i + 1) * md.mesh.numberofvertices2d), :] = vector * (1.0 - (1.0 - i / (md.mesh.numberoflayers - 1.0))**degree)
            else:
                projected_vector[((layer - 1) * md.mesh.numberofvertices2d):(layer * md.mesh.numberofvertices2d), :] = vector * (1.0 - (1.0 - layer / (md.mesh.numberoflayers - 1.0))**degree)
    else:
        raise TypeError("pyissm.model.mesh.project_3d: Unknown projection type")

    return projected_vector


def project_2d(md, value, layer):
    """
    Extract a 3D field value from a specified layer onto the 2D mesh.

    This function projects a field value from a given layer of an extruded 3D mesh
    back onto the original 2D mesh. It is commonly used to compare values across
    different layers or to analyze vertical variations in a 3D model.

    Parameters
    ----------
    md : ISSM Model object
        ISSM Model object containing a 3D mesh. Must have valid 3D mesh
        attributes (e.g., md.mesh.numberofvertices, md.mesh.numberoflayers).
    value : ndarray or array_like
        3D field values to be projected onto 2D. Can be:
        - (md.mesh.numberofvertices,) for node-based data
        - (md.mesh.numberofvertices + 1,) for node data with extra element
        - (md.mesh.numberofelements,) for element-based data
        - (md.mesh.numberofelements + 1,) for element data with extra element
        - 2D array with shape (n, 1) which will be reshaped to 1D
    layer : int
        Layer number (1-indexed) from which to extract data. Must be between
        1 and md.mesh.numberoflayers (where 1 is the base layer).

    Returns
    -------
    projection_value : ndarray
        2D field values extracted from the specified layer. Shape depends on input:
        - (md.mesh.numberofvertices2d,) for node-based input
        - (md.mesh.numberofvertices2d + 1,) for node input with extra element
        - (md.mesh.numberofelements2d,) for element-based input
        - (md.mesh.numberofelements2d + 1,) for element input with extra element

    Raises
    ------
    TypeError
        If md is not provided or does not contain a 3D mesh.
    ValueError
        If layer is not between 1 and md.mesh.numberoflayers.

    Examples
    --------
    Extract velocity from the second layer (1 = base):

    >>> vel_layer2 = project_2d(md3d, md3d.initialization.vel, 2)

    Extract temperature from the surface layer:

    >>> temp_surface = project_2d(md3d, md3d.thermal.temperature, md3d.mesh.numberoflayers)

    Notes
    -----
    Layer numbering is 1-indexed, where layer 1 corresponds to the base of the mesh
    and layer md.mesh.numberoflayers corresponds to the surface.

    This function is useful for:
    - Comparing field values across different vertical layers
    - Extracting surface or basal values from 3D simulations
    - Analyzing vertical variations in model results
    """
    
    ## NOTE: This function is taken directly from $ISSM_DIR/src/m/extrusion/project3d.py with only minor modifications for pyISSM integration.

    # Error checks
    if not md:
        raise TypeError("pyissm.model.mesh.project_2d: md must be provided")
    if md.mesh.domain_type().lower() != '3d':
        raise TypeError("pyissm.model.mesh.project_2d: md must contain a 3D mesh")

    if layer < 1 or layer > md.mesh.numberoflayers:
        raise ValueError(f"pyissm.model.mesh.project_2d: Layer must be between 1 and {md.mesh.numberoflayers}")

    # coerce to array in case float is passed
    if type(value) not in [np.ndarray, np.ma.core.MaskedArray]:
        value = np.array(value)

    # Check if 2D array
    vec2d = False
    if value.ndim == 2 and value.shape[1] == 1:
        value = value.reshape(-1, )
        vec2d = True

    # Project
    if value.size == 1:
        projection_value = value[(layer - 1) * md.mesh.numberofelements2d:layer * md.mesh.numberofelements2d]
    elif value.shape[0] == md.mesh.numberofvertices:
        projection_value = value[(layer - 1) * md.mesh.numberofvertices2d:layer * md.mesh.numberofvertices2d]
    elif value.shape[0] == md.mesh.numberofvertices + 1:
        if np.ndim(value) == 1:
            projection_value = np.hstack((value[(layer - 1) * md.mesh.numberofvertices2d:layer * md.mesh.numberofvertices2d], value[-1]))
        else:
            projection_value = np.vstack((value[(layer - 1) * md.mesh.numberofvertices2d:layer * md.mesh.numberofvertices2d], value[-1]))
    else:
        projection_value = value[(layer - 1) * md.mesh.numberofelements2d:layer * md.mesh.numberofelements2d]

    if vec2d:
        projection_value = projection_value.reshape(-1, )

    return projection_value

def depth_average(md, vector):
    """
    Compute depth average of a 3D vector using the trapezoidal rule.

    This function computes the depth-averaged value of a 3D field defined on an
    extruded mesh and returns the result projected onto the corresponding 2D mesh.
    The depth averaging is performed using the trapezoidal integration rule.

    Parameters
    ----------
    md : ISSM Model object
        ISSM Model object containing a 3D mesh. Must have valid 3D mesh
        attributes (e.g., md.mesh.numberofvertices, md.mesh.numberoflayers,
        md.mesh.z) and geometry information (md.geometry.thickness).
    vector : ndarray or array_like
        3D field values to be depth-averaged. Can be:
        - (md.mesh.numberofvertices,) for node-based data
        - (md.mesh.numberofelements,) for element-based data

    Returns
    -------
    vector_average : ndarray
        Depth-averaged field values projected onto the 2D mesh. Shape depends on input:
        - (md.mesh.numberofvertices2d,) for node-based input
        - (md.mesh.numberofelements2d,) for element-based input

    Raises
    ------
    TypeError
        If md is not provided or does not contain a 3D mesh.
    ValueError
        If vector size does not match expected dimensions.

    Notes
    -----
    The function uses the trapezoidal rule for integration:
    
    .. math::
        \\bar{v} = \\frac{1}{H} \\int_0^H v(z) dz
    
    where H is the total thickness and the integral is approximated using
    layer-by-layer trapezoidal integration.

    Examples
    --------
    Compute depth-averaged velocity:

    >>> vel_bar = depth_average(md, md.initialization.vel)

    Compute depth-averaged temperature:

    >>> temp_bar = depth_average(md, md.thermal.temperature)
    """

    # Error checks 
    if not md:
        raise TypeError("pyissm.model.mesh.depth_average: md must be provided")
    if md.mesh.domain_type().lower() != '3d':
        raise TypeError("pyissm.model.mesh.depth_average: md must contain a 3D mesh")

    # coerce to array in case float is passed
    if type(vector) not in [np.ndarray, np.ma.core.MaskedArray]:
        vector = np.array(vector)

    # Check if 2D array
    vec2d = False
    if vector.ndim == 2:
        vec2d = True
        vector = vector.reshape(-1, )

    # Handle node data
    if vector.shape[0] == md.mesh.numberofvertices:
        vector_average = np.zeros(md.mesh.numberofvertices2d)
        for i in range(1, md.mesh.numberoflayers):
            vector_average = vector_average + (project_2d(md, vector, i) + project_2d(md, vector, i + 1)) / 2. * (project_2d(md, md.mesh.z, i + 1) - project_2d(md, md.mesh.z, i))
        vector_average = vector_average / project_2d(md, md.geometry.thickness, 1)

    # Handle element data
    elif vector.shape[0] == md.mesh.numberofelements:
        vector_average = np.zeros(md.mesh.numberofelements2d)
        for i in range(1, md.mesh.numberoflayers):
            vertices_dz = (project_2d(md, md.mesh.z, i + 1) - project_2d(md, md.mesh.z, i))
            elements_dz = vertices_dz.mean(1)
            vector_average = vector_average + project_2d(md, vector, i) * elements_dz
    #vector_average = vector_average + project2d(md, vector, i) * (project2d(md, md.mesh.z, i + 1) - project2d(md, md.mesh.z, i))
        vertices_thickness = project_2d(md, md.geometry.thickness, 1)
        elements_thickness = vertices_thickness.mean(1)
        vector_average = vector_average / elements_thickness
    #vector_average = vector_average / project2d(md, md.geometry.thickness, 1)

    else:
        raise ValueError('pyissm.model.mesh.depth_average: vector size not supported')

    if vec2d:
        vector_average = vector_average.reshape(-1, )

    return vector_average
