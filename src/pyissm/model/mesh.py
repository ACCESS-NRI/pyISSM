"""

Functions for building and interacting with an ISSM model mesh.

"""
import numpy as np
import matplotlib.tri as tri
from scipy.interpolate import griddata
import scipy.sparse
import warnings
from .. import utils

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
    if utils.general.has_nested_attr(md, 'mesh', 'elements2d'):

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
            if utils.wrappers.check_wrappers_installed():
                flag = utils.wrappers.ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y, region, 'element', 1).astype(bool)
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
