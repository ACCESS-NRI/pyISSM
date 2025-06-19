"""

Functions for building and interacting with an ISSM model mesh.

"""
import numpy as np
import matplotlib.tri as tri
from scipy.interpolate import griddata
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
    """

    ## Set default values
    is3d = False

    ## Process a 3D model
    if utils.has_nested_attr(md, 'mesh', 'elements2d'):

        # Create mesh object
        mesh = get_mesh(md.mesh.x2d, md.mesh.y2d, md.mesh.elements2d)

        # Extract X/Y Coordinates for 2D mesh
        mesh_x = md.mesh.x2d
        mesh_y = md.mesh.y2d

        # Extract and adjust 2D element numbering to be 0-indexed
        mesh_elements = md.mesh.elements2d - 1

        # Return is3d and display warning for 3D mesh
        is3d = True
        warnings.warn('3D model found. Processing as 2D mesh.')

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
    """

    ## Process model mesh
    mesh, mesh_x, mesh_y, mesh_elements, is3d = process_mesh(md)

    ## Identify ice/ocean nodes from surface layer of 3D model
    if is3d:
        ice_nodes = ice_levelset[md.mesh.vertexonsurface == 1] < 0
        ice_front_nodes = ice_levelset[md.mesh.vertexonsurface == 1] == 0
        ocean_nodes = ocean_levelset[md.mesh.vertexonsurface == 1] < 0

        warnings.warn('3D model found. Processing surface layer only.')

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
        raise ValueError('model_field must be defined on vertices or elements')

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
        raise ValueError('domain_mask should be defined on grid_x / grid_y.')
    elif domain_mask.dtype != bool:
        # If a custom domain_mask is supplied, it must be boolean
        raise TypeError('domain_mask should be boolean')

    # Apply mask
    gridded_model_field[:, ~domain_mask] = fill_value

    ## Squeeze output to remove time dimension if it's static
    ## -------------------------------------
    gridded_model_field = gridded_model_field.squeeze()

    return gridded_model_field
