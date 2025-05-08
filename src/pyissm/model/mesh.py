"""

Functions for building and interacting with an ISSM model mesh.

"""
import numpy as np
import matplotlib.tri as tri
from scipy.interpolate import griddata

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
