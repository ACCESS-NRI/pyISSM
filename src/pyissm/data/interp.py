"""
Tools to interpolate data to/from ISSM model mesh.

This module contains various interpolation functions that can be used in conjunction with ISSM models.
"""

import xarray as xr
import numpy as np
from scipy.interpolate import RegularGridInterpolator, griddata
from pyissm import tools

def xr_to_mesh(data,
               var_name,
               mesh_x,
               mesh_y,
               x_var = 'x',
               y_var = 'y',
               default_value = np.nan,
               interpolation_type = 'bilinear',
               issm_wrapper = True,
               crop_to_mesh = True,
               crop_buffer = 5000,
               fill_nan = False,
               fill_nan_interpolation_type = 'nearest'):
    
    """
    Interpolate a variable from an xarray dataset onto mesh nodes.
    
    Assumes rectilinear (structured) grid in the xarray dataset.
    
    Parameters
    ----------
    data : str or xr.Dataset
        Path to a netCDF file or an xarray Dataset containing the gridded data.
    var_name : str
        Name of the variable to interpolate.
    mesh_x : ndarray
        X-coordinates of mesh nodes.
    mesh_y : ndarray
        Y-coordinates of mesh nodes.
    x_var : str, optional
        Name of the x-coordinate variable in the dataset. Default is 'x'.
    y_var : str, optional
        Name of the y-coordinate variable in the dataset. Default is 'y'.
    default_value : float, optional
        Value to assign to points outside the grid domain. When ``fill_nan`` is True this
        acts as a last resort, applied only to nodes the fill pass could not reach.
        Default is np.nan.
    interpolation_type : str, optional
        Type of interpolation method. For ISSM wrapper: 'bilinear', 'nearest', etc.
        For scipy: 'linear', 'nearest', 'slinear', 'cubic', 'quintic', 'pchip'.
        Default is 'bilinear'.
    issm_wrapper : bool, optional
        If True, use ISSM wrapper functions for interpolation. If False, use scipy.
        Default is True.
    crop_to_mesh: bool, optional
        If True, the xarray Dataset is cropped to the extent of the mesh prior to
        interpolation. Only applied when the dataset has 1D coordinates. Default is True.
    crop_buffer: float or int, optional
        Buffer distance (m) applied to mesh bounding box to prevent edge effects.
        Default is 5000 m
    fill_nan: bool, optional
        If True, mesh nodes left unresolved by the initial interpolation are filled using the
        ``fill_nan_interpolation_type`` interpolation method. A node is unresolved if it falls
        outside the source domain, or if the source data is NaN there. ``default_value`` is then
        applied only to nodes the fill pass could not reach. Default is False.
    fill_nan_interpolation_type: str, optional
        Interpolation type used to fill unresolved nodes when ``fill_nan`` is True. Default
        is 'nearest'.
    
    Returns
    -------
    ndarray
        Interpolated variable values at mesh nodes.
    
    Raises
    ------
    TypeError
        If data is neither a file path nor an xarray Dataset.
    ValueError
        If variable is not 2D, coordinates are inconsistent, or grid is not rectilinear.
    ImportError
        If issm_wrapper is True but ISSM wrappers are not installed.
    """

    # Load xarray dataset if a filepath was given
    if isinstance(data, str):
        data = xr.open_dataset(data)
        close = True
    elif isinstance(data, xr.Dataset):
        close = False
    else:
        raise TypeError("pyissm.data.interp.xr_to_mesh: data must be a file path or an xarray Dataset")

    # Crop data to mesh extent
    # NOTE: label-based cropping requires 1D dimension coordinates. Datasets
    # carrying 2D (broadcast) coordinates are left uncropped and are reduced to
    # 1D further below.
    if crop_to_mesh and data[x_var].ndim == 1 and data[y_var].ndim == 1:

        ## Compute mesh bounding box, including buffer distance
        x_min = np.min(mesh_x) - crop_buffer
        x_max = np.max(mesh_x) + crop_buffer
        y_min = np.min(mesh_y) - crop_buffer
        y_max = np.max(mesh_y) + crop_buffer   

        ## Get coordinate arrays
        x = data[x_var]
        y = data[y_var]    

        ## Determine coordinate ordering
        x_ascending = (x[1] - x[0]) > 0
        y_ascending = (y[1] - y[0]) > 0

        ## Build slices based on coordinate ordering
        x_slice = slice(x_min, x_max) if x_ascending else slice(x_max, x_min)
        y_slice = slice(y_min, y_max) if y_ascending else slice(y_max, y_min)

        ## Crop data
        data = data.sel({x_var: x_slice,
                         y_var: y_slice})

    # Extract and squeeze arrays
    x = np.asarray(data[x_var].values).squeeze()
    y = np.asarray(data[y_var].values).squeeze()
    var_data = np.asarray(data[var_name].values).squeeze()

    if close:
        data.close()

    # Convert everything to float64 (but keep shapes)
    x = x.astype(np.float64, copy = False)
    y = y.astype(np.float64, copy = False)
    var_data = var_data.astype(np.float64, copy = False)
    mesh_x = mesh_x.astype(np.float64, copy = False)
    mesh_y = mesh_y.astype(np.float64, copy = False)

    # Check for rectilinear grid
    if var_data.ndim != 2:
        raise ValueError(f"pyissm.data.interp.xr_to_mesh: variable '{var_name}' must be 2D on a rectilinear grid")

    # If coordinates are 1D arrays, check shapes
    if x.ndim == 1 and y.ndim == 1:
        if var_data.shape != (y.size, x.size):
            raise ValueError(f"pyissm.data.interp.xr_to_mesh: variable '{var_name}' has shape {var_data.shape}, "
                             f"expected ({y.size}, {x.size})")

    # If coordinates are 2D arrays, check shapes and rectilinearity
    # If 2D, they should be repeated 1D arrays
    elif x.ndim == 2 and y.ndim == 2:
        if x.shape != var_data.shape or y.shape != var_data.shape:
            raise ValueError("pyissm.data.interp.xr_to_mesh: x, y, and variable must have identical shapes "
                             "for 2D coordinate grids")

        # Check rectilinearity
        if not np.allclose(x, x[0, :][None, :]):
            raise ValueError("pyissm.data.interp.xr_to_mesh: 2D x-coordinate is not rectilinear")

        if not np.allclose(y, y[:, 0][:, None]):
            raise ValueError("pyissm.data.interp.xr_to_mesh: 2D y-coordinate is not rectilinear")

        # Reduce to 1D
        x = x[0, :]
        y = y[:, 0]

    else:
        raise ValueError(
            "pyissm.data.interp.xr_to_mesh: x and y must both be 1D or both be 2D"
        )


    # Ensure monotonic increasing coordinates
    if np.any(np.diff(x) < 0):
        x = x[::-1]
        var_data = var_data[:, ::-1]

    if np.any(np.diff(y) < 0):
        y = y[::-1]
        var_data = var_data[::-1, :]

    # When a fallback fill is requested, the primary interpolation must leave
    # unresolved nodes as NaN so that they can be identified afterwards. Filling
    # them with default_value up front would make them indistinguishable from
    # genuinely interpolated values. default_value is applied at the end instead,
    # to whatever the fallback pass could not reach.
    primary_default_value = np.nan if fill_nan else default_value

    # If use_wrapper is True, use ISSM wrappers for interpolation
    if issm_wrapper:
        
        # Check that wrappers are installed
        if not tools.wrappers.check_wrappers_installed():
            raise ImportError("pyissm.data.interp.xr_to_mesh: ISSM wrappers are not installed. Please install them or set issm_wrapper = False to use scipy interpolation.")
        
        # Interpolate using ISSM wrapper
        var_on_mesh = tools.wrappers.InterpFromGridToMesh(
            x,
            y,
            var_data,
            mesh_x,
            mesh_y,
            primary_default_value,
            interpolation_type
        )
    
    # Otherwise, use scipy for interpolation        
    else:

        scipy_method_list = ['linear', 'nearest', 'slinear', 'cubic', 'quintic', 'pchip']
        if interpolation_type not in scipy_method_list:
            raise ValueError(f"pyissm.data.interp.xr_to_mesh: interpolation_type '{interpolation_type}' is not supported by scipy. Choose from {scipy_method_list}.")

        interp = RegularGridInterpolator(
                (y, x),
                var_data,
                method = interpolation_type,
                bounds_error = False,
                fill_value = primary_default_value,
            )

        # scipy expects points as (y, x)
        mesh_points = np.column_stack((mesh_y, mesh_x))

        # Extract variable on mesh
        var_on_mesh = interp(mesh_points)

    # If fill_nan is True, fill unresolved nodes from their nearest resolved
    # neighbour rather than leaving them at a constant. Nodes are unresolved
    # either because they fall outside the source domain, or because the source
    # data is NaN there (e.g. a product that is undefined over open ocean).
    if fill_nan:

        # Identify points that are non-nan
        valid_mask = np.isfinite(var_on_mesh)

        # Nothing to do if every node was resolved by the primary interpolation
        if not np.all(valid_mask):

            # Fill nan values using fill_nan_interpolation_type interpolation
            filled_points = points_to_mesh(
                data_x = mesh_x[valid_mask],
                data_y = mesh_y[valid_mask],
                data_values = var_on_mesh[valid_mask],
                mesh_x = mesh_x,
                mesh_y = mesh_y,
                interpolation_type = fill_nan_interpolation_type
            )

            # Replace nan values
            var_on_mesh[~valid_mask] = filled_points[~valid_mask]

            # Apply default_value to any node the fallback pass could not reach
            unfilled_mask = ~np.isfinite(var_on_mesh)
            if np.any(unfilled_mask):
                var_on_mesh[unfilled_mask] = default_value

    return var_on_mesh

def points_to_mesh(data_x,
                   data_y,
                   data_values,
                   mesh_x,
                   mesh_y,
                   default_value = np.nan,
                   interpolation_type = 'linear',
                   fill_nan = False,
                   fill_nan_interpolation_type = 'nearest'):
    """
    Interpolate scattered points onto mesh node coordinates using scipy.

    Parameters
    ----------
    data_x : array_like
        X-coordinates of the scattered data points. Can be 1D or 2D; if 2D it will be flattened.
    data_y : array_like
        Y-coordinates of the scattered data points. Must have the same shape as `data_x`.
    data_values : array_like
        Values at the scattered data points. Must have the same shape as `data_x`.
    mesh_x : array_like
        X-coordinates of mesh nodes. 1D array.
    mesh_y : array_like
        Y-coordinates of mesh nodes. 1D array.
    default_value : float, optional
        Value used to fill points outside the convex hull of the input data. Default is `np.nan`.
    interpolation_type : str, optional
        Interpolation method passed to `scipy.interpolate.griddata`. Supported options: 'linear', 'nearest', 'cubic'.
        Default is 'linear'.
    fill_nan : bool, optional
        If True, mesh nodes left unresolved by the initial interpolation (i.e. those outside the convex hull of the
        input data) are filled from the input data using the ``fill_nan_interpolation_type`` method, instead of being
        assigned ``default_value``. ``default_value`` is then applied only to nodes that remain unresolved.
        Default is False.
    fill_nan_interpolation_type : str, optional
        Interpolation type used to fill unresolved nodes when ``fill_nan`` is True. Default is 'nearest'.

    Returns
    -------
    ndarray
        1D array of interpolated values at the mesh nodes (shape equals `mesh_x`/`mesh_y`).
        Points outside the convex hull of the input data are assigned `default_value`, unless `fill_nan` is True.

    Raises
    ------
    ValueError
        If `interpolation_type` is not supported by scipy, if input shapes are inconsistent, or if no valid
        data points remain after removing NaNs/Infs.
    """

    # Check interpolation type    
    scipy_method_list = ['linear', 'nearest', 'cubic']
    if interpolation_type not in scipy_method_list:
        raise ValueError(f"pyissm.data.interp.points_to_mesh: interpolation_type '{interpolation_type}' is not supported by scipy. Choose from {scipy_method_list}.")
    
    if fill_nan and fill_nan_interpolation_type not in scipy_method_list:
        raise ValueError(f"pyissm.data.interp.points_to_mesh: fill_nan_interpolation_type '{fill_nan_interpolation_type}' is not supported by scipy. Choose from {scipy_method_list}.")

    # Validate input shapes
    if data_x.shape != data_y.shape or data_x.shape != data_values.shape:
        raise ValueError("pyissm.data.interp.points_to_mesh: data_x, data_y, and data_values must have the same shape.")
    
    # Flatten if needed:
    if data_values.ndim == 2:
        data_x = data_x.flatten()
        data_y = data_y.flatten()
        data_values = data_values.flatten()
    elif data_values.ndim == 1:
        pass
    else:
        raise ValueError("pyissm.data.interp.points_to_mesh: data_values must be 1D or 2D array.")
    
    # Remove NaNs/Inf from input data
    mask = np.isfinite(data_x) & np.isfinite(data_y) & np.isfinite(data_values)

    if not np.any(mask):
        raise ValueError("pyissm.data.interp.points_to_mesh: no valid data points available for interpolation")

    data_x = data_x[mask]
    data_y = data_y[mask]
    data_values = data_values[mask]

    # Prepare points for interpolation
    data_points = np.column_stack((data_y, data_x))
    interp_points = np.column_stack((mesh_y, mesh_x))
    
    # As above, unresolved nodes must stay NaN when a fallback fill is requested
    # so they remain distinguishable from interpolated values.
    primary_default_value = np.nan if fill_nan else default_value

    # Perform interpolation    
    data_on_mesh = griddata(data_points,
                            data_values,
                            interp_points,
                            method = interpolation_type,
                            fill_value = primary_default_value)

    # Fill nodes left unresolved by the initial interpolation (those outside the
    # convex hull of the input data) from the input data itself, rather than
    # assigning them a constant.
    if fill_nan:

        valid_mask = np.isfinite(data_on_mesh)

        # Nothing to do if every node was resolved by the initial interpolation
        if not np.all(valid_mask):

            # Only the unresolved nodes need re-interpolating
            data_on_mesh[~valid_mask] = griddata(data_points,
                                                 data_values,
                                                 interp_points[~valid_mask],
                                                 method = fill_nan_interpolation_type)

            # Apply default_value to any node the fallback pass could not reach
            unfilled_mask = ~np.isfinite(data_on_mesh)
            if np.any(unfilled_mask):
                data_on_mesh[unfilled_mask] = default_value

    return data_on_mesh

def mesh_to_xr(md,
               model_field,
               spacing = None,
               x = None,
               y = None,
               interpolation_type = 'linear',
               domain_mask = None,
               fill_value = np.nan,
               time = None,
               crs = None,
               name = None,
               attrs = None):
    """
    Interpolate a model field from the mesh onto a regular grid, as an xarray DataArray.

    This is the counterpart to :func:`xr_to_mesh`. The interpolation itself is performed by
    :func:`pyissm.model.mesh.grid_model_field`, which masks the result to the model domain rather
    than to the convex hull of the mesh nodes. This function adds the grid construction and the
    coordinate, time and CRS metadata needed to use the result as data (e.g. to write to netCDF,
    or to compare against gridded observations).

    Parameters
    ----------
    md : :class:`pyissm.model.Model`
        Model object providing the mesh that ``model_field`` is defined on.
    model_field : ndarray
        Field to interpolate. Either ``(npoints,)`` for a static field or ``(nt, npoints)`` for a
        time-varying one, where ``npoints`` matches the number of mesh vertices or elements.
    spacing : float, optional
        Grid spacing, in mesh coordinate units. The grid extent is taken from the mesh bounding box.
        Ignored if both ``x`` and ``y`` are given; one of ``spacing`` or ``x``/``y`` is required.
    x : ndarray, optional
        1D array of x-coordinates defining the output grid. Must be given together with ``y``.
    y : ndarray, optional
        1D array of y-coordinates defining the output grid. Must be given together with ``x``.
    interpolation_type : str, optional
        Interpolation method passed to `scipy.interpolate.griddata`. Supported options: 'linear',
        'nearest', 'cubic'. Default is 'linear'.
    domain_mask : ndarray of bool, optional
        Mask of valid grid cells, shaped ``(len(y), len(x))``. If not given, a mask is generated
        from the model mesh. Cells where the mask is False are set to ``fill_value``.
    fill_value : float, optional
        Value assigned to grid cells outside the model domain. Default is np.nan.
    time : array_like, optional
        Values for the time coordinate of a time-varying field. Defaults to a simple integer index.
        Ignored for static fields.
    crs : int or str or pyproj.CRS, optional
        Coordinate reference system of the mesh (e.g. 3031 for Antarctic Polar Stereographic, or
        'EPSG:3413' for Greenland). If given, it is attached as a CF-style ``spatial_ref`` scalar
        coordinate. The mesh itself carries no projection information, so no CRS is assumed when
        this is not given.
    name : str, optional
        Name of the returned DataArray.
    attrs : dict, optional
        Additional attributes (e.g. ``units``, ``long_name``) to attach to the returned DataArray.

    Returns
    -------
    xr.DataArray
        Gridded field with dimensions ``(y, x)`` for a static field, or ``(time, y, x)`` for a
        time-varying one.

    Raises
    ------
    ValueError
        If neither ``spacing`` nor both ``x`` and ``y`` are given, if ``x``/``y`` are not 1D, or if
        ``model_field`` is not defined on mesh vertices or elements.

    Examples
    --------
    .. code-block:: python

        >>> # Static field on a 5 km grid, tagged as Antarctic Polar Stereographic
        >>> thickness = pyissm.data.interp.mesh_to_xr(md, md.geometry.thickness,
        ...                                           spacing = 5000, crs = 3031,
        ...                                           name = 'thickness',
        ...                                           attrs = {'units': 'm'})
        >>> thickness.to_netcdf('thickness.nc')

        >>> # Time-varying field, onto a grid the caller defines
        >>> vel = pyissm.data.interp.mesh_to_xr(md, vel_over_time,
        ...                                     x = np.arange(-3e6, 3e6, 10000),
        ...                                     y = np.arange(-3e6, 3e6, 10000),
        ...                                     time = solution_times)

    See Also
    --------
    xr_to_mesh : Interpolate gridded data onto mesh nodes (the reverse operation).
    pyissm.model.mesh.grid_model_field : Underlying mesh-to-grid interpolation.
    """

    # Imported here rather than at module level: pyissm.data is initialised before pyissm.model,
    # so a module-level import of pyissm.model would be circular.
    from pyissm.model import mesh as model_mesh

    model_field = np.asarray(model_field)

    ## Build the output grid
    ## -------------------------------------
    if x is not None and y is not None:
        x = np.asarray(x, dtype = np.float64)
        y = np.asarray(y, dtype = np.float64)

        if x.ndim != 1 or y.ndim != 1:
            raise ValueError('pyissm.data.interp.mesh_to_xr: x and y must be 1D arrays')

    elif spacing is not None:
        # Derive the grid extent from the mesh bounding box. np.arange can stop just short of the
        # maximum, so the stop is nudged by half a cell to keep the far edge of the mesh covered.
        x = np.arange(np.min(md.mesh.x), np.max(md.mesh.x) + spacing / 2, spacing)
        y = np.arange(np.min(md.mesh.y), np.max(md.mesh.y) + spacing / 2, spacing)

    else:
        raise ValueError('pyissm.data.interp.mesh_to_xr: either spacing, or both x and y, must be provided')

    grid_x, grid_y = np.meshgrid(x, y)

    ## Interpolate mesh -> grid
    ## -------------------------------------
    gridded_field = model_mesh.grid_model_field(md,
                                                model_field,
                                                grid_x,
                                                grid_y,
                                                method = interpolation_type,
                                                domain_mask = domain_mask,
                                                fill_value = fill_value)

    ## Assemble the DataArray
    ## -------------------------------------
    # grid_model_field squeezes its output, which would also drop a length-1 spatial axis, so the
    # expected shape is restored explicitly here.
    is_static = (model_field.ndim == 1)

    if is_static:
        gridded_field = np.reshape(gridded_field, (y.size, x.size))
        dims = ('y', 'x')
        coords = {'y': y, 'x': x}
    else:
        nt = model_field.shape[0]
        gridded_field = np.reshape(gridded_field, (nt, y.size, x.size))
        dims = ('time', 'y', 'x')
        coords = {'time': np.arange(nt) if time is None else np.asarray(time),
                  'y': y,
                  'x': x}

    data_array = xr.DataArray(gridded_field,
                              dims = dims,
                              coords = coords,
                              name = name,
                              attrs = dict(attrs) if attrs else {})

    ## Attach the CRS, if one was given
    ## -------------------------------------
    # Follows the CF grid mapping convention: a scalar coordinate holding the CRS definition, which
    # the data variable points at via a 'grid_mapping' attribute. This is what rioxarray, GDAL and
    # QGIS look for.
    if crs is not None:
        try:
            from pyproj import CRS

            crs_attrs = CRS.from_user_input(crs).to_cf()

        except ImportError:
            # pyproj is not a hard dependency of pyissm, so fall back to recording the CRS verbatim
            crs_attrs = {'crs_wkt': str(crs)}

        data_array = data_array.assign_coords(spatial_ref = 0)
        data_array.coords['spatial_ref'].attrs = crs_attrs
        data_array.attrs['grid_mapping'] = 'spatial_ref'

    return data_array

