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
               crop_buffer = 5000):
    
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
        Value to assign to points outside the grid domain. Default is np.nan.
    interpolation_type : str, optional
        Type of interpolation method. For ISSM wrapper: 'bilinear', 'nearest', etc.
        For scipy: 'linear', 'nearest', 'slinear', 'cubic', 'quintic', 'pchip'.
        Default is 'bilinear'.
    issm_wrapper : bool, optional
        If True, use ISSM wrapper functions for interpolation. If False, use scipy.
        Default is True.
    crop_to_mesh: bool, optional
        If True, the xarraxy Dataset is cropped to the extent of the mesh prior to
        interpolation. Default is True.
    crop_buffer: float or int, optional
        Buffer distance (m) applied to mesh bounding box to prevent edge effects.
        Default is 5000 m
    
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
    if crop_to_mesh:

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
            default_value,
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
                fill_value = default_value,
            )

        # scipy expects points as (y, x)
        mesh_points = np.column_stack((mesh_y, mesh_x))

        # Extract variable on mesh
        var_on_mesh = interp(mesh_points)

    return var_on_mesh

def points_to_mesh(data_x,
                   data_y,
                   data_values,
                   mesh_x,
                   mesh_y,
                   default_value = np.nan,
                   interpolation_type = 'linear'):
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

    Returns
    -------
    ndarray
        1D array of interpolated values at the mesh nodes (shape equals `mesh_x`/`mesh_y`).
        Points outside the convex hull of the input data are assigned `default_value`.

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
    
    # Perform interpolation    
    data_on_mesh = griddata(data_points,
                            data_values,
                            interp_points,
                            method = interpolation_type,
                            fill_value = default_value)

    return data_on_mesh
    
