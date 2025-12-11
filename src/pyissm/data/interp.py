"""
Tools to interpolate data to/from ISSM model mesh.

This module contains various interpolation functions that can be used in conjunction with ISSM models.
"""

import xarray as xr
import numpy as np
from pyissm import tools

def xr_to_mesh(data,
               var_name,
               mesh_x,
               mesh_y,
               x_var = 'x',
               y_var = 'y',
               default_value = np.nan,
               interpolation_type = 'bilinear'):
    
    """Interpolate a variable from an xarray dataset onto the mesh nodes."""

    # Load xarray dataset if a filepath was given
    if isinstance(data, str):
        data = xr.open_dataset(data)
        close = True
    elif isinstance(data, xr.Dataset):
        close = False
    else:
        raise TypeError("pyissm.data.interp.xr_to_mesh: data must be a file path or an xarray Dataset")

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

    # Ensure Y increases (ISSM expects increasing y-axis)
    if y.ndim == 1:
        if np.any(np.diff(y) < 0):
            y = y[::-1]
            var_data = var_data[::-1, :]
    else:
        # 2-D grid
        if np.diff(y[:, 0]).mean() < 0:
            y = y[::-1, :]
            var_data = var_data[::-1, :]

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

    return var_on_mesh
