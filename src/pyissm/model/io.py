"""
Functions for reading and writing ISSM models.

This module contains functions for reading and writing ISSM models to and from files.
"""

from types import SimpleNamespace
from ..core import Model
from .. import analysis
from .. import model
from .. import utils
import netCDF4 as nc
import numpy as np
import pandas as pd
import os

def load_model_group(nc_file,
                     group_name,
                     sub_group_name = None):
    """
    Load variables from a specified group (and optional subgroup) in a NetCDF file.

    Extracts all variables from the given group (or subgroup if specified) in an open
    NetCDF file and returns them as attributes of a `SimpleNamespace` for convenient access.

    Parameters
    ----------
    nc_file : netCDF4.Dataset
        An open NetCDF file object.
    group_name : str
        Name of the group from which to load variables.
    sub_group_name : str, optional
        Name of the subgroup within the specified group. If not provided, variables are
        loaded from the main group.

    Returns
    -------
    types.SimpleNamespace
        A namespace object containing the variables from the specified group or subgroup.
        Each attribute corresponds to a variable, stored as a NumPy array.
    """

    ## Initialise empty variable_dict
    variable_dict = {}

    ## Extract the group (or sub_group, if defined)
    group = nc_file.groups[group_name]
    if sub_group_name:
        group = group.groups[sub_group_name]
    else:
        pass

    ## List all variables within the specified group
    variable_names = group.variables

    ## Loop over variables and add to variable_dict
    for variable in variable_names:
        variable_dict[variable] = group.variables[variable][:]

    ## Convert variable_dict to SimpleNamespace to return
    variables = SimpleNamespace(**variable_dict)

    return variables

def load_model(filepath):

    """
    Load an ISSM model from a NetCDF file.

    This function reads a NetCDF file containing an ISSM model, reconstructs the hierarchical group and subgroup
    structure, and populates an instance of the `Model` class with the data.

    Parameters
    ----------
    filepath : str
        Path to the NetCDF file containing the model data.

    Returns
    -------
    Model
        An instance of the `Model` class populated with data from the NetCDF file.

    Notes
    -----
    - Each top-level group in the NetCDF file is added as an attribute of the `Model` instance.
    - If a group contains nested subgroups, those are recursively added as attributes of their parent group.
    - Uses `load_model_group()` to extract group and subgroup data.
    - The function expects the NetCDF file to follow the standard ISSM model structure.
    """

    ## Initialise empty model to add information to
    md = Model()

    ## Open netCDF file (read-only)
    nc_file = nc.Dataset(filepath, mode = 'r')

    ## List available groups
    model_groups = list(nc_file.groups.keys())

    ## Loop through available groups and add them to md
    for group_name in model_groups:

        # Load data for given group
        group_data = load_model_group(nc_file, group_name)

        # Add group_data to md
        setattr(md, group_name, group_data)

        ## Check for subgroups
        if len(nc_file.groups[group_name].groups) > 0:

            # List all subgroup_names
            subgroup_names = list(nc_file.groups[group_name].groups.keys())

            for subgroup in subgroup_names:

                # Get the subgroup data
                subgroup_data = load_model_group(nc_file, group_name, subgroup)

                # Get the group so the subgroup can be added
                group = getattr(md, group_name)

                # Add the subgroup_data
                setattr(group, subgroup, subgroup_data)
        else:
            pass

    return md

def export_gridded_model(md,
                         out_file,
                         grid_x,
                         grid_y,
                         variable_map = None,
                         method = 'linear',
                         domain_mask = None,
                         fill_value = np.nan):
    """
    Export gridded model variables to a NetCDF file based on a variable mapping specification.

    This function interpolates ISSM model output variables onto a regular 2D grid and writes the
    results to a NetCDF file. Variables are defined in a CSV variable map which specifies the
    desired output name, source location in the model, and optional unit conversions.

    Parameters
    ----------
    md : object
        The ISSM model object containing simulation results and mesh information.
    out_file : str
        Path to the output NetCDF file.
    grid_x : ndarray
        2D array of X coordinates for the regular grid.
    grid_y : ndarray
        2D array of Y coordinates for the regular grid.
    variable_map : str, optional
        Path to a CSV file mapping model variables to output variable names and metadata.
        If not provided, a default variable map located in `../files/default_variable_map.csv`
        relative to this script will be used.
    method : str, optional
        Interpolation method used to grid model data. Specific variables override this option.
        Options are `'linear'`, `'nearest'`, etc. Default is `'linear'`.
    domain_mask : ndarray, optional
        A boolean or numerical mask to apply to the output grid. If provided, values outside
        the domain are masked or set to `fill_value`.
    fill_value : float, optional
        Value to use for missing or masked data in the NetCDF output. Default is `np.nan`.

    Raises
    ------
    ValueError
        If the variable map file contains duplicate output variable names.
    Exception
        If any unexpected error occurs during export. The output file is removed in this case.

    Notes
    -----
    - Supports both static and time-dependent model fields.
    - Handles custom unit conversions as defined in the variable map.
    - Includes ISMIP-specific variable logic where applicable.
    - Output NetCDF will contain `grid_x`, `grid_y`, and optionally `time` dimensions.

    """

    ## Define list of variable to force nearest-neighbour interpolation
    nn_interp_list = ['ice_levelset', 'ocean_levelset', 'MaskOceanLevelset', 'MaskIceLevelset']

    ## Get variable map
    # If not defined, use default map in ../files
    if variable_map is None:
        variable_map = os.path.join(os.path.dirname(__file__), '../files/default_variable_map.csv')
        variable_map = os.path.abspath(variable_map)

    var_map = pd.read_csv(variable_map)

    ## Error Checks
    # Check that all outputVariableNames are unique
    if any(var_map['outputVariableName'].duplicated()):
        raise ValueError(f"export_gridded_model: Duplicate outputVariableName found in {variable_map}.")

    ## Wrap in try so that file can be removed if any error occurs allowing easy re-try
    try:
        ## Create NetCDF file & dimensions
        ny, nx = grid_x.shape
        nc_file = nc.Dataset(out_file, 'w', format = 'NETCDF4')
        nc_file.createDimension('grid_x', nx)
        nc_file.createDimension('grid_y', ny)
        # Add variables
        var_x = nc_file.createVariable('grid_x', 'f4', ('grid_x'), fill_value=fill_value)
        var_y = nc_file.createVariable('grid_y', 'f4', ('grid_y'), fill_value=fill_value)
        var_x[:] = grid_x[0,:]
        var_y[:] = grid_y[:,0]

        # If TransientSolution is in the request and exists in the model, create the time dimension
        if 'TransientSolution' in var_map['issmModelSubgroup'].values and utils.has_nested_attr(md, 'results', 'TransientSolution'):
            time = getattr(md.results.TransientSolution, 'time')
            nc_file.createDimension('time', len(time))
            var_time = nc_file.createVariable('time', 'f4', ('time'), fill_value=fill_value)
            var_time[:] = time

        ## Loop over each row of the var_map dataframe.
        for _, row in var_map.iterrows():

            ## Extract relevant group / sub-group / variable information
            issm_group = row['issmModelGroup']
            issm_subgroup = row['issmModelSubgroup']
            issm_variable = row['issmVariableName']

            ## Check for ISMIP6 specific requirements
            if issm_group == 'ISMIP6':

                # Try to get variable
                ismip_variable = row['outputVariableName']
                variable = analysis.ismip.get_ismip_variable(md, ismip_variable)

                # If ISMIP specific requirements were not met, continue
                if variable is None:
                    continue
                else:
                    print(f"Computing ISMIP6 Variable: \033[1m{ismip_variable}\033[0m")

            ## Otherwise, try to get the variable directly from the model
            else:
                ## Check that issm_group exists in the model. If not, continue
                if not hasattr(md, issm_group):
                    print(f"The following group is missing and will be skipped: \033[1m{issm_group}\033[0m")
                    continue

                ## Extract the issm_group
                group = getattr(md, issm_group)

                ## Check that issm_subgroup exists. If not, continue
                if pd.isnull(issm_subgroup):

                    ## If no subgroup is defined, it's not a nested result. Use the parent group
                    sub_group = group
                else:
                    if not hasattr(group, issm_subgroup):
                        print(f"\033[1m{issm_group}.{issm_subgroup}\033[0m is missing and will be skipped.")
                        continue

                    ## Extract the issm_subgroup
                    sub_group = getattr(group, issm_subgroup)

                ## Check that the variable exists
                if not hasattr(sub_group, issm_variable):
                    print(f"\033[1m{issm_variable}\033[0m is missing in \033[1m{issm_group}.{issm_subgroup}\033[0m and will be skipped. ")
                    continue

                ## Extract the variable
                variable = getattr(sub_group, issm_variable)

                ## Check if unit conversion is required
                if pd.isnull(row['issmVariableUnit']) and pd.isnull(row['outputVariableUnit']):
                    pass
                elif row['issmVariableUnit'] == row['outputVariableUnit']:
                    pass
                else:
                    variable = utils.convert_units(row['issmVariableUnit'], row['outputVariableUnit'], variable)


            ## ------------------------------------------------
            ## At this point, the variable exists and it should be added to the NetCDF file

            ## CASE 1 - Transient 2D field
            if variable.ndim == 2:

                ## Grid the variable
                # If mask/levelset variable requested, force nearest-neighbour interpolation method
                if issm_variable in nn_interp_list:
                    print(f"Gridding: \033[1m{issm_variable}\033[0m using NN interpolation")
                    variable_grid = model.mesh.grid_model_field(md, variable, grid_x, grid_y, method = 'nearest', domain_mask = domain_mask)

                # Otherwise interpolate using the specified method
                else:
                    variable_grid = model.mesh.grid_model_field(md, variable, grid_x, grid_y, method = method, domain_mask = domain_mask)
                    print(f"Gridding: \033[1m{issm_variable}\033[0m")

                ## Create variable in nc_file with t/y/x dimensions
                nc_var = nc_file.createVariable(row['outputVariableName'], 'f4', ('time', 'grid_y', 'grid_x'), fill_value = fill_value)
                nc_var[:] = variable_grid
                nc_var.long_name = row['outputVariableLongName']
                nc_var.units = row['outputVariableUnit']


            ## Static 2D fields & timeseries data
            if variable.ndim == 1:

                ## CASE 2 - Static 2D field (defined on vertices OR elements)
                if (len(variable) == md.mesh.numberofvertices) or (len(variable) == md.mesh.numberofelements):

                    ## Grid the variable
                    # If mask/levelset variable requested, force nearest-neighbour interpolation method
                    if issm_variable in nn_interp_list:
                        print(f"Gridding: \033[1m{issm_variable}\033[0m using NN interpolation")
                        variable_grid = model.mesh.grid_model_field(md, variable, grid_x, grid_y, method='nearest', domain_mask=domain_mask)

                    # Otherwise interpolate using the specified method
                    else:
                        variable_grid = model.mesh.grid_model_field(md, variable, grid_x, grid_y, method=method, domain_mask=domain_mask)
                        print(f"Gridding: \033[1m{issm_variable}\033[0m")

                    ## Create variable in nc_file with y/x dimensions
                    nc_var = nc_file.createVariable(row['outputVariableName'], 'f4', ('grid_y', 'grid_x'), fill_value = fill_value)
                    nc_var[:] = variable_grid
                    nc_var.long_name = row['outputVariableLongName']
                    nc_var.units = row['outputVariableUnit']

                ## CASE 3 - Timeseries
                # NOTE: time is not defined if it's not requested or doesn't exist in the model
                try:
                    if len(variable) == len(time):

                        if not issm_group == 'ISMIP6':
                            print(f"Found: \033[1m{issm_variable}\033[0m")

                        ## Create variable in nc_file with t dimensions
                        nc_var = nc_file.createVariable(row['outputVariableName'], 'f4', ('time'), fill_value=fill_value)
                        nc_var[:] = variable
                        nc_var.long_name = row['outputVariableLongName']
                        nc_var.units = row['outputVariableUnit']

                except NameError:
                    pass


        nc_file.close()
        return

    ## If there is any error, remove the out_file and raise the error
    except Exception as e:
        if os.path.exists(out_file):
            os.remove(out_file)
        print(f"export_grided_model: Export failed -- {e}\n\033[91mModel not written to file.\033[0m")
        raise
