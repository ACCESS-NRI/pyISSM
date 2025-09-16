"""
Functions for reading and writing ISSM models.

This module contains functions for reading and writing ISSM models to and from files.
"""

import netCDF4 as nc
import numpy as np
import pandas as pd
import os
import math

from .. import core, analysis, model, param, utils

def load_model(path: str) -> core.Model:

    ## Helper function to load different variables
    def get_variables(state, group, group_name):
        for var_name, var in group.variables.items():
            try:
                data = var[:]
                # If it's a scalar, extract the item, otherwise make numpy array
                if isinstance(data, np.ndarray) and data.shape == ():
                    state[var_name] = data.item()
                else:
                    state[var_name] = np.array(data)
            except Exception as e:
                print(f"⚠️ Failed to read variable '{group_name}.{var_name}': {e}")
                continue
        return state

    ## Helper function to load attributes
    def get_attributes(state, group):
        for attr in group.ncattrs():
            if attr != "classtype":
                state[attr] = group.getncattr(attr)
        return state

    ## Helper function to retrieve classtype and create new instance object
    def get_class(group):
        classtype = group.getncattr("classtype")
        obj = param.class_registry.create_instance(classtype)
        return classtype, obj

    ## Helper function to normalise NaN values (convert all NaN to np.nan)
    def normalize_nans(obj):
        if isinstance(obj, dict):
            return {k: normalize_nans(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [normalize_nans(item) for item in obj]
        elif isinstance(obj, tuple):
            return tuple(normalize_nans(item) for item in obj)
        elif isinstance(obj, np.ndarray):
            if np.issubdtype(obj.dtype, np.floating):
                obj = np.where(np.isnan(obj), np.nan, obj)
            return obj
        elif isinstance(obj, float) and math.isnan(obj):
            return np.nan
        else:
            return obj

    ## Initialise empty model class
    md = core.Model()

    ## Open the model netcdf file
    with nc.Dataset(path, 'r') as ds:

        ## Iterate through all top-level groups in the NetCDF file...
        for grp_name in ds.groups:
            grp = ds.groups[grp_name]

            ## --------------------------------------------------------
            ## Process the 'results' group
            ## - 'results' doesn't have a classtype attribute, but it contains
            ## various subgroups, like 'TransientSolution', 'StressbalanceSolution' etc.
            ## --------------------------------------------------------

            if grp_name == "results":

                ## Process indivdiual subgroups...
                for sub_grp_name, sub_grp in grp.groups.items():
                    print(f"ℹ️ Processing results group: {sub_grp_name}")

                    ## Check that a valid classtype exists:
                    if "classtype" in sub_grp.ncattrs():

                        ## Get the classtype for the subgroup & create new instance
                        classtype, obj = get_class(sub_grp)
                        # If obj is None, carry on (a warning is printed by create_instance in get_class)
                        if obj is None:
                            continue

                        ## Create empty state
                        state = {}

                        ## Get scalar attributes (those that are not stored as variables) and add to state
                        state = get_attributes(state, sub_grp)

                        ## Get variables
                        state = get_variables(state, sub_grp, sub_grp_name)

                        ## Convert all NaN values to np.nan
                        state = normalize_nans(state)

                        ## Set the state for the model class
                        try:
                            obj.__setstate__(state)
                        except Exception as e:
                            print(f"⚠️ Failed to set state for '{sub_grp_name}': {e}")
                            continue

                        ## Assign the object to the model (e.g., md.results.TransientSolution)
                        setattr(md.results, sub_grp_name, obj)

                    else:
                        print(f"⚠️️ classtype does not exist for group {grp_name}. Skipping...")
                        continue

            else:
                ## --------------------------------------------------------
                ## Process other model groups
                ## --------------------------------------------------------

                ## Check that a valid classtype exists:
                if "classtype" in grp.ncattrs():
                    ## Get the classtype for the group & create new instance
                    classtype, obj = get_class(grp)
                    # If obj is None, carry on (a warning is printed by create_instance in get_class)
                    if obj is None:
                        continue

                    ## Create empty state
                    state = {}

                    ## Get scalar attributes (those that are not stored as variables) and add to state
                    state = get_attributes(state, grp)

                    ## Get variables
                    state = get_variables(state, grp, grp_name)

                    ## Convert all NaN values to np.nan
                    state = normalize_nans(state)

                    ## Set the state for the model class
                    try:
                        obj.__setstate__(state)
                    except Exception as e:
                        print(f"⚠️ Failed to set state for '{grp_name}': {e}")
                        continue

                    ## Assign the object to the model (e.g., md.mesh)
                    setattr(md, grp_name, obj)

                else:
                    print(f"⚠️️ classtype does not exist for group {grp_name}. Skipping...")
                    continue

    return md


def save_model(path: str, md):

    ## Helper function to convert character array to string for NetCDF writing
    def char_array_to_strings(arr):
        arr = np.asarray(arr)  # Ensure it's a NumPy array
        if arr.ndim == 1:
            # Convert 1D string to single byte string array
            return np.array(["".join(arr.astype(str))], dtype='S')
        elif arr.ndim == 2:
            # Convert 2D strings to multiple byte string array
            return np.array(["".join(row.astype(str)) for row in arr], dtype='S')
        else:
            raise ValueError("Input must be a 1D or 2D char array with dtype='S1'")

    # Helper function to serialize an object's state
    def serialize_object(obj, group):
        """
        Serializes an object's state, including nested objects.
        """

        ## Get state from the object
        state = obj.__getstate__()

        ## For each item, write attributes and variables...
        for attr_name, value in state.items():

            # Handle scalars
            if isinstance(value, (int, float, str, bool)):
                group.setncattr(attr_name, value)

            # Handle arrays and lists (convert lists to arrays)
            elif isinstance(value, np.ndarray) or isinstance(value, list):

                # If it's a list, convert to an array
                if isinstance(value, list):
                    value = np.array(value, dtype='S')
                # Otherwise, check the array type
                else:
                    # If it's an object array, convert to a string array (object arrays can't be written to NetCDF)
                    if value.dtype == object:
                        value = np.array(value, dtype='S')
                    # Special handling for 'S1' datatype (these come from NetCDF Char variables when output from MATLAB)
                    elif value.dtype.kind == 'S':
                        value = char_array_to_strings(value)
                    else:
                        value = value

                # Handle the dimensions -- define a name from the size. If it doesn't already exist, create it.
                dim_names = []
                for i, size in enumerate(value.shape):
                    dim_name = f"dim_{i}_{size}"
                    if dim_name not in defined_dimensions:
                        ds.createDimension(dim_name, size)
                        defined_dimensions[dim_name] = size
                    dim_names.append(dim_name)

                # Create variable.
                # If it's a string array, turn off zlib compression
                if value.dtype.kind == 'S':
                    var = group.createVariable(attr_name, value.dtype, dimensions=dim_names, zlib=False)
                else:
                    var = group.createVariable(attr_name, value.dtype, dimensions=dim_names, zlib=True)

                ## Add data to variable
                var[:] = value

            ## Handle nested objects (recursively serialize them)
            elif isinstance(value, object):
                # If the value is a class instance, treat it as a group and recurse
                if hasattr(value, '__getstate__'):
                    nested_group = group.createGroup(attr_name)
                    serialize_object(value, nested_group)

            else:
                print(f"⚠️ Skipping unsupported field: {attr_name} ({type(value).__name__})")
                continue

    def get_registered_name(obj):
        classname = obj.__class__
        matching_keys = [k for k, v in param.class_registry.CLASS_REGISTRY.items() if v is classname]
        if not matching_keys:
            raise ValueError(f"Class {classname} is not registered.")
        registered_name = min(matching_keys, key=len)
        return registered_name

    ## Create the NetCDF file
    with nc.Dataset(path, 'w', format='NETCDF4') as ds:

        # Initialise dictionary to track existing dimension sizes & names
        defined_dimensions = {}

        # Loop through model attributes (top-level groups)
        for name, obj in vars(md).items():
            # Handle 'results' group specially
            if name == "results":
                results_group = ds.createGroup("results")

                # Loop through each solution type in md.results
                for solution_name, solution_obj in vars(md.results).items():
                    if solution_obj is None:
                        print(f"⚠️ Skipping solution type: {solution_name})")
                        continue

                    # Create subgroup for this solution (e.g., TransientSolution)
                    solution_group = results_group.createGroup(solution_name)

                    # Attach class type metadata
                    classname = get_registered_name(solution_obj)
                    solution_group.setncattr("classtype", classname)

                    # Serialize the solution state
                    serialize_object(solution_obj, solution_group)

            else:
                # For regular model components (e.g., mesh, materials, geometry)
                if obj is None:
                    continue

                # Create group for the model component
                group = ds.createGroup(name)

                # Attach class type metadata
                classname = get_registered_name(obj)
                group.setncattr("classtype", classname)

                # Serialize the component state
                serialize_object(obj, group)


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
    results to a NetCDF file. Variables are defined in a variable map which specifies the
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
    variable_map : str or pd.DataFrame, optional
        Custom variable mapping specification. If a string, it should be the path to a CSV file
        mapping model variables to output variable names and metadata. If a DataFrame, it should
        contain the same columns as the CSV file. If not provided, a default variable map will
        be used located in `../files/default_variable_map.csv` relative to this script.
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
    FileNotFoundError
        If the variable map CSV file does not exist.
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

    # If it's a string, read the CSV file
    if isinstance(variable_map, str):
        if not os.path.exists(variable_map):
            raise FileNotFoundError(f"export_gridded_model: Variable map file {variable_map} does not exist.")
        
        # Read the variable map CSV file
        var_map = pd.read_csv(variable_map)
    
    # If it's a DataFrame, use it directly
    elif isinstance(variable_map, pd.DataFrame):
        var_map = variable_map

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
        if 'TransientSolution' in var_map['issmModelSubgroup'].values and utils.general.has_nested_attr(md, 'results', 'TransientSolution'):
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

                ## If the variable is empty, skip it and print a warning
                if np.isnan(variable).all() or variable.shape[0] == 0:
                    if pd.isna(issm_subgroup):
                        print(f"\033[1m{issm_variable}\033[0m is empty in \033[1m{issm_group}\033[0m and will be skipped.")
                        continue
                    else:
                        print(f"\033[1m{issm_variable}\033[0m is empty in \033[1m{issm_group}.{issm_subgroup}\033[0m and will be skipped.")
                        continue   

                ## Check if unit conversion is required
                if pd.isnull(row['issmVariableUnit']) and pd.isnull(row['outputVariableUnit']):
                    pass
                elif row['issmVariableUnit'] == row['outputVariableUnit']:
                    pass
                else:
                    variable = utils.general.convert_units(row['issmVariableUnit'], row['outputVariableUnit'], variable)


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
