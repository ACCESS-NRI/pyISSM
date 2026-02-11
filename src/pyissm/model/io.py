"""
Functions for reading and writing ISSM models.

This module contains functions for reading and writing ISSM models to and from files.
"""

import netCDF4 as nc
import numpy as np
import pandas as pd
import os
import math
import subprocess
import warnings
import sys
import shutil
import collections

from pyissm import analysis, model, tools

def load_model(path):
    """
    Load an ISSM model from a NetCDF file.
    
    This function reads an ISSM model that has been saved in NetCDF format,
    reconstructing all model components including mesh, materials, geometry,
    boundary conditions, and results. The function handles nested objects
    and properly deserializes the complete model state.
    
    Parameters
    ----------
    path : str
        Path to the NetCDF file containing the ISSM model data.
    
    Returns
    -------
    model.Model
        The reconstructed ISSM model object with all components loaded
        from the NetCDF file.
    
    Raises
    ------
    FileNotFoundError
        If the specified NetCDF file does not exist.
    ValueError
        If the NetCDF file is corrupted or missing required metadata.
    
    Notes
    -----
    - The function automatically handles different data types and converts
      NaN values to numpy.nan for consistency
    - TransientSolution results are automatically expanded from collapsed
      format back to individual timesteps
    - Missing or corrupted groups are skipped with warning messages
    
    See Also
    --------
    save_model : Save an ISSM model to NetCDF format
    
    Examples
    --------
    >>> md = load_model('my_model.nc')
    >>> print(f"Model loaded with {md.mesh.numberofvertices} vertices")
    """

    # Helper function to load different variables
    def _get_variables(state, group, group_name):
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

    # Helper function to load attributes
    def _get_attributes(state, group):
        for attr in group.ncattrs():
            if attr != "classtype":
                val = group.getncattr(attr)
                if isinstance(val, str) and val == "__EMPTY_LIST__":
                    state[attr] = []
                    continue
                else:
                    state[attr] = val
        return state

    # Helper function to retrieve classtype and create new instance object
    def _get_class(group):
        classtype = group.getncattr("classtype")
        obj = model.classes.class_registry.create_instance(classtype)
        return classtype, obj

    # Helper function to normalise NaN values (convert all NaN to np.nan)
    def _normalize_nans(obj):
        if isinstance(obj, dict):
            return {k: _normalize_nans(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [_normalize_nans(item) for item in obj]
        elif isinstance(obj, tuple):
            return tuple(_normalize_nans(item) for item in obj)
        elif isinstance(obj, np.ndarray):
            if np.issubdtype(obj.dtype, np.floating):
                obj = np.where(np.isnan(obj), np.nan, obj)
            return obj
        elif isinstance(obj, float) and math.isnan(obj):
            return np.nan
        else:
            return obj
        
    def _normalize_loaded_attributes(obj):
        """
        Normalize a loaded model object by converting data types for consistency.
        
        This function performs post-loading normalization on model objects to ensure
        data types are consistent with Python conventions. It handles various data
        type conversions that may be necessary after loading from NetCDF format.
        
        Parameters
        ----------
        obj : object
            The model object to normalize. Must have a __dict__ attribute containing
            the attributes to be normalized.
        
        Notes
        -----
        The function performs the following normalizations:
        
        - NumPy integers are converted to Python int for consistency
        - 1D object arrays containing strings are converted to Python lists
        - Legacy character arrays are converted to lists of strings:
          
          - 1D char arrays become single-element lists: ['string']
          - 2D char arrays (MATLAB-style) become lists of strings
          
        - Nested objects with __dict__ attributes are recursively processed
        
        The function modifies the object in-place by updating its attributes
        using setattr().
        """
        
        ## Helper function to unpack MATLAB-style 2D char arrays
        def _unpack_char_cell(char_arr):
            """
            Convert MATLAB-style 2D char array (num_strings x max_len) to list of Python strings
            """
            num_strings, max_len = char_arr.shape
            str_list = []
            for i in range(num_strings):
                row = char_arr[i, :]
                # decode bytes if needed
                if row.dtype.kind == 'S':
                    s = b"".join(row.astype('S1')).decode('utf-8').rstrip()
                else:
                    s = "".join(row).rstrip()
                if s:
                    str_list.append(s)
            return str_list

        ## Iterate through all attributes of the object
        for attr_name, value in vars(obj).items():

            ### Convert NumPy integers to Python integers
            if isinstance(value, np.integer):
                setattr(obj, attr_name, int(value))

            ### Convert 1D object arrays to Python lists
            elif isinstance(value, np.ndarray) and value.dtype == object and value.ndim == 1:
                setattr(obj, attr_name, value.tolist())

            ### Convert legacy Char arrays to lists of strings
            elif isinstance(value, np.ndarray) and value.dtype.kind in ['S','U']:

                if value.ndim == 2:
                    #### MATLAB-style 2D char array: convert to list of strings
                    str_list = _unpack_char_cell(value)
                    setattr(obj, attr_name, str_list)

                elif value.ndim == 1:
                    #### 1D char array: single string, put in a list
                    s = b"".join(value.astype('S1')).decode('utf-8') if value.dtype.kind == 'S' else "".join(value)
                    setattr(obj, attr_name, [s])

            ### Recursively normalize nested objects
            elif hasattr(value, "__dict__"):
                _normalize_loaded_attributes(value)

    # Initialise empty model class
    md = model.Model()

    # Open the model netcdf file
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
                        classtype, obj = _get_class(sub_grp)
                        # If obj is None, carry on (a warning is printed by create_instance in get_class)
                        if obj is None:
                            continue

                        ## Create empty state
                        state = {}

                        ## Get scalar attributes (those that are not stored as variables) and add to state
                        state = _get_attributes(state, sub_grp)

                        ## Get variables
                        state = _get_variables(state, sub_grp, sub_grp_name)

                        ## Convert all NaN values to np.nan
                        state = _normalize_nans(state)

                        ## Set the state for the model class
                        try:
                            obj.__setstate__(state)
                        except Exception as e:
                            print(f"⚠️ Failed to set state for '{sub_grp_name}': {e}")
                            continue

                        ## If the object is a collapsed solutionstep from TransientSolution, expand back into solution for consistency
                        if isinstance(obj, model.classes.results.solutionstep) and sub_grp_name == "TransientSolution":
                            obj = _expand_step_to_solution(obj)

                        ## Normalize data types
                        _normalize_loaded_attributes(obj)

                        ## Assign the object to the model                        
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
                    classtype, obj = _get_class(grp)

                    # If obj is None, carry on (a warning is printed by create_instance in get_class)
                    if obj is None:
                        continue

                    ## Create empty state
                    state = {}

                    ## Get scalar attributes (those that are not stored as variables) and add to state
                    state = _get_attributes(state, grp)

                    ## Get variables
                    state = _get_variables(state, grp, grp_name)

                    ## Convert all NaN values to np.nan
                    state = _normalize_nans(state)

                    ## Set the state for the model class
                    try:
                        obj.__setstate__(state)
                    except Exception as e:
                        print(f"⚠️ Failed to set state for '{grp_name}': {e}")
                        continue

                    ## Normalize data types
                    _normalize_loaded_attributes(obj)

                    ## Assign the object to the model (e.g., md.mesh)
                    setattr(md, grp_name, obj)

                else:
                    print(f"⚠️️ classtype does not exist for group {grp_name}. Skipping...")
                    continue

    return md


def save_model(md, path):
    """
    Save an ISSM model to a NetCDF file.

    This function serializes an ISSM model and all its components to NetCDF format,
    including mesh, materials, geometry, boundary conditions, and results. The function
    handles nested objects and properly serializes the complete model state for later
    reconstruction.

    Parameters
    ----------
    md : model.Model
        The ISSM model object to be saved, containing all model components
        such as mesh, geometry, materials, boundary conditions, and results.

    path : str
        Path to the output NetCDF file where the model will be saved.

    Raises
    ------
    OSError
        If the output file cannot be created or written to.
    ValueError
        If the model contains unregistered classes that cannot be serialized.

    Notes
    -----
    - The function automatically handles different data types including scalars,
      arrays, lists, and nested objects
    - TransientSolution results are automatically collapsed from individual
      timesteps to a compact format for efficient storage
    - Boolean values are converted to integers for NetCDF compatibility
    - Object arrays and string arrays are properly handled with appropriate
      data type conversions
    - Compression is enabled for most variables to reduce file size

    See Also
    --------
    load_model : Load an ISSM model from NetCDF format

    Examples
    --------
    >>> save_model(md, 'my_model.nc')
    """

    # Helper function to convert character array to string for NetCDF writing
    def _char_array_to_strings(arr):
        arr = np.asarray(arr)  # Ensure it's a NumPy array
        if arr.ndim == 1:
            # Convert 1D string to single byte string array
            return np.array(["".join(arr.astype(str))], dtype='S')
        elif arr.ndim == 2:
            # Convert 2D strings to multiple byte string array
            return np.array(["".join(row.astype(str)) for row in arr], dtype='S')
        else:
            raise ValueError("Input must be a 1D or 2D char array with dtype='S1'")
        
    def _serialize_scalars(value, group, attr_name):
        if isinstance(value, (int, float, str, bool)):
            ## If it's a boolean, convert to int for NetCDF writing
            if isinstance(value, (bool)):
                value = int(value)
            group.setncattr(attr_name, value)
            return True
        return False
    
    def _serialize_arrays_lists(value, group, attr_name):
        # Handle arrays and lists (convert lists to arrays)
        if isinstance(value, np.ndarray) or isinstance(value, list):

            # If it's a list, convert to an array
            if isinstance(value, list) and len(value) == 0:
                # Handle empty list case
                group.setncattr(attr_name, "__EMPTY_LIST__")
                    
            if isinstance(value, list):
                value = np.array(value, dtype='S')
            # Otherwise, check the array type
            else:
                # If it's an object array, convert to a string array (object arrays can't be written to NetCDF)
                if value.dtype == object:
                    value = np.array(value, dtype='S')
                # Special handling for 'S1' datatype (these come from NetCDF Char variables when output from MATLAB)
                elif value.dtype.kind == 'S':
                    value = _char_array_to_strings(value)
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

            return True
        return False
    
    # Helper function to serialize an object's state
    def _serialize_object(obj, group, path = ''):
        """
        Serializes an object's state, including nested objects.
        """
        
        # If the object is None, skip it.
        if obj is None:
            print(f"⚠️  {path} is NoneType. This will be skipped from saving and re-initialised (default) when re-loading the NetCDF file.")
            return

        # Get state from the object
        state = None
        if hasattr(obj, '__getstate__'):
            state = obj.__getstate__()
        elif hasattr(obj, '__dict__'):
            state = obj.__dict__

        ## Is state is still None (failed above conditions), or not a dictionary, skip it
        if state is None or not isinstance(state, dict):
            print(f"⚠️  {path} cannot be serialized. This will be skipped from saving and re-initialised (default) when re-loading the NetCDF file.")
            return

        ## For each item, write attributes and variables...
        for attr_name, value in state.items():

            # Build path of item
            full_path = f"{path}.{attr_name}" if path else attr_name

            # Handle scalars
            if _serialize_scalars(value, group, attr_name):
                continue

            # Handle arrays and lists
            if _serialize_arrays_lists(value, group, attr_name):
                continue

            # Handle dict / OrderedDict
            if isinstance(value, collections.abc.Mapping):
                nested_group = group.createGroup(attr_name)

                for k, v in value.items():
                    key_name = str(k)
                    key_path = f"{full_path}.{key_name}"

                    if _serialize_scalars(v, nested_group, key_name):
                        continue

                    if _serialize_arrays_lists(v, nested_group, key_name):
                        continue

                    if isinstance(v, collections.abc.Mapping) or hasattr(v, '__getstate__') or hasattr(v, '__dict__'):
                        sub_group = nested_group.createGroup(key_name)
                        _serialize_object(v, sub_group, key_path)
                        continue

                    if v is None:                        
                        print(f"⚠️  {key_path} is NoneType. This will be skipped from saving and re-initialised (default) when re-loading the NetCDF file.")
                        continue

                    print(f"⚠️  Skipping unsupported field in Mapping: {key_path} ({type(v).__name__}). This will be skipped from saving and re-initialised (default) when re-loading the NetCDF file.")
                continue
                
            # Nested objects
            if hasattr(value, '__getstate__') or hasattr(value, '__dict__'):
                nested_group = group.createGroup(attr_name)
                _serialize_object(value, nested_group, full_path)
                continue
                
            # None / unsupported
            if value is None:
                print(f"⚠️  {full_path} is NoneType. This will be skipped from saving and re-initialised (default) when re-loading the NetCDF file.")
                continue
            
            print(f"⚠️  Skipping unsupported field: {full_path} ({type(value).__name__}). This will be skipped from saving and re-initialised (default) when re-loading the NetCDF file.")
            
    def _get_registered_name(obj):
        classname = obj.__class__
        matching_keys = [k for k, v in model.classes.class_registry.CLASS_REGISTRY.items() if v is classname]
        if not matching_keys:
            raise ValueError(f"Class {classname} is not registered.")
        registered_name = min(matching_keys, key=len)
        return registered_name

    # Create the NetCDF file
    with nc.Dataset(path, 'w', format='NETCDF4') as ds:

        ## Initialise dictionary to track existing dimension sizes & names
        defined_dimensions = {}

        ## Loop through model attributes (top-level groups)
        for name, obj in vars(md).items():
            ## Handle 'results' group specially
            if name == "results":
                results_group = ds.createGroup("results")

                ## Loop through each solution type in md.results
                for solution_name, solution_obj in vars(md.results).items():
                    if solution_obj is None:
                        print(f"⚠️  Skipping solution type: {solution_name})")
                        continue

                    ## Create subgroup for this solution (e.g., TransientSolution)
                    solution_group = results_group.createGroup(solution_name)

                    ## If it's a TransientSolution, collapse solution to a single step for writing
                    if isinstance(solution_obj, model.classes.results.solution) and solution_name == "TransientSolution":
                        solution_obj = _collapse_solution_to_step(solution_obj)

                    ## Attach class type metadata
                    classname = _get_registered_name(solution_obj)
                    solution_group.setncattr("classtype", classname)

                    _serialize_object(solution_obj, solution_group)

            else:
                ## For regular model components (e.g., mesh, materials, geometry)
                if obj is None:
                    continue

                ## Create group for the model component
                group = ds.createGroup(name)

                ## Attach class type metadata
                classname = _get_registered_name(obj)
                group.setncattr("classtype", classname)

                ## Serialize the component state
                _serialize_object(obj, group, path = name)

def _collapse_solution_to_step(solution):
    """
    Collapse a solution object with multiple timesteps into a single solutionstep for NetCDF storage.
    
    This function combines multiple solutionstep objects from a solution into a single
    solutionstep object with time-stacked arrays and consolidated metadata. This is used
    to efficiently store transient solutions in NetCDF format by reducing the complexity
    of nested timestep structures.
    
    Parameters
    ----------
    solution : model.classes.results.solution
        A solution object containing a list of solutionstep instances, one for each
        timestep, with individual arrays and metadata that need to be consolidated.
    
    Returns
    -------
    model.classes.results.solutionstep
        A collapsed solutionstep object containing time-stacked arrays and metadata
        from all timesteps, suitable for NetCDF storage.
    
    Notes
    -----
    - Arrays from individual timesteps are stacked along axis 0 (time dimension)
    - Time and step metadata are collected into arrays
    - Scalar values that are consistent across timesteps are stored as single values
    - Arrays with consistent shapes across timesteps are stacked into higher-dimensional arrays
    - Inconsistent data falls back to list storage
    - Empty solutions return an empty solutionstep
    - This function is the inverse operation of `_expand_step_to_solution`
    
    See Also
    --------
    _expand_step_to_solution : The inverse operation that expands a collapsed 
        solutionstep back into a solution with multiple timesteps
    """
    
    # Create empty solutionstep
    step = model.classes.results.solutionstep()
    steps = solution.steps or []

    # Skip if there are no steps
    if not steps:
        return step

    # Collect all fields across all steps
    all_fields = set()
    for s in steps:
        all_fields.update(s.__dict__.keys())

    # Loop over all fields and consolidate
    for field in sorted(all_fields):
        values = [getattr(s, field, None) for s in steps]
        non_none = [v for v in values if v is not None]

        # Skip if all values are none
        if not non_none:
            continue

        # Handle time and steps
        if field.lower() in ['time', 'step']:
            arr = np.array(non_none, dtype = float if field.lower() == 'time' else int)
            setattr(step, field, arr)
            continue

        # Handle single consistent scalar value
        if len(non_none) == 1:
            setattr(step, field, non_none[0])
            continue

        # Handle arrays with consistent shapes
        if all(isinstance(v, np.ndarray) for v in non_none):
            shapes = {v.shape for v in non_none}
            if len(shapes) == 1:
                stacked = np.stack(non_none, axis=0)
                setattr(step, field, stacked.squeeze())
                continue

        # Fallback: store as list
        setattr(step, field, values)

    return step

def _expand_step_to_solution(step_obj):
    """
    Expand a collapsed solutionstep back into a solution object with multiple timesteps.
    
    This function reverses the collapse operation performed by `_collapse_solution_to_step`,
    converting a single solutionstep object (with time-stacked arrays) back into a solution
    object containing a list of individual solutionstep instances for each timestep.
    
    Parameters
    ----------
    step_obj : model.classes.results.solutionstep
        A collapsed solutionstep object containing time-stacked arrays and metadata
        from multiple timesteps that were previously collapsed for NetCDF storage.
    
    Returns
    -------
    model.classes.results.solution
        A solution object containing a list of solutionstep instances, one for each
        timestep, with arrays split along the time axis and metadata appropriately 
        distributed.
    
    Notes
    -----
    - Time-varying arrays are split along axis 0 (time dimension)
    - Scalar values are propagated to all timesteps, except for metadata fields
      like 'SolutionType' which are assigned only to the first timestep
    - Arrays with singleton dimensions are flattened to match original format
    - If no time information is available, creates a single-timestep solution
    - This function is the inverse operation of `_collapse_solution_to_step`
    
    See Also
    --------
    _collapse_solution_to_step : The inverse operation that collapses multiple 
        timesteps into a single solutionstep for NetCDF storage
    """

    # Determine number of timesteps
    nt = None
    if hasattr(step_obj, 'time'):
        if isinstance(step_obj.time, (np.ndarray, list)):
            nt = len(step_obj.time)
        else:
            # Scalar time -- single step
            nt = 1
    else:
        # Fallback: find first array with first dimension >= 1
        for v in step_obj.__dict__.values():
            if isinstance(v, np.ndarray) and v.ndim >= 1:
                nt = v.shape[0]
                break

    if nt is None:
        # No time dimension -- set on step
        nt = 1

    # Create empty solution and steps
    sol = model.classes.results.solution([])
    sol.steps = [model.classes.results.solutionstep() for _ in range(nt)]

    # Loop over all fields in step_obj
    for field, value in step_obj.__dict__.items():
        if value is None:
            continue

        # Handle time / step fields
        if field.lower() in ['time', 'step']:
            if isinstance(value, (np.ndarray, list)):
                for t in range(nt):
                    val = value[t]
                    # Flatten single arrays
                    if isinstance(val, np.ndarray) and val.shape == (1,):
                        val = val[0]
                    sol.steps[t].__dict__[field] = val
            else:
                # Single value -- assign to first step
                sol.steps[0].__dict__[field] = value
            continue

        # Handle arrays of shape (nt, ...)
        if isinstance(value, np.ndarray):
            if value.shape[0] == nt:
                for t in range(nt):
                    val_t = value[t]
                    # Flatten single dimensions (0D arrays or [1,] arrays)
                    if isinstance(val_t, np.ndarray):
                        if val_t.ndim == 0:
                            val_t = val_t.item()
                        elif val_t.shape == (1,):
                            val_t = val_t[0]
                    sol.steps[t].__dict__[field] = val_t
                continue

        # Handle lists of length nt
        if isinstance(value, list) and len(value) == nt:
            for t in range(nt):
                val_t = value[t]
                # Flatten single arrays
                if isinstance(val_t, np.ndarray) and val_t.shape == (1,):
                    val_t = val_t[0]
                sol.steps[t].__dict__[field] = val_t
            continue

        # Handle scalar values
        if isinstance(value, (int, float, str, bool)):
            # Metadata like SolutionType: assign only to first step
            if field.lower() == 'solutiontype':
                sol.steps[0].__dict__[field] = value
            else:
                for t in range(nt):
                    sol.steps[t].__dict__[field] = value
            continue

        # Fallback: assign entire value to first step
        sol.steps[0].__dict__[field] = value

    return sol

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
        be used located in `../assets/default_variable_map.csv` relative to this script.
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
        variable_map = os.path.join(os.path.dirname(__file__), '../assets/default_variable_map.csv')
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
        if 'TransientSolution' in var_map['issmModelSubgroup'].values and tools.general.has_nested_attr(md, 'results', 'TransientSolution'):
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
                    variable = tools.general.convert_units(row['issmVariableUnit'], row['outputVariableUnit'], variable)


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

def issm_scp_out(host,
                 path,
                 login,
                 port,
                 packages):
    """
    Copy files to a remote host using SCP or create symbolic links for local transfers.
    This function transfers files either by creating symbolic links (for local transfers 
    where host matches hostname) or by using the SCP protocol for remote transfers. 
    For remote transfers, it attempts standard SCP first and falls back to legacy SSH 
    options if the initial attempt fails.

    Parameters
    ----------
    host : str
        Target hostname or IP address for file transfer.
    path : str
        Destination directory path on the target host.
    login : str
        Username for authentication on the remote host.
    port : int or None
        SSH port number for connection. If None, uses default port 22.
    packages : list of str
        List of file/package names to transfer from current working directory.

    Raises
    ------
    Exception
        If SCP transfer fails after attempting both standard and legacy options.

    Warnings
    --------
    UserWarning
        When a package file does not exist and will be skipped during local transfer.

    Notes
    -----
    For local transfers (same hostname), the function:
    - Creates symbolic links in the destination path
    - Removes existing files with same names before linking
    - Skips non-existent packages with warnings
    For remote transfers, the function:
    - Attempts standard SCP first
    - Falls back to SCP with legacy SSH options (-OT) if standard fails
    - Supports custom port specification with -P flag
    Examples
    --------
    >>> # Local transfer
    >>> issm_scp_out('localhost', '/tmp/dest', 'user', None, ['file1.txt', 'file2.dat'])
    >>> # Remote transfer with custom port
    >>> issm_scp_out('remote.server.com', '/home/user/data', 'username', 2222, ['data.bin'])
    """
    
    # Get hostname
    hostname = tools.config.get_hostname()

    # If host and hostname are the same, do a simple copy
    if host.lower() == hostname.lower():
        for package in packages:
            
            ## Check package exists
            if os.path.exists(package):

                ### Get current working directory
                pwd = os.getcwd()
                
                ### Change the current working directory to the target path
                os.chdir(path)
                try:
                    ### Remove any existing file with the same name in the target directory
                    os.remove(package)
                except OSError:
                    ### Ignore errors if the file doesn't exist or can't be removed
                    pass
                ### Create a symbolic link from the package in the original directory to the target path
                subprocess.call('ln -s %s %s' % (os.path.join(pwd, package), path), shell=True)
                ### Change back to the original working directory
                os.chdir(pwd)

            else:
                ## If the package does not exist, print a warning
                warnings.warn(f'pyissm.model.io.issm_scp_out: {package} does not exist and will be skipped')

    # If this is not a local machine, use scp to transfer the files
    else:
        ## Get current working directory and build full paths to package files
        pwd = os.getcwd()
        file_list = [os.path.join(pwd, x) for x in packages]
        file_list_str = ' '.join([str(x) for x in file_list])
        
        ## Handle scp with custom port
        if port:
            ## First attempt: try scp with custom port
            subproc_cmd = 'scp -P {} {} {}@{}:{}'.format(port, file_list_str, login, host, path)
            subproc = subprocess.Popen(subproc_cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE, universal_newlines = True)
            outs, errs = subproc.communicate()
            
            ## If first attempt failed, try with legacy SSH options
            if errs != '':
                subproc_cmd = 'scp -OT -P {} {} {}@{}:{}'.format(port, file_list_str, login, host, path)
                subproc = subprocess.Popen(subproc_cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE, universal_newlines = True)
                outs, errs = subproc.communicate()
        else:
            ## Handle scp with default port (22)
            ## First attempt: try standard scp
            subproc_cmd = 'scp {} {}@{}:{}'.format(file_list_str, login, host, path)
            subproc = subprocess.Popen(subproc_cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE, universal_newlines = True)
            outs, errs = subproc.communicate()
            
            ## If first attempt failed, try with legacy SSH options
            if errs != '':
                subproc_cmd = 'scp -OT {} {}@{}:{}'.format(file_list_str, login, host, path)
                subproc = subprocess.Popen(subproc_cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE, universal_newlines = True)
                outs, errs = subproc.communicate()
            
        ## Check scp worked
        if errs != '':
            raise Exception(f'pyissm.model.io.issm_scp_out: scp failed with the following error: {errs}')

def issm_ssh(host,
             login,
             port,
             command):

    """
    Execute a command on a remote host via SSH or locally if on the same machine.
    
    This function determines whether to run a command locally or remotely based on
    hostname comparison. For remote execution, it uses platform-specific SSH clients:
    plink.exe on Windows and standard ssh on Mac/Linux. Includes a workaround for
    macOS file descriptor blocking issues.
    
    Parameters
    ----------
    host : str
        The hostname or IP address of the target machine.
    login : str
        The username for SSH authentication.
    port : int or None
        The SSH port number. If None, uses default SSH port (22).
    command : str
        The command to execute on the target machine.
    
    Notes
    -----
    - On Windows, requires plink.exe in the ISSM external packages directory
    - Prompts for username and password interactively on Windows
    - On macOS, applies file descriptor flags to prevent "Resource temporarily 
      unavailable" errors
    - Uses shell=True for subprocess calls, which may have security implications
    """
        
    # Get hostname
    hostname = tools.config.get_hostname()

    # If host and hostname are the same, just run the command
    if host.lower() == hostname.lower():
        subprocess.call(command, shell = True)
    
    # If this is not a local machine, use ssh to run the command
    else:
        ## Windows requires plink.exe for ssh
        if tools.config.is_pc():
            issm_dir = tools.config.get_issm_dir()

            username = eval(input('Enter your username: '))
            key = eval(input('Enter your key: '))

            subprocess.call('%s/externalpackages/ssh/plink.exe-ssh -l "%s" -pw "%s" %s "%s"' % (issm_dir, username, key, host, command), shell = True)
        ## Mac/Linux use standard ssh
        else:
            if port:
                subprocess.call('ssh -l %s -p %d localhost "%s"' % (login, port, command), shell=True)
            else:
                subprocess.call('ssh -l {} {} "{}"'.format(login, host, command), shell=True)

    # "IOError: [Errno 35] Resource temporarily unavailable"
    # on the Mac when trying to display md after the solution.
    # (from http://code.google.com/p/robotframework/issues/detail?id=995)
    if sys.platform == 'darwin':
        import fcntl

        fd = sys.stdin.fileno()
        flags = fcntl.fcntl(fd, fcntl.F_GETFL)
        fcntl.fcntl(fd, fcntl.F_SETFL, flags & ~os.O_NONBLOCK)

        fd = sys.stdout.fileno()
        flags = fcntl.fcntl(fd, fcntl.F_GETFL)
        fcntl.fcntl(fd, fcntl.F_SETFL, flags & ~os.O_NONBLOCK)

def issm_scp_in(host,
                login,
                port,
                path,
                packages):
    """
    Transfer files from a remote host to the current working directory using SCP or local copy.

    This function transfers specified packages (files) from a remote path to the current
    working directory. If the host is the same as the local hostname, it performs a local
    copy operation. Otherwise, it uses SCP (Secure Copy Protocol) to transfer files from
    the remote host.
    
    Parameters
    ----------
    host : str
        The hostname or IP address of the remote host.
    login : str
        The username for authentication on the remote host.
    port : int or None
        The SSH port number for the remote connection. If None, uses default port 22.
    path : str
        The remote directory path where the packages are located.
    packages : list of str
        List of filenames to transfer from the remote host.
    
    Raises
    ------
    Exception
        If the SCP command fails with an error.
    OSError
        If a package does not exist after the transfer operation.
    
    Warnings
    --------
    UserWarning
        If a package does not exist on the local host during local copy operation.
    
    Notes
    -----
    The function first attempts standard SCP commands. If those fail, it retries with
    legacy SSH options (-OT flags) to handle compatibility issues with different SSH
    configurations.
    For local transfers (when host matches local hostname), the function uses shutil.copy
    and ignores any OSError exceptions that may occur during the copy operation.
    """
      
    # Get hostname
    hostname = tools.config.get_hostname()

    # If host and hostname are the same, do a simple copy
    if host.lower() == hostname.lower():
        
        for package in packages:

            ## Check package exists
            if os.path.exists(os.path.join(path, package)):

                ### Get current working directory
                pwd = os.getcwd()
                
                try:
                    shutil.copy(os.path.join(path, package), pwd)
                except OSError:
                    # Ignore errors
                    pass

            else:
                ## If the package does not exist, print a warning
                warnings.warn(f'pyissm.model.io.issm_scp_in: {package} does not exist and will be skipped')

    # If this is not a local machine, use scp to transfer the files
    else:

        ## Get current working directory
        pwd = os.getcwd()
        
        ## Transfer files individually to handle missing errlog files gracefully
        for package in packages:
            ## Build remote file path with proper host prefix
            remote_file = f'{login}@{host}:{os.path.join(path, package)}'
            
            ## Handle scp with custom port
            if port:
                ## First attempt: try scp with custom port
                subproc_cmd = 'scp -P {} {} {}'.format(port, remote_file, pwd)
                subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                outs, errs = subproc.communicate()
                
                ## If first attempt failed, try with legacy SSH options
                if errs != '':
                    subproc_cmd = 'scp -OT -P {} {} {}'.format(port, remote_file, pwd)
                    subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                    outs, errs = subproc.communicate()
            else:
                ## Handle scp with default port (22)
                ## First attempt: try standard scp
                subproc_cmd = 'scp {} {}'.format(remote_file, pwd)
                subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                outs, errs = subproc.communicate()
                
                ## If first attempt failed, try with legacy SSH options
                if errs != '':
                    subproc_cmd = 'scp -OT {} {}'.format(remote_file, pwd)
                    subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                    outs, errs = subproc.communicate()
            
            ## Check if transfer succeeded
            if errs != '':
                ## Check if this is an errlog file - these are optional
                if package.endswith('.errlog'):
                    warnings.warn(f'pyissm.model.io.issm_scp_in: {package} was not transferred (likely no errors occurred)')
                else:
                    raise Exception(f'pyissm.model.io.issm_scp_in: scp failed for {package} with error: {errs}')
            
            ## Verify file was transferred (only for successful transfers)
            if errs == '' and not os.path.exists(os.path.join('.', package)):
                if package.endswith('.errlog'):
                    warnings.warn(f'pyissm.model.io.issm_scp_in: {package} was not transferred (likely no errors occurred)')
                else:
                    raise OSError(f'pyissm.model.io.issm_scp_in: {package} does not exist after transfer')