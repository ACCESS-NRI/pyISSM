"""
Functions for reading and writing data.

This module contains functions for reading and writing data to and from files.
"""

from types import SimpleNamespace
from ..model import Model
import netCDF4 as nc

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