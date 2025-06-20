"""
ISMIP-related functions for ISSM

This module contains functions for handling ISMIP-related data and operations.
"""


from .. import utils

def calc_perc_ice_cover(total_ice_area, ice_area):
    """
    Calculate the percentage ice cover, relative to the of total ice area.

    This function computes the fraction of the ice area that is grounded (i.e.,
    not floating). The result is returned as a percentage.

    Parameters
    ----------
    total_ice_area : ndarray
        Array representing the total ice area.

    ice_area : ndarray
        Array representing the ice area for percentage coverage to be computed.

    Returns
    -------
    ice_area_perc : ndarray
        Array with the same shape as `ice_area`, where each
        element represents the percentage of the total_ice_area that is covered by the ice_area.

    """

    ## Calculate ice area fraction
    ice_area_perc = (ice_area / total_ice_area) * 100

    return ice_area_perc


def get_ismip_variable(md, ismip_variable_name):
    """
    Retrieve a specific ISMIP variable from a model object.

    This function extracts and calculates a requested ISMIP (Ice Sheet Model Intercomparison Project)
    variable from a given model object, if the necessary data is available. It supports area-based
    variables calculated from grounded and floating ice components stored in the model's
    transient solution.

    Parameters
    ----------
    md : object
        The model object containing transient solution results, typically with
        attributes like `md.results.TransientSolution.GroundedArea` and `FloatingArea`.

    ismip_variable_name : str
        Name of the ISMIP variable to compute. Must be one of:
        - 'land_ice_area_fraction' (assumes all ice in ISSM is land ice and returns 100 %)
        - 'floating_ice_shelf_area_fraction'
        - 'grounded_ice_sheet_area_fraction'

    Returns
    -------
    variable : array-like or None
        The computed variable as a percentage, or `None` if the required data is
        not available in the model object.

    Raises
    ------
    None explicitly, but prints a warning if the requested variable cannot be computed.

    Notes
    -----
    - 'land_ice_area_fraction' assumes all ice is land ice (returns 100% where data is available).
    - 'floating_ice_shelf_area_fraction' computes the fraction of total ice that is floating.
    - 'grounded_ice_sheet_area_fraction' computes the fraction of total ice that is grounded.
    - Missing required attributes in `md` will result in the variable being skipped.
    """

    variable = None

    if ismip_variable_name == 'land_ice_area_fraction':
        if utils.has_nested_attr(md, 'results', 'TransientSolution', 'GroundedArea') and utils.has_nested_attr(md, 'results', 'TransientSolution', 'FloatingArea'):
            total_area = md.results.TransientSolution.GroundedArea + md.results.TransientSolution.FloatingArea
            variable = calc_perc_ice_cover(total_area, total_area)
        else:
            print("\033[1mland_ice_area_fraction\033[0m cannot be computed and will be skipped.")

    elif ismip_variable_name == 'floating_ice_shelf_area_fraction':
        if utils.has_nested_attr(md, 'results', 'TransientSolution', 'FloatingArea'):
            total_area = md.results.TransientSolution.GroundedArea + md.results.TransientSolution.FloatingArea
            variable = calc_perc_ice_cover(total_area, md.results.TransientSolution.FloatingArea)
        else:
            print("\033[1mfloating_ice_shelf_area_fraction\033[0m cannot be computed and will be skipped.")

    elif ismip_variable_name == 'grounded_ice_sheet_area_fraction':
        if utils.has_nested_attr(md, 'results', 'TransientSolution', 'GroundedArea'):
            total_area = md.results.TransientSolution.GroundedArea + md.results.TransientSolution.FloatingArea
            variable = calc_perc_ice_cover(total_area, md.results.TransientSolution.GroundedArea)
        else:
            print("\033[1mgrounded_ice_sheet_area_fraction\033[0m cannot be computed and will be skipped.")

    else:
        print(f"\033[1m{ismip_variable_name}\033[0m is not recognized and will be skipped.")

    return variable
