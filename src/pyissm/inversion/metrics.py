"""
Inversion metrics for pyISSM.

This module contains various functions for computing metrics related to inversion processes in ISSM models.
"""

from .. import model, tools
import numpy as np
import pandas as pd
import itertools


def cost_function_values(md):
    """
    Extract final inversion cost-function values.

    Parameters
    ----------
    md : :class: `pyissm.Model`
        ISSM model object containing inversion results.

    Returns
    -------
    :class:`dict`
        Dictionary containing final cost-function values for each cost function and total cost.
    """

    # Get cost function history from model results
    J = md.results.StressbalanceSolution.J

    # Get final cost function values (last entry in J)
    final_costs = J[-1]

    # Create output dictionary
    output = {}

    # Loop over cost functions and add to output
    for i, cf in enumerate(md.inversion.cost_functions):

        output[f"cost_{cf}"] = final_costs[i]

    # Get total cost (last entry in J)
    output["cost_total"] = final_costs[-1]

    return output

def cost_function_ratios(md):
    """
    Compute ratios between all inversion cost functions.

    Parameters
    ----------
    md : :class:`pyissm.Model`
        ISSM model object containing inversion results.

    Returns
    -------
    :class:`dict`
        Dictionary of pairwise cost-function ratios.
    """

    # Get individual cost values
    costs = cost_function_values(md)

    # Remove total cost
    cost_keys = [
        key for key in costs.keys()
        if key != "cost_total"
    ]

    # Create output dictionary
    output = {}

    # Generate all unique cost-function pairs
    for key1, key2 in itertools.combinations(cost_keys, 2):

        # Extract cost values for the pair
        value1 = costs[key1]
        value2 = costs[key2]

        # Extract numeric IDs (e.g. "cost_101" -> "101") for naming ratios
        cf1 = key1.replace("cost_", "")
        cf2 = key2.replace("cost_", "")

        # Construct ratio name (e.g. "ratio_101_103")
        ratio_name = f"ratio_{cf1}_{cf2}"

        # Calculate ratio and avoid divide-by-zero
        if value2 == 0:
            output[ratio_name] = np.nan
        else:
            output[ratio_name] = value1 / value2

    return output

def extract_convergence_history(md):
    """
    Extract convergence history of cost functions.

    Parameters
    ----------
    md : :class:`pyissm.Model`

    Returns
    -------
    :class:`pandas.DataFrame`
        DataFrame containing cost function history.
    """

    # Get cost function history from model results
    J = md.results.StressbalanceSolution.J

    # Create output dictionary
    output = {}

    # Loop over cost functions and add to output
    for i, cf in enumerate(md.inversion.cost_functions):

        # Add cost history for this cost function (column in J) to output dictionary
        output[f"cost_{cf}"] = J[:, i]

    # Add total cost history (last column in J)
    output["cost_total"] = J[:, -1]

    # Convert to DataFrame for easier analysis
    output = pd.DataFrame(output)

    # Add iteration numbers as first column
    output.insert(0, "iteration", np.arange(1, len(output) + 1))

    return output

def velocity_residual_area_metrics(md):
    """
    Compute positive/negative residual area metrics.

    Parameters
    ----------
    md : :class:`pyissm.Model`
        ISSM model object containing mesh and solution data.

    Returns
    -------
    :class:`dict`
        Dictionary containing the computed metrics.
    """

    # Compute velocity residuals at vertices
    residuals_vertex = tools.diagnostics.velocity_residuals(md)

    # Convert residuals to elements
    residuals_element = tools.interp.vertex_to_element(md, residuals_vertex)

    # Calculate elemtent areas
    element_areas = model.mesh.get_element_areas_volumes(md.mesh.elements, md.mesh.x, md.mesh.y)

    # Create masks for positive and negative residuals
    positive_mask = residuals_element > 0
    negative_mask = residuals_element < 0

    # Calculate total area of positive and negative residuals
    positive_area = np.sum(element_areas[positive_mask])
    negative_area = np.sum(element_areas[negative_mask])

    # Calculate total area for normalization
    total_area = np.sum(element_areas)

    # Calculate area fractions
    positive_fraction = positive_area / total_area
    negative_fraction = negative_area / total_area

    # Return metrics in a dictionary
    return {
        "positive_residual_area": positive_area,
        "negative_residual_area": negative_area,
        "positive_residual_fraction": positive_fraction,
        "negative_residual_fraction": negative_fraction,
    }

def field_smoothness_metrics(md,
                            field = None):
        """
        Compute field smoothness metrics.
    
        Parameters
        ----------
        md : :class:`pyissm.Model`
            ISSM model object containing mesh and (optionally) solution data.
        field : :class:`numpy.ndarray`
            Field values, defined on vertices. If None, will attempt to extract from model results based on control parameter.
    
        Returns
        -------
        :class:`dict`
            Dictionary containing smoothness metrics.
        """

        # Extract field from model if not provided
        if field is None:
            if md.inversion.control_parameters == ['MaterialsRheologyBbar']:
                field = md.results.StressbalanceSolution.MaterialsRheologyBbar
            elif md.inversion.control_parameters == ['FrictionCoefficient']:
                field = md.results.StressbalanceSolution.FrictionCoefficient
            elif md.inversion.control_parameters == ['FrictionC']:
                field = md.results.StressbalanceSolution.FrictionC
            else:
                raise ValueError("pyissm.inversion.metrics.field_smoothness_metrics: Unsupported control parameter for smoothness metrics")
        
        # Compute slope (magnitude of gradient)
        _, _, slope = tools.geometry.slope(md, field)

        # Compute area-weighted mean slope
        area_weighted_mean_slope = tools.diagnostics.weighted_mean(md, slope)

        # Compute area-weighted RMSE of slope
        area_weighted_rmse_slope = tools.diagnostics.weighted_rmse(md, slope)
    
        return {
            "mean_gradient_magnitude": area_weighted_mean_slope,
            "rmse_gradient_magnitude": area_weighted_rmse_slope,
        }