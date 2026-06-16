"""
Diagnostic functions for ISSM models

This module contains various functions that support diagnostic analysis of ISSM models.
"""

import numpy as np
from . import interp
from .. import model

def weighted_mean(md,
                  values,
                  weights = None):
    """
    Compute weighted mean.

    Parameters
    ----------
    md : :class:`pyissm.Model`
        ISSM model.

    values : :class:`numpy.ndarray`
        Values to compute weighted mean. Either vertex-based or element-based.

    weights : :class:`numpy.ndarray`, optional
        Custom weights. Default = None, in which case element areas will be used.

    Returns
    -------
        :class:`float`
            Weighted mean value.
    """

    # Ensure values is a numpy array
    values = np.asarray(values)

    # If values are vertex-based, convert to element-based
    if values.shape[0] == md.mesh.numberofvertices:
        values = interp.vertex_to_element(md, values)

    # If weights are not provided, use element areas
    if weights is None:
        weights = model.mesh.get_element_areas_volumes(md.mesh.elements, md.mesh.x, md.mesh.y)
    
    # Ensure weights is a numpy array
    weights = np.asarray(weights)

    # Mask out non-finite values and corresponding weights
    mask = np.isfinite(values) & np.isfinite(weights)

    # Compute area-weighted mean
    wm = np.sum(values[mask] * weights[mask]) / np.sum(weights[mask])
    
    return wm 

def weighted_rmse(md,
                  values,
                  weights = None):
    """
    Compute weighted RMSE.

    Parameters
    ----------
    md : :class:`pyissm.Model`
        ISSM model.

    values : :class:`numpy.ndarray`
        Values to compute weighted RMSE.

    weights : :class:`numpy.ndarray`, optional.
        Custom weights. Default = None, in which case element areas will be used.

    Returns
    -------
    :class:`float`
        Weighted RMSE value.
    """

    # Ensure values is a numpy array
    values = np.asarray(values)

    # Compute weighted mean square error
    wmse = weighted_mean(md, values**2, weights = weights)

    # Return RMSE
    return np.sqrt(wmse)

def velocity_residuals(md,
                       obs_vel = None,
                       model_vel = None):
    """
    Compute velocity residuals.

    Parameters
    ----------
    md : :class:`pyissm.Model`
        ISSM model.
        
    obs_vel : :class:`numpy.ndarray`, optional
        Observed velocities. If None, will use md.inversion.vel_obs.

    model_vel : :class:`numpy.ndarray`, optional
        Model velocities. If None, will use md.results.StressbalanceSolution.Vel

    Returns
    -------
    :class:`numpy.ndarray`
        Velocity residuals.
    """

    # If obs_vel is not provided, extract from model inversion data
    if obs_vel is None:
        obs_vel = md.inversion.vel_obs

    # If model_vel is not provided, extract from model results
    if model_vel is None:
        model_vel = md.results.StressbalanceSolution.Vel

    # Calculate residuals
    residuals = model_vel - obs_vel

    return residuals

def velocity_rmse(md,
                  obs_vel = None,
                  model_vel = None):
    """
    Compute weighted velocity RMSE.

    Parameters
    ----------
    md : :class:`pyissm.Model`
        ISSM model.

    obs_vel : :class:`numpy.ndarray`, optional
        Observed velocities. If None, will use md.inversion.vel_obs.

    model_vel : :class:`numpy.ndarray`, optional
        Model velocities. If None, will use md.results.StressbalanceSolution.Vel

    Returns
    -------
    :class:`float`
        Weighted velocity RMSE.
    """

    # Compute velocity residuals
    residuals = velocity_residuals(md, obs_vel, model_vel)

    # Compute and return weighted RMSE of velocity residuals
    return weighted_rmse(md, residuals)