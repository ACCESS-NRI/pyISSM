"""
Functions for parameterising ISSM models.
"""

import numpy as np
import os
from .. import model
from .. import utils

def set_mask(md,
             floating_ice_name = '',
             grounded_ice_name = 'all',
             ice_domain = None,
             **kwargs):

    """
    Set the mask of a model based on a floating ice region and a grounded ice region.

    Parameters
    ----------
    md : ISSM model object
        The model to set the mask on.
    floating_ice_name : str, optional
        The name of the floating ice region. Default is '', which means no 
        floating ice region is set.
    grounded_ice_name : str, optional
        The name of the grounded ice region. Default is 'all', which means all 
        grounded ice is set.
    ice_domain : str, optional
        Path to file defining the ice domain contour. If provided, used to define
        the ice levelset field. Default is None.
    **kwargs : dict
        Additional keyword arguments to pass to the flag_elements function.

    Returns
    -------
    model
        The modified model with updated mask fields.

    Raises
    ------
    FileNotFoundError
        If ice_domain file path is provided but the file does not exist.
    """

    # Get mesh information
    elements = md.mesh.elements
    x = md.mesh.x
    y = md.mesh.y

    # Flag elements on floating and grounded ice
    elements_on_floating_ice = model.mesh.flag_elements(md, region = floating_ice_name, **kwargs)
    elements_on_grounded_ice = model.mesh.flag_elements(md, region = grounded_ice_name, **kwargs)

    # Ensure mutual exclusivity of floating and grounded ice elements (contours can intersect).
    # When there are overlaps, grounded wins; grounded = not floating ice.
    elements_on_floating_ice = elements_on_floating_ice & ~elements_on_grounded_ice
    elements_on_grounded_ice = ~elements_on_floating_ice    

    # Get vertices for floating and grounded ice elements
    vertex_on_floating_ice = np.zeros(md.mesh.numberofvertices, dtype = bool)
    vertex_on_grounded_ice = np.zeros(md.mesh.numberofvertices, dtype = bool)
    vertex_on_grounded_ice[md.mesh.elements[elements_on_grounded_ice].ravel() - 1] = True
    vertex_on_floating_ice[md.mesh.elements[elements_on_floating_ice].ravel() - 1] = True

    # Populate the ocean levelset field
    md.mask.ocean_levelset = -1. * np.ones(md.mesh.numberofvertices)
    md.mask.ocean_levelset[md.mesh.elements[elements_on_grounded_ice].ravel() - 1] = 1.

    # Populate the ice levelset field
    ## If an ice domain is provided, use it to define the levelset
    if ice_domain is not None:

        ### Check that the file exists
        if not os.path.isfile(ice_domain):
            raise FileNotFoundError(f"Ice domain file {ice_domain} not found.")
        
        ### Initialize ice_levelset to 1 (no ice)
        md.mask.ice_levelset = np.ones(md.mesh.numberofvertices)
        
        ### Get vertices inside the ice domain contour and set ice_levelset to -1 (ice)
        vertex_inside = utils.wrappers.ContourToMesh(elements, x, y, ice_domain, 'node', 1)
        print('type of vertex_inside:', type(vertex_inside))
        print(vertex_inside)
        md.mask.ice_levelset[vertex_inside.astype(bool)] = -1.
    
    ## Otherwise, set ice_levelset to -1 (all ice)
    else:
        md.mask.ice_levelset = -1. * np.ones(md.mesh.numberofvertices)

    return md



