"""
Functions for parameterising ISSM models.
"""

import numpy as np
import os
import warnings
from pyissm import model, tools

def set_mask(md,
             floating_ice_name = None,
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
        if tools.wrappers.check_wrappers_installed():
            vertex_inside = tools.wrappers.ContourToMesh(elements, x, y, ice_domain, 'node', 1)
            md.mask.ice_levelset[vertex_inside.astype(bool)] = -1.
        else:
            raise ImportError("pyissm.tools.param.set_mask: Wrappers are not installed. Cannot use ice_domain option.")
    
    ## Otherwise, set ice_levelset to -1 (all ice)
    else:
        md.mask.ice_levelset = -1. * np.ones(md.mesh.numberofvertices)

    return md

def set_flow_equation(md,
                      SIA = None,
                      SSA = None,
                      HO = None,
                      L1L2 = None,
                      MOLHO = None,
                      FS = None,
                      fill = None,
                      coupling = 'tiling',
                      **kwargs):

    """
    Set flow equation for model    
    """

    # Error checks
    if coupling.lower() not in ['tiling', 'penalties']:
        raise ValueError("pyissm.tools.param.set_flow_equation: Coupling must be 'tiling' or 'penalties'.")
    
    if fill is not None:
        if fill.lower() not in ['sia', 'ssa', 'ho', 'l1l2', 'molho', 'fs']:
            raise ValueError("pyissm.tools.param.set_flow_equation: Fill must be one of: None, 'SIA', 'SSA', 'HO', 'L1L2', 'MOLHO', or 'FS'.")
    
    # Get flags for each flow equation
    sia_flag = model.mesh.flag_elements(md, SIA)
    ssa_flag = model.mesh.flag_elements(md, SSA)
    ho_flag = model.mesh.flag_elements(md, HO)
    l1l2_flag = model.mesh.flag_elements(md, L1L2)
    molho_flag = model.mesh.flag_elements(md, MOLHO)
    fs_flag = model.mesh.flag_elements(md, FS)
    none_flag = np.zeros(md.mesh.numberofelements, dtype = bool)

    # If fill is specified, fill unassigned elements with the specified flow equation
    if fill is not None:
        if fill.lower() == 'sia':
            sia_flag = ~ssa_flag & ~ho_flag
        elif fill.lower() == 'ssa':
            ssa_flag = ~sia_flag & ~ho_flag & ~fs_flag
        elif fill.lower() == 'ho':
            ho_flag = ~sia_flag & ~ssa_flag & ~fs_flag
    
    # Check that all elements only have one (compatible) flow equation assigned
    flag = [sia_flag, ssa_flag, ho_flag, l1l2_flag, molho_flag, fs_flag]
    
    ## Check all elements have been assigned a flow equation
    if not np.all(np.logical_or.reduce(flag)):
        raise ValueError("pyissm.tools.param.set_flow_equation: One or more elements have not been assigned a flow equation.")    
    
    ## Check no elements have been assigned multiple flow equations
    if np.any(np.sum(flag, axis=0) > 1):
        raise ValueError("pyissm.tools.param.set_flow_equation: One or more elements have been assigned multiple flow equations.")

    ## Check that L1L2, MOLHO, and FS are not coupled to any other model for now
    if np.any(l1l2_flag) and np.any(sia_flag | ssa_flag | ho_flag | fs_flag | molho_flag):
        raise ValueError('pyissm.tools.param.set_flow_equation: L1L2 cannot be coupled to other flow equations.')
    if np.any(molho_flag) and np.any(sia_flag | ssa_flag | ho_flag | fs_flag | l1l2_flag):
        raise ValueError('pyissm.tools.param.set_flow_equation: MOLHO cannot be coupled to other flow equations.')
    if np.any(fs_flag) & np.any(sia_flag | ssa_flag | ho_flag | l1l2_flag | molho_flag):
        raise ValueError('pyissm.tools.param.set_flow_equation: FS cannot be coupled to other flow equations.')
        
    ##  Check HO or FS are not used for a 2D mesh
    if md.mesh.domain_type() == '2Dhorizontal':
        if np.any(ho_flag):
            raise ValueError('pyissm.tools.param.set_flow_equation: HO cannot be used for a 2D mesh. Extrude it first.')
        if np.any(fs_flag):
            raise ValueError('pyissm.tools.param.set_flow_equation: FS cannot be used for a 2D mesh. Extrude it first.')
    
    # Initialize flow equation fields
    sia_node = np.zeros(md.mesh.numberofvertices, dtype = bool)
    ssa_node = np.zeros(md.mesh.numberofvertices, dtype = bool)
    ho_node = np.zeros(md.mesh.numberofvertices, dtype = bool)
    l1l2_node = np.zeros(md.mesh.numberofvertices, dtype = bool)
    molho_node = np.zeros(md.mesh.numberofvertices, dtype = bool)
    fs_node = np.zeros(md.mesh.numberofvertices, dtype = bool)

    # Populate flow equation fields
    sia_node[md.mesh.elements[sia_flag].ravel() - 1] = True
    ssa_node[md.mesh.elements[ssa_flag].ravel() - 1] = True
    ho_node[md.mesh.elements[ho_flag].ravel() - 1] = True
    l1l2_node[md.mesh.elements[l1l2_flag].ravel() - 1] = True
    molho_node[md.mesh.elements[molho_flag].ravel() - 1] = True

    # Handle FS separately
    ## Modify fs_flag to remove elements that are constrained everywhere (spc + border with HO or SSA)
    if np.any(fs_flag):
        
        ## Find all the nodes on the boundary of the domain without icefront
        full_spc_nodes = np.logical_or(~np.isnan(md.stressbalance.spcvx)
                                       & ~np.isnan(md.stressbalance.spcvy)
                                       & ~np.isnan(md.stressbalance.spcvz),
                                       np.logical_and(ho_node, fs_node))
        
        ## Find all the elements on the boundary of the domain without icefront
        full_spc_elements = np.sum(full_spc_nodes[md.mesh.elements - 1], axis=1) == 6

        fs_flag[np.where(full_spc_elements.reshape(-1))] = False

        fs_node[md.mesh.elements[fs_flag].ravel() - 1] = True

    ## Complete with NoneApproximateion or the other model if there is no FS
    if any(fs_flag):

        ### Fill with HO if possible
        if any(ho_flag):
            ho_flag[~fs_flag] = True
            ho_node[md.mesh.elements[ho_flag].ravel() - 1] = True
        elif any(ssa_flag):
            ssa_flag[~fs_flag] = True
            ssa_node[md.mesh.elements[ssa_flag].ravel() - 1] = True
        else:
            none_flag[~fs_flag] = True

    # Complete coupling
    if coupling.lower() == 'tiling':
        md.stressbalance.vertex_pairing = np.array([])
    
    ssaho_node = np.zeros(md.mesh.numberofvertices, dtype = bool)
    hofs_node = np.zeros(md.mesh.numberofvertices, dtype = bool)
    ssafs_node = np.zeros(md.mesh.numberofvertices, dtype = bool)
    ssaho_flag = np.zeros(md.mesh.numberofelements, dtype = bool)
    hofs_flag = np.zeros(md.mesh.numberofelements, dtype = bool)
    ssafs_flag = np.zeros(md.mesh.numberofelements, dtype = bool)

    if coupling.lower() == 'penalties':
        ## Create border between HO and SSA
        num_nodes_2d = md.mesh.numberofvertices2d
        num_layers = md.mesh.numberoflayers

        ### Nodes connected to two different types of elements
        border_nodes_2d = np.where(np.logical_and(ho_node[0:num_nodes_2d], ssa_node[0:num_nodes_2d]))[0] + 1

        ## Initialise and populate penalties structure
        if np.all(np.logical_not(np.isnan(border_nodes_2d))):
            penalties = np.zeros((0, 2))
            for i in range(1, num_layers):
                penalties = np.vstack((penalties, np.vstack((border_nodes_2d, border_nodes_2d + md.mesh.numberofvertices2d * (i))).T))
            md.stressbalance.vertex_pairing = penalties

    elif coupling.lower() == 'tiling':

        ## Couple SSA and HO
        if any(ssa_flag) and any(ho_flag):
            ### Find border nodes
            ssaho_node[ssa_node & ho_node] = True

            ### SSA elements in contact with this layer become SSAHO elements
            matrix_elements = ssaho_node[md.mesh.elements - 1]
            common_elements = np.sum(matrix_elements, axis=1) != 0
            common_elements[ho_flag] = False
            ssa_flag[common_elements] = False
            ssaho_flag[common_elements] = True
            ssa_node[:] = False
            ssa_node[md.mesh.elements[ssa_flag].ravel() - 1] = True
            
            ### Rule out elements that don't touch the 2 boundaries
            pos = hofs_flag.nonzero()[0]
            elist = (
                np.sum(fs_node[md.mesh.elements[pos, :] - 1], axis=1).astype(int)
                - np.sum(ho_node[md.mesh.elements[pos, :] - 1], axis=1).astype(int)
            )

            ### Update flags based on elist values
            pos1 = pos[elist == 1]
            fs_flag[pos1] = True
            hofs_flag[pos1] = False

            pos2 = pos[elist == -1]
            ho_flag[pos2] = True
            hofs_flag[pos2] = False

            ### Recompute nodes associated to these elements
            fs_node[:] = False
            fs_node[md.mesh.elements[fs_flag, :] - 1] = True
            ho_node[:] = False
            ho_node[md.mesh.elements[ho_flag, :] - 1] = True
            hofs_node[:] = False
            hofs_node[md.mesh.elements[hofs_flag, :] - 1] = True

        ## Couple FS and SSA
        elif any(fs_flag) and any(ssa_flag):

            ### Find border nodes
            ssafs_node[ssa_node & fs_node] = True

            ### FS elements in contact with this layer become SSAFS elements
            matrix_elements = ssafs_node[md.mesh.elements - 1]
            common_elements = np.sum(matrix_elements, axis=1) != 0
            common_elements[ssa_flag] = False
            fs_flag[common_elements] = False
            ssafs_flag[common_elements] = True
            fs_node = np.zeros(md.mesh.numberofvertices, dtype = bool)
            fs_node[md.mesh.elements[fs_flag].ravel() - 1] = True

            ### Rule out elements that don't touch the 2 boundaries
            pos = ssafs_flag.nonzero()[0]
            elist = (
                np.sum(ssa_node[md.mesh.elements[pos, :] - 1], axis=1).astype(int)
                - np.sum(fs_node[md.mesh.elements[pos, :] - 1], axis=1).astype(int)
            )

            ### Update flags based on elist values
            pos1 = pos[elist == 1]
            ssa_flag[pos1] = True
            ssafs_flag[pos1] = False

            pos2 = pos[elist == -1]
            fs_flag[pos2] = True
            ssafs_flag[pos2] = False

            ### Recompute nodes associated to these elements
            ssa_node[:] = False
            ssa_node[md.mesh.elements[ssa_flag, :] - 1] = True
            fs_node[:] = False
            fs_node[md.mesh.elements[fs_flag, :] - 1] = True
            ssafs_node[:] = False
            ssafs_node[md.mesh.elements[ssafs_flag, :] - 1] = True
    
        elif any(fs_flag) and any(sia_flag):
            raise TypeError('pyissm.tools.param.set_flow_equation: Type of coupling not supported yet.')

    # Create SSAHOApproximation where needed
    md.flowequation.element_equation = np.zeros(md.mesh.numberofelements, int)
    md.flowequation.element_equation[none_flag] = 0
    md.flowequation.element_equation[sia_flag] = 1
    md.flowequation.element_equation[ssa_flag] = 2
    md.flowequation.element_equation[l1l2_flag] = 3
    md.flowequation.element_equation[molho_flag] = 4
    md.flowequation.element_equation[ho_flag] = 5
    md.flowequation.element_equation[fs_flag] = 6
    md.flowequation.element_equation[ssaho_flag] = 7
    md.flowequation.element_equation[ssafs_flag] = 8
    md.flowequation.element_equation[hofs_flag] = 9

    # Define border
    md.flowequation.borderHO = ho_node.astype(int)
    md.flowequation.borderSSA = ssa_node.astype(int)
    md.flowequation.borderFS = fs_node.astype(int)

    # Create vertices_type
    md.flowequation.vertex_equation = np.zeros(md.mesh.numberofvertices, int)
    md.flowequation.vertex_equation[ssa_node] = 2
    md.flowequation.vertex_equation[l1l2_node] = 3
    md.flowequation.vertex_equation[molho_node] = 4
    md.flowequation.vertex_equation[ho_node] = 5
    md.flowequation.vertex_equation[fs_node] = 6
    ## Do SIA last so spcs are setup correctly (SIA has priority )
    md.flowequation.vertex_equation[sia_node] = 1
    if any(fs_flag):
        if not (any(ho_flag) or any(ssa_flag)):
            md.flowequation.vertex_equation[~fs_node] = 0
    md.flowequation.vertex_equation[ssaho_node] = 7
    md.flowequation.vertex_equation[hofs_node] = 8
    md.flowequation.vertex_equation[ssafs_node] = 9
    
    # Define solution types
    md.flowequation.isSIA = int(any(md.flowequation.element_equation == 1))
    md.flowequation.isSSA = int(any(md.flowequation.element_equation == 2))
    md.flowequation.isL1L2= int(any(md.flowequation.element_equation == 3))
    md.flowequation.isMOLHO= int(any(md.flowequation.element_equation == 4))
    md.flowequation.isHO = int(any(md.flowequation.element_equation == 5))
    md.flowequation.isFS = int(any(md.flowequation.element_equation == 6))

    return md