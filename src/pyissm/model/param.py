"""
Functions for parameterising ISSM models.
"""

import numpy as np
import os
from pathlib import Path
from datetime import datetime
from pyissm import model, tools
import copy

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
            raise ImportError("pyissm.model.param.set_mask: Wrappers are not installed. Cannot use ice_domain option.")
    
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
        raise ValueError("pyissm.model.param.set_flow_equation: Coupling must be 'tiling' or 'penalties'.")
    
    if fill is not None:
        if fill.lower() not in ['sia', 'ssa', 'ho', 'l1l2', 'molho', 'fs']:
            raise ValueError("pyissm.model.param.set_flow_equation: Fill must be one of: None, 'SIA', 'SSA', 'HO', 'L1L2', 'MOLHO', or 'FS'.")
    
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
        raise ValueError("pyissm.model.param.set_flow_equation: One or more elements have not been assigned a flow equation.")    
    
    ## Check no elements have been assigned multiple flow equations
    if np.any(np.sum(flag, axis=0) > 1):
        raise ValueError("pyissm.model.param.set_flow_equation: One or more elements have been assigned multiple flow equations.")

    ## Check that L1L2, MOLHO, and FS are not coupled to any other model for now
    if np.any(l1l2_flag) and np.any(sia_flag | ssa_flag | ho_flag | fs_flag | molho_flag):
        raise ValueError('pyissm.model.param.set_flow_equation: L1L2 cannot be coupled to other flow equations.')
    if np.any(molho_flag) and np.any(sia_flag | ssa_flag | ho_flag | fs_flag | l1l2_flag):
        raise ValueError('pyissm.model.param.set_flow_equation: MOLHO cannot be coupled to other flow equations.')
    if np.any(fs_flag) & np.any(sia_flag | l1l2_flag | molho_flag):
        raise ValueError('pyissm.model.param.set_flow_equation: FS cannot be coupled to SAI, L1L2, or MOLHO.')
        
    ##  Check HO or FS are not used for a 2D mesh
    if md.mesh.domain_type() == '2Dhorizontal':
        if np.any(ho_flag):
            raise ValueError('pyissm.model.param.set_flow_equation: HO cannot be used for a 2D mesh. Extrude it first.')
        if np.any(fs_flag):
            raise ValueError('pyissm.model.param.set_flow_equation: FS cannot be used for a 2D mesh. Extrude it first.')
    
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

            ### Recompute nodes associated to these elements
            ssa_node[:] = False
            ssa_node[md.mesh.elements[ssa_flag].ravel() - 1] = True

            ### Rule out elements that don't touch the 2 boundaries
            pos = np.where(ssaho_flag)[0]
            elist = (
                np.sum(ho_node[md.mesh.elements[pos, :] - 1], axis=1).astype(int)
                - np.sum(ssa_node[md.mesh.elements[pos, :] - 1], axis=1).astype(int)
            )
            
            pos1 = pos[elist == 1]
            ho_flag[pos1] = True
            ssaho_flag[pos1] = False

            pos2 = pos[elist == -1]
            ssa_flag[pos2] = True
            ssaho_flag[pos2] = False

            ### Recompute nodes associated to these elements
            ssa_node[:] = False
            ssa_node[md.mesh.elements[ssa_flag].ravel() - 1] = True
            ho_node[:] = False
            ho_node[md.mesh.elements[ho_flag].ravel() - 1] = True
            ssaho_node[:] = False
            ssaho_node[md.mesh.elements[ssaho_flag].ravel() - 1] = True

        ## Couple HO and FS
        elif any(ho_flag) and any(fs_flag):
            hofs_node[ho_node & fs_node] = True

            ### FS elements in contact with this layer become HOFS elements
            matrix_elements = hofs_node[md.mesh.elements - 1]
            common_elements = np.sum(matrix_elements, axis=1) != 0
            common_elements[ho_flag] = False

            fs_flag[common_elements] = False
            hofs_flag[common_elements] = True

            ### Recompute nodes associated to these elements
            fs_node[:] = False
            fs_node[md.mesh.elements[fs_flag].ravel() - 1] = True

            ### Rule out elements that don't touch the 2 boundaries
            pos = np.where(hofs_flag)[0]
            elist = (
                np.sum(fs_node[md.mesh.elements[pos, :] - 1], axis=1).astype(int)
                - np.sum(ho_node[md.mesh.elements[pos, :] - 1], axis=1).astype(int)
            )
            
            pos1 = pos[elist == 1]
            fs_flag[pos1] = True
            hofs_flag[pos1] = False

            pos2 = pos[elist == -1]
            ho_flag[pos2] = True
            hofs_flag[pos2] = False

            ## Recompute nodes again
            fs_node[:] = False
            fs_node[md.mesh.elements[fs_flag].ravel() - 1] = True
            ho_node[:] = False
            ho_node[md.mesh.elements[ho_flag].ravel() - 1] = True
            hofs_node[:] = False
            hofs_node[md.mesh.elements[hofs_flag].ravel() - 1] = True

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
            raise TypeError('pyissm.model.param.set_flow_equation: Type of coupling not supported yet.')

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

def parameterize(md, parameter_file):
    """
    Parameterize an ISSM model from a Python parameter file.

    Parameters
    ----------
    md : object
        ISSM model instance to populate. The parameter file is expected to
        mutate this object in-place and may refer to it as ``md``.
    parameter_file : str or pathlib.Path
        Path to a Python file containing parameter-setting code. This file is
        executed with an empty global namespace and a local namespace where
        ``md`` is pre-defined.

    Returns
    -------
    object
        The modified model instance (same object passed in).

    Raises
    ------
    FileNotFoundError
        If the provided parameter_file does not exist.
    Exception
        Any exception raised while executing the parameter file is propagated.

    Notes
    -----
    The parameter file must be a valid Python script that assigns values or
    calls functions that modify the provided ``md`` object. For security,
    be cautious when executing untrusted parameter files since they are run
    with ``exec()`` and may perform arbitrary operations.

    Examples
    --------
    >>> parameterize(md, "parameters.py")
    """

    # Get path
    path = Path(parameter_file)

    # Error checks
    if not path.exists():
        raise FileNotFoundError(f"pyissm.model.param.parameterize: Parameter file {parameter_file} not found.")
    
    # Execute the parameter file
    # NOTE: Local execution namespace (parameter files expect "md" to exist)
    local_env = {"md": md}

    with path.open("r") as f:
        code = compile(f.read(), str(path), "exec")
        exec(code, {}, local_env)

    # Set name if missing
    if not getattr(md.miscellaneous, "name", None):
        md.miscellaneous.name = path.stem

    # Set timestamp note
    timestamp = datetime.now().strftime("%c")
    md.miscellaneous.notes = (
        f"Model created using parameter file '{parameter_file}' on {timestamp}."
    )

    return md

def contour_envelope(mh, flags = None):
    """
    Build a set of segments enveloping a contour.
    This function computes segments that form the boundary envelope of a contour
    within a given mesh. It identifies elements on the domain boundary and extracts
    the segments that define the contour envelope.
    Parameters
    ----------
    mh : ISSM model mesh object containing vertex and element information.

    flags : {str, int, float, bool}, optional
        Currently not supported. Reserved for future use. If provided, raises
        NotImplementedError. Default is None.

    Returns
    -------
    segments : ndarray
        Array of shape (n_segments, 3) containing segment information where:
        - Column 0: First node index of the segment
        - Column 1: Second node index of the segment
        - Column 2: Associated element index (1-indexed)

    Raises
    ------
    NotImplementedError
        If `flags` argument is provided (not yet supported).
    ImportError
        If ISSM wrappers are not installed.

    Notes
    -----
    This function requires ISSM wrappers to be installed. It computes node and
    element connectivity tables and identifies boundary elements that separate
    interior from exterior regions.

    Examples
    --------
    >>> segments = contour_envelope(mh)
    """

    # Some checks
    if flags is not None:

        ## NOTE: This function is taken from the original ISSM Python code: $ISSM_DIR/src/m/mesh/contour_envelope.py with only minor modifications for pyISSM integration.
        ## It differs with MATLAB when Flags is provided as an array or file. For now, raise error if flags argument is used (not yet supported)
        if flags is not None:
            raise NotImplementedError("pyissm.model.param.contour_envelope: The `flags` argument is not yet supported. Contact ACCESS-NRI for assistance")

        # Original checks below
        if isinstance(flags, str):
            file = flags
            if not os.path.exists(file):
                raise IOError(f"pyissm.model.param.contour_envelope: file {file} not found")
            isfile = 1
        elif isinstance(flags, (bool, int, float)):
            #do nothing for now
            isfile = 0
        else:
            raise TypeError("pyissm.model.param.contour_envelope: second argument should be a file or an elements flag")
        
    # Check that wrappers are installed
    if not tools.wrappers.check_wrappers_installed():
        raise ImportError("pyissm.model.param.contour_envelope: This function requires ISSM wrappers to be installed.")

    # Now, build the connectivity tables for this mesh
    # Computing connectivity
    if np.size(mh.vertexconnectivity, axis=0) != mh.numberofvertices and np.size(mh.vertexconnectivity, axis=0) != mh.numberofvertices2d:
        mh.vertexconnectivity = tools.wrappers.NodeConnectivity(mh.elements, mh.numberofvertices)
    if np.size(mh.elementconnectivity, axis=0) != mh.numberofelements and np.size(mh.elementconnectivity, axis=0) != mh.numberofelements2d:
        mh.elementconnectivity = tools.wrappers.ElementConnectivity(mh.elements, mh.vertexconnectivity)

    #get nodes inside profile
    elementconnectivity = copy.deepcopy(mh.elementconnectivity)
    if mh.dimension() == 2:
        elements = copy.deepcopy(mh.elements)
        x = copy.deepcopy(mh.x)
        y = copy.deepcopy(mh.y)
        numberofvertices = copy.deepcopy(mh.numberofvertices)
        numberofelements = copy.deepcopy(mh.numberofelements)
    else:
        elements = copy.deepcopy(mh.elements2d)
        x = copy.deepcopy(mh.x2d)
        y = copy.deepcopy(mh.y2d)
        numberofvertices = copy.deepcopy(mh.numberofvertices2d)
        numberofelements = copy.deepcopy(mh.numberofelements2d)

    if flags is not None:
        if isfile:
            # Get flag list of elements and nodes inside the contour
            nodein = tools.wrappers.ContourToMesh(elements, x, y, file, 'node', 1)
            elemin = (np.sum(nodein(elements), axis=1) == np.size(elements, axis=1))
            # Modify element connectivity
            elemout = np.nonzero(np.logical_not(elemin))[0]
            elementconnectivity[elemout, :] = 0
            elementconnectivity[np.nonzero(np.isin(elementconnectivity, elemout + 1))] = 0
        else:
            # Get flag list of elements and nodes inside the contour
            nodein = np.zeros(numberofvertices)
            elemin = np.zeros(numberofelements)

            pos = np.nonzero(flags)
            elemin[pos] = 1
            nodein[elements[pos, :] - 1] = 1

            # Modify element connectivity
            elemout = np.nonzero(np.logical_not(elemin))[0]
            elementconnectivity[elemout, :] = 0
            elementconnectivity[np.nonzero(np.isin(elementconnectivity, elemout + 1))] = 0

    # Find element on boundary
    # First: find elements on the boundary of the domain
    flag = copy.deepcopy(elementconnectivity)
    if flags is not None:
        flag[np.nonzero(flag)] = elemin[flag[np.nonzero(flag)]]
    elementonboundary = np.logical_and(np.prod(flag, axis=1) == 0, np.sum(flag, axis=1) > 0)

    # Find segments on boundary
    pos = np.nonzero(elementonboundary)[0]
    num_segments = np.size(pos)
    segments = np.zeros((num_segments * 3, 3), int)
    count = 0

    for el1 in pos:
        els2 = elementconnectivity[el1, np.nonzero(elementconnectivity[el1, :])[0]] - 1
        if np.size(els2) > 1:
            flag = np.intersect1d(np.intersect1d(elements[els2[0], :], elements[els2[1], :]), elements[el1, :])
            nods1 = elements[el1, :]
            nods1 = np.delete(nods1, np.nonzero(nods1 == flag))
            segments[count, :] = [nods1[0], nods1[1], el1 + 1]

            ord1 = np.nonzero(nods1[0] == elements[el1, :])[0][0]
            ord2 = np.nonzero(nods1[1] == elements[el1, :])[0][0]

    #swap segment nodes if necessary
            if ((ord1 == 0 and ord2 == 1) or (ord1 == 1 and ord2 == 2) or (ord1 == 2 and ord2 == 0)):
                temp = segments[count, 0]
                segments[count, 0] = segments[count, 1]
                segments[count, 1] = temp
            segments[count, 0:2] = np.flipud(segments[count, 0:2])
            count += 1
        else:
            nods1 = elements[el1, :]
            flag = np.setdiff1d(nods1, elements[els2, :])
            for j in range(0, 3):
                nods = np.delete(nods1, j)
                if np.any(np.isin(flag, nods)):
                    segments[count, :] = [nods[0], nods[1], el1 + 1]
                    ord1 = np.nonzero(nods[0] == elements[el1, :])[0][0]
                    ord2 = np.nonzero(nods[1] == elements[el1, :])[0][0]
                    if ((ord1 == 0 and ord2 == 1) or (ord1 == 1 and ord2 == 2) or (ord1 == 2 and ord2 == 0)):
                        temp = segments[count, 0]
                        segments[count, 0] = segments[count, 1]
                        segments[count, 1] = temp
                    segments[count, 0:2] = np.flipud(segments[count, 0:2])
                    count += 1
    segments = segments[0:count, :]

    return segments

def kill_icebergs(md):
    """
    Remove isolated floating ice patches (icebergs) by setting their ice_levelset values to +1.

    Parameters
    ----------
    md : ISSM model object
        The model containing the mesh and mask information.

    Returns
    -------
    ndarray
        Modified ice_levelset array with isolated icebergs set to +1.
        If no icebergs are found, returns a copy of the original ice_levelset.

    Notes
    -----
    The algorithm works in three main steps:

    1. Mark elements without ice as done
    2. Initialize mask from grounded ice elements
    3. Use flood-fill algorithm to propagate connectivity from grounded ice through 
       floating ice patches
    4. Set vertices that remain unconnected and have ice to +1 (removing icebergs)

    Isolated floating ice patches are identified as ice vertices that cannot be 
    reached through the flood-fill algorithm starting from grounded ice regions.
    """

    elements = md.mesh.elements - 1
    ice_ls   = md.mask.ice_levelset
    ocean_ls = md.mask.ocean_levelset

    nverts = md.mesh.numberofvertices
    nelems = md.mesh.numberofelements

    mask = np.zeros(nverts, dtype=np.int8)
    element_flag = np.zeros(nelems, dtype=np.int8)

    print("Looking for isolated patches of floating ice (icebergs)")

    # Identify and mark elements with ice
    isice = np.min(ice_ls[elements], axis=1) < 0
    element_flag[~isice] = 1

    # Identify grounded ice elements
    isgrounded = np.sum(ocean_ls[elements] > 0, axis=1) > 2
    grounded_idx = np.where(isgrounded)[0]

    if grounded_idx.size > 0:
        mask[elements[grounded_idx].ravel()] = 1

    # Ice-free vertices should stay 0
    mask[ice_ls >= 0] = 0

    # Flood-fill from grounded ice through floating ice
    iteration = 1
    more = True

    while more:
        print(f"   -- iteration {iteration}")
        more = False

        remaining_elems = np.where(element_flag == 0)[0]

        for i in remaining_elems:
            idx = elements[i]
            # If at least 2 vertices are already connected → activate element
            if np.sum(mask[idx] > 0) > 1:
                element_flag[i] = 1
                mask[idx] = 1
                more = True

        iteration += 1

    # Set isolated floating ice to +1
    pos = np.where((mask == 0) & (ice_ls < 0))[0]

    if pos.size > 0:
        print(f"REMOVING {pos.size} vertex{'es' if pos.size > 1 else ''} on icebergs")
        new_ice_ls = ice_ls.copy()
        new_ice_ls[pos] = +1
        return new_ice_ls

    print("No iceberg found!")
    return ice_ls.copy()