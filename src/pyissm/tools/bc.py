"""
Tools for setting boundary conditions on ISSM models.
"""

import warnings
import numpy as np
import os
from .. import utils

def get_ice_front_nodes(md, ice_front_exp):
    """Identify nodes on the ice front from the provided contour file."""

    # Error checks
    if not os.path.exists(ice_front_exp):
        raise IOError(f"SetIceSheetBC: Ice front file '{ice_front_exp}' not found. Please provide a valid file path.")
    
    # Identify nodes inside the contour
    ## If Python wrappers are installed, use them to identify nodes inside the contour
    if utils.wrappers.check_wrappers_installed():
        node_inside_ice_front = utils.wrappers.ContourToMesh(md.mesh.elements,
                                                                    md.mesh.x,
                                                                    md.mesh.y,
                                                                    ice_front_exp,
                                                                    interp_type = 'node',
                                                                    edge_value = 2)
        
        # Identify nodes on the mesh boundary AND inside the ice front contour
        node_on_ice_front = np.logical_and(md.mesh.vertexonboundary, node_inside_ice_front.reshape(-1))

    else:
        ## If Python wrappers are not installed, issue a warning and return all False (No ice front)
        warnings.warn('pyissm.tools.bc.get_ice_front_nodes: Python wrappers not installed. Cannot identify ice front nodes.\n'
                      'Returning all False array (no ice front).')
        node_on_ice_front = np.zeros((md.mesh.numberofvertices), bool)
    
    return node_on_ice_front

def set_neumann_bc(md, node_on_ice_front):
    """Set Neumann boundary conditions for the ice front."""

    # Find indices of nodes on the ice front
    ice_front_node_indices = np.nonzero(node_on_ice_front)[0]

    # Mark ice front position in the level set
    md.mask.ice_levelset[ice_front_node_indices] = 0

def set_sb_dirichlet_bc(md,
                     node_on_ice_front):
    """Set Dirichlet boundary conditions on boundary (excluding ice front)."""

    # Set empty spc arrays
    md.stressbalance.spcvx = np.nan * np.ones(md.mesh.numberofvertices)
    md.stressbalance.spcvy = np.nan * np.ones(md.mesh.numberofvertices)
    md.stressbalance.spcvz = np.nan * np.ones(md.mesh.numberofvertices)
    md.stressbalance.referential = np.nan * np.ones((md.mesh.numberofvertices, 6))
    md.stressbalance.loadingforce = 0 * np.ones((md.mesh.numberofvertices, 3))

    # Identify mesh type and number of nodes on ice front segments
    if md.mesh.element_type() == 'Penta':
        num_nodes_on_front = 4
    elif md.mesh.element_type() == 'Tria':
        num_nodes_on_front = 2
    else:
        raise NameError(f"set_dirichlet_bc: Unsupported mesh type '{md.mesh.element_type()}'.")
    
    # Find segnments not completely on the ice front
    if np.any(md.mask.ice_levelset <= 0):
        ## Get ice_levelset values for each boundary segment (exclude last column)
        ice_levelset_values = md.mask.ice_levelset[md.mesh.segments[:, :-1] - 1]
        
        ## Convert to binary
        ice_front_segments = 1 - ice_levelset_values
        
        ## Identify segments not completely on the ice front -- A segment is not on the ice front if the sum of its ice_front_segments values is not equal to num_nodes_on_front
        segments_not_on_front = np.nonzero(np.sum(ice_front_segments, axis=1) != num_nodes_on_front)[0]

        ## Get all nodes not on ice-front
        boundary_node_indices = md.mesh.segments[segments_not_on_front, :-1] - 1
    else:
        ## If no nodes are ice, select all boundary nodes
        boundary_node_indices = np.nonzero(md.mesh.vertexonboundary)[0]
    
    # Apply Dirichlet conditions
    ## If observed velocities are available, use them as Dirichlet values
    if (isinstance(md.inversion.vx_obs, np.ndarray) and
        md.inversion.vx_obs.shape[0] == md.mesh.numberofvertices and
        isinstance(md.inversion.vy_obs, np.ndarray) and
        md.inversion.vy_obs.shape[0] == md.mesh.numberofvertices):
        
        # Reshape to rank-2 if necessary to match spc arrays
        if md.inversion.vx_obs.ndim == 1:
            md.inversion.vx_obs = md.inversion.vx_obs.reshape(-1, )
        if md.inversion.vy_obs.ndim == 1:
            md.inversion.vy_obs = md.inversion.vy_obs.reshape(-1, )

        warnings.warn('pyissm.tools.bc.set_sb_dirichlet_bc: Using observed velocities for stressbalance model boundary conditions.')
        md.stressbalance.spcvx[boundary_node_indices] = md.inversion.vx_obs[boundary_node_indices]
        md.stressbalance.spcvy[boundary_node_indices] = md.inversion.vy_obs[boundary_node_indices]
    
    ## If observed velocities are not available, set Dirichlet values to zero
    else:
        warnings.warn('pyissm.tools.bc.set_sb_dirichlet_bc: No observed velocities found. Setting stressbalance model boundary conditions as 0.')
        md.stressbalance.spcvx[boundary_node_indices] = 0
        md.stressbalance.spcvy[boundary_node_indices] = 0
        md.stressbalance.spcvz[boundary_node_indices] = 0

def set_ice_shelf_bc(md,
                     ice_front_exp = None):
    """
    Set ice shelf boundary conditions for the ISSM model.

    Parameters:
    md : Model
        The ISSM model object.
    ice_front_exp : str, optional
        Path to the ice front contour file. If None, no ice front is assumed.

    Returns:
    node_on_ice_front : np.ndarray
        Boolean array indicating nodes on the ice front.
    """

    # Identify ice front 
    if not ice_front_exp is None:
        ## If an ice front file is provided, identify nodes on the ice front
        node_on_ice_front = get_ice_front_nodes(md, ice_front_exp)
    else: 
        ## If no ice front file is provided, assume no ice front
        node_on_ice_front = np.zeros((md.mesh.numberofvertices), bool)

    # Set neumann BC on ice front
    set_neumann_bc(md, node_on_ice_front)

    # Set dirichlet BC on boundary (excluding ice front)
    set_sb_dirichlet_bc(md, node_on_ice_front)

    # Define other boundary conditions
    ## Initialize smb and basalforcings
    md.smb.initialise(md)
    md.basalforcings.initialise(md)

    ## Define balancethickness BCs
    if np.all(np.isnan(md.balancethickness.thickening_rate)):
        md.balancethickness.thickening_rate = np.zeros((md.mesh.numberofvertices))
        warnings.warn('pyissm.tools.bc.set_ice_shelf_bc: no balancethickness.thickening_rate specified -- values set as 0.')
    md.masstransport.spcthickness = np.nan * np.ones((md.mesh.numberofvertices))
    md.balancethickness.spcthickness = np.nan * np.ones((md.mesh.numberofvertices))
    md.damage.spcdamage = np.nan * np.ones((md.mesh.numberofvertices))

    ## Define thermal BCs
    if (isinstance(md.initialization.temperature, np.ndarray) and 
        md.initialization.temperature.shape[0] == md.mesh.numberofvertices):
        
        md.thermal.spctemperature = np.nan * np.ones((md.mesh.numberofvertices))
        
        if hasattr(md.mesh, 'vertexonsurface'):
            vertex_on_surface = np.nonzero(md.mesh.vertexonsurface)[0]
            md.thermal.spctemperature[vertex_on_surface] = md.initialization.temperature[vertex_on_surface]  #impose observed temperature on surface
        
        if (not isinstance(md.basalforcings.geothermalflux, np.ndarray) or 
            md.basalforcings.geothermalflux.shape[0] != md.mesh.numberofvertices):
            md.basalforcings.geothermalflux = np.zeros((md.mesh.numberofvertices))

    else:
        warnings.warn('pyissm.tools.bc.set_ice_shelf_bc: No observed temperature found. No thermal boundary conditions created.')
    
    return md



