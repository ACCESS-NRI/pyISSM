"""
Tools for setting boundary conditions on ISSM models.
"""

import warnings
import numpy as np
import os
from pyissm import tools

def get_ice_front_nodes(md, ice_front_exp):
    """Identify nodes on the ice front from the provided contour file."""

    # Error checks
    if not os.path.exists(ice_front_exp):
        raise IOError(f"SetIceSheetBC: Ice front file '{ice_front_exp}' not found. Please provide a valid file path.")
    
    # Identify nodes inside the contour
    ## If Python wrappers are installed, use them to identify nodes inside the contour
    if tools.wrappers.check_wrappers_installed():
        node_inside_ice_front = tools.wrappers.ContourToMesh(md.mesh.elements,
                                                                    md.mesh.x,
                                                                    md.mesh.y,
                                                                    ice_front_exp,
                                                                    interp_type = 'node',
                                                                    edge_value = 2)
        
        # Identify nodes on the mesh boundary AND inside the ice front contour
        node_on_ice_front = np.logical_and(md.mesh.vertexonboundary, node_inside_ice_front.reshape(-1))

    else:
        ## If Python wrappers are not installed, issue a warning and return all False (No ice front)
        warnings.warn('pyissm.model.bc.get_ice_front_nodes: Python wrappers not installed. Cannot identify ice front nodes.\n'
                      'Returning all False array (no ice front).')
        node_on_ice_front = np.zeros((md.mesh.numberofvertices), bool)
    
    return node_on_ice_front

def set_neumann_bc(md, node_on_ice_front):
    """Set Neumann boundary conditions for the ice front."""

    # Find indices of nodes on the ice front
    ice_front_nodes = np.nonzero(node_on_ice_front)[0]

    # Mark ice front position in the level set
    md.mask.ice_levelset[ice_front_nodes] = 0

def set_sb_dirichlet_bc(md):
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

        warnings.warn('pyissm.model.bc.set_sb_dirichlet_bc: Using observed velocities for vx and vy stressbalance model boundary conditions. vz boundary conditions are set to 0.')
        md.stressbalance.spcvx[boundary_node_indices] = md.inversion.vx_obs[boundary_node_indices]
        md.stressbalance.spcvy[boundary_node_indices] = md.inversion.vy_obs[boundary_node_indices]
        md.stressbalance.spcvz[boundary_node_indices] = 0
    
    ## If observed velocities are not available, set Dirichlet values to zero
    else:
        warnings.warn('pyissm.model.bc.set_sb_dirichlet_bc: No observed velocities found. Setting stressbalance model boundary conditions as 0.')
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
    md : Model
        The ISSM model object with updated boundary conditions.
    """

    # Identify ice front 
    if not ice_front_exp is None:
        ## If an ice front file is provided, identify nodes on the ice front
        node_on_ice_front = get_ice_front_nodes(md, ice_front_exp)
    else: 
        ## If no ice front file is provided, assume no ice front
        warnings.warn('pyissm.model.bc.set_ice_shelf_bc: No ice front file provided. Assuming no ice front.')
        node_on_ice_front = np.zeros((md.mesh.numberofvertices), bool)

    # Set neumann BC on ice front
    set_neumann_bc(md, node_on_ice_front)

    # Set dirichlet BC on boundary (excluding ice front)
    set_sb_dirichlet_bc(md)

    # Define other boundary conditions
    ## Initialize smb and basalforcings
    md.smb.initialise(md)
    md.basalforcings.initialise(md)

    ## Define balancethickness BCs
    if np.all(np.isnan(md.balancethickness.thickening_rate)):
        md.balancethickness.thickening_rate = np.zeros((md.mesh.numberofvertices))
        warnings.warn('pyissm.model.bc.set_ice_shelf_bc: no balancethickness.thickening_rate specified -- values set as 0.')
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
        warnings.warn('pyissm.model.bc.set_ice_shelf_bc: No observed temperature found. No thermal boundary conditions created.')
    
    return md

def set_molho_bc(md):
    """
    Set MOLHO boundary conditions for ISSM model.

    Parameters:
    md : Model
        The ISSM model object.
    Returns:
    md : Model
        The ISSM model object with updated boundary conditions.
    
    """

    md.stressbalance.spcvx_base = md.stressbalance.spcvx
    md.stressbalance.spcvy_base = md.stressbalance.spcvy
    md.stressbalance.spcvx_shear = np.full(md.mesh.numberofvertices, np.nan)
    md.stressbalance.spcvy_shear = np.full(md.mesh.numberofvertices, np.nan)

    return md

def set_ice_sheet_bc(md):
    """
    Set ice sheet boundary conditions for the ISSM model.

    Parameters:
    md : Model
        The ISSM model object.

    Returns:
    md : Model
        The ISSM model object with updated boundary conditions.
    """

    # Set dirichlet BC on boundary (excluding ice front)
    set_sb_dirichlet_bc(md)

    # Initialise surface and basal forcings
    md.smb.initialise(md)
    md.basalforcings.initialise(md)

    # Initialise ocean forcings and sea level
    md.dsl.initialise(md)

    ## Define balancethickness BCs
    if np.all(np.isnan(md.balancethickness.thickening_rate)):
        md.balancethickness.thickening_rate = np.zeros((md.mesh.numberofvertices))
        warnings.warn('pyissm.model.bc.set_ice_sheet_bc: no balancethickness.thickening_rate specified -- values set as 0.')
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
            md.basalforcings.geothermalflux = 50.0 * pow(10, -3) * np.ones((md.mesh.numberofvertices))  # 50 mW/m^2

    else:
        warnings.warn('pyissm.model.bc.set_ice_sheet_bc: no observed temperature found -- no thermal boundary conditions created.')

    return md

def set_marine_ice_sheet_bc(md,
                            ice_front_exp = None):
    """
    Set marine ice sheet boundary conditions for the ISSM model.

    Parameters:
    md : Model
        The ISSM model object.
    ice_front_exp : str, optional
        Path to the ice front contour file. If None, the ice front is guessed from ocean levelset.

    Returns:
    md : Model
        The ISSM model object with updated boundary conditions.
    """

    # Identify ice front 
    if not ice_front_exp is None:
        ## If an ice front file is provided, identify nodes on the ice front
        node_on_ice_front = get_ice_front_nodes(md, ice_front_exp)
        vertex_on_ice_front = md.mesh.vertexonboundary & np.asarray(node_on_ice_front, dtype=bool).ravel()
    else: 
        # Guess where the ice front is
        vertex_on_floating_ice = np.zeros((md.mesh.numberofvertices))
        pos = np.nonzero(np.sum(md.mask.ocean_levelset[md.mesh.elements - 1] < 0.0, axis=1) > 0.0)[0]
        vertex_on_floating_ice[md.mesh.elements[pos].astype(int) - 1] = 1.
        vertex_on_ice_front = np.logical_and(np.reshape(md.mesh.vertexonboundary, (-1, )), vertex_on_floating_ice > 0.0)

    # Set neumann BC on ice front
    pos = np.nonzero(np.logical_and(md.mesh.vertexonboundary, np.logical_not(vertex_on_ice_front)))[0]
    if not np.size(pos):
        warnings.warn('pyissm.model.bc.set_marine_ice_sheet_bc: Ice front all around the glacier, no dirichlet found. Dirichlet must be added manually.')

    md.stressbalance.spcvx = np.nan * np.ones(md.mesh.numberofvertices)
    md.stressbalance.spcvy = np.nan * np.ones(md.mesh.numberofvertices)
    md.stressbalance.spcvz = np.nan * np.ones(md.mesh.numberofvertices)
    md.stressbalance.referential = np.nan * np.ones((md.mesh.numberofvertices, 6))
    md.stressbalance.loadingforce = 0 * np.ones((md.mesh.numberofvertices, 3))

    # Position of ice front
    pos = np.nonzero(vertex_on_ice_front)[0]
    md.mask.ice_levelset[pos] = 0
    
    # Set Dirichlet BCs
    set_sb_dirichlet_bc(md)

    # Define other boundary conditions
    ## Initialize smb and basalforcings
    md.smb.initialise(md)
    md.basalforcings.initialise(md)

    ## Define hydrology BCs
    md.hydrology.spcwatercolumn = np.zeros((md.mesh.numberofvertices, 2))
    pos = np.nonzero(md.mesh.vertexonboundary)[0]
    md.hydrology.spcwatercolumn[pos, 0] = 1

    ## Define balancethickness BCs
    if np.all(np.isnan(md.balancethickness.thickening_rate)):
        md.balancethickness.thickening_rate = np.zeros((md.mesh.numberofvertices))
        warnings.warn('pyissm.model.bc.set_ice_shelf_bc: no balancethickness.thickening_rate specified -- values set as 0.')
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
            md.basalforcings.geothermalflux[md.mask.ocean_levelset > 0.0] = 50e-3  # 50 mW/m²

    else:
        warnings.warn('pyissm.model.bc.set_ice_shelf_bc: No observed temperature found. No thermal boundary conditions created.')
    
    return md