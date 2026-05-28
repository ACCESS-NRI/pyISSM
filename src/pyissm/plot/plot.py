"""
Functions to visualize ISSM Models
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from pyissm import model, tools
from pyissm.plot.colormaps import (
    seacolor, landcolor, ibcao, demmap, get_colormap, truncate_colormap
)
import warnings

## ------------------------------------------------------------------------------------
## MESH PLOTTING
## ------------------------------------------------------------------------------------
def plot_mesh2d(md,
                ax = None,
                color = 'k',
                linewidth = 0.1,
                xlabel = 'X (m)',
                ylabel = 'Y (m)',
                figsize = (6.4, 4.8),
                constrained_layout = True,
                show_nodes = False,
                node_kwargs = {},
                **kwargs):
    """
    Plot a 2D triangular mesh using matplotlib.

    Parameters
    ----------
    md : ISSM Model object
        ISSM Model object containing mesh. Must be compatible with process_mesh().
    ax : matplotlib.axes.Axes, optional
        An existing matplotlib axes object to plot on. If None, a new figure and axes are created.
    color : str, optional
        Color of the triangle edges. Default is 'k' (black).
    linewidth : float, optional
        Width of the triangle edges. Default is 0.1.
    xlabel : str, optional
        Label for the x-axis. Default is 'X (m)'.
    ylabel : str, optional
        Label for the y-axis. Default is 'Y (m)'.
    figsize : tuple of float, optional
        Figure size in inches if a new figure is created. Default is (6.4, 4.8).
    constrained_layout : bool, optional
        Whether to use constrained layout when creating a new figure. Default is 'True'.
    show_nodes : bool, optional
        If True, overlay nodes on the mesh using ax.scatter(). Default is 'False'.
    node_kwargs : dict, optional
        Keyword arguments passed to ax.scatter() for node display. Defaults to
        {'marker': '.', 'color': 'k', 's': 5}.
    **kwargs
        Additional keyword arguments passed to ax.triplot().

    Returns
    -------
    matplotlib.figure.Figure or matplotlib.axes.Axes
        If 'ax' is None, returns '(fig, ax)' of the created figure and axes.
        If 'ax' is provided, returns the modified 'ax'.

    Example
    -------
    fig, ax = plot_mesh2d(md)
    fig, ax = plot_mesh2d(md, color = 'blue', linewidth = 0.5)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1 = plot_mesh2d(md, ax = ax1)
    ax2 = plot_mesh2d(md, ax = ax2, show_nodes = True, node_kwargs = {'color': 'red'})
    """

    ## Is an ax passed?
    ax_defined = ax is not None

    ## Define default node kwargs and update node_kwargs with passed arguments
    default_node_kwargs = {'marker': '.',
                         'color': 'k',
                         's': 5,
                         'zorder': 3}
    default_node_kwargs.update(**node_kwargs)


    ## If no ax is defined, create new figure
    ## otherwise, plot on defined ax
    if ax is None:
        fig, ax = plt.subplots(figsize = figsize, constrained_layout = constrained_layout)
    else:
        fig = ax.get_figure()

    ## Process mesh for plotting
    mesh, mesh_x, mesh_y, mesh_elements, is3d = model.mesh.process_mesh(md)

    ## Make plot
    ax.triplot(mesh,
               color = color,
               linewidth = linewidth,
               zorder = 2,
               **kwargs)

    ## Add optional nodes
    if show_nodes:
        ax.scatter(mesh.x, mesh.y, **default_node_kwargs)

    ## Add axis labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if not ax_defined:
        return fig, ax
    else:
        return ax

## ------------------------------------------------------------------------------------
## NODE / ELEMENT TYPE PLOTTING
## ------------------------------------------------------------------------------------
def plot_model_nodes(md,
                     ice_levelset,
                     ocean_levelset,
                     type = 'ice_nodes',
                     ax = None,
                     s = 5,
                     colors = ['r', 'k'],
                     marker = 'o',
                     xlabel = 'X (m)',
                     ylabel = 'Y (m)',
                     figsize = (6.4, 4.8),
                     constrained_layout = True,
                     show_mesh = True,
                     mesh_kwargs = {},
                     show_legend = True,
                     legend_kwargs = {}):

    """
    Plot model nodes by type (ice, ice-front, ocean, floating, grounded) on a 2D mesh.

    This function uses level set fields to classify mesh nodes and visualize them
    in a scatter plot. Optionally overlays the finite element mesh and includes
    a custom legend. Supports plotting in existing Matplotlib axes or creating
    a new figure.

    Parameters
    ----------
    md : ISSM Model object
        ISSM Model object containing mesh. Must be compatible with process_mesh().
    ice_levelset : ndarray
        Array of ice level set values. Negative values indicate ice-covered nodes; zero indicates the ice front.
    ocean_levelset : ndarray
        Array of ocean level set values. Negative values indicate ocean-covered nodes.
    type : {'ice_nodes', 'ice_front_nodes', 'ocean_nodes', 'floating_ice_nodes', 'grounded_ice_nodes'}, optional
        The node type to visualize. Default is 'ice_nodes'.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw on. If None, a new figure and axes are created.
    s : float, optional
        Marker size for scatter plot. Default is 5.
    colors : list of str, optional
        Colors used for binary classification. First color is for "False", second for "True". Default is ['r', 'k'].
    marker : str, optional
        Marker style for nodes. Default is 'o'.
    figsize : tuple of float, optional
        Size of the figure in inches (width, height). Default is (6.4, 4.8).
    constrained_layout : bool, optional
        Whether to use constrained layout in figure. Default is True.
    show_mesh : bool, optional
        Whether to overlay the triangular mesh. Default is True.
    mesh_kwargs : dict, optional
        Additional keyword arguments passed to plot_mesh2d() for customizing mesh appearance.
    show_legend : bool, optional
        Whether to display a legend identifying node types. Default is True.
    legend_kwargs : dict, optional
        Additional keyword arguments passed to ax.legend().

    Returns
    -------
    matplotlib.figure.Figure or matplotlib.axes.Axes
        If 'ax' is None, returns a tuple (fig, ax) with the created figure and axes.
        If 'ax' is provided, returns the updated axes.

    Example
    -------
    fig, ax = plot_model_nodes(md, md.mask.ice_levelset, md.mask.ocean_levelset)
    fig, ax = plot_model_nodes(md, md.mask.ice_levelset, md.mask.ocean_levelset, show_mesh = False)
    fig, ax = plot_model_nodes(md, md.mask.ice_levelset, md.mask.ocean_levelset, type = 'ice_front_nodes', mesh_kwargs = {'linewidth': 0.5, 'color': 'cyan'})
    fig, ax = plot_model_nodes(md, md.mask.ice_levelset, md.mask.ocean_levelset, type = 'floating_ice_nodes', colors = ['blue', 'magenta'])
    """

    ## Set defaults
    ax_defined = ax is not None

    default_mesh_kwargs = {'alpha': 0.5}
    default_mesh_kwargs.update(**mesh_kwargs)

    default_legend_kwargs = {'title': 'Node types',
                           'fontsize': 10}
    default_legend_kwargs.update(**legend_kwargs)

    ## Define binary colormap
    cmap = matplotlib.colors.ListedColormap(colors)

    ## Process model mesh
    mesh, mesh_x, mesh_y, mesh_elements, is3d = model.mesh.process_mesh(md)

    ## Find node types
    node_types = model.mesh.find_node_types(md,
                                            ice_levelset,
                                            ocean_levelset)

    ## Set-up (or retrieve) figure
    if ax is None:
        fig, ax = plt.subplots(figsize = figsize, constrained_layout = constrained_layout)
    else:
        fig = ax.get_figure()

    ## Plot nodes
    ax.scatter(mesh_x,
               mesh_y,
               c = node_types[type],
               s = s,
               cmap = cmap,
               vmin = 0,
               vmax = 1,
               marker = marker)

    ## Add mesh (optional) with specific arguments
    if show_mesh:
        plot_mesh2d(md,
                    ax=ax,
                    **default_mesh_kwargs)

    ## Add legend
    if show_legend:
        labels_dict = {
            'ice_nodes': ['Ice', 'No ice'],
            'ice_front_nodes': ['Ice front', 'No ice front'],
            'ocean_nodes': ['Ocean', 'No ocean'],
            'floating_ice_nodes': ['Floating ice', 'No floating ice'],
            'grounded_ice_nodes': ['Grounded ice', 'No grounded ice']
        }

        legend_elements = [
            matplotlib.lines.Line2D([0], [0], marker = marker, linestyle = 'none', color = 'none', markerfacecolor = colors[1], markeredgecolor='none', label = labels_dict[type][0]),
            matplotlib.lines.Line2D([0], [0], marker = marker, linestyle = 'none', color = 'none', markerfacecolor = colors[0], markeredgecolor='none', label = labels_dict[type][1])
        ]
        ax.legend(handles = legend_elements, **default_legend_kwargs)

    ## Add axis labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ## Return
    if not ax_defined:
        return fig, ax
    else:
        return ax

def plot_model_elements(md,
                        ice_levelset,
                        ocean_levelset,
                        type = 'ice_elements',
                        ax = None,
                        color = 'blue',
                        xlabel = 'X (m)',
                        ylabel = 'Y (m)',
                        figsize = (6.4, 4.8),
                        constrained_layout = True,
                        show_mesh = True,
                        mesh_kwargs = {},
                        show_legend = True,
                        legend_kwargs={}):

    """
    Plot model elements by type (ice, ice-front, ocean, floating, grounded, grounding line) on a 2D mesh.

    This function uses level set fields to classify mesh elements and visualize them
    in a triplot visualisation. Optionally overlays the finite element mesh and includes
    a custom legend. Supports plotting in existing Matplotlib axes or creating
    a new figure.

    Parameters
    ----------
    md : ISSM Model object
        ISSM Model object containing mesh. Must be compatible with process_mesh().
    ice_levelset : ndarray
        Array of ice level set values. Negative values indicate ice-covered nodes; zero indicates the ice front.
    ocean_levelset : ndarray
        Array of ocean level set values. Negative values indicate ocean-covered nodes.
    type : {'ice_elements', 'ice_front_elements', 'ocean_elements', 'floating_ice_elements',
            'grounded_ice_elements', 'grounding_line_elements'}, optional
        The element type to visualize. Default is 'ice_elements'.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw on. If None, a new figure and axes are created.
    color : str, optional
        Color used for elements. Default is "blue".
    figsize : tuple of float, optional
        Size of the figure in inches (width, height). Default is (6.4, 4.8).
    constrained_layout : bool, optional
        Whether to use constrained layout in figure. Default is True.
    show_mesh : bool, optional
        Whether to overlay the triangular mesh. Default is True.
    mesh_kwargs : dict, optional
        Additional keyword arguments passed to plot_mesh2d() for customizing mesh appearance.
    show_legend : bool, optional
        Whether to display a legend identifying node types. Default is True.
    legend_kwargs : dict, optional
        Additional keyword arguments passed to ax.legend().

    Returns
    -------
    matplotlib.figure.Figure or matplotlib.axes.Axes
        If 'ax' is None, returns a tuple (fig, ax) with the created figure and axes.
        If 'ax' is provided, returns the updated axes.

    Example
    -------
    fig, ax = plot_model_elements(md, md.mask.ice_levelset, md.mask.ocean_levelset)
    fig, ax = plot_model_elements(md, md.mask.ice_levelset, md.mask.ocean_levelset, show_mesh = False)
    fig, ax = plot_model_elements(md, md.mask.ice_levelset, md.mask.ocean_levelset, type = 'ice_front_elements', mesh_kwargs = {'linewidth': 0.5, 'color': 'cyan'})
    fig, ax = plot_model_elements(md, md.mask.ice_levelset, md.mask.ocean_levelset, type = 'floating_ice_elements', color = 'red')
    """

    ## Set defaults
    ax_defined = ax is not None

    default_mesh_kwargs = {'alpha': 0.5}
    default_mesh_kwargs.update(**mesh_kwargs)

    default_legend_kwargs = {'title': 'Element type',
                           'fontsize': 10,
                           'loc': 'upper right'}
    default_legend_kwargs.update(**legend_kwargs)

    ## Process model mesh
    mesh, mesh_x, mesh_y, mesh_elements, is3d = model.mesh.process_mesh(md)

    ## Find element types
    element_types = model.mesh.find_element_types(md,
                                                  ice_levelset,
                                                  ocean_levelset)
    # Isolate requested elements
    select_elements = element_types[type]

    # If selected_elements is all False, no elements exist
    if not np.any(select_elements):
        raise ValueError(f'plot_model_elements: No {type} elements exist in the model.')

    ## Get position of elements > 0
    element_pos = np.where(select_elements > 0)

    ## Set-up (or retrieve) figure
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=constrained_layout)
    else:
        fig = ax.get_figure()

    ## Make plot
    # Create colors of required length & static cmap of given colour
    colors = np.ones(np.shape(mesh_elements[element_pos])[0])
    cmap = matplotlib.colors.ListedColormap(color)

    ## Plot elements (shading = 'flat' [default] is required when data are defined on elements)
    ax.tripcolor(mesh_x, mesh_y, mesh_elements[element_pos], facecolors=colors, cmap = cmap, edgecolors = 'none', shading = 'flat')

    ## Add mesh (optional) with specific arguments
    if show_mesh:
        plot_mesh2d(md, ax = ax, **default_mesh_kwargs)

    ## Add legend
    if show_legend:
        labels_dict = {
            'ice_elements': 'Ice',
            'ocean_elements': 'Ocean',
            'floating_ice_elements': 'Floating ice',
            'grounded_ice_elements': 'Grounded ice',
            'grounding_line_elements': 'Grounding line',
            'ice_front_elements': 'Ice front'}

        legend_elements = [
            matplotlib.patches.Patch(color = color, label = labels_dict[type]),
        ]
        ax.legend(handles = legend_elements, **default_legend_kwargs)

    ## Add axis labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ## Return
    if not ax_defined:
        return fig, ax
    else:
        return ax


def plot_model_field(md,
                     field,
                     layer = None,
                     ax = None,
                     depth_average = False,
                     plot_data_on = 'points',
                     xlabel = 'X (m)',
                     ylabel = 'Y (m)',
                     edgecolors = 'face',
                     figsize = (6.4, 4.8),
                     log = False,
                     constrained_layout = True,
                     show_mesh = False,
                     show_cbar = False,
                     mesh_kwargs = {},
                     cbar_kwargs = {},
                     colormap = None,
                     **kwargs):

    """
    Plot a 2D scalar field defined on a 2D or 3D model mesh.

    This function visualizes a field defined on the mesh of an ISSM model.
    It works both 2D and 3D meshes by extracting a layer from 3D models. If no layer is specified for a 3D model,
    the surface values are plotted by default.

    Parameters
    ----------
    md : ISSM Model object
        ISSM Model object containing mesh. Must be compatible with process_mesh()
    field : np.ndarray
        A 1D array representing the scalar field to plot. Must be defined on vertices or elements. For 3D models, may be defined across all layers.
    layer : int, optional
        Index of the horizontal layer to extract for 3D models (1-based indexing). Ignored for 2D models.
        If not provided, the surface layer is used by default for vertex- and element-based 3D fields.
    ax : matplotlib.axes.Axes, optional
        An existing matplotlib Axes object to plot on. If `None`, a new figure and axes are created.
    depth_average : bool, optional
        If True and the model is 3D, the field is depth-averaged before plotting. Cannot be True if `layer` is specified. Default is False.
    plot_data_on: str, optional
        Should data be plotted on points or elements? Default is 'points'. These options are converted to 'shading' used by plt.tripcolor(), as follows:
        plot_data_on = 'points': shading = 'gouraud' / plot_data_on = 'elements': shading = 'flat'. When data are defined on elements, plot_data_on = 'elements'
        is selected automatically internally.
    xlabel : str, optional
        Label for the x-axis. Default is 'X (m)'.
    ylabel : str, optional
        Label for the y-axis. Default is 'Y (m)'.
    edgecolors : str, optional
        Color of triangle edges in the plot. Passed to 'tripcolor'. Default is 'face'.
    figsize : tuple, optional
        Figure size in inches when creating a new figure. Default is (6.4, 4.8).
    log : bool, optional
        If True, applies a logarithmic color normalization. Default is 'False'.
    constrained_layout : bool, optional
        Whether to use constrained layout for the figure. Default is 'True'.
    show_mesh : bool, optional
        If True, overlays the 2D mesh lines on the plot. Default is 'False'.
    show_cbar : bool, optional
        If True, adds a colorbar to the plot. Default is 'False'.
    mesh_kwargs : dict, optional
        Keyword arguments passed to the mesh plotting function 'plot_mesh2d()'.
    cbar_kwargs : dict, optional
        Keyword arguments passed to 'fig.colorbar'.
    **kwargs : dict, optional
        Additional keyword arguments passed to 'ax.tripcolor'.

    Returns
    -------
    matplotlib.figure.Figure or matplotlib.axes.Axes
        If 'ax' is not provided, returns a tuple '(fig, ax)'. If 'ax' is provided, returns only 'ax'.

    Example
    -------
    fig, ax = plot_model_field(md, md.inversion.vel_obs)
    fig, ax = plot_model_field(md, md.results.TransientSolution.Vel[0], log = True)
    fig, ax = plot_model_field(md, md.results.TransientSolution.Temperature[-1], layer = 3)
    fig, ax = plot_model_field(md, md.results.TransientSolution.Temperature[-1], layer = 3, show_cbar = True, show_mesh = True)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, constrained_layout = True)
    ax1 = plot_model_field(md, md.inversion.vel_obs, ax = ax1, log = True, cmap = 'viridis')
    ax2 = plot_model_field(md, md.results.StressbalanceSolution.Vel, ax = ax2, log = True, cmap = 'viridis')
    ax3 = issm.plot.plot_model_field(md, md.inversion.vel_obs - md.results.StressbalanceSolution.Vel, ax=ax3, cmap='RdBu', show_cbar = True)
    ax1.set_title('Observed Velocities'); ax2.set_title('Modelled Velocities'); ax3.set_title('Velocity Error')
    """

    ## Set defaults
    ax_defined = ax is not None
    norm = matplotlib.colors.LogNorm() if log else None
    shading = 'gouraud' # Consistent with plot_data_on = 'points'

    # Accept 'colormap' as an alias for 'cmap' (consistent with other pyissm functions)
    if colormap is not None and 'cmap' not in kwargs:
        kwargs['cmap'] = colormap

    ## Process mesh
    mesh, mesh_x, mesh_y, mesh_elements, is3d = model.mesh.process_mesh(md)

    ## Error check
    if layer is not None and depth_average:
        raise ValueError('plot_model_field: Cannot specify both a layer and depth_average = True. Please choose one or the other.')

    ## Dimension checks
    if is3d:
        # If a 3D model is provided, extract the layer (if provided), or continue to default extraction of surface layer
        if layer is not None:
            # Extract the specified layer
            field = model.mesh._project_2d(md, field, layer)

        elif depth_average:
            field = model.mesh.depth_average(md, field)

        else:
            # Default behaviour for 3D model with no layer specified
            if field.shape[0] == md.mesh.numberofvertices:
                warnings.warn('plot_model_field: 3D model found with no layer specified. Plotting surface vertices layer only.')
                field = field[md.mesh.vertexonsurface == 1]
            elif field.shape[0] == md.mesh.numberofelements:
                warnings.warn('plot_model_field: 3D model found with no layer specified. Plotting surface elements layer only.')
                field = field[np.isnan(md.mesh.upperelements) == 1]
            elif field.shape[0] == md.mesh.numberofelements2d:
                pass # Field is defined on 2D mesh elements. Already 2D compatible. Continue
            elif field.shape[0] == md.mesh.numberofvertices2d:
                pass # Field is defined on 2D mesh vertices. Already 2D compatible. Continue
            else:
                raise Exception('plot_model_field: The provided field is an unexpected shape.')
    else:
        # If layer is defined, raise warning to explicitly state that it isn't used
        if layer is not None:
            warnings.warn('plot_model_field: 2D model found. Layer definition is ignored.')
        
        # If a 2D model is provided, the field should be defined on vertices or elements.
        if field.shape[0] not in (md.mesh.numberofvertices, md.mesh.numberofelements):
            
            # If the field has one extra row, it is a timestep. Remove this row for plotting.
            if (field.shape[0] == md.mesh.numberofvertices + 1) or (field.shape[0] == md.mesh.numberofelements + 1):
                warnings.warn(f'plot_model_field: Ignoring the timestep value {field[-1]}.')
                field = field[:-1]
            else:
                raise Exception('plot_model_field: The provided field is an unexpected shape.')

    ## Update shading, if necessary. When field is defined on elements, shading = 'flat' is required.
    if (is3d and field.shape[0] == md.mesh.numberofelements2d and plot_data_on == 'points') or (not is3d and field.shape[0] == md.mesh.numberofelements and plot_data_on == 'points'):
        shading = 'flat'
        warnings.warn("Using plot_data_on = 'elements'. Data are defined on elements")
    if plot_data_on == 'elements':
        shading = 'flat'

    ## If no ax is defined, create new figure
    ## otherwise, plot on defined ax
    if ax is None:
        fig, ax = plt.subplots(figsize = figsize, constrained_layout = constrained_layout)
    else:
        fig = ax.get_figure()

    ## Plot field
    trip = ax.tripcolor(mesh, field, edgecolors = edgecolors, norm = norm, **kwargs, shading = shading)

    ## Add optional mesh
    if show_mesh:
        plot_mesh2d(md, ax = ax, **mesh_kwargs)

    ## Add optional colorbar
    if show_cbar:
        fig.colorbar(trip, ax = ax, **cbar_kwargs)

    ## Assign x/y labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ## Return — always yield (fig, ax, trip) so callers can access the mappable
    return fig, ax, trip

def plot_model_bc(md,
                  type = 'stressbalance',
                  ax = None,
                  layer = None,
                  scale = 10,
                  xlabel = 'X (m)',
                  ylabel = 'Y (m)',
                  figsize = (6.4, 4.8),
                  constrained_layout = True,
                  show_mesh = True,
                  mesh_kwargs = {},
                  show_legend = True,
                  legend_kwargs = {}):
    """
    Plot Dirichlet and Neumann boundary conditions from an ISSM model.

    This function visualizes boundary conditions for a specified model
    component (e.g., `stressbalance`, `masstransport`, `thermal`, etc.) on
    a 2D or 3D mesh. Dirichlet conditions are plotted as colored markers,
    and Neumann boundaries (e.g. ice-front) plotted as coloured elements.

    Parameters
    ----------
    md : ISSM Model object
        ISSM Model object containing mesh. Must be compatible with process_mesh()
    type : str, optional
        The boundary condition type to plot. Must be one of:
        'stressbalance', 'masstransport', 'thermal',
        'balancethickness', 'hydrology', 'debris', or 'levelset'.
        Default is 'stressbalance'.
    ax : matplotlib.axes.Axes, optional
        Existing matplotlib Axes object. If not provided, a new figure and
        axes are created.
    layer : int, optional
        Index of the horizontal layer to extract for 3D models (1-based indexing). Ignored for 2D models.
        If not provided, the surface layer is used by default.
    scale : float, optional
        Scaling factor for Dirichlet marker sizes. Default is 10.
    figsize : tuple of float, optional
        Size of the figure in inches (width, height). Default is (6.4, 4.8).
    constrained_layout : bool, optional
        Whether to use constrained layout for figure spacing. Default is True.
    show_mesh : bool, optional
        Whether to display the model mesh beneath boundary markers. Default is True.
    mesh_kwargs : dict, optional
        Additional keyword arguments passed to the mesh plotting function.
        Overrides default {'alpha': 0.5}.
    show_legend : bool, optional
        Whether to display a legend showing boundary condition types. Default is True.
    legend_kwargs : dict, optional
        Additional keyword arguments passed to `ax.legend()`. Overrides default
        {'title': 'Boundary condition', 'fontsize': 10, 'loc': 'upper right'}.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The matplotlib Figure object (only returned if `ax` is not provided).
    ax : matplotlib.axes.Axes
        The matplotlib Axes object containing the plot.

    Notes
    -----
    - For 3D models, only surface boundary conditions are plotted.
    - Neumann (ice-front) elements are included by default as blue elements.
    - If no constraints are found for a given boundary condition, a message
      is printed and nothing is plotted for that type.

    Examples
    --------
    fig, ax = plot_model_bc(md)
    fig, ax = plot_model_bc(md, scale = 1)
    fig, ax = plot_model_bc(md, type='thermal', mesh_kwargs = {'color': 'grey'}, legend_kwargs = {'title': 'Model BCs'})

    See Also
    --------
    plot_model_elements : Used internally to visualize Neumann conditions and mesh.
    """

    ## Set defaults
    ax_defined = ax is not None

    default_mesh_kwargs = {'alpha': 0.5}
    default_mesh_kwargs.update(**mesh_kwargs)

    default_legend_kwargs = {'title': 'Boundary condition',
                           'fontsize': 10,
                           'loc': 'upper right'}
    default_legend_kwargs.update(**legend_kwargs)

    ## Process model mesh
    mesh, mesh_x, mesh_y, mesh_elements, is3d = model.mesh.process_mesh(md)

    ## Get SPC boundaries
    ## -------------------------------------
    
    # Initialize empty dict (used to check supported BCs are present)
    spc_dict = {}

    ## STESSBALANCE
    if type == 'stressbalance':
        spc_dict = {'spcvx': {'data': md.stressbalance.spcvx,
                              'label': 'vx Dirichlet',
                              'col': 'red',
                              'marker': 'o',
                              'size': 10 * scale},
                    'spcvy': {'data': md.stressbalance.spcvy,
                              'label': 'vy Dirichlet',
                              'col': 'blue',
                              'marker': 'o',
                              'size': 6 * scale},
                    'spcvz': {'data': md.stressbalance.spcvz,
                              'label': 'vz Dirichlet',
                              'col': 'yellow',
                              'marker': 'o',
                              'size': 2 * scale}
                    }

    ## MASSTRANSPORT
    if type == 'masstransport':
        spc_dict = {'spcthickness': {'data': md.masstransport.spcthickness,
                                     'label': 'Thickness Dirichlet',
                                     'col': 'red',
                                     'marker': 'o',
                                     'size': 5 * scale}
                    }

    ## THERMAL
    if type == 'thermal':
        spc_dict = {'spctemperature': {'data': md.thermal.spctemperature,
                                       'label': 'Temperature Dirichlet',
                                       'col': 'red',
                                       'marker': 'o',
                                       'size': 5 * scale}
                    }

    ## BALANCE THICKNESS
    if type == 'balancethickness':
        spc_dict = {'spcthickness': {'data': md.balancethickness.spcthickness,
                                     'label': 'Thickness Dirichlet',
                                     'col': 'red',
                                     'marker': 'o',
                                     'size': 5 * scale}
                    }
    
    ## HYDROLOGY
    if type == 'hydrology':

        if isinstance(md.hydrology, model.classes.hydrology.dc):
            spc_dict = {'spcsediment_head': {'data': md.hydrology.spcsediment_head,
                                             'label': 'Sediment water head Dirichlet',
                                             'col': 'red',
                                             'marker': 'o',
                                             'size': 10 * scale},
                             'spcepl_head': {'data': md.hydrology.spcepl_head,
                                             'label': 'EPL water head Dirichlet',
                                             'col': 'blue',
                                             'marker': 'o',
                                             'size': 6 * scale}
                        }
        
        if isinstance(md.hydrology, model.classes.hydrology.glads):
            spc_dict = {'spcphi': {'data': md.hydrology.spcphi,
                                   'label': 'Hydraulic potential Dirichlet',
                                   'col': 'red',
                                   'marker': 'o',
                                   'size': 10 * scale}
                        }
            
        if isinstance(md.hydrology, model.classes.hydrology.shakti):
            spc_dict = {'spchead': {'data': md.hydrology.spchead,
                                    'label': 'Water head Dirichlet',
                                    'col': 'red',
                                    'marker': 'o',
                                    'size': 10 * scale}
                        }
        if isinstance(md.hydrology, model.classes.hydrology.shreve) or isinstance(md.hydrology, model.classes.hydrology.tws):
            spc_dict = {'spcwatercolumn': {'data': md.hydrology.spcwatercolumn,
                                           'label': 'Water column Dirichlet',
                                           'col': 'red',
                                           'marker': 'o',
                                           'size': 5 * scale}
                        }
        
        # If no supported hydrology model is found, spc_dict remains empty and a warning is raised.
        if spc_dict == {}:
            warnings.warn('No supported hydrology model found for plotting. Supported models include: dc, glads, shakti, shreve, and tws.')

    ## DEBRIS
    if type == 'debris':
        spc_dict = {'spcthickness': {'data': md.debris.spcthickness,
                                     'label': 'Thickness Dirichlet',
                                     'col': 'red',
                                     'marker': 'o',
                                     'size': 5 * scale}
                    }

    ## LEVELSET
    if type == 'levelset':
        spc_dict = {'spclevelset': {'data': md.levelset.spclevelset,
                                    'label': 'Levelset Dirichlet',
                                    'col': 'red',
                                    'marker': 'o',
                                    'size': 5 * scale}
                    }
        
    if spc_dict == {}:
        raise ValueError(f'pyissm.plot.plot_model_bc: No supported boundary conditions found for type "{type}". Please check the model and type definition.')

    ## Set-up (or retrieve) figure
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=constrained_layout)
    else:
        fig = ax.get_figure()

    ## Initiate plot with Neumann BCs (ice-front)
    try:
        ax = plot_model_elements(md,
                                md.mask.ice_levelset,
                                md.mask.ocean_levelset,
                                ax = ax,
                                type = 'ice_front_elements',
                                show_mesh = show_mesh,
                                show_legend = False,
                                mesh_kwargs = default_mesh_kwargs)
        neumann_legend = True
    except ValueError:
        print('No ice-front (Neumann) elements found to plot.')
        ax = plot_mesh2d(md, ax = ax, **default_mesh_kwargs)
        neumann_legend = False

    ## Add Dirichlet BCs
    for key, spc in spc_dict.items():
        data = spc['data']

        # If 2D (e.g. spcwatercolumn), use first column
        if data.ndim == 2:
            data = data[:, 0]

        # If the data are all NaN, there are no constraints
        if np.isnan(data).all():
            print(f'No constraints found in {key}')
            pass
        else:
            # If model is 3D, extract the BCs on the specified layer (or surface if no layer specified)
            if is3d:
                if layer is None:
                    data = data[md.mesh.vertexonsurface == 1]
                    warnings.warn(f'3D model found. Plotting surface BCs only.')
                else:
                    data = model.mesh._project_2d(md, data, layer)
                    warnings.warn(f'3D model found. Plotting BCs  on layer {layer}.')

            # Make plot
            ax.scatter(mesh_x[~np.isnan(data)],
                       mesh_y[~np.isnan(data)],
                       c = spc['col'],
                       marker = spc['marker'],
                       s = spc['size'],
                       label = spc['label'])

    ## Add optional legend
    if show_legend:
        handles, labels = ax.get_legend_handles_labels()

        # Add Neumann entry only if it exists
        if neumann_legend:
            ice_front = matplotlib.patches.Patch(color = 'blue',
                                                 label ='Neumann (ice-front)')
            handles = [ice_front] + handles
        
        ax.legend(handles = handles, **default_legend_kwargs)

    ## Add axis labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ## Return
    if not ax_defined:
        return fig, ax
    else:
        return ax
    
def plot_model_ts(md,
                  data_list = None,
                  variable_list = None,
                  unit_list = None,
                  time = (0, np.inf),
                  xlabel = 'Time (years)',
                  figsize = (6.4, 4.8),
                  fontsize = 8,
                  color = 'black',
                  show_grid = True,
                  constrained_layout = True,
                  sharex = True,
                  **kwargs):
    """
    Plot time series data from an ISSM model's TransientSolution.

    Automatically extracts 1D time series arrays from `md.results.TransientSolution` if `data_list` is not provided.
    Plots each variable as a separate subplot, sharing the x-axis (time).
    
    Parameters
    ----------
    md : object
        ISSM model object containing results with a `TransientSolution` attribute.
    data_list : list of np.ndarray, optional
        List of 1D numpy arrays to plot. If None, all 1D arrays in `md.results.TransientSolution` (except 'time' and 'step') are used.
    variable_list : list of str, optional
        List of variable names corresponding to `data_list`. If None, names are auto-generated or extracted from `TransientSolution`.
    unit_list : list of str, optional
        List of units for each variable. If None, defaults to 'Units unspecified' or 'ISSM Units'.
    time : tuple of float, optional
        Time range (start, end) to plot. Defaults to (0, np.inf).
    figsize : tuple of int, optional
        Figure size in inches. Defaults to (6.4, 4.8).
    fontsize : int, optional
        Font size for labels and titles. Defaults to 8.
    color : str, optional
        Line color for plots. Defaults to 'black'.
    show_grid : bool, optional
        Whether to show grid lines on plots. Defaults to True.
    constrained_layout : bool, optional
        Whether to use matplotlib's constrained layout. Defaults to True.
    sharex : bool, optional
        Whether subplots share the x-axis. Defaults to True.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The matplotlib figure object containing the plots.
    axs : np.ndarray of matplotlib.axes.Axes
        Array of axes objects for each subplot.
    
    Raises
    ------
    ValueError
        If no time array is found, no 1D time series are available, input lists are mismatched, or data arrays are not 1D.

    Examples
    --------
    >>> fig, axs = plot_model_ts(md)
    >>> fig, axs = plot_model_ts(md, data_list=[arr1, arr2], variable_list=['Var1', 'Var2'], unit_list=['m', 'kg'])
    """

    # Check if md has the required attributes
    if hasattr(md, 'results') is False or not hasattr(md.results, 'TransientSolution'):
        raise ValueError("The provided model does not have a 'results.TransientSolution' attribute.")

    # Extract model time from the TransientSolution
    model_time = getattr(md.results.TransientSolution, 'time', None)
    if model_time is None:
        raise ValueError("No 'time' array found in md.results.TransientSolution.")

    # Auto-extract 1D arrays if no data_list provided
    if data_list is None:
        variables = {}
        for k, v in md.results.TransientSolution._field_major_cache.items():
            # Skip specific entries
            if k in ('step', 'time', 'StressbalanceConvergenceNumSteps', 'outlog', 'errlog'):
                continue
            # If it's a 1D array, add to variables dict
            if isinstance(v, np.ndarray) and v.ndim == 1:
                variables[k] = v
        
        # If no 1D time series found, raise an error
        if not variables:
            raise ValueError("No 1D time series found in md.results.TransientSolution.")

        # Construct data_list, variable_list, and unit_list from the variables dictionary. Use 'ISSM units'.    
        data_list = list(variables.values())
        variable_list = list(variables.keys())
        unit_list = ['ISSM Units'] * len(data_list)

    # Check that data_list is not empty
    n_vars = len(data_list)
    if n_vars == 0:
        raise ValueError("data_list must contain at least one 1D numpy array.")

    # Fill missing variable names/units
    if variable_list is None:
        variable_list = [f'Variable {i + 1}' for i in range(n_vars)]
        
    if unit_list is None:
        unit_list = ['Units unspecified'] * n_vars

    # Validate input lengths
    if not (len(data_list) == len(variable_list) == len(unit_list)):
        raise ValueError("data_list, variable_list, and unit_list must all be the same length.")

    # Check that all data arrays are 1D
    if not all(arr.ndim == 1 for arr in data_list):
        raise ValueError("All data in data_list must be 1D numpy arrays.")

    # Filter time range
    time_idx = np.where((model_time >= time[0]) & (model_time <= time[1]))[0]
    model_time = model_time[time_idx]

    # Create subplots & ensure axs is iterable even if n_vars == 1
    fig, axs = plt.subplots(n_vars, 1, figsize=figsize, constrained_layout=constrained_layout, sharex=sharex)
    axs = np.atleast_1d(axs)

    # Plot each variable iteratively
    for data, variable, unit, ax in zip(data_list, variable_list, unit_list, axs):
        ax.plot(model_time, data[time_idx], color = color, **kwargs)
        ax.set_title(variable, fontsize = fontsize)
        ax.set_ylabel(f'{variable}\n({unit})', fontsize = fontsize)
        ax.tick_params(axis = 'both', labelsize = fontsize)
        if show_grid:
            ax.grid(alpha = 0.4)

    # Add a common x-label for the last subplot
    axs[-1].set_xlabel(xlabel, fontsize = fontsize)

    return fig, axs

## ------------------------------------------------------------------------------------
## CONTOUR, STREAMLINES, QUIVER

def plot_contour(md, field, levels=10, ax=None, colors='y', linestyles='-',
                 linewidths=1.0, label=False, label_fmt='%.2f', label_fontsize=10,
                 figsize=(6.4, 4.8), constrained_layout=True, **kwargs):
    """Overlay contour lines of a scalar field on a triangular mesh.

    Parameters
    ----------
    md : model
        ISSM model object.
    field : array-like, shape (nvertices,) or (nvertices, ntimesteps)
        Scalar field to contour. If 2-D the first time-step is used.
    levels : int or array-like
        Number of levels or explicit level values.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None a new figure is created.
    colors : str or list
        Contour line colour(s). Default 'y'.
    linestyles : str or list
        Contour line style(s). Default '-'.
    linewidths : float or list
        Contour line width(s). Default 1.0.
    label : bool
        Add inline contour labels. Default False.
    label_fmt : str
        Format string for contour labels. Default '%.2f'.
    label_fontsize : float
        Font size for contour labels. Default 10.
    figsize : tuple
        Figure size when ax is None. Default (6.4, 4.8).
    constrained_layout : bool
        Use constrained layout when creating a new figure. Default True.

    Returns
    -------
    (fig, ax) when ax was None, otherwise ax.
    """
    field = np.asarray(field, dtype=float)
    if field.ndim == 2:
        field = field[:, 0]

    x = md.mesh.x
    y = md.mesh.y
    elements = md.mesh.elements - 1  # 1-based → 0-based

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=constrained_layout)

    cs = ax.tricontour(x, y, elements, field, levels=levels, colors=colors,
                       linestyles=linestyles, linewidths=linewidths, **kwargs)
    if label:
        ax.clabel(cs, fmt=label_fmt, fontsize=label_fontsize, inline=True)

    return (fig, ax) if own_fig else ax


def plot_streamlines(md, vx=None, vy=None, ax=None, grid_n=200, color='k',
                     linewidth=1.0, density=1.0, arrowsize=1.0, speed_scale=False,
                     speed_scale_factor=1.0, figsize=(6.4, 4.8),
                     constrained_layout=True, **kwargs):
    """Draw streamlines of a 2-D velocity field on the mesh.

    Requires *scipy*.

    Parameters
    ----------
    md : model
        ISSM model object.
    vx, vy : array-like, shape (nvertices,), optional
        X / Y velocity components. Defaults to md.initialization.vx / vy.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None a new figure is created.
    grid_n : int
        Resolution of the regular grid for streamplot. Default 200.
    color : str
        Streamline colour. Default 'k'.
    linewidth : float or array-like
        Streamline width. Default 1.0. Set to None to auto-scale by speed
        when *speed_scale* is True.
    density : float
        Streamline density passed to ``ax.streamplot``. Default 1.0.
    arrowsize : float
        Arrow size multiplier. Default 1.0.
    speed_scale : bool
        If True, scale linewidth by speed magnitude. Default False.
    speed_scale_factor : float
        Multiplier for speed-scaled linewidths. Default 1.0.
    figsize : tuple
        Figure size when ax is None. Default (6.4, 4.8).
    constrained_layout : bool
        Use constrained layout when creating a new figure. Default True.

    Returns
    -------
    (fig, ax) when ax was None, otherwise ax.
    """
    try:
        from scipy.interpolate import griddata
    except ImportError as err:
        raise ImportError("plot_streamlines requires scipy") from err

    if vx is None:
        vx = md.initialization.vx
    if vy is None:
        vy = md.initialization.vy

    vx = np.asarray(vx, dtype=float)
    vy = np.asarray(vy, dtype=float)
    if vx.ndim == 2:
        vx = vx[:, 0]
    if vy.ndim == 2:
        vy = vy[:, 0]

    x = md.mesh.x
    y = md.mesh.y

    xg = np.linspace(x.min(), x.max(), grid_n)
    yg = np.linspace(y.min(), y.max(), grid_n)
    Xg, Yg = np.meshgrid(xg, yg)
    pts = np.column_stack((x, y))
    ug = griddata(pts, vx, (Xg, Yg), method='linear')
    vg = griddata(pts, vy, (Xg, Yg), method='linear')

    # Mask NaNs so streamplot won't crash
    ug = np.where(np.isnan(ug), 0.0, ug)
    vg = np.where(np.isnan(vg), 0.0, vg)

    if speed_scale and linewidth is None:
        speed = np.sqrt(ug ** 2 + vg ** 2)
        smax = speed.max()
        lw = speed_scale_factor * (speed / smax) if smax > 0 else np.ones_like(speed)
    else:
        lw = linewidth

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=constrained_layout)

    ax.streamplot(xg, yg, ug, vg, color=color, linewidth=lw, density=density,
                  arrowsize=arrowsize, **kwargs)

    return (fig, ax) if own_fig else ax


def plot_quiver(md, vx=None, vy=None, ax=None, color='k', scale=None,
                width=5e-3, headwidth=6.0, headlength=None,
                normalize=False, figsize=(6.4, 4.8),
                constrained_layout=True, **kwargs):
    """Draw a quiver (arrow) plot of a 2-D velocity field on the mesh.

    Parameters
    ----------
    md : model
        ISSM model object.
    vx, vy : array-like, shape (nvertices,), optional
        X / Y velocity components. Defaults to md.initialization.vx / vy.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None a new figure is created.
    color : str or array-like
        Arrow colour or 'speed' to colour by velocity magnitude.
    scale : float, optional
        Scaling factor for arrows. Auto-computed from domain size if None.
    width : float
        Arrow shaft width in axes fraction. Default 5e-3.
    headwidth : float
        Head width relative to shaft. Default 6.
    headlength : float, optional
        Head length. Defaults to headwidth * 1.5.
    normalize : bool
        If True all arrows have equal length (unit vectors). Default False.
    figsize : tuple
        Figure size when ax is None. Default (6.4, 4.8).
    constrained_layout : bool
        Use constrained layout when creating a new figure. Default True.

    Returns
    -------
    (fig, ax) when ax was None, otherwise ax.
    """
    if vx is None:
        vx = md.initialization.vx
    if vy is None:
        vy = md.initialization.vy

    vx = np.asarray(vx, dtype=float)
    vy = np.asarray(vy, dtype=float)
    if vx.ndim == 2:
        vx = vx[:, 0]
    if vy.ndim == 2:
        vy = vy[:, 0]

    x = md.mesh.x
    y = md.mesh.y

    # Auto-scale: scale arrows so max arrow is ~sqrt(domain area / nvertices)
    if scale is None:
        speed = np.sqrt(vx ** 2 + vy ** 2)
        max_speed = speed.max()
        if max_speed > 0:
            xdist = x.max() - x.min()
            ydist = y.max() - y.min()
            ref_len = np.sqrt(xdist * ydist / len(x))
            scale = max_speed / ref_len
        else:
            scale = 1.0

    if normalize:
        speed = np.sqrt(vx ** 2 + vy ** 2)
        speed = np.where(speed == 0, 1.0, speed)
        vx = vx / speed
        vy = vy / speed
        scale = scale / speed.max() if not isinstance(scale, float) else scale

    if headlength is None:
        headlength = headwidth * 1.5

    use_array_color = color == 'speed'
    if use_array_color:
        c = np.sqrt(vx ** 2 + vy ** 2)
    else:
        c = color

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=constrained_layout)

    if use_array_color:
        # Pass speed as the positional C array for array-based colouring
        ax.quiver(x, y, vx, vy, c, scale=scale, scale_units='xy', width=width,
                  headwidth=headwidth, headlength=headlength, **kwargs)
    else:
        # Pass colour as a keyword so matplotlib doesn't try np.isfinite on a string
        ax.quiver(x, y, vx, vy, scale=scale, scale_units='xy', width=width, color=c,
                  headwidth=headwidth, headlength=headlength, **kwargs)

    return (fig, ax) if own_fig else ax


## ------------------------------------------------------------------------------------
## EDGE OVERLAY, ICE FRONT, COASTLINES

def plot_edgeoverlay(md, data, ax=None, edgemin=None, edgetype='thickness',
                     edgeranges=2.0, colormap='jet',
                     figsize=(6.4, 4.8), constrained_layout=True, **kwargs):
    """Overlay coloured or thickness-weighted mesh edges on a plot.

    Requires md.mesh.edges (edge connectivity array).

    Parameters
    ----------
    md : model
        ISSM model object.
    data : array-like, shape (nedges,) or (nvertices,)
        Scalar values to colour / weight the edges. If per-vertex, the
        edge value is taken as the mean of its two endpoint values.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None a new figure is created.
    edgemin : float, optional
        Minimum quantile threshold (0–1). Edges below are not drawn.
        Default 0 (draw all edges).
    edgetype : str
        'color'     — colour edges by data value (uses ax.quiver colour hack).
        'thickness' — vary line thickness by data value (default).
    edgeranges : float
        Maximum line width or colour magnitude for the thickest/brightest edges.
        Default 2.
    colormap : str or Colormap
        Colormap for 'thickness' mode. Default 'jet'.
    figsize : tuple
        Figure size when ax is None. Default (6.4, 4.8).
    constrained_layout : bool
        Use constrained layout when creating a new figure. Default True.

    Returns
    -------
    (fig, ax) when ax was None, otherwise ax.
    """
    from matplotlib.collections import LineCollection
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors

    x = md.mesh.x
    y = md.mesh.y

    if not hasattr(md.mesh, 'edges') or md.mesh.edges is None:
        raise AttributeError("md.mesh.edges is required for plot_edgeoverlay")

    edges = np.asarray(md.mesh.edges, dtype=int) - 1  # 1-based → 0-based
    # edges shape: (nedges, 2+) — first two columns are vertex indices
    e0 = edges[:, 0]
    e1 = edges[:, 1]

    data = np.asarray(data, dtype=float)
    if data.ndim == 2:
        data = data[:, 0]

    # Per-vertex → per-edge
    if data.shape[0] == len(x):
        edge_vals = 0.5 * (data[e0] + data[e1])
    elif data.shape[0] == len(e0):
        edge_vals = data
    else:
        raise ValueError("data length must match number of edges or vertices")

    if edgemin is not None:
        threshold = np.nanquantile(edge_vals, edgemin)
        mask = edge_vals >= threshold
        e0, e1, edge_vals = e0[mask], e1[mask], edge_vals[mask]

    segs = np.stack((np.column_stack((x[e0], y[e0])),
                     np.column_stack((x[e1], y[e1]))), axis=1)

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=constrained_layout)

    ev_min, ev_max = edge_vals.min(), edge_vals.max()
    ev_range = ev_max - ev_min if ev_max != ev_min else 1.0
    ev_norm = (edge_vals - ev_min) / ev_range  # 0-1

    if edgetype == 'color':
        cmap = get_colormap(colormap)
        lc = LineCollection(segs, cmap=cmap, norm=mcolors.Normalize(ev_min, ev_max),
                            **kwargs)
        lc.set_array(edge_vals)
        ax.add_collection(lc)
    else:  # thickness
        lw = ev_norm * edgeranges
        lw = np.clip(lw, 0.1, None)
        cmap = get_colormap(colormap)
        rgba = cmap(ev_norm)
        lc = LineCollection(segs, linewidths=lw, colors=rgba, **kwargs)
        ax.add_collection(lc)

    ax.autoscale()
    return (fig, ax) if own_fig else ax


def plot_icefront(md, ax=None, ice_color='tab:blue', gl_color='tab:green',
                  alpha=0.6, figsize=(6.4, 4.8), constrained_layout=True, **kwargs):
    """Highlight ice-front and grounding-line elements on the mesh.

    Ice-front elements (elements at the ocean–ice boundary) are coloured
    *ice_color*. Grounding-line elements are coloured *gl_color*.

    Parameters
    ----------
    md : model
        ISSM model object (must have md.mask.ice_levelset and
        md.mask.ocean_levelset).
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None a new figure is created.
    ice_color : str
        Colour for ice-front elements. Default 'tab:blue'.
    gl_color : str
        Colour for grounding-line elements. Default 'tab:green'.
    alpha : float
        Alpha transparency. Default 0.6.
    figsize : tuple
        Figure size when ax is None. Default (6.4, 4.8).
    constrained_layout : bool
        Use constrained layout when creating a new figure. Default True.

    Returns
    -------
    (fig, ax) when ax was None, otherwise ax.
    """
    x = md.mesh.x
    y = md.mesh.y
    elements = md.mesh.elements - 1  # 0-based

    ice_ls = np.asarray(md.mask.ice_levelset, dtype=float)
    ocean_ls = np.asarray(md.mask.ocean_levelset, dtype=float)

    # Element-mean levelset values
    ice_elem = ice_ls[elements].mean(axis=1)
    ocean_elem = ocean_ls[elements].mean(axis=1)

    # Ice-front: element spans ice/ocean boundary (sign change in ocean levelset)
    # A sign change means not all vertices have the same sign
    def _has_sign_change(ls, elems):
        mins = ls[elems].min(axis=1)
        maxs = ls[elems].max(axis=1)
        return (mins < 0) & (maxs > 0)

    ice_front_mask = _has_sign_change(ocean_ls, elements)
    gl_mask = _has_sign_change(ice_ls, elements)

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=constrained_layout)

    # Draw base mesh
    ax.triplot(x, y, elements, color='0.8', linewidth=0.3, zorder=1)

    if ice_front_mask.any():
        ice_values = np.zeros(len(elements))
        ice_values[ice_front_mask] = 1.0
        ax.tripcolor(x, y, elements, facecolors=ice_values,
                     cmap=matplotlib.colors.ListedColormap(['none', ice_color]),
                     vmin=0, vmax=1, alpha=alpha, zorder=2)

    if gl_mask.any():
        gl_values = np.zeros(len(elements))
        gl_values[gl_mask] = 1.0
        ax.tripcolor(x, y, elements, facecolors=gl_values,
                     cmap=matplotlib.colors.ListedColormap(['none', gl_color]),
                     vmin=0, vmax=1, alpha=alpha, zorder=3)

    import matplotlib.patches as mpatches
    handles = []
    if ice_front_mask.any():
        handles.append(mpatches.Patch(color=ice_color, alpha=alpha, label='Ice front'))
    if gl_mask.any():
        handles.append(mpatches.Patch(color=gl_color, alpha=alpha, label='Grounding line'))
    if handles:
        ax.legend(handles=handles)

    return (fig, ax) if own_fig else ax


def plot_coastlines(md, ax=None, color='k', linewidth=0.5, resolution='10m',
                    figsize=(6.4, 4.8), constrained_layout=True, **kwargs):
    """Overlay coastline geometries on the model domain.

    Uses *cartopy* (preferred) when available. Falls back to *shapely* / Natural
    Earth shapefiles. If neither is installed a warning is emitted and the axes
    are returned unchanged.

    Parameters
    ----------
    md : model
        ISSM model object.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None a new figure is created.
    color : str
        Coastline colour. Default 'k'.
    linewidth : float
        Coastline line width. Default 0.5.
    resolution : str
        Natural Earth resolution: '10m', '50m', or '110m'. Default '10m'.
    figsize : tuple
        Figure size when ax is None. Default (6.4, 4.8).
    constrained_layout : bool
        Use constrained layout when creating a new figure. Default True.

    Returns
    -------
    (fig, ax) when ax was None, otherwise ax.
    """
    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=constrained_layout)

    try:
        import cartopy.feature as cfeature
        import cartopy.crs as ccrs
        # If the axes is a GeoAxes, add the feature directly
        if hasattr(ax, 'add_feature'):
            res_map = {'10m': '10m', '50m': '50m', '110m': '110m'}
            scale = res_map.get(resolution, '50m')
            coast = cfeature.NaturalEarthFeature('physical', 'coastline', scale,
                                                 edgecolor=color, facecolor='none',
                                                 linewidth=linewidth)
            ax.add_feature(coast, **kwargs)
            return (fig, ax) if own_fig else ax
    except ImportError:
        pass

    try:
        import shapefile
    except ImportError:
        warnings.warn(
            "plot_coastlines: neither cartopy nor pyshp (shapefile) is installed. "
            "Coastlines cannot be drawn.",
            stacklevel=2,
        )
        return (fig, ax) if own_fig else ax

    # Fallback: draw model mesh boundary as a rough coastline proxy
    try:
        from pyissm.tools.general import mesh_coastline  # may not exist
        coords = mesh_coastline(md)
        for segment in coords:
            ax.plot(segment[:, 0], segment[:, 1], color=color,
                    linewidth=linewidth, **kwargs)
    except Exception:
        warnings.warn(
            "plot_coastlines: could not load coastline data. "
            "Install cartopy for full functionality.",
            stacklevel=2,
        )

    return (fig, ax) if own_fig else ax


## ------------------------------------------------------------------------------------
## GRIDDED DATA AND GEOTIFF OVERLAY

def plot_gridded(md, field, ax=None, xlim=None, ylim=None, posting=None,
                 caxis=None, colormap='viridis', log=False, show_mesh=False,
                 figsize=(6.4, 4.8), constrained_layout=True, **kwargs):
    """Interpolate a vertex field to a regular grid and display as an image.

    Requires *scipy*.

    Parameters
    ----------
    md : model
        ISSM model object.
    field : array-like, shape (nvertices,) or (nvertices, ntimesteps)
        Scalar field to display. If 2-D the first time-step is used.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None a new figure is created.
    xlim, ylim : array-like [min, max], optional
        Axis limits. Defaults to mesh extent.
    posting : float, optional
        Grid spacing. Defaults to ~300 grid points along the longer axis.
    caxis : array-like [vmin, vmax], optional
        Colour axis limits.
    colormap : str or Colormap
        Colormap. Default 'viridis'.
    log : bool
        If True apply log10 to data before plotting. Default False.
    show_mesh : bool
        If True overlay the mesh edges. Default False.
    figsize : tuple
        Figure size when ax is None. Default (6.4, 4.8).
    constrained_layout : bool
        Use constrained layout when creating a new figure. Default True.

    Returns
    -------
    (fig, ax, im) when ax was None, otherwise (ax, im).
    """
    try:
        from scipy.interpolate import griddata
    except ImportError as err:
        raise ImportError("plot_gridded requires scipy") from err

    field = np.asarray(field, dtype=float)
    if field.ndim == 2:
        field = field[:, 0]

    x = md.mesh.x
    y = md.mesh.y

    if xlim is None:
        xlim = [x.min(), x.max()]
    if ylim is None:
        ylim = [y.min(), y.max()]

    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]

    if posting is None:
        posting = max(xrange, yrange) / 300.0

    nx = max(2, int(np.round(xrange / posting)))
    ny = max(2, int(np.round(yrange / posting)))

    xg = np.linspace(xlim[0], xlim[1], nx)
    yg = np.linspace(ylim[0], ylim[1], ny)
    Xg, Yg = np.meshgrid(xg, yg)

    pts = np.column_stack((x, y))
    Zg = griddata(pts, field, (Xg, Yg), method='linear')

    if log:
        Zg = np.log10(np.abs(Zg))

    vmin, vmax = (None, None) if caxis is None else (caxis[0], caxis[1])
    cmap = get_colormap(colormap)

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=constrained_layout)

    im = ax.imshow(Zg, extent=[xlim[0], xlim[1], ylim[0], ylim[1]],
                   origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, **kwargs)

    if show_mesh:
        elements = md.mesh.elements - 1
        ax.triplot(x, y, elements, color='k', linewidth=0.2, alpha=0.3)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    return (fig, ax, im) if own_fig else (ax, im)


def plot_overlay(md, geotiff_path, ax=None, xlim=None, ylim=None,
                 overlaylims=None, greyscale=True, alpha=1.0,
                 figsize=(6.4, 4.8), constrained_layout=True, **kwargs):
    """Overlay a GeoTIFF raster image on the model domain.

    Requires *osgeo.gdal*.

    Parameters
    ----------
    md : model
        ISSM model object (used for default xlim/ylim from mesh extent).
    geotiff_path : str
        Path to the GeoTIFF file.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None a new figure is created.
    xlim, ylim : array-like [min, max], optional
        Domain extent in model coordinates. Defaults to mesh extent.
    overlaylims : array-like [min, max], optional
        Intensity clipping limits for greyscale images.
    greyscale : bool
        Convert to greyscale. Default True.
    alpha : float
        Transparency of the overlay. Default 1.0.
    figsize : tuple
        Figure size when ax is None. Default (6.4, 4.8).
    constrained_layout : bool
        Use constrained layout when creating a new figure. Default True.

    Returns
    -------
    (fig, ax) when ax was None, otherwise ax.
    """
    try:
        from osgeo import gdal
    except ImportError as err:
        raise ImportError("plot_overlay requires osgeo.gdal (GDAL)") from err

    ds = gdal.Open(geotiff_path)
    if ds is None:
        raise IOError(f"plot_overlay: could not open '{geotiff_path}'")

    gt = ds.GetGeoTransform()
    # gt: (x_min, pixel_w, rotation, y_max, rotation, pixel_h) — pixel_h is negative
    xmin_r = gt[0]
    xmax_r = gt[0] + ds.RasterXSize * gt[1]
    ymax_r = gt[3]
    ymin_r = gt[3] + ds.RasterYSize * gt[5]

    band = ds.GetRasterBand(1)
    data = band.ReadAsArray().astype(float)
    nodata = band.GetNoDataValue()
    if nodata is not None:
        data[data == nodata] = np.nan

    if greyscale:
        if overlaylims is not None:
            data = np.clip(data, overlaylims[0], overlaylims[1])
            dmin, dmax = overlaylims
        else:
            finite = data[np.isfinite(data)]
            dmin, dmax = finite.min(), finite.max()
        if dmax != dmin:
            img = (data - dmin) / (dmax - dmin)
        else:
            img = np.zeros_like(data)
        cmap_img = 'gray'
    else:
        img = data
        cmap_img = 'viridis'

    if xlim is None:
        xlim = [md.mesh.x.min(), md.mesh.x.max()]
    if ylim is None:
        ylim = [md.mesh.y.min(), md.mesh.y.max()]

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=constrained_layout)

    ax.imshow(img, extent=[xmin_r, xmax_r, ymin_r, ymax_r],
              origin='upper', cmap=cmap_img, alpha=alpha,
              interpolation='bilinear', **kwargs)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    return (fig, ax) if own_fig else ax


## ------------------------------------------------------------------------------------
## MESH ELEMENT AND VERTEX NUMBERING

def plot_elementnumbering(md, ax=None, highlight=None, fontsize=6,
                          figsize=(6.4, 4.8), constrained_layout=True, **kwargs):
    """Plot the mesh with element numbers labelled at centroids.

    Parameters
    ----------
    md : model
        ISSM model object.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None a new figure is created.
    highlight : list of int, optional
        1-based element indices to highlight in red. Others shown in black.
    fontsize : float
        Font size for labels. Default 6.
    figsize : tuple
        Figure size when ax is None. Default (6.4, 4.8).
    constrained_layout : bool
        Use constrained layout when creating a new figure. Default True.

    Returns
    -------
    (fig, ax) when ax was None, otherwise ax.
    """
    x = md.mesh.x
    y = md.mesh.y
    elements = md.mesh.elements - 1  # 0-based

    highlight_set = set()
    if highlight is not None:
        highlight_set = {int(h) - 1 for h in highlight}  # 0-based

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=constrained_layout)

    ax.triplot(x, y, elements, color='0.5', linewidth=0.5)

    cx = x[elements].mean(axis=1)
    cy = y[elements].mean(axis=1)
    for i, (xi, yi) in enumerate(zip(cx, cy)):
        c = 'red' if i in highlight_set else 'black'
        ax.text(xi, yi, str(i + 1), ha='center', va='center',
                fontsize=fontsize, color=c, **kwargs)

    return (fig, ax) if own_fig else ax


def plot_vertexnumbering(md, ax=None, highlight=None, fontsize=6,
                         figsize=(6.4, 4.8), constrained_layout=True, **kwargs):
    """Plot the mesh with vertex numbers labelled at each node.

    Parameters
    ----------
    md : model
        ISSM model object.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None a new figure is created.
    highlight : list of int, optional
        1-based vertex indices to highlight in red. Others shown in white.
    fontsize : float
        Font size for labels. Default 6.
    figsize : tuple
        Figure size when ax is None. Default (6.4, 4.8).
    constrained_layout : bool
        Use constrained layout when creating a new figure. Default True.

    Returns
    -------
    (fig, ax) when ax was None, otherwise ax.
    """
    x = md.mesh.x
    y = md.mesh.y
    elements = md.mesh.elements - 1  # 0-based

    highlight_set = set()
    if highlight is not None:
        highlight_set = {int(h) - 1 for h in highlight}  # 0-based

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=constrained_layout)

    ax.triplot(x, y, elements, color='0.5', linewidth=0.5)

    for i, (xi, yi) in enumerate(zip(x, y)):
        bg = 'red' if i in highlight_set else 'white'
        ax.text(xi, yi, str(i + 1), ha='center', va='center',
                fontsize=fontsize,
                bbox=dict(boxstyle='circle,pad=0.1', facecolor=bg,
                          edgecolor='0.5', linewidth=0.3),
                **kwargs)

    return (fig, ax) if own_fig else ax
