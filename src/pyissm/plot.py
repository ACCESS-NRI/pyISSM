"""
Functions to visualize ISSM Models
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from . import model, utils
import warnings

## ------------------------------------------------------------------------------------
## MESH PLOTTING
## ------------------------------------------------------------------------------------
def plot_mesh2d(mesh,
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
    mesh : matplotlib.tri.Triangulation or similar
        A 2D mesh object containing 'x', 'y', and triangle connectivity information.
        Can be created for ISSM models with get_mesh() or process_mesh().
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
    fig, ax = plot_mesh2d(mesh)
    fig, ax = plot_mesh2d(mesh, color = 'blue', linewidth = 0.5)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1 = plot_mesh2d(mesh, ax = ax1)
    ax2 = plot_mesh2d(mesh, ax = ax2, show_nodes = True, node_kwargs = {'color': 'red'})
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
        plot_mesh2d(mesh,
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
        plot_mesh2d(mesh, ax = ax, **default_mesh_kwargs)

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

    ## Process mesh
    mesh, mesh_x, mesh_y, mesh_elements, is3d = model.mesh.process_mesh(md)

    ## Dimension checks
    if is3d:
        # If a 3D model is provided, extract the layer (if provided), or continue to default extraction of surface layer
        if layer is not None:
            # Extract the specified layer
            field, _ = utils.general.extract_field_layer(md, field, layer)
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
    if (is3d and field.shape == md.mesh.numberofelements2d and plot_data_on == 'points') or (not is3d and field.shape == md.mesh.numberofelements and plot_data_on == 'points'):
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
        plot_mesh2d(mesh, ax = ax, **mesh_kwargs)

    ## Add optional colorbar
    if show_cbar:
        fig.colorbar(trip, ax = ax, **cbar_kwargs)

    ## Assign x/y labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ## Return
    if not ax_defined:
        return fig, ax
    else:
        return ax

def plot_model_bc(md,
                  type = 'stressbalance',
                  ax = None,
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

    if type == 'masstransport':
        spc_dict = {'spcthickness': {'data': md.masstransport.spcthickness,
                                     'label': 'Thickness Dirichlet',
                                     'col': 'red',
                                     'marker': 'o',
                                     'size': 5 * scale}
                    }

    if type == 'thermal':
        spc_dict = {'spctemperature': {'data': md.thermal.spctemperature,
                                       'label': 'Temperature Dirichlet',
                                       'col': 'red',
                                       'marker': 'o',
                                       'size': 5 * scale}
                    }

    if type == 'balancethickness':
        spc_dict = {'spcthickness': {'data': md.balancethickness.spcthickness,
                                     'label': 'Thickness Dirichlet',
                                     'col': 'red',
                                     'marker': 'o',
                                     'size': 5 * scale}
                    }

    if type == 'hydrology':
        if md.hydrology.__class__.__name__.lower() == 'glads':
            spc_dict = {'spcphi': {'data': md.hydrology.spcphi,
                    'label': 'phi Dirichlet',
                    'col': 'green',
                    'marker': 's',
                    'size': 8 * scale}
                    }
        else:
           spc_dict = {'spcwatercolumn': {'data': md.hydrology.spcwatercolumn,
                        'label': 'Water column Dirichlet',
                        'col': 'red',
                        'marker': 'o',
                        'size': 5 * scale}
                        }

    if type == 'debris':
        spc_dict = {'spcthickness': {'data': md.debris.spcthickness,
                                     'label': 'Thickness Dirichlet',
                                     'col': 'red',
                                     'marker': 'o',
                                     'size': 5 * scale}
                    }

    if type == 'levelset':
        spc_dict = {'spclevelset': {'data': md.levelset.spclevelset,
                                    'label': 'Levelset Dirichlet',
                                    'col': 'red',
                                    'marker': 'o',
                                    'size': 5 * scale}
                    }

    ## Set-up (or retrieve) figure
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=constrained_layout)
    else:
        fig = ax.get_figure()

    ## Initiate plot with Neumann BCs (ice-front)
    ax = plot_model_elements(md,
                             md.mask.ice_levelset,
                             md.mask.ocean_levelset,
                             ax = ax,
                             type = 'ice_front_elements',
                             show_mesh = show_mesh,
                             show_legend = False,
                             mesh_kwargs = default_mesh_kwargs)

    ## Add Dirichlet BCs
    for key, spc in spc_dict.items():
        data = spc['data']

        # If the data are all NaN, there are no constraints
        if np.isnan(data).all():
            print(f'No constraints found in {key}')
            pass
        else:
            # If model is 3D, extract the BCs on the surface layer
            if is3d:
                data = data[md.mesh.vertexonsurface == 1]
                warnings.warn(f'3D model found. Plotting surface BCs only.')

            # Make plot
            ax.scatter(mesh_x[~np.isnan(data)],
                       mesh_y[~np.isnan(data)],
                       c = spc['col'],
                       marker = spc['marker'],
                       s = spc['size'],
                       label = spc['label'])

    ## Add optional legend (including manual entry for Neumann ice-front)
    if show_legend:
        ice_front = matplotlib.patches.Patch(color = 'blue', label ='Neumann (ice-front)')
        ax.legend(handles=[ice_front] + ax.get_legend_handles_labels()[0], **default_legend_kwargs)

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
        for k, v in vars(md.results.TransientSolution).items():
            # If the variable is a 1D numpy array and not 'time' or 'step' variables, add it to variables dictionary
            if isinstance(v, np.ndarray) and v.ndim == 1 and k not in ('step', 'time'):
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