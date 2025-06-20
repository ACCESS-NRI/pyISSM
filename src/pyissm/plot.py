"""
Functions to visualize ISSM Models
"""

import matplotlib.pyplot as plt
import matplotlib
from .model.mesh import process_mesh, find_node_types

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
                node_args = {},
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
    node_args : dict, optional
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
    ax2 = plot_mesh2d(mesh, ax = ax2, show_nodes = True, node_args = {'color': 'red'})
    """

    ## Is an ax passed?
    ax_defined = ax is not None

    ## Define default node argus and update node_args with passed arguments
    default_node_args = {'marker': '.',
                         'color': 'k',
                         's': 5}
    default_node_args.update(**node_args)


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
               **kwargs)

    ## Add optional nodes
    if show_nodes:
        ax.scatter(mesh.x, mesh.y, **default_node_args)

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
                     figsize = (6.4, 4.8),
                     constrained_layout = True,
                     show_mesh = True,
                     mesh_args = {},
                     show_legend = True,
                     legend_args = {}):

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
    mesh_args : dict, optional
        Additional keyword arguments passed to plot_mesh2d() for customizing mesh appearance.
    show_legend : bool, optional
        Whether to display a legend identifying node types. Default is True.
    legend_args : dict, optional
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
    fig, ax = plot_model_nodes(md, md.mask.ice_levelset, md.mask.ocean_levelset, type = 'ice_front_nodes', mesh_args = {'linewidth': 0.5, 'color': 'cyan'})
    fig, ax = plot_model_nodes(md, md.mask.ice_levelset, md.mask.ocean_levelset, type = 'floating_ice_nodes', colors = ['blue', 'magenta'])
    """

    ## Set defaults
    ax_defined = ax is not None

    default_mesh_args = {'alpha': 0.5}
    default_mesh_args.update(**mesh_args)

    default_legend_args = {'title': 'Node types',
                           'fontsize': 10}
    default_legend_args.update(**legend_args)

    ## Define binary colormap
    cmap = matplotlib.colors.ListedColormap(colors)

    ## Process model mesh
    mesh, mesh_x, mesh_y, mesh_elements, is3d = process_mesh(md)

    ## Find node types
    node_types = find_node_types(md,
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
                    **default_mesh_args)

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
        ax.legend(handles = legend_elements, **default_legend_args)

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
                        figsize = (6.4, 4.8),
                        constrained_layout = True,
                        show_mesh = True,
                        mesh_args = {},
                        show_legend = True,
                        legend_args={}):

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
    mesh_args : dict, optional
        Additional keyword arguments passed to plot_mesh2d() for customizing mesh appearance.
    show_legend : bool, optional
        Whether to display a legend identifying node types. Default is True.
    legend_args : dict, optional
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
    fig, ax = plot_model_elements(md, md.mask.ice_levelset, md.mask.ocean_levelset, type = 'ice_front_elements', mesh_args = {'linewidth': 0.5, 'color': 'cyan'})
    fig, ax = plot_model_elements(md, md.mask.ice_levelset, md.mask.ocean_levelset, type = 'floating_ice_elements', color = 'red')
    """

    ## Set defaults
    ax_defined = ax is not None

    default_mesh_args = {'alpha': 0.5}
    default_mesh_args.update(**mesh_args)

    default_legend_args = {'title': 'Element type',
                           'fontsize': 10,
                           'loc': 'upper right'}
    default_legend_args.update(**legend_args)

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
        raise ValueError(f'No {type} elements exist in the model.')

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

    ## Plot elements
    ax.tripcolor(mesh_x, mesh_y, mesh_elements[element_pos], facecolors=colors, cmap = cmap, edgecolors = 'none')

    ## Add mesh (optional) with specific arguments
    if show_mesh:
        plot_mesh2d(mesh, ax = ax, **default_mesh_args)

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
        ax.legend(handles = legend_elements, **default_legend_args)

    ## Return
    if not ax_defined:
        return fig, ax
    else:
        return ax