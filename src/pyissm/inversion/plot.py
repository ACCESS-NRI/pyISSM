import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def plot_convergence_history(convergence_history,
                             metrics = None,
                             log_scale = True,
                             ax = None,
                             figsize = (6.4, 4.8),
                             constrained_layout = True):
    """
    Plot optimization convergence metrics over iterations.

    Parameters
    ----------
    convergence_history : :class:`pandas.DataFrame`
        Table containing an ``"iteration"`` column and one or more metric
        columns to plot against iteration.

    metrics : :class:`list` of :class:`str` or ``None``, optional
        Metric column names to plot. If ``None`` (default), all columns except
        ``"iteration"`` are plotted.

    log_scale : :class:`bool`, optional
        If ``True`` (default), use a logarithmic scale for the y-axis.

    ax : matplotlib.axes.Axes, optional
        Axes object to draw on. If None, a new figure and axes are created.

    figsize : tuple of float, optional
        Size of the figure in inches (width, height). Default is (6.4, 4.8).

    constrained_layout : bool, optional
        Whether to use constrained layout in figure. Default is True.

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
        The created Matplotlib figure. Returned only if a new figure was created.

    ax : :class:`matplotlib.axes.Axes`
        The axes containing the heatmap. If an axes was passed in, this is the same object; if no axes was passed, this is the newly created axes.

    Notes
    -----
        This function does not validate metric names. A missing metric column will
        raise a pandas ``KeyError`` during plotting.
    """

    ## Is an ax passed?
    ax_defined = ax is not None

    ## If no ax is defined, create new figure
    ## otherwise, plot on defined ax
    if ax is None:
        fig, ax = plt.subplots(figsize = figsize, constrained_layout = constrained_layout)
    else:
        fig = ax.get_figure()

    # If metrics is not provided, plot all metrics in convergence history (except iteration)
    if metrics is None:

        metrics = [
            col for col in convergence_history.columns
            if col != "iteration"
            ]
    
    # Loop over metrics and plot each one
    for metric in metrics:
        ax.plot(convergence_history["iteration"],
                convergence_history[metric],
                label = metric)
    
    # Set log scale if requested
    if log_scale:
        ax.set_yscale('log')
    
    # Set labels and legend
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Metric Value")
    ax.legend(title = "Metric")
    
    if not ax_defined:
        return fig, ax
    else:
        return ax

def plot_sensitivity_heatmap(diagnostics,
                             x,
                             y,
                             value,
                             cmap = "viridis",
                             ax = None,
                             figsize = (6.4, 4.8),
                             constrained_layout = True,
                             **kwargs):
    """
    Plot a sensitivity heatmap from tabular diagnostics data.

    Parameters
    ----------
    diagnostics : :class:`pandas.DataFrame`
        Input table containing at least the columns specified by ``x``, ``y``,
        and ``value``.

    x : :class:`str`
        Column name used for heatmap x-axis categories.

    y : :class:`str`
        Column name used for heatmap y-axis categories.

    value : :class:`str`
        Column name whose values populate the heatmap cells.

    cmap : :class:`str` or :class:`matplotlib.colors.Colormap`, optional
        Colormap used to render the heatmap. Defaults to ``"viridis"``.

    ax : matplotlib.axes.Axes, optional
        Axes object to draw on. If None, a new figure and axes are created.

    figsize : :class:`tuple` of :class:`float`, optional
        Figure size in inches as ``(width, height)``. Defaults to
        ``(6.4, 4.8)``.

    constrained_layout : :class:`bool`, optional
        Whether to use Matplotlib constrained layout. Defaults to ``True``.

    **kwargs
        Additional keyword arguments forwarded to :meth:`matplotlib.axes.Axes.imshow`.

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
        The created Matplotlib figure. Returned only if a new figure was created.

    ax : :class:`matplotlib.axes.Axes`
        The axes containing the heatmap. If an axes was passed in, this is the same object; if no axes was passed, this is the newly created axes.

    Notes
    -----
        The heatmap is generated via :meth:`pandas.DataFrame.pivot_table` using
        ``y`` as the index, ``x`` as the columns, and ``value`` as cell values.
    """

    ## Is an ax passed?
    ax_defined = ax is not None

    ## If no ax is defined, create new figure
    ## otherwise, plot on defined ax
    if ax is None:
        fig, ax = plt.subplots(figsize = figsize, constrained_layout = constrained_layout)
    else:
        fig = ax.get_figure()

    # Pivot table
    heatmap = diagnostics.pivot_table(index = y, columns = x, values = value)

    # Populate heatmap plot
    im = ax.imshow(
        heatmap.values,
        aspect = "auto",
        origin = "lower",
        cmap = cmap,
        **kwargs
    )

    # Tick labels
    ax.set_xticks(range(len(heatmap.columns)))
    ax.set_xticklabels(heatmap.columns)

    ax.set_yticks(range(len(heatmap.index)))
    ax.set_yticklabels(heatmap.index)

    # Labels and title
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title(value)

    # Add colorbar
    plt.colorbar(im, ax = ax, label=value)

    if not ax_defined:
        return fig, ax
    else:
        return ax