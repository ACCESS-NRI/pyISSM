import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from . import metrics
from .. import tools, plot
from matplotlib.backends.backend_pdf import PdfPages

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

def plot_run_summary(md,
                     output_pdf,
                     diagnostics = None,
                     plot_kwargs = None):
    """
    Generate PDF summary report for inversion run.

    Saves a multi-page PDF report to the specified path, containing visualizations and diagnostics of the inversion run.
    The first page includes spatial plots of observed velocity, modelled velocity, velocity residuals, and the inverted
    control parameter field. The second page includes a convergence history plot, histograms of velocity residuals and
    inverted field slopes, and a text summary of key diagnostics. Plot appearance can be customized via the `plot_kwargs` parameter.

    Parameters
    ----------
    md : :class:`pyissm.Model`
        ISSM model object containing inversion results.
    output_pdf : :class:`str`
        Path to output PDF file.
    diagnostics : :class:`dict`, optional
        Diagnostics configuration. Default is None.
    plot_kwargs : :class:`dict`, optional
        Plot keyword arguments. Default is None.

    Returns
    -------
    None
    """

    # Compute velocity residuals for diagnostics and plotting
    residuals = tools.diagnostics.velocity_residuals(md)
    
    # Determine if velocity residuals exceed colorbar limits, and set colorbar extension accordingly
    extend_min = False
    extend_max = False
    if np.max(residuals) > 50:
        extend_max = True
        extend = 'max'
    if np.min(residuals) < -50:
        extend_min = True
        extend = 'min'
    if extend_max and extend_min:
        extend = 'both'
    else:
        extend = 'neither'

    # Extract control parameter field and name for plotting
    inversion_field, inversion_field_name = metrics._find_control_parameter_field(md)

    # Define default plot kwargs for each plot type, which can be overridden by user-provided plot_kwargs
    default_plot_kwargs = {

        # PAGE 1 - SPATIAL FIELDS
        ## Observed velocity
        "observed": {
            "show_cbar": True,
            "xlabel": "",
            "cbar_kwargs": {
                "label": "Observed velocity (m/yr)",
            },
        },
        ## Modelled velocity
        "modelled": {
            "show_cbar": True,
            "xlabel": "",
            "ylabel": "",
            "cbar_kwargs": {
                "label": "Modelled velocity (m/yr)",
            },
        },
        ## Velocity residuals
        "residuals": {
            "cmap": "RdBu_r",
            "show_cbar": True,
            "vmin": -50,
            "vmax": 50,
            "cbar_kwargs": {
                "label": "Velocity residual (m/yr)",
                "extend": extend,
            },
        },
        ## Inversion field
        "inversion": {
            "show_cbar": True,
            "ylabel": "",
            "cbar_kwargs": {
                "label": f"Inverted {inversion_field_name}",
            },
        },

        # PAGE 2 - CONVERGENCE, DISTRIBUTIONS & TEXT SUMMARY
        ## Convergence history
        "convergence": {
            "log_scale": True,
        },
        ## Residual histogram
        "residual_hist": {
            "bins": 100,
        },
        ## Inversion field slope histogram
        "slope_hist": {
            "bins": 100,
        },

        # FIGURE SETTINGS (PER PAGE)
        ## Page 1 - spatial fields
        "page1_figure": {
            "figsize": (15, 10),
            "sharex": True,
            "sharey": True,
            "constrained_layout": True,
        },
        ## Page 2 - convergence, distributions & text summary
        "page2_figure": {
            "figsize": (15, 10),
            "constrained_layout": True,
        }
    }

    # Merge user-provided plot kwargs with defaults
    plot_kwargs = {} if plot_kwargs is None else plot_kwargs
    
    merged_plot_kwargs = {}

    for key, defaults in default_plot_kwargs.items():
        merged_plot_kwargs[key] = {
            **defaults,
            **plot_kwargs.get(key, {}),
        }

    # Generate PDF report with two pages: (1) spatial fields and (2) convergence, distributions, and text summary
    with PdfPages(output_pdf) as pdf:

        # --------------------------------------------------
        # PAGE 1 — SPATIAL FIELDS
        # --------------------------------------------------

        # Set up 2x2 figure for spatial plots
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, **merged_plot_kwargs["page1_figure"])

        # Observed velocity
        plot.plot_model_field(
            md,
            md.inversion.vel_obs,
            ax = ax1,
            **merged_plot_kwargs["observed"]
        )
        ax1.set_title("Observed velocity")

        # Modelled velocity
        plot.plot_model_field(
            md,
            md.results.StressbalanceSolution.Vel,
            ax = ax2,
            **merged_plot_kwargs["modelled"]
        )
        ax2.set_title("Modelled velocity")

        # Residuals
        plot.plot_model_field(
            md,
            residuals,
            ax = ax3,
            **merged_plot_kwargs["residuals"]
        )
        ax3.set_title("Velocity residual (mod - obs)")

        # Inversion field
        plot.plot_model_field(
            md,
            inversion_field,
            ax = ax4,
            **merged_plot_kwargs["inversion"]
        )
        ax4.set_title(f"Inverted {inversion_field_name}")

        # Save page 1 to PDF
        pdf.savefig(fig)
        plt.close(fig)

        # --------------------------------------------------
        # PAGE 2 — CONVERGENCE, DISTRIBUTIONS & TEXT SUMMARY
        # --------------------------------------------------
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, **merged_plot_kwargs["page2_figure"])

        # Convergence history
        convergence_history = metrics.extract_convergence_history(md)

        plot_convergence_history(convergence_history,
                                 ax = ax1,
                                 **merged_plot_kwargs["convergence"]
                                 )
        ax1.set_title("Cost function convergence history")
        ax1.grid()

        # Residual histogram
        finite = residuals[np.isfinite(residuals)]

        ax2.hist(finite,
                 **merged_plot_kwargs["residual_hist"]
                 )

        ax2.set_title("Velocity residual distribution")
        ax2.set_xlabel("Velocity residual (m/yr)")
        ax2.set_ylabel("Count")

        # Slope histogram
        _, _, slope = tools.geometry.slope(md, inversion_field)

        finite = slope[np.isfinite(slope)]

        ax3.hist(finite,
                 **merged_plot_kwargs["slope_hist"]
                 )

        ax3.set_title(f"Inverted {inversion_field_name} distribution")
        ax3.set_xlabel(f"{inversion_field_name} slope")
        ax3.set_ylabel("Count")

        # Text Summary
        ## If diagnostics are not provided, compute a default set of diagnostics to display
        if diagnostics is None:
            diagnostics = {
                "vel_rmse": tools.diagnostics.velocity_rmse(md),

                # Expand dictionaries into flat metrics
                **metrics.cost_function_values(md),
                **metrics.cost_function_ratios(md),
                **metrics.velocity_residual_area_metrics(md),
                **metrics.field_smoothness_metrics(md)
            }
        ## Turn off axis
        ax4.axis("off")

        ## Format diagnostics into lines of text for display
        lines = []
        for key, value in diagnostics.items():

            if isinstance(value, (np.integer, int)):
                lines.append(f"{key}: {value}")

            elif isinstance(value, (np.floating, float)):
                lines.append(f"{key}: {value:.3g}")

            else:
                lines.append(f"{key}: {value}")
        
        ## Display diagnostics as text in the axis
        ax4.text(
            0.01,
            0.99,
            "\n".join(lines),
            va = "top",
            family = "monospace",
        )

        # Save page 2 to PDF
        pdf.savefig(fig)
        plt.close(fig)