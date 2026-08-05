"""
Sensitivity analysis tools for ISSM inversions.

This module provides functions to systematically vary inversion parameters, execute runs, and compute diagnostics for analyzing sensitivity of inversion results to parameter choices.
"""

import copy
import itertools
from pathlib import Path
import numpy as np
import pandas as pd
from .. import model, tools
from . import metrics, plot

def build_parameter_grid(coefficients, initial_run_id = 1):
    """
    Build parameter sensitivity combinations.

    Parameters
    ----------
    coefficients : :class:`dict`
        Dictionary mapping ISSM cost function IDs to coefficient lists.

        Example
        -------
        {101: [0.1, 1, 10],
        103: [1, 10],
        502: [1e-19]}

    initial_run_id : :class:`int`, default = 1
        Starting index for run_id assignment. Subsequent runs will have IDs incremented by 1 from this value.
        Useful for tracking runs across multiple calls to this function without overwriting previous run_ids.

    Returns
    -------
    :class:`pandas.DataFrame`
        Parameter combination table.

    Notes
    -----
    Returned DataFrame contains:
        - run_id
        - cfXXX columns for each cost function
        - run_name

    Example
    -------
    : code-block:: python

        >>> coefficients = {101: [1, 10, 100, 1000],
        ...                 103: [1, 10, 100, 1000]}
        >>> param_grid = pyissm.inversion.sensitivity.build_parameter_grid(coefficients)
        >>> param_grid = pyissm.inversion.sensitivity.build_parameter_grid(coefficients, initial_run_id = 101)
    """

    # Sort cost functions for consistent ordering
    cost_functions = sorted(coefficients.keys())

    # Build cartesian product of parameter combinations
    combinations = list(
        itertools.product(
            *[coefficients[cf] for cf in cost_functions]
        )
    )

    # Store output rows before converting to DataFrame
    records = []

    # Loop through all parameter combinations
    for run_id, combo in enumerate(combinations, start = initial_run_id):

        # Initialize record with run_id
        record = {
            "run_id": run_id,
        }

        # Add coefficient columns
        for cf, value in zip(cost_functions, combo):
            record[f"cf{cf}"] = value

        # Generate unique run name based on parameter values
        pair_string = "_".join(
            [f"{v:g}" for v in combo]
        )

        # Assign run_name to record
        record["run_name"] = f"run_{run_id:03d}_{pair_string}"

        # Append record to list
        records.append(record)

    # Convert list of records to DataFrame
    return pd.DataFrame(records)


def assign_cost_functions(md,
                          coefficients,
                          global_mask = None,
                          coeff_masks = None):
    """
    Assign inversion cost function coefficients.

    Parameters
    ----------
    md : :class:`pyissm.Model`
        ISSM model object.

    coefficients : :class:`dict`
        Dictionary mapping cost function IDs to scalar coefficients.

        Example
        -------
        {101: 1000,
        103: 25,
        502: 1e-19}

    global_mask : :class:`numpy.ndarray` of :class:`bool`, optional
        Global boolean mask applied to all cost functions.
            - True  -> coefficient active
            - False -> coefficient set to zero

    coeff_masks : :class:`dict`, optional
        Dictionary mapping cost-function IDs to boolean masks.
        Masks are applied in the same way as the global_mask, but only affect their respective cost functions.

        Example
        -------
        {101: vel_mask,
        103: grounded_mask}

    Returns
    -------
    md : :class:`pyissm.Model`
        Updated ISSM model.

    Example
    -------
    : code-block:: python

        >>> coefficients = {101: 1000,
        ...                 103: 25,
        ...                 502: 1e-19}
        >>> global_mask = md.inversion.vel_obs > 0
        >>> md = pyissm.inversion.sensitivity.assign_cost_functions(md, coefficients, global_mask = global_mask)

        >>> coeff_masks = {101: md.inversion.vel_obs > 0,
        ...                103: md.mask.ice_levelset > 0}
        >>> md = pyissm.inversion.sensitivity.assign_cost_functions(md, coefficients, coeff_masks = coeff_masks)
    """

    # Sort cost functions for consistent ordering
    cost_functions = sorted(coefficients.keys())

    # Determine number of vertices
    nverts = md.mesh.numberofvertices

    # Allocate coefficient matrix (default = 0 [inactive])
    coeff_matrix = np.zeros((nverts, len(cost_functions)), dtype=float)

    # Validate global_mask
    if global_mask is None:
        global_mask = np.ones(nverts, dtype = bool)
    else:

        # Ensure global_mask is a boolean array of the correct shape. Enforce boolean to prevent potential unintentional partial weighting with floats
        global_mask = np.asarray(global_mask)

        if global_mask.dtype != bool:
            raise TypeError("pyissm.inversion.sensitivity.assign_cost_functions: global_mask must be boolean.")

        if global_mask.shape != (nverts,):
            raise ValueError(f"pyissm.inversion.sensitivity.assign_cost_functions: global_mask must have shape ({nverts},)")
    

    # Validate coeff_masks
    if coeff_masks is not None:

        # Ensure coeff_masks is a dictionary mapping cost function IDs to boolean masks
        for cf, cf_mask in coeff_masks.items():

            cf_mask = np.asarray(cf_mask)

            # Ensure cf_mask is a boolean array of the correct shape. Enforce boolean to prevent potential unintentional partial weighting with floats
            if cf_mask.dtype != bool:
                raise TypeError(f"pyissm.inversion.sensitivity.assign_cost_functions: Mask for cost function {cf} must be boolean.")

            if cf_mask.shape != (nverts,):
                raise ValueError(f"pyissm.inversion.sensitivity.assign_cost_functions: Mask for cost function {cf} must have shape ({nverts},)")

            coeff_masks[cf] = cf_mask

    # Populate coefficient matrix
    for i, cf in enumerate(cost_functions):

        # Start from global mask
        active = global_mask.copy()

        # Apply per-cost-function mask
        if coeff_masks is not None and cf in coeff_masks:
            active &= coeff_masks[cf]

        # Assign coefficient only where active
        coeff_matrix[active, i] = coefficients[cf]

    # Assign to ISSM model
    md.inversion.cost_functions = cost_functions
    md.inversion.cost_functions_coefficients = coeff_matrix

    return md

def parameter_sensitivity(md,
                          parameter_grid,
                          output_dir,
                          run = True,
                          load_only = False,
                          global_mask = None,
                          coeff_masks = None):
    """
    Run inversion sensitivity experiments.

    Parameters
    ----------
    md : :class:`pyissm.Model`
        ISSM model object.

    parameter_grid : :class:`pandas.DataFrame`
        Output from :func:`pyissm.inversion.sensitivity.build_parameter_grid`.

    output_dir : :class:`str` or :class:`pathlib.Path`
        Directory used to store sensitivity outputs, manifests, and saved models.
        
        **NOTE:**
        This is intentionally separate from: md.cluster.executionpath:
            - `md.cluster.executionpath` controls where ISSM stages runtime files for execution.
            - `output_dir` controls where pyISSM stores organized experiment outputs.

    run : :class:`bool`, default = True
        Submit inversion runs.

    load_only : :class:`bool`, default = False
        Load completed inversion results.

    global_mask : :class:`numpy.ndarray` of :class:`bool`, optional
        Global boolean mask applied to all cost functions. Passed to :func:`assign_cost_functions`.
            - True  -> coefficient active
            - False -> coefficient set to zero

    coeff_masks : :class:`dict`, optional
        Dictionary mapping cost-function IDs to boolean masks. Passed to :func:`assign_cost_functions`.

        Example
        -------
        {101: vel_mask,
        103: grounded_mask}
    
    Returns
    -------
    :class:`pandas.DataFrame`
        Experiment manifest.

    Example
    -------
    : code-block:: python

        >>> coefficients = {101: [1, 10, 100, 1000],
        ...                 103: [1, 10, 100, 1000]}
        >>> parameter_grid = pyissm.inversion.sensitivity.build_parameter_grid(coefficients)
        >>> manifest = pyissm.inversion.sensitivity.parameter_sensitivity(md, parameter_grid, output_dir = './sensitivity_runs', run = True, load_only = False)
    """

    # Prevent ambiguous execution mode
    if run and load_only:
        raise ValueError("pyissm.inversion.sensitivity.parameter_sensitivity: Cannot set both run and load_only to True.")

    if not run and not load_only:
        raise ValueError("pyissm.inversion.sensitivity.parameter_sensitivity: Either run or load_only must be True.")

    # Create output directory for sensitivity study
    output_dir = Path(output_dir)
    output_dir.mkdir(parents = True, exist_ok = True)

    # Initialize manifest records
    manifest = []

    # Identify cost-function columns (e.g., cf101, cf102, cf502)
    cf_columns = [
        col for col in parameter_grid.columns
        if col.startswith("cf")
        ]

    # Loop through all parameter combinations
    for _, row in parameter_grid.iterrows():

        # Extract run metadata
        run_id = int(row["run_id"])
        run_name = row["run_name"]

        print(f"Running sensitivity experiment: {run_name}")

        # Create independent model copy
        mdi = copy.deepcopy(md)

        # Build coefficient dictionary
        coefficients = {}

        # Extract coefficient values from the parameter grid row
        for col in cf_columns:

            # Extract cost function ID from column name (e.g., cf101 -> 101)
            cf = int(col.replace("cf", ""))

            # Store coefficient value for this cost function
            coefficients[cf] = row[col]

        # Assign inversion coefficients
        mdi = assign_cost_functions(mdi, coefficients, global_mask = global_mask, coeff_masks = coeff_masks)

        # Assign unique ISSM runtime name
        mdi.miscellaneous.name = run_name

        # Create dedicated run directory to store saved model and convergence history for this run
        run_dir = output_dir / run_name
        run_dir.mkdir(exist_ok = True)

        # Default status
        status = "pending"

        # Try to execute or load the run, and catch any exceptions to record failure in the manifest
        try:

            # Submit inversion run
            if run:
                
                mdi = model.execute.solve(mdi, "Stressbalance", runtime_name = False, load_only = False)

                # Distinguish between submitted vs completed based on waitonlock setting
                if mdi.settings.waitonlock == 0:
                    # If waitonlock is 0 the run is submitted, but we do not wait for completion, so we mark as "submitted"
                    status = "submitted"
                else:
                    # If waitonlock is not 0, we wait for the run to complete before proceeding, so we can mark as "completed"
                    status = "completed"
                    
                    # Save completed model
                    model.io.save_model(mdi, run_dir / f"{run_name}.nc")

            # Load previously completed results
            elif load_only:
                
                # If the model file already exists, we can assume the run was completed and just load the model
                if (run_dir / f"{run_name}.nc").exists():
                    print(f"Model already exists for {run_name}, loading model.")
                    status = "completed"
                else:
                    print(f"Model does not exist for {run_name}, attempting to load model from execution path.")

                    # Load the model from the execution path (this will raise an error if the model is not found, which we catch to mark the run as failed)
                    mdi = model.execute.solve(mdi, "Stressbalance", runtime_name = False, load_only = True)

                    # Save loaded model
                    model.io.save_model(mdi, run_dir / f"{run_name}.nc")

                    status = "completed"

        # If exceptions occur during run submission or loading, catch them and record failure in the manifest
        except Exception as exc:

            print(f"FAILED: {run_name}")
            print(exc)

            status = "failed"

        # Build manifest entry
        manifest_entry = {
            "run_id": run_id,
            "run_name": run_name,
            "run_dir": str(run_dir),
            "execution_path": str(mdi.cluster.executionpath),
            "run_status": status,
        }

        # Store coefficient metadata (e.g. cf101, cf103) in manifest for easy reference when analyzing diagnostics later
        manifest_entry.update({f"cf{cf}": value for cf, value in coefficients.items()})

        # Append manifest entry
        manifest.append(manifest_entry)

    # Convert manifest to DataFrame
    manifest = pd.DataFrame(manifest)

    # Save manifest
    manifest.to_csv(output_dir / "manifest.csv", index = False)

    return manifest

def load_parameter_sensitivity_run(manifest_row):
    """
    Load a parameter sensitivity run from a manifest row.

    Parameters
    ----------
    manifest_row : :class:`pandas.Series`

    Returns
    -------
    md : :class:`pyissm.Model`
        Loaded ISSM model for the specified run.

    Example
    -------
    : code-block:: python

        >>> manifest = pd.read_csv('./sensitivity_runs/manifest.csv')
        >>> md = pyissm.inversion.sensitivity.load_parameter_sensitivity_run(manifest.iloc[0])
    """

    # Extract run directory and name from manifest row
    run_dir = Path(manifest_row["run_dir"])
    run_name = manifest_row["run_name"]

    # Construct model path
    model_path = run_dir / f"{run_name}.nc"

    # Check if model file exists
    if not model_path.exists():
        raise FileNotFoundError(f"pyissm.inversion.sensitivity.load_parameter_sensitivity_run: Model file not found: {model_path}")

    # Load model
    md = model.io.load_model(model_path)

    return md

def compute_sensitivity_diagnostics(manifest,
                                    output_dir = None,
                                    plot_run_summaries = True):
    """
    Compute diagnostics for parameter sensitivity runs from a manifest.

    Parameters
    ----------
    manifest : :class:`pandas.DataFrame`
        Output from :func:`pyissm.inversion.sensitivity.parameter_sensitivity`.

    output_dir : :class:`str` or :class:`pathlib.Path`, optional
        Directory used to store computed diagnostics. If None, diagnostics are not saved to disk.

    plot_run_summaries : :class:`bool`, default = True
        Whether to generate summary PDF reports for each run, including diagnostic plots and convergence history. If True, summary PDFs are saved in each run's directory.

    Returns
    -------
    :class:`pandas.DataFrame`
        Diagnostics table.

    Example
    -------
    : code-block:: python

        >>> manifest = pd.read_csv('./sensitivity_runs/manifest.csv')
        >>> diagnostics = pyissm.inversion.sensitivity.compute_sensitivity_diagnostics(manifest, output_dir = './sensitivity_runs/diagnostics')
    
    """

    # Copy manifest so coefficient metadata is preserved
    diagnostics = manifest.copy()

    # Loop through all runs
    for idx, row in manifest.iterrows():

        # Extract run name for logging
        run_name = row["run_name"]
        print(f"Processing diagnostics: {run_name}")

        # Try to load the model and compute diagnostics, and catch any exceptions to record failure in the diagnostics table
        try:

            # Load model
            mdi = load_parameter_sensitivity_run(row)

            # -----------------------------------------
            # Compute diagnostics
            # -----------------------------------------

            # Compute all diagnostics for this run and store in a dictionary
            run_metrics = {
                "vel_rmse": tools.diagnostics.velocity_rmse(mdi),

                # Expand dictionaries into flat metrics
                **metrics.cost_function_values(mdi),
                **metrics.cost_function_ratios(mdi),
                **metrics.velocity_residual_area_metrics(mdi),
                **metrics.field_smoothness_metrics(mdi)
            }

            # Dynamically assign all metrics to the diagnostics DataFrame for this run
            for key, value in run_metrics.items():

                diagnostics.at[idx, key] = value

            diagnostics.at[idx, "diagnostics_status"] = "complete"

            if output_dir is not None:
                # Save intermediate diagnostics to CSV after each run in case of long-running computations and to preserve results even if some runs fail
                diagnostics.to_csv(Path(output_dir) / "diagnostics.csv", index = False)

            # -----------------------------------------
            # Extract & save convergence history
            # -----------------------------------------
            
            # Extract convergence history DataFrame from model results
            convergence_history = metrics.extract_convergence_history(mdi)
    
            # Save convergence history to CSV in the run directory
            convergence_history.to_csv(Path(row["run_dir"]) / f"{run_name}_convergence_history.csv", index = False)
    
            # -----------------------------------------
            # Generate summary PDF for this run
            # -----------------------------------------
            if plot_run_summaries:
                plot.plot_run_summary(mdi,
                                      output_pdf = Path(row["run_dir"]) / f"{run_name}_summary.pdf",
                                      diagnostics = diagnostics.loc[idx])

        # If exceptions occur during loading or diagnostics computation, catch them and record failure in the diagnostics table
        except Exception as exc:

            # Record failure in diagnostics table & print error message
            diagnostics.at[idx, "diagnostics_status"] = "failed"
            print(f"Error processing run {run_name}: {exc}")

    return diagnostics

def normalize_diagnostics(diagnostics,
                          columns = None,
                          higher_is_better = None,
                          suffix = "_norm"):
    """
    Normalize diagnostic columns to [0, 1].

    Parameters
    ----------
    diagnostics : :class:`pandas.DataFrame`
        DataFrame containing model diagnostics.

    columns : :class:`list` of :class:`str`, optional
        Columns to normalize. Defaults to all numeric columns.

    higher_is_better : :class:`dict`, optional
        Dictionary mapping column name -> bool.
        - True  = higher values are better
        - False = lower values are better

        Example
        -------
        {"velocity_rmse": False,
        "thickness_rmse": False,
        "correlation": True}

        Any unspecified column defaults to False.

    suffix : :class:`str`, optional
        Suffix for normalized columns.

    Returns
    -------
    :class:`pandas.DataFrame`
        Copy of dataframe with normalized columns added.

    Example
    -------
    : code-block:: python
    
        >>> manifest = pd.read_csv('./sensitivity_runs/manifest.csv')
        >>> diagnostics = pyissm.inversion.sensitivity.compute_sensitivity_diagnostics(manifest, output_dir = './sensitivity_runs/diagnostics')
        >>> diagnostics_norm = pyissm.inversion.sensitivity.normalize_diagnostics(diagnostics, columns = ['vel_rmse', 'thickness_rmse'])
    """

    # Create a copy of the diagnostics DataFrame to avoid modifying the original
    out = diagnostics.copy()

    # If columns to normalize are not specified, default to all numeric columns
    if columns is None:
        columns = out.select_dtypes(include=[np.number]).columns.tolist()

    # If higher_is_better is not provided, default to all False (lower is better)
    if higher_is_better is None:
        higher_is_better = {}

    # Loop through specified columns and normalize each one
    for col in columns:

        # Extract values and compute min/max for normalization
        values = out[col].astype(float)

        vmin = values.min()
        vmax = values.max()

        # Avoid divide-by-zero if all values identical
        if np.isclose(vmin, vmax):
            out[f"{col}{suffix}"] = 1.0
            continue
        
        # Normalize based on whether higher values are better or not
        if higher_is_better.get(col, False):
            # Higher is better
            norm = (values - vmin) / (vmax - vmin)

        else:
            # Lower is better
            norm = (vmax - values) / (vmax - vmin)

        # Assign normalized values to new column in output DataFrame
        out[f"{col}{suffix}"] = norm

    return out