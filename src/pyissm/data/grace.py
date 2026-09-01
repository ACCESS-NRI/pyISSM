"""
GRACE liquid-water-equivalent thickness loading and mesh interpolation.

This module loads a GRACE mascon NetCDF product and interpolates its
liquid-water-equivalent thickness onto an ISSM model mesh, replacing the
upstream MATLAB tutorial's `examples/Functions/grace.m` (hand-rolled
longitude-seam gap-filling, hardcoded fill-value thresholds, `griddata`).
"""

import numpy as np
import pandas as pd
import xarray as xr

from pyissm.data.interp import xr_to_mesh

# Name of the liquid-water-equivalent thickness variable in the GRACE mascon
# product this module targets (units: cm).
LWE_THICKNESS_VAR = 'lwe_thickness'

# ---------------------------------------------------------------------------
# Dataset loading
# ---------------------------------------------------------------------------

def _open_grace_dataset(filename):

    """
    Open a GRACE NetCDF product with CF decoding applied.

    Relies on xarray's native CF handling (`mask_and_scale=True`,
    `decode_times=True`) so this product's `_FillValue` sentinel becomes NaN
    and the `time` axis decodes to real datetimes automatically, rather than
    hardcoding the MATLAB tutorial's fill-value thresholds or reimplementing
    its manual date arithmetic.

    Parameters
    ----------
    filename : str
        Path to the GRACE NetCDF file.

    Returns
    -------
    xr.Dataset
        CF-decoded dataset, with `lwe_thickness` guaranteed to be ordered
        `(time, lat, lon, ...)`.

    Raises
    ------
    ValueError
        If the dataset has no `lwe_thickness` variable, or that variable is
        missing a `time`, `lat`, or `lon` dimension.
    """

    ds = xr.open_dataset(filename, mask_and_scale = True, decode_times = True)

    if LWE_THICKNESS_VAR not in ds:
        ds.close()
        raise ValueError(f"pyissm.data.grace._open_grace_dataset: dataset is missing '{LWE_THICKNESS_VAR}'")

    missing_dims = {'time', 'lat', 'lon'} - set(ds[LWE_THICKNESS_VAR].dims)
    if missing_dims:
        ds.close()
        raise ValueError(f"pyissm.data.grace._open_grace_dataset: "
                          f"'{LWE_THICKNESS_VAR}' is missing dimension(s) {sorted(missing_dims)}")

    # Confirmed against the real product (GRCTellus.JPL.200204_201701.LND.
    # RL05_1.DSTvSCS1411.nc, 2026-09-01): dims already arrive as
    # (time, lat, lon). Transpose defensively rather than assume every GRACE
    # product orders them the same way.
    return ds.transpose('time', 'lat', 'lon', ...)

# ---------------------------------------------------------------------------
# Epoch selection
# ---------------------------------------------------------------------------

def _to_decimal_year(times):

    """
    Convert datetimes to decimal years, accounting for year length.

    Each timestamp's fractional part is its true elapsed proportion of its
    own calendar year (365 or 366 days), computed from the actual duration
    between that year's start and the next — not a fixed `days/365` divisor,
    which silently mis-converts every date in a leap year.

    Parameters
    ----------
    times : array-like of datetime64
        Decoded timestamps (e.g. `ds['time'].values`).

    Returns
    -------
    ndarray of float
        Decimal-year value for each input timestamp.
    """

    times = pd.DatetimeIndex(times)
    year_start = pd.to_datetime({'year': times.year, 'month': 1, 'day': 1})
    year_end = pd.to_datetime({'year': times.year + 1, 'month': 1, 'day': 1})
    fraction_of_year = (times - year_start) / (year_end - year_start)

    return times.year.to_numpy(dtype = float) + fraction_of_year.to_numpy()

def _select_epoch_indices(decimal_year, tmin, tmax):

    """
    Select the contiguous run of time indices nearest a decimal-year range.

    Parameters
    ----------
    decimal_year : ndarray of float
        Decimal year for each time index, ascending (as produced by
        `_to_decimal_year` on a chronologically-ordered `time` axis).
    tmin, tmax : float
        Decimal-year bounds. The nearest available epoch to each bound is
        selected; `tmin == tmax` selects a single epoch.

    Returns
    -------
    ndarray of int
        Indices into `decimal_year`, inclusive of both endpoints.

    Raises
    ------
    ValueError
        If the nearest index to `tmax` precedes the nearest index to `tmin`
        (i.e. `tmax < tmin`).
    """

    decimal_year = np.asarray(decimal_year, dtype = float)
    idx_min = int(np.argmin(np.abs(decimal_year - tmin)))
    idx_max = int(np.argmin(np.abs(decimal_year - tmax)))

    if idx_max < idx_min:
        raise ValueError(f"pyissm.data.grace._select_epoch_indices: tmax ({tmax}) resolves to an "
                          f"earlier epoch than tmin ({tmin})")

    return np.arange(idx_min, idx_max + 1)

# ---------------------------------------------------------------------------
# Longitude-seam handling
# ---------------------------------------------------------------------------

def _to_longitude_0_360(long):

    """
    Convert longitudes to the GRACE product's `[0, 360)` convention.

    The mesh's `long` is `[-180, 180]`; rather than reprojecting the grid,
    the query points are converted to match it.

    Parameters
    ----------
    long : array-like
        Longitudes in any convention (e.g. `[-180, 180]`).

    Returns
    -------
    ndarray of float
        Longitudes wrapped into `[0, 360)`.
    """

    return np.mod(np.asarray(long, dtype = float), 360.0)

def _pad_longitude_seam(ds):

    """
    Pad the longitude axis by one cell on each side across the 0/360 seam.

    Without this, a query point just west of 0° (e.g. `359.9`) or just east
    of the last grid column has no bracketing cell on one side and
    interpolates to NaN, even though the data wraps continuously around the
    globe. Wrapping the data (`mode='wrap'`) reintroduces cell 0's field
    values next to the last column and vice versa, but `xr.Dataset.pad`
    wraps the `lon` *coordinate* values along with the data, which would
    leave the two new edge cells labelled with the wrong (already-used)
    coordinate — so those two coordinate values are overwritten with the
    true out-of-range values (`lon[0] - spacing`, `lon[-1] + spacing`),
    keeping the padded axis monotonically increasing and safe to
    interpolate against directly.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with a `lon` dimension/coordinate in `[0, 360)`.

    Returns
    -------
    xr.Dataset
        `ds` padded by one cell at each longitude edge.

    Raises
    ------
    ValueError
        If `lon` is not uniformly spaced.
    """

    lon = ds['lon'].values
    spacing = np.diff(lon)
    if not np.allclose(spacing, spacing[0]):
        raise ValueError("pyissm.data.grace._pad_longitude_seam: 'lon' must be uniformly spaced")

    padded = ds.pad(lon = (1, 1), mode = 'wrap')
    new_lon = np.concatenate(([lon[0] - spacing[0]], lon, [lon[-1] + spacing[0]]))

    return padded.assign_coords(lon = new_lon)

# ---------------------------------------------------------------------------
# Element centroids
# ---------------------------------------------------------------------------

def _centroids(lat, long, elements):

    """
    Compute element centroids as spherical (unit-vector) averages.

    Averaging lat/long directly breaks down near the poles and across the
    longitude seam (e.g. a triangle spanning 359° and 1° would naively
    average to 180°, the opposite side of the globe). Averaging the
    vertices' unit-sphere Cartesian coordinates instead, then converting the
    resultant vector back to lat/long, sidesteps both cases without any
    pole/seam special-casing.

    Parameters
    ----------
    lat, long : ndarray, shape (nv,)
        Vertex latitudes and longitudes, in degrees.
    elements : ndarray of int, shape (ne, 3)
        1-based vertex indices for each triangular element (pyISSM
        convention).

    Returns
    -------
    centroid_lat, centroid_long : ndarray, shape (ne,)
        Element centroid latitudes and longitudes, in degrees.
    """

    vertex_index = elements - 1
    lat_rad = np.radians(lat[vertex_index])
    long_rad = np.radians(long[vertex_index])

    x = (np.cos(lat_rad) * np.cos(long_rad)).mean(axis = 1)
    y = (np.cos(lat_rad) * np.sin(long_rad)).mean(axis = 1)
    z = np.sin(lat_rad).mean(axis = 1)

    centroid_lat = np.degrees(np.arctan2(z, np.hypot(x, y)))
    centroid_long = np.degrees(np.arctan2(y, x))

    return centroid_lat, centroid_long

# ---------------------------------------------------------------------------
# Per-month interpolation
# ---------------------------------------------------------------------------

def _interpolate_months(ds_padded, target_lat, target_long_0_360, time_indices):

    """
    Interpolate `lwe_thickness` onto target points, one month at a time.

    `xr_to_mesh` only accepts a single 2-D field per call, so each selected
    month is sliced out and interpolated separately.

    Parameters
    ----------
    ds_padded : xr.Dataset
        Dataset already padded across the longitude seam by
        `_pad_longitude_seam`.
    target_lat, target_long_0_360 : ndarray, shape (n,)
        Query latitudes and `[0, 360)`-convention longitudes (vertices or
        element centroids).
    time_indices : ndarray of int
        Time-axis indices to interpolate, as selected by
        `_select_epoch_indices`.

    Returns
    -------
    ndarray, shape (n, len(time_indices))
        Interpolated `lwe_thickness` (cm), one column per selected month.
    """

    water_load_cm = np.empty((len(target_lat), len(time_indices)))

    for month_index, time_index in enumerate(time_indices):
        ds_month = ds_padded.isel(time = int(time_index))
        water_load_cm[:, month_index] = xr_to_mesh(
            ds_month, LWE_THICKNESS_VAR,
            mesh_x = target_long_0_360, mesh_y = target_lat,
            x_var = 'lon', y_var = 'lat',
            interpolation_type = 'linear',
            issm_wrapper = False, crop_to_mesh = False)

    return water_load_cm

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def grace(md, tmin, tmax, filename, onvertex = True):

    """
    Interpolate GRACE liquid-water-equivalent thickness onto the model mesh.

    Reimplements the upstream MATLAB tutorial's `examples/Functions/grace.m`
    from scratch on top of `xarray` and `pyissm.data.interp.xr_to_mesh`,
    rather than porting its flattened-grid indexing, hand-rolled longitude-
    seam gap-filling, hardcoded fill-value thresholds, and `griddata` call.

    Parameters
    ----------
    md : pyissm.model.Model
        Model whose mesh (`md.mesh.lat`/`long`/`elements`) the data is
        interpolated onto.
    tmin, tmax : float
        Decimal-year bounds; the nearest available months are selected
        (inclusive). `tmin == tmax` selects a single epoch.
    filename : str
        Path to the GRACE NetCDF file. Required — the default data location
        is user-specific, so the caller (e.g. the tutorial notebook)
        supplies it; this function never hardcodes a path.
    onvertex : bool, optional
        True (default) returns per-vertex loads; False returns per-element
        loads, evaluated at each element's spherical centroid.

    Returns
    -------
    ndarray, shape (nv or ne, n_months)
        Water load in METRES of water equivalent, one column per selected
        month. NaNs (points outside GRACE's land-only coverage, or genuinely
        missing data) are replaced with 0.
    """

    ds = _open_grace_dataset(filename)
    try:
        decimal_year = _to_decimal_year(ds['time'].values)
        time_indices = _select_epoch_indices(decimal_year, tmin, tmax)
        ds_padded = _pad_longitude_seam(ds)

        if onvertex:
            target_lat, target_long = md.mesh.lat, md.mesh.long
        else:
            target_lat, target_long = _centroids(md.mesh.lat, md.mesh.long, md.mesh.elements)

        target_long_0_360 = _to_longitude_0_360(target_long)
        water_load_cm = _interpolate_months(ds_padded, target_lat, target_long_0_360, time_indices)
    finally:
        ds.close()

    water_load = water_load_cm / 100.0
    water_load[np.isnan(water_load)] = 0.0

    return water_load
