"""
Ocean/land classification for ISSM model mesh points.

This module classifies geographic points as ocean or land using Natural
Earth land polygons and a vectorised point-in-polygon test, replacing the
upstream ISSM tutorial's shell-out to the GMT CLI (`gmt select`).
"""

import numpy as np
import geopandas as gpd
from shapely.ops import unary_union

# Natural Earth land-polygon resolutions this module knows how to fetch.
VALID_RESOLUTIONS = ('10m', '50m', '110m')

def _load_natural_earth_land(resolution = '50m'):

    """
    Load Natural Earth land polygons as a list of shapely geometries.

    Parameters
    ----------
    resolution : {'10m', '50m', '110m'}, optional
        Natural Earth land-polygon resolution. Default '50m'.

    Returns
    -------
    list of shapely.Geometry
        Land polygon/multipolygon geometries in EPSG:4326.

    Raises
    ------
    ValueError
        If `resolution` is not one of '10m', '50m', '110m'.
    RuntimeError
        If the `cartopy` package is not installed.
    """

    if resolution not in VALID_RESOLUTIONS:
        raise ValueError("pyissm.data.ocean_mask._load_natural_earth_land: "
                          f"resolution must be one of {VALID_RESOLUTIONS}, got {resolution!r}")

    try:
        import cartopy.feature as cfeature
    except ImportError:
        raise RuntimeError("pyissm.data.ocean_mask._load_natural_earth_land: the 'cartopy' package is "
                            "required (conda install -c conda-forge cartopy)")

    return list(cfeature.NaturalEarthFeature('physical', 'land', resolution).geometries())

def gmtmask(lat, long, resolution = '50m'):

    """
    Classify points as ocean or land using Natural Earth land polygons.

    Parameters
    ----------
    lat : array-like
        Point latitudes, in degrees [-90, 90].
    long : array-like
        Point longitudes, in degrees [-180, 180].
    resolution : {'10m', '50m', '110m'}, optional
        Natural Earth land-polygon resolution. Default '50m'.

    Returns
    -------
    ndarray of int, shape (n,)
        1 where the point is over ocean, 0 where over land (matches the
        upstream MATLAB `gmtmask` convention: `land = (ocean_mask == 0)`).

    Raises
    ------
    ValueError
        If `resolution` is not one of '10m', '50m', '110m'.
    RuntimeError
        If the `cartopy` package is not installed.
    """

    lat = np.asarray(lat, dtype = float)
    long = np.asarray(long, dtype = float)

    land = unary_union(_load_natural_earth_land(resolution))
    points = gpd.GeoSeries(gpd.points_from_xy(long, lat), crs = 'EPSG:4326')
    is_land = points.within(land).to_numpy()

    return (~is_land).astype(int)
