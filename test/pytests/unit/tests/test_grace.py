"""
Unit tests for pyissm.data.grace.

Uses a small synthetic in-memory GRACE-like NetCDF (never the real ~40 MB
product) so these tests run in CI without any external data file. A
separate, real-file test is gated on the PYISSM_GRACE_NC environment
variable and skips cleanly when it is unset.
"""

import os
from importlib import import_module
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

try:
    import xarray as xr
    XARRAY_AVAILABLE = True
except ImportError:
    XARRAY_AVAILABLE = False

pytestmark = pytest.mark.skipif(not XARRAY_AVAILABLE, reason="xarray not available")

if XARRAY_AVAILABLE:
    grace_module = import_module('pyissm.data.grace')

# Synthetic grid: 1-degree spacing, matching the real product's convention,
# but a small 8x6 patch (not global) to keep the fixture lightweight.
LON = np.arange(0.5, 8.5, 1.0)   # 0.5 ... 7.5 (8 points)
LAT = np.arange(-2.5, 3.5, 1.0)  # -2.5 ... 2.5 (6 points)
MONTHS = pd.date_range('2010-01-15', periods=12, freq='MS') + pd.Timedelta(days=14)
FILL_VALUE = 32767.0

# Grid indices of the single deliberately-missing cell, used by the NaN->0 test.
NAN_TIME_INDEX = 0
NAN_LAT_INDEX = 2
NAN_LON_INDEX = 3


def _make_synthetic_grace_dataset():

    """
    Build a small synthetic GRACE-like dataset with a known, deterministic
    `lwe_thickness` field: value(t, lat, lon) = 100*t_index + 10*lat_index +
    lon_index (cm), so exact-recovery checks can compare against a formula
    rather than a stored expected array. One cell is set to the raw
    `_FillValue` sentinel to exercise CF masking.
    """

    n_time, n_lat, n_lon = len(MONTHS), len(LAT), len(LON)
    time_idx, lat_idx, lon_idx = np.meshgrid(
        np.arange(n_time), np.arange(n_lat), np.arange(n_lon), indexing = 'ij')
    lwe = (100 * time_idx + 10 * lat_idx + lon_idx).astype(np.float32)
    lwe[NAN_TIME_INDEX, NAN_LAT_INDEX, NAN_LON_INDEX] = FILL_VALUE

    ds = xr.Dataset(
        data_vars = {
            'lwe_thickness': (('time', 'lat', 'lon'), lwe,
                               {'units': 'cm', '_FillValue': FILL_VALUE}),
        },
        coords = {
            'time': MONTHS,
            'lat': LAT,
            'lon': LON,
        },
    )
    return ds


@pytest.fixture
def grace_nc_path(tmp_path):
    """Write the synthetic dataset to a temporary NetCDF file and return its path."""
    ds = _make_synthetic_grace_dataset()
    path = tmp_path / 'synthetic_grace.nc'
    ds.to_netcdf(path)
    return str(path)


@pytest.fixture
def synthetic_md():
    """
    A minimal `md` stand-in (`SimpleNamespace`, not a real `pyissm.model.Model`)
    exposing just the `md.mesh.lat`/`long`/`elements` fields `grace()` reads.
    Avoids depending on a native `Model()`, which `grace()` itself never
    constructs or touches.
    """
    lat = np.array([0.0, 1.0, -1.0])
    long = np.array([2.5, 5.5, 0.5])
    elements = np.array([[1, 2, 3]])
    return SimpleNamespace(mesh = SimpleNamespace(lat = lat, long = long, elements = elements))


class TestOpenGraceDataset:

    def test_dims_and_fill_masking(self, grace_nc_path):
        ds = grace_module._open_grace_dataset(grace_nc_path)
        assert ds['lwe_thickness'].dims == ('time', 'lat', 'lon')
        assert np.isnan(ds['lwe_thickness'].values[NAN_TIME_INDEX, NAN_LAT_INDEX, NAN_LON_INDEX])
        ds.close()

    def test_missing_variable_raises(self, tmp_path):
        path = tmp_path / 'no_lwe.nc'
        xr.Dataset({'foo': (('x',), np.arange(3))}).to_netcdf(path)
        with pytest.raises(ValueError):
            grace_module._open_grace_dataset(str(path))


class TestDecimalYear:

    def test_matches_hand_computed_value(self):
        times = np.array(['2002-04-16T12:00:00'], dtype = 'datetime64[ns]')
        decimal_year = grace_module._to_decimal_year(times)[0]
        expected = 2002 + (pd.Timestamp('2002-04-16T12:00:00') - pd.Timestamp('2002-01-01')) / \
                   (pd.Timestamp('2003-01-01') - pd.Timestamp('2002-01-01'))
        assert decimal_year == pytest.approx(expected)

    def test_leap_year_end_of_year_stays_below_next_year(self):
        # 2016 is a leap year (366 days); Dec 31 should be just under 2017,
        # not overshoot it the way a fixed /365 divisor would.
        times = np.array(['2016-12-31T12:00:00'], dtype = 'datetime64[ns]')
        decimal_year = grace_module._to_decimal_year(times)[0]
        assert 2016.99 < decimal_year < 2017.0

    def test_monotonic_for_ascending_input(self):
        decimal_year = grace_module._to_decimal_year(MONTHS.values)
        assert np.all(np.diff(decimal_year) > 0)


class TestSelectEpochIndices:

    def setup_method(self):
        self.decimal_year = np.array([2010.0, 2010.2, 2010.4, 2010.6, 2010.8])

    def test_single_epoch_nearest_match(self):
        idx = grace_module._select_epoch_indices(self.decimal_year, 2010.19, 2010.19)
        assert list(idx) == [1]

    def test_inclusive_range(self):
        idx = grace_module._select_epoch_indices(self.decimal_year, 2010.15, 2010.65)
        assert list(idx) == [1, 2, 3]

    def test_tmax_before_tmin_raises(self):
        with pytest.raises(ValueError):
            grace_module._select_epoch_indices(self.decimal_year, 2010.8, 2010.0)


class TestPadLongitudeSeam:

    def test_padded_shape_and_monotonic(self, grace_nc_path):
        ds = grace_module._open_grace_dataset(grace_nc_path)
        padded = grace_module._pad_longitude_seam(ds)
        assert padded.sizes['lon'] == ds.sizes['lon'] + 2
        assert np.all(np.diff(padded['lon'].values) > 0)
        ds.close()

    def test_wraps_data_across_seam(self, grace_nc_path):
        ds = grace_module._open_grace_dataset(grace_nc_path)
        padded = grace_module._pad_longitude_seam(ds)
        first_month = ds['lwe_thickness'].isel(time = 1).values  # avoid the NaN month
        padded_first_month = padded['lwe_thickness'].isel(time = 1).values
        assert np.array_equal(padded_first_month[:, 0], first_month[:, -1])
        assert np.array_equal(padded_first_month[:, -1], first_month[:, 0])
        ds.close()

    def test_non_uniform_spacing_raises(self, grace_nc_path):
        ds = grace_module._open_grace_dataset(grace_nc_path)
        irregular_lon = np.concatenate([ds['lon'].values[:-1], [ds['lon'].values[-1] + 5.0]])
        bad = ds.assign_coords(lon = irregular_lon)
        with pytest.raises(ValueError):
            grace_module._pad_longitude_seam(bad)
        ds.close()


class TestCentroids:

    def test_planar_triangle_matches_simple_average(self):
        lat = np.array([0.0, 0.0, 10.0])
        long = np.array([0.0, 10.0, 5.0])
        elements = np.array([[1, 2, 3]])
        centroid_lat, centroid_long = grace_module._centroids(lat, long, elements)
        assert centroid_lat[0] == pytest.approx(10.0 / 3.0, abs = 1e-2)
        assert centroid_long[0] == pytest.approx(5.0, abs = 1e-6)

    def test_seam_straddling_triangle_centres_near_zero_not_180(self):
        lat = np.array([0.0, 0.0, 1.0])
        long = np.array([359.0, 1.0, 0.0])
        elements = np.array([[1, 2, 3]])
        _, centroid_long = grace_module._centroids(lat, long, elements)
        # A naive arithmetic mean of 359/1/0 would wrongly give 120; the
        # spherical average must stay near the seam itself.
        assert abs(centroid_long[0]) < 1.0 or abs(centroid_long[0] - 360) < 1.0


class TestGrace:

    def test_onvertex_shape_single_epoch(self, grace_nc_path, synthetic_md):
        decimal_year = grace_module._to_decimal_year(MONTHS.values)[1]
        water_load = grace_module.grace(synthetic_md, decimal_year, decimal_year, grace_nc_path, onvertex = True)
        assert water_load.shape == (3, 1)

    def test_onvertex_shape_multi_epoch(self, grace_nc_path, synthetic_md):
        dy = grace_module._to_decimal_year(MONTHS.values)
        water_load = grace_module.grace(synthetic_md, dy[0], dy[4], grace_nc_path, onvertex = True)
        assert water_load.shape == (3, 5)

    def test_onvertex_false_returns_one_row_per_element(self, grace_nc_path, synthetic_md):
        dy = grace_module._to_decimal_year(MONTHS.values)[1]
        water_load = grace_module.grace(synthetic_md, dy, dy, grace_nc_path, onvertex = False)
        assert water_load.shape == (1, 1)

    def test_cm_to_m_conversion_and_known_point_recovery(self, grace_nc_path):
        # Query exactly at a grid vertex: linear interpolation must recover
        # the source value exactly, so the /100 conversion can be checked
        # precisely rather than approximately.
        dy = grace_module._to_decimal_year(MONTHS.values)[1]
        lat_idx, lon_idx = 3, 4
        md = SimpleNamespace(mesh = SimpleNamespace(
            lat = np.array([LAT[lat_idx]]), long = np.array([LON[lon_idx]]), elements = np.array([[1, 1, 1]])))
        water_load = grace_module.grace(md, dy, dy, grace_nc_path, onvertex = True)
        expected_cm = 100 * 1 + 10 * lat_idx + lon_idx  # matches the fixture's formula, time index 1
        assert water_load[0, 0] == pytest.approx(expected_cm / 100.0, abs = 1e-6)

    def test_nan_filled_to_zero(self, grace_nc_path):
        # Query exactly at the deliberately-missing cell's location, at the
        # same month it was blanked out.
        dy = grace_module._to_decimal_year(MONTHS.values)[NAN_TIME_INDEX]
        md = SimpleNamespace(mesh = SimpleNamespace(
            lat = np.array([LAT[NAN_LAT_INDEX]]), long = np.array([LON[NAN_LON_INDEX]]), elements = np.array([[1, 1, 1]])))
        water_load = grace_module.grace(md, dy, dy, grace_nc_path, onvertex = True)
        assert water_load[0, 0] == 0.0


class TestGraceRealFile:
    """
    Gated on the real ~40 MB GRACE product, via PYISSM_GRACE_NC. Never
    fabricated or bundled in the repo; skips cleanly when the env var is
    unset. Not a golden-fixture comparison (see grace_golden/, Phase 6) —
    just a basic real-data sanity check.
    """

    @pytest.mark.skipif('PYISSM_GRACE_NC' not in os.environ, reason = "PYISSM_GRACE_NC not set")
    def test_runs_against_real_file(self):
        filename = os.environ['PYISSM_GRACE_NC']
        md = SimpleNamespace(mesh = SimpleNamespace(
            lat = np.array([23.0, -25.0]), long = np.array([13.0, 133.0]), elements = np.array([[1, 2, 1]])))
        water_load = grace_module.grace(md, 2010.0, 2010.0, filename, onvertex = True)
        assert water_load.shape == (2, 1)
        assert np.all(np.isfinite(water_load))
