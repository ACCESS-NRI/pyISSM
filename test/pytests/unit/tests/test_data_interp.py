"""
Unit tests for the pyissm.data.interp module.

Tests cover:
- xr_to_mesh: gridded (rectilinear) data -> mesh nodes
- points_to_mesh: scattered data -> mesh nodes

Both functions are pure Python (scipy backend), so these tests do not
require the compiled ISSM wrappers. The ISSM wrapper code path
(issm_wrapper=True) is not exercised here.
"""

import numpy as np
import pytest
import xarray as xr

from pyissm.data.interp import xr_to_mesh, points_to_mesh


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def make_grid_dataset(nx = 11,
                      ny = 11,
                      x_range = (0.0, 10.0),
                      y_range = (0.0, 10.0),
                      var_name = 'var',
                      descending_x = False,
                      descending_y = False,
                      coords_2d = False):
    """
    Build a rectilinear xarray Dataset holding the plane var = x + y.

    A linear field means bilinear interpolation is exact, so interpolated
    values can be compared against the analytic result.
    """
    x = np.linspace(*x_range, nx)
    y = np.linspace(*y_range, ny)

    if descending_x:
        x = x[::-1]
    if descending_y:
        y = y[::-1]

    grid_x, grid_y = np.meshgrid(x, y)
    values = grid_x + grid_y

    if coords_2d:
        return xr.Dataset(
            {var_name: (('y', 'x'), values)},
            coords = {'x': (('y', 'x'), grid_x),
                      'y': (('y', 'x'), grid_y)},
        )

    return xr.Dataset({var_name: (('y', 'x'), values)},
                      coords = {'x': x, 'y': y})


# Mesh nodes wholly inside the [0, 10] x [0, 10] source domain
INSIDE_X = np.array([2.0, 5.0, 8.0])
INSIDE_Y = np.array([2.0, 5.0, 8.0])

# As above, plus one node well outside the source domain
OUTSIDE_X = np.array([2.0, 5.0, 8.0, 20.0])
OUTSIDE_Y = np.array([2.0, 5.0, 8.0, 20.0])


# ==================================================================
# xr_to_mesh -- characterisation of existing behaviour
# ==================================================================

class TestXrToMeshBasics:
    """Behaviour that must not change."""

    def test_linear_field_is_interpolated_exactly(self):
        data = make_grid_dataset()
        result = xr_to_mesh(data, 'var', INSIDE_X, INSIDE_Y,
                            interpolation_type = 'linear',
                            issm_wrapper = False)

        np.testing.assert_allclose(result, INSIDE_X + INSIDE_Y)

    def test_accepts_a_netcdf_path(self, tmp_path):
        data = make_grid_dataset()
        path = tmp_path / 'grid.nc'
        data.to_netcdf(path)

        from_path = xr_to_mesh(str(path), 'var', INSIDE_X, INSIDE_Y,
                               interpolation_type = 'linear',
                               issm_wrapper = False)
        from_dataset = xr_to_mesh(data, 'var', INSIDE_X, INSIDE_Y,
                                  interpolation_type = 'linear',
                                  issm_wrapper = False)

        np.testing.assert_allclose(from_path, from_dataset)

    @pytest.mark.parametrize('descending_x, descending_y',
                             [(True, False), (False, True), (True, True)])
    def test_descending_axes_are_flipped(self, descending_x, descending_y):
        data = make_grid_dataset(descending_x = descending_x,
                                 descending_y = descending_y)
        result = xr_to_mesh(data, 'var', INSIDE_X, INSIDE_Y,
                            interpolation_type = 'linear',
                            issm_wrapper = False)

        np.testing.assert_allclose(result, INSIDE_X + INSIDE_Y)

    def test_accepts_2d_rectilinear_coordinates(self):
        data = make_grid_dataset(coords_2d = True)
        result = xr_to_mesh(data, 'var', INSIDE_X, INSIDE_Y,
                            interpolation_type = 'linear',
                            issm_wrapper = False)

        np.testing.assert_allclose(result, INSIDE_X + INSIDE_Y)

    def test_cropping_does_not_change_the_result(self):
        data = make_grid_dataset(nx = 51, ny = 51)

        cropped = xr_to_mesh(data, 'var', INSIDE_X, INSIDE_Y,
                             interpolation_type = 'linear',
                             issm_wrapper = False,
                             crop_to_mesh = True, crop_buffer = 1.0)
        uncropped = xr_to_mesh(data, 'var', INSIDE_X, INSIDE_Y,
                               interpolation_type = 'linear',
                               issm_wrapper = False,
                               crop_to_mesh = False)

        np.testing.assert_allclose(cropped, uncropped)

    def test_nearest_interpolation_is_supported(self):
        data = make_grid_dataset()
        result = xr_to_mesh(data, 'var', INSIDE_X, INSIDE_Y,
                            interpolation_type = 'nearest',
                            issm_wrapper = False)

        np.testing.assert_allclose(result, INSIDE_X + INSIDE_Y)


class TestXrToMeshErrors:
    """Error paths that must keep raising."""

    def test_rejects_non_dataset_input(self):
        with pytest.raises(TypeError):
            xr_to_mesh(42, 'var', INSIDE_X, INSIDE_Y, issm_wrapper = False)

    def test_rejects_a_variable_that_is_not_2d(self):
        data = make_grid_dataset()
        data['var3d'] = (('t', 'y', 'x'),
                         np.repeat(data['var'].values[None, :, :], 3, axis = 0))

        with pytest.raises(ValueError, match = 'must be 2D'):
            xr_to_mesh(data, 'var3d', INSIDE_X, INSIDE_Y,
                       issm_wrapper = False, crop_to_mesh = False)

    def test_rejects_an_unsupported_interpolation_type(self):
        data = make_grid_dataset()
        with pytest.raises(ValueError, match = 'not supported by scipy'):
            xr_to_mesh(data, 'var', INSIDE_X, INSIDE_Y,
                       interpolation_type = 'bilinear',
                       issm_wrapper = False)


class TestXrToMeshDefaultValue:
    """Constant fill outside the source domain (fill_nan disabled)."""

    def test_outside_nodes_receive_default_value(self):
        data = make_grid_dataset()
        result = xr_to_mesh(data, 'var', OUTSIDE_X, OUTSIDE_Y,
                            default_value = -999.0,
                            interpolation_type = 'linear',
                            issm_wrapper = False)

        np.testing.assert_allclose(result[:3], INSIDE_X + INSIDE_Y)
        assert result[3] == -999.0

    def test_outside_nodes_are_nan_by_default(self):
        data = make_grid_dataset()
        result = xr_to_mesh(data, 'var', OUTSIDE_X, OUTSIDE_Y,
                            interpolation_type = 'linear',
                            issm_wrapper = False)

        assert np.isnan(result[3])

    @pytest.mark.parametrize('default_value', [np.nan, 0.0, -999.0])
    def test_fill_nan_disabled_is_a_plain_constant_fill(self, default_value):
        """The fill_nan=False path must stay byte-identical."""
        data = make_grid_dataset()
        result = xr_to_mesh(data, 'var', OUTSIDE_X, OUTSIDE_Y,
                            default_value = default_value,
                            interpolation_type = 'linear',
                            issm_wrapper = False,
                            fill_nan = False)

        np.testing.assert_allclose(result[:3], INSIDE_X + INSIDE_Y)
        if np.isnan(default_value):
            assert np.isnan(result[3])
        else:
            assert result[3] == default_value


# ==================================================================
# xr_to_mesh -- nearest-neighbour fallback fill
# ==================================================================

class TestXrToMeshFillNan:
    """
    The second (fallback) interpolation pass.

    Nodes that the primary interpolation could not resolve -- because they lie
    outside the source domain, or because the source is NaN there -- are given
    the value of the nearest resolved node instead of a constant.
    """

    def test_nan_default_value_is_filled(self):
        data = make_grid_dataset()
        result = xr_to_mesh(data, 'var', OUTSIDE_X, OUTSIDE_Y,
                            default_value = np.nan,
                            interpolation_type = 'linear',
                            issm_wrapper = False,
                            fill_nan = True)

        np.testing.assert_allclose(result[:3], INSIDE_X + INSIDE_Y)
        # Nearest resolved node to (20, 20) is (8, 8) -> 16
        assert result[3] == pytest.approx(16.0)

    def test_finite_default_value_is_still_filled(self):
        """
        Regression: the fallback pass used to key off np.isfinite, so a finite
        default_value made every node look resolved and the pass silently did
        nothing. default_value is a last resort, not a blocker.
        """
        data = make_grid_dataset()
        result = xr_to_mesh(data, 'var', OUTSIDE_X, OUTSIDE_Y,
                            default_value = 0.0,
                            interpolation_type = 'linear',
                            issm_wrapper = False,
                            fill_nan = True)

        np.testing.assert_allclose(result[:3], INSIDE_X + INSIDE_Y)
        assert result[3] != 0.0
        assert result[3] == pytest.approx(16.0)

    def test_interior_source_nans_are_filled(self):
        """
        The shelf-edge case: the source product is undefined (NaN) over part of
        its own domain, so nodes there are unresolved even though they are
        inside the grid.
        """
        data = make_grid_dataset(nx = 21, ny = 21)
        values = data['var'].values.copy()
        values[15:, 15:] = np.nan          # blank the top-right corner
        data['var'] = (('y', 'x'), values)

        mesh_x = np.array([2.0, 5.0, 9.5])
        mesh_y = np.array([2.0, 5.0, 9.5])

        result = xr_to_mesh(data, 'var', mesh_x, mesh_y,
                            default_value = 0.0,
                            interpolation_type = 'linear',
                            issm_wrapper = False,
                            fill_nan = True)

        assert np.all(np.isfinite(result))
        # The blanked node takes the nearest resolved value, not 0
        assert result[2] == pytest.approx(10.0)

    def test_mesh_entirely_outside_the_source_domain_is_filled(self):
        data = make_grid_dataset()
        mesh_x = np.array([100.0, 200.0])
        mesh_y = np.array([100.0, 200.0])

        with pytest.raises(ValueError, match = 'no valid data points'):
            xr_to_mesh(data, 'var', mesh_x, mesh_y,
                       default_value = 0.0,
                       interpolation_type = 'linear',
                       issm_wrapper = False,
                       fill_nan = True,
                       crop_to_mesh = False)

    def test_fully_resolved_mesh_is_untouched(self):
        """With nothing to fill, enabling the pass must change nothing."""
        data = make_grid_dataset()

        without = xr_to_mesh(data, 'var', INSIDE_X, INSIDE_Y,
                             interpolation_type = 'linear',
                             issm_wrapper = False, fill_nan = False)
        with_fill = xr_to_mesh(data, 'var', INSIDE_X, INSIDE_Y,
                               interpolation_type = 'linear',
                               issm_wrapper = False, fill_nan = True)

        np.testing.assert_array_equal(without, with_fill)

    def test_rejects_an_unsupported_fill_interpolation_type(self):
        data = make_grid_dataset()
        with pytest.raises(ValueError, match = 'not supported by scipy'):
            xr_to_mesh(data, 'var', OUTSIDE_X, OUTSIDE_Y,
                       issm_wrapper = False,
                       fill_nan = True,
                       fill_nan_interpolation_type = 'bilinear')


# ==================================================================
# points_to_mesh -- characterisation of existing behaviour
# ==================================================================

def scattered_source(n = 21):
    """Scattered samples of var = x + y over [0, 10] x [0, 10]."""
    x = np.linspace(0.0, 10.0, n)
    y = np.linspace(0.0, 10.0, n)
    grid_x, grid_y = np.meshgrid(x, y)
    return grid_x.ravel(), grid_y.ravel(), (grid_x + grid_y).ravel()


class TestPointsToMeshBasics:

    def test_linear_field_is_interpolated_exactly(self):
        data_x, data_y, data_values = scattered_source()
        result = points_to_mesh(data_x, data_y, data_values,
                                INSIDE_X, INSIDE_Y)

        np.testing.assert_allclose(result, INSIDE_X + INSIDE_Y, atol = 1e-10)

    def test_2d_inputs_are_flattened(self):
        x = np.linspace(0.0, 10.0, 21)
        y = np.linspace(0.0, 10.0, 21)
        grid_x, grid_y = np.meshgrid(x, y)

        result = points_to_mesh(grid_x, grid_y, grid_x + grid_y,
                                INSIDE_X, INSIDE_Y)

        np.testing.assert_allclose(result, INSIDE_X + INSIDE_Y, atol = 1e-10)

    def test_non_finite_source_points_are_dropped(self):
        data_x, data_y, data_values = scattered_source()
        data_values = data_values.copy()
        data_values[0] = np.nan
        data_values[1] = np.inf

        result = points_to_mesh(data_x, data_y, data_values,
                                INSIDE_X, INSIDE_Y)

        assert np.all(np.isfinite(result))

    def test_outside_the_convex_hull_receives_default_value(self):
        data_x, data_y, data_values = scattered_source()
        result = points_to_mesh(data_x, data_y, data_values,
                                OUTSIDE_X, OUTSIDE_Y,
                                default_value = -999.0)

        assert result[3] == -999.0

    def test_rejects_an_unsupported_interpolation_type(self):
        data_x, data_y, data_values = scattered_source()
        with pytest.raises(ValueError, match = 'not supported by scipy'):
            points_to_mesh(data_x, data_y, data_values,
                           INSIDE_X, INSIDE_Y,
                           interpolation_type = 'quintic')

    def test_rejects_mismatched_input_shapes(self):
        data_x, data_y, data_values = scattered_source()
        with pytest.raises(ValueError, match = 'same shape'):
            points_to_mesh(data_x, data_y[:-1], data_values,
                           INSIDE_X, INSIDE_Y)

    def test_rejects_an_all_nan_source(self):
        data_x, data_y, _ = scattered_source()
        data_values = np.full_like(data_x, np.nan)

        with pytest.raises(ValueError, match = 'no valid data points'):
            points_to_mesh(data_x, data_y, data_values,
                           INSIDE_X, INSIDE_Y)


class TestPointsToMeshFillNan:
    """The fallback pass, absent from points_to_mesh until now."""

    def test_outside_the_convex_hull_is_filled(self):
        data_x, data_y, data_values = scattered_source()
        result = points_to_mesh(data_x, data_y, data_values,
                                OUTSIDE_X, OUTSIDE_Y,
                                default_value = 0.0,
                                fill_nan = True)

        np.testing.assert_allclose(result[:3], INSIDE_X + INSIDE_Y, atol = 1e-10)
        # Nearest source sample to (20, 20) is the corner (10, 10) -> 20
        assert result[3] == pytest.approx(20.0)

    def test_fill_nan_disabled_is_unchanged(self):
        data_x, data_y, data_values = scattered_source()

        without = points_to_mesh(data_x, data_y, data_values,
                                 OUTSIDE_X, OUTSIDE_Y,
                                 default_value = -999.0)
        with_flag_off = points_to_mesh(data_x, data_y, data_values,
                                       OUTSIDE_X, OUTSIDE_Y,
                                       default_value = -999.0,
                                       fill_nan = False)

        np.testing.assert_array_equal(without, with_flag_off)

    def test_fully_resolved_mesh_is_untouched(self):
        data_x, data_y, data_values = scattered_source()

        without = points_to_mesh(data_x, data_y, data_values,
                                 INSIDE_X, INSIDE_Y)
        with_fill = points_to_mesh(data_x, data_y, data_values,
                                   INSIDE_X, INSIDE_Y,
                                   fill_nan = True)

        np.testing.assert_array_equal(without, with_fill)

    def test_nearest_primary_interpolation_needs_no_fill(self):
        data_x, data_y, data_values = scattered_source()
        result = points_to_mesh(data_x, data_y, data_values,
                                OUTSIDE_X, OUTSIDE_Y,
                                interpolation_type = 'nearest',
                                fill_nan = True)

        assert np.all(np.isfinite(result))

    def test_rejects_an_unsupported_fill_interpolation_type(self):
        data_x, data_y, data_values = scattered_source()
        with pytest.raises(ValueError, match = 'not supported by scipy'):
            points_to_mesh(data_x, data_y, data_values,
                           INSIDE_X, INSIDE_Y,
                           fill_nan = True,
                           fill_nan_interpolation_type = 'bilinear')
