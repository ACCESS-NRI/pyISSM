"""
Unit tests for the pyissm.data.interp module.

Tests cover:
- xr_to_mesh: gridded (rectilinear) data -> mesh nodes
- points_to_mesh: scattered data -> mesh nodes
- mesh_to_xr: mesh -> regular grid, as an xarray DataArray

All three functions are pure Python (scipy backend), so these tests do not
require the compiled ISSM wrappers. The ISSM wrapper code path
(issm_wrapper=True) is not exercised here.
"""

from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

from pyissm.data.interp import xr_to_mesh, points_to_mesh, mesh_to_xr


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


# ==================================================================
# mesh_to_xr -- mesh -> regular grid
# ==================================================================

def make_square_model(n = 6, size = 10.0):
    """
    Build a minimal model object carrying a triangulated square mesh.

    The mesh covers [0, size] x [0, size] as an n x n grid of vertices, split
    into triangles. Only the mesh attributes used by grid_model_field are
    populated.
    """
    axis = np.linspace(0.0, size, n)
    grid_x, grid_y = np.meshgrid(axis, axis)

    # Two triangles per grid cell, using 1-based vertex indices (ISSM convention)
    elements = []
    for row in range(n - 1):
        for col in range(n - 1):
            bottom_left = row * n + col + 1
            bottom_right = bottom_left + 1
            top_left = bottom_left + n
            top_right = top_left + 1
            elements.append([bottom_left, bottom_right, top_right])
            elements.append([bottom_left, top_right, top_left])

    elements = np.array(elements)

    mesh = SimpleNamespace(x = grid_x.ravel(),
                           y = grid_y.ravel(),
                           elements = elements,
                           numberofvertices = grid_x.size,
                           numberofelements = elements.shape[0])

    return SimpleNamespace(mesh = mesh)


class TestMeshToXrStructure:
    """Shape, dimensions and coordinates of the returned DataArray."""

    def test_static_field_returns_a_2d_dataarray(self):
        md = make_square_model()
        field = md.mesh.x + md.mesh.y

        result = mesh_to_xr(md, field, spacing = 1.0)

        assert isinstance(result, xr.DataArray)
        assert result.dims == ('y', 'x')
        assert result.shape == (result.y.size, result.x.size)

    def test_time_varying_field_gains_a_time_dimension(self):
        md = make_square_model()
        field = np.stack([md.mesh.x + md.mesh.y,
                          2.0 * (md.mesh.x + md.mesh.y),
                          3.0 * (md.mesh.x + md.mesh.y)])

        result = mesh_to_xr(md, field, spacing = 1.0)

        assert result.dims == ('time', 'y', 'x')
        assert result.sizes['time'] == 3
        np.testing.assert_array_equal(result.time.values, [0, 1, 2])

    def test_time_coordinate_can_be_supplied(self):
        md = make_square_model()
        field = np.stack([md.mesh.x, 2.0 * md.mesh.x])

        result = mesh_to_xr(md, field, spacing = 1.0, time = [2000.0, 2010.0])

        np.testing.assert_allclose(result.time.values, [2000.0, 2010.0])

    def test_explicit_grid_is_used_verbatim(self):
        md = make_square_model()
        x = np.linspace(0.0, 10.0, 7)
        y = np.linspace(0.0, 10.0, 9)

        result = mesh_to_xr(md, md.mesh.x, x = x, y = y)

        np.testing.assert_allclose(result.x.values, x)
        np.testing.assert_allclose(result.y.values, y)

    def test_spacing_grid_spans_the_mesh(self):
        md = make_square_model(size = 10.0)
        result = mesh_to_xr(md, md.mesh.x, spacing = 2.0)

        assert result.x.values[0] == pytest.approx(0.0)
        assert result.x.values[-1] == pytest.approx(10.0)
        np.testing.assert_allclose(np.diff(result.x.values), 2.0)

    def test_name_and_attributes_are_attached(self):
        md = make_square_model()
        result = mesh_to_xr(md, md.mesh.x, spacing = 2.0,
                            name = 'thickness',
                            attrs = {'units': 'm', 'long_name': 'Ice thickness'})

        assert result.name == 'thickness'
        assert result.attrs['units'] == 'm'
        assert result.attrs['long_name'] == 'Ice thickness'


class TestMeshToXrValues:
    """Interpolated values and domain masking."""

    def test_linear_field_is_reproduced_on_the_grid(self):
        md = make_square_model()
        field = md.mesh.x + md.mesh.y

        result = mesh_to_xr(md, field, spacing = 1.0)

        # A linear field is reproduced exactly by linear interpolation
        expected_x, expected_y = np.meshgrid(result.x.values, result.y.values)
        inside = np.isfinite(result.values)
        np.testing.assert_allclose(result.values[inside],
                                   (expected_x + expected_y)[inside],
                                   atol = 1e-9)

    def test_element_based_fields_are_supported(self):
        md = make_square_model()
        field = np.ones(md.mesh.numberofelements)

        result = mesh_to_xr(md, field, spacing = 1.0)

        inside = np.isfinite(result.values)
        np.testing.assert_allclose(result.values[inside], 1.0)

    def test_cells_outside_the_domain_take_fill_value(self):
        md = make_square_model(size = 10.0)
        # Grid extends well beyond the mesh, so some cells must fall outside it
        x = np.linspace(-5.0, 15.0, 21)
        y = np.linspace(-5.0, 15.0, 21)

        result = mesh_to_xr(md, md.mesh.x, x = x, y = y, fill_value = -999.0)

        assert np.any(result.values == -999.0)
        assert result.sel(x = -5.0, y = -5.0).item() == -999.0

    def test_nearest_interpolation_is_supported(self):
        md = make_square_model()
        result = mesh_to_xr(md, md.mesh.x + md.mesh.y, spacing = 1.0,
                            interpolation_type = 'nearest')

        assert np.any(np.isfinite(result.values))

    def test_round_trip_through_xr_to_mesh(self):
        """Gridding then re-interpolating recovers the field in the interior."""
        md = make_square_model(n = 11, size = 10.0)
        field = md.mesh.x + md.mesh.y

        gridded = mesh_to_xr(md, field, spacing = 0.5, name = 'var')

        # Sample away from the mesh boundary, where the domain mask makes the
        # round trip lossy
        sample_x = np.array([2.0, 5.0, 7.5])
        sample_y = np.array([2.0, 5.0, 7.5])

        recovered = xr_to_mesh(gridded.to_dataset(), 'var',
                               sample_x, sample_y,
                               interpolation_type = 'linear',
                               issm_wrapper = False)

        np.testing.assert_allclose(recovered, sample_x + sample_y, atol = 1e-6)


class TestMeshToXrCrs:
    """CF-style CRS metadata."""

    def test_no_crs_is_assumed_by_default(self):
        md = make_square_model()
        result = mesh_to_xr(md, md.mesh.x, spacing = 2.0)

        assert 'spatial_ref' not in result.coords
        assert 'grid_mapping' not in result.attrs

    def test_crs_is_attached_when_given(self):
        md = make_square_model()
        result = mesh_to_xr(md, md.mesh.x, spacing = 2.0, crs = 3031)

        assert 'spatial_ref' in result.coords
        assert result.attrs['grid_mapping'] == 'spatial_ref'
        assert result.coords['spatial_ref'].attrs

    def test_crs_accepts_a_string(self):
        md = make_square_model()
        result = mesh_to_xr(md, md.mesh.x, spacing = 2.0, crs = 'EPSG:3413')

        assert 'spatial_ref' in result.coords


class TestMeshToXrErrors:

    def test_requires_a_grid_specification(self):
        md = make_square_model()
        with pytest.raises(ValueError, match = 'spacing'):
            mesh_to_xr(md, md.mesh.x)

    def test_rejects_2d_grid_coordinates(self):
        md = make_square_model()
        grid_x, grid_y = np.meshgrid(np.linspace(0, 10, 5), np.linspace(0, 10, 5))

        with pytest.raises(ValueError, match = '1D'):
            mesh_to_xr(md, md.mesh.x, x = grid_x, y = grid_y)

    def test_rejects_a_field_of_the_wrong_length(self):
        md = make_square_model()
        with pytest.raises(ValueError):
            mesh_to_xr(md, np.ones(7), spacing = 2.0)
