"""
Tests for pyissm.model.mesh.grid_model_field and its reusable interpolator (issue #23: gridding is slow).

The interpolator computes the triangulation and weights once per mesh/grid; results must match
scipy.interpolate.griddata for every method.
"""

import numpy as np
import pytest
from types import SimpleNamespace
from scipy.interpolate import griddata
from scipy.spatial import Delaunay

from pyissm.model import mesh


@pytest.fixture(scope = 'module')
def setup():
    """Irregular triangulated mesh on [0, 1e5]^2 (corners included) and a grid extending past it."""
    rng = np.random.default_rng(42)
    pts = np.vstack((rng.uniform(0, 1e5, (400, 2)), [[0, 0], [1e5, 0], [0, 1e5], [1e5, 1e5]]))
    elements = Delaunay(pts).simplices + 1
    md = SimpleNamespace(mesh = SimpleNamespace(x = pts[:, 0], y = pts[:, 1], elements = elements,
                                                numberofvertices = len(pts), numberofelements = len(elements)))
    grid_x, grid_y = np.meshgrid(np.linspace(-1e4, 1.1e5, 37), np.linspace(-1e4, 1.1e5, 41))
    return md, pts, grid_x, grid_y, rng


def _vertex_field(pts, scale = 1.):
    return scale * (np.sin(pts[:, 0] / 2e4) + np.cos(pts[:, 1] / 3e4))


class TestGridInterpolator:

    @pytest.mark.parametrize('method', ['linear', 'nearest', 'cubic'])
    def test_matches_griddata(self, setup, method):
        md, pts, grid_x, grid_y, _ = setup
        values = _vertex_field(pts)
        expected = griddata(pts, values, (grid_x, grid_y), method = method)
        result = mesh._grid_interpolator(pts, grid_x, grid_y, method)(values)
        assert result.shape == grid_x.shape
        np.testing.assert_array_equal(np.isnan(result), np.isnan(expected))
        np.testing.assert_allclose(result, expected, rtol = 0, atol = 1e-12, equal_nan = True)

    def test_reused_for_many_fields(self, setup):
        md, pts, grid_x, grid_y, rng = setup
        interpolate = mesh._grid_interpolator(pts, grid_x, grid_y, 'linear')
        for _ in range(3):
            values = rng.random(len(pts))
            np.testing.assert_allclose(interpolate(values), griddata(pts, values, (grid_x, grid_y)),
                                       rtol = 0, atol = 1e-12, equal_nan = True)

    def test_nan_values_propagate(self, setup):
        md, pts, grid_x, grid_y, _ = setup
        values = _vertex_field(pts)
        values[0] = np.nan
        expected = griddata(pts, values, (grid_x, grid_y))
        result = mesh._grid_interpolator(pts, grid_x, grid_y, 'linear')(values)
        np.testing.assert_array_equal(np.isnan(result), np.isnan(expected))

    def test_unknown_method(self, setup):
        md, pts, grid_x, grid_y, _ = setup
        with pytest.raises(ValueError, match = 'unknown interpolation method'):
            mesh._grid_interpolator(pts, grid_x, grid_y, 'spline')


class TestGridModelField:

    def _griddata_reference(self, md, field, grid_x, grid_y, method):
        """grid_model_field as implemented before #23: griddata per step, masked to the mesh."""
        el = md.mesh.elements - 1
        if field.shape[-1] == md.mesh.numberofvertices:
            points = np.column_stack((md.mesh.x, md.mesh.y))
        else:
            points = np.column_stack((md.mesh.x[el].mean(axis = 1), md.mesh.y[el].mean(axis = 1)))
        out = np.stack([griddata(points, f, (grid_x, grid_y), method = method) for f in np.atleast_2d(field)])
        mask = mesh.make_gridded_domain_mask(md.mesh.x, md.mesh.y, el, grid_x, grid_y)
        out[:, ~mask] = np.nan
        return out.squeeze()

    @pytest.mark.parametrize('method', ['linear', 'nearest', 'cubic'])
    def test_transient_vertex_field(self, setup, method):
        md, pts, grid_x, grid_y, _ = setup
        field = np.stack([_vertex_field(pts, scale) for scale in (1., 2., 3.)])
        result = mesh.grid_model_field(md, field, grid_x, grid_y, method = method)
        assert result.shape == (3,) + grid_x.shape
        np.testing.assert_allclose(result, self._griddata_reference(md, field, grid_x, grid_y, method),
                                   rtol = 0, atol = 1e-12, equal_nan = True)

    @pytest.mark.parametrize('method', ['linear', 'nearest'])
    def test_element_field(self, setup, method):
        md, pts, grid_x, grid_y, rng = setup
        field = rng.random((2, md.mesh.numberofelements))
        result = mesh.grid_model_field(md, field, grid_x, grid_y, method = method)
        np.testing.assert_allclose(result, self._griddata_reference(md, field, grid_x, grid_y, method),
                                   rtol = 0, atol = 1e-12, equal_nan = True)

    def test_static_field(self, setup):
        md, pts, grid_x, grid_y, _ = setup
        result = mesh.grid_model_field(md, _vertex_field(pts), grid_x, grid_y)
        assert result.shape == grid_x.shape
        assert np.isnan(result[0, 0]) and np.isfinite(result[20, 18])  # outside / inside the mesh

    def test_domain_mask_and_fill_value(self, setup):
        md, pts, grid_x, grid_y, _ = setup
        domain_mask = grid_x < 5e4
        result = mesh.grid_model_field(md, _vertex_field(pts), grid_x, grid_y, domain_mask = domain_mask, fill_value = -9999.)
        assert np.all(result[~domain_mask] == -9999.)
        assert np.all(result[domain_mask & (grid_x > 1e4) & (grid_y > 1e4) & (grid_y < 9e4)] != -9999.)
