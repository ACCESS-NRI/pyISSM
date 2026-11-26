"""
Tests for pyissm.plot.plot_model_field colour normalisation (issue #183: vmin/vmax with log = True).
"""

import numpy as np
import pytest

matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

from pyissm.tools import wrappers


def _wrappers_installed():
    try:
        return bool(wrappers.check_wrappers_installed())
    except Exception:
        return False


# square_mesh calls the compiled ISSM wrappers (NodeConnectivity, ElementConnectivity)
pytestmark = pytest.mark.skipif(not _wrappers_installed(), reason = 'ISSM Python wrappers not installed')


@pytest.fixture
def md_and_field():
    from pyissm.model import Model, mesh
    md = mesh.square_mesh(Model(), 1.e5, 1.e5, 6, 6)
    return md, np.logspace(0, 3, md.mesh.numberofvertices)


@pytest.fixture(autouse = True)
def _close_figures():
    yield
    plt.close('all')


def _plot(md, field, **kwargs):
    from pyissm.plot import plot_model_field
    fig, ax, trip = plot_model_field(md, field, **kwargs)
    fig.canvas.draw()
    return trip


class TestLogNormLimits:

    def test_log_without_limits(self, md_and_field):
        trip = _plot(*md_and_field, log = True)
        assert isinstance(trip.norm, mcolors.LogNorm)
        assert (trip.norm.vmin, trip.norm.vmax) == (1., 1000.)

    def test_log_with_vmin_vmax(self, md_and_field):
        trip = _plot(*md_and_field, log = True, vmin = 10, vmax = 100)
        assert isinstance(trip.norm, mcolors.LogNorm)
        assert (trip.norm.vmin, trip.norm.vmax) == (10., 100.)

    def test_log_with_vmin_only(self, md_and_field):
        trip = _plot(*md_and_field, log = True, vmin = 10)
        assert (trip.norm.vmin, trip.norm.vmax) == (10., 1000.)

    @pytest.mark.parametrize('limits', [dict(vmin = 0, vmax = 100), dict(vmin = -1), dict(vmax = 0)])
    def test_log_rejects_non_positive_limits(self, md_and_field, limits):
        with pytest.raises(ValueError, match = 'must be positive'):
            _plot(*md_and_field, log = True, **limits)

    def test_log_and_norm_together_rejected(self, md_and_field):
        with pytest.raises(ValueError, match = 'either log = True or norm'):
            _plot(*md_and_field, log = True, norm = mcolors.LogNorm())


class TestLinearLimits:

    def test_linear_vmin_vmax(self, md_and_field):
        trip = _plot(*md_and_field, vmin = 5, vmax = 50)
        assert not isinstance(trip.norm, mcolors.LogNorm)
        assert (trip.norm.vmin, trip.norm.vmax) == (5., 50.)

    def test_user_norm(self, md_and_field):
        trip = _plot(*md_and_field, norm = mcolors.LogNorm(5, 50))
        assert isinstance(trip.norm, mcolors.LogNorm)
        assert (trip.norm.vmin, trip.norm.vmax) == (5., 50.)
