"""
Tests for pyissm.plot.animate_model_field (issues #40, #50: animate transient 2D results).
"""

import numpy as np
import pytest

matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')
import matplotlib.animation
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

NSTEPS = 4


@pytest.fixture
def md():
    """Square-mesh model with a fake 4-step TransientSolution (Vel on vertices, as ISSM stores it: (nv, 1))."""
    from pyissm.model import Model, mesh
    from pyissm.model.classes import results
    md = mesh.square_mesh(Model(), 1.e5, 1.e5, 6, 6)
    nv = md.mesh.numberofvertices
    steps = []
    for i in range(NSTEPS):
        step = results.solutionstep()
        step.time = 10. * i
        step.step = i + 1
        step.Vel = (1. + i) * np.linspace(1., 100., nv).reshape(-1, 1)
        step.IceVolume = 1.e12 + i
        steps.append(step)
    md.results.TransientSolution = results.solution(steps)
    return md


@pytest.fixture(autouse = True)
def _close_figures():
    yield
    plt.close('all')


def _frames(anim):
    return list(anim.new_frame_seq())


def _render(anim, frame):
    # The first canvas draw starts the animation, which redraws its first frame; do that before drawing `frame`
    anim._fig.canvas.draw()
    anim._func(frame)
    anim._fig.canvas.draw()


class TestAnimateModelField:

    def test_from_result_name(self, md):
        from pyissm.plot import animate_model_field
        fig, ax, anim = animate_model_field(md, 'Vel')
        assert isinstance(anim, matplotlib.animation.FuncAnimation)
        assert _frames(anim) == list(range(NSTEPS))
        assert ax.get_title() == 'Time = 0'
        _render(anim, 2)
        assert ax.get_title() == 'Time = 20'

    def test_colour_limits_span_all_frames(self, md):
        from pyissm.plot import animate_model_field
        fig, ax, anim = animate_model_field(md, 'Vel')
        trip = ax.collections[-1]
        assert (trip.norm.vmin, trip.norm.vmax) == (1., 100. * NSTEPS)
        _render(anim, NSTEPS - 1)
        trip = ax.collections[-1]
        assert (trip.norm.vmin, trip.norm.vmax) == (1., 100. * NSTEPS)

    def test_one_artist_per_frame(self, md):
        from pyissm.plot import animate_model_field
        fig, ax, anim = animate_model_field(md, 'Vel')
        n = len(ax.collections)
        for frame in range(NSTEPS):
            _render(anim, frame)
        assert len(ax.collections) == n

    def test_log_with_limits(self, md):
        from pyissm.plot import animate_model_field
        fig, ax, anim = animate_model_field(md, 'Vel', log = True, vmin = 10, vmax = 200)
        trip = ax.collections[-1]
        assert isinstance(trip.norm, mcolors.LogNorm)
        assert (trip.norm.vmin, trip.norm.vmax) == (10., 200.)

    def test_array_and_list_input(self, md):
        from pyissm.plot import animate_model_field
        vel = md.results.TransientSolution.Vel
        fig, ax, anim = animate_model_field(md, vel, time = [0., 1., 2., 3.], title = 'yr {time:.0f} step {step}')
        assert ax.get_title() == 'yr 0 step 0'
        fig, ax, anim = animate_model_field(md, [v[:, 0] for v in vel])
        assert ax.get_title() == 'Step 0'

    def test_frames_subset(self, md):
        from pyissm.plot import animate_model_field
        fig, ax, anim = animate_model_field(md, 'Vel', frames = [1, 3])
        assert _frames(anim) == [1, 3]
        assert ax.get_title() == 'Time = 10'

    def test_element_field(self, md):
        from pyissm.plot import animate_model_field
        ne = md.mesh.numberofelements
        data = np.vstack([np.full(ne, i + 1.) for i in range(NSTEPS)])
        fig, ax, anim = animate_model_field(md, data)
        _render(anim, NSTEPS - 1)

    def test_save_gif(self, md, tmp_path):
        from pyissm.plot import animate_model_field
        out = tmp_path / 'vel.gif'
        animate_model_field(md, 'Vel', save = str(out), fps = 2)
        assert out.stat().st_size > 0

    def test_no_title(self, md):
        from pyissm.plot import animate_model_field
        fig, ax, anim = animate_model_field(md, 'Vel', title = '', show_cbar = False)
        assert ax.get_title() == ''
        assert len(fig.axes) == 1


class TestAnimateModelFieldErrors:

    def test_unknown_result(self, md):
        from pyissm.plot import animate_model_field
        with pytest.raises(ValueError, match = "'Nope' not found"):
            animate_model_field(md, 'Nope')

    def test_unknown_solution(self, md):
        from pyissm.plot import animate_model_field
        with pytest.raises(ValueError, match = 'has no'):
            animate_model_field(md, 'Vel', solution = 'StressbalanceSolution')

    def test_scalar_series_rejected(self, md):
        from pyissm.plot import animate_model_field
        with pytest.raises(ValueError, match = 'one row per time step'):
            animate_model_field(md, 'IceVolume')

    def test_time_length_mismatch(self, md):
        from pyissm.plot import animate_model_field
        with pytest.raises(ValueError, match = 'time has 2 values'):
            animate_model_field(md, 'Vel', time = [0., 1.])

    def test_bad_frames(self, md):
        from pyissm.plot import animate_model_field
        with pytest.raises(ValueError, match = 'frames must be indices'):
            animate_model_field(md, 'Vel', frames = [0, NSTEPS])
