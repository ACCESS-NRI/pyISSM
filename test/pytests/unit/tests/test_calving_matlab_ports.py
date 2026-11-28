"""
Tests for the calving laws ported from MATLAB-only upstream classes:
calvingvonmisesAD (law 11), calvinghab (5), calvingdev2 (7), calvingtest (8), calvingpollard (10)
and calvingcalvingmip (12).

Law numbers and record formats follow $ISSM_DIR/src/m/classes/calving*.m and must match what
ISSM reads in LevelsetAnalysis.cpp (IoCodeToEnumCalving).
"""

import os
import tempfile
from types import SimpleNamespace

import numpy as np
import pytest

from pyissm.model.classes import calving, class_registry

YTS = 365.25 * 24 * 3600
CODES = {2: 'Integer', 3: 'Double', 6: 'IntMat', 7: 'DoubleMat'}


class _FakeModel(SimpleNamespace):
    """Minimal model: records consistency failures like Model.check_message."""

    def check_message(self, message):
        self.private.isconsistent = False
        self.messages.append(message)
        return self


def _fake_md(nv = 6, ne = 4, calving_obj = None):
    return _FakeModel(mesh = SimpleNamespace(numberofvertices = nv, numberofelements = ne),
                      constants = SimpleNamespace(yts = YTS),
                      transient = SimpleNamespace(ismovingfront = 1),
                      private = SimpleNamespace(isconsistent = True),
                      calving = calving_obj, messages = [])


def _records(obj, md):
    """Marshall obj and return {name: (format, mattype, values)}."""
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, 'out.bin')
        with open(path, 'wb') as fid:
            obj.marshall_class(fid, 'md.calving', md)
        buf = open(path, 'rb').read()
    out, pos = {}, 0
    while pos < len(buf):
        n = int(np.frombuffer(buf, np.int32, 1, pos)[0]); pos += 4
        name = buf[pos:pos + n].decode(); pos += n
        length = int(np.frombuffer(buf, np.int64, 1, pos)[0]); pos += 8
        payload = buf[pos:pos + length]; pos += length
        code = int(np.frombuffer(payload, np.int32, 1, 0)[0])
        fmt = CODES.get(code, code)
        if fmt == 'Integer':
            out[name] = (fmt, None, int(np.frombuffer(payload, np.int32, 1, 4)[0]))
        elif fmt == 'Double':
            out[name] = (fmt, None, float(np.frombuffer(payload, np.float64, 1, 4)[0]))
        else:
            mattype, rows, cols = (int(v) for v in np.frombuffer(payload, np.int32, 3, 4))
            out[name] = (fmt, mattype, np.frombuffer(payload, np.float64, rows * cols, 16).reshape(rows, cols))
    return out


def _check(md, solution = 'TransientSolution'):
    md.private.isconsistent = True
    md.calving.check_consistency(md, solution, [])
    return md.private.isconsistent


class TestRegistration:

    @pytest.mark.parametrize('name', ['vonmisesAD', 'hab', 'dev2', 'test', 'pollard', 'calvingmip'])
    def test_registered_and_listed(self, name):
        assert class_registry.CLASS_REGISTRY[f'calving.{name}'] is getattr(calving, name)
        assert name in calving.available()

    @pytest.mark.parametrize('matlab, python', [
        ('calvingvonmisesAD.calvingvonmisesAD', 'calving.vonmisesAD'), ('calvinghab.calvinghab', 'calving.hab'),
        ('calvingdev2.calvingdev2', 'calving.dev2'), ('calvingtest.calvingtest', 'calving.test'),
        ('calvingpollard.calvingpollard', 'calving.pollard'), ('calvingcalvingmip.calvingcalvingmip', 'calving.calvingmip'),
    ])
    def test_legacy_matlab_names(self, matlab, python):
        assert class_registry.map_classtype(matlab) == python


class TestDefaults:

    def test_defaults_match_upstream(self):
        ad = calving.vonmisesAD()
        assert np.isnan(ad.basin_id) and ad.num_basins == 0
        assert (ad.stress_threshold_groundedice, ad.stress_threshold_floatingice, ad.min_thickness) == (1e6, 150e3, 0.)
        assert calving.hab().flotation_fraction == 0.15
        d = calving.dev2()
        assert (d.stress_threshold_groundedice, d.stress_threshold_floatingice, d.height_above_floatation) == (1e6, 150e3, 0.)
        t = calving.test()
        assert (t.speedfactor, t.independentrate) == (1, 0)
        assert calving.pollard().rc == 0.75
        m = calving.calvingmip()
        assert (m.min_thickness, m.experiment) == (0., 1)

    @pytest.mark.parametrize('cls', [calving.vonmisesAD, calving.hab, calving.dev2, calving.test, calving.pollard, calving.calvingmip])
    def test_repr_and_str(self, cls):
        obj = cls()
        assert str(obj) == f'ISSM - calving.{cls.__name__} Class'
        assert repr(obj)


class TestMarshall:

    def test_vonmisesAD(self):
        md = _fake_md()
        c = calving.vonmisesAD()
        c.num_basins = 2
        c.basin_id = np.array([1, 2, 2, 1])
        c.stress_threshold_groundedice = np.array([1e6, 8e5])
        c.stress_threshold_floatingice = np.array([150e3, 2e5])
        c.min_thickness = 10.
        r = _records(c, md)
        assert r['md.calving.law'] == ('Integer', None, 11)
        assert r['md.calving.num_basins'] == ('Integer', None, 2)
        fmt, mattype, basin = r['md.calving.basin_id']
        assert (fmt, mattype) == ('IntMat', 2)
        np.testing.assert_array_equal(basin.ravel(), [0, 1, 1, 0])  # 0-indexed
        assert r['md.calving.stress_threshold_groundedice'][0] == 'DoubleMat'
        np.testing.assert_array_equal(r['md.calving.stress_threshold_groundedice'][2].ravel(), [1e6, 8e5])
        np.testing.assert_array_equal(r['md.calving.stress_threshold_floatingice'][2].ravel(), [150e3, 2e5])
        assert r['md.calving.min_thickness'] == ('Double', None, 10.)

    def test_hab(self):
        r = _records(calving.hab(), _fake_md())
        assert r['md.calving.law'] == ('Integer', None, 5)
        fmt, mattype, value = r['md.calving.flotation_fraction']
        assert (fmt, mattype, value.item()) == ('DoubleMat', 1, 0.15)

    def test_dev2(self):
        md = _fake_md()
        c = calving.dev2()
        c.height_above_floatation = 5.
        r = _records(c, md)
        assert r['md.calving.law'] == ('Integer', None, 7)
        assert r['md.calving.stress_threshold_groundedice'][:2] == ('DoubleMat', 1)
        assert r['md.calving.stress_threshold_floatingice'][:2] == ('DoubleMat', 1)
        assert r['md.calving.height_above_floatation'] == ('Double', None, 5.)

    def test_test_scalar_and_rate_scaling(self):
        c = calving.test()
        c.speedfactor = 0.5
        c.independentrate = 100.  # m/yr
        r = _records(c, _fake_md())
        assert r['md.calving.law'] == ('Integer', None, 8)
        assert r['md.calving.speedfactor'][2].item() == 0.5
        assert r['md.calving.independentrate'][2].item() == pytest.approx(100. / YTS)  # written in m/s

    def test_test_time_series(self):
        c = calving.test()
        c.independentrate = np.array([[0., 100.], [0., 10.]])  # values row, then times [yr]
        r = _records(c, _fake_md())
        written = r['md.calving.independentrate'][2]
        np.testing.assert_allclose(written[0], [0., 100. / YTS])
        np.testing.assert_allclose(written[1], [0., 10. * YTS])  # times converted to seconds

    def test_pollard(self):
        r = _records(calving.pollard(), _fake_md())
        assert r == {'md.calving.law': ('Integer', None, 10), 'md.calving.rc': ('Double', None, 0.75)}

    def test_calvingmip(self):
        c = calving.calvingmip()
        c.experiment = 3
        r = _records(c, _fake_md())
        assert r == {'md.calving.law': ('Integer', None, 12), 'md.calving.min_thickness': ('Double', None, 0.),
                     'md.calving.experiment': ('Integer', None, 3)}


class TestCheckConsistency:

    def test_skipped_outside_moving_front_transient(self):
        md = _fake_md(calving_obj = calving.pollard())
        md.calving.rc = 2.
        assert _check(md, 'StressbalanceSolution')
        md.transient.ismovingfront = 0
        assert _check(md)

    def test_pollard_rc_range(self):
        md = _fake_md(calving_obj = calving.pollard())
        assert _check(md)
        for bad in (0., 1.5):
            md.calving.rc = bad
            assert not _check(md)

    def test_calvingmip_experiment_values(self):
        md = _fake_md(calving_obj = calving.calvingmip())
        assert _check(md)
        md.calving.experiment = 6
        assert not _check(md)

    def test_hab_fraction(self):
        md = _fake_md(calving_obj = calving.hab())
        assert _check(md)
        md.calving.flotation_fraction = -0.1
        assert not _check(md)

    def test_dev2_thresholds(self):
        md = _fake_md(calving_obj = calving.dev2())
        assert _check(md)
        md.calving.stress_threshold_floatingice = 0.
        assert not _check(md)

    def test_test_rates(self):
        md = _fake_md(calving_obj = calving.test())
        assert _check(md)
        md.calving.speedfactor = -1.
        assert not _check(md)

    def test_vonmisesAD_basins(self):
        md = _fake_md(calving_obj = calving.vonmisesAD())
        md.calving.num_basins = 2
        md.calving.basin_id = np.array([1, 2, 2, 1])
        md.calving.stress_threshold_groundedice = np.array([1e6, 8e5])
        md.calving.stress_threshold_floatingice = np.array([150e3, 2e5])
        assert _check(md)
        md.calving.basin_id = np.array([1, 3, 2, 1])  # basin 3 > num_basins
        assert not _check(md)
        md.calving.basin_id = np.array([1, 2, 2, 1])
        md.calving.stress_threshold_floatingice = np.array([150e3])  # one threshold for two basins
        assert not _check(md)


def _requires_wrappers():
    from pyissm.tools import wrappers
    try:
        return not bool(wrappers.check_wrappers_installed())
    except Exception:
        return True


@pytest.mark.skipif(_requires_wrappers(), reason = 'ISSM Python wrappers not installed')
class TestExtrudeCollapse:

    def _model(self, calving_obj):
        from pyissm.model import Model, mesh
        md = mesh.square_mesh(Model(), 1.e5, 1.e5, 5, 5)
        nv, ne = md.mesh.numberofvertices, md.mesh.numberofelements
        md.geometry.thickness = 1000. * np.ones(nv)
        md.geometry.base = -500. * np.ones(nv)
        md.geometry.surface = md.geometry.base + md.geometry.thickness
        md.geometry.bed = -600. * np.ones(nv)
        md.mask.ice_levelset = -np.ones(nv)
        md.mask.ocean_levelset = np.ones(nv)
        md.materials.rheology_B = 1.e8 * np.ones(nv)
        md.materials.rheology_n = 3. * np.ones(ne)
        md.calving = calving_obj(nv, ne)
        return md, nv, ne

    def test_hab_vertex_field(self):
        def make(nv, ne):
            c = calving.hab()
            c.flotation_fraction = 0.1 + 0.01 * np.arange(nv)
            return c
        md, nv, ne = self._model(make)
        md = md.extrude(3, 1.)
        assert np.size(md.calving.flotation_fraction) == md.mesh.numberofvertices
        md = md.collapse()
        np.testing.assert_allclose(md.calving.flotation_fraction, 0.1 + 0.01 * np.arange(nv))

    def test_dev2_scalar_and_vertex_thresholds(self):
        def make(nv, ne):
            c = calving.dev2()
            c.stress_threshold_floatingice = 1e5 + np.arange(nv)
            return c
        md, nv, ne = self._model(make)
        md = md.extrude(3, 1.)
        assert md.calving.stress_threshold_groundedice == 1e6
        assert np.size(md.calving.stress_threshold_floatingice) == md.mesh.numberofvertices
        md = md.collapse()
        np.testing.assert_allclose(md.calving.stress_threshold_floatingice, 1e5 + np.arange(nv))

    def test_vonmisesAD_basin_id(self):
        def make(nv, ne):
            c = calving.vonmisesAD()
            c.num_basins = 2
            c.basin_id = 1 + (np.arange(ne) % 2)
            c.stress_threshold_groundedice = np.array([1e6, 8e5])
            c.stress_threshold_floatingice = np.array([150e3, 2e5])
            return c
        md, nv, ne = self._model(make)
        md = md.extrude(3, 1.)
        assert np.size(md.calving.basin_id) == md.mesh.numberofelements
        np.testing.assert_array_equal(md.calving.stress_threshold_groundedice, [1e6, 8e5])
        md = md.collapse()
        np.testing.assert_array_equal(md.calving.basin_id, 1 + (np.arange(ne) % 2))

    @pytest.mark.parametrize('cls', [calving.test, calving.pollard, calving.calvingmip])
    def test_scalar_laws_unchanged(self, cls):
        md, nv, ne = self._model(lambda nv, ne: cls())
        before = dict(vars(md.calving))
        md = md.extrude(3, 1.)
        assert vars(md.calving) == before
