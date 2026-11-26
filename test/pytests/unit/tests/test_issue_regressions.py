"""
Regression tests for fixed GitHub issues:

- #141: friction.coulomb2 must implement ISSM friction law 13 (C, m), not duplicate friction.coulomb (law 7).
- #224: Model.collapse() must handle material classes without rheology_n (estar) and collapse
  rheology_E / rheology_Ec / rheology_Es.
- #128: a dependent with fos_reverse_index must not require nods (MATLAB scripts omit it; it is never marshalled).
- #170: load_results_from_cluster() must download from the runtime_name it is given, and
  load_results_from_disk() must return the model when the output file is missing.
"""

import os
import tempfile

import numpy as np
import pytest
from types import SimpleNamespace

from pyissm.model import execute
from pyissm.model.classes import friction

from pyissm.tools import wrappers


def _wrappers_installed():
    try:
        return bool(wrappers.check_wrappers_installed())
    except Exception:
        return False


# Mesh construction and collapse() call the compiled ISSM wrappers (NodeConnectivity, ElementConnectivity)
requires_wrappers = pytest.mark.skipif(not _wrappers_installed(), reason = 'ISSM Python wrappers not installed')


def _marshall_records(obj, md, prefix):
    """Marshall obj to a temp file and return {record name: payload bytes}."""
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, 'out.bin')
        with open(path, 'wb') as fid:
            obj.marshall_class(fid, prefix, md)
        buf = open(path, 'rb').read()
    records, pos = {}, 0
    while pos < len(buf):
        n = int(np.frombuffer(buf, np.int32, 1, pos)[0]); pos += 4
        name = buf[pos:pos + n].decode(); pos += n
        length = int(np.frombuffer(buf, np.int64, 1, pos)[0]); pos += 8
        records[name] = buf[pos:pos + length]; pos += length
    return records


def _integer_record(payload):
    # Integer record payload: int32 format code, then int32 value
    return int(np.frombuffer(payload, np.int32, 1, 4)[0])


def _fake_md(nv = 4, ne = 2):
    return SimpleNamespace(mesh = SimpleNamespace(numberofvertices = nv, numberofelements = ne),
                           constants = SimpleNamespace(yts = 365.25 * 24 * 3600))


class TestCoulomb2Law13:
    """friction.coulomb2 is the Coulomb-limited law 13, with fields C and m (issue #141)."""

    def test_fields(self):
        f = friction.coulomb2()
        for field in ['C', 'm', 'coupling', 'effective_pressure', 'effective_pressure_limit']:
            assert hasattr(f, field)
        for field in ['coefficient', 'coefficientcoulomb', 'p', 'q']:
            assert not hasattr(f, field)

    def test_marshall_writes_law_13(self):
        md = _fake_md()
        f = friction.coulomb2()
        f.C = np.ones(md.mesh.numberofvertices)
        f.m = np.ones(md.mesh.numberofelements) / 3.
        records = _marshall_records(f, md, 'md.friction')
        assert _integer_record(records['md.friction.law']) == 13
        assert {'md.friction.C', 'md.friction.m', 'md.friction.coupling', 'md.friction.effective_pressure_limit'} <= set(records)
        assert 'md.friction.effective_pressure' not in records

    @pytest.mark.parametrize('coupling', [3, 4])
    def test_marshall_effective_pressure_for_coupling_3_4(self, coupling):
        md = _fake_md()
        f = friction.coulomb2()
        f.C = np.ones(md.mesh.numberofvertices)
        f.m = np.ones(md.mesh.numberofelements) / 3.
        f.coupling = coupling
        f.effective_pressure = np.ones(md.mesh.numberofvertices)
        records = _marshall_records(f, md, 'md.friction')
        assert 'md.friction.effective_pressure' in records

    def test_differs_from_coulomb(self):
        md = _fake_md()
        c = friction.coulomb()
        c.coefficient = np.ones(md.mesh.numberofvertices)
        c.coefficientcoulomb = np.ones(md.mesh.numberofvertices)
        c.p = np.ones(md.mesh.numberofelements)
        c.q = np.ones(md.mesh.numberofelements)
        assert _integer_record(_marshall_records(c, md, 'md.friction')['md.friction.law']) == 7


class TestDependentNods:
    """A fos_reverse dependent does not need nods (issue #128)."""

    def test_fos_reverse_dependent_without_nods(self):
        from pyissm.model.classes import dependent
        dep = dependent()
        dep.name = 'Outputdefinition1'
        dep.type = 'scalar'
        dep.fos_reverse_index = 1
        md = SimpleNamespace()
        assert dep.check_consistency(md, 'TransientSolution', []) is md


class TestLoadResults:
    """Result loading respects runtime_name and survives a missing output file (issue #170)."""

    def test_load_results_from_disk_missing_file_returns_md(self, tmp_path):
        md = SimpleNamespace(qmu = SimpleNamespace(isdakota = 0))
        assert execute.load_results_from_disk(md, str(tmp_path / 'missing.outbin')) is md

    def test_load_results_from_cluster_uses_runtime_name(self, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        downloaded = []
        cluster = SimpleNamespace(name = 'not-this-host', download = lambda d, files: downloaded.append(d))
        md = SimpleNamespace(cluster = cluster,
                             qmu = SimpleNamespace(isdakota = 0),
                             miscellaneous = SimpleNamespace(name = 'model'),
                             private = SimpleNamespace(runtimename = 'stale-runtime'))
        monkeypatch.setattr(execute, 'load_results_from_disk', lambda md, filename: md)
        with pytest.warns(UserWarning):
            execute.load_results_from_cluster(md, runtime_name = 'model-11-24-2026-10-00-00-1234')
        assert downloaded == ['model-11-24-2026-10-00-00-1234']

    def test_load_results_from_cluster_defaults_to_private_runtimename(self, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        downloaded = []
        cluster = SimpleNamespace(name = 'not-this-host', download = lambda d, files: downloaded.append(d))
        md = SimpleNamespace(cluster = cluster,
                             qmu = SimpleNamespace(isdakota = 0),
                             miscellaneous = SimpleNamespace(name = 'model'),
                             private = SimpleNamespace(runtimename = 'model-runtime'))
        monkeypatch.setattr(execute, 'load_results_from_disk', lambda md, filename: md)
        with pytest.warns(UserWarning):
            execute.load_results_from_cluster(md)
        assert downloaded == ['model-runtime']


def _model_3d(materials, friction_obj = None, numlayers = 3):
    """Small extruded square model with the given (2D-sized) materials and friction."""
    from pyissm.model import Model, mesh
    md = mesh.square_mesh(Model(), 1.e5, 1.e5, 5, 5)
    nv, ne = md.mesh.numberofvertices, md.mesh.numberofelements
    md.geometry.thickness = 1000. * np.ones(nv)
    md.geometry.base = -500. * np.ones(nv)
    md.geometry.surface = md.geometry.base + md.geometry.thickness
    md.geometry.bed = -600. * np.ones(nv)
    md.mask.ice_levelset = -np.ones(nv)
    md.mask.ocean_levelset = np.ones(nv)
    md.materials = materials(nv, ne)
    if friction_obj is not None:
        md.friction = friction_obj(nv, ne)
    md = md.extrude(numlayers, 1.)
    return md, nv, ne


@requires_wrappers
class TestCollapseMaterials:
    """Model.collapse() works for every ice material class (issue #224)."""

    def test_collapse_estar(self):
        from pyissm.model.classes import materials

        def make(nv, ne):
            m = materials.estar()
            m.rheology_B = 1.e8 * np.ones(nv)
            m.rheology_Ec = 1. + np.arange(nv, dtype = float)
            m.rheology_Es = 2. + np.arange(nv, dtype = float)
            return m

        md, nv, ne = _model_3d(make)
        assert np.size(md.materials.rheology_Ec) == md.mesh.numberofvertices
        md = md.collapse()
        assert md.mesh.numberofvertices == nv
        assert not hasattr(md.materials, 'rheology_n')
        np.testing.assert_array_equal(md.materials.rheology_Ec, 1. + np.arange(nv))
        np.testing.assert_array_equal(md.materials.rheology_Es, 2. + np.arange(nv))
        assert np.size(md.materials.rheology_B) == nv

    def test_collapse_enhancedice(self):
        from pyissm.model.classes import materials

        def make(nv, ne):
            m = materials.enhancedice()
            m.rheology_B = 1.e8 * np.ones(nv)
            m.rheology_E = 1. + np.arange(nv, dtype = float)
            m.rheology_n = 3. * np.ones(ne)
            return m

        md, nv, ne = _model_3d(make)
        md = md.collapse()
        np.testing.assert_array_equal(md.materials.rheology_E, 1. + np.arange(nv))
        np.testing.assert_array_equal(md.materials.rheology_n, 3. * np.ones(ne))

    def test_collapse_coulomb2(self):
        from pyissm.model.classes import materials

        def make_materials(nv, ne):
            m = materials.ice()
            m.rheology_B = 1.e8 * np.ones(nv)
            m.rheology_n = 3. * np.ones(ne)
            return m

        def make_friction(nv, ne):
            f = friction.coulomb2()
            f.C = 1. + np.arange(nv, dtype = float)
            f.m = np.ones(ne) / 3.
            return f

        md, nv, ne = _model_3d(make_materials, make_friction)
        assert np.size(md.friction.C) == md.mesh.numberofvertices
        md = md.collapse()
        np.testing.assert_array_equal(md.friction.C, 1. + np.arange(nv))
        np.testing.assert_array_equal(md.friction.m, np.ones(ne) / 3.)
