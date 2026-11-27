"""
Tests for pyissm.tools.general.compare_models (issue #67).
"""

import collections

import numpy as np
import pytest

from pyissm.model import Model
from pyissm.model.classes import friction
from pyissm.tools.general import compare_models


def _model(nv = 10):
    md = Model()
    md.geometry.thickness = 1000. * np.ones(nv)
    md.geometry.base = -500. * np.ones(nv)
    md.geometry.surface = md.geometry.base + md.geometry.thickness
    md.friction.coefficient = 30. * np.ones(nv)
    md.miscellaneous.name = 'test'
    return md


def _paths(differences):
    return [path for path, _ in differences]


class TestCompareModels:

    def test_identical(self, capsys):
        assert compare_models(_model(), _model()) == []
        assert 'no differences found' in capsys.readouterr().out

    def test_default_models_identical(self):
        assert compare_models(Model(), Model(), verbose = False) == []

    def test_value_difference(self, capsys):
        md2 = _model()
        md2.geometry.thickness[:3] += 0.5
        diffs = compare_models(_model(), md2)
        assert diffs == [('geometry.thickness', 'values differ (3 of 10 entries, max abs diff 0.5)')]
        assert 'geometry.thickness: values differ' in capsys.readouterr().out

    def test_tolerance(self):
        md2 = _model()
        md2.geometry.thickness *= 1 + 1e-12
        assert _paths(compare_models(_model(), md2, verbose = False)) == ['geometry.thickness']
        assert compare_models(_model(), md2, rtol = 1e-10, verbose = False) == []
        assert compare_models(_model(), md2, atol = 1e-6, verbose = False) == []

    def test_nan_equal_to_nan(self):
        md1, md2 = _model(), _model()
        md1.geometry.bed = np.full(10, np.nan)
        md2.geometry.bed = np.full(10, np.nan)
        md1.geometry.thickness[2] = md2.geometry.thickness[2] = np.nan
        assert compare_models(md1, md2, verbose = False) == []

    def test_nan_positions_differ(self):
        md2 = _model()
        md2.geometry.thickness[0] = np.nan
        diffs = compare_models(_model(), md2, verbose = False)
        assert diffs == [('geometry.thickness', 'NaN positions differ (0 vs 1 NaN)')]

    def test_shape_difference(self):
        md2 = _model()
        md2.geometry.thickness = np.ones(12)
        (path, description), = compare_models(_model(), md2, fields = ['geometry'], ignore = ['geometry.surface'], verbose = False)
        assert path == 'geometry.thickness' and description.startswith('shape differs')

    def test_column_vectors_match_by_default(self):
        md2 = _model()
        md2.geometry.thickness = md2.geometry.thickness.reshape(-1, 1)
        assert compare_models(_model(), md2, verbose = False) == []
        diffs = compare_models(_model(), md2, ignore_singleton_dims = False, verbose = False)
        assert diffs == [('geometry.thickness', 'shape differs ((10,) vs (10, 1))')]

    def test_class_difference_not_recursed(self):
        md2 = _model()
        md2.friction = friction.schoof()
        diffs = compare_models(_model(), md2, verbose = False)
        assert diffs == [('friction', 'class differs (friction.default vs friction.schoof), not compared further')]

    def test_string_and_type_differences(self):
        md2 = _model()
        md2.miscellaneous.name = 'other'
        md2.inversion.cost_functions = ['101']
        diffs = dict(compare_models(_model(), md2, verbose = False))
        assert diffs['miscellaneous.name'] == "differs ('test' vs 'other')"
        assert diffs['inversion.cost_functions[0]'] == 'type differs (int vs str)'

    def test_list_length_difference(self):
        md1, md2 = _model(), _model()
        # Numeric lists compare as arrays; other lists element by element
        md2.inversion.cost_functions = [101, 501]
        assert dict(compare_models(md1, md2, verbose = False)) == {'inversion.cost_functions': 'shape differs ((1,) vs (2,))'}
        md1.transient.requested_outputs = ['default']
        md2.transient.requested_outputs = ['default', 'IceVolume']
        assert dict(compare_models(md1, md2, fields = ['transient'], verbose = False)) == {'transient.requested_outputs': 'length differs (1 vs 2)'}

    def test_missing_field(self):
        md2 = _model()
        md2.geometry.extra = 1.
        assert compare_models(_model(), md2, verbose = False) == [('geometry.extra', 'only in md2')]
        assert compare_models(md2, _model(), verbose = False) == [('geometry.extra', 'only in md1')]

    def test_dict_and_ordereddict_compare_by_content(self):
        md1, md2 = _model(), _model()
        md1.miscellaneous.dummy = collections.OrderedDict(a = 1., b = 'x')
        md2.miscellaneous.dummy = {'a': 1., 'b': 'x'}
        assert compare_models(md1, md2, verbose = False) == []
        md2.miscellaneous.dummy['a'] = 2.
        assert _paths(compare_models(md1, md2, verbose = False)) == ["miscellaneous.dummy['a']"]

    def test_bool_vs_int(self):
        md1, md2 = _model(), _model()
        md1.transient.isstressbalance = True
        md2.transient.isstressbalance = 1
        assert compare_models(md1, md2, verbose = False) == []

    def test_fields_and_ignore(self):
        md2 = _model()
        md2.geometry.thickness[0] += 1.
        md2.friction.coefficient[0] += 1.
        assert _paths(compare_models(_model(), md2, fields = ['friction'], verbose = False)) == ['friction.coefficient']
        assert _paths(compare_models(_model(), md2, ignore = ['geometry'], verbose = False)) == ['friction.coefficient']
        assert compare_models(_model(), md2, ignore = ['geometry.thickness', 'friction'], verbose = False) == []

    def test_nested_results(self):
        from pyissm.model.classes import results
        md1, md2 = _model(), _model()
        for md, vel in ((md1, 1.), (md2, 2.)):
            step = results.solutionstep()
            step.Vel = vel * np.ones(10)
            md.results.TransientSolution = results.solution([step])
        assert _paths(compare_models(md1, md2, verbose = False)) == ['results.TransientSolution.steps[0].Vel']
        assert compare_models(md1, md2, ignore = ['results'], verbose = False) == []

    def test_netcdf_paths(self, tmp_path):
        from pyissm.model import io
        md = _model()
        io.save_model(md, str(tmp_path / 'a.nc'))
        io.save_model(md, str(tmp_path / 'b.nc'))
        assert compare_models(str(tmp_path / 'a.nc'), str(tmp_path / 'b.nc'), verbose = False) == []
