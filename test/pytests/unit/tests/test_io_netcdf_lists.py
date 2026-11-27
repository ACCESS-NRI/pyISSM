"""
save_model / load_model round trip for list fields.

Lists of numbers (e.g. md.inversion.cost_functions = [101]) used to be written as strings and came
back as ['101'], which then failed validation when marshalling the reloaded model.
"""

import numpy as np
import pytest

netCDF4 = pytest.importorskip('netCDF4')

from pyissm.model import Model, io
from pyissm.model.classes import class_utils


def _round_trip(md, tmp_path):
    path = str(tmp_path / 'model.nc')
    io.save_model(md, path)
    return io.load_model(path), path


class TestNumericLists:

    def test_default_cost_functions(self, tmp_path):
        md = Model()
        assert md.inversion.cost_functions == [101]
        md2, _ = _round_trip(md, tmp_path)
        assert md2.inversion.cost_functions == [101]
        assert type(md2.inversion.cost_functions[0]) is int

    def test_reloaded_cost_functions_marshall(self, tmp_path):
        md = Model()
        md.inversion.cost_functions = [101, 103, 501]
        md2, _ = _round_trip(md, tmp_path)
        assert class_utils.marshall_inversion_cost_functions(md2.inversion.cost_functions) == \
            ['SurfaceAbsVelMisfit', 'SurfaceLogVelMisfit', 'DragCoefficientAbsGradient']

    @pytest.mark.parametrize('values, kind', [([1, 2, 3], int), ([1.5, 2.], float), ([1, 2.5], float), ([True, False], bool)])
    def test_list_types(self, tmp_path, values, kind):
        md = Model()
        md.miscellaneous.dummy = {'values': values}
        md2, _ = _round_trip(md, tmp_path)
        loaded = md2.miscellaneous.dummy['values']
        assert isinstance(loaded, list)
        assert loaded == values
        assert all(type(v) is kind for v in loaded)

    def test_stored_as_numbers(self, tmp_path):
        md = Model()
        md.inversion.cost_functions = [101, 501]
        _, path = _round_trip(md, tmp_path)
        with netCDF4.Dataset(path) as ds:
            var = ds['inversion'].variables['cost_functions']
            assert var.dtype.kind == 'i'
            assert var.getncattr('__LIST_TYPE__') == 'int'

    def test_string_lists_unchanged(self, tmp_path):
        md = Model()
        md.transient.requested_outputs = ['default', 'IceVolume']
        md2, path = _round_trip(md, tmp_path)
        assert md2.transient.requested_outputs == ['default', 'IceVolume']
        with netCDF4.Dataset(path) as ds:
            assert '__LIST_TYPE__' not in ds['transient'].variables['requested_outputs'].ncattrs()

    def test_numeric_arrays_stay_arrays(self, tmp_path):
        md = Model()
        md.inversion.cost_functions_coefficients = np.arange(8.).reshape(4, 2)
        md2, _ = _round_trip(md, tmp_path)
        assert isinstance(md2.inversion.cost_functions_coefficients, np.ndarray)
        np.testing.assert_array_equal(md2.inversion.cost_functions_coefficients, np.arange(8.).reshape(4, 2))

    def test_old_files_still_load(self, tmp_path, monkeypatch):
        """Files written before this change hold numeric lists as untagged strings; they load as before."""
        monkeypatch.setattr(io, '_numeric_list_type', lambda values: None)  # old behaviour: every list as strings
        md = Model()
        _, path = _round_trip(md, tmp_path)
        monkeypatch.undo()
        with netCDF4.Dataset(path) as ds:
            assert '__LIST_TYPE__' not in ds['inversion'].variables['cost_functions'].ncattrs()
        assert io.load_model(path).inversion.cost_functions == ['101']


class TestNumericListType:

    @pytest.mark.parametrize('values, expected', [
        ([1, 2], 'int'), ([np.int64(1)], 'int'), ([1.5, 2], 'float'), ([np.float32(1.)], 'float'),
        ([True], 'bool'), ([np.bool_(False), True], 'bool'),
        ([True, 1], None), (['a'], None), ([1, 'a'], None), ([{'a': 1}], None), ([], None),
    ])
    def test_cases(self, values, expected):
        assert io._numeric_list_type(values) == expected
