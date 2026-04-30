"""
Unit tests for pyissm.model.execute binary write functions and format_to_code.
Tests use io.BytesIO as an in-memory file to avoid disk I/O.
"""

import pytest
import io
import struct
import numpy as np
from types import SimpleNamespace

try:
    from pyissm.model.execute import (
        format_to_code,
        _write_field_name,
        _write_boolean,
        _write_integer,
        _write_double,
        _write_string,
        _write_string_array,
        _apply_scaling_and_time_series,
    )
    EXECUTE_AVAILABLE = True
except ImportError:
    EXECUTE_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not EXECUTE_AVAILABLE,
    reason="pyissm.model.execute not available"
)


# ============== FORMAT_TO_CODE ==============

class TestFormatToCode:
    """Tests for format_to_code()."""

    def test_boolean(self):
        assert format_to_code('Boolean') == 1

    def test_integer(self):
        assert format_to_code('Integer') == 2

    def test_double(self):
        assert format_to_code('Double') == 3

    def test_string(self):
        assert format_to_code('String') == 4

    def test_boolean_mat(self):
        assert format_to_code('BooleanMat') == 5

    def test_int_mat(self):
        assert format_to_code('IntMat') == 6

    def test_double_mat(self):
        assert format_to_code('DoubleMat') == 7

    def test_mat_array(self):
        assert format_to_code('MatArray') == 8

    def test_string_array(self):
        assert format_to_code('StringArray') == 9

    def test_compressed_mat(self):
        assert format_to_code('CompressedMat') == 10

    def test_unsupported_raises_ioerror(self):
        with pytest.raises(IOError):
            format_to_code('Unsupported')

    def test_case_sensitive(self):
        """format_to_code should be case-sensitive."""
        with pytest.raises(IOError):
            format_to_code('double')

    def test_empty_string_raises(self):
        with pytest.raises(IOError):
            format_to_code('')


# ============== _WRITE_FIELD_NAME ==============

class TestWriteFieldName:
    """Tests for _write_field_name()."""

    def test_writes_name_length_as_int32(self):
        fid = io.BytesIO()
        _write_field_name(fid, 'hello')
        fid.seek(0)
        length = struct.unpack('<i', fid.read(4))[0]
        assert length == 5  # len('hello')

    def test_writes_name_bytes(self):
        fid = io.BytesIO()
        _write_field_name(fid, 'test')
        fid.seek(4)
        assert fid.read(4) == b'test'

    def test_empty_name(self):
        fid = io.BytesIO()
        _write_field_name(fid, '')
        fid.seek(0)
        length = struct.unpack('<i', fid.read(4))[0]
        assert length == 0

    def test_total_bytes_written(self):
        fid = io.BytesIO()
        name = 'myfield'
        _write_field_name(fid, name)
        assert fid.tell() == 4 + len(name)

    def test_unicode_ascii(self):
        fid = io.BytesIO()
        _write_field_name(fid, 'abc')
        fid.seek(0)
        length = struct.unpack('<i', fid.read(4))[0]
        assert length == 3
        assert fid.read(3) == b'abc'


# ============== _WRITE_BOOLEAN ==============

class TestWriteBoolean:
    """Tests for _write_boolean()."""

    def test_writes_correct_record_length(self):
        fid = io.BytesIO()
        _write_boolean(fid, True, 'flag')
        fid.seek(0)
        record_length = struct.unpack('<q', fid.read(8))[0]
        assert record_length == 8  # 4 bytes code + 4 bytes data

    def test_writes_correct_format_code(self):
        fid = io.BytesIO()
        _write_boolean(fid, True, 'flag')
        fid.seek(8)
        code = struct.unpack('<i', fid.read(4))[0]
        assert code == format_to_code('Boolean')

    def test_true_written_as_1(self):
        fid = io.BytesIO()
        _write_boolean(fid, True, 'flag')
        fid.seek(12)
        value = struct.unpack('<i', fid.read(4))[0]
        assert value == 1

    def test_false_written_as_0(self):
        fid = io.BytesIO()
        _write_boolean(fid, False, 'flag')
        fid.seek(12)
        value = struct.unpack('<i', fid.read(4))[0]
        assert value == 0

    def test_total_bytes(self):
        fid = io.BytesIO()
        _write_boolean(fid, True, 'flag')
        assert fid.tell() == 8 + 4 + 4


# ============== _WRITE_INTEGER ==============

class TestWriteInteger:
    """Tests for _write_integer()."""

    def test_writes_correct_record_length(self):
        fid = io.BytesIO()
        _write_integer(fid, 42, 'num')
        fid.seek(0)
        record_length = struct.unpack('<q', fid.read(8))[0]
        assert record_length == 8

    def test_writes_correct_format_code(self):
        fid = io.BytesIO()
        _write_integer(fid, 42, 'num')
        fid.seek(8)
        code = struct.unpack('<i', fid.read(4))[0]
        assert code == format_to_code('Integer')

    def test_positive_value(self):
        fid = io.BytesIO()
        _write_integer(fid, 42, 'num')
        fid.seek(12)
        value = struct.unpack('<i', fid.read(4))[0]
        assert value == 42

    def test_zero(self):
        fid = io.BytesIO()
        _write_integer(fid, 0, 'num')
        fid.seek(12)
        value = struct.unpack('<i', fid.read(4))[0]
        assert value == 0

    def test_negative_value(self):
        fid = io.BytesIO()
        _write_integer(fid, -7, 'num')
        fid.seek(12)
        value = struct.unpack('<i', fid.read(4))[0]
        assert value == -7


# ============== _WRITE_DOUBLE ==============

class TestWriteDouble:
    """Tests for _write_double()."""

    def test_writes_correct_record_length(self):
        fid = io.BytesIO()
        _write_double(fid, 3.14, 'val')
        fid.seek(0)
        record_length = struct.unpack('<q', fid.read(8))[0]
        assert record_length == 12  # 4 bytes code + 8 bytes float64

    def test_writes_correct_format_code(self):
        fid = io.BytesIO()
        _write_double(fid, 3.14, 'val')
        fid.seek(8)
        code = struct.unpack('<i', fid.read(4))[0]
        assert code == format_to_code('Double')

    def test_writes_float_value(self):
        fid = io.BytesIO()
        _write_double(fid, 1.5, 'val')
        fid.seek(12)
        value = struct.unpack('<d', fid.read(8))[0]
        assert np.isclose(value, 1.5)

    def test_writes_pi(self):
        fid = io.BytesIO()
        _write_double(fid, np.pi, 'pi')
        fid.seek(12)
        value = struct.unpack('<d', fid.read(8))[0]
        assert np.isclose(value, np.pi)

    def test_writes_zero(self):
        fid = io.BytesIO()
        _write_double(fid, 0.0, 'zero')
        fid.seek(12)
        value = struct.unpack('<d', fid.read(8))[0]
        assert value == 0.0

    def test_total_bytes(self):
        fid = io.BytesIO()
        _write_double(fid, 1.0, 'val')
        assert fid.tell() == 8 + 4 + 8


# ============== _WRITE_STRING ==============

class TestWriteString:
    """Tests for _write_string()."""

    def test_writes_correct_record_length(self):
        fid = io.BytesIO()
        s = 'hello'
        _write_string(fid, s)
        fid.seek(0)
        record_length = struct.unpack('<q', fid.read(8))[0]
        # len(s) + 4 (code) + 4 (str len field)
        assert record_length == len(s) + 8

    def test_writes_format_code(self):
        fid = io.BytesIO()
        _write_string(fid, 'test')
        fid.seek(8)
        code = struct.unpack('<i', fid.read(4))[0]
        assert code == format_to_code('String')

    def test_writes_string_length(self):
        fid = io.BytesIO()
        _write_string(fid, 'hello')
        fid.seek(12)
        str_len = struct.unpack('<i', fid.read(4))[0]
        assert str_len == 5

    def test_writes_string_bytes(self):
        fid = io.BytesIO()
        _write_string(fid, 'hello')
        fid.seek(16)
        assert fid.read(5) == b'hello'

    def test_empty_string(self):
        fid = io.BytesIO()
        _write_string(fid, '')
        fid.seek(0)
        record_length = struct.unpack('<q', fid.read(8))[0]
        assert record_length == 8  # 4 + 4 + 0

    def test_total_bytes(self):
        fid = io.BytesIO()
        s = 'abc'
        _write_string(fid, s)
        assert fid.tell() == 8 + 4 + 4 + len(s)


# ============== _WRITE_STRING_ARRAY ==============

class TestWriteStringArray:
    """Tests for _write_string_array()."""

    def test_writes_format_code(self):
        fid = io.BytesIO()
        _write_string_array(fid, ['a', 'bb'])
        fid.seek(8)
        code = struct.unpack('<i', fid.read(4))[0]
        assert code == format_to_code('StringArray')

    def test_writes_array_length(self):
        fid = io.BytesIO()
        _write_string_array(fid, ['a', 'bb', 'ccc'])
        fid.seek(12)
        arr_len = struct.unpack('<i', fid.read(4))[0]
        assert arr_len == 3

    def test_single_string(self):
        fid = io.BytesIO()
        _write_string_array(fid, ['only'])
        fid.seek(12)
        arr_len = struct.unpack('<i', fid.read(4))[0]
        assert arr_len == 1

    def test_empty_array(self):
        fid = io.BytesIO()
        _write_string_array(fid, [])
        fid.seek(12)
        arr_len = struct.unpack('<i', fid.read(4))[0]
        assert arr_len == 0

    def test_string_contents_written(self):
        fid = io.BytesIO()
        _write_string_array(fid, ['hi'])
        fid.seek(16)  # past: 8 (reclen) + 4 (code) + 4 (arr_len)
        str_len = struct.unpack('<i', fid.read(4))[0]
        assert str_len == 2
        assert fid.read(2) == b'hi'


# ============== _APPLY_SCALING_AND_TIME_SERIES ==============

class TestApplyScalingAndTimeSeries:
    """Tests for _apply_scaling_and_time_series()."""

    def test_no_scale_no_yts_unchanged(self):
        data = np.array([1.0, 2.0, 3.0])
        result = _apply_scaling_and_time_series(data, 'DoubleMat', 5, None, None)
        np.testing.assert_array_equal(result, [1.0, 2.0, 3.0])

    def test_scale_applied_to_1d(self):
        data = np.array([1.0, 2.0, 3.0])
        result = _apply_scaling_and_time_series(data, 'DoubleMat', 5, 2.0, None)
        np.testing.assert_allclose(result, [2.0, 4.0, 6.0])

    def test_scale_applied_to_2d_non_timeseries(self):
        data = np.array([[1.0, 2.0], [3.0, 4.0]])
        result = _apply_scaling_and_time_series(data, 'DoubleMat', 5, 3.0, None)
        np.testing.assert_allclose(result, [[3.0, 6.0], [9.0, 12.0]])

    def test_scale_on_timeseries_skips_last_row(self):
        """For timeseries (nrows==timeserieslength), scale only rows[:-1]."""
        data = np.array([[1.0, 2.0], [3.0, 4.0], [10.0, 20.0]])
        timeserieslength = 3
        result = _apply_scaling_and_time_series(data, 'DoubleMat', timeserieslength, 2.0, None)
        np.testing.assert_allclose(result[:-1, :], [[2.0, 4.0], [6.0, 8.0]])
        np.testing.assert_allclose(result[-1, :], [10.0, 20.0])  # unchanged

    def test_yts_applied_to_last_row_only(self):
        data = np.array([[1.0, 2.0], [10.0, 20.0]])
        timeserieslength = 2
        result = _apply_scaling_and_time_series(data, 'DoubleMat', timeserieslength, None, 365.25)
        np.testing.assert_allclose(result[0, :], [1.0, 2.0])  # unchanged
        np.testing.assert_allclose(result[-1, :], [10.0 * 365.25, 20.0 * 365.25])

    def test_mat_array_scale(self):
        data = [np.array([1.0, 2.0]), np.array([3.0, 4.0])]
        result = _apply_scaling_and_time_series(data, 'MatArray', 5, 2.0, None)
        np.testing.assert_allclose(result[0], [2.0, 4.0])
        np.testing.assert_allclose(result[1], [6.0, 8.0])

    def test_mat_array_no_change_with_none(self):
        data = [np.array([1.0]), np.array([2.0])]
        result = _apply_scaling_and_time_series(data, 'MatArray', 5, None, None)
        np.testing.assert_allclose(result[0], [1.0])
        np.testing.assert_allclose(result[1], [2.0])

    def test_list_data_scaled(self):
        """Lists should have each element scaled."""
        data = [1.0, 2.0, 3.0]
        result = _apply_scaling_and_time_series(data, 'DoubleMat', 5, 2.0, None)
        assert result == [2.0, 4.0, 6.0]
