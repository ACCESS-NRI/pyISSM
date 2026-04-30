"""
Unit tests for pyissm.tools.exp - exp_write and exp_read functions.
Uses tempfiles to avoid disk persistence.
"""

import pytest
import tempfile
import os
import numpy as np

try:
    from pyissm.tools.exp import exp_write, exp_read
    EXP_AVAILABLE = True
except ImportError:
    EXP_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not EXP_AVAILABLE,
    reason="pyissm.tools.exp not available"
)


# ============== EXP_WRITE ==============

class TestExpWrite:
    """Tests for exp_write()."""

    def test_single_contour_creates_file(self):
        contour = {'x': [0.0, 1.0, 0.0], 'y': [0.0, 0.0, 1.0], 'name': 'triangle'}
        with tempfile.NamedTemporaryFile(suffix='.exp', delete=False) as f:
            fname = f.name
        try:
            exp_write(contour, fname)
            assert os.path.exists(fname)
        finally:
            os.unlink(fname)

    def test_file_contains_name(self):
        contour = {'x': [0.0, 1.0], 'y': [0.0, 1.0], 'name': 'mycontour'}
        with tempfile.NamedTemporaryFile(suffix='.exp', delete=False, mode='w') as f:
            fname = f.name
        try:
            exp_write(contour, fname)
            with open(fname) as f:
                content = f.read()
            assert 'mycontour' in content
        finally:
            os.unlink(fname)

    def test_file_contains_coordinates(self):
        contour = {'x': [1.5, 2.5], 'y': [3.5, 4.5], 'name': 'line'}
        with tempfile.NamedTemporaryFile(suffix='.exp', delete=False, mode='w') as f:
            fname = f.name
        try:
            exp_write(contour, fname)
            with open(fname) as f:
                content = f.read()
            assert '1.5' in content
            assert '3.5' in content
        finally:
            os.unlink(fname)

    def test_contour_uses_filename_when_no_name(self):
        contour = {'x': [0.0, 1.0], 'y': [0.0, 1.0]}
        with tempfile.NamedTemporaryFile(suffix='.exp', delete=False, mode='w') as f:
            fname = f.name
        try:
            exp_write(contour, fname)
            with open(fname) as f:
                content = f.read()
            assert fname in content
        finally:
            os.unlink(fname)

    def test_list_of_contours(self):
        contours = [
            {'x': [0.0, 1.0], 'y': [0.0, 0.0], 'name': 'c1'},
            {'x': [2.0, 3.0], 'y': [0.0, 0.0], 'name': 'c2'},
        ]
        with tempfile.NamedTemporaryFile(suffix='.exp', delete=False, mode='w') as f:
            fname = f.name
        try:
            exp_write(contours, fname)
            with open(fname) as f:
                content = f.read()
            assert 'c1' in content
            assert 'c2' in content
        finally:
            os.unlink(fname)

    def test_contour_with_density(self):
        contour = {'x': [0.0, 1.0], 'y': [0.0, 1.0], 'density': 2}
        with tempfile.NamedTemporaryFile(suffix='.exp', delete=False, mode='w') as f:
            fname = f.name
        try:
            exp_write(contour, fname)
            with open(fname) as f:
                content = f.read()
            assert '2' in content
        finally:
            os.unlink(fname)

    def test_x_y_length_mismatch_raises(self):
        contour = {'x': [0.0, 1.0, 2.0], 'y': [0.0, 1.0]}  # different lengths
        with tempfile.NamedTemporaryFile(suffix='.exp', delete=False, mode='w') as f:
            fname = f.name
        try:
            with pytest.raises(RuntimeError, match='identical size'):
                exp_write(contour, fname)
        finally:
            if os.path.exists(fname):
                os.unlink(fname)

    def test_numpy_arrays_as_xy(self):
        contour = {
            'x': np.array([0.0, 1.0, 0.0]),
            'y': np.array([0.0, 0.0, 1.0]),
            'name': 'numpy_tri'
        }
        with tempfile.NamedTemporaryFile(suffix='.exp', delete=False, mode='w') as f:
            fname = f.name
        try:
            exp_write(contour, fname)
            assert os.path.exists(fname)
        finally:
            os.unlink(fname)


# ============== EXP_READ ==============

class TestExpRead:
    """Tests for exp_read()."""

    def _write_and_read(self, contour):
        """Helper: write then read back a contour."""
        with tempfile.NamedTemporaryFile(suffix='.exp', delete=False, mode='w') as f:
            fname = f.name
        try:
            exp_write(contour, fname)
            result = exp_read(fname)
        finally:
            os.unlink(fname)
        return result

    def test_returns_list(self):
        contour = {'x': [0.0, 1.0], 'y': [0.0, 1.0], 'name': 'line'}
        result = self._write_and_read(contour)
        assert isinstance(result, list)

    def test_single_contour_length(self):
        contour = {'x': [0.0, 1.0], 'y': [0.0, 1.0], 'name': 'line'}
        result = self._write_and_read(contour)
        assert len(result) == 1

    def test_name_round_trip(self):
        contour = {'x': [0.0, 1.0], 'y': [0.0, 1.0], 'name': 'myshape'}
        result = self._write_and_read(contour)
        assert result[0]['name'] == 'myshape'

    def test_xy_round_trip(self):
        contour = {'x': [1.5, 2.5, 3.5], 'y': [4.5, 5.5, 6.5], 'name': 'test'}
        result = self._write_and_read(contour)
        np.testing.assert_allclose(result[0]['x'], [1.5, 2.5, 3.5], rtol=1e-6)
        np.testing.assert_allclose(result[0]['y'], [4.5, 5.5, 6.5], rtol=1e-6)

    def test_nods_correct(self):
        contour = {'x': [0.0, 1.0, 2.0], 'y': [0.0, 1.0, 2.0], 'name': 'test'}
        result = self._write_and_read(contour)
        assert result[0]['nods'] == 3

    def test_density_round_trip(self):
        contour = {'x': [0.0, 1.0], 'y': [0.0, 1.0], 'name': 'test', 'density': 3}
        result = self._write_and_read(contour)
        assert result[0]['density'] == 3.0

    def test_multiple_contours_round_trip(self):
        contours = [
            {'x': [0.0, 1.0], 'y': [0.0, 0.0], 'name': 'first'},
            {'x': [2.0, 3.0], 'y': [0.0, 0.0], 'name': 'second'},
        ]
        with tempfile.NamedTemporaryFile(suffix='.exp', delete=False, mode='w') as f:
            fname = f.name
        try:
            exp_write(contours, fname)
            result = exp_read(fname)
        finally:
            os.unlink(fname)
        assert len(result) == 2
        assert result[0]['name'] == 'first'
        assert result[1]['name'] == 'second'

    def test_closed_contour_detection(self):
        """A contour whose first and last points match should be detected as closed."""
        contour = {
            'x': [0.0, 1.0, 0.0],  # last == first
            'y': [0.0, 0.0, 0.0],
            'name': 'closed'
        }
        result = self._write_and_read(contour)
        assert result[0]['closed'] is True

    def test_open_contour_detection(self):
        """A contour whose endpoints differ should be detected as open."""
        contour = {
            'x': [0.0, 1.0, 2.0],  # last != first
            'y': [0.0, 0.0, 0.0],
            'name': 'open'
        }
        result = self._write_and_read(contour)
        assert result[0]['closed'] is False

    def test_file_not_found_raises(self):
        with pytest.raises(IOError, match='does not exist'):
            exp_read('/nonexistent/path/to/file.exp')

    def test_x_y_arrays_are_numpy(self):
        contour = {'x': [0.0, 1.0], 'y': [0.0, 1.0], 'name': 'test'}
        result = self._write_and_read(contour)
        assert isinstance(result[0]['x'], np.ndarray)
        assert isinstance(result[0]['y'], np.ndarray)
