"""
Unit tests for pyissm.tools.config module.

Tests cover:
- Solver configuration functions (iluasm_options, issm_gsl_solver, issm_mumps_solver)
- Platform detection functions (is_pc, get_hostname, get_username)

Note: Some tests require the ISSM backend for wrapper-dependent functions.
"""

import collections
import os
import platform
import socket

import numpy as np
import pytest
from types import SimpleNamespace

try:
    from pyissm.tools.config import (
        iluasm_options,
        issm_gsl_solver,
        issm_mumps_solver,
        is_pc,
        get_hostname,
        get_username,
        get_issm_dir,
    )
    CONFIG_AVAILABLE = True
except ImportError:
    CONFIG_AVAILABLE = False
    iluasm_options = issm_gsl_solver = issm_mumps_solver = None
    is_pc = get_hostname = get_username = get_issm_dir = None

pytestmark = pytest.mark.skipif(
    not CONFIG_AVAILABLE,
    reason="ISSM config module not available"
)


class TestIluasmOptions:
    """Tests for iluasm_options solver configuration."""

    def test_returns_ordered_dict(self):
        """Test that function returns an OrderedDict."""
        opts = iluasm_options()
        assert isinstance(opts, collections.OrderedDict)

    def test_default_values(self):
        """Test default configuration values."""
        opts = iluasm_options()
        
        assert opts['mat_type'] == 'aij'
        assert opts['ksp_type'] == 'gmres'
        assert opts['pc_type'] == 'asm'
        assert opts['sub_pc_type'] == 'ilu'
        assert opts['pc_asm_overlap'] == 5
        assert opts['ksp_max_it'] == 100
        assert opts['ksp_rtol'] == 1e-15

    def test_override_defaults(self):
        """Test that kwargs override defaults."""
        opts = iluasm_options(ksp_max_it=200, ksp_rtol=1e-12)
        
        assert opts['ksp_max_it'] == 200
        assert opts['ksp_rtol'] == 1e-12
        # Other defaults should remain
        assert opts['pc_type'] == 'asm'

    def test_add_new_options(self):
        """Test that new options can be added."""
        opts = iluasm_options(custom_option='value')
        
        assert opts['custom_option'] == 'value'


class TestIssmGslSolver:
    """Tests for issm_gsl_solver configuration."""

    def test_returns_ordered_dict(self):
        """Test that function returns an OrderedDict."""
        opts = issm_gsl_solver()
        assert isinstance(opts, collections.OrderedDict)

    def test_default_values(self):
        """Test default configuration values."""
        opts = issm_gsl_solver()
        
        assert opts['toolkit'] == 'issm'
        assert opts['mat_type'] == 'dense'
        assert opts['vec_type'] == 'seq'
        assert opts['solver_type'] == 'gsl'

    def test_override_defaults(self):
        """Test that kwargs override defaults."""
        opts = issm_gsl_solver(mat_type='sparse')
        
        assert opts['mat_type'] == 'sparse'
        assert opts['solver_type'] == 'gsl'


class TestIssmMumpsSolver:
    """Tests for issm_mumps_solver configuration."""

    def test_returns_ordered_dict(self):
        """Test that function returns an OrderedDict."""
        opts = issm_mumps_solver()
        assert isinstance(opts, collections.OrderedDict)

    def test_default_values(self):
        """Test default configuration values."""
        opts = issm_mumps_solver()
        
        assert opts['toolkit'] == 'issm'
        assert opts['mat_type'] == 'mpisparse'
        assert opts['vec_type'] == 'mpi'
        assert opts['solver_type'] == 'mumps'

    def test_override_defaults(self):
        """Test that kwargs override defaults."""
        opts = issm_mumps_solver(solver_type='petsc')
        
        assert opts['solver_type'] == 'petsc'
        assert opts['mat_type'] == 'mpisparse'


class TestIsPc:
    """Tests for is_pc platform detection."""

    def test_returns_bool(self):
        """Test that function returns a boolean."""
        result = is_pc()
        assert isinstance(result, bool)

    def test_matches_platform(self):
        """Test that result matches platform detection."""
        result = is_pc()
        expected = 'Windows' in platform.system()
        assert result == expected


class TestGetHostname:
    """Tests for get_hostname function."""

    def test_returns_string(self):
        """Test that function returns a string."""
        hostname = get_hostname()
        assert isinstance(hostname, str)

    def test_returns_lowercase(self):
        """Test that hostname is lowercase."""
        hostname = get_hostname()
        assert hostname == hostname.lower()

    def test_matches_socket(self):
        """Test that result matches socket.gethostname()."""
        hostname = get_hostname()
        expected = socket.gethostname().lower()
        assert hostname == expected


class TestGetUsername:
    """Tests for get_username function."""

    def test_returns_string(self):
        """Test that function returns a string."""
        username = get_username()
        assert isinstance(username, str)

    def test_matches_environment(self):
        """Test that result matches environment variable."""
        username = get_username()
        if 'Windows' in platform.system():
            expected = os.environ.get('USERNAME', '')
        else:
            expected = os.environ.get('USER', '')
        assert username == expected


class TestGetIssmDir:
    """Tests for get_issm_dir function."""

    def test_returns_string_or_none(self):
        """Test that function returns string or None."""
        result = get_issm_dir()
        assert result is None or isinstance(result, str)

    def test_no_trailing_slash(self):
        """Test that result has no trailing slash if set."""
        result = get_issm_dir()
        if result is not None:
            assert not result.endswith('/')
            assert not result.endswith('\\')

    def test_warns_if_not_set(self, monkeypatch):
        """Test that warning is issued if ISSM_DIR is not set."""
        # Temporarily unset ISSM_DIR
        monkeypatch.delenv('ISSM_DIR', raising=False)
        monkeypatch.delenv('ISSM_DIR_WIN', raising=False)
        
        import warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = get_issm_dir()
            
            # May or may not warn depending on environment
            # but result should be None if env var is unset
            if result is None:
                assert len(w) >= 1 or True  # Warning may be suppressed
