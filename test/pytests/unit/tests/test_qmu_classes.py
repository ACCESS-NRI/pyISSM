"""
Unit tests for pyissm.model.classes.qmu - QMU Dakota model classes.
"""

import pytest
import numpy as np

try:
    from pyissm.model.classes import qmu
    QMU_AVAILABLE = True
except ImportError:
    QMU_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not QMU_AVAILABLE,
    reason="pyissm.model.classes.qmu not available"
)


# ============== QMU.DEFAULT TESTS ==============

class TestQmuDefault:
    """Tests for qmu.default class."""

    def test_init(self):
        q = qmu.default()
        assert q is not None

    def test_has_isdakota(self):
        q = qmu.default()
        assert hasattr(q, 'isdakota')
        assert q.isdakota == 0

    def test_has_output(self):
        q = qmu.default()
        assert hasattr(q, 'output')
        assert q.output == 0

    def test_has_variables(self):
        q = qmu.default()
        assert hasattr(q, 'variables')

    def test_has_correlation_matrix(self):
        q = qmu.default()
        assert hasattr(q, 'correlation_matrix')

    def test_has_responses(self):
        q = qmu.default()
        assert hasattr(q, 'responses')

    def test_has_method(self):
        q = qmu.default()
        assert hasattr(q, 'method')
        assert isinstance(q.method, dict)

    def test_has_params(self):
        q = qmu.default()
        assert hasattr(q, 'params')

    def test_has_statistics(self):
        q = qmu.default()
        assert hasattr(q, 'statistics')

    def test_has_results(self):
        q = qmu.default()
        assert hasattr(q, 'results')

    def test_has_numberofresponses(self):
        q = qmu.default()
        assert hasattr(q, 'numberofresponses')
        assert q.numberofresponses == 0

    def test_has_variabledescriptors(self):
        q = qmu.default()
        assert hasattr(q, 'variabledescriptors')

    def test_has_variablepartitions(self):
        q = qmu.default()
        assert hasattr(q, 'variablepartitions')

    def test_has_variablepartitions_npart(self):
        q = qmu.default()
        assert hasattr(q, 'variablepartitions_npart')

    def test_has_variablepartitions_nt(self):
        q = qmu.default()
        assert hasattr(q, 'variablepartitions_nt')

    def test_has_responsedescriptors(self):
        q = qmu.default()
        assert hasattr(q, 'responsedescriptors')

    def test_has_responsepartitions(self):
        q = qmu.default()
        assert hasattr(q, 'responsepartitions')

    def test_has_responsepartitions_npart(self):
        q = qmu.default()
        assert hasattr(q, 'responsepartitions_npart')

    def test_has_responsepartitions_nt(self):
        q = qmu.default()
        assert hasattr(q, 'responsepartitions_nt')

    def test_has_mass_flux_profile_directory(self):
        q = qmu.default()
        assert hasattr(q, 'mass_flux_profile_directory')

    def test_has_mass_flux_profiles(self):
        q = qmu.default()
        assert hasattr(q, 'mass_flux_profiles')

    def test_has_mass_flux_segments(self):
        q = qmu.default()
        assert hasattr(q, 'mass_flux_segments')
        assert isinstance(q.mass_flux_segments, list)

    def test_has_adjacency(self):
        q = qmu.default()
        assert hasattr(q, 'adjacency')

    def test_has_vertex_weight(self):
        q = qmu.default()
        assert hasattr(q, 'vertex_weight')

    def test_repr(self):
        q = qmu.default()
        r = repr(q)
        assert isinstance(r, str)

    def test_str(self):
        q = qmu.default()
        assert isinstance(str(q), str)

    def test_init_from_other(self):
        q1 = qmu.default()
        q1.numberofresponses = 5
        q2 = qmu.default(other=q1)
        assert q2.numberofresponses == 5


# ============== QMU.STATISTICS TESTS ==============

class TestQmuStatistics:
    """Tests for qmu.statistics class."""

    def test_init(self):
        s = qmu.statistics()
        assert s is not None

    def test_has_nfiles_per_directory(self):
        s = qmu.statistics()
        assert hasattr(s, 'nfiles_per_directory')
        assert s.nfiles_per_directory == 5

    def test_has_ndirectories(self):
        s = qmu.statistics()
        assert hasattr(s, 'ndirectories')
        assert s.ndirectories == 50

    def test_has_method(self):
        s = qmu.statistics()
        assert hasattr(s, 'method')
        assert isinstance(s.method, list)

    def test_method_has_name(self):
        s = qmu.statistics()
        assert 'name' in s.method[0]
        assert s.method[0]['name'] == 'None'

    def test_method_has_fields(self):
        s = qmu.statistics()
        assert 'fields' in s.method[0]
        assert isinstance(s.method[0]['fields'], list)

    def test_method_has_steps(self):
        s = qmu.statistics()
        assert 'steps' in s.method[0]
        assert isinstance(s.method[0]['steps'], list)

    def test_method_has_nbins(self):
        s = qmu.statistics()
        assert 'nbins' in s.method[0]

    def test_method_has_indices(self):
        s = qmu.statistics()
        assert 'indices' in s.method[0]
        assert isinstance(s.method[0]['indices'], list)

    def test_repr(self):
        s = qmu.statistics()
        r = repr(s)
        assert isinstance(r, str)

    def test_str(self):
        s = qmu.statistics()
        assert isinstance(str(s), str)

    def test_init_from_other(self):
        s1 = qmu.statistics()
        s1.nfiles_per_directory = 10
        s2 = qmu.statistics(other=s1)
        assert s2.nfiles_per_directory == 10
