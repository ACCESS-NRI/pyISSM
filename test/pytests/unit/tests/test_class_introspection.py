"""
Tests for class introspection helpers on model classes:

- #76: get_fieldnames() on every model class (inherited from manage_state).
- #94: available() lists the parameterisations of a model component, from the module,
  from a class, or from a model field.
"""

import pickle

import pytest

from pyissm.model.classes import class_registry, friction, results, smb


def _registered():
    return sorted(class_registry.CLASS_REGISTRY.items())


class TestGetFieldnames:

    def test_friction_fields_in_definition_order(self):
        f = friction.default()
        assert f.get_fieldnames() == list(vars(f))
        assert f.get_fieldnames()[0] == 'coefficient'

    def test_returns_list_of_public_names(self):
        f = friction.default()
        f._scratch = 1
        names = f.get_fieldnames()
        assert isinstance(names, list)
        assert '_scratch' not in names

    @pytest.mark.parametrize('key, cls', _registered())
    def test_every_registered_class_has_it(self, key, cls):
        assert callable(getattr(cls, 'get_fieldnames', None)), key

    def test_solutionstep(self):
        step = results.solutionstep()
        step.Vel = 1.
        step.time = 0.
        assert step.get_fieldnames() == ['Vel', 'time']

    def test_solution_lists_result_fields_across_steps(self):
        sol = results.solution([results.solutionstep(), results.solutionstep()])
        sol[0].Vel, sol[0].time = 1., 0.
        sol[1].Vel, sol[1].Thickness, sol[1].time = 2., 3., 1.
        assert sol.get_fieldnames() == ['Vel', 'time', 'Thickness']
        assert 'Thickness' in repr(sol)

    def test_empty_solution(self):
        assert results.solution().get_fieldnames() == []

    def test_not_part_of_saved_state(self):
        f = pickle.loads(pickle.dumps(friction.default()))
        assert 'get_fieldnames' not in vars(f) and 'available' not in vars(f)


class TestAvailable:

    def test_module_function(self):
        options = friction.available()
        assert isinstance(options, dict)
        assert list(options)[:3] == ['default', 'coulomb', 'coulomb2']
        assert {'schoof', 'weertman', 'regcoulomb'} <= set(options)
        assert options['default'] == 'Default friction class for ISSM.'

    def test_class_and_instance_methods_match_module(self):
        assert friction.schoof.available() == friction.available()
        assert friction.schoof().available() == friction.available()
        assert smb.default().available() == smb.available()

    def test_model_field(self):
        from pyissm.model import Model
        md = Model()
        assert md.friction.available() == friction.available()
        assert md.smb.available() == smb.available()

    def test_single_class_component(self):
        from pyissm.model.classes import geometry
        assert list(geometry.available()) == ['geometry']

    def test_every_registered_class_is_listed(self):
        for key, cls in _registered():
            category, name = key.rsplit('.', 1)
            assert class_registry.available(category)[name], key

    def test_descriptions_are_single_line(self):
        for key, cls in _registered():
            category, name = key.rsplit('.', 1)
            assert '\n' not in class_registry.available(category)[name]

    def test_unknown_component(self):
        with pytest.raises(ValueError, match = 'no classes registered'):
            class_registry.available('not_a_component')
