"""
Unit tests for pyissm.model.mesh.gmshplanet.

Property-based only (no old-vs-new comparison against the upstream MATLAB
gmshplanet): the algorithm here is intentionally different (Gmsh OCC sphere
primitive vs. a manually-constructed .geo file), so there is no meaningful
element-by-element or point-by-point "old vs. new" comparison to make. See
`goals.md` / `Implementation plan.md` for the rationale.
"""
from importlib import import_module
from unittest.mock import patch

import numpy as np
import pytest

try:
    import gmsh
    GMSH_AVAILABLE = True
except ImportError:
    GMSH_AVAILABLE = False

pytestmark = pytest.mark.skipif(not GMSH_AVAILABLE, reason="gmsh not available")

if GMSH_AVAILABLE:
    mesh_module = import_module('pyissm.model.mesh')
    from pyissm.model.Model import Model

RADIUS_KM = 6371.0
RESOLUTION_KM = 2000.0


def _assert_valid_planet_mesh(md, radius_km):
    """Shared property checks for a populated mesh3dsurface from gmshplanet."""
    m = md.mesh
    radius_m = radius_km * 1.e3

    assert m.numberofvertices > 0
    assert m.numberofelements > 0

    assert m.elements.shape == (m.numberofelements, 3)
    assert m.elements.min() >= 1
    assert m.elements.max() <= m.numberofvertices
    # No orphan nodes: every vertex is referenced by at least one element.
    nodes = np.arange(1, m.numberofvertices + 1)
    assert np.all(np.isin(nodes, m.elements))

    for field in ('x', 'y', 'z', 'lat', 'long', 'r'):
        array = getattr(m, field)
        assert array.shape == (m.numberofvertices,)
        assert np.all(np.isfinite(array))

    r = np.sqrt(m.x**2 + m.y**2 + m.z**2)
    assert np.allclose(r, radius_m, atol=1.)
    assert np.allclose(m.r, radius_m, atol=1.)
    assert np.all(np.abs(m.lat) < 90.)


def _mean_edge_length_km(md):
    """Mean triangle edge length [km] of a populated gmshplanet mesh."""
    m = md.mesh
    c = np.column_stack([m.x, m.y, m.z])
    els = m.elements - 1
    p0, p1, p2 = c[els[:, 0]], c[els[:, 1]], c[els[:, 2]]
    return np.concatenate([
        np.linalg.norm(p0 - p1, axis=1),
        np.linalg.norm(p1 - p2, axis=1),
        np.linalg.norm(p2 - p0, axis=1),
    ]).mean() / 1.e3


class TestGmshplanetBasic:
    """Tests for the initial (uniform-resolution) mesh path."""

    def test_produces_valid_populated_mesh(self):
        md = Model()
        md = mesh_module.gmshplanet(md, radius=RADIUS_KM, resolution=RESOLUTION_KM)
        _assert_valid_planet_mesh(md, RADIUS_KM)

    def test_consistency_check_passes(self):
        md = Model()
        md = mesh_module.gmshplanet(md, radius=RADIUS_KM, resolution=RESOLUTION_KM)
        md.mesh.check_consistency(md, solution='', analyses=['MasstransportAnalysis'])
        assert md.private.isconsistent

    def test_raises_if_mesh_not_empty(self):
        md = Model()
        md = mesh_module.gmshplanet(md, radius=RADIUS_KM, resolution=RESOLUTION_KM)
        with pytest.raises(RuntimeError):
            mesh_module.gmshplanet(md, radius=RADIUS_KM, resolution=RESOLUTION_KM)

    def test_finer_resolution_yields_more_vertices(self):
        md_coarse = mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=3000.0)
        md_fine = mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=1000.0)
        assert md_fine.mesh.numberofvertices > md_coarse.mesh.numberofvertices

    def test_preexisting_session_survives_and_is_left_undisturbed(self):
        """gmshplanet() must not finalize a Gmsh session it did not create,
        must remove only its own temporary model(s), and must restore every
        session-global option (algorithm, mesh-size bounds/factor) it changed."""
        gmsh.initialize()
        try:
            gmsh.model.add('caller_model')
            gmsh.option.setNumber('Mesh.Algorithm', 8)
            gmsh.option.setNumber('Mesh.MeshSizeMin', 42.)
            gmsh.option.setNumber('Mesh.MeshSizeMax', 4242.)
            gmsh.option.setNumber('Mesh.MeshSizeFactor', 0.1)
            before_model_names = set(gmsh.model.list())

            md = mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM)

            assert gmsh.isInitialized()
            assert gmsh.model.getCurrent() == 'caller_model'
            assert set(gmsh.model.list()) == before_model_names
            assert gmsh.option.getNumber('Mesh.Algorithm') == 8
            assert gmsh.option.getNumber('Mesh.MeshSizeMin') == 42.
            assert gmsh.option.getNumber('Mesh.MeshSizeMax') == 4242.
            assert gmsh.option.getNumber('Mesh.MeshSizeFactor') == 0.1
        finally:
            if gmsh.isInitialized():
                gmsh.finalize()

        _assert_valid_planet_mesh(md, RADIUS_KM)

    def test_preexisting_meshsizefactor_does_not_scale_output_and_is_restored(self):
        """A non-default Mesh.MeshSizeFactor inherited from a caller-owned
        session must not silently scale the requested resolution (Gmsh
        applies this factor after Mesh.MeshSizeMin/Max), and must be
        restored to the caller's value afterward."""
        gmsh.initialize()
        try:
            gmsh.model.add('caller_model')
            gmsh.option.setNumber('Mesh.MeshSizeFactor', 0.1)

            md = mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM)

            assert gmsh.option.getNumber('Mesh.MeshSizeFactor') == 0.1
        finally:
            if gmsh.isInitialized():
                gmsh.finalize()

        _assert_valid_planet_mesh(md, RADIUS_KM)
        # A factor of 0.1 left in effect would produce elements roughly an
        # order of magnitude smaller than RESOLUTION_KM; confirm the mesh
        # is sized to the requested resolution instead.
        mean_edge_km = _mean_edge_length_km(md)
        assert mean_edge_km > 0.5 * RESOLUTION_KM, (
            f'mean edge length {mean_edge_km:.1f} km suggests a pre-existing '
            f'Mesh.MeshSizeFactor scaled down the requested resolution')

    def test_preexisting_models_named_planet_do_not_collide(self):
        """A caller-owned session may already have a model literally named
        'planet' (the name gmshplanet() itself used to hardcode); it must
        not be touched, mistaken for gmshplanet's own temporary model, or
        cause a name collision."""
        gmsh.initialize()
        try:
            gmsh.model.add('planet')
            gmsh.model.occ.addSphere(0, 0, 0, 1.)
            gmsh.model.occ.synchronize()
            caller_planet_entities = gmsh.model.getEntities()

            md = mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM)

            assert gmsh.isInitialized()
            assert gmsh.model.getCurrent() == 'planet'
            assert gmsh.model.getEntities() == caller_planet_entities
        finally:
            if gmsh.isInitialized():
                gmsh.finalize()

        _assert_valid_planet_mesh(md, RADIUS_KM)

    def test_failure_during_generation_cleans_up_preexisting_session(self):
        """If mesh generation fails partway through, gmshplanet() must still
        leave a caller-owned session's models and options undisturbed, and
        the original exception must still propagate."""
        gmsh.initialize()
        try:
            gmsh.model.add('caller_model')
            gmsh.option.setNumber('Mesh.Algorithm', 8)
            before_model_names = set(gmsh.model.list())

            with patch.object(gmsh.model.mesh, 'generate', side_effect=RuntimeError('forced failure')):
                with pytest.raises(RuntimeError, match='forced failure'):
                    mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM)

            assert gmsh.isInitialized()
            assert gmsh.model.getCurrent() == 'caller_model'
            assert set(gmsh.model.list()) == before_model_names
            assert gmsh.option.getNumber('Mesh.Algorithm') == 8
        finally:
            if gmsh.isInitialized():
                gmsh.finalize()


class TestGmshplanetRefine:
    """Tests for the refine=/refinemetric= adaptive-remesh path."""

    def _build_prior(self):
        md = mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM)
        return md.mesh

    def test_refine_requires_both_args(self):
        prior = self._build_prior()
        with pytest.raises(ValueError):
            mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM, refine=prior)
        with pytest.raises(ValueError):
            mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM,
                                    refinemetric=np.ones(prior.numberofvertices))

    def test_refinemetric_length_mismatch_raises(self):
        prior = self._build_prior()
        with pytest.raises(ValueError):
            mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM,
                                    refine=prior, refinemetric=np.ones(prior.numberofvertices - 1))

    def test_refinemetric_scalar_raises_valueerror_not_typeerror(self):
        prior = self._build_prior()
        with pytest.raises(ValueError):
            mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM,
                                    refine=prior, refinemetric=500.e3)

    def test_refinemetric_2d_array_raises(self):
        prior = self._build_prior()
        metric = np.full((prior.numberofvertices, 2), 500.e3)
        with pytest.raises(ValueError):
            mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM,
                                    refine=prior, refinemetric=metric)

    def test_refinemetric_non_finite_raises(self):
        prior = self._build_prior()
        metric = np.full(prior.numberofvertices, 500.e3)
        metric[0] = np.nan
        with pytest.raises(ValueError):
            mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM,
                                    refine=prior, refinemetric=metric)

    def test_refinemetric_non_positive_raises(self):
        prior = self._build_prior()
        metric = np.full(prior.numberofvertices, 500.e3)
        metric[0] = 0.
        with pytest.raises(ValueError):
            mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM,
                                    refine=prior, refinemetric=metric)

    def test_refinemetric_above_1000km_is_not_clamped(self):
        """A large uniform metric should produce elements sized to match
        it, not silently capped at ~1,000 km by a leftover Mesh.MeshSizeMax."""
        prior = self._build_prior()
        metric = np.full(prior.numberofvertices, 3000.e3)
        md = mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM,
                                     refine=prior, refinemetric=metric)
        _assert_valid_planet_mesh(md, RADIUS_KM)

        # A clamp at 1,000 km would keep this well under 1,000 km; an
        # honored 3,000 km metric should land close to that instead.
        mean_edge_km = _mean_edge_length_km(md)
        assert mean_edge_km > 1500., (
            f'mean edge length {mean_edge_km:.1f} km suggests the >1000 km '
            f'metric was clamped rather than honored')

    def test_refine_ignores_preexisting_narrow_meshsizemax(self):
        """A narrow Mesh.MeshSizeMax already set on a Gmsh session that was
        live *before* gmshplanet() was called must not silently clamp the
        refine metric. gmsh.initialize() is a no-op (not a reset) when a
        session is already active, so a stale option can otherwise survive
        into gmshplanet()'s call. gmshplanet() must also leave that
        caller-owned session running and undisturbed afterward - it doesn't
        own it, so it must not finalize it, leak its own temporary
        model(s)/view into it, or leave any option it changed unrestored."""
        prior = self._build_prior()
        metric = np.full(prior.numberofvertices, 3000.e3)

        gmsh.initialize()
        try:
            gmsh.model.add('caller_model')
            gmsh.option.setNumber('Mesh.MeshSizeMax', 200.e3)
            gmsh.option.setNumber('Mesh.MeshSizeFactor', 0.1)
            before_model_names = set(gmsh.model.list())

            md = mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM,
                                         refine=prior, refinemetric=metric)

            # The caller's session, current model, and own models must all
            # have survived, with none of gmshplanet's temporary models left
            # behind and the caller's own option values restored.
            assert gmsh.isInitialized()
            assert gmsh.model.getCurrent() == 'caller_model'
            assert set(gmsh.model.list()) == before_model_names
            assert gmsh.option.getNumber('Mesh.MeshSizeMax') == 200.e3
            assert gmsh.option.getNumber('Mesh.MeshSizeFactor') == 0.1
        finally:
            if gmsh.isInitialized():
                gmsh.finalize()

        _assert_valid_planet_mesh(md, RADIUS_KM)
        mean_edge_km = _mean_edge_length_km(md)
        assert mean_edge_km > 1500., (
            f'mean edge length {mean_edge_km:.1f} km suggests a pre-existing '
            f'Mesh.MeshSizeMax clamped the metric')

    def test_refine_preexisting_models_named_refine_source_do_not_collide(self):
        """A caller-owned session may already have models literally named
        'planet' and/or 'refine_source' (the names this refine path itself
        used to hardcode); they must not be touched or cause a collision."""
        prior = self._build_prior()
        metric = np.full(prior.numberofvertices, 500.e3)

        gmsh.initialize()
        try:
            gmsh.model.add('refine_source')
            gmsh.model.occ.addSphere(0, 0, 0, 1.)
            gmsh.model.occ.synchronize()
            caller_refine_source_entities = gmsh.model.getEntities()
            gmsh.model.add('planet')
            before_model_names = set(gmsh.model.list())

            md = mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM,
                                         refine=prior, refinemetric=metric)

            assert gmsh.isInitialized()
            assert gmsh.model.getCurrent() == 'planet'
            assert set(gmsh.model.list()) == before_model_names
            gmsh.model.setCurrent('refine_source')
            assert gmsh.model.getEntities() == caller_refine_source_entities
        finally:
            if gmsh.isInitialized():
                gmsh.finalize()

        _assert_valid_planet_mesh(md, RADIUS_KM)

    def test_refine_produces_valid_mesh(self):
        prior = self._build_prior()
        metric = np.full(prior.numberofvertices, 500.e3)
        md = mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM,
                                     refine=prior, refinemetric=metric)
        _assert_valid_planet_mesh(md, RADIUS_KM)

    def test_refine_produces_spatially_varying_resolution(self):
        """Elements near a small-metric target region come out meaningfully
        smaller than elements near a large-metric region, and the mesh
        stays valid throughout."""
        prior = self._build_prior()
        target = np.array([RADIUS_KM * 1.e3, 0., 0.])
        p = np.column_stack([prior.x, prior.y, prior.z])
        dist = np.linalg.norm(p - target, axis=1)
        metric = np.clip(dist * 0.15, 50.e3, 2000.e3)

        md = mesh_module.gmshplanet(Model(), radius=RADIUS_KM, resolution=RESOLUTION_KM,
                                     refine=prior, refinemetric=metric)
        _assert_valid_planet_mesh(md, RADIUS_KM)

        m = md.mesh
        c = np.column_stack([m.x, m.y, m.z])
        els = m.elements - 1
        centroids = (c[els[:, 0]] + c[els[:, 1]] + c[els[:, 2]]) / 3.
        dist_c = np.linalg.norm(centroids - target, axis=1)
        near = dist_c < 1000.e3
        far = dist_c > 4000.e3
        assert near.sum() > 0 and far.sum() > 0

        def _mean_edge(mask):
            p0, p1, p2 = c[els[mask, 0]], c[els[mask, 1]], c[els[mask, 2]]
            return np.concatenate([
                np.linalg.norm(p0 - p1, axis=1),
                np.linalg.norm(p1 - p2, axis=1),
                np.linalg.norm(p2 - p0, axis=1),
            ]).mean()

        assert _mean_edge(near) < _mean_edge(far)
