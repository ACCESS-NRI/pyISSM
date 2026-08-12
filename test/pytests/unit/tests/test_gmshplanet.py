"""
Unit tests for pyissm.model.mesh.gmshplanet.

Property-based only (no old-vs-new comparison against the upstream MATLAB
gmshplanet): the algorithm here is intentionally different (Gmsh OCC sphere
primitive vs. a manually-constructed .geo file), so there is no meaningful
element-by-element or point-by-point "old vs. new" comparison to make. See
`goals.md` / `Implementation plan.md` for the rationale.
"""
from importlib import import_module

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
