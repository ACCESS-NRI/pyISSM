"""
Integration test for the discrete-model + View transfer mechanism in
scripts/gmsh_refine_transfer.py: proves it produces spatially varying mesh
sizes from a prior mesh's per-vertex metric, before that mechanism is
folded into gmshplanet's refine= mode.

This tests the resolved Gmsh API mechanism directly, not pyissm code -
gmshplanet() itself does not exist yet. Once it is implemented, this can be
superseded by (or merged into) test_gmshplanet.py's own refine-mode
coverage.
"""
import importlib.util
import pathlib

import numpy as np
import pytest

try:
    import gmsh
    GMSH_AVAILABLE = True
except ImportError:
    GMSH_AVAILABLE = False

pytestmark = pytest.mark.skipif(not GMSH_AVAILABLE, reason="gmsh not available")

if GMSH_AVAILABLE:
    _script_path = (pathlib.Path(__file__).resolve().parents[4]
                     / 'scripts' / 'gmsh_refine_transfer.py')
    _spec = importlib.util.spec_from_file_location('gmsh_refine_transfer', _script_path)
    _mod = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_mod)
    attach_refine_background_field = _mod.attach_refine_background_field


RADIUS = 1000.0


def _build_prior_mesh():
    gmsh.model.add('prior')
    gmsh.model.occ.addSphere(0, 0, 0, RADIUS)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber('Mesh.MeshSizeMin', 300)
    gmsh.option.setNumber('Mesh.MeshSizeMax', 300)
    gmsh.model.mesh.generate(2)

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    coords = node_coords.reshape(-1, 3)
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(dim=2)
    tri_idx = list(elem_types).index(2)

    dist_from_pole = np.linalg.norm(coords - np.array([RADIUS, 0, 0]), axis=1)
    metric = np.clip(dist_from_pole * 0.15, 20.0, 300.0)
    return node_tags, node_coords, elem_tags[tri_idx], elem_node_tags[tri_idx], metric


def test_refine_transfer_produces_spatially_varying_mesh_and_stays_valid():
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)
    try:
        (prior_node_tags, prior_node_coords,
         prior_tri_tags, prior_tri_nodes, metric) = _build_prior_mesh()

        gmsh.model.add('planet')
        gmsh.model.setCurrent('planet')
        gmsh.model.occ.addSphere(0, 0, 0, RADIUS)
        gmsh.model.occ.synchronize()

        field_tag = attach_refine_background_field(
            gmsh, prior_node_tags, prior_node_coords,
            prior_tri_tags, prior_tri_nodes, metric)
        gmsh.model.mesh.field.setAsBackgroundMesh(field_tag)
        gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
        gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
        gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)
        gmsh.option.setNumber('Mesh.MeshSizeMin', 1)
        gmsh.option.setNumber('Mesh.MeshSizeMax', 1e6)
        gmsh.model.mesh.generate(2)

        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        coords = node_coords.reshape(-1, 3)
        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(dim=2)
        tri_idx = list(elem_types).index(2)
        tris = np.asarray(elem_node_tags[tri_idx]).reshape(-1, 3)

        # Validity: nonempty, well-formed mesh; every referenced node tag
        # exists (no dangling connectivity).
        assert len(node_tags) > 0
        assert len(tris) > 0
        assert set(np.unique(tris)).issubset(set(node_tags))

        tag_to_idx = {t: i for i, t in enumerate(node_tags)}
        idx_tris = np.vectorize(tag_to_idx.get)(tris)
        p0 = coords[idx_tris[:, 0]]
        p1 = coords[idx_tris[:, 1]]
        p2 = coords[idx_tris[:, 2]]
        centroids = (p0 + p1 + p2) / 3.0
        edge_len = np.concatenate([
            np.linalg.norm(p0 - p1, axis=1),
            np.linalg.norm(p1 - p2, axis=1),
            np.linalg.norm(p2 - p0, axis=1),
        ])
        assert np.all(np.isfinite(edge_len))
        assert np.all(edge_len > 0)

        # Refinement effect: elements near the pole (small metric) must be
        # meaningfully smaller than elements far from it (large metric).
        dist_c = np.linalg.norm(centroids - np.array([RADIUS, 0, 0]), axis=1)
        near_pole = dist_c < 200
        far_pole = dist_c > 800
        assert near_pole.sum() > 0 and far_pole.sum() > 0

        def _mean_edge(mask):
            return np.concatenate([
                np.linalg.norm(p0 - p1, axis=1)[mask],
                np.linalg.norm(p1 - p2, axis=1)[mask],
                np.linalg.norm(p2 - p0, axis=1)[mask],
            ]).mean()

        mean_near = _mean_edge(near_pole)
        mean_far = _mean_edge(far_pole)
        assert mean_near < 0.5 * mean_far, (
            f'refinement had no meaningful effect: mean edge near pole '
            f'({mean_near:.1f}) not much smaller than far from pole '
            f'({mean_far:.1f})')
    finally:
        gmsh.finalize()
