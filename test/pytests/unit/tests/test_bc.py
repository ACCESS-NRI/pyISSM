"""Tests for boundary-condition helpers."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyissm.model import bc


def _model(mesh, ice_levelset):
    return SimpleNamespace(
        mesh=mesh,
        mask=SimpleNamespace(ice_levelset=np.asarray(ice_levelset)),
        stressbalance=SimpleNamespace(),
        inversion=SimpleNamespace(vx_obs=np.nan, vy_obs=np.nan),
    )


def test_set_sb_dirichlet_bc_uses_2d_segments_for_triangles():
    mesh = SimpleNamespace(
        numberofvertices=5,
        segments=np.array([
            [1, 2, 1],
            [2, 3, 1],
            [3, 4, 2],
            [4, 1, 2],
        ]),
        element_type=lambda: 'Tria',
    )
    md = _model(mesh, -np.ones(mesh.numberofvertices))

    with pytest.warns(UserWarning, match='No observed velocities found'):
        bc._set_sb_dirichlet_bc(md)

    np.testing.assert_array_equal(md.stressbalance.spcvx[:4], 0)
    assert np.isnan(md.stressbalance.spcvx[4])


def test_set_sb_dirichlet_bc_expands_3d_segments_through_every_layer():
    numberofvertices2d = 6
    numberoflayers = 3
    mesh = SimpleNamespace(
        numberofvertices=numberofvertices2d * numberoflayers,
        numberofvertices2d=numberofvertices2d,
        numberoflayers=numberoflayers,
        segments2d=np.array([
            [1, 2, 1],
            [2, 3, 1],
            [3, 4, 2],
            [4, 5, 3],
            [5, 6, 4],
            [6, 1, 4],
        ]),
        element_type=lambda: 'Penta',
    )
    levelset2d = np.array([0, 0, 0, -1, -1, -1])
    md = _model(mesh, np.tile(levelset2d, numberoflayers))

    with pytest.warns(UserWarning, match='No observed velocities found'):
        bc._set_sb_dirichlet_bc(md)

    # Node 2 belongs only to ice-front edges and remains unconstrained on all
    # layers. Every other boundary node also belongs to a non-front face.
    unconstrained = np.array([1, 7, 13])
    constrained = np.setdiff1d(np.arange(mesh.numberofvertices), unconstrained)
    assert np.all(np.isnan(md.stressbalance.spcvx[unconstrained]))
    np.testing.assert_array_equal(md.stressbalance.spcvx[constrained], 0)
    np.testing.assert_array_equal(md.stressbalance.spcvy[constrained], 0)
    np.testing.assert_array_equal(md.stressbalance.spcvz[constrained], 0)
