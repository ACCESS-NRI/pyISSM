"""
Shared pytest fixtures for pyISSM unit tests.
"""

import numpy as np
import pytest
from types import SimpleNamespace


@pytest.fixture
def simple_mesh_2d():
    """
    Create a simple 2D triangular mesh fixture.
    
    Returns a SimpleNamespace with mesh attributes for a unit square
    divided into 2 triangles.
    """
    mesh = SimpleNamespace()
    
    # 4 vertices forming a unit square
    mesh.x = np.array([0.0, 1.0, 1.0, 0.0])
    mesh.y = np.array([0.0, 0.0, 1.0, 1.0])
    
    # 2 triangles (1-based indexing)
    mesh.elements = np.array([
        [1, 2, 3],
        [1, 3, 4],
    ])
    
    mesh.numberofvertices = 4
    mesh.numberofelements = 2
    
    def dimension():
        return 2
    mesh.dimension = dimension
    
    return mesh


@pytest.fixture
def simple_mesh_3d():
    """
    Create a simple 3D mesh fixture with 2D surface layer.
    
    Returns a SimpleNamespace mimicking a 3D ISSM mesh structure.
    """
    mesh = SimpleNamespace()
    
    # 2D base mesh (4 vertices, 2 triangles)
    mesh.x2d = np.array([0.0, 1.0, 1.0, 0.0])
    mesh.y2d = np.array([0.0, 0.0, 1.0, 1.0])
    mesh.elements2d = np.array([
        [1, 2, 3],
        [1, 3, 4],
    ])
    mesh.numberofvertices2d = 4
    mesh.numberofelements2d = 2
    
    # 3D mesh with 3 layers (12 vertices total)
    mesh.numberoflayers = 3
    mesh.x = np.tile(mesh.x2d, 3)
    mesh.y = np.tile(mesh.y2d, 3)
    mesh.z = np.concatenate([np.zeros(4), np.ones(4) * 0.5, np.ones(4)])
    mesh.numberofvertices = 12
    mesh.numberofelements = 4  # 2 elements per layer gap = 4 prisms
    
    # Surface marker
    mesh.vertexonsurface = np.concatenate([np.zeros(4), np.zeros(4), np.ones(4)])
    
    def dimension():
        return 3
    mesh.dimension = dimension
    
    return mesh


@pytest.fixture
def fake_model_2d(simple_mesh_2d):
    """
    Create a fake 2D ISSM model object with minimal required attributes.
    """
    md = SimpleNamespace()
    md.mesh = simple_mesh_2d
    
    # Geometry
    md.geometry = SimpleNamespace()
    md.geometry.surface = np.array([100.0, 100.0, 100.0, 100.0])
    md.geometry.bed = np.array([0.0, 0.0, 0.0, 0.0])
    md.geometry.thickness = md.geometry.surface - md.geometry.bed
    
    # Mask
    md.mask = SimpleNamespace()
    md.mask.ice_levelset = np.array([-1.0, -1.0, -1.0, -1.0])  # All ice
    md.mask.ocean_levelset = np.array([1.0, 1.0, 1.0, 1.0])    # No ocean
    
    return md


@pytest.fixture
def fake_model_3d(simple_mesh_3d):
    """
    Create a fake 3D ISSM model object with minimal required attributes.
    """
    md = SimpleNamespace()
    md.mesh = simple_mesh_3d
    
    # Geometry (surface values)
    md.geometry = SimpleNamespace()
    md.geometry.surface = np.ones(12) * 100.0
    md.geometry.bed = np.zeros(12)
    md.geometry.thickness = md.geometry.surface - md.geometry.bed
    
    # Mask
    md.mask = SimpleNamespace()
    md.mask.ice_levelset = -1.0 * np.ones(12)
    md.mask.ocean_levelset = np.ones(12)
    
    return md


@pytest.fixture
def regular_grid():
    """
    Create a regular grid fixture for interpolation tests.
    """
    x = np.linspace(0, 1, 11)
    y = np.linspace(0, 1, 11)
    xx, yy = np.meshgrid(x, y)
    return SimpleNamespace(x=x, y=y, xx=xx, yy=yy)


@pytest.fixture
def sample_field_values():
    """
    Create sample field values for testing.
    """
    return SimpleNamespace(
        uniform=np.array([1.0, 1.0, 1.0, 1.0]),
        linear_x=np.array([0.0, 1.0, 1.0, 0.0]),
        linear_y=np.array([0.0, 0.0, 1.0, 1.0]),
        zero_crossing=np.array([-1.0, 1.0, 1.0, -1.0]),
    )
