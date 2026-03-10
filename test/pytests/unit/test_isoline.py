from types import SimpleNamespace

import numpy as np
import pytest

from pyissm.tools.exp import isoline


def make_fake_md():
    """
    Build a tiny fake model object with only the attributes isoline needs.
    """
    md = SimpleNamespace()
    md.mesh = SimpleNamespace()

    # Square split into 2 triangles
    md.mesh.x = np.array([0.0, 1.0, 1.0, 0.0])
    md.mesh.y = np.array([0.0, 0.0, 1.0, 1.0])

    # 1-based indexing if MATLAB-style conventions are used internally
    md.mesh.elements = np.array([
        [1, 2, 3],
        [1, 3, 4],
    ])

    return md


def test_isoline_returns_struct_output():
    md = make_fake_md()

    # Field crosses zero through the domain
    field = np.array([-1.0, 1.0, 1.0, -1.0])

    contours, edges_tria = isoline(md, field, value=0.0, output="struct")

    assert contours is not None
    assert isinstance(contours, list)
    assert len(contours) >= 1

    assert edges_tria is not None
    assert edges_tria.shape[1] == 3


def test_isoline_matrix_output_has_two_columns():
    md = make_fake_md()
    field = np.array([-1.0, 1.0, 1.0, -1.0])

    contour_matrix, edges_tria = isoline(md, field, value=0.0, output="matrix")

    assert contour_matrix is not None
    assert contour_matrix.ndim == 2
    assert contour_matrix.shape[1] == 2
    assert edges_tria is not None


def test_isoline_longest_output_not_empty():
    md = make_fake_md()
    field = np.array([-1.0, 1.0, 1.0, -1.0])

    longest, _ = isoline(md, field, value=0.0, output="longest")

    assert longest is not None
    assert isinstance(longest, np.ndarray)
    assert longest.ndim == 2
    assert longest.shape[1] == 2


def test_isoline_rejects_bad_field_length():
    md = make_fake_md()
    bad_field = np.array([1.0, 2.0])

    with pytest.raises((ValueError, IndexError, AssertionError)):
        isoline(md, bad_field, value=0.0, output="struct")