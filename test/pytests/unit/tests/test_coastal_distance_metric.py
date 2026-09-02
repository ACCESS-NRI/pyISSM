"""
Unit tests for pyissm.model.mesh.coastal_distance_metric.

Property-based only: this function has no upstream equivalent to compare
against (the MATLAB SlrGRACE tutorial's mesh refinement uses a hand-rolled
O(N^2) nested loop, not a reusable metric function), so there is no
meaningful old-vs-new comparison to make. See `goals.md` / `Implementation
plan.md` for the never-port rationale. Needs no optional dependency
(gmsh/cartopy) or a native Model() - it operates on plain lat/long/ocean
arrays.
"""
from importlib import import_module

import numpy as np
import pytest

mesh_module = import_module('pyissm.model.mesh')

RADIUS_M = 6371.012e3


def _lon_split_ocean_land(n_per_side=50, band_deg=5.0):
    """
    Synthetic points along the equator straddling a long=0 boundary: ocean
    (long < 0) on one side, land (long >= 0) on the other, at a range of
    distances from the boundary.
    """
    offsets = np.linspace(0.01, band_deg, n_per_side)
    long = np.concatenate([-offsets[::-1], offsets])
    lat = np.zeros_like(long)
    ocean = (long < 0).astype(int)
    return lat, long, ocean


class TestCoastalDistanceMetric:

    def test_shape_matches_input(self):
        lat, long, ocean = _lon_split_ocean_land()
        metric = mesh_module.coastal_distance_metric(lat, long, ocean,
                                                      mindist_coast=1e3, mindist_land=2e3, maxdist=1e6)
        assert metric.shape == (len(lat), )

    def test_all_finite_and_positive(self):
        lat, long, ocean = _lon_split_ocean_land()
        metric = mesh_module.coastal_distance_metric(lat, long, ocean,
                                                      mindist_coast=1e3, mindist_land=2e3, maxdist=1e6)
        assert np.all(np.isfinite(metric))
        assert np.all(metric > 0)

    def test_floor_is_respected_per_class(self):
        # Points essentially on the boundary should clamp to each class's
        # own floor, not the true (near-zero) distance.
        lat = np.array([0.0, 0.0])
        long = np.array([-1e-6, 1e-6])
        ocean = np.array([1, 0])
        metric = mesh_module.coastal_distance_metric(lat, long, ocean,
                                                      mindist_coast=150e3, mindist_land=300e3, maxdist=600e3)
        assert metric[0] == pytest.approx(150e3)
        assert metric[1] == pytest.approx(300e3)

    def test_ceiling_is_respected_far_from_any_opposite_class_vertex(self):
        lat = np.array([0.0, 89.0])
        long = np.array([-90.0, 90.0])
        ocean = np.array([1, 0])
        metric = mesh_module.coastal_distance_metric(lat, long, ocean,
                                                      mindist_coast=1e3, mindist_land=1e3, maxdist=50e3)
        assert np.all(metric == pytest.approx(50e3))

    def test_metric_increases_away_from_boundary_on_ocean_side(self):
        lat, long, ocean = _lon_split_ocean_land(n_per_side=50, band_deg=5.0)
        is_ocean = ocean.astype(bool)
        metric = mesh_module.coastal_distance_metric(lat, long, ocean,
                                                      mindist_coast=1.0, mindist_land=1.0, maxdist=1e9)
        order = np.argsort(np.abs(long[is_ocean]))
        sorted_metric = metric[is_ocean][order]
        assert np.all(np.diff(sorted_metric) >= -1e-6)

    def test_matches_known_great_circle_separation(self):
        # Two points at a known angular separation - one ocean, one land -
        # so each is the other's only opposite-class candidate. Verifies
        # the chord -> great-circle conversion against an analytic value.
        lat = np.array([0.0, 0.0])
        long = np.array([-1.0, 1.0])
        ocean = np.array([1, 0])
        metric = mesh_module.coastal_distance_metric(lat, long, ocean,
                                                      mindist_coast=1.0, mindist_land=1.0, maxdist=1e9)
        expected = np.radians(2.0) * RADIUS_M
        assert metric[0] == pytest.approx(expected, rel=1e-3)
        assert metric[1] == pytest.approx(expected, rel=1e-3)

    def test_custom_radius_is_honoured(self):
        lat = np.array([0.0, 0.0])
        long = np.array([-1.0, 1.0])
        ocean = np.array([1, 0])
        metric_default = mesh_module.coastal_distance_metric(lat, long, ocean,
                                                              mindist_coast=1.0, mindist_land=1.0, maxdist=1e9)
        metric_half_radius = mesh_module.coastal_distance_metric(lat, long, ocean,
                                                                  mindist_coast=1.0, mindist_land=1.0, maxdist=1e9,
                                                                  radius=RADIUS_M / 2.0)
        assert metric_half_radius[0] == pytest.approx(metric_default[0] / 2.0, rel=1e-3)
