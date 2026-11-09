"""Automated geometric tests for tower_morphometry.py (synthetic shapes with
known answers).

Run from the project root:
    .venv\\Scripts\\python.exe -m pytest sensitivity_analysis\\tests

Author: Anonymous (blinded for peer review)
Project: 3D Geospatial Approaches for Geoheritage Documentation and Analysis
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import tower_morphometry as tm  # noqa: E402


def circle(r=10.0, n=720, cx=0.0, cy=0.0):
    """Points on a circle of radius r centred on (cx, cy)."""
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return cx + r * np.cos(t), cy + r * np.sin(t)


def ellipse(a, b, azimuth_deg=None, n=720):
    """Ellipse with semi-axis a along the azimuth direction (clockwise from
    north)."""
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    u = a * np.cos(t)  # along major axis
    v = b * np.sin(t)  # along minor axis
    if azimuth_deg is None:
        return u, v  # major along easting (azimuth 90)
    az = np.radians(azimuth_deg)
    d = np.array([np.sin(az), np.cos(az)])  # (E, N) unit vector of azimuth
    p = np.array([-d[1], d[0]])  # perpendicular
    x = u * d[0] + v * p[0]
    y = u * d[1] + v * p[1]
    return x, y


def rectangle(w, h, n_per_edge=60):
    """Points along the outline of a w by h rectangle centred on the origin."""
    xs = np.linspace(-w / 2, w / 2, n_per_edge)
    ys = np.linspace(-h / 2, h / 2, n_per_edge)
    edge = np.concatenate(
        [
            np.column_stack([xs, np.full_like(xs, -h / 2)]),
            np.column_stack([xs, np.full_like(xs, h / 2)]),
            np.column_stack([np.full_like(ys, -w / 2), ys]),
            np.column_stack([np.full_like(ys, w / 2), ys]),
        ]
    )
    return edge[:, 0], edge[:, 1]


def fill(x, y, n=4000, seed=0):
    """Uniform random points inside the convex hull of an outline."""
    from scipy.spatial import Delaunay

    rng = np.random.default_rng(seed)
    px = rng.uniform(x.min(), x.max(), n * 4)
    py = rng.uniform(y.min(), y.max(), n * 4)
    tri = Delaunay(np.column_stack([x, y]))
    inside = tri.find_simplex(np.column_stack([px, py])) >= 0
    return px[inside][:n], py[inside][:n]


def test_circle_area_and_isotropy():
    x, y = fill(*circle(10.0), n=20000)
    m = tm.band_metrics(x, y, 0, 2, min_points=1, min_coverage_deg=0)
    assert m.hull_area_m2 == pytest.approx(np.pi * 100, rel=0.02)
    assert m.elongation_poly == pytest.approx(1.0, abs=0.02)
    assert m.elongation_mbr == pytest.approx(1.0, abs=0.02)


def test_axis_aligned_ellipse_azimuth_is_east_west():
    x, y = fill(*ellipse(10.0, 5.0), n=20000)
    m = tm.band_metrics(x, y, 0, 2, min_points=1, min_coverage_deg=0)
    assert m.major_span_m == pytest.approx(20.0, rel=0.03)
    assert m.minor_span_m == pytest.approx(10.0, rel=0.03)
    assert m.elongation_mbr == pytest.approx(2.0, rel=0.05)
    assert m.azimuth_poly_deg == pytest.approx(90.0, abs=2.0)
    assert m.azimuth_mbr_deg == pytest.approx(90.0, abs=2.0)
    assert m.elongation_poly == pytest.approx(2.0, rel=0.05)


@pytest.mark.parametrize("az", [0.0, 30.0, 60.0, 120.0, 150.0, 175.0])
def test_rotated_ellipse_azimuth_recovered(az):
    x, y = fill(*ellipse(12.0, 4.0, azimuth_deg=az), n=20000)
    m = tm.band_metrics(x, y, 0, 2, min_points=1, min_coverage_deg=0)
    diff = abs(((m.azimuth_poly_deg - az) + 90) % 180 - 90)  # axial difference
    assert diff < 2.0, (az, m.azimuth_poly_deg)
    diff_mbr = abs(((m.azimuth_mbr_deg - az) + 90) % 180 - 90)
    assert diff_mbr < 3.0, (az, m.azimuth_mbr_deg)


def test_azimuth_is_axial_periodic():
    x1, y1 = fill(*ellipse(12.0, 4.0, azimuth_deg=40.0), n=8000, seed=1)
    x2, y2 = fill(*ellipse(12.0, 4.0, azimuth_deg=220.0), n=8000, seed=1)
    m1 = tm.band_metrics(x1, y1, 0, 2, min_points=1, min_coverage_deg=0)
    m2 = tm.band_metrics(x2, y2, 0, 2, min_points=1, min_coverage_deg=0)
    diff = abs(((m1.azimuth_poly_deg - m2.azimuth_poly_deg) + 90) % 180 - 90)
    assert diff < 2.0
    assert 0.0 <= m1.azimuth_poly_deg < 180.0


def test_rectangle_spans_and_elongation():
    x, y = fill(*rectangle(20.0, 4.0), n=20000)
    m = tm.band_metrics(x, y, 0, 2, min_points=1, min_coverage_deg=0)
    assert m.major_span_m == pytest.approx(20.0, rel=0.02)
    assert m.minor_span_m == pytest.approx(4.0, rel=0.02)
    assert m.hull_area_m2 == pytest.approx(80.0, rel=0.02)
    # polygon moments for a uniform rectangle: sqrt(lambda1/lambda2) = 5
    assert m.elongation_poly == pytest.approx(5.0, rel=0.03)


def test_projected_coordinates_do_not_change_orientation():
    # Regression test: UTM-scale offsets must not corrupt the polygon moments.
    x, y = fill(*ellipse(12.0, 4.0, azimuth_deg=30.0), n=8000, seed=3)
    m0 = tm.band_metrics(x, y, 0, 2, min_points=1, min_coverage_deg=0)
    m1 = tm.band_metrics(
        x + 444750.0, y + 3245840.0, 0, 2, min_points=1, min_coverage_deg=0
    )
    assert m1.elongation_poly == pytest.approx(m0.elongation_poly, rel=1e-6)
    diff = abs(((m1.azimuth_poly_deg - m0.azimuth_poly_deg) + 90) % 180 - 90)
    assert diff < 0.5
    assert m1.elongation_poly == pytest.approx(3.0, rel=0.05)


def test_hull_perimeter_of_square():
    x, y = rectangle(10.0, 10.0)
    m = tm.band_metrics(x, y, 0, 2, min_points=1, min_coverage_deg=0)
    assert m.hull_perimeter_m == pytest.approx(40.0, rel=1e-6)
    assert m.hull_area_m2 == pytest.approx(100.0, rel=1e-6)


def test_empty_band_is_flagged_not_crashed():
    m = tm.band_metrics(np.array([]), np.array([]), 0, 2)
    assert "empty" in m.flags and m.n_points == 0
    assert np.isnan(m.hull_area_m2)


def test_degenerate_collinear_points_flagged():
    x = np.linspace(0, 10, 50)
    y = np.zeros_like(x)
    m = tm.band_metrics(x, y, 0, 2, min_points=1, min_coverage_deg=0)
    assert "degenerate" in m.flags
    assert np.isnan(m.hull_area_m2)


def test_too_few_points_flagged_low_points():
    x, y = fill(*circle(10.0), n=50)
    m = tm.band_metrics(x, y, 0, 2, min_points=200)
    assert "low_points" in m.flags


def test_duplicate_xy_do_not_break_hull():
    x, y = fill(*circle(10.0), n=2000)
    x = np.concatenate([x, x[:500]])
    y = np.concatenate([y, y[:500]])
    m = tm.band_metrics(x, y, 0, 2, min_points=1, min_coverage_deg=0)
    assert np.isfinite(m.hull_area_m2)


def test_coverage_full_ring_vs_one_sided_surface():
    # Surface samples (ring), as in real photogrammetric bands.
    x, y = circle(10.0, n=720)
    _, cov_full = tm.sector_coverage(x, y)
    assert cov_full == pytest.approx(360.0)
    # One visible side only: the front half of the ring. The gap behind the
    # centroid should reduce coverage well below 360 degrees.
    keep = x > 0
    _, cov_front = tm.sector_coverage(x[keep], y[keep])
    assert 200.0 <= cov_front <= 280.0


def test_coverage_of_filled_disk_is_full_not_a_completeness_measure():
    # Documents a known limit: for filled point sets every direction is
    # occupied,
    # so sector coverage does not detect a missing side. Use it on
    # surface samples.
    x, y = fill(*circle(10.0), n=20000)
    keep = x > x.mean()
    _, cov_half_filled = tm.sector_coverage(x[keep], y[keep])
    assert cov_half_filled == pytest.approx(360.0)


def test_two_separated_cylinders_are_two_major_components():
    x1, y1 = fill(*circle(5.0), n=2000)
    x2, y2 = fill(*circle(5.0, cx=40.0), n=2000, seed=2)
    x = np.concatenate([x1, x2])
    y = np.concatenate([y1, y2])
    assert tm.n_major_components(x, y, link_m=0.5) == 2
    m = tm.band_metrics(x, y, 0, 2, min_points=1, min_coverage_deg=0)
    assert "multiple_components" in m.flags


def test_two_touching_cylinders_are_one_major_component():
    x1, y1 = fill(*circle(5.0), n=2000)
    x2, y2 = fill(*circle(5.0, cx=9.0), n=2000, seed=2)
    x = np.concatenate([x1, x2])
    y = np.concatenate([y1, y2])
    assert tm.n_major_components(x, y, link_m=0.5) == 1
    m = tm.band_metrics(x, y, 0, 2, min_points=1, min_coverage_deg=0)
    assert "multiple_components" not in m.flags


def test_small_fragments_do_not_count_as_objects():
    x1, y1 = fill(*circle(5.0), n=4000)
    # ten isolated 5-point fragments far away: reconstruction debris, not
    # objects
    fx = np.arange(10) * 3.0 + 100.0
    x = np.concatenate([x1, fx.repeat(5)])
    y = np.concatenate([y1, np.zeros(50)])
    assert tm.n_components(x, y, link_m=0.5) > 1
    assert tm.n_major_components(x, y, link_m=0.5) == 1


def test_band_edges_absolute_multiples():
    e = tm.band_edges(807.65, 949.22, 2.0)
    assert e[0] == 806.0 and e[-1] >= 949.22
    assert np.allclose(np.diff(e), 2.0)
    e5 = tm.band_edges(807.65, 949.22, 5.0)
    assert np.allclose(e5 % 5.0, 0.0)


def test_slice_indices_covers_every_point_once():
    rng = np.random.default_rng(4)
    z = rng.uniform(100, 110, 5000)
    edges = tm.band_edges(z.min(), z.max(), 1.0)
    seen = np.zeros(z.size, dtype=int)
    for _, _, idx in tm.slice_indices(z, edges):
        seen[idx] += 1
    assert np.all(seen == 1)


def test_stratified_subsample_size_and_reproducibility():
    x, y = fill(*circle(10.0), n=5000)
    a = tm.stratified_subsample(x, y, 500, np.random.default_rng(7))
    b = tm.stratified_subsample(x, y, 500, np.random.default_rng(7))
    assert np.array_equal(a, b)
    assert abs(len(a) - 500) <= 25


def test_axial_sd_of_constant_is_zero_and_of_uniform_is_large():
    assert tm.axial_sd_deg(np.full(50, 42.0)) == pytest.approx(0.0, abs=1e-6)
    rng = np.random.default_rng(5)
    assert tm.axial_sd_deg(rng.uniform(0, 180, 5000)) > 40.0
