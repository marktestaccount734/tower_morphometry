"""Horizontal cross-section morphometry for tower point clouds.

Conventions
- Coordinates: projected metres, EPSG:32649 (WGS 84 / UTM 49N).
- Elevation: actual Z in metres. Band edges are integer multiples of the band
  thickness in absolute Z, so bands line up across towers and thickness runs.
- Azimuth: axial, degrees clockwise from north, in [0, 180). Direction (dx, dy)
  with dx = easting, dy = northing gives azimuth atan2(dx, dy) mod 180.
- Convex-hull area is the area of the envelope of sampled points. It is not the
  solid rock cross-section and must not be used for volume or mass.

Methods
- Hull area and perimeter: scipy ConvexHull (2-D volume and area).
- Spans and elongation (primary for spans): oriented minimum-area bounding
  rectangle, found by testing each hull edge direction (exact for convex
  polygons).
- Azimuth and polygon elongation (primary for orientation): second moments of
the
  hull polygon treated as a uniform area (Green's-theorem formulas), not PCA on
  hull vertices. Elongation = sqrt(lambda1 / lambda2).
- Angular coverage: number of occupied 10-degree sectors around the point
centroid,
  reported in degrees (0 to 360).
- Multi-object check: connected components of points linked within a distance
  threshold (default 0.5 m).

QC flags are recorded, never used to delete points.

Author: Anonymous (blinded for peer review)
Project: 3D Geospatial Approaches for Geoheritage Documentation and Analysis
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import ConvexHull, cKDTree

try:  # scipy >= 1.11 exposes QhullError in the public namespace
    from scipy.spatial import QhullError
except ImportError:  # pragma: no cover
    from scipy.spatial.qhull import QhullError

SECTOR_DEG = 10.0
N_SECTORS = 36
TOL = 1e-12


# ---------------------------------------------------------------- I/O


def read_tower(path: str):
    """Read X, Y, Z and 8-bit RGB from a LAS file (read-only)."""
    import laspy

    las = laspy.read(path)
    x = np.asarray(las.x, dtype=np.float64)
    y = np.asarray(las.y, dtype=np.float64)
    z = np.asarray(las.z, dtype=np.float64)
    rgb = np.column_stack([las.red, las.green, las.blue]).astype(np.uint16)
    rgb8 = (rgb >> 8).astype(np.float64)  # 16-bit RGB holds 8-bit values x 256
    return x, y, z, rgb8


# ---------------------------------------------------------------- bands


def band_edges(zmin: float, zmax: float, thickness: float) -> np.ndarray:
    """Band edges at integer multiples of thickness covering [zmin, zmax]."""
    if thickness <= 0:
        raise ValueError("thickness must be positive")
    lo = np.floor(zmin / thickness) * thickness
    hi = np.ceil(zmax / thickness) * thickness
    n = int(round((hi - lo) / thickness)) + 1
    return lo + thickness * np.arange(n)


def slice_indices(z: np.ndarray, edges: np.ndarray):
    """Yield (lower, upper, index array) for each band. Upper edge is
    exclusive,
    except the last band, which includes its upper edge."""
    order = np.argsort(z, kind="stable")
    zs = z[order]
    cuts = np.searchsorted(zs, edges, side="left")
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        a, b = cuts[i], cuts[i + 1]
        idx = order[a:b]
        if i == len(edges) - 2:
            # include points exactly on the top edge in the last band
            extra = np.searchsorted(zs, hi, side="right")
            idx = order[a:extra]
        yield lo, hi, idx


# ---------------------------------------------------------------- geometry


def _polygon_moments(poly: np.ndarray):
    """Area, centroid and covariance of a simple polygon with uniform density.
    poly: (n, 2) vertices, counter-clockwise, first vertex not repeated.

    The polygon is translated to its vertex mean before the moments are summed.
    Projected coordinates are of order 1e5 to 1e6 m, so unshifted sums of cross
    products cancel catastrophically in float64. Covariance is translation
    invariant,
    so the shift changes only numerical error, not the result."""
    origin = poly.mean(axis=0)
    shifted = poly - origin
    x = shifted[:, 0]
    y = shifted[:, 1]
    x1 = np.roll(x, -1)
    y1 = np.roll(y, -1)
    cross = x * y1 - x1 * y
    A2 = np.sum(cross)  # twice the signed area
    A = 0.5 * A2
    if abs(A) < TOL:
        return None
    cx = np.sum((x + x1) * cross) / (3.0 * A2)
    cy = np.sum((y + y1) * cross) / (3.0 * A2)
    Sxx = np.sum(cross * (x * x + x * x1 + x1 * x1)) / 12.0
    Syy = np.sum(cross * (y * y + y * y1 + y1 * y1)) / 12.0
    Sxy = np.sum(cross * (x * y1 + 2 * x * y + 2 * x1 * y1 + x1 * y)) / 24.0
    # central second moments per unit area
    cxx = Sxx / A - cx * cx
    cyy = Syy / A - cy * cy
    cxy = Sxy / A - cx * cy
    cov = np.array([[cxx, cxy], [cxy, cyy]])
    return A, (cx + origin[0], cy + origin[1]), cov


def _axial_azimuth(dx: float, dy: float) -> float:
    """Axial azimuth of direction (dx, dy): degrees clockwise from north,
    in [0, 180).
    dx is easting and dy is northing, so (0, 1) gives 0 and (1, 0) gives 90."""
    return float(np.degrees(np.arctan2(dx, dy)) % 180.0)


def _min_bounding_rectangle(poly: np.ndarray):
    """Oriented minimum-area bounding rectangle of a convex polygon.
    Returns (major, minor, azimuth_of_major) or None."""
    best = None
    n = len(poly)
    for i in range(n):
        p0 = poly[i]
        p1 = poly[(i + 1) % n]
        e = p1 - p0
        length = np.hypot(*e)
        if length < TOL:
            continue
        u = e / length
        v = np.array([-u[1], u[0]])
        pu = poly @ u
        pv = poly @ v
        wu = pu.max() - pu.min()
        wv = pv.max() - pv.min()
        area = wu * wv
        if best is None or area < best[0]:
            best = (area, wu, wv, u)
    if best is None:
        return None
    _, wu, wv, u = best
    if wu >= wv:
        major, minor, d = wu, wv, u
    else:
        major, minor, d = wv, wu, np.array([-u[1], u[0]])
    return major, minor, _axial_azimuth(d[0], d[1])


def sector_coverage(x: np.ndarray, y: np.ndarray):
    """Occupied 10-degree sectors around the point centroid."""
    if len(x) == 0:
        return 0, 0.0
    dx = x - x.mean()
    dy = y - y.mean()
    ang = np.degrees(np.arctan2(dx, dy)) % 360.0
    sec = np.floor(ang / SECTOR_DEG).astype(int) % N_SECTORS
    occupied = int(np.unique(sec).size)
    return occupied, occupied * SECTOR_DEG


def component_sizes(
    x: np.ndarray, y: np.ndarray, link_m: float = 0.5
) -> np.ndarray:
    """Sizes of connected components when points closer than link_m are
    linked."""
    n = len(x)
    if n == 0:
        return np.array([], dtype=int)
    pts = np.column_stack([x, y])
    tree = cKDTree(pts)
    pairs = tree.query_pairs(link_m, output_type="ndarray")
    parent = np.arange(n)

    # Union-find: walk parent links to the root of a's component,
    # halving the path as we go.
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    roots = np.array([find(i) for i in range(n)])
    _, counts = np.unique(roots, return_counts=True)
    return np.sort(counts)[::-1]


def n_components(x: np.ndarray, y: np.ndarray, link_m: float = 0.5) -> int:
    """Raw number of connected components (sensitive to photogrammetric
    holes)."""
    return int(component_sizes(x, y, link_m).size)


def n_major_components(
    x, y, link_m: float = 0.5, min_frac: float = 0.05, min_pts: int = 50
) -> int:
    """Components holding at least min_frac of the band's points and at
    least min_pts points.
    Used for the multiple-object flag. Small fragments caused by
    reconstruction holes are ignored.
    """
    sizes = component_sizes(x, y, link_m)
    if sizes.size == 0:
        return 0
    thr = max(min_pts, min_frac * len(x))
    return int(np.sum(sizes >= thr))


# ---------------------------------------------------------------- band metrics


@dataclass
class BandResult:
    """Morphometry of one horizontal band of a tower.

    Missing or undefined values are NaN. ``flags`` lists QC issues (for example
    'low_points', 'low_coverage', 'multiple_components', 'azimuth_unreliable')
    and is never used to delete points."""

    lower_m: float
    upper_m: float
    midpoint_m: float
    n_points: int
    hull_area_m2: float = np.nan
    hull_perimeter_m: float = np.nan
    major_span_m: float = np.nan
    minor_span_m: float = np.nan
    elongation_mbr: float = np.nan
    azimuth_mbr_deg: float = np.nan
    elongation_poly: float = np.nan
    azimuth_poly_deg: float = np.nan
    poly_centroid_x: float = np.nan
    poly_centroid_y: float = np.nan
    point_centroid_x: float = np.nan
    point_centroid_y: float = np.nan
    angular_coverage_deg: float = np.nan
    n_components: int = (
        0  # major components (>= 5% of points and >= 50 points)
    )
    n_components_raw: int = 0  # all connected components at the link distance
    flags: list = field(default_factory=list)


def band_metrics(
    x: np.ndarray,
    y: np.ndarray,
    lower: float,
    upper: float,
    min_points: int = 200,
    min_coverage_deg: float = 180.0,
    elong_threshold: float = 1.2,
    link_m: float = 0.5,
) -> BandResult:
    """Morphometry of one horizontal band. Failures are recorded as flags."""
    n = int(len(x))
    res = BandResult(lower, upper, 0.5 * (lower + upper), n)
    if n == 0:
        res.flags.append("empty")
        return res
    res.point_centroid_x = float(x.mean())
    res.point_centroid_y = float(y.mean())
    res.n_components_raw = n_components(x, y, link_m)
    res.n_components = n_major_components(x, y, link_m)
    occ, cov = sector_coverage(x, y)
    res.angular_coverage_deg = cov
    if res.n_components > 1:
        res.flags.append("multiple_components")
    if n < min_points:
        res.flags.append("low_points")
    if cov < min_coverage_deg:
        res.flags.append("low_coverage")
    if n < 3:
        res.flags.append("degenerate")
        return res
    pts = np.column_stack([x, y])
    try:
        hull = ConvexHull(pts)
    except (QhullError, ValueError):
        res.flags.append("degenerate")
        return res
    poly = pts[hull.vertices]  # scipy returns counter-clockwise order in 2-D
    res.hull_area_m2 = float(hull.volume)
    res.hull_perimeter_m = float(hull.area)

    mbr = _min_bounding_rectangle(poly)
    if mbr is not None and mbr[1] > TOL:
        res.major_span_m, res.minor_span_m = float(mbr[0]), float(mbr[1])
        res.elongation_mbr = float(mbr[0] / mbr[1])
        res.azimuth_mbr_deg = float(mbr[2])

    pm = _polygon_moments(poly)
    if pm is not None:
        A, (cx, cy), cov = pm
        res.poly_centroid_x, res.poly_centroid_y = float(cx), float(cy)
        w, vecs = np.linalg.eigh(cov)
        if w[0] > TOL:
            # eigh sorts eigenvalues ascending, so w[1] is the major axis.
            res.elongation_poly = float(np.sqrt(w[1] / w[0]))
            v1 = vecs[:, 1]
            res.azimuth_poly_deg = _axial_azimuth(v1[0], v1[1])
    if (
        np.isfinite(res.elongation_poly)
        and res.elongation_poly < elong_threshold
    ):
        res.flags.append("azimuth_unreliable")
    if not np.isfinite(res.azimuth_poly_deg):
        res.flags.append("azimuth_undefined")
    return res


def slice_tower(
    x,
    y,
    z,
    thickness: float = 2.0,
    min_points: int = 200,
    min_coverage_deg: float = 180.0,
    elong_threshold: float = 1.2,
    link_m: float = 0.5,
):
    """Morphometry for every band of one tower. Returns a list of
    BandResult."""
    edges = band_edges(float(z.min()), float(z.max()), thickness)
    out = []
    for lo, hi, idx in slice_indices(z, edges):
        out.append(
            band_metrics(
                x[idx],
                y[idx],
                lo,
                hi,
                min_points=min_points,
                min_coverage_deg=min_coverage_deg,
                elong_threshold=elong_threshold,
                link_m=link_m,
            )
        )
    return out


# ---------------------------------------------------------------- subsampling


def stratified_subsample(
    x, y, n_target: int, rng: np.random.Generator, grid: int = 4
):
    """Spatially stratified subsample: the same fraction is drawn from each
    cell of a
    grid x grid partition of the band's bounding box. Returns sorted indices.

    Counts per cell are rounded, so the realised total can differ slightly from
    n_target. Cells with no points, or with fewer than one selected point after
    rounding, contribute nothing. If n_target >= len(x), all indices are
    returned.
    """
    n = len(x)
    if n_target >= n:
        return np.arange(n)
    frac = n_target / n
    x0, x1 = x.min(), x.max()
    y0, y1 = y.min(), y.max()
    ix = np.clip(
        ((x - x0) / max(x1 - x0, TOL) * grid).astype(int), 0, grid - 1
    )
    iy = np.clip(
        ((y - y0) / max(y1 - y0, TOL) * grid).astype(int), 0, grid - 1
    )
    cell = ix * grid + iy
    chosen = []
    for c in np.unique(cell):
        members = np.flatnonzero(cell == c)
        k = int(round(frac * len(members)))
        if k > 0:
            chosen.append(
                # Draw k members of this cell without replacement.
                rng.choice(members, size=min(k, len(members)), replace=False)
            )
    idx = np.concatenate(chosen) if chosen else np.array([], dtype=int)
    return np.sort(idx)


def axial_circular_stats(az_deg: np.ndarray):
    """Mean axial direction and resultant length for 180-degree periodic
    angles."""
    # Double the angles so that an axis at 10 and one at 190 degrees are
    # treated as the same direction.
    a = np.radians(np.asarray(az_deg, dtype=float)) * 2.0
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan, np.nan
    c, s = np.cos(a).mean(), np.sin(a).mean()
    r = np.hypot(c, s)
    # Halve the doubled mean to return to the axial 0-180 degree scale.
    mean = (np.degrees(np.arctan2(s, c)) / 2.0) % 180.0
    return float(mean), float(r)


def axial_sd_deg(az_deg: np.ndarray) -> float:
    """Circular standard deviation for axial data, in degrees."""
    a = np.asarray(az_deg, dtype=float)
    a = a[np.isfinite(a)]
    if a.size < 2:
        return np.nan
    _, r = axial_circular_stats(a)
    r = min(max(r, 1e-12), 1.0)
    # Circular standard deviation (Fisher, 1993), halved to the axial scale.
    return float(np.degrees(np.sqrt(-2.0 * np.log(r))) / 2.0)
