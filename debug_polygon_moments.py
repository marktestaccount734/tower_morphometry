"""Debug: compare polygon-moment covariance with the covariance of uniformly
sampled points
inside the same hull, for one J_T01 band. Read-only on the cache.

Diagnostic script, not part of the analysis pipeline.

Author: Anonymous (blinded for peer review)
Project: 3D Geospatial Approaches for Geoheritage Documentation and Analysis
"""

import sys
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull, Delaunay

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import tower_morphometry as tm  # noqa: E402

d = np.load(HERE / "cache" / "J_T01.npz")
x, y, z = d["x"], d["y"], d["z"]
edges = tm.band_edges(z.min(), z.max(), 2.0)
for lo, hi, idx in tm.slice_indices(z, edges):
    if abs(lo - 940.0) < 1e-6:
        break
pts = np.column_stack([x[idx], y[idx]])
hull = ConvexHull(pts)
poly = pts[hull.vertices]
print("hull vertices", len(poly), "area", hull.volume)
print(
    "signed area via shoelace",
    0.5
    * np.sum(
        poly[:, 0] * np.roll(poly[:, 1], -1)
        - np.roll(poly[:, 0], -1) * poly[:, 1]
    ),
)
A, c, cov = tm._polygon_moments(poly)
print("polygon centroid", c, "point mean", pts.mean(axis=0))
print("polygon cov\n", cov)
print("polygon eig", np.linalg.eigh(cov)[0])

# sampled reference
tri = Delaunay(poly)
rng = np.random.default_rng(0)
u = rng.uniform(poly.min(axis=0), poly.max(axis=0), size=(400000, 2))
u = u[tri.find_simplex(u) >= 0]
cov_s = np.cov(u.T)
print("sampled cov\n", cov_s)
print("sampled eig", np.linalg.eigh(cov_s)[0])
print("mbr", tm._min_bounding_rectangle(poly))
