"""Debug: list the largest connected components of a Q_T01 band, to see whether
the
multiple-component flag reflects real separate objects or sparse-sampling
fragments.
Read-only on the cache.

Diagnostic script, not part of the analysis pipeline.

Author: Anonymous (blinded for peer review)
Project: 3D Geospatial Approaches for Geoheritage Documentation and Analysis
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import tower_morphometry as tm  # noqa: E402

d = np.load(HERE / "cache" / "Q_T01.npz")
x, y, z = d["x"], d["y"], d["z"]
for band_lo in (760.0, 840.0, 900.0):
    sel = (z >= band_lo) & (z < band_lo + 2.0)
    xs, ys = x[sel], y[sel]
    pts = np.column_stack([xs, ys])
    from scipy.spatial import cKDTree

    tree = cKDTree(pts)
    pairs = tree.query_pairs(0.5, output_type="ndarray")
    n = len(pts)
    parent = np.arange(n)

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
    labs, counts = np.unique(roots, return_counts=True)
    order = np.argsort(counts)[::-1][:6]
    print(
        f"band {band_lo:.0f}-{band_lo + 2:.0f} m: "
        f"n={n}, components={len(labs)}, "
        f"major(>=5%,>=50)={tm.n_major_components(xs, ys)}"
    )
    for k in order:
        m = roots == labs[k]
        print(
            f"   size {counts[k]:6d}  "
            f"centroid ({xs[m].mean():.1f}, {ys[m].mean():.1f})  "
            f"extent {np.ptp(xs[m]):.1f} x {np.ptp(ys[m]):.1f} m"
        )
