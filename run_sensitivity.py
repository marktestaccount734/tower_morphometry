"""Run sensitivity phases 4-10 in order. Requires the baseline
(sensitivity_analysis.py baseline).

    python sensitivity_analysis\\run_sensitivity.py [phase ...]

Phases: thickness, pointcount, coverage, azimuth, rockbase, crosstower, all

Author: Anonymous (blinded for peer review)
Project: 3D Geospatial Approaches for Geoheritage Documentation and Analysis
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sensitivity_analysis as base  # noqa: E402
import sensitivity_phases as ph  # noqa: E402

RESULTS = HERE / "results"


def main(argv):
    """Run the requested sensitivity phases in dependency order.

    Phases that need subsampling results (azimuth reliability, cross-tower and
    leave-one-tower-out) cause the subsampling phase to run in the same call,
    because nothing is cached between calls."""
    phases = argv[1:] or ["all"]
    run_all = "all" in phases
    cfg = base.load_config()
    RESULTS.mkdir(parents=True, exist_ok=True)
    clouds = {}
    for tower in cfg["towers"]:
        clouds[tower] = base.load_cloud(tower, cfg)
    baseline = pd.read_csv(RESULTS / "baseline_measurements.csv")
    base.log_event("sensitivity_start", phases=phases)

    def timed(name, fn, *a):
        t0 = time.time()
        out = fn(*a)
        base.log_event(
            "phase_done", phase=name, seconds=round(time.time() - t0, 1)
        )
        print(f"{name}: done ({time.time() - t0:.1f} s)", flush=True)
        return out

    point_ret = None
    sub_summary = pd.DataFrame()
    if run_all or "thickness" in phases:
        timed("band_thickness", ph.phase_band_thickness, clouds, cfg)
    if run_all or "pointcount" in phases or "crosstower" in phases:
        point_ret = timed(
            "point_count_retention",
            ph.phase_point_count_retention,
            baseline,
            cfg,
        )
    if (
        run_all
        or "pointcount" in phases
        or "crosstower" in phases
        or "azimuth" in phases
    ):
        _, sub_summary = timed(
            "subsampling", ph.phase_subsampling, clouds, baseline, cfg
        )
    if run_all or "coverage" in phases:
        timed(
            "angular_coverage",
            ph.phase_angular_coverage,
            clouds,
            baseline,
            cfg,
        )
    if run_all or "azimuth" in phases:
        timed(
            "azimuth_reliability",
            ph.phase_azimuth_reliability,
            baseline,
            sub_summary,
            cfg,
        )
    if run_all or "rockbase" in phases:
        timed("rock_base", ph.phase_rock_base, clouds, cfg)
    if run_all or "crosstower" in phases:
        if point_ret is None:
            point_ret = ph.phase_point_count_retention(baseline, cfg)
        cross = timed(
            "cross_tower",
            ph.phase_cross_tower,
            baseline,
            point_ret,
            sub_summary,
            cfg,
        )
        timed("leave_one_tower_out", ph.phase_leave_one_tower_out, cross, cfg)
    base.log_event("sensitivity_end", phases=phases)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
