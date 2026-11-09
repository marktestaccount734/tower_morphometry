"""Sensitivity phases 4-10 for tower cross-sectional morphometry.

Each function writes a CSV into sensitivity_analysis/results and returns a
DataFrame.
Inputs are read from the cache (written by sensitivity_analysis.py from the LAS
files).

Author: Anonymous (blinded for peer review)
Project: 3D Geospatial Approaches for Geoheritage Documentation and Analysis
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import tower_morphometry as tm

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


# ---------------------------------------------------------------- shared
# helpers


def polygon_metrics(x: np.ndarray, y: np.ndarray) -> dict:
    """Light metrics for replicate runs: area, polygon elongation and axial
    azimuth.
    Skips component analysis, which is expensive and not needed per replicate.
    """
    out = {"area": np.nan, "elong": np.nan, "az": np.nan, "mbr_elong": np.nan}
    if len(x) < 3:
        return out
    from scipy.spatial import ConvexHull
    from scipy.spatial import QhullError

    pts = np.column_stack([x, y])
    try:
        hull = ConvexHull(pts)
    except (QhullError, ValueError):
        return out
    poly = pts[hull.vertices]
    out["area"] = float(hull.volume)
    pm = tm._polygon_moments(poly)
    if pm is not None:
        _, _, cov = pm
        w, vecs = np.linalg.eigh(cov)
        if w[0] > 1e-12:
            # Same axis-ratio definition as tower_morphometry.band_metrics.
            out["elong"] = float(np.sqrt(w[1] / w[0]))
            v1 = vecs[:, 1]
            out["az"] = tm._axial_azimuth(v1[0], v1[1])
    mbr = tm._min_bounding_rectangle(poly)
    if mbr is not None and mbr[1] > 1e-12:
        out["mbr_elong"] = float(mbr[0] / mbr[1])
    return out


def band_slices(x, y, z, thickness: float):
    """List of (lower, upper, index array) for all bands, with
    absolute-origin edges."""
    edges = tm.band_edges(float(z.min()), float(z.max()), thickness)
    return [(lo, hi, idx) for lo, hi, idx in tm.slice_indices(z, edges)]


def band_table(
    tower: str,
    x,
    y,
    z,
    thickness: float,
    prov: dict,
    rep_cfg: dict,
    base_floor=None,
) -> pd.DataFrame:
    """Morphometry table for every band of one tower at one band thickness.

    Representative bands must meet the minimum points, minimum angular coverage
    and single-component rule from rep_cfg, and lie at or above base_floor
    (the visual base) when one is given."""
    rows = []
    for lo, hi, idx in band_slices(x, y, z, thickness):
        r = tm.band_metrics(
            x[idx],
            y[idx],
            lo,
            hi,
            min_points=prov["min_points_per_band"],
            min_coverage_deg=prov["min_angular_coverage_deg"],
            elong_threshold=prov["elongation_threshold"],
        )
        rows.append(
            {
                "tower": tower,
                "band_thickness_m": thickness,
                "lower_m": lo,
                "upper_m": hi,
                "n_points": r.n_points,
                "hull_area_m2": r.hull_area_m2,
                "elongation_poly": r.elongation_poly,
                "azimuth_poly_deg": r.azimuth_poly_deg,
                "elongation_mbr": r.elongation_mbr,
                "azimuth_mbr_deg": r.azimuth_mbr_deg,
                "angular_coverage_deg": r.angular_coverage_deg,
                "n_major_components": r.n_components,
                "flags": ";".join(r.flags),
                "representative": bool(
                    r.n_points >= rep_cfg["min_points"]
                    and r.angular_coverage_deg
                    >= rep_cfg["min_angular_coverage_deg"]
                    and r.n_components == 1
                    and not np.isnan(r.hull_area_m2)
                    and (base_floor is None or r.lower_m >= base_floor)
                ),
            }
        )
    return pd.DataFrame(rows)


def select_sections(df: pd.DataFrame, max_per_tower: int = 6) -> pd.DataFrame:
    """Representative sections, spread evenly over elevation
    (deterministic)."""
    out = []
    for tower, g in df[df.representative].groupby("tower"):
        g = g.sort_values("lower_m").reset_index(drop=True)
        k = min(max_per_tower, len(g))
        if k == 0:
            continue
        # Pick up to k sections spread evenly across the elevation range.
        pick = np.unique(np.linspace(0, len(g) - 1, k).round().astype(int))
        out.append(g.iloc[pick])
    if not out:
        return pd.DataFrame()
    sel = pd.concat(out, ignore_index=True)
    sel["section_id"] = np.arange(len(sel))
    return sel


def rock_share_per_band(r, g, b, slices) -> list:
    """Share of rock-coloured points per band, using the Excess Green rule
    as a heuristic.
    Rock here means 'not flagged as vegetation'. Not a validated classifier."""
    exg = 2 * g - r - b
    veg = (exg > 20) & (g > 1.08 * r) & (g > 1.08 * b)
    return [
        float(1.0 - veg[idx].mean()) if idx.size else np.nan
        for _, _, idx in slices
    ]


# ---------------------------------------------------------------- phase 4:
# band thickness


def phase_band_thickness(clouds: dict, cfg: dict) -> pd.DataFrame:
    """Phase 4. Repeat the band analysis at each candidate band thickness
    and write
    band_thickness_results.csv and band_thickness_summary.csv."""
    prov = dict(cfg["provisional_values"])
    frames = []
    for thick in cfg["candidates"]["band_thickness_m"]:
        for tower, (x, y, z, _) in clouds.items():
            frames.append(
                band_table(
                    tower,
                    x,
                    y,
                    z,
                    float(thick),
                    prov,
                    cfg["representative_section"],
                    base_floor=cfg["towers"][tower].get(
                        "exclude_below_base_m"
                    ),
                )
            )
    out = pd.concat(frames, ignore_index=True)
    out.to_csv(RESULTS / "band_thickness_results.csv", index=False)
    summ = (
        out.groupby(["tower", "band_thickness_m"])
        .agg(
            bands=("n_points", "size"),
            representative_bands=("representative", "sum"),
            median_area_representative=(
                "hull_area_m2",
                lambda s: (
                    float(np.nanmedian(s[out.loc[s.index, "representative"]]))
                    if out.loc[s.index, "representative"].any()
                    else np.nan
                ),
            ),
            median_elong_representative=(
                "elongation_poly",
                lambda s: (
                    float(np.nanmedian(s[out.loc[s.index, "representative"]]))
                    if out.loc[s.index, "representative"].any()
                    else np.nan
                ),
            ),
        )
        .reset_index()
    )
    summ.to_csv(RESULTS / "band_thickness_summary.csv", index=False)
    return out


# ---------------------------------------------------------------- phase 5:
# point count


def phase_point_count_retention(
    baseline: pd.DataFrame, cfg: dict
) -> pd.DataFrame:
    """Phase 5a. For each candidate point count, the fraction of all bands
    and of
    representative bands that keep at least that many points. Writes
    point_count_retention.csv."""
    rows = []
    rep = baseline[baseline.representative]
    for thr in cfg["candidates"]["min_points_per_band"]:
        for tower, g in baseline.groupby("tower"):
            gr = rep[rep.tower == tower]
            rows.append(
                {
                    "tower": tower,
                    "min_points": thr,
                    "bands_total": len(g),
                    "bands_retained_all": int((g.n_points >= thr).sum()),
                    "retained_fraction_all": float((g.n_points >= thr).mean()),
                    "representative_total": len(gr),
                    "representative_retained": int((gr.n_points >= thr).sum()),
                    "representative_retained_fraction": (
                        float((gr.n_points >= thr).mean())
                        if len(gr)
                        else np.nan
                    ),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "point_count_retention.csv", index=False)
    return out


def phase_subsampling(clouds: dict, baseline: pd.DataFrame, cfg: dict):
    """Phase 5b. Repeated stratified subsampling of selected representative
    sections.

    For each selected section and candidate point count, draws nrep
    replicates with
    stratified_subsample, each with a seed derived from the base seed, the
    section,
    the size and the replicate number. Compares each replicate's hull area and
    axial azimuth with the full-resolution section. Writes
    subsampling_replicates.csv
    and subsampling_results.csv. Returns (replicates, summary)."""
    prov = cfg["provisional_values"]
    thick = float(prov["band_thickness_m"])
    sel = select_sections(baseline[baseline.representative].copy())
    if sel.empty:
        return pd.DataFrame(), pd.DataFrame()
    nrep = int(cfg["project"]["n_subsample_reps"])
    base_seed = int(cfg["project"]["random_seed"])
    rep_rows, sum_rows = [], []
    for _, s in sel.iterrows():
        tower = s["tower"]
        x, y, z, _ = clouds[tower]
        lo = float(s["lower_m"])
        idx = None
        for l2, h2, ii in band_slices(x, y, z, thick):
            if abs(l2 - lo) < 1e-9:
                idx = ii
                break
        if idx is None:
            continue
        xs, ys = x[idx], y[idx]
        ref = polygon_metrics(xs, ys)
        for size in cfg["candidates"]["subsample_sizes"]:
            if size >= len(xs):
                continue
            for rep in range(nrep):
                seed = (
                    base_seed
                    + 100000 * int(s["section_id"])
                    + 1000 * size
                    + rep
                )
                rng = np.random.default_rng(seed)
                sub = tm.stratified_subsample(xs, ys, size, rng)
                m = polygon_metrics(xs[sub], ys[sub])
                rep_rows.append(
                    {
                        "tower": tower,
                        "section_lower_m": lo,
                        "section_id": int(s["section_id"]),
                        "n_full": len(xs),
                        "requested_n": int(size),
                        "subsample_n": int(sub.size),
                        "rep": rep,
                        "seed": seed,
                        "area": m["area"],
                        "elong": m["elong"],
                        "az": m["az"],
                        "ref_area": ref["area"],
                        "ref_elong": ref["elong"],
                        "ref_az": ref["az"],
                    }
                )
    reps = pd.DataFrame(rep_rows)
    reps.to_csv(RESULTS / "subsampling_replicates.csv", index=False)
    if reps.empty:
        return reps, pd.DataFrame()
    for (tower, sec, size), g in reps.groupby(
        ["tower", "section_id", "requested_n"]
    ):
        ref_area = g.ref_area.iloc[0]
        ref_elong = g.ref_elong.iloc[0]
        ref_az = g.ref_az.iloc[0]
        sum_rows.append(
            {
                "tower": tower,
                "section_id": sec,
                "section_lower_m": g.section_lower_m.iloc[0],
                "n_full": g.n_full.iloc[0],
                # requested points per band (candidate value)
                "subsample_n": size,
                "mean_realised_n": float(g.subsample_n.mean()),
                "reps": len(g),
                "ref_area_m2": ref_area,
                "area_mean_m2": g.area.mean(),
                "area_bias_fraction": (
                    (g.area.mean() - ref_area) / ref_area
                    if ref_area
                    else np.nan
                ),
                "area_cv_fraction": (
                    g.area.std(ddof=1) / g.area.mean()
                    if g.area.mean()
                    else np.nan
                ),
                "elong_bias": g.elong.mean() - ref_elong,
                "elong_cv_fraction": (
                    g.elong.std(ddof=1) / g.elong.mean()
                    if g.elong.mean()
                    else np.nan
                ),
                "azimuth_axial_sd_deg": tm.axial_sd_deg(g.az.values),
                "ref_azimuth_deg": ref_az,
            }
        )
    summ = pd.DataFrame(sum_rows)
    summ.to_csv(RESULTS / "subsampling_results.csv", index=False)
    return reps, summ


# ---------------------------------------------------------------- phase 6:
# angular coverage


def phase_angular_coverage(
    clouds: dict, baseline: pd.DataFrame, cfg: dict
) -> pd.DataFrame:
    """Phase 6. Two parts. (1) Retention of bands at each candidate
    angular-coverage
    threshold. (2) Simulated missing sectors: on representative sections,
    removes a
    random angular gap and measures the error in area, elongation and azimuth.
    Writes angular_coverage_results.csv, angular_coverage_simulation.csv and
    angular_coverage_error_summary.csv. Surface samples only; this does not
    detect
    missing sides of filled point sets."""
    prov = cfg["provisional_values"]
    thick = float(prov["band_thickness_m"])
    rows = []
    for thr in cfg["candidates"]["min_angular_coverage_deg"]:
        for tower, g in baseline.groupby("tower"):
            rows.append(
                {
                    "tower": tower,
                    "min_angular_coverage_deg": thr,
                    "bands_total": len(g),
                    "bands_retained": int(
                        (g.angular_coverage_deg >= thr).sum()
                    ),
                    "retained_fraction": float(
                        (g.angular_coverage_deg >= thr).mean()
                    ),
                }
            )
    retention = pd.DataFrame(rows)

    # Simulated missing sectors on representative sections (surface
    # samples, see tests).
    sel = select_sections(baseline[baseline.representative].copy())
    base_seed = int(cfg["project"]["random_seed"])
    nrep = int(cfg["project"]["n_subsample_reps"])
    sim = []
    for _, s in sel.iterrows():
        tower = s["tower"]
        x, y, z, _ = clouds[tower]
        lo = float(s["lower_m"])
        idx = None
        for l2, h2, ii in band_slices(x, y, z, thick):
            if abs(l2 - lo) < 1e-9:
                idx = ii
                break
        if idx is None:
            continue
        xs, ys = x[idx], y[idx]
        ref = polygon_metrics(xs, ys)
        cx, cy = xs.mean(), ys.mean()
        ang = np.degrees(np.arctan2(xs - cx, ys - cy)) % 360.0
        for kept in cfg["candidates"]["min_angular_coverage_deg"]:
            gap = 360.0 - kept
            if gap <= 0:
                continue
            for rep in range(nrep):
                rng = np.random.default_rng(
                    base_seed + 7 * int(s["section_id"]) + 13 * kept + rep
                )
                start = rng.uniform(0, 360)
                rel = (ang - start) % 360.0
                keep = rel >= gap
                if keep.sum() < 3:
                    continue
                m = polygon_metrics(xs[keep], ys[keep])
                occ, cov = tm.sector_coverage(xs[keep], ys[keep])
                sim.append(
                    {
                        "tower": tower,
                        "section_id": int(s["section_id"]),
                        "section_lower_m": lo,
                        "kept_coverage_target_deg": kept,
                        "rep": rep,
                        "measured_coverage_deg": cov,
                        "n_kept": int(keep.sum()),
                        "area_error_fraction": (m["area"] - ref["area"])
                        / ref["area"],
                        "elong_error": m["elong"] - ref["elong"],
                        "az_error_axial_deg": abs(
                            ((m["az"] - ref["az"]) + 90) % 180 - 90
                        ),
                    }
                )
    sim_df = pd.DataFrame(sim)
    sim_df.to_csv(RESULTS / "angular_coverage_simulation.csv", index=False)
    if not sim_df.empty:
        agg = (
            sim_df.groupby("kept_coverage_target_deg")
            .agg(
                reps=("rep", "size"),
                median_abs_area_error=(
                    "area_error_fraction",
                    lambda s: float(np.median(np.abs(s))),
                ),
                p90_abs_area_error=(
                    "area_error_fraction",
                    lambda s: float(np.percentile(np.abs(s), 90)),
                ),
                median_abs_elong_error=(
                    "elong_error",
                    lambda s: float(np.median(np.abs(s))),
                ),
                median_az_error_deg=("az_error_axial_deg", "median"),
            )
            .reset_index()
        )
        agg.to_csv(RESULTS / "angular_coverage_error_summary.csv", index=False)
    retention.to_csv(RESULTS / "angular_coverage_results.csv", index=False)
    return retention


# ---------------------------------------------------------------- phase 7:
# elongation and azimuth


def phase_azimuth_reliability(
    baseline: pd.DataFrame, subsampling_summary: pd.DataFrame, cfg: dict
) -> pd.DataFrame:
    """Phase 7. Fraction of bands at or above each candidate elongation
    threshold, and
    the axial azimuth SD across subsamples grouped by reference elongation bin.
    Writes azimuth_reliability_results.csv and azimuth_sd_by_elongation.csv."""
    rows = []
    for thr in cfg["candidates"]["elongation_threshold"]:
        for tower, g in baseline.groupby("tower"):
            gb = g[np.isfinite(g.elongation_poly)]
            rows.append(
                {
                    "tower": tower,
                    "elongation_threshold": thr,
                    "bands_with_elongation": len(gb),
                    "bands_at_or_above_threshold": int(
                        (gb.elongation_poly >= thr).sum()
                    ),
                    "fraction_at_or_above": (
                        float((gb.elongation_poly >= thr).mean())
                        if len(gb)
                        else np.nan
                    ),
                }
            )
    frac = pd.DataFrame(rows)
    # Azimuth stability versus elongation, from the subsampling
    # replicates (all sizes).
    bins = [1.0, 1.05, 1.10, 1.20, 1.30, 1.50, 2.0, 1e9]
    labels = [
        "<1.05",
        "1.05-1.10",
        "1.10-1.20",
        "1.20-1.30",
        "1.30-1.50",
        "1.50-2.00",
        ">=2.00",
    ]
    if not subsampling_summary.empty:
        ss = subsampling_summary.copy()
        ss["ref_elong"] = ss["section_id"].map(
            baseline.reset_index()
            .drop_duplicates("tower")
            .set_index("tower")["elongation_poly"]
        )
        # reference elongation of each section, from the section table
        sel = select_sections(baseline[baseline.representative].copy())
        ref_map = sel.set_index("section_id")["elongation_poly"]
        ss["ref_elong"] = ss["section_id"].map(ref_map)
        ss["elong_bin"] = pd.cut(
            ss["ref_elong"], bins=bins, labels=labels, right=False
        )
        binned = (
            ss.groupby("elong_bin", observed=False)
            .agg(
                sections=("section_id", "nunique"),
                median_az_axial_sd_deg=("azimuth_axial_sd_deg", "median"),
                max_az_axial_sd_deg=("azimuth_axial_sd_deg", "max"),
            )
            .reset_index()
        )
        binned.to_csv(RESULTS / "azimuth_sd_by_elongation.csv", index=False)
    frac.to_csv(RESULTS / "azimuth_reliability_results.csv", index=False)
    return frac


# ---------------------------------------------------------------- phase 8:
# rock base persistence


def phase_rock_base(clouds: dict, cfg: dict) -> pd.DataFrame:
    """Phase 8. Lowest and highest excess-green detections for each band
    thickness,
    rock-share threshold and number of consecutive bands, compared with
    each tower's
    visual base reference. The colour rule is a heuristic, so results are
    reported as
    persistence, not validation. Writes rock_base_sensitivity.csv."""
    rb = cfg["rock_base"]
    rows = []
    for thick in cfg["candidates"]["band_thickness_m"]:
        for tower, (x, y, z, rgb) in clouds.items():
            slices = band_slices(x, y, z, float(thick))
            share = rock_share_per_band(
                rgb[:, 0], rgb[:, 1], rgb[:, 2], slices
            )
            lows = np.array([lo for lo, _, _ in slices])
            ref = cfg["towers"][tower]
            for thr_share in rb["rock_share_thresholds"]:
                for k in cfg["candidates"]["consecutive_bands"]:
                    ok = np.array([s >= thr_share for s in share])
                    lowest = highest = np.nan
                    n = len(ok)
                    # Lowest band starting k consecutive rock-like bands.
                    for i in range(n - k + 1):
                        stop = i + k
                        if ok[i:stop].all():
                            lowest = float(lows[i])
                            break
                    # Highest band ending k consecutive rock-like bands.
                    for j in range(n - 1, k - 2, -1):
                        start = j - k + 1
                        stop = j + 1
                        if ok[start:stop].all():
                            highest = float(lows[j] + thick)
                            break
                    rows.append(
                        {
                            "tower": tower,
                            "band_thickness_m": thick,
                            "rock_share_threshold": thr_share,
                            "consecutive_bands": k,
                            "lowest_detection_m": lowest,
                            "highest_detection_m": highest,
                            "detected": bool(np.isfinite(lowest)),
                            "reference_m": ref["base_reference_m"],
                            "reference_type": ref["base_reference_type"],
                            "reference_status": ref["base_reference_status"],
                            "lowest_minus_reference_m": (
                                (lowest - ref["base_reference_m"])
                                if np.isfinite(lowest)
                                else np.nan
                            ),
                        }
                    )
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "rock_base_sensitivity.csv", index=False)
    return out


# ---------------------------------------------------------------- phase
# 9-10: cross-tower and LOTO


def _criteria_met(sub: pd.DataFrame, crit: dict) -> bool:
    """True if a subsampling summary meets the stability criteria for one
    tower."""
    if sub.empty:
        return False
    # Median rule (chosen 9 Nov 2026): the median over sections must meet
    # each limit.
    bias_ok = (
        sub.area_bias_fraction.abs().median() <= crit["max_abs_bias_fraction"]
    )
    cv_ok = sub.area_cv_fraction.median() <= crit["max_cv_fraction"]
    az_ok = (
        sub.azimuth_axial_sd_deg.median() <= crit["max_azimuth_axial_sd_deg"]
    )
    return bool(bias_ok and cv_ok and az_ok)


def phase_cross_tower(baseline, point_ret, sub_summary, cfg):
    """Phase 9. For each candidate point count and tower, retention and
    stability
    summaries, and whether that tower meets the stability criteria under
    the median
    rule. Writes cross_tower_comparison.csv."""
    crit = cfg["stability_criteria"]
    rows = []
    for size in cfg["candidates"]["min_points_per_band"]:
        for tower in sorted(baseline.tower.unique()):
            pr = point_ret[
                (point_ret.tower == tower) & (point_ret.min_points == size)
            ]
            sub = (
                sub_summary[
                    (sub_summary.tower == tower)
                    & (sub_summary.subsample_n == size)
                ]
                if not sub_summary.empty
                else pd.DataFrame()
            )
            rows.append(
                {
                    "parameter": "min_points_per_band",
                    "value": size,
                    "tower": tower,
                    "retention_fraction": (
                        float(pr.retained_fraction_all.iloc[0])
                        if len(pr)
                        else np.nan
                    ),
                    "median_abs_area_bias": (
                        float(sub.area_bias_fraction.abs().median())
                        if len(sub)
                        else np.nan
                    ),
                    "median_area_cv": (
                        float(sub.area_cv_fraction.median())
                        if len(sub)
                        else np.nan
                    ),
                    "median_az_sd_deg": (
                        float(sub.azimuth_axial_sd_deg.median())
                        if len(sub)
                        else np.nan
                    ),
                    "criteria_met": _criteria_met(sub, crit),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "cross_tower_comparison.csv", index=False)
    return out


def phase_leave_one_tower_out(cross: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Phase 10. Leave-one-tower-out check: four folds, one per tower. For
    each held-out
    tower, select a value using the other three towers, then test that
    value on the
    held-out tower. Selection rule: among values that meet the criteria on
    at least
    two training towers, take the one with the highest mean retention. Writes
    leave_one_tower_out_results.csv. With only four towers, each fold
    leaves three
    towers to select from, so results are sensitive to which tower is held out.
    """
    rows = []
    towers = sorted(cross.tower.unique())
    for held in towers:
        train = cross[cross.tower != held]
        test = cross[cross.tower == held]
        for param in sorted(cross.parameter.unique()):
            tp = train[train.parameter == param]
            vals = sorted(tp.value.unique())
            best, best_ret = None, -1.0
            for v in vals:
                tv = tp[tp.value == v]
                if tv.criteria_met.sum() >= 2:
                    ret = float(tv.retention_fraction.mean())
                    if ret > best_ret:
                        best, best_ret = v, ret
            if best is None:
                rows.append(
                    {
                        "held_out_tower": held,
                        "parameter": param,
                        "selected_value": np.nan,
                        "held_out_criteria_met": np.nan,
                        "note": "no value met criteria on >=2 training towers",
                    }
                )
                continue
            tvt = test[test.value == best]
            rows.append(
                {
                    "held_out_tower": held,
                    "parameter": param,
                    "selected_value": best,
                    "held_out_criteria_met": (
                        bool(tvt.criteria_met.iloc[0]) if len(tvt) else np.nan
                    ),
                    "held_out_retention": (
                        float(tvt.retention_fraction.iloc[0])
                        if len(tvt)
                        else np.nan
                    ),
                    "note": "",
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "leave_one_tower_out_results.csv", index=False)
    return out
