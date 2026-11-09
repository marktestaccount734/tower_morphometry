"""Figures for the sensitivity analysis (300 dpi PNG and vector PDF). Reads
results CSVs only.

Run from the project root:
    .venv\\Scripts\\python.exe sensitivity_analysis\\make_figures.py

Author: Anonymous (blinded for peer review)
Project: 3D Geospatial Approaches for Geoheritage Documentation and Analysis
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
R = HERE / "results"
F = HERE / "figures"
F.mkdir(exist_ok=True)
TOWERS = ["J_T01", "J_T02", "Q_T01", "Y_T01"]
COL = {
    "J_T01": "#1b9e77",
    "J_T02": "#d95f02",
    "Q_T01": "#7570b3",
    "Y_T01": "#e7298a",
}
plt.rcParams.update(
    {"font.size": 9, "axes.spines.top": False, "axes.spines.right": False}
)


def save(fig, name):
    fig.savefig(F / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(F / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_area_elevation():
    """Convex-hull envelope area against band midpoint elevation, one panel
    per tower,
    coloured by band thickness (1, 2 and 5 m). Representative bands only."""
    d = pd.read_csv(R / "band_thickness_results.csv")
    d = d[d.representative]
    fig, axes = plt.subplots(1, 4, figsize=(10, 3.4), sharey=False)
    cmap = plt.get_cmap("viridis")
    thick_colour = {
        1.0: cmap(0.15),
        2.0: cmap(0.5),
        5.0: cmap(0.85),
    }  # viridis ramp, one colour per thickness
    for ax, t in zip(axes, TOWERS):
        g = d[d.tower == t]
        for thick, style in [(1.0, "o"), (2.0, "s"), (5.0, "^")]:
            h = g[g.band_thickness_m == thick]
            ax.scatter(
                h.lower_m + thick / 2,
                h.hull_area_m2,
                s=14,
                marker=style,
                color=thick_colour[thick],
                label=f"{thick:g} m",
                alpha=0.8,
            )
        ax.set_title(t)
    axes[0].set_ylabel("Convex-hull envelope area (m²)")
    axes[0].legend(title="Band", frameon=False)
    fig.subplots_adjust(wspace=0.35, bottom=0.22)
    fig.supxlabel(
        "Band midpoint elevation (m)", y=0.02
    )  # one shared x-label for all four panels
    save(fig, "area_vs_elevation_by_band_thickness")


def fig_elongation_elevation():
    """Polygon elongation against band midpoint elevation, one panel per
    tower, by band
    thickness, with the 1.2 provisional threshold marked. Representative
    bands only.
    """
    d = pd.read_csv(R / "band_thickness_results.csv")
    d = d[d.representative & np.isfinite(d.elongation_poly)]
    fig, axes = plt.subplots(1, 4, figsize=(10, 3.4), sharey=True)
    cmap = plt.get_cmap("viridis")
    thick_colour = {
        1.0: cmap(0.15),
        2.0: cmap(0.5),
        5.0: cmap(0.85),
    }  # viridis ramp, one colour per thickness
    for ax, t in zip(axes, TOWERS):
        g = d[d.tower == t]
        for thick, style in [(1.0, "o"), (2.0, "s"), (5.0, "^")]:
            h = g[g.band_thickness_m == thick]
            ax.scatter(
                h.lower_m + thick / 2,
                h.elongation_poly,
                s=14,
                marker=style,
                color=thick_colour[thick],
                label=f"{thick:g} m",
                alpha=0.8,
            )
        ax.axhline(1.2, color="grey", lw=0.8, ls="--")
        ax.set_title(t)
    axes[0].set_ylabel("Polygon elongation sqrt(λ1/λ2)")
    axes[0].legend(title="Band", frameon=False)
    fig.subplots_adjust(wspace=0.2, bottom=0.22)
    fig.supxlabel(
        "Band midpoint elevation (m)", y=0.02
    )  # one shared x-label for all four panels
    save(fig, "elongation_vs_elevation_by_band_thickness")


def fig_stability_vs_size():
    """Median area bias, area CV and axial azimuth SD against requested
    points per band,
    with the stability limits marked. Medians are taken over representative
    sections.
    Levels are placed at equal spacing, not in proportion to point count."""
    s = pd.read_csv(R / "subsampling_results.csv")
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2))
    specs = [
        ("area_bias_fraction", "Absolute area bias (fraction)", True),
        ("area_cv_fraction", "Area coefficient of variation", False),
        ("azimuth_axial_sd_deg", "Axial azimuth SD (deg)", False),
    ]
    sizes = [50, 100, 200, 500, 1000]
    pos = list(
        range(len(sizes))
    )  # tested levels placed at equal spacing (categorical, not to scale)
    cmap = plt.get_cmap("viridis")
    tower_colour = {
        t: cmap(v) for t, v in zip(TOWERS, [0.1, 0.4, 0.7, 0.95])
    }  # viridis, one colour per tower
    for ax, (col, label, absval) in zip(axes, specs):
        for t in TOWERS:
            g = s[s.tower == t]
            if g.empty:
                continue
            vals = g[col].abs() if absval else g[col]
            agg = vals.groupby(g.subsample_n).median().reindex(sizes)
            ax.plot(
                pos, agg.values, "o-", color=tower_colour[t], label=t, ms=4
            )
        ax.set_xticks(pos)
        ax.set_xticklabels([f"{n:,}" for n in sizes])
        ax.set_xlim(-0.3, len(sizes) - 0.7)
        ax.set_ylabel(label)
    axes[0].axhline(0.05, color="grey", ls="--", lw=0.8)
    axes[1].axhline(0.10, color="grey", ls="--", lw=0.8)
    axes[2].axhline(15, color="grey", ls="--", lw=0.8)
    axes[0].legend(frameon=False)
    fig.subplots_adjust(
        wspace=0.45, bottom=0.22
    )  # gap between panels; room below for the shared label
    fig.supxlabel(
        "Points per band", y=0.02
    )  # one shared x-label for all three panels
    save(fig, "stability_vs_point_sample_size")


def fig_coverage_error():
    """Area error, elongation error and azimuth error against kept angular
    coverage,
    from simulated missing sectors on surface samples."""
    e = pd.read_csv(R / "angular_coverage_error_summary.csv")
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2))
    cmap = plt.get_cmap("viridis")
    c_median, c_p90 = cmap(0.15), cmap(
        0.7
    )  # viridis: median and 90th percentile
    axes[0].plot(
        e.kept_coverage_target_deg,
        e.median_abs_area_error,
        "o-",
        color=c_median,
        label="median",
    )
    axes[0].plot(
        e.kept_coverage_target_deg,
        e.p90_abs_area_error,
        "s--",
        color=c_p90,
        label="90th percentile",
    )
    axes[0].set_ylabel("Absolute area error (fraction)")
    axes[1].plot(
        e.kept_coverage_target_deg,
        e.median_abs_elong_error,
        "o-",
        color=c_median,
    )
    axes[1].set_ylabel("Absolute elongation error")
    axes[2].plot(
        e.kept_coverage_target_deg, e.median_az_error_deg, "o-", color=c_median
    )
    axes[2].set_ylabel("Median axial azimuth error (deg)")
    for ax in axes:
        ax.set_xticks(e.kept_coverage_target_deg)
    axes[0].legend(frameon=False)
    fig.subplots_adjust(bottom=0.22)
    fig.supxlabel(
        "Kept angular coverage (deg)", y=0.02
    )  # one shared x-label for all three panels
    save(fig, "angular_coverage_error")


def fig_azimuth_vs_elongation():
    """Axial azimuth SD across subsamples against reference section
    elongation, with the
    15-degree limit and the 1.2 elongation threshold marked."""
    s = pd.read_csv(R / "subsampling_results.csv")
    sel = pd.read_csv(R / "baseline_measurements.csv")
    sel = sel[sel.representative][["tower", "lower_m", "elongation_poly"]]
    # section reference elongation via the subsampling table (section_lower_m)
    s = s.merge(
        sel,
        left_on=["tower", "section_lower_m"],
        right_on=["tower", "lower_m"],
        how="left",
    )
    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    cmap = plt.get_cmap("viridis")
    tower_colour = {
        t: cmap(v) for t, v in zip(TOWERS, [0.1, 0.4, 0.7, 0.95])
    }  # same viridis mapping as the other tower figures
    for t in TOWERS:
        g = s[(s.tower == t) & s.elongation_poly.notna()]
        if g.empty:
            continue
        ax.scatter(
            g.elongation_poly,
            g.azimuth_axial_sd_deg,
            s=18,
            color=tower_colour[t],
            label=t,
            alpha=0.8,
        )
    ax.axhline(15, color="grey", ls="--", lw=0.8)
    ax.axvline(1.2, color="grey", ls=":", lw=0.8)
    ax.set_xlabel("Section elongation (polygon moments)")
    ax.set_ylabel("Axial azimuth SD across subsamples (deg)")
    ax.legend(frameon=False, fontsize=8)
    save(fig, "azimuth_sd_vs_elongation")


def fig_rock_base():
    """Lowest and highest colour-based rock-base detections against
    consecutive bands
    required (2 m bands, 50% rock share), with each tower's visual
    reference as a
    dotted line. The colour rule is a heuristic; the figure shows persistence.
    """
    r = pd.read_csv(R / "rock_base_sensitivity.csv")
    r = r[(r.band_thickness_m == 2.0) & (r.rock_share_threshold == 0.5)]
    fig, axes = plt.subplots(1, 4, figsize=(10, 3.4), sharey=False)
    cmap = plt.get_cmap("viridis")
    tower_colour = {
        t: cmap(v) for t, v in zip(TOWERS, [0.1, 0.4, 0.7, 0.95])
    }  # same viridis mapping as the stability figure
    for ax, t in zip(axes, TOWERS):
        g = r[r.tower == t].sort_values("consecutive_bands")
        ax.plot(
            g.consecutive_bands,
            g.lowest_detection_m,
            "o-",
            color=tower_colour[t],
            label="lowest",
        )
        ax.plot(
            g.consecutive_bands,
            g.highest_detection_m,
            "s--",
            color=tower_colour[t],
            alpha=0.6,
            label="highest",
        )
        ref = g.reference_m.iloc[0]
        ax.axhline(ref, color="black", lw=0.8, ls=":")
        status = g.reference_status.iloc[0]
        ax.set_title(t)
        ax.set_xticks([1, 2, 3, 4, 5])
    axes[0].set_ylabel("Elevation (m); dotted line = reference")
    axes[0].legend(frameon=False)
    fig.subplots_adjust(bottom=0.22)
    fig.supxlabel(
        "Consecutive bands required", y=0.02
    )  # one shared x-label for all four panels
    save(fig, "rock_base_vs_persistence")


def fig_heatmap():
    """Median absolute area bias for each tower and candidate point count,
    annotated
    with pass or fail under the median rule."""
    c = pd.read_csv(R / "cross_tower_comparison.csv")
    pivot_bias = c.pivot(
        index="tower", columns="value", values="median_abs_area_bias"
    )
    pivot_ok = c.pivot(
        index="tower", columns="value", values="criteria_met"
    ).astype(float)
    fig, ax = plt.subplots(figsize=(6.2, 2.8))
    im = ax.imshow(pivot_bias.values, cmap="viridis_r", aspect="auto")
    for i in range(pivot_bias.shape[0]):
        for j in range(pivot_bias.shape[1]):
            ok = pivot_ok.values[i, j] == 1.0
            ax.text(
                j,
                i,
                f"{pivot_bias.values[i, j]:.3f}\n{'pass' if ok else 'fail'}",
                ha="center",
                va="center",
                color="white" if pivot_bias.values[i, j] > 0.06 else "black",
                fontsize=7,
            )
    ax.set_xticks(range(pivot_bias.shape[1]))
    ax.set_xticklabels(pivot_bias.columns)
    ax.set_yticks(range(pivot_bias.shape[0]))
    ax.set_yticklabels(pivot_bias.index)
    ax.set_xlabel("Points per band (requested)")
    fig.colorbar(im, ax=ax, label="Median |area bias|")
    save(fig, "cross_tower_sensitivity_heatmap")


def fig_loto():
    """Leave-one-tower-out result: whether the value selected from
    the other three towers meets the criteria on the held-out
    tower."""
    lot = pd.read_csv(R / "leave_one_tower_out_results.csv")
    fig, ax = plt.subplots(figsize=(4.8, 2.8))
    vals = lot.held_out_criteria_met.astype(str)
    ax.bar(
        lot.held_out_tower,
        [1 if v == "True" else 0 for v in vals],
        color=[COL[t] for t in lot.held_out_tower],
    )
    ax.set_ylim(0, 1.2)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["not met", "met"])
    save(fig, "leave_one_tower_out")


if __name__ == "__main__":
    fig_area_elevation()
    fig_elongation_elevation()
    fig_stability_vs_size()
    fig_coverage_error()
    fig_azimuth_vs_elongation()
    fig_rock_base()
    fig_heatmap()
    fig_loto()
    print("figures written to", F)
