# Sensitivity analysis of cross-sectional morphometry: Zhangjiajie sandstone towers

Status: analysis complete for points per band and rock base. Cross-tower and leave-one-tower-out checks for the other parameters are not done (section 6). Four towers are used (J_T01, J_T02, Q_T01, Y_T01); M_T01 is excluded by user decision (6 Nov 2026).

## 1. Purpose and scope

This analysis tests how sensitive the vertical cross-section measurements are to analysis choices: points per band, band thickness, angular coverage, elongation threshold, consecutive-band persistence, and the rock-base reference.

**Stability is not accuracy.** Results show whether measurements change when sampling changes. They do not show that the outlines match the true rock geometry. Accuracy has been checked against one visual reference only (J_T01).

## 2. Data and definitions

- Point clouds: five LAS files, EPSG:32649 (UTM 49N), projected metres. Q_T01 is a 20% seeded subsample for all metrics.
- Bands: horizontal slices by actual Z elevation, default 2 m; sensitivity at 1 m and 5 m.
- Envelope area: convex hull of the sampled points in a band. A sampled-surface envelope, not solid rock area. Not for volume or mass.
- Elongation: sqrt(lambda1 / lambda2) from polygon second moments of the hull (uniform density).
- Azimuth: axial direction of the major axis, 0 to 180 degrees clockwise from north. Polygon moments are computed with the origin shifted to the vertex mean to avoid UTM-scale numerical error.
- Angular coverage: number of occupied 10-degree sectors around the centroid. Valid for surface samples only.
- Representative section: at least 1,000 points, at least 180 degrees of angular coverage, a single object component (approved 6 Nov 2026), and, for J_T02 only, a lower edge at or above the visual base of 888 m (user decision 9 Nov 2026: bands below the base treated as canopy or skirt).

## 3. Stability criteria and the rule used

| Criterion | Limit |
|---|---|
| Bias: \|subsample area - full area\| / full area | median over sections < 5% |
| Variability: CV across 100 seeded subsamples | median < 10% |
| Axial azimuth SD across subsamples | median < 15 degrees |
| Band retention at the chosen setting | >= 80% |
| Cross-tower consistency | Criteria met on at least 3 of 4 towers |

**Rule.** A tower passes a setting if the median over its sections meets each limit (the median rule, chosen 9 Nov 2026). Stage A did not specify whether the limit applies to the median or to a majority of sections. An earlier version of the code used a majority-of-sections rule; it was changed to the median rule before the results below were produced. Under the majority rule, 500 points would have passed Q_T01. Under the median rule it does not.

## 4. Recommended setting: 1,000 points per band

Selected by the user on 9 Nov 2026.

| Tower | Retention | Median abs. area bias | Median area CV | Median axial azimuth SD |
|---|---|---|---|---|
| J_T01 | 0.96 | 0.015 | 0.004 | 0.21 deg |
| J_T02 | 0.97 | 0.013 | 0.004 | 0.53 deg |
| Q_T01 | 0.99 | 0.034 | 0.009 | 1.17 deg |
| Y_T01 | 0.98 | 0.024 | 0.010 | 0.36 deg |

All four towers pass all three stability limits at 1,000 points.

Median absolute area bias at other point counts (median rule):

| Points per band | J_T01 | J_T02 | Q_T01 | Y_T01 | Towers passing |
|---|---|---|---|---|---|
| 50 | 0.100 | 0.089 | 0.174 | 0.134 | 0 |
| 100 | 0.071 | 0.062 | 0.126 | 0.101 | 0 |
| 200 | 0.048 | 0.041 | 0.094 | 0.066 | 2 (J_T01, J_T02) |
| 500 | 0.027 | 0.023 | 0.057 | 0.038 | 3 (Q_T01 fails) |
| 1000 | 0.015 | 0.013 | 0.034 | 0.024 | 4 |

The 1,000-point setting equals the representative-section threshold, so it does not change which sections are reported.

## 5. Other parameters

**Band thickness.** Representative bands (at least 1,000 points) at 1 / 2 / 5 m: J_T01 28 / 23 / 15; J_T02 24 / 17 / 13; Q_T01 0 / 2 / 6; Y_T01 15 / 9 / 9. Median envelope area at 2 m: J_T01 493 m2, J_T02 601 m2, Q_T01 877 m2 (two sections), Y_T01 424 m2. Median envelope area increases with thickness for J_T01 (464 / 493 / 601 m2) and J_T02 (592 / 601 / 631 m2). Median elongation falls with thickness for J_T01 and Y_T01 and is roughly flat for J_T02 (1.49 to 1.53 across the three thicknesses). Q_T01 has no representative section at 1 m. Band thickness is reported as a sensitivity variable; it has no stability criterion and was not selected. The 2 m default is kept.

**Angular coverage.** On simulated missing-sector surface samples, median absolute area error is 0.86 (90 degrees), 0.74 (120), 0.51 (180), 0.23 (240) and 0.11 (270). Median axial azimuth error is 23, 27, 24, 13 and 6 degrees. The approved 180-degree threshold admits sections whose area could be about half wrong. Section-level shape claims should be made with the 240 or 270 degree figures in mind.

**Elongation threshold.** Azimuth SD by elongation bin: 1.10-1.20 median 7.6 degrees (maximum 16.9; one section); 1.20-1.30 median 1.8 (maximum 7.9; four sections); 1.30-1.50 median 1.5 (maximum 3.9); 1.50-2.00 median 0.7; 2.00 and above median 0.55 (maximum 2.4). The provisional 1.20 is the lowest threshold where the azimuth spread stays within 15 degrees in every bin tested. It is a lower bound, not an optimum.

**Consecutive bands (rock-base persistence).** At 2 m bands and 50% rock share, the lowest colour-based detection is: J_T01 806 m for all 1 to 5 bands (visual base 807.65 m); J_T02 888 m for 1 to 3 bands and 902 m for 4 to 5 (visual base 888 m); Q_T01 758 m for 1 to 4 bands and 774 m at 5 (visual base 736.96 m); Y_T01 936 m for all (visual base 938.74 m). Persistence changes the detected lower limit only for J_T02 and Q_T01. The provisional value of 3 is not validated.

## 6. Cross-tower and leave-one-tower-out checks

- **Cross-tower (points per band):** the table in section 4 and `results/cross_tower_comparison.csv`.
- **Leave-one-tower-out (points per band):** on the other three towers, select the value that passes on at least two of them, with the highest retention; then test it on the held-out tower. Results under the median rule:

| Held-out tower | Value selected from the other three | Held-out criteria met |
|---|---|---|
| J_T01 | 500 points | Yes |
| J_T02 | 500 points | Yes |
| Q_T01 | 200 points | No |
| Y_T01 | 200 points | No |

Two of four folds fail. The selection rule chooses 200 points when both 200 and 500 pass on the training towers and 200 has the marginally higher retention. The check therefore shows that the rule cannot be relied on to choose a passing value for every tower, which supports fixing the point count by inspection of the full stability results rather than by an automated rule. It does not validate 1,000 points, and it should be reported as a limitation.

- **Not done:** cross-tower and leave-one-tower-out checks for band thickness, angular coverage, elongation threshold and consecutive bands. These parameters have no stability criteria yet, so leave-one-tower-out cannot be applied to them as it stands. A descriptive cross-tower comparison can be built from existing outputs if needed.

## 7. Rock-base reference

- Reference values used: J_T01 807.65 m (visual clip floor, usable). J_T02 888 m, Q_T01 736.96 m and Y_T01 938.74 m (approximate side-view estimates from one observer).
- `results/rock_base_sensitivity.csv` uses these visual references. The original cloud-minimum values are kept in the column `reference_cloud_min_m`.
- The J_T02 visual estimate (888 m) matches the lowest detection at 1 to 3 consecutive bands.
- Q_T01 colour-based detections begin about 21 m above its visual estimate. The colour rule does not reproduce the Q_T01 base, so the Q_T01 reference is uncertain.
- The colour rule (excess green) is a heuristic. Detections are reported as stability, not validation, except for J_T01.
- The cloud minimum is never labelled a rock base in any output.

## 8. J_T02 lowest band

The lowest representative J_T02 band at 2 m (882 to 884 m; 3,320 points; elongation 6.6) lay below the 888 m visual base. It is excluded from representative sections as canopy or skirt (user decision, 9 Nov 2026). The 2 m representative count for J_T02 falls from 18 to 17. Median envelope area and elongation for J_T02 are reported after exclusion.

## 9. Limitations

1. Stability is not accuracy. Only J_T01 has a usable independent reference.
2. Q_T01 is the weakest tower: a 20% subsample, two representative sections at 2 m, none at 1 m.
3. The 180-degree coverage threshold is generous for surface samples (section 5).
4. Leave-one-tower-out covers point count only. Under the median rule it fails in two of four folds, so the automated selection rule should not be used to set the threshold.
5. The 1,000-point threshold already exceeds all lower point-count candidates, so point-count sensitivity mainly governs non-representative bands.
6. Visual rock-base estimates are approximate, from one observer.
7. Vegetation removal is manual clipping. The excess-green rule is used only for rock-base detection, and no all-points versus rock-only comparison has been done for the morphometry (analysis plan, Phase 2).
8. Convex-hull envelopes include canopy where it remains.

## 10. Suggested wording for the paper

"At 1,000 points per 2 m band, median absolute area bias was below 3.5% and median axial azimuth SD below 1.2 degrees across four towers, with 96 to 100% of bands retained. At 500 points, the held-out Qiankunzhu tower failed the stability criteria. An automated leave-one-tower-out selection failed in two of four folds, so the threshold was fixed on the full results rather than by automated selection. These results indicate that the metrics do not depend strongly on the random selection or density of sampled points. They do not validate the accuracy of the reconstructed outlines, which has been checked against one independent visual reference only."

## 11. Files

- `results/`: all CSV outputs, including `parameter_recommendations.csv`.
- `figures/`: PNG (300 dpi) and PDF versions, generated by `make_figures.py`.
- `logs/run_log.jsonl`: timestamped run log.
- `README.md`: commands and environment.
