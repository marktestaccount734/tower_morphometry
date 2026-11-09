# Stage A: dataset audit and design (cross-sectional morphometry and sensitivity)

Status: complete. Stage B (implementation and analysis) waits for approval.

## 1. Existing code

- `analyse_tower_sections.py` (earlier prototype): not present anywhere in the project, Desktop or Documents. Nothing to reuse; a new module is required.
- `zhangjiajie_python_package/zhangjiajie_facets_analysis.py`: facet statistics only. No cross-section, segmentation, vegetation filter or QC code.
- CloudCompare outputs (`*_metrics.txt`, `facets.csv`): usable inputs, not analysis code.

## 2. Dataset inventory (all five, source files under `data/segmented/`)

| Tower | Points (file) | Format | CRS | Bbox X x Y (m) | Z range (m) | Density (pts/m², bbox) | Attributes | Notes |
|---|---|---|---|---|---|---|---|---|
| J_T01 | 535,224 | LAS 1.2, fmt 2 | EPSG:32649 | 44.8 x 65.1 | 807.65–949.22 | ~140 | RGB (16-bit), confidence, normal x/y/z (stored) | Re-saved file; provenance noted in methods |
| J_T02 | 438,885 | LAS 1.2, fmt 2 | EPSG:32649 | 32.0 x 73.8 | 881.35–995.42 | ~186 | RGB, confidence, normals (stored) | Canopy skirt included |
| M_T01 | 1,638,144 | LAS 1.2, fmt 2 | EPSG:32649 | 127.5 x 92.1 | 581.20–857.48 | ~140 | RGB, confidence, normals (stored) | Footprint suggests more than one pillar; to check |
| Q_T01 | 3,435,314 (full); 687,064 (subsample) | LAS 1.2, fmt 2 | EPSG:32649 | 56.5 x 72.9 | 735.56–937.82 | ~835 full; ~167 subsample | RGB, confidence, normals (stored) | Use subsample for all metrics (agreed) |
| Y_T01 | 499,978 | LAS 1.2, fmt 2 | EPSG:32649 | 39.4 x 27.7 | 937.33–1,025.47 | ~458 | RGB, confidence, normals (stored) | Densest of the full-file towers after Q_T01 |

Common checks:
- Coordinates are projected metres (UTM 49N, GeoKey 32649 and WKT both confirmed). No degree-based coordinates.
- Vertical datum: not recorded in any file. Treated as unknown.
- Colour: 16-bit RGB, maximum about 65,280, i.e. 8-bit values scaled by 256. Divide by 256 before any colour rule.
- Classification: all points are class 0 (unclassified). No segmentation or vegetation labels exist.
- Confidence: Metashape values, integer-valued, range differs by tower (J_T01 1–153, J_T02 1–95, M_T01 1–108, Q_T01 1–134, Y_T01 1–36). Retain as supplementary QC only, as the spec says.
- Stored normal x/y/z fields are source-file fields, not the computed normals. Excluded from analysis.
- Duplicates: about 4–9% of points share exact coordinates. Retained, per the agreed rule.

## 3. Geometric validity

- **Individual tower:** J_T01, J_T02, Q_T01 and Y_T01 are single towers on the evidence of their footprints and earlier visual checks. M_T01's footprint (127.5 x 92.1 m) suggests more than one pillar; it needs a visual check before its sections are used.
- **Vegetation:** manual clipping only, as far as practicable. Residual canopy remains, especially on the lower bands (J_T02, Y_T01 and Q_T01 in particular). Convex hulls of these bands describe canopy as well as rock.
- **Convex hulls:** a hull is a geometric envelope of sampled points. It bridges concavities, gaps and any canopy or neighbouring material. It is not a solid rock cross-section, and it must not be used for volume or mass. Hull area is therefore a diagnostic of the sampled envelope, and it is most suspect where canopy dominates or where the outline is concave.
- **Completeness:** the surfaces visible in the scans are one-sided for some towers. Angular coverage must be measured, not assumed.

## 4. Proposed design (Stage B)

Code (new, in `sensitivity_analysis/`):

- `tower_morphometry.py`: read LAS with laspy, slice on actual Z, convex hull (scipy), hull area and perimeter, polygon second-moment orientation with 180° axial periodicity, minimum-area bounding rectangle as the span and elongation check, sector-based angular coverage (36 sectors of 10°), QC flags. Records failures explicitly.
- `sensitivity_analysis.py`: the experiments in Phases 4–9.
- `sensitivity_config.yaml`: the candidate values from the spec, the seeds and the stability criteria.

Verification (automated tests before any real data):

- circle, ellipse, rectangle, rotated ellipse with known axes and azimuth, two touching cylinders, empty and degenerate bands.
- Azimuth tests must check 180° periodicity.

Measurement definitions to fix before the analysis:

- **Major/minor spans:** from the oriented minimum bounding rectangle. PCA on hull vertices is not used.
- **Azimuth:** from the second moments of the hull polygon, reported as an axial angle in [0, 180).
- **Centroid:** the area centroid of the hull polygon, with the point centroid kept as a separate column.
- **Angular coverage:** the number of occupied 10° sectors around the centroid, in degrees.

## 5. Stability criteria (to be fixed before selection)

- **Bias:** |mean subsample area − full-resolution area| / full-resolution area < 5%.
- **Variability:** coefficient of variation across 100 seeded subsamples < 10%.
- **Convergence:** the bias and variability above are met at or below the chosen point count.
- **Azimuth reliability:** axial SD across subsamples < 15°, and elongation above the chosen threshold.
- **Retention:** at least 80% of bands retained at the chosen settings, reported per tower.
- **Cross-tower consistency:** the same setting meets the criteria for at least four of five towers.

## 6. Computational estimate

- Slicing and hulls: a few seconds per tower for full-resolution bands at 2 m, and a few minutes for 1 m bands on Q_T01 and M_T01.
- Subsampling: 5 point counts x 100 seeds x about 30 sections = 15,000 hull evaluations, each on at most a few thousand points. A few minutes in total.
- Memory: the largest cloud (Q_T01 full, 3.4 M points) fits in memory at 16 bytes per point for the XYZ columns.

## 7. Decisions needed before Stage B

1. **Representative sections:** the spec asks for well-reconstructed sections. Proposed criterion: at least 1,000 points, at least 180° of angular coverage, and no flagged neighbour contamination. Approve or change.
2. **Rock-base reference (Phase 8):** the spec asks for comparison with independently inspected rock–vegetation boundaries. We have only the fallback cloud-minimum bases, plus the J_T01 clip floor. We need side-view visual base heights for at least J_T02, M_T01 and Q_T01 before persistence sensitivity can be judged.
3. **M_T01:** confirm whether it is one tower or more. If more than one, it needs splitting before its sections are used.
4. **Q_T01:** all Q_T01 sensitivity runs use the subsample, as agreed. The full-density file is not used.
5. **Output location:** `sensitivity_analysis/` at the project root, as the spec says, with results CSVs there.
6. **Candidate values:** adopt the spec's candidate table as written, and replace the provisional values in `config/default.yaml` for this analysis.

## 8. Decisions recorded (user, 6 Nov 2026)

1. Representative sections: approved (at least 1,000 points, at least 180° coverage, no neighbour contamination).
2. Rock-base reference (Phase 8): OPEN. Side-view base heights for J_T02 and Q_T01 (and Y_T01 if wanted) are needed, or Phase 8 is reported as uncertain against the cloud-minimum fallback bases.
3. M_T01: EXCLUDED from all analyses. Four towers remain (J_T01, J_T02, Q_T01, Y_T01). Leave-one-tower-out is therefore four folds.
4. Output location (`sensitivity_analysis/` at the project root) and the spec's candidate values: approved.
5. Q_T01: subsample only, as already agreed.

## 9. Phase 8 reference table (agreed 6 Nov 2026)

| Tower | Reference elevation (m) | Reference type | Label in outputs | Rock-base status |
|---|---|---|---|---|
| J_T01 | 807.65 | Clip floor set visually at lowest visible rock above canopy skirt | Visual rock-base judgement | Usable reference |
| J_T02 | 881.35 | Cloud minimum (canopy-skirt minimum) | Lower extent of reconstruction, not rock base | Uncertain |
| Q_T01 | 735.56 | Cloud minimum (canopy-skirt minimum) | Lower extent of reconstruction, not rock base | Uncertain |
| Y_T01 | 937.33 | Cloud minimum (canopy-skirt minimum) | Lower extent of reconstruction, not rock base | Uncertain |

M_T01 is excluded from all analyses.

Interpretation rules:

- Persistence results are reported as stability of detections, not validation of the rock base, except for J_T01.
- The cloud minimum is never labelled as a rock base in any output, figure or table.
- Detections below a cloud-minimum reference are flagged as outside the reconstruction.
