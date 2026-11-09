# Sensitivity analysis: cross-sectional morphometry of Zhangjiajie towers

This folder contains the Stage B sensitivity analysis. Read `sensitivity_analysis_report.md` for the methods, results and limitations, and `results/parameter_recommendations.csv` for the recommended settings.

## Environment

Run all commands from the project root, `C:\tmp\China\Zhangjiajie_Geomorphology`, using the project virtual environment:

- Python environment: `.venv` (created from the QGIS Python; see `requirements.txt`).
- Required packages: numpy, pandas, scipy, laspy, matplotlib, pyyaml, pytest.

## Commands

Run the baseline first, then the phases. Each command writes to `results/` and `logs/run_log.jsonl`.

```bash
.venv\Scripts\python.exe sensitivity_analysis\sensitivity_analysis.py baseline
.venv\Scripts\python.exe sensitivity_analysis\run_sensitivity.py all
.venv\Scripts\python.exe sensitivity_analysis\make_figures.py
```

To run one phase at a time:

```bash
.venv\Scripts\python.exe sensitivity_analysis\run_sensitivity.py thickness
.venv\Scripts\python.exe sensitivity_analysis\run_sensitivity.py pointcount
.venv\Scripts\python.exe sensitivity_analysis\run_sensitivity.py coverage
.venv\Scripts\python.exe sensitivity_analysis\run_sensitivity.py azimuth
.venv\Scripts\python.exe sensitivity_analysis\run_sensitivity.py rockbase
.venv\Scripts\python.exe sensitivity_analysis\run_sensitivity.py crosstower
```

Run the automated geometry tests before any new analysis:

```bash
.venv\Scripts\python.exe -m pytest sensitivity_analysis\tests
```

## Configuration

All candidate values, stability limits, the representative-section rule, and the rock-base references are in `sensitivity_config.yaml`. Every output records the values it used. Change the configuration there, not in the scripts.

Random seed is fixed at 42. Subsampling uses 100 seeded replicates per setting.

## Files

| File or folder | Purpose |
|---|---|
| `stage_a_audit.md` | Stage A: dataset audit, design, stability criteria and approved decisions |
| `sensitivity_config.yaml` | Candidate values, criteria, references |
| `tower_morphometry.py` | Band slicing, hull area, polygon-moment azimuth and elongation, bounding rectangle, angular coverage, component checks |
| `sensitivity_analysis.py` | Config loading, cloud loading (with cache), inventory, baseline |
| `sensitivity_phases.py` | The sensitivity phases and the cross-tower and leave-one-tower-out logic |
| `run_sensitivity.py` | Runs the phases in order |
| `make_figures.py` | Figures from the results CSVs |
| `tests/test_morphometry.py` | Automated geometry tests on synthetic shapes |
| `debug_polygon_moments.py`, `debug_components.py` | Diagnostics used to find and check the polygon-moment and component issues. Read-only. Kept for reproducibility |
| `cache/` | Cached point-cloud arrays (`.npz`). Can be regenerated from the source files; safe to delete |
| `results/` | All CSV outputs, including `parameter_recommendations.csv` |
| `figures/` | Figures as PNG (300 dpi) and PDF |
| `logs/run_log.jsonl` | Timestamped run log |

## Known issues

- The leave-one-tower-out check (median rule) fails in two of four folds (Q_T01 and Y_T01 held out, value selected 200). See report section 6.
- `results/rock_base_sensitivity.csv` uses the visual references. The original cloud-minimum values are in `reference_cloud_min_m`.
- The lowest representative band for J_T02 (882 m) lies below the visual base (888 m). It is kept pending a decision (report section 8).
- The cross-tower and leave-one-tower-out checks cover the point-count parameter only.
