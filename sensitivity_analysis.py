"""Sensitivity analysis for tower cross-sectional morphometry (Stage B).

Run from the project root:
    python sensitivity_analysis\\sensitivity_analysis.py baseline
    python sensitivity_analysis\\sensitivity_analysis.py all

Every result table records the parameter values used. Inputs are read-only.

Author: Anonymous (blinded for peer review)
Project: 3D Geospatial Approaches for Geoheritage Documentation and Analysis
"""

from __future__ import annotations

import json
import platform
import sys
import time
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import tower_morphometry as tm  # noqa: E402

CFG_PATH = HERE / "sensitivity_config.yaml"
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
CACHE = HERE / "cache"
LOG = HERE / "logs" / "run_log.jsonl"


def load_config() -> dict:
    """Read sensitivity_config.yaml (candidate values, criteria and tower
    references)."""
    with open(CFG_PATH, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def log_event(event: str, **kw) -> None:
    """Append one timestamped JSON line to logs/run_log.jsonl."""
    LOG.parent.mkdir(parents=True, exist_ok=True)
    rec = {
        "time": datetime.now().isoformat(timespec="seconds"),
        "event": event,
        **kw,
    }
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, default=str) + "\n")


def load_cloud(tower: str, cfg: dict):
    """Load X, Y, Z and 8-bit RGB for a tower, caching to .npz (read-only
    on inputs)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    src = ROOT / cfg["towers"][tower]["cloud"]
    cache = CACHE / f"{tower}.npz"
    if cache.exists() and cache.stat().st_mtime >= src.stat().st_mtime:
        d = np.load(cache)
        return d["x"], d["y"], d["z"], d["rgb"]
    t0 = time.time()
    x, y, z, rgb = tm.read_tower(str(src))
    np.savez(cache, x=x, y=y, z=z, rgb=rgb)
    log_event(
        "cloud_loaded",
        tower=tower,
        source=str(src.relative_to(ROOT)),
        points=int(x.size),
        seconds=round(time.time() - t0, 2),
    )
    return x, y, z, rgb


def inventory(cfg: dict) -> pd.DataFrame:
    """Dataset inventory: source, CRS, extent, elevation range, density,
    attributes."""
    rows = []
    for tower in cfg["towers"]:
        x, y, z, rgb = load_cloud(tower, cfg)
        area = np.ptp(x) * np.ptp(y)
        rows.append(
            {
                "tower": tower,
                "source_file": cfg["towers"][tower]["cloud"],
                "points": int(x.size),
                "crs": cfg["project"]["crs"],
                "x_extent_m": round(float(np.ptp(x)), 2),
                "y_extent_m": round(float(np.ptp(y)), 2),
                "z_min_m": round(float(z.min()), 2),
                "z_max_m": round(float(z.max()), 2),
                "bbox_points_per_m2": round(float(x.size / area), 1),
                "rgb_bit_depth": "16 (stored); 8-bit after >> 8",
                "classification": "all unclassified (0)",
                "confidence_attribute": (
                    "Metashape confidence (supplementary QC)"
                ),
                "stored_normals": "normal x/y/z (source fields, excluded)",
                "vertical_datum": "unknown",
            }
        )
    return pd.DataFrame(rows)


def bands_to_frame(tower: str, results, prov: dict, cfg: dict) -> pd.DataFrame:
    """Convert BandResult objects to a table and mark each band as
    representative.

    A band is representative if it has at least the minimum points, the minimum
    angular coverage and (when required) a single object component, has a valid
    hull, and lies at or above the tower's exclude_below_base_m if one is set.
    The provisional parameter values are written into the table as param_*
    columns.
    """
    rows = []
    rep = cfg["representative_section"]
    for r in results:
        d = asdict(r)
        d.pop("flags")
        d["flags"] = ";".join(r.flags)
        d["tower"] = tower
        d["representative"] = bool(
            r.n_points >= rep["min_points"]
            and r.angular_coverage_deg >= rep["min_angular_coverage_deg"]
            and (r.n_components == 1 or not rep["require_single_component"])
            and not np.isnan(r.hull_area_m2)
            and (
                cfg["towers"][tower].get("exclude_below_base_m") is None
                or r.lower_m >= cfg["towers"][tower]["exclude_below_base_m"]
            )
        )
        d.update({f"param_{k}": v for k, v in prov.items()})
        rows.append(d)
    return pd.DataFrame(rows)


def run_baseline(cfg: dict) -> pd.DataFrame:
    """Slice every tower at the provisional settings and write
    results/baseline_measurements.csv."""
    prov = cfg["provisional_values"]
    frames = []
    for tower in cfg["towers"]:
        x, y, z, _ = load_cloud(tower, cfg)
        t0 = time.time()
        res = tm.slice_tower(
            x,
            y,
            z,
            thickness=prov["band_thickness_m"],
            min_points=prov["min_points_per_band"],
            min_coverage_deg=prov["min_angular_coverage_deg"],
            elong_threshold=prov["elongation_threshold"],
        )
        df = bands_to_frame(tower, res, prov, cfg)
        frames.append(df)
        log_event(
            "baseline_tower",
            tower=tower,
            bands=len(df),
            representative=int(df.representative.sum()),
            seconds=round(time.time() - t0, 2),
            params=prov,
        )
        print(
            f"{tower}: {len(df)} bands, "
            f"{int(df.representative.sum())} representative "
            f"({time.time() - t0:.1f} s)",
            flush=True,
        )
    out = pd.concat(frames, ignore_index=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    out.to_csv(RESULTS / "baseline_measurements.csv", index=False)
    return out


def main(argv: list[str]) -> int:
    """Command-line entry point. Always writes the dataset inventory; runs the
    baseline for the 'baseline' and 'all' commands."""
    cmd = argv[1] if len(argv) > 1 else "baseline"
    cfg = load_config()
    RESULTS.mkdir(parents=True, exist_ok=True)
    log_event(
        "start",
        command=cmd,
        python=platform.python_version(),
        numpy=np.__version__,
        pandas=pd.__version__,
    )
    inv = inventory(cfg)
    inv.to_csv(RESULTS / "dataset_inventory.csv", index=False)
    if cmd in ("baseline", "all"):
        run_baseline(cfg)
    log_event("end", command=cmd)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
