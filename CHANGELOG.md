# Changelog

All notable changes are recorded here. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Fixed

- C-factor cells replaced by the default value are counted and written to the
  report warnings instead of being substituted silently.

## [0.1.0] - 2026-11-17

Initial public release.

### Added

- `rusle run` computes A = R·K·LS·C·P on a single-band metric GeoTIFF DEM and writes `erosion.tif`, `erosion-map.png`, `slope-class-stats.csv`, `summary.json` and `report.md`.
- `rusle demo` generates a deterministic synthetic 300×300 @ 30 m terraced south-facing DEM, a parameter file and a `dem_truth.json` truth record, then executes the full chain.
- `rusle factors` computes R and K and prints the P table for a slope series without a DEM.
- DEM validation: rejects missing CRS, geographic CRS, non-metre linear units, non-square pixels, multi-band rasters and unsupported dtypes.
- D8 flow direction and topological flow accumulation with cycle fallback; LS from unit contributing area with m = 0.5, n = 1.3 and a 0.01 floor on flats.
- Optional C-factor raster aligned to the DEM, and `--p-factor` scalar override.
- Severity classes following a simplified SL190-2007 banding, plus per-slope-class area shares, means and P90.
- UTF-8 BOM CSV and JSON summary with explicit synthetic/external data labels and a warnings block.
- Test suite (69 tests) covering factor equations, analytic LS benchmarks, D8 routing, NoData, CRS rejection, output schemas and CLI exit codes; ruff lint clean.
- CI workflow on Ubuntu and Windows for Python 3.11 and 3.12.

### Known limitations

- No depression filling; D8 stops on flats and inside sinks.
- C and P are literature values, not locally calibrated.
- LS uses single-flow D8 without upslope-area cutoff.
- 30 m grids cannot resolve terrace risers.

[0.1.0]: https://github.com/uppercrusteve/rusle-terrace/releases/tag/v0.1.0
