"""Command-line interface for the RUSLE terrace screening tool."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path

import numpy as np

from . import __version__
from .factors import PRACTICES, ParameterError, k_factor, load_params, p_factor, r_factor
from .report import plot_map, summarize, write_report_md, write_stats_csv, write_summary_json
from .runoff import DemError, analyze, write_erosion
from .synth import generate_demo, is_synthetic_dem

PROG = "rusle"


def _warnings(support_practice: str, synthetic: bool) -> list[str]:
    notes = [
        "未填洼：D8 在洼地和平坦区可能截断或改道汇流，运行前应先填洼。",
        "C 与 P 是最大不确定性来源，必须用小区实测或本地文献率定。",
        "结果为规划级相对筛查，不是小区实测替代。",
    ]
    if support_practice == "terraces":
        notes.append("P 值假设梯田完好且全覆盖；梯田损毁会显著提高 P。")
    if synthetic:
        notes.append("DEM 为合成演示数据，不代表任何真实地块。")
    return notes


def _run(
    dem_path: Path,
    params_path: Path,
    outdir: Path,
    c_raster: Path | None,
    p_override: float | None,
    title: str,
) -> int:
    params = load_params(params_path)
    if c_raster is not None and not c_raster.is_file():
        raise DemError(f"C-factor raster not found: {c_raster}")
    synthetic = is_synthetic_dem(dem_path)
    analysis = analyze(dem_path, params, c_raster=c_raster, p_override=p_override)
    outdir.mkdir(parents=True, exist_ok=True)
    warnings = _warnings(params.support_practice, synthetic)
    warnings.extend(analysis.notes)
    if p_override is not None:
        warnings.append(f"P 已被命令行覆盖为 {p_override:g}，未使用 {params.support_practice} 表。")

    erosion_path = write_erosion(outdir / "erosion.tif", analysis)
    summary = summarize(analysis, dem_path, synthetic, warnings)
    write_summary_json(outdir / "summary.json", summary)
    write_stats_csv(outdir / "slope-class-stats.csv", analysis.stats)
    write_report_md(outdir / "report.md", analysis, dem_path, synthetic, warnings)
    plot_map(analysis, outdir / "erosion-map.png", title=title, synthetic=synthetic)

    values = analysis.valid_erosion
    print(f"DEM: {dem_path}  ({analysis.dem.shape[0]}x{analysis.dem.shape[1]} @ "
          f"{analysis.dem.cellsize:g} m, {analysis.dem.crs})")
    print(f"Data: {'SYNTHETIC demo' if synthetic else 'external/observed'}")
    print(f"R = {analysis.r:.2f}   K = {analysis.k:.4f}   "
          f"C mean = {np.nanmean(analysis.c_factor):.3f}   "
          f"P mean = {np.nanmean(analysis.p_factor):.3f}")
    if values.size:
        print(f"A (t/ha/yr): mean {values.mean():.3f}  median {np.median(values):.3f}  "
              f"P90 {np.percentile(values, 90):.3f}  max {values.max():.3f}")
    print(f"Outputs written to: {outdir.resolve()}")
    print(f"  {erosion_path.name}, erosion-map.png, slope-class-stats.csv, summary.json, report.md")
    if synthetic:
        print("NOTE: synthetic DEM — not observed data.")
    return 0


def _cmd_run(args: argparse.Namespace) -> int:
    return _run(args.dem, args.params, args.out, args.c_factor_raster, args.p_factor, args.title)


def _cmd_demo(args: argparse.Namespace) -> int:
    outdir = args.out
    files = generate_demo(outdir)
    truth = json.loads(files["truth"].read_text(encoding="utf-8"))
    slope = truth["analytic_slope_deg"]
    print("Generated synthetic 300x300 @ 30 m south-facing terraced DEM "
          f"(analytic slope {slope['min']}–{slope['max']} deg).")
    return _run(files["dem"], files["params"], outdir, None, None, args.title)


def _cmd_factors(args: argparse.Namespace) -> int:
    r = r_factor(args.precip)
    k = k_factor(args.sand, args.silt, args.clay, args.oc)
    payload = {
        "R_MJ_mm_ha_h_yr": round(r, 4),
        "K_t_ha_h_ha_MJ_mm": round(k, 6),
        "SNI": round(1.0 - args.sand / 100.0, 4),
        "organic_carbon_fraction": round(args.oc / 100.0, 5),
        "P_by_slope_deg": {
            f"{deg}": float(p_factor(deg, args.practice)) for deg in (3.0, 7.0, 12.0, 20.0)
        },
        "practice": args.practice,
    }
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=PROG,
        description="RUSLE soil-erosion screening for terraced hillslope vineyards "
        "(planning level, metric DEM required). / 坡地梯田葡萄园侵蚀筛查。",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="run RUSLE for a DEM and parameter file",
                         description="Compute A = R*K*LS*C*P and write rasters and reports.")
    run.add_argument("dem", type=Path, help="single-band metric GeoTIFF DEM")
    run.add_argument("params", type=Path, help="YAML parameter file")
    run.add_argument("-o", "--out", type=Path, required=True, help="output directory")
    run.add_argument("--c-factor-raster", type=Path, default=None,
                     help="optional C raster aligned to the DEM (overrides cover_factor)")
    run.add_argument("--p-factor", type=float, default=None,
                     help="override P with a scalar in 0-1")
    run.add_argument("--title", default="坡地葡萄园土壤侵蚀估算（RUSLE）",
                     help="map title")
    run.set_defaults(func=_cmd_run)

    demo = sub.add_parser(
        "demo",
        help="generate a synthetic DEM and run end-to-end",
        description="Create a synthetic terraced DEM, parameters and full outputs.",
    )
    demo.add_argument("-o", "--out", type=Path, default=Path("demo_out"),
                      help="output directory (default: demo_out)")
    demo.add_argument("--title", default="坡地葡萄园土壤侵蚀估算 · 合成 DEM 演示",
                      help="map title")
    demo.set_defaults(func=_cmd_demo)

    factors = sub.add_parser("factors", help="print R, K and P factors",
                             description="Compute factors without a DEM.")
    factors.add_argument("--precip", type=float, required=True, help="annual precipitation (mm)")
    factors.add_argument("--sand", type=float, required=True, help="sand (%)")
    factors.add_argument("--silt", type=float, required=True, help="silt (%)")
    factors.add_argument("--clay", type=float, required=True, help="clay (%)")
    factors.add_argument("--oc", type=float, required=True, help="organic carbon (%)")
    factors.add_argument("--practice", choices=PRACTICES, default="none",
                         help="support practice for the P example (default: none)")
    factors.set_defaults(func=_cmd_factors)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (ParameterError, DemError, ValueError) as exc:
        print(f"{PROG}: error: {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"{PROG}: I/O error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
