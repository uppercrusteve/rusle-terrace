"""Map rendering and tabular/JSON/Markdown report writers."""

from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LightSource, LogNorm, to_rgba  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

from . import __version__  # noqa: E402
from .runoff import SEVERITY_BANDS, Analysis  # noqa: E402

SEVERITY_COLORS = ("#fff7bc", "#fec44f", "#fd8d3c", "#e34a33", "#99000d")
MAP_CMAP = "YlOrRd"
_FOOTER = "规划级筛查，非小区实测替代 | synthetic demo data"


def _setup_fonts() -> None:
    plt.rcParams["font.sans-serif"] = [
        "Microsoft YaHei",
        "SimHei",
        "Source Han Serif SC",
        "Noto Sans CJK SC",
        "DejaVu Sans",
    ]
    plt.rcParams["axes.unicode_minus"] = False


def _format_km(value: float) -> str:
    return f"{value / 1000.0:.0f}"


def plot_map(
    analysis: Analysis,
    out_png: str | Path,
    title: str = "坡地葡萄园土壤侵蚀估算（RUSLE）",
    synthetic: bool = False,
) -> Path:
    """Render hillshade + semi-transparent erosion overlay with legend and scale bar."""
    _setup_fonts()
    out = Path(out_png)
    dem = analysis.dem
    z = np.asarray(dem.data, dtype="float64")
    finite = np.isfinite(z)
    z_plot = np.where(finite, z, np.nanmin(z[finite]))
    valid = np.isfinite(analysis.erosion)
    values = analysis.erosion[valid]
    vmax = float(np.percentile(values, 99.5)) if values.size else 10.0
    vmax = max(vmax, 5.0)
    norm = LogNorm(vmin=0.5, vmax=vmax)
    cmap = plt.get_cmap(MAP_CMAP).with_extremes(bad=(0.0, 0.0, 0.0, 0.0))

    left, bottom = dem.transform.c, dem.transform.f
    right = left + dem.transform.a * dem.shape[1]
    top = bottom + dem.transform.e * dem.shape[0]
    extent = (left, right, bottom, top)

    fig = plt.figure(figsize=(10.0, 9.8), dpi=180)
    ax = fig.add_axes((0.075, 0.165, 0.87, 0.705))
    light = LightSource(azdeg=315.0, altdeg=38.0)
    shade = light.hillshade(
        z_plot, vert_exag=4.0, dx=dem.cellsize, dy=dem.cellsize, fraction=1.4
    )
    lo, hi = np.percentile(shade, [2.0, 98.0])
    shade = np.clip((shade - lo) / max(hi - lo, 1e-6), 0.0, 1.0)
    shade_cmap = plt.get_cmap("gray").with_extremes(bad=(1.0, 1.0, 1.0, 1.0))
    ax.imshow(
        np.ma.masked_where(~valid, shade),
        extent=extent,
        cmap=shade_cmap,
        vmin=0.0,
        vmax=1.0,
        origin="upper",
    )
    overlay = np.ma.masked_where(~valid, analysis.erosion)
    image = ax.imshow(
        overlay,
        extent=extent,
        cmap=cmap,
        norm=norm,
        origin="upper",
        alpha=0.55,
        interpolation="nearest",
    )

    # Contours keep the terrace geometry legible under the overlay.
    if np.isfinite(z).sum() > 0 and np.nanmax(z) > np.nanmin(z):
        levels = np.linspace(np.nanmin(z), np.nanmax(z), 14)
        ax.contour(z_plot, levels=levels, colors="#222222", linewidths=0.3, alpha=0.5,
                   extent=extent, origin="upper")

    ax.set_xlim(left, right)
    ax.set_ylim(bottom, top)
    ax.set_aspect("equal")
    ax.set_xlabel("东向坐标 Easting (km)", fontsize=10)
    ax.set_ylabel("北向坐标 Northing (km)", fontsize=10)
    x_ticks = np.linspace(left, right, 5)
    y_ticks = np.linspace(bottom, top, 5)
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.set_xticklabels([_format_km(v) for v in x_ticks], fontsize=8)
    ax.set_yticklabels([_format_km(v) for v in y_ticks], fontsize=8)
    ax.grid(color="white", linewidth=0.4, alpha=0.35)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

    mean_a = float(values.mean()) if values.size else float("nan")
    subtitle = "合成 DEM 演示（非实测）" if synthetic else "实测/外部 DEM"
    rows, cols = analysis.dem.shape
    fig.text(
        0.075,
        0.945,
        title,
        fontsize=17,
        fontweight="bold",
        ha="left",
        va="top",
    )
    fig.text(
        0.075,
        0.905,
        f"{subtitle} · {rows}×{cols} 格网 @ {analysis.dem.cellsize:g} m · "
        f"均值 A = {mean_a:.2f} t/ha/yr · R = {analysis.r:.0f} · K = {analysis.k:.3f}",
        fontsize=9.5,
        color="#444444",
        ha="left",
        va="top",
    )

    # Colorbar (log scale so treads and risers are both legible).
    cax = fig.add_axes((0.075, 0.083, 0.46, 0.017))
    cbar = fig.colorbar(image, cax=cax, orientation="horizontal", extend="both")
    ticks = [1.0, 2.0, 5.0, 10.0, 25.0, 50.0, 100.0, 250.0, 500.0]
    cbar.set_ticks([t for t in ticks if t < vmax])
    cbar.ax.set_xticklabels([f"{t:g}" for t in ticks if t < vmax])
    cbar.ax.xaxis.set_label_position("top")
    cbar.set_label("年均土壤流失量 A = R·K·LS·C·P  (t/ha/yr，对数色标)", fontsize=9, labelpad=3)
    cbar.ax.tick_params(labelsize=8)

    # Severity legend with class thresholds.
    handles = []
    labels = []
    for (label, lo, hi), color in zip(SEVERITY_BANDS, SEVERITY_COLORS, strict=True):
        handles.append(Patch(facecolor=color, edgecolor="#555555", linewidth=0.5))
        if np.isinf(hi):
            labels.append(f"{label}  >{lo:g}")
        else:
            labels.append(f"{label}  {lo:g}–{hi:g}")
    ax.legend(
        handles,
        labels,
        title="侵蚀强度分级 (t/ha/yr, SL190-2007 简化)",
        loc="upper right",
        bbox_to_anchor=(0.985, 0.852),
        fontsize=8,
        title_fontsize=8.5,
        framealpha=0.92,
        borderpad=0.6,
        handlelength=1.4,
    )

    # Scale bar (1 km), alternating segments.
    bar_frac = 1000.0 / (right - left)
    bar_x, bar_y, bar_h = 0.045, 0.035, 0.014
    ax.add_patch(
        Rectangle((bar_x, bar_y), bar_frac / 2, bar_h, transform=ax.transAxes,
                  facecolor="black", edgecolor="black", linewidth=0.6, zorder=6)
    )
    ax.add_patch(
        Rectangle((bar_x + bar_frac / 2, bar_y), bar_frac / 2, bar_h, transform=ax.transAxes,
                  facecolor="white", edgecolor="black", linewidth=0.6, zorder=6)
    )
    ax.text(bar_x + bar_frac / 2, bar_y + bar_h + 0.008, "1 km", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=8.5, zorder=6)
    ax.text(bar_x, bar_y - 0.006, "0", transform=ax.transAxes, ha="center", va="top",
            fontsize=7, zorder=6)

    # North arrow (grid north, north-up raster).
    ax.annotate(
        "",
        xy=(0.93, 0.975),
        xytext=(0.93, 0.888),
        xycoords="axes fraction",
        arrowprops={"facecolor": "black", "edgecolor": "black", "width": 2.2, "headwidth": 8.0},
        zorder=6,
    )
    ax.text(0.93, 0.982, "N", transform=ax.transAxes, ha="center", va="bottom",
            fontsize=11, fontweight="bold", zorder=6)

    # NoData patch in the map corner (top-left, clear of legend and scale bar).
    if not valid.all():
        ax.add_patch(
            Rectangle((0.102, 0.950), 0.024, 0.020, transform=ax.transAxes,
                      facecolor="#dddddd", edgecolor="#777777", linewidth=0.5, zorder=6)
        )
        ax.text(0.132, 0.960, "NoData", transform=ax.transAxes, va="center", fontsize=7.5,
                zorder=6)

    fig.text(0.075, 0.022, _FOOTER, fontsize=8, color="#666666", ha="left")
    fig.text(0.945, 0.022, "RUSLE terrace screening · v" + __version__,
             fontsize=8, color="#666666", ha="right")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def write_stats_csv(path: str | Path, stats: list[dict[str, Any]]) -> Path:
    out = Path(path)
    fields = [
        "class_index",
        "class_lo_deg",
        "class_hi_deg",
        "cell_count",
        "area_ha",
        "area_pct",
        "mean_erosion_t_ha_yr",
        "p90_erosion_t_ha_yr",
        "severity_by_mean",
    ]
    with out.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in stats:
            writer.writerow(row)
    return out


def summarize(
    analysis: Analysis,
    dem_path: str | Path,
    synthetic: bool,
    warnings: list[str],
) -> dict[str, Any]:
    values = analysis.valid_erosion
    cell_area_ha = analysis.dem.cellsize**2 / 10_000.0
    valid_area = float(values.size * cell_area_ha)
    total_area = float(analysis.dem.valid.sum() * cell_area_ha)
    severity = []
    for band in _severity_rows(values, valid_area):
        severity.append(band)
    return {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "software": {"name": "rusle-terrace", "version": __version__},
        "data": {
            "dem": str(dem_path),
            "synthetic": synthetic,
            "label": "合成演示 DEM（非实测）" if synthetic else "外部/实测 DEM",
            "crs": str(analysis.dem.crs),
            "cell_size_m": analysis.dem.cellsize,
            "shape": list(analysis.dem.shape),
            "total_area_ha": round(total_area, 4),
            "valid_area_ha": round(valid_area, 4),
        },
        "parameters": analysis.params.to_dict(),
        "factors": {
            "R_MJ_mm_ha_h_yr": round(analysis.r, 4),
            "K_t_ha_h_ha_MJ_mm": round(analysis.k, 6),
            "C_mean": float(np.nanmean(analysis.c_factor)),
            "P_mean": float(np.nanmean(analysis.p_factor)),
        },
        "results_t_ha_yr": {
            "mean": round(float(values.mean()), 4) if values.size else None,
            "median": round(float(np.median(values)), 4) if values.size else None,
            "p90": round(float(np.percentile(values, 90)), 4) if values.size else None,
            "max": round(float(values.max()), 4) if values.size else None,
            "valid_cells": int(values.size),
        },
        "severity_area_pct": severity,
        "warnings": warnings,
    }


def _severity_rows(values: np.ndarray, valid_area: float) -> list[dict[str, Any]]:
    rows = []
    lower = 0.0
    for label, _, upper in SEVERITY_BANDS:
        mask = (values >= lower) & (values < upper)
        area = float(mask.sum() * valid_area / values.size) if values.size else 0.0
        pct = 100.0 * mask.sum() / values.size if values.size else 0.0
        rows.append(
            {
                "severity": label,
                "lower_t_ha_yr": lower,
                "upper_t_ha_yr": None if np.isinf(upper) else upper,
                "area_pct": round(float(pct), 3),
                "area_ha": round(area, 4),
            }
        )
        lower = upper
    return rows


def write_summary_json(path: str | Path, summary: dict[str, Any]) -> Path:
    out = Path(path)
    out.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return out


def write_report_md(
    path: str | Path,
    analysis: Analysis,
    dem_path: str | Path,
    synthetic: bool,
    warnings: list[str],
) -> Path:
    values = analysis.valid_erosion
    cell_area_ha = analysis.dem.cellsize**2 / 10_000.0
    lines = [
        "# RUSLE 土壤侵蚀估算报告",
        "",
        f"- 生成时间 (UTC): {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"- 数据: `{dem_path}` — {'**合成演示 DEM，非实测**' if synthetic else '外部/实测 DEM'}",
        "- 方法: RUSLE 规划级筛查；公式、假设与局限见 `docs/methods.md`。",
        "",
        "> 结论仅用于坡位/地块间相对比较和措施排序，不能替代小区实测。",
        "",
        "## 参数",
        "",
        "| 参数 | 数值 |",
        "| --- | --- |",
    ]
    for key, value in analysis.params.to_dict().items():
        lines.append(f"| {key} | {value} |")
    lines += [
        "",
        "## 因子",
        "",
        f"- R = {analysis.r:.2f} MJ mm ha⁻¹ h⁻¹ yr⁻¹",
        f"- K = {analysis.k:.4f} t ha h ha⁻¹ MJ⁻¹ mm⁻¹",
        f"- C 均值 = {np.nanmean(analysis.c_factor):.3f}",
        f"- P 均值 = {np.nanmean(analysis.p_factor):.3f}",
        "",
        "## 结果摘要",
        "",
        f"- 有效面积: {values.size * cell_area_ha:.2f} ha",
    ]
    if values.size:
        lines += [
            f"- 平均 A: {values.mean():.2f} t/ha/yr",
            f"- 中位数 A: {np.median(values):.2f} t/ha/yr",
            f"- P90 A: {np.percentile(values, 90):.2f} t/ha/yr",
            f"- 最大 A: {values.max():.2f} t/ha/yr",
        ]
    lines += [
        "",
        "## 坡度分级统计",
        "",
        "| 坡度 (°) | 面积占比 (%) | 均值 (t/ha/yr) | P90 (t/ha/yr) | 分级 |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in analysis.stats:
        lines.append(
            f"| {row['class_lo_deg']:g}–{row['class_hi_deg']:g} "
            f"| {row['area_pct']:.2f} | {row['mean_erosion_t_ha_yr']:.2f} "
            f"| {row['p90_erosion_t_ha_yr']:.2f} | {row['severity_by_mean']} |"
        )
    lines += ["", "## 警告", ""]
    lines += [f"- {text}" for text in warnings] if warnings else ["- 无"]
    lines += ["", "## 复现", "", "```bash", "rusle run DEM.tif params.yaml -o out/", "```", ""]
    out = Path(path)
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def legend_color(label: str) -> tuple[float, float, float, float]:
    """Public helper returning the severity colour, used by the map tests."""
    for (band_label, _, _), color in zip(SEVERITY_BANDS, SEVERITY_COLORS, strict=True):
        if band_label == label:
            return to_rgba(color)
    raise KeyError(label)
