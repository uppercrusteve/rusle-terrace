# Contributing

## 开发环境

需要 Python 3.10+，以及可用的 GDAL/rasterio。Linux 先装系统 GDAL 库；Windows 与 macOS 直接装 wheel 即可。

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .\.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
rusle --help
```

## 提交前

```bash
python -m ruff check .
python -m pytest -q
```

两条命令必须零告警、全绿。新增或改动公式、单位、分级阈值时，同步更新 `docs/methods.md` 与对应手算夹具测试；只改代码不改口径文档的 PR 不合并。

## 数据口径

- DEM 必须是**单波段、米制投影 CRS、方形像素**的 GeoTIFF。经纬度或英尺制 DEM 先重投影，程序不做内部重投影。
- 参数文件是 UTF-8 YAML，六个必填键：`annual_precip_mm`、`sand_pct`、`silt_pct`、`clay_pct`、`organic_carbon_pct`、`cover_factor`。
- `sand_pct + silt_pct + clay_pct` 必须为 100 ± 1；`cover_factor` 与 `--p-factor` 在 0–1。
- 面积单位统一 ha，侵蚀量统一 t/ha/yr；网格坐标与高程单位 m。
- 示例、夹具、截图必须注明“合成演示”，不得伪装成观测数据。合成 DEM 由 `rusle demo` 解析式生成，`dem_truth.json` 记录解析坡度的真值。
- 没有本地小区实测的 C、P 不得在文档里写成“已率定”。

## Issue

请附：DEM 元数据（CRS、分辨率、范围，可脱敏）、`params.yaml` 关键字段、完整命令、退出码、期望与实际输出、Python 与 rasterio 版本。数值差异请先按 `docs/methods.md` 给出独立手算，便于定位是公式还是数据问题。

## PR

- 一个 PR 只做一件事；描述里贴 `ruff` 与 `pytest` 输出。
- 依赖保持精简，不引入 Web 框架或 GIS 桌面栈。
- 新的分级或因子表必须在 `docs/methods.md` 给出出处和适用坡度范围。
- 提交信息用 `feat:` / `fix:` / `test:` / `docs:` 前缀。

## License

贡献即同意以 MIT 许可发布。
