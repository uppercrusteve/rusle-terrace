# RUSLE 土壤侵蚀估算报告

- 生成时间 (UTC): 2026-11-16T18:00:04+00:00
- 数据: `examples\demo\dem.tif` — **合成演示 DEM，非实测**
- 方法: RUSLE 规划级筛查；公式、假设与局限见 `docs/methods.md`。

> 结论仅用于坡位/地块间相对比较和措施排序，不能替代小区实测。

## 参数

| 参数 | 数值 |
| --- | --- |
| annual_precip_mm | 620.0 |
| sand_pct | 45.0 |
| silt_pct | 35.0 |
| clay_pct | 20.0 |
| organic_carbon_pct | 1.2 |
| cover_factor | 0.25 |
| support_practice | terraces |
| slope_classes_deg | [0.0, 5.0, 10.0, 15.0, 25.0, 90.0] |

## 因子

- R = 1512.51 MJ mm ha⁻¹ h⁻¹ yr⁻¹
- K = 0.0426 t ha h ha⁻¹ MJ⁻¹ mm⁻¹
- C 均值 = 0.250
- P 均值 = 0.185

## 结果摘要

- 有效面积: 8061.12 ha
- 平均 A: 103.91 t/ha/yr
- 中位数 A: 8.22 t/ha/yr
- P90 A: 253.55 t/ha/yr
- 最大 A: 6931.33 t/ha/yr

## 坡度分级统计

| 坡度 (°) | 面积占比 (%) | 均值 (t/ha/yr) | P90 (t/ha/yr) | 分级 |
| --- | --- | --- | --- | --- |
| 0–5 | 76.20 | 9.10 | 16.13 | 中度 (moderate) |
| 5–10 | 0.45 | 145.33 | 187.20 | 极强烈 (very strong) |
| 10–15 | 1.25 | 174.96 | 372.91 | 极强烈 (very strong) |
| 15–25 | 22.11 | 425.85 | 707.88 | 极强烈 (very strong) |
| 25–90 | 0.00 | nan | nan | 无数据 (no data) |

## 警告

- 未填洼：D8 在洼地和平坦区可能截断或改道汇流，运行前应先填洼。
- C 与 P 是最大不确定性来源，必须用小区实测或本地文献率定。
- 结果为规划级相对筛查，不是小区实测替代。
- P 值假设梯田完好且全覆盖；梯田损毁会显著提高 P。
- DEM 为合成演示数据，不代表任何真实地块。

## 复现

```bash
rusle run DEM.tif params.yaml -o out/
```
