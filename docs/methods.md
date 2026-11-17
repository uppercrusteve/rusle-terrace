# 方法与局限 / Methods and limitations

本文档记录 `rusle-terrace` 的全部公式、单位、实现选择、文献出处和已知偏差。代码里出现的每个常数都能在这里找到对应说明。**这是规划级筛查工具，不是小区实测的替代品。**

---

## 1. 总体框架

采用修正通用土壤流失方程 RUSLE（Renard et al., 1997）的乘积形式：

```
A = R · K · LS · C · P
```

| 符号 | 含义 | 单位 |
| --- | --- | --- |
| A | 年均土壤流失量 | t ha⁻¹ yr⁻¹ |
| R | 降雨侵蚀力 | MJ mm ha⁻¹ h⁻¹ yr⁻¹ |
| K | 土壤可蚀性 | t ha h ha⁻¹ MJ⁻¹ mm⁻¹ |
| LS | 坡长-坡度复合因子 | 无量纲 |
| C | 植被与地表覆盖因子 | 无量纲，0–1 |
| P | 水土保持措施因子 | 无量纲，0–1 |

RUSLE 假设各因子独立可乘、长期平均、片蚀与细沟蚀为主；不计算沟蚀、崩岗、风蚀和耕作位移。对单场暴雨和年际波动没有解释力。

实现位置：`src/rusle_terrace/runoff.py::analyze` 中的 `erosion = r * k * ls * c * p`。无效像元在乘积后统一置 NaN，写出栅格时改为 `-9999`。

---

## 2. 降雨侵蚀力 R

```
R = 0.0483 · P^1.61
```

- `P` 为年降水量，mm；`R` 的单位是 MJ mm ha⁻¹ h⁻¹ yr⁻¹。
- 出处：Renard & Freimund (1994), equation 2，基于美国西部站点回归。适用于年降水约 250–1500 mm 的干旱半干旱区；湿润区偏高，热带季风区需谨慎。
- 实现：`factors.py::r_factor`。要求有限且 > 0，否则抛 `ParameterError`。

**局限**：只用年降水量，不反映雨强、雨滴动能和季节分配；不能代替基于自记雨量计或日雨量模型的 R 值计算。云南高原的干湿季极端降雨可能被低估。

---

## 3. 土壤可蚀性 K

采用 EPIC 模型（Williams, 1995）的土壤可蚀性近似：

```
K_US = 0.2 + 0.3 · exp[−0.0256 · SAN · (1 − SIL/100)]
            · (SIL / (CLA + SIL))^0.3
            · [1 − 0.25C / (C + exp(3.72 − 2.95C))]
            · [1 − 0.7 · SNI / (SNI + exp(−5.51 + 22.9 · SNI))]

K = 0.1317 · K_US
```

- `SAN`、`SIL`、`CLA` 为砂、粉、黏粒的质量百分数（0–100，三者合计 100 ± 1）。
- `C` 为有机碳**质量分数**（`organic_carbon_pct / 100`，输入按 % 给出）。
- `SNI = 1 − SAN/100`。
- `0.1317` 是把美制单位（ton acre h / (100 ft ton acre in)）换算到 SI 的系数。
- 出处：Williams (1995)；单位换算见 Renard et al. (1997) 附录。
- 实现：`factors.py::k_factor`。校验三项质地均在 0–100、有机碳在 0–100、`SIL + CLA > 0`。

**局限**：EPIC 是经验公式，对火山灰土壤和强团聚体土壤偏差较大（常高估可蚀性）；有机碳项单调压低 K，高有机质土壤可能反向低估；`0.2` 截距意味着式子在美制单位下恒 > 0.2。有实测 K（小区或基准公式）时应改用实测值。本工具只有单一 K 标量，没有 K 的空间栅格入口。

---

## 4. 坡长-坡度因子 LS

### 4.1 公式

```
LS = (m + 1) · (As / 22.13)^m · (sin β / 0.0896)^n
```

- `As` 为单位宽度上游汇水面积，m²/m（数值等于米）。`As = acc · cellsize`，`acc` 是 D8 汇流累积的像元数（含自身），`cellsize` 为像元边长 m。
- `β` 为坡度角（弧度）。`22.13` 是标准小区坡长 22.13 m；`0.0896 = sin(5.14°)` 是标准小区坡度。
- `m = 0.5`，`n = 1.3`，取自 Wischmeier & Smith (1978) 与 Renard et al. (1997) 的常用取值，**不随坡度或细沟/片蚀比调整**。
- 坡度 ≤ 0 时 `sin β = 0` 会让 LS 归零，平地统一取 `LS = 0.01`。
- 实现：`runoff.py::ls_factor`，常数 `M_EXP`、`N_EXP`、`LS_SLOPE_ZERO`。

### 4.2 为什么用单位汇水面积

把坡长改写为上游汇水面积，可以把二维汇流信息引入 LS。这是 GIS-RUSLE 的常规做法（Desmet & Govers, 1996；Van Remortel et al., 2004），代价是假设汇流路径即坡长方向。发生淤积、分流的坡面会得到偏大的 `As`。

### 4.3 实现细节与已知偏差

- 在计算坡度与流向之前，所有无效像元（含边界与内部 NoData）用最近邻有效值填充，避免边缘出现虚假陡坡、内部空洞产生假洼地；最终输出仍把真实无效区置 NaN 并掩膜。
- 只做单流向 D8，不累积多向分配，`As` 对流向误差敏感。
- 设固定 m、n，不按坡度分段（陡坡段 n 通常应更高），陡坡可能低估。
- 没有上坡面积上限（cutoff），长坡会持续累积到很大的 LS。
- 30 m 栅格上的 LS 是像元平均意义的陡坡带，不是垄沟尺度的真实坡长。

### 4.4 验证

`tests/test_runoff.py` 用解析平面（均匀 9% 坡降）对照解析 LS，要求相对误差 < 5%，并单独验证平地取 0.01、`As` 单调性和 NaN 掩膜。这是数值一致性验证，不是野外验证。

---

## 5. 坡度与 D8 流向

### 5.1 坡度

中心差分：

```
slope = atan( hypot(∂z/∂x, ∂z/∂y) ) · 180/π
```

实现：`flow.py::slope_degrees`，`cellsize` 必须为正。中心差分在台坎边缘会把陡坎平滑成 2–3 个像元的过渡带，真实坎脚坡度被低估。

### 5.2 D8 流向

- 8 邻域，取 `(z − z_neighbor) / (distance · cellsize)` 最大的方向，即最陡下降；对角距离乘 √2，保证平面坡面走正交方向。
- 等坡降时按方向顺序（东、东南、南、西南、西、西北、北、东北）取第一个，偏置向东。
- 无下坡邻域（洼地、平地、出水口）记 `-1`。
- 实现：`flow.py::d8_directions`。

### 5.3 汇流累积

`flow.py::flow_accumulation` 在 D8 图上做拓扑队列，每个像元计 1（自身）再加全部上游像元数。出水口手动入队；因未填洼而产生的环路在队列耗尽后按索引反序补推一次，保证水量不丢。**没有填洼（depression filling）**：洼地与平地都记为无方向（`-1`），上游汇流在洼地/平地边缘被截断；等坡降平局则按方向顺序取东，天然偏向东西向路由。

**结论**：运行前应先在外部 GIS 完成填洼与（必要时）平地抬升；本工具不修改输入 DEM。

---

## 6. 覆盖因子 C

两种入口：

1. `cover_factor` 标量（YAML 必填，0–1），全场均一。适合整片园地同一管理方式的粗筛。
2. `--c-factor-raster`：单波段栅格，形状、仿射变换必须与 DEM 完全一致，值域 0–1，无效像元回退到 `cover_factor`。适合由 NDVI 或物候曲线生成的月/季 C 图。

`C` 是 RUSLE 中随管理剧烈变化的项。葡萄园行间生草、行内清耕、冬季裸土、覆盖物和树龄都会造成数量级差异。

**局限**：C 的文献值多为国外作物或通用土地覆盖表。**必须用本地小区实测（天然降雨或模拟降雨）或至少本地文献率定**；直接套用默认值只能得到相对排序。单一标量也无法表达行间/行内的空间结构。

---

## 7. 措施因子 P

按坡度分箱的查表值（`factors.py::P_TABLE`）：

| support_practice | 0–5° | 5–10° | 10–15° | ≥15° |
| --- | --- | --- | --- | --- |
| `none` | 1.00 | 1.00 | 1.00 | 1.00 |
| `contour` | 0.55 | 0.60 | 0.65 | 0.70 |
| `terraces` | 0.15 | 0.20 | 0.25 | 0.30 |

出处：Wischmeier & Smith (1978) 与 Renard et al. (1997) 的等高耕作/梯田表，按坡度内插为阶梯值。

`--p-factor` 可用 0–1 的标量覆盖整张表（例如梯田大面积损毁时统一提高 P）。覆盖时报告与 `summary.json` 会写入警告。

**局限**：

- 表值假设措施**完好、连续、全覆盖**；梯田垮塌、断行、排水口冲刷都会显著提高实际 P。
- 分箱是阶梯函数，5° 与 5.0001° 会跳级。
- 没有区分石坎/土坎、坎高、田面宽度和坡度。云南坡地梯田的规格差异很大，P 的本地率定优先于任何查表。

---

## 8. 坡度分级与侵蚀强度

默认坡度分级（`slope_classes_deg`，可改，但必须以 0 开始、严格递增、以 ≥90 结束）：

```
0–5, 5–10, 10–15, 15–25, 25–90°
```

每个分级输出像元数、面积 ha、面积占比、A 均值与 P90，并按**分级均值**给出强度等级。空分级记 NaN 与“无数据”。

侵蚀强度分级（`runoff.py::SEVERITY_BANDS`）为 SL190-2007《土壤侵蚀分类分级标准》水力侵蚀面蚀分级的简化：

| 等级 | A (t ha⁻¹ yr⁻¹) |
| --- | --- |
| 微度 very slight | 0–2 |
| 轻度 slight | 2–5 |
| 中度 moderate | 5–10 |
| 强烈 strong | 10–25 |
| 极强烈 very strong | > 25 |

上界左闭右开。SL190-2007 原表按 5 个侵蚀类型区分别定阈值，这里只取一套并集式的通用带，**不是标准原文的逐区阈值**，跨区引用前应核对原文。

---

## 9. 输入、坐标与单位

- DEM：单波段 GeoTIFF，必须是**投影 CRS 且线性单位为米**，方形、正向分辨率像素。经纬度 CRS、英尺制 CRS、多波段、非方形像素都会被拒绝并返回退出码 2。
- 不做内部重投影。拿到 EPSG:4326 的 DEM 先在 QGIS/GDAL 里 `gdalwarp -t_srs <metric UTM>` 并重采样。
- `cellsize` 只读取 `transform.a`；非方形像素已被排除。
- 高程单位假定为米。使用英尺高程的 DEM 需先换算，否则坡度与 `As` 都会错。
- 面积由像元面积 `cellsize²` 换算为 ha，除以 10 000。

---

## 10. DEM 分辨率的影响

- **汇水面积随分辨率平方变化**：像元变粗，`As` 以 `cellsize` 线性增大，`LS` 随 `(As)^0.5` 增大。同一坡面在 90 m 与 10 m 上得到的 A 不可直接比较。
- **台坎不可分辨**：30 m 像元大于多数梯田坎宽（演示 DEM 的坎面仅 7 m，田面 23 m），坎坡被平滑成缓坡，坎脚侵蚀带被摊薄。要刻画坎脚，建议 ≤ 10 m；30 m 只适合区域筛查与地块间排序。
- **流向不确定性**：平地与微洼在粗分辨率下更普遍，D8 截断更多。
- **建议**：同一项目内固定分辨率；跨分辨率比较前先重采样到统一网格并重算，而不是按面积加权旧结果。

演示数据为 300×300 @ 30 m 的解析合成地表（南向、三级台地、两条横向波），`dem_truth.json` 记录解析坡度真值（1.72–17.88°），**不是实测 DEM**。

---

## 11. 输出口径

| 文件 | 内容 |
| --- | --- |
| `erosion.tif` | float32，deflate，NoData = −9999，波段描述 “RUSLE A (t/ha/yr), planning screening” |
| `erosion-map.png` | 山体阴影 + 对数色标侵蚀叠加 + 强度图例 + 比例尺 + 指北针 |
| `slope-class-stats.csv` | 逐坡度分级：像元数、面积 ha、占比、均值、P90、按均值定的等级 |
| `summary.json` | schema 1.0：数据标签（合成/外部）、CRS、面积、参数、R/K/C/P、分位数、强度面积占比、警告 |
| `report.md` | 中文报告，含参数表、因子、坡度分级表、警告与复现命令 |

`summary.json` 的 `warnings` 和报告固定包含“未填洼”“C/P 需本地率定”“非小区实测替代”三条；梯田模式与合成 DEM 会追加对应警告。合成判定靠 DEM 同目录的 `dem_truth.json` 中 `synthetic: true`，不靠文件名。

---

## 12. 与其他实现的差异

- 不含 USLE 的 P 值以外措施（带状耕作、等高缓冲）细分。
- 不含坡度分段 m/n（USLE-M、RUSLE2 的细沟/片蚀比）。
- 不含汇流面积上限、沉积联结、输沙比。
- 不含 terraces 的几何改造（把台地修平），只在 P 上打折。

因此结果适合：地块内坡位排序、措施情景的相对比较、复垦/整地优先级。不适合：绝对流失量申报、小区试验替代、沟蚀量估算、单场暴雨设计。

---

## 13. 参考文献

- Desmet, P.J.J., Govers, G. (1996). A GIS procedure for automatically calculating the USLE LS factor on topographically complex landscape units. *Journal of Soil and Water Conservation*, 51(5), 427–433.
- Renard, K.G., Foster, G.R., Weesies, G.A., McCool, D.K., Yoder, D.C. (1997). *Predicting Soil Erosion by Water: A Guide to Conservation Planning with the Revised Universal Soil Loss Equation (RUSLE)*. USDA Agriculture Handbook 703.
- Renard, K.G., Freimund, J.R. (1994). Using monthly precipitation data to estimate the R-factor in the revised USLE. *Journal of Hydrology*, 157(1–4), 287–306. https://doi.org/10.1016/0022-1694(94)90110-4
- Van Remortel, R.D., Maichle, R.W., Hickey, R.J. (2004). Computing the LS factor for the Revised Universal Soil Loss Equation through array-based slope processing of digital elevation data using a C++ executable. *Computers & Geosciences*, 30(9–10), 1043–1053. https://doi.org/10.1016/j.cageo.2004.08.001
- Williams, J.R. (1995). The EPIC model. In: Singh, V.P. (ed.) *Computer Models of Watershed Hydrology*, 909–1000. Water Resources Publications.
- Wischmeier, W.H., Smith, D.D. (1978). *Predicting Rainfall Erosion Losses: A Guide to Conservation Planning*. USDA Agriculture Handbook 537.
- 中华人民共和国水利部 (2007). *土壤侵蚀分类分级标准* SL190-2007. 中国水利水电出版社.

---

## 14. 一句话局限

**R、K、LS 是公式近似，C、P 是未经本地率定的文献值；30 m 合成演示 DEM 只证明流程可跑，所有绝对数值在真实地块上都应重算并重新率定 C/P。**
