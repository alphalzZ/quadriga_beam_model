# 14 个 sheet 波束的 XLSX 方向图

根据 `beam_model_sheet.mat` 的 14 组 Pitch/Azimuth，使用当前时域链路相同的全数据拟合 QuaDRiGa 天线模型导出。载频 1.97 GHz、502 个三角栅格接收阵子、理想 RHCP；角度直接使用 θ=ele、φ=azi。**这是拟合模型产生的数据，不是新增的实测数据。** 仅 C9 的 50.37° 扫描需要从 50° 边界状态外推 0.37°。

每个文件遵循原始 `pattern_Pitch_0_Azimuth_0.txt` 的表格布局：

- 文件名 `pattern_Pitch_<角度>_Azimuth_<角度>.xlsx`；对应 C1～C14 的关系见 `beam_index.csv`。
- 单工作表 `Sheet1`，共 92 行、361 列；A1 为 `theta\phi`。
- 第一行 B1～MW1 为 φ=0～359°；第一列 A2～A92 为 θ=0～90°，均为 1° 步进。
- B2～MW92 为 91×360 个归一化增益 dB 数值；每份完整方向图的采样峰值归一化为 0 dB，保留两位小数。
- 为匹配原始数据的输出动态范围，导出时将低于 −65 dB 的值截为 −65 dB。不会人为复制原文件中的 −60 dB 钳位片区或固定奇异点；不会改变拟合模型本身。

数值为 QuaDRiGa 复电场插值后的 `|Fa|²+|Fb|²`。与前一版极坐标切面共用同一个天线模型；XLSX 按完整半球网格峰值归一化，极坐标图按两条切面采样峰值归一化，两者参考值可能略有差异。原始 1° 网格不一定恰好经过文件名中的小数指向，因此网格峰值位置可与指向相差少量角度。数据只有归一化功率；绝对增益、复相位、极化纯度仍不是实测量。

每份表的 A1 批注及 `manifest.json` 保存模型来源、指向、归一化参考、截底点数和外推标志。验证程序重新打开全部文件，逐一检查所有 458,640 个增益数值、轴标签、表格尺寸和 0 dB 峰值，同时验证 ZIP。

## 重建

从仓库根目录运行：

```bash
python -m pip install -r experiments/measured_array_1970/pattern_xlsx/requirements.txt
octave --quiet --eval "addpath('experiments/measured_array_1970/pattern_xlsx'); export_sheet_pattern_grid('/path/to/QuaDRiGa');"
python experiments/measured_array_1970/pattern_xlsx/export_xlsx.py
```

所有输出位于本模块 `results/`。`14_sheet_beam_patterns_xlsx.zip` 包含 14 个 XLSX、索引和说明；`sheet_pattern_grid.mat` 保留归一化、截底和取整前的完整线性功率，便于改变导出规则。
