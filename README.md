# QuaDRiGa 波束建模与链路仿真

新手请先阅读 [入门 Handbook](doc/handbook.md)：从基础波束仿真、实测数据拟合到方向图导出和链路级仿真。原理与设计文档统一收录在 [doc/](doc/README.md)。

按功能组织代码；每个模块维护自己的入口、说明和 `results/`，原始数据与共享函数单独保存。

| 模块 | 用途 | 运行入口 |
|---|---|---|
| [矩形阵波束码本](src/beam_codebook/README.md) | 3.5 GHz、32×16 阵列，photo/sheet 角度缓存 | `run_beam_model.m` |
| [三角栅格与圆极化](src/beam_model_triangular_cp/README.md) | 可选栅格、LHCP/RHCP/双极化及 QuaDRiGa 适配 | `run_beam_model_triangular_cp.m` |
| [频域链路](src/link_frequency_domain/README.md) | 原始频域信道与多用户 OFDM 示例 | `iface_demo.m` |
| [时域链路](src/link_time_domain/README.md) | 512 阵子时域卷积与逐采样波束合并 | `iface_demo_time_domain.m` |
| [实测物理基线](experiments/measured_array_1970/baseline_fit/README.md) | 1.97 GHz 排布搜索及留出比较 | `fit_model.py`、`fit_ellipse.py` |
| [全数据拟合](experiments/measured_array_1970/all_data_fit/README.md) | 全部 78 个扫描状态、全部方向图采样点拟合 | `fit_all_data.py` |
| [sheet XLSX 方向图](experiments/measured_array_1970/pattern_xlsx/README.md) | 14 个完整半球方向图 Excel 导出与 ZIP | `export_sheet_pattern_grid.m`、`export_xlsx.py` |
| [sheet 极坐标方向图](experiments/measured_array_1970/beam_patterns/README.md) | 14 个拟合波束的离轴/方位切面与 PDF | `export_sheet_pattern_cuts.m`、`plot_sheet_polar.py` |
| [sheet 时域链路](experiments/measured_array_1970/time_domain/README.md) | 全数据拟合的 502 阵子、sheet 离轴 Pitch/Azimuth | `iface_allfit_time_domain.m` |

```text
doc/                                  入门手册、原理与设计说明
data/measured_patterns/                原始方向图及数据说明
src/
  common/                              通用方向角、阵列权重、阵元增益函数
  beam_codebook/                       矩形阵码本生成与 results/
  beam_model_triangular_cp/             三角栅格/圆极化示例与 results/
  link_frequency_domain/               频域链路与 results/
  link_time_domain/                    时域链路（运行时写入 results/）
experiments/measured_array_1970/
  common/                              实测数据、拟合数学、QuaDRiGa 适配函数
  baseline_fit/                        物理基线及独立 results/
  all_data_fit/                        全数据拟合及独立 results/
  time_domain/                         sheet 时域链路及独立 results/
  beam_patterns/                       sheet 极坐标切面及独立 results/
  pattern_xlsx/                        sheet 完整方向图 XLSX 导出及独立 results/
```

本次 sheet 链路直接使用用户确认的 `theta=ele_tab`、`phi=azi_tab`；13 个波束位于 0～50° 扫描范围内，仅 C9 的 50.37° 需要外推。完整 4×14 符号验证及逐波束结果见该模块说明。

MATLAB/Octave 入口通过各自目录加入路径；QuaDRiGa 使用 `QD_PATH`、函数的 `qd_path` 参数或预先加入的源码路径。Python 依赖见 [requirements.txt](experiments/measured_array_1970/requirements.txt)。数据路径和默认输出路径相对模块文件解析，可从其他工作目录运行。

## 目录迁移

| 原位置 | 当前位置 |
|---|---|
| `src/beam_model_photo/run_beam_model.m` | `src/beam_codebook/run_beam_model.m` |
| `src/beam_model_photo/beam_model_*.mat` | 对应码本模块的 `results/` |
| `src/beam_model_photo/iface_demo.m` | `src/link_frequency_domain/iface_demo.m` |
| `src/beam_model_photo/iface_demo_time_domain.m` | `src/link_time_domain/iface_demo_time_domain.m` |
| `src/beam_model_photo/*triangular*`、圆极化适配/绘图 | `src/beam_model_triangular_cp/` |
| `experiments/measured_array_1970/model.py`、阵列适配函数 | `experiments/measured_array_1970/common/` |
| 实测项目根目录的早期拟合脚本及 `results/` | `experiments/measured_array_1970/baseline_fit/` |
| 原模块 `iface_out.mat` | `src/link_frequency_domain/results/iface_out_module_legacy.mat` |
| 根目录 `iface_out.mat` | `src/link_frequency_domain/results/iface_out_root_legacy.mat` |

已有调用需按表更新 `addpath` 或脚本路径。原始数据和历史结果均保留；两个历史 `iface_out.mat` 内容不同，因此分别保留。新增需求以独立模块目录实现，仅将真正通用的计算抽取到 `common/`，入口自动加载必要依赖。
