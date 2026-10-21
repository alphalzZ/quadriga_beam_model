# 1.97 GHz 实测方向图建模

本项目使用同一副天线的 78 个扫描状态；上行接收端为约 500 个阵子，当前拟合候选为 502 个三角栅格阵子。

| 模块 | 职责 | 入口 |
|---|---|---|
| [common](common/) | 数据读取、几何与指标、QuaDRiGa 物理阵列/实测等效波束适配 | `model.py`、`build_physical_array1970.m` |
| [baseline_fit](baseline_fit/README.md) | 早期物理排布搜索、训练/留出比较 | `fit_model.py`、`fit_ellipse.py` |
| [all_data_fit](all_data_fit/README.md) | 全部 78 份数据和全部采样点拟合 | `fit_all_data.py`、`build_allfit_beam1970.m` |
| [beam_patterns](beam_patterns/README.md) | 14 个 sheet 拟合波束极坐标切面与 PDF | `export_sheet_pattern_cuts.m`、`plot_sheet_polar.py` |
| [time_domain](time_domain/README.md) | sheet 离轴 Pitch/Azimuth 驱动的时域 OFDM 链路 | `iface_allfit_time_domain.m` |

原始方向图保存在 [data/measured_patterns](../../data/measured_patterns/)，数据说明见[波束方向图数据说明](../../data/measured_patterns/beam_pattern_data_notes.md)。各功能模块的结果仅写入自己的 `results/`；Python 依赖集中在本目录 `requirements.txt`。新增功能应新建同级模块，共享计算函数放入 `common/`，不依赖其他模块的当前工作目录。
