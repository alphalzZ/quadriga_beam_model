# 原始频域 OFDM 链路

入口 `iface_demo.m` 使用 `../beam_codebook/results/beam_model_photo.mat`，3.5 GHz、32×16 阵子。QuaDRiGa 路径通过环境变量提供，结果保存到本模块 `results/iface_out.mat`。

```bash
QD_PATH=/path/to/QuaDRiGa octave --quiet --eval "run('src/link_frequency_domain/iface_demo.m');"
```

本模块保留原始频域传播方案；逐采样传播和合并的版本见 [link_time_domain](../link_time_domain/README.md)，全数据拟合天线版本见 [sheet 时域链路](../../experiments/measured_array_1970/time_domain/README.md)。

目录迁移时同步修正了旧脚本写死的 QuaDRiGa 路径、依赖当前目录的缓存/结果路径，以及与现有缓存不一致的 16×32 阵列排列。
