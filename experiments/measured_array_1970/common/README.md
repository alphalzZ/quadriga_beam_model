# 实测天线共享模型

`model.py` 负责读取原始数据、坐标、几何、阵列因子和比较指标。`build_physical_array1970.m` 构造独立阵子端口或物理合并阵列，默认读取 `../baseline_fit/results/ellipse_physical_model.mat`，也可显式传入全数据拟合模型。`build_measured_beam1970.m` 将一个已测扫描状态构造成等效波束端口。

本目录仅保存共享计算与适配代码；基线搜索、全数据拟合和时域链路各自管理配置与 `results/`。
