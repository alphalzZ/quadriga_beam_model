# 1.97 GHz 实测扫描方向图建模

本目录把 `../../data/measured_patterns/` 的 78 份文件作为**同一副天线的扫描状态**，搜索三角栅格阵列，并与全部方向图对比。详细结果见 [建模报告](results/REPORT.md)。

当前推荐候选为 **41×17 候选外框、椭圆裁剪后 502 个接收阵子**。这不是 41×17 的满阵，也不能据此确认真实硬件结构。按用户澄清，`2×500` 是发射、接收各约 500 个阵子；**上行接收端约 500 个阵子，不因圆极化选项翻倍**。

支持为一套阵列选择 LHCP 或 RHCP 理想极化。数据没有复电场、极化分量和绝对增益，所以这三者均为模型假设，未由实测验证。

## QuaDRiGa 使用

```matlab
addpath('/path/to/quadriga_beam_model/experiments/measured_array_1970');
addpath(genpath('/path/to/QuaDRiGa/quadriga_src'));

% 一个扫描状态下的接收阵列：502 个物理阵子，合并到 1 个 RHCP 信号端口。
[a, meta] = build_physical_array1970(31);

% 保留 502 个独立阵子端口，自行做时域/频域合并。
[a_elements, meta] = build_physical_array1970([]);

% 选择 LHCP；不会把阵子数量加倍。
a_lhcp = build_physical_array1970(31, [], 1, 'LHCP');

% 使用同一副天线另一种扫描状态。
M = load(fullfile(fileparts(which('build_physical_array1970')), ...
    'results', 'ellipse_physical_model.mat'));
disp(M.commands); % 每一行是 [Pitch, Azimuth]，对应 beam_index
```

将 `a` 赋给 `qd_layout.rx_array` 可用于上行接收，将另行构造的阵列赋给 `tx_array` 可用于发射。`coupling=conj(W_receive)`，QuaDRiGa 在接收端使用 coupling 的转置，因此等价于 `W_receive' * antenna_samples`，不是再次取错共轭。

阵面位于局部 x-y 平面，法向为 +z；观察角映射为 QuaDRiGa `azimuth=phi`、`elevation=90-theta`。`commands` 中的 Pitch 是离阵面法向的扫描角，不直接等同于 QuaDRiGa elevation。实际安装朝向应通过布局/天线姿态另行设置。

前述物理模型仅近似实测旁瓣。若要在已有扫描状态上保留原始功率方向图，可使用：

```matlab
% 选择实测 Pitch=20、Azimuth=30 的整副天线等效波束，1 个 RHCP 端口。
[beam, meta] = build_measured_beam1970(20, 30);
beam_lhcp = build_measured_beam1970(20, 30, 0, 'LHCP');
```

它逐网格保留实测功率，但不是 502 阵子的几何识别结果；只表示一个选中的扫描状态，不能把 78 个状态当作 78 个物理阵子。相位和极化关系仍是合成的，默认 0 dB 标度没有绝对增益校准意义。未测后半球填零，不表示真实后向辐射为零。

## 完整复现

需要 Python 3.12、NumPy 2.3.5、SciPy 1.17.0、Matplotlib 3.10.8，以及 MATLAB/Octave 和 QuaDRiGa 2.8.1。QuaDRiGa 源码版本为 `277866650eb115adb5b3e8ac252b0d1df073596d`；本次验证使用 GNU Octave 9.4.0。

从本目录运行：

```bash
python -m pip install -r requirements.txt
# 控制并行数，避免小矩阵运算因 BLAS 线程过多变慢。
export OPENBLAS_NUM_THREADS=2
python fit_model.py       # 满外框三角阵、矩形对照、审计及留出比较
python fit_ellipse.py     # 椭圆裁剪三角栅格搜索、精化及全方向图比较
python verify1970.py      # 数据、几何、权重、解析 Jacobian 与缓存检查
```

随后在 MATLAB/Octave 中运行：

```matlab
report = verify_quadriga1970('/path/to/QuaDRiGa');
```

再运行 `python plot_report.py` 重新生成对比图与波束宽度摘要。`results/REPORT.md` 记录本次建模的参数和结论；若改变搜索范围、数据或模型族，需结合新生成的 CSV/JSON 更新报告。

## 文件说明

| 文件 | 作用 |
|---|---|
| `model.py` | 数据解析、坐标、阵列因子、几何、指标与大圆切片 |
| `fit_model.py` | 完整外框三角阵的训练搜索和矩形基线 |
| `fit_ellipse.py` | 椭圆裁剪三角阵的离散排布与连续参数搜索 |
| `build_physical_array1970.m` | 推荐的 502 阵子 QuaDRiGa 阵列 |
| `build_measured_beam1970.m` | 一种实测扫描状态的功率等效模型 |
| `verify1970.py` | 数学与数据检查 |
| `verify_quadriga1970.m` | 全部实测功率导入、物理场和上行合并检查 |
| `results/ellipse_physical_model.mat` | 推荐候选的阵子坐标、幅度加权、72 个扫描命令及接收权重 |
| `results/physical_model.mat` | 满外框三角阵对照模型 |
| `results/*per_beam_metrics.csv` | 每份文件的主瓣、较强旁瓣、全域误差和宽度 |
| `results/*predicted_patterns.npz` | 全部 78 个扫描状态的预测功率图，顺序见 pitch/azimuth 数组 |
| `results/data_audit.csv` | 原始文件 SHA-256、实测峰值、坐标假设和训练/留出归属 |
| `results/*verification.*` | 当前运行的验证结果 |

几何参数、幅度加权和单元方向图由训练集共享估计；每个扫描状态仅改变相位权重。没有把实测功率开平方后冒充阵子级复信道，也没有把实测波束重放当作排布识别精度。
