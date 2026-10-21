# 三角栅格与圆极化波束模型

`run_beam_model_triangular_cp.m` 参考 `run_beam_model.m` 的 14 个波束、全角域扫描、增益与串扰计算流程，增加可选阵列栅格和圆极化端口模型。无需 QuaDRiGa 即可生成波束模型；将结果转换为 QuaDRiGa 天线、运行完整验证时需要 QuaDRiGa 2.8.1。

## 运行与配置

```matlab
addpath('/path/to/quadriga_beam_model/src/beam_model_triangular_cp');

% 默认：3.5 GHz，纵向 32 行、横向 16 列，三角栅格、RHCP。
model = run_beam_model_triangular_cp();

% 左旋圆极化。
model = run_beam_model_triangular_cp(struct('polarization','LHCP'));

% 同一空间位置上设置 LHCP/RHCP 两个独立端口。
model = run_beam_model_triangular_cp(struct('polarization','DUAL'));

% 矩形栅格对照。
model = run_beam_model_triangular_cp(struct( ...
    'lattice','rectangular','polarization','RHCP'));
```

| 参数 | 默认值 | 含义 |
|---|---|---|
| `fc` | `3.5e9` | 载频，Hz |
| `Mv`、`Mh` | `32`、`16` | 行数与列数；三角栅格的相邻行横向错开 |
| `spacing_lambda` | `0.5` | 最近邻阵子间距，以波长为单位 |
| `lattice` | `'triangular'` | `'triangular'` 或 `'rectangular'` |
| `polarization` | `'RHCP'` | `'LHCP'`、`'RHCP'` 或 `'DUAL'` |
| `angle_source` | `'photo'` | `'photo'` 或原表格的 `'sheet'` 方位角 |
| `grid_step_deg` | `1` | 方向图扫描步长，需整除 180 且不超过 10° |
| `output_dir` | 本模块 `results/` | 输出位置；`''` 表示仅返回结果、不写文件 |

默认输出前缀为 `beam_model_triangular_rhcp`，包括 MAT 缓存、波束统计 `_summary.csv`、方位切片 `_cuts.csv`、串扰 `_crosstalk.csv` 和空间坐标 `_sites.csv`。其他栅格、极化使用独立前缀。不同频率、行列数或角度来源使用相同栅格/极化时，应指定不同 `output_dir` 以保留多份结果。

绘图：

```bash
python plot_beam_model_cp.py /path/to/beam_model_triangular_rhcp
```

输出 `_overview.png`，包含空间排布、目标方向、串扰矩阵及方位增益切片。图中的空间点是阵子位置数；双极化在相同位置有两个端口，不会画成两套分开的阵子。

## 几何与极化模型

阵列位于局部 `y-z` 平面，法向为 `+x`。三角栅格在中心化前的坐标为：

```text
y(row,col) = d * (col + 0.5*mod(row,2))
z(row,col) = sqrt(3)/2 * d * row
```

因此同一行内的间距为 `d`，相邻行高度差为 `sqrt(3)/2*d`，相邻行横向错开 `d/2`，最近邻距离仍为 `d`。默认 3.5 GHz 下 `d≈0.0428275 m`。相同阵子数与最近邻间距下，三角阵列和矩形阵列的物理孔径不同，不能将其差异全部解释为栅格的性能优势。

每个有效圆极化端口使用同一个 3GPP 标量单元功率方向图 `elem_gain_dB`，将电场分成两个正交、等幅、相差 90° 的分量。在 QuaDRiGa 的局部 `[Fa;Fb]` 基底与官方馈电符号约定下：

```text
LHCP Jones = [1; +j] / sqrt(2)
RHCP Jones = [1; -j] / sqrt(2)
```

端口电场等于 `sqrt(element_gain_linear) * Jones`。两个分量的总功率等于原单元增益，未额外增加 3 dB。这是合成的理想圆极化方向图：没有实测轴比、幅相误差、互耦和扫描引起的极化劣化，也不是交叉偶极子在全角域的真实方向图。

| 配置 | 空间位置数 `Nsite` | 端口数 `N` | 波束数 `B` | `W` 维度 |
|---|---:|---:|---:|---|
| LHCP 或 RHCP | 512 | 512 | 14 | 512 × 14 |
| DUAL | 512 | 1024 | 28 | 1024 × 28 |

双极化端口顺序：前 512 个为 LHCP，后 512 个为 RHCP；波束顺序同样先是 14 个 LHCP，再是 14 个 RHCP。每个波束仅使用其对应极化的 512 个端口，权重范数为 1，目标方向的阵列因子峰值为 512，未将两种极化的端口数错误地算成 1024 的相干增益。

主要缓存变量：

- `p_site`：独立空间位置，3 × `Nsite`。
- `p`：每个端口的位置，3 × `N`；双极化位置重复。
- `port_site`、`port_pol`、`port_jones`：端口对应的位置、极化编号和 Jones 向量。
- `beam_pol`、`beam_jones`、`beam_name`：每个波束的极化与标签。
- `W`、`U`：接收权重及波束方向。
- `Gx`、`XT`：考虑极化匹配的阵列因子和归一化串扰。

串扰计算使用共同局部横向基底中的接收匹配关系：

```text
A(port,j) = exp(j*k*p(:,port)'*U(:,j)) ...
            * (port_jones(:,port)' * beam_jones(:,j))
Gx = abs(W' * A).^2
```

理想相反圆极化在这个基底中的 Jones 内积为零，因此相反极化块的串扰为 `-Inf dB`。这表示理想单元的极化正交性，不表示卫星传播链路始终具有无限隔离。真实传播还涉及收发坐标转换、观察方向、反射和去极化；QuaDRiGa 信道生成器会处理其传播模型中的几何与极化效应。

## 与 QuaDRiGa 对接

```matlab
addpath(genpath('/path/to/QuaDRiGa/quadriga_src'));
model = load('beam_model_triangular_rhcp.mat');
array = make_cp_arrayant(model);  % 默认使用 5° 的方向图采样
l = qd_layout;
l.simpar.center_frequency = model.fc;
l.tx_array = array;
```

`make_cp_arrayant` 将每个有效 CP 端口存为一个阵列元素，设置其真实位置和复数 `Fa/Fb`，并采用单位耦合矩阵；无额外的双端口馈电层。因此 `array` 的信道端口顺序与 `model.W` 完全一致。

适配器的 5° 方向图采样独立于波束模型的 1° 扫描。QuaDRiGa 在采样点之间插值，可使用 `make_cp_arrayant(model,1)` 提高分辨率，但全尺寸双极化方向图的内存需求会显著增加。

接入时域链路时，需要读取这份新缓存，用 `make_cp_arrayant(model)` 创建天线，并根据 `N`、`B` 和用户配置设置维度；原 `iface_demo_time_domain` 中针对矩形单极化缓存的检查不适用于三角/双圆极化模型，不能仅替换 MAT 文件名。

## 验证

```matlab
report = verify_beam_model_triangular_cp('/path/to/QuaDRiGa');
```

已在 GNU Octave 9.4.0、QuaDRiGa 2.8.1 上通过：

1. 全尺寸三角栅格的行错位、行列间距、最近邻距离与几何中心检查。
2. LHCP/RHCP 权重归一化、匹配增益一致性、双圆极化交叉块零响应检查。
3. 矩形模式的 `p`、`W` 与现有 3.5 GHz 矩形缓存回归对照。
4. 全尺寸 1024 端口 QuaDRiGa 复电场、空间投影距离及总单元增益检查。
5. 与 QuaDRiGa 官方 H/V 模板及圆极化馈电矩阵的正前方场分量对照。
6. 4×3 三角阵列、24 个双圆极化端口与双圆极化接收天线的真实 QuaDRiGa LOS 信道生成，检查信道维度及有限非零系数。
7. MAT/CSV 写出与 MAT 重载检查。

仓库附带默认完整尺寸 RHCP 的生成结果与概览图。验证中的实际信道为小阵列集成检查，不是完整卫星场景、全尺寸链路 BER 或实测天线性能验证。
