# 时域信道与波束合并示例

入口是 `iface_demo_time_domain.m`，采用带配置参数的 MATLAB/Octave 函数形式。

```matlab
addpath('/path/to/quadriga_beam_model/src/beam_model_photo');
result = iface_demo_time_domain(struct( ...
    'qd_path', '/path/to/QuaDRiGa', ...
    'nSlot', 4, ...
    'output_file', '/tmp/iface_time_domain_out.mat'));
```

也可以设置环境变量 `QD_PATH`，或提前将 QuaDRiGa 加入 MATLAB/Octave 路径，然后直接运行 `result = iface_demo_time_domain();`。波束文件始终从函数所在目录读取；结果默认保存到当前工作目录的 `iface_time_domain_out.mat`。无需修改原来的 `iface_demo.m`。

本云环境可运行：

```bash
cd /workspace/.cloud-setup/validation
/workspace/.cloud-setup/bin/octave --quiet --eval "addpath('/workspace/quadriga_beam_model/src/beam_model_photo'); iface_demo_time_domain(struct('qd_path','/workspace/.cloud-setup/deps/QuaDRiGa'));"
```

## 链路处理

1. 为每个用户生成 QPSK 符号，映射到零频两侧的 288 个子载波，DC 留空。
2. IFFT 生成时域采样并添加循环前缀。
3. 将 QuaDRiGa 的多径系数、时延转换为每个用户到每个天线端口的离散 FIR。小数时延采用加窗 sinc 插值，默认每条路径使用 33 个插值抽头。
4. 在时域对每个天线端口执行线性卷积，并叠加各用户信号及天线噪声。
5. 对每一个接收采样执行 `y_beam = Wq' * r_ant`，完成时域波束合并。
6. 接收端去 CP、FFT，使用理想 CSI 做逐子载波标量 MMSE 均衡和硬判决，统计 BER。

传播和波束合并过程对应：

```text
r_ant(a,n) = sum_u sum_t h_ant{u}(a,t+1) * tx(u,n-t) + noise(a,n)
y_beam(k,n) = sum_a conj(Wq(a,k)) * r_ant(a,n)
```

`h_ant{u}` 为 `512 × L`，`Wq` 为 `512 × 14`，`r_ant` 为 `512 × 每时隙采样数`，`y_beam` 为 `14 × 每时隙采样数`。权重是随时间不变的复数系数，直接逐采样作用，无需对权重做 IFFT。当前模型使用一个静态信道快照，不模拟多普勒或时变波束。

## 与原示例的区别及边界

- 使用与 `run_beam_model.m` 相同的配置：3.5 GHz、纵向 32/横向 16，共 512 个端口。启动时检查波束缓存的载频、阵列尺寸、端口位置及映射唯一性；缓存不一致时提示重新运行生成脚本。现有 `beam_model_photo.mat` 及配套 CSV、PNG 已按此配置重新生成。
- 默认 `nFFT=512`、子载波间隔 120 kHz，采样率为 61.44 MHz。CP 根据离散 FIR 长度自动设置，保证线性卷积后的 OFDM 符号不受前一符号干扰；这不是固定的 3GPP 标准 CP 配置。
- 每个用户移除一个公共最小时延，模拟理想定时提前；保留多径超额时延及天线间相对时延。插值 FIR 另引入 `fd_half` 个采样的公共因果延迟，包含在 CP 和接收 CSI 中。
- 使用有限长度 sinc 近似连续时延，其误差随插值长度和有效频带而变。`Heff` 由实际使用的离散 FIR 计算，保证均衡器与该离散信道一致；它不表示与连续时延模型完全相同。
- 采用单位酉 IFFT/FFT，因此归一化时域噪声方差与每个子载波上的复噪声方差均为 `N0`。噪声在天线端口加入，保留非正交波束之间的噪声相关性。
- FFT 得到的信道仅用于接收 CSI、SINR 预测和无噪声参考检查。发射波形的传播使用时域线性卷积，波束合并发生在去 CP 和 FFT 之前。
- 同原示例，利用信道互易性将 QuaDRiGa 生成的基站到用户端口系数用于上行；端口转置不取共轭。
- BER 预测使用高斯干扰近似；有限用户 QPSK 干扰及有限统计样本会使测量值与预测存在差异。

可配置参数：`qd_path`、`nSlot`（默认 4）、`nSym`（默认 14）、`nFFT`（默认 512）、`fd_half`（默认 16）、`seed`（默认 1）、`output_file`。将 `output_file` 设为空字符串可关闭保存。时域仿真比原来的频域链路耗时更多，可先用 `nSlot=1,nSym=2` 检查运行，再增加时隙数。

结果中保存 FIR、权重、接收 CSI、SINR、BER、误比特数和总比特数，以及最后一个时隙的发送采样 `tx_last` 和合并后接收采样 `y_beam_last`。`consistency_error` 记录时域链路与相同离散 FIR 的频域参考之间的相对误差，代码要求其小于 `1e-10`。
