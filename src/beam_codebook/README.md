# 矩形阵波束码本

3.5 GHz、纵向 32/横向 16、半波长间距的矩形阵示例。公共计算函数位于 `../common/`，由入口自动加载；全部缓存、CSV 和图片位于本模块 `results/`。

从仓库根目录运行：

```bash
octave --quiet --eval "run('src/beam_codebook/run_beam_model.m');"
BEAM_ANGLE_SOURCE=sheet octave --quiet --eval "run('src/beam_codebook/run_beam_model.m');"
python src/beam_codebook/plot_beam_model.py photo
python src/beam_codebook/plot_beam_model.py sheet
```

默认 `photo`；`sheet` 保留用户给出的原始 14 组 `ele_tab/azi_tab`。本历史矩形阵示例将 ele 作为水平面下倾角，其权重及总增益属于该旧坐标模型。实测天线的 sheet 时域链路**仅读取该缓存的角度表与名称**，并按用户确认的固定法线 Pitch/Azimuth 重新构造 1.97 GHz、502 阵子权重和校准波束，不复用旧权重、旧几何或旧增益。
