"""Plot 14 calibrated sheet beams, using the exported QuaDRiGa field cuts."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.io import loadmat

OUT = Path(__file__).resolve().parent / "results"
FLOOR = -60
families = {f.name for f in font_manager.fontManager.ttflist}
font = next((name for name in ("Noto Sans CJK SC", "Noto Sans CJK JP", "Microsoft YaHei", "SimHei")
             if name in families), "DejaVu Sans")
plt.rcParams.update({"font.family": font, "axes.unicode_minus": False,
                     "font.size": 10, "pdf.fonttype": 42})


def polar(ax, angle, gain, color, title, command, theta_cut=False):
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    if theta_cut:
        ax.set_thetamin(0)
        ax.set_thetamax(90)
        ax.set_xticks(np.deg2rad([0, 15, 30, 45, 60, 75, 90]))
    else:
        ax.set_xticks(np.deg2rad(np.arange(0, 360, 45)))
    ax.set_ylim(0, -FLOOR)
    ticks = np.arange(10, 61, 10)
    ax.set_yticks(ticks, [f"{t + FLOOR}" for t in ticks])
    ax.set_rlabel_position(225 if not theta_cut else 0)
    ax.tick_params(labelsize=8)
    ax.grid(color="#dce3ec", linewidth=0.6)
    ax.spines["polar"].set_color("#a4afbd")
    ax.plot(np.deg2rad(angle), np.clip(gain, FLOOR, 0) - FLOOR,
            color=color, linewidth=1.1)
    ax.plot([np.deg2rad(command)] * 2, [0, -FLOOR],
            "--", color="#e78a23", linewidth=0.9, alpha=0.9)
    ax.set_title(title, fontsize=11, pad=18, color="#203047")


def main():
    d = loadmat(OUT / "sheet_pattern_cuts.mat", simplify_cells=True)
    names = list(d["beam_name"])
    theta, phi = d["theta_deg"], d["phi_deg"]
    tg, pg = d["theta_normalized_db"], d["phi_normalized_db"]
    assert tg.shape[0] == pg.shape[0] == len(names) == 14
    assert np.max(np.abs(pg[:, 0] - pg[:, -1])) < 1e-7
    assert np.max(tg) < 1e-10 and np.max(pg) < 1e-10
    extrap = np.asarray(d["extrapolated"], dtype=bool)
    summary = []
    footer = ("每个波束的两切面共用峰值归一化；外圈 0 dB，内圈 −60 dB，低于 −60 dB 仅显示截底。\n"
              "曲线来自时域链路使用的全数据拟合 QuaDRiGa 波束；复相位及圆极化为模型假设。橙色虚线为指定指向。")
    with PdfPages(OUT / "14_beams_polar_cuts.pdf") as pdf:
        for k, name in enumerate(names):
            pitch, azi = d["pitch_deg"][k], d["azimuth_deg"][k]
            fig, axes = plt.subplots(1, 2, figsize=(12, 6.3), subplot_kw={"projection": "polar"})
            polar(axes[0], theta, tg[k], "#1768b5",
                  f"离轴切面：固定 φ = {azi:.2f}°\nθ = 0°～90°", pitch, True)
            polar(axes[1], phi, pg[k], "#00836b",
                  f"方位切面：固定 θ = {pitch:.2f}°\nφ = 0°～360°（周期闭合）", azi)
            note = "  |  C9：50° 边界方向图外推 0.37°" if extrap[k] else ""
            fig.suptitle(f"{name}  |  Pitch = {pitch:.2f}°，Azimuth = {azi:.2f}°{note}\n"
                         "1.97 GHz · 502 阵子 · 三角栅格 · 全数据拟合 · RHCP", fontsize=13, y=.99)
            fig.text(.5, .025, footer, ha="center", fontsize=9, color="#586477")
            fig.subplots_adjust(left=.05, right=.95, bottom=.19, top=.78, wspace=.35)
            fig.savefig(OUT / f"{name}_polar_cuts.png", dpi=160)
            pdf.savefig(fig)
            plt.close(fig)
            summary.append({"beam": str(name), "pitch_deg": float(pitch), "azimuth_deg": float(azi),
                            "extrapolated": bool(extrap[k]), "normalization_db": float(d["normalization_db"][k]),
                            "theta_cut_peak_error_deg": float(d["peak_error_deg"][k, 0]),
                            "phi_cut_peak_angular_error_deg": float(d["peak_error_deg"][k, 1])})

    for mode, angles, gains, color, suffix in [
        ("离轴切面（各固定本波束 φ）", theta, tg, "#1768b5", "theta"),
        ("方位切面（各固定本波束 θ）", phi, pg, "#00836b", "phi"),
    ]:
        fig, axes = plt.subplots(4, 4, figsize=(16, 17), subplot_kw={"projection": "polar"})
        for k, ax in enumerate(axes.flat):
            if k >= len(names):
                ax.set_visible(False)
                continue
            command = d["pitch_deg"][k] if suffix == "theta" else d["azimuth_deg"][k]
            polar(ax, angles, gains[k], color,
                  f'{names[k]}  θ={d["pitch_deg"][k]:.2f}°  φ={d["azimuth_deg"][k]:.2f}°'
                  + (" *" if extrap[k] else ""), command, suffix == "theta")
        fig.suptitle(f"14 个 sheet 波束 · {mode}\n1.97 GHz / 502 阵子 / 全数据拟合", fontsize=18, y=.99)
        fig.text(.5, .015, "共同显示范围 −60～0 dB；各波束两切面共用同一归一化系数。* C9 外推 0.37°。", ha="center")
        fig.subplots_adjust(top=.89, bottom=.07, hspace=.62, wspace=.48)
        fig.savefig(OUT / f"14_beams_{suffix}_polar.png", dpi=160)
        plt.close(fig)

    for suffix, angles, gain, raw in [("theta", theta, tg, d["theta_gain_db"]),
                                      ("phi", phi, pg, d["phi_gain_db"])]:
        with (OUT / f"{suffix}_cuts.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["beam", "angle_deg", "model_gain_db", "normalized_gain_db"])
            for k, name in enumerate(names):
                writer.writerows((name, f"{a:.6f}", f"{r:.8f}", f"{g:.8f}")
                                 for a, r, g in zip(angles, raw[k], gain[k]))
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("PASS: 14 paired PNGs, 2 polar overviews, 14-page PDF, CSV cuts and summary")


if __name__ == "__main__":
    main()
