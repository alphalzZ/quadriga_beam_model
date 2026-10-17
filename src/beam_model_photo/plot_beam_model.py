import csv, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

tag = sys.argv[1] if len(sys.argv) > 1 else "sheet"
rows = list(csv.DictReader(open(f"beam_summary_{tag}.csv")))
names = [r["beam"] for r in rows]
ele = np.array([float(r["ele_deg"]) for r in rows])          # downtilt
azi = np.array([float(r["azi_deg"]) for r in rows])
dist = np.array([float(r["dist_m"]) for r in rows])
gcmd = np.array([float(r["gain_at_cmd_dBi"]) for r in rows])
gmax = np.array([float(r["gmax_dBi"]) for r in rows])
mirror = np.array([float(r["mirror_az_deg"]) for r in rows])
cuts = np.loadtxt(f"beam_cuts_az_{tag}.csv", delimiter=",")         # K x 361, total gain dBi
xt = np.loadtxt(f"beam_crosstalk_{tag}.csv", delimiter=",")         # K x K
az_grid = np.arange(-180, 181, 1.0)
wrap = lambda a: (a + 180) % 360 - 180
az_w = wrap(azi)

fig = plt.figure(figsize=(16, 13))
gs = GridSpec(4, 7, figure=fig, height_ratios=[1.25, 1.25, 1, 1], hspace=0.55, wspace=0.35)

# (a) direction map: solid = commanded, hollow = front-back twin
ax = fig.add_subplot(gs[0:2, 0:3])
sc = ax.scatter(az_w, -ele, c=np.log10(dist), cmap="viridis", s=180, edgecolor="k", zorder=3)
back = np.abs(az_w) > 90
if back.any():
    ax.scatter(wrap(mirror[back]), -ele[back], facecolor="none", edgecolor="crimson", s=180, linewidth=1.8, zorder=2,
               label="front-back twin of back-facing beams")
for i, n in enumerate(names):
    ax.annotate(n, (az_w[i], -ele[i]), xytext=(6, 5), textcoords="offset points", fontsize=10)
    if back[i]:
        ax.plot([az_w[i], wrap(mirror[i])], [-ele[i], -ele[i]], ls=":", color="crimson", lw=1)
ax.axvspan(-90, 90, color="0.93", zorder=0)
ax.text(0, -57, "front hemisphere (|az|<90)", ha="center", color="0.4")
ax.set_xlim(-180, 180); ax.set_ylim(-60, 0)
ax.set_xlabel("azimuth [deg]"); ax.set_ylabel("elevation = -downtilt [deg]")
ax.set_title("(a) 14 beams: commanded direction vs front-back twin")
if back.any():
    ax.legend(loc="upper left", fontsize=9)
cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
cb.set_label("log10(implied distance / m)")

# (b) crosstalk matrix, total gain relative to own target
ax = fig.add_subplot(gs[0:2, 3:7])
im = ax.imshow(np.clip(xt, -40, 25), cmap="RdBu_r", vmin=-40, vmax=25)
ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=90, fontsize=8)
ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=8)
ax.set_xlabel("interfering direction j"); ax.set_ylabel("beam k")
for k in range(len(names)):
    for j in range(len(names)):
        if k != j and xt[k, j] > -10:
            ax.text(j, k, f"{xt[k,j]:.0f}", ha="center", va="center", fontsize=6.5)
ax.set_title("(b) crosstalk XT(k,j) [dB]: UE at dir j seen by beam k, relative to own target\n(>0 dB = interferer stronger than the target)")
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("dB")

# (c) azimuth cuts of total gain at each beam's elevation
for k in range(len(names)):
    r, c = divmod(k, 7)
    ax = fig.add_subplot(gs[2 + r, c])
    ax.plot(az_grid, cuts[k], lw=1.1)
    ax.axvline(az_w[k], color="g", ls="--", lw=1, label="commanded")
    if abs(az_w[k]) > 90:
        ax.axvline(wrap(mirror[k]), color="crimson", ls=":", lw=1, label="twin")
    ax.set_title(f"{names[k]}  el={-ele[k]:.1f}", fontsize=9)
    ax.set_xlim(-180, 180); ax.set_ylim(-5, 40)
    ax.set_xticks([-180, -90, 0, 90, 180]); ax.tick_params(labelsize=7)
    ax.grid(alpha=0.3)
    if k == 0:
        ax.legend(fontsize=6, loc="lower center")
fig.text(0.5, 0.005, "(c) total gain [dBi] vs azimuth at each beam's elevation (array factor x 3GPP element pattern)",
         ha="center", fontsize=11)
title_src = "spreadsheet azimuths" if tag == "sheet" else "azimuths corrected from panorama (120 deg FOV, centre = boresight)"
fig.suptitle(f"512-element (16x32, lambda/2) receive beam model, 28 GHz, {title_src}", fontsize=13, y=0.995)
fig.savefig(f"beam_model_overview_{tag}.png", dpi=120, bbox_inches="tight")
print(f"saved beam_model_overview_{tag}.png")
