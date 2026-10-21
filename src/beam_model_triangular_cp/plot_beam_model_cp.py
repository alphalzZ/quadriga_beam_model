"""Plot CP model CSV exports: python plot_beam_model_cp.py /path/to/prefix."""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    prefix = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent/"results"/"beam_model_triangular_rhcp"
    with Path(f"{prefix}_summary.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    sites = np.atleast_2d(np.loadtxt(f"{prefix}_sites.csv", delimiter=","))
    cuts = np.atleast_2d(np.loadtxt(f"{prefix}_cuts.csv", delimiter=","))
    xt = np.atleast_2d(np.loadtxt(f"{prefix}_crosstalk.csv", delimiter=","))
    assert cuts.shape[0] == len(rows) and xt.shape == (len(rows), len(rows))
    names = [r["beam"] for r in rows]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    ax = axes[0, 0]
    ax.scatter(sites[:, 1], sites[:, 2], s=9)
    ax.set(xlabel="y [m]", ylabel="z [m]", title=f"Array geometry: {len(sites)} sites")
    ax.set_aspect("equal")
    ax = axes[0, 1]
    for pol, marker in (("LHCP", "o"), ("RHCP", "x")):
        group = [r for r in rows if r["polarization"] == pol]
        if group:
            ax.scatter([float(r["azi_deg"]) for r in group],
                       [-float(r["ele_deg"]) for r in group], marker=marker, label=pol)
    ax.set(xlabel="azimuth [deg]", ylabel="elevation [deg]", title="Commanded directions")
    ax.legend()
    ax = axes[1, 0]
    image = ax.imshow(np.clip(xt, -60, 10), vmin=-60, vmax=10, cmap="RdBu_r")
    ticks = range(len(names))
    ax.set_xticks(ticks, names, rotation=90, fontsize=6)
    ax.set_yticks(ticks, names, fontsize=6)
    polarizations = {r["polarization"] for r in rows}
    title = ("Crosstalk [dB], ideal opposite CP = -inf" if len(polarizations) > 1
             else f"Crosstalk [dB], matched {rows[0]['polarization']} illumination")
    ax.set(xlabel="Incident beam / polarization", ylabel="Receive beam", title=title)
    fig.colorbar(image, ax=ax, label="dB")
    ax = axes[1, 1]
    az = np.linspace(-180, 180, cuts.shape[1])
    # LHCP/RHCP have identical matched cuts in this ideal model.
    for index, row in enumerate(rows[:14]):
        ax.plot(az, cuts[index], lw=0.9, label=row["beam"].split("_")[0])
    ax.set(xlabel="azimuth [deg]", ylabel="Matched total gain [dBi]",
           title="Azimuth cuts at each beam's nearest elevation sample", ylim=(-35, 40))
    ax.legend(ncol=4, fontsize=7)
    fig.suptitle(prefix.name)
    target = Path(f"{prefix}_overview.png")
    fig.savefig(target, dpi=140)
    plt.close(fig)
    print(f"saved {target}")


if __name__ == "__main__":
    main()
