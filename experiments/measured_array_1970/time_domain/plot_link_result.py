"""Plot the verified default all-fit time-domain link results."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

OUT=Path(__file__).resolve().parent/"results"


def main():
    with (OUT/"beam_metrics.csv").open() as f:
        rows=list(csv.DictReader(f))
    report=json.loads((OUT/"validation.json").read_text())
    data=loadmat(OUT/"receiver_reference.mat")
    assert report["slots"]==4 and report["symbols_per_slot"]==14
    assert len(rows)==14 and data["Heff"].shape==(14,288,14)
    names=[r["beam"] for r in rows]
    x=np.arange(14)
    fig,axes=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
    axes[0].semilogy(x,[float(r["ber_pred"]) for r in rows],"o-",label="Gaussian-interference prediction")
    axes[0].semilogy(x,[float(r["ber_meas"]) for r in rows],"x--",label="Time-domain measured BER")
    axes[0].set_xticks(x,names,rotation=45)
    axes[0].set(xlabel="User / assigned beam",ylabel="BER",title="1.97 GHz, calibrated array, 4 x 14 OFDM symbols")
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    axes[1].bar(x,[float(r["mean_sinr_db"]) for r in rows])
    axes[1].set_xticks(x,names,rotation=45)
    axes[1].set(xlabel="User / assigned beam",ylabel="Mean SINR [dB]",title="Same transmit PSD and receiver noise as original demo")
    fig.savefig(OUT/"ber_and_sinr.png",dpi=150);plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
    for extrapolated,marker,label in ((0,"o","Interpolated measured scans"),(1,"x","Boundary-pattern extrapolation")):
        group=[r for r in rows if int(r["extrapolated"])==extrapolated]
        axes[0].scatter([float(r["fit_scan_phi_deg"]) for r in group],
                        [float(r["fit_scan_theta_deg"]) for r in group],marker=marker,label=label)
        for r in group:
            axes[0].annotate(r["beam"],(float(r["fit_scan_phi_deg"]),float(r["fit_scan_theta_deg"])),xytext=(3,4),textcoords="offset points",fontsize=8)
    axes[0].axhline(50,color="gray",ls=":")
    axes[0].set(xlabel="Canonical fitted scan azimuth [deg]",ylabel="Off-normal scan angle [deg]",title="Original world-frame directions are retained")
    axes[0].legend(fontsize=8)
    covariance=data["noise_covariance"]
    normalized=abs(covariance)/np.sqrt(np.outer(np.real(np.diag(covariance)),np.real(np.diag(covariance))))
    im=axes[1].imshow(normalized,vmin=0,vmax=1)
    axes[1].set_xticks(x,names,rotation=90,fontsize=7);axes[1].set_yticks(x,names,fontsize=7)
    axes[1].set(title="Retained cross-beam noise correlation",xlabel="Beam",ylabel="Beam")
    fig.colorbar(im,ax=axes[1])
    fig.savefig(OUT/"directions_and_noise.png",dpi=150);plt.close(fig)
    print("Saved default link BER/SINR, direction and noise-correlation plots")


if __name__=="__main__":main()
