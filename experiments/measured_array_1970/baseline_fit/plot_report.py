"""Reproduce all physical-versus-measured comparison figures."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'common'))
from model import geometry, great_circle_cut, load_data, split

OUT=Path(__file__).resolve().parent/"results"


def read_csv(name):
    with (OUT/name).open() as f:
        return list(csv.DictReader(f))


def save(fig,name):
    fig.savefig(OUT/name,dpi=140,bbox_inches="tight")
    plt.close(fig)


def main():
    records=load_data();lookup={(r["pitch"],r["azimuth"]):i for i,r in enumerate(records)}
    block=json.loads((OUT/"physical_config.json").read_text())
    ellipse=json.loads((OUT/"ellipse_config.json").read_text())
    eb=np.load(OUT/"ellipse_predicted_patterns.npz")["gains"]
    rb=np.load(OUT/"predicted_patterns.npz")["gains"]
    em=read_csv("ellipse_per_beam_metrics.csv")
    ranking=read_csv("ellipse_candidate_ranking.csv")
    fig,axes=plt.subplots(1,3,figsize=(15,5),constrained_layout=True)
    for ax,config,title in zip(axes[:2],(block,ellipse),("Full triangular block","Elliptically clipped triangular lattice")):
        p,_=geometry(config)
        ax.scatter(p[0],p[1],s=9)
        ax.set(xlabel="x [m]",ylabel="y [m]",title=f"{title}\n{p.shape[1]} active RX elements; one port/element")
        ax.set_aspect("equal")
    axes[2].scatter([int(r["sites"]) for r in ranking],[float(r["training_rmse_db"]) for r in ranking],s=12)
    axes[2].set(xlabel="Active elements per RX/TX aperture",ylabel="Subset training RMSE [dB]",title="Discrete ellipse search (not confidence bounds)",ylim=(0,2))
    save(fig,"geometry_search.png")

    examples=[(0,0),(10,60),(20,30),(30,150),(40,90),(50,330)]
    fig,axes=plt.subplots(6,2,figsize=(14,21),constrained_layout=True)
    for row,command in enumerate(examples):
        index=lookup[command];r=records[index];g=r["gain"]-r["gain"].max()
        for col,tangent in enumerate(("radial","tangential")):
            ax=axes[row,col]
            for gain,label,color in ((g,"Measured","black"),(eb[index],"502-element ellipse","tab:blue"),(rb[index],"Full block","tab:orange")):
                x,y=great_circle_cut(gain,*command,tangent)
                ax.plot(x,y,label=label,color=color,lw=1.1)
            ax.axhline(-3,color="gray",ls=":",lw=.7)
            ax.set(xlabel="Great-circle offset from command [deg]",ylabel="Normalized gain [dB]",ylim=(-65,1),xlim=(-15,15),
                title=f"Pitch={command[0]}, Azimuth={command[1]} | {split(r)} | {tangent}")
            ax.grid(alpha=.2)
            if row==0:
                ax.legend(fontsize=8)
    save(fig,"cut_comparisons.png")

    fig,axes=plt.subplots(2,3,figsize=(16,9),constrained_layout=True)
    for row,command in enumerate(((20,30),(40,90))):
        index=lookup[command];g=records[index]["gain"];g=g-g.max()
        for col,(values,title) in enumerate(((g,"Measured"),(eb[index],"502-element ellipse"),(eb[index]-g,"Model - measured"))):
            kw=dict(vmin=-60,vmax=0,cmap="viridis") if col<2 else dict(vmin=-15,vmax=15,cmap="RdBu_r")
            im=axes[row,col].imshow(values,origin="lower",aspect="auto",extent=(0,359,0,90),**kw)
            axes[row,col].set(xlabel="Fixed-frame phi [deg]",ylabel="Theta from array normal [deg]",title=f"P{command[0]} A{command[1]} | {title}")
            fig.colorbar(im,ax=axes[row,col],label="dB")
    save(fig,"map_comparisons.png")

    fig,axes=plt.subplots(1,3,figsize=(16,5),constrained_layout=True)
    for ax,key,title,limit in zip(axes,("high_gain_rmse_db","above_minus20_rmse_db","full_clip65_rmse_db"),
        ("Measured >= -10 dB","Measured >= -20 dB","All samples, clipped at -65 dB"),(1,4,15)):
        values=np.zeros((6,12))
        for r in em:
            if int(r["azimuth"])<360:
                values[int(r["pitch"])//10,int(r["azimuth"])//30]=float(r[key])
        im=ax.imshow(values,origin="lower",aspect="auto",vmin=0,vmax=limit)
        ax.set_xticks(range(12),range(0,360,30),rotation=90)
        ax.set_yticks(range(6),range(0,60,10))
        for p in range(6):
            for a in range(12):
                label="T" if (p==0 or (p in (1,3,5) and a%2==0)) else "H"
                ax.text(a,p,label,ha="center",va="center",color="white",fontsize=7)
        ax.set(xlabel="Command azimuth [deg]",ylabel="Command pitch [deg]",title=title+"\nT=train, H=held out")
        fig.colorbar(im,ax=ax,label="RMSE [dB]")
    save(fig,"error_by_command.png")

    fig,axes=plt.subplots(1,2,figsize=(11,5),constrained_layout=True)
    held=[r for r in em if r["split"].startswith("holdout")]
    stats={}
    for ax,tangent in zip(axes,("radial","tangential")):
        x=np.array([float(r[f"measured_{tangent}_width_deg"]) for r in held])
        y=np.array([float(r[f"model_{tangent}_width_deg"]) for r in held])
        stats[f"{tangent}_beamwidth_mae_deg"]=float(np.mean(abs(x-y)))
        ax.scatter(x,y,s=18);lim=(min(x.min(),y.min())-.3,max(x.max(),y.max())+.3)
        ax.plot(lim,lim,color="black",ls=":")
        ax.set(xlabel="Measured 3 dB width [deg]",ylabel="Model 3 dB width [deg]",title=f"Held-out {tangent} cuts",xlim=lim,ylim=lim)
    stats["mean_peak_error_deg"]=float(np.mean([float(r["peak_error_deg"]) for r in held]))
    (OUT/"beamwidth_summary.json").write_text(json.dumps(stats,indent=2)+"\n")
    save(fig,"beamwidth_comparison.png")

    fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
    for col,(title,threshold) in enumerate((("High gain (measured >= -10 dB)",-10),("Full pattern (both clipped at -65 dB)",-1e9))):
        for values,label in ((eb,"Ellipse"),(rb,"Full block")):
            errors=[]
            for i,r in enumerate(records):
                if not split(r).startswith("holdout"):
                    continue
                g=r["gain"]-r["gain"].max();mask=g[1:]>=threshold
                errors.extend(abs(np.maximum(values[i,1:],-65)-np.maximum(g[1:],-65))[mask])
            data=np.sort(errors)
            axes[col].plot(data,np.linspace(0,1,len(data)),label=label)
        axes[col].set(xlabel="Absolute gain error [dB]",ylabel="Fraction of held-out samples",title=title,xlim=(0,2 if col==0 else 40))
        axes[col].legend();axes[col].grid(alpha=.2)
    save(fig,"heldout_error_cdf.png")
    print("Saved six comparison figures and beamwidth_summary.json")


if __name__=="__main__":
    main()
