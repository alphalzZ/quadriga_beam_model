"""Compare a physically sampled elliptical triangular aperture to full blocks.

This is a second predeclared FOOTPRINT family, motivated by training sidelobes.
The 42 held-out files are never used in parameter fitting or shape ranking.
"""
import json
from pathlib import Path

import numpy as np
from scipy.io import savemat
from scipy.optimize import least_squares

from fit_model import training_samples, write_csv
from model import C0, FC, direction, geometry, load_data, metrics, predict, split

OUT = Path(__file__).resolve().parent/"results"


class DiscreteObjective:
    def __init__(self, rows, cols, samples):
        cfg = dict(rows=rows, cols=cols, footprint="ellipse", parameters=[1,0,0,0,0])
        self.xy = geometry(cfg)[0][:2]*FC/C0
        m,n=np.meshgrid(np.arange(rows),np.arange(cols),indexing="ij")
        x=np.sqrt(3)/2*m; y=n+.5*(m%2);x=x-x.mean();y=y-y.mean()
        mask=(x/(np.sqrt(3)/2*(rows-1)/2))**2+(y/((cols-1)/2))**2<=1+1e-12
        self.cr=-np.cos(2*np.pi*m/(rows-1))[mask]
        self.cc=-np.cos(2*np.pi*n/(cols-1))[mask]
        self.duv,theta,pitch,self.target=samples
        self.eq=10*np.log10(np.maximum(np.cos(np.deg2rad(theta)),1e-12)/np.cos(np.deg2rad(pitch)))
        self.last=None

    def evaluate(self,v):
        if self.last is not None and np.array_equal(v,self.last):
            return self.residual,self.jacobian
        d,psi,ar,ac,q=v;c,s=np.cos(np.deg2rad(psi)),np.sin(np.deg2rad(psi))
        u=self.duv
        dot=self.xy[0,:,None]*(c*u[0]+s*u[1])+self.xy[1,:,None]*(-s*u[0]+c*u[1])
        dotpsi=(self.xy[0,:,None]*(-s*u[0]+c*u[1])+self.xy[1,:,None]*(-c*u[0]-s*u[1]))*np.pi/180
        e=np.exp(2j*np.pi*d*dot)
        w=(1+ar*self.cr)*(1+ac*self.cc);norm=w.sum()
        a=w@e/norm
        d_ar=self.cr*(1+ac*self.cc);d_ac=self.cc*(1+ar*self.cr)
        da=[(w[:,None]*e*(2j*np.pi*dot)).sum(axis=0)/norm,
            (w[:,None]*e*(2j*np.pi*d*dotpsi)).sum(axis=0)/norm,
            (d_ar@e-a*d_ar.sum())/norm,
            (d_ac@e-a*d_ac.sum())/norm]
        self.residual=20*np.log10(np.maximum(abs(a),1e-15))+q*self.eq-self.target
        self.jacobian=np.column_stack([20/np.log(10)*np.real(z/(a+1e-20)) for z in da]+[self.eq])
        self.last=v.copy()
        return self.residual,self.jacobian


def fit(rows,cols,samples,max_nfev=35,initial=None):
    objective=DiscreteObjective(rows,cols,samples)
    orientation=0 if np.std(objective.xy[0])>np.std(objective.xy[1]) else 90
    if initial is None:
        # Airy-equivalent broadside main beam gives a scale initializer only.
        largest=np.sqrt(np.max(np.var(objective.xy,axis=1)))
        initial=[np.clip(4.65/largest,.32,.88),orientation,.02,.02,3]
    result=least_squares(lambda v:objective.evaluate(v)[0],initial,
        jac=lambda v:objective.evaluate(v)[1],
        bounds=([.3,orientation-15,0,0,0],[.9,orientation+15,.8,.8,8]),
        max_nfev=max_nfev,ftol=1e-7,xtol=1e-7,gtol=1e-7)
    return dict(rows=rows,cols=cols,sites=objective.xy.shape[1],ports=objective.xy.shape[1],
        footprint="ellipse",lattice="triangular",parameters=result.x.tolist(),
        training_rmse_db=float(np.sqrt(np.mean(result.fun**2))),converged=bool(result.success),
        evaluations=int(result.nfev))


def main():
    records=load_data(); samples=training_samples(records)
    # Balanced deterministic subset for the discrete sweep; refine leading
    # hypotheses with every training sample before making the final selection.
    ids=np.linspace(0,len(samples[-1])-1,420).astype(int)
    sparse=(samples[0][:,ids],samples[1][ids],samples[2][ids],samples[3][ids])
    candidates=[]
    for rows in range(18,57):
        for cols in range(12,39):
            cfg=dict(rows=rows,cols=cols,footprint="ellipse",parameters=[1,0,0,0,0])
            count=geometry(cfg)[0].shape[1]
            if 450<=count<=550:
                candidates.append(fit(rows,cols,sparse))
        if rows%5==0:
            print(f"Ellipse rows through {rows}: {len(candidates)} candidates",flush=True)
    candidates.sort(key=lambda c:c["training_rmse_db"])
    refined=[]
    for c in candidates[:12]:
        refined.append(fit(c["rows"],c["cols"],samples,100,c["parameters"]))
    refined.sort(key=lambda c:c["training_rmse_db"])
    selected=refined[0]
    selected |= dict(frequency_hz=FC,parameter_names=["spacing_lambda","rotation_deg","row_taper","col_taper","element_cos_power"],
        coordinate_convention="theta from +z; phi from +x toward +y",
        counting="~500 elements per RX or TX direction; polarization does not double count",
        polarization="selectable ideal LHCP or RHCP; not measured",site_count_prior=[450,550],
        element_peak_gain_dbi_assumed=0,selection="training-only ellipse sweep then full-training refinement")
    (OUT/"ellipse_config.json").write_text(json.dumps(selected,indent=2)+"\n")
    fields=("rows","cols","sites","ports","training_rmse_db","converged","evaluations")
    write_csv(OUT/"ellipse_candidate_ranking.csv",[
        {k:c[k] for k in fields}|dict(zip(selected["parameter_names"],c["parameters"])) for c in candidates])
    (OUT/"ellipse_refined_candidates.json").write_text(json.dumps(refined,indent=2)+"\n")
    results=[];predictions=[]
    for index,r in enumerate(records):
        prediction=predict(selected,r["pitch"],r["azimuth"])
        predictions.append(prediction.astype(np.float32))
        results.append({k:r[k] for k in ("filename","pitch","azimuth")}|{"split":split(r)}|
            metrics(r["gain"],prediction,r["pitch"],r["azimuth"]))
        if index%13==0:print(f"Ellipse validation file {index+1}/78",flush=True)
    write_csv(OUT/"ellipse_per_beam_metrics.csv",results)
    np.savez_compressed(OUT/"ellipse_predicted_patterns.npz",gains=np.stack(predictions),
                       pitch=[r["pitch"] for r in records],azimuth=[r["azimuth"] for r in records])
    summary=dict(candidate_count=len(candidates),refined_count=len(refined))
    for label in ("train","holdout_pitch","holdout_azimuth","duplicate_360"):
        group=[r for r in results if r["split"]==label]
        names=[k for k in group[0] if k.endswith("_db") or k.endswith("_deg")]
        summary[label]={"files":len(group)}|{k:float(np.nanmean([r[k] for r in group])) for k in names}
    held=[r for r in results if r["split"].startswith("holdout")]
    summary["holdout_high_gain_rmse_db"]=float(np.mean([r["high_gain_rmse_db"] for r in held]))
    summary["holdout_above_minus20_rmse_db"]=float(np.mean([r["above_minus20_rmse_db"] for r in held]))
    summary["holdout_full_clip65_rmse_db"]=float(np.mean([r["full_clip65_rmse_db"] for r in held]))
    summary["quality"]="finite triangular ellipse hypothesis; deep sidelobes/nulls remain uncertain"
    (OUT/"ellipse_summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    p,taper=geometry(selected)
    commands=np.array([(r["pitch"],r["azimuth"]) for r in records if r["azimuth"]<360])
    U=direction(commands[:,0],commands[:,1])
    W=taper[:,None]*np.exp(2j*np.pi*FC/C0*(p.T@U))/np.linalg.norm(taper)
    savemat(OUT/"ellipse_physical_model.mat",dict(fc=FC,p_site=p,taper=taper,commands=commands,W_receive=W,
        element_cos_power=selected["parameters"][4],rows=selected["rows"],cols=selected["cols"],
        parameters=selected["parameters"],Nsite=p.shape[1],Nports=p.shape[1],
        cp_jones=np.array([[1,1],[1j,-1j]])/np.sqrt(2)),do_compression=True)
    print(json.dumps(selected,indent=2),flush=True)
    print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__":main()
