"""Fit every measured value, with no holdout or regularization.

A 502-element physical baseline is refined using every command's high-gain
samples. A free 2-D nodal residual per command then fits ALL 2,555,280 input
values. The residual least-squares problem has the identity design matrix at
measurement nodes, so coefficients are solved analytically. This is empirical
beam calibration, not unique identification of 502 element-level patterns.
"""
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.io import savemat

PARENT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PARENT/'common'))
sys.path.insert(0,str(PARENT/'baseline_fit'))
from fit_ellipse import fit
from model import C0, FC, THETA, UV, direction, geometry, load_data

OUT=Path(__file__).resolve().parent/"results"


def fields(config,pitch,azimuth):
    p,taper=geometry(config)
    delta=(UV-direction(pitch,azimuth)[:2,None,None]).reshape(2,-1)
    af=np.empty(THETA.size,dtype=complex)
    for start in range(0,delta.shape[1],2048):
        e=np.exp(2j*np.pi*FC/C0*(p[:2].T@delta[:,start:start+2048]))
        af[start:start+2048]=taper@e/taper.sum()
    af=af.reshape(THETA.shape)
    db=20*np.log10(np.maximum(abs(af),1e-15))+10*config['parameters'][4]*np.log10(
        np.maximum(np.cos(np.deg2rad(THETA)),1e-12)/np.cos(np.deg2rad(pitch)))
    db=db-db.max()
    # A finite reference permits nonzero measured gains at physical nulls.
    db=np.maximum(db,-80.0)
    phase=np.angle(af)
    phase[abs(af)<1e-12]=0  # phase at a true model null is unobservable
    return db,phase


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    records=load_data()
    initial=json.loads((PARENT/"baseline_fit/results/ellipse_config.json").read_text())
    xs=[];ts=[];ps=[];ys=[]
    for r in records:
        # No command is excluded, including 360 aliases and broadside repeats.
        mask=(r['gain']>=-15)&(THETA>0)
        xs.append((UV-direction(r['pitch'],r['azimuth'])[:2,None,None])[:,mask])
        ts.append(THETA[mask]);ps.append(np.full(np.sum(mask),r['pitch']))
        ys.append(r['gain'][mask]-r['gain'].max())
    samples=(np.concatenate(xs,axis=1),np.concatenate(ts),np.concatenate(ps),np.concatenate(ys))
    print(f"Refining baseline: 78/78 commands, {len(samples[-1])} high-gain samples",flush=True)
    config=fit(initial['rows'],initial['cols'],samples,100,initial['parameters'])
    config.update(frequency_hz=FC,selection="all 78 commands; no holdout",counting="one RX/TX aperture has 502 elements, not 1004",
                  baseline_high_gain_samples=len(samples[-1]),calibration_floor_db=-80,
                  final_fit="independent 2D nodal dB residuals; all samples; zero regularization")
    (OUT/"all_data_config.json").write_text(json.dumps(config,indent=2)+"\n")
    p,taper=geometry(config);commands=np.array([(r['pitch'],r['azimuth']) for r in records])
    U=direction(commands[:,0],commands[:,1])
    W=taper[:,None]*np.exp(2j*np.pi*FC/C0*(p.T@U))/np.linalg.norm(taper)
    savemat(OUT/"physical_baseline.mat",dict(fc=FC,p_site=p,taper=taper,commands=commands,
        W_receive=W,element_cos_power=config['parameters'][4],Nsite=len(taper),Nports=len(taper),
        rows=config['rows'],cols=config['cols'],parameters=config['parameters']),do_compression=True)

    baselines=[];coefficients=[];phases=[];fitted=[];metrics=[]
    for index,r in enumerate(records):
        base,phase=fields(config,r['pitch'],r['azimuth'])
        residual=r['gain']-base
        fit_db=base+residual
        error=fit_db-r['gain']
        assert np.max(abs(error))<1e-10
        baselines.append(base);coefficients.append(residual);phases.append(phase.astype(np.float32));fitted.append(fit_db)
        row={k:r[k] for k in ('filename','pitch','azimuth','sha256')}
        row.update(role="fit",samples=r['gain'].size,
            baseline_rmse_db=float(np.sqrt(np.mean((base-r['gain'])**2))),
            calibrated_rmse_db=float(np.sqrt(np.mean(error**2))),
            calibrated_max_error_db=float(np.max(abs(error))),
            max_abs_calibration_db=float(np.max(abs(residual))))
        metrics.append(row)
        if index%13==12:print(f"Fitted full angular grid: {index+1}/78 commands",flush=True)
    savemat(OUT/"calibration_tables.mat",dict(fc=FC,commands=commands,theta_deg=np.arange(91),phi_deg=np.arange(360),
        baseline_db=np.stack(baselines,axis=2),residual_db=np.stack(coefficients,axis=2),
        assumed_phase_rad=np.stack(phases,axis=2),reference_floor_db=-80),do_compression=True)
    np.savez_compressed(OUT/"fitted_patterns.npz",gains=np.stack(fitted),commands=commands)
    with (OUT/"all_beam_metrics.csv").open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(metrics[0]),lineterminator='\n');writer.writeheader();writer.writerows(metrics)
    summary=dict(input_files=78,fit_files=78,holdout_files=0,angular_samples_per_file=91*360,
        total_fitted_samples=sum(r['samples'] for r in metrics),receive_elements=len(taper),
        calibrated_rmse_db=float(np.sqrt(np.mean([(r['calibrated_rmse_db'])**2 for r in metrics]))),
        calibrated_max_error_db=max(r['calibrated_max_error_db'] for r in metrics),
        phase="from physical baseline; not measured",polarization="selectable ideal LHCP/RHCP; not measured",
        model="all-data empirical beam-state residual calibration on physical baseline",
        geometry_identified_from_calibration=False)
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__":main()
