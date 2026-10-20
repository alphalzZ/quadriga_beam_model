"""Plot fitted measured-grid powers for all 78 commands (no holdout)."""
from pathlib import Path
import sys
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from model import load_data, great_circle_cut


def main():
    out=HERE/'results'; records=load_data()
    data=np.load(out/'fitted_patterns.npz'); fitted=data['gains']
    cache=loadmat(out/'calibration_tables.mat',variable_names=['baseline_db','residual_db'])
    baseline=np.moveaxis(cache['baseline_db'],2,0)
    coefficients=np.moveaxis(cache['residual_db'],2,0)
    assert fitted.shape==(78,91,360)
    for i,r in enumerate(records):
        assert np.array_equal(data['commands'][i],[r['pitch'],r['azimuth']])
        assert np.max(abs(fitted[i]-r['gain']))<1e-10
        assert np.max(abs(baseline[i]+coefficients[i]-r['gain']))<1e-10
    fig,axes=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
    with (out/'all_beam_metrics.csv').open() as f: metrics=list(csv.DictReader(f))
    axes[0].plot([float(r['baseline_rmse_db']) for r in metrics],label='Physical baseline')
    axes[0].plot([float(r['calibrated_rmse_db']) for r in metrics],label='Calibrated all-data fit')
    axes[0].set(xlabel='Command index (all 78 files)',ylabel='Full-grid RMSE [dB]',title='No command or angular sample held out')
    axes[0].legend()
    maximum=np.array([float(r['calibrated_max_error_db']) for r in metrics])
    axes[1].plot(maximum,'o',ms=3)
    axes[1].set(xlabel='Command index',ylabel='Maximum absolute error [dB]',title='Calibrated measurement nodes: floating-point residual')
    fig.savefig(out/'all_commands_error.png',dpi=150);plt.close(fig)
    lookup={(r['pitch'],r['azimuth']):i for i,r in enumerate(records)}
    examples=[(0,0),(10,60),(20,30),(30,150),(40,90),(50,330)]
    fig,axes=plt.subplots(6,2,figsize=(14,21),constrained_layout=True)
    for row,command in enumerate(examples):
        i=lookup[command]
        for col,tangent in enumerate(('radial','tangential')):
            ax=axes[row,col]
            for gain,label,color,style in ((baseline[i],'Physical baseline','tab:orange','-'),
                    (records[i]['gain'],'Measured','black','-'),(fitted[i],'All-data calibrated','tab:blue','--')):
                x,y=great_circle_cut(gain,*command,tangent)
                ax.plot(x,y,color=color,ls=style,lw=1.2,label=label)
            ax.set(xlabel='Great-circle offset [deg]',ylabel='Gain [dB]',title=f'P{command[0]} A{command[1]} | {tangent}',ylim=(-65,1),xlim=(-15,15))
            if row==0:ax.legend(fontsize=8)
            ax.grid(alpha=.2)
    fig.savefig(out/'calibrated_cut_comparisons.png',dpi=140);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(15,9),constrained_layout=True)
    for row,command in enumerate(((20,30),(50,330))):
        i=lookup[command]
        for col,(g,title) in enumerate(((records[i]['gain'],'Measured'),(fitted[i],'Calibrated fit'),(fitted[i]-records[i]['gain'],'Numerical difference'))):
            options=dict(vmin=-65,vmax=0,cmap='viridis') if col<2 else dict(vmin=-1e-12,vmax=1e-12,cmap='RdBu_r')
            im=axes[row,col].imshow(g,origin='lower',aspect='auto',extent=(0,359,0,90),**options)
            axes[row,col].set(xlabel='phi [deg]',ylabel='theta [deg]',title=f'P{command[0]} A{command[1]} | {title}')
            fig.colorbar(im,ax=axes[row,col],label='dB')
    fig.savefig(out/'calibrated_map_comparisons.png',dpi=140);plt.close(fig)
    result=dict(verified_files=78,verified_samples=78*91*360,all_grid_powers_match=True,
                residual_coefficients_round_trip=True)
    (out/'python_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print('All 78 complete angular grids verified; saved three comparison figures.')


if __name__=='__main__':main()
