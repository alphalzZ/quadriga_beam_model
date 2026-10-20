"""Independent mathematical/data checks for the fitted models."""
import json
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist
from scipy.io import loadmat

from fit_ellipse import DiscreteObjective
from fit_model import training_samples
from model import C0, FC, array_factor, direction, geometry, load_data, split


def main():
    out=Path(__file__).resolve().parent/"results"
    records=load_data()
    train=[r for r in records if split(r)=="train"]
    held=[r for r in records if split(r).startswith("holdout")]
    key=lambda r:(r["pitch"],r["azimuth"]%360 if r["pitch"] else 0)
    assert {key(r) for r in train}.isdisjoint({key(r) for r in held})
    assert len(train)==30 and len(held)==42 and len({key(r) for r in train})==19
    assert len({key(r) for r in records})==61
    cfg=json.loads((out/"physical_config.json").read_text())
    p,w=geometry(cfg)
    block_mat=loadmat(out/"physical_model.mat")
    assert int(block_mat["Nports"].item())==p.shape[1]
    assert np.array_equal(block_mat["p_site"],p)
    rng=np.random.default_rng(1970)
    delta=rng.uniform(-.6,.6,(2,140))
    direct=w@np.exp(2j*np.pi*FC/C0*(p[:2].T@delta))/w.sum()
    closed=array_factor(cfg["parameters"],cfg["rows"],cfg["cols"],delta)
    # Centering removes a global phase; power, not arbitrary phase origin, must agree.
    assert np.max(abs(abs(direct)**2-abs(closed)**2))<1e-11
    assert abs(pdist(p.T).min()-cfg["parameters"][0]*C0/FC)<1e-12
    ellipse=json.loads((out/"ellipse_config.json").read_text())
    pe,we=geometry(ellipse)
    assert 450<=pe.shape[1]<=550 and np.max(abs(pe[2]))==0
    assert np.linalg.norm(pe.mean(axis=1))<1e-12
    assert abs(pdist(pe.T).min()-ellipse["parameters"][0]*C0/FC)<1e-12
    mat=loadmat(out/"ellipse_physical_model.mat")
    assert np.array_equal(mat["p_site"],pe)
    assert int(mat["Nports"].item())==pe.shape[1], "RX count must not double for CP"
    assert np.max(abs(np.sum(abs(mat["W_receive"])**2,axis=0)-1))<1e-12
    # Check steering at several nonzero commands against direct summation.
    for index in (0,19,37,71):
        pitch,az=mat["commands"][index]
        u=direction(pitch,az)
        response=np.vdot(mat["W_receive"][:,index],np.exp(2j*np.pi*FC/C0*pe.T@u))
        expected=we.sum()/np.linalg.norm(we)
        assert abs(response-expected)<1e-10
    samples=training_samples(records)
    idx=np.arange(0,100,2)
    small=(samples[0][:,idx],samples[1][idx],samples[2][idx],samples[3][idx])
    objective=DiscreteObjective(ellipse["rows"],ellipse["cols"],small)
    v=np.array(ellipse["parameters"]); _,jac=objective.evaluate(v)
    fd=np.empty(jac.shape)
    for k in range(5):
        dv=np.zeros(5);dv[k]=1e-6
        fd[:,k]=(objective.evaluate(v+dv)[0]-objective.evaluate(v-dv)[0])/(2e-6)
    rel=np.linalg.norm(jac-fd)/max(np.linalg.norm(fd),1)
    assert rel<1e-5,rel
    # All output files must refer to the actual current data set/order.
    for name in ("predicted_patterns.npz","ellipse_predicted_patterns.npz"):
        data=np.load(out/name)
        assert data["gains"].shape==(78,91,360)
        assert np.isfinite(data["gains"]).all()
        assert np.max(abs(data["gains"].max(axis=(1,2))))<1e-5
        for i,r in enumerate(records):
            assert data["pitch"][i]==r["pitch"] and data["azimuth"][i]==r["azimuth"]
    result=dict(data_grid=True,alias_and_split_isolation=True,lattice_geometry=True,
        direct_vs_closed_AF=True,receive_steering=True,rx_element_count=True,
        optimizer_jacobian_relative_error=float(rel),cache_consistency=True)
    (out/"python_verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
