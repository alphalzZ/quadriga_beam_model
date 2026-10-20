"""Reproduce the audit, training-only lattice search and held-out comparison."""
import argparse
import csv
import json
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.io import savemat
from scipy.optimize import least_squares

from model import (FC, THETA, UV, direction, geometry, load_data, metrics,
                   predict, relative_db, split)


def write_csv(path, rows):
    with Path(path).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def audit(records):
    rows = []
    for r in records:
        error = np.rad2deg(np.arccos(np.clip(np.dot(
            direction(r["pitch"], r["azimuth"]), direction(r["peak_theta"], r["peak_phi"])), -1, 1)))
        rows.append({k: r[k] for k in ("filename", "pitch", "azimuth", "sha256", "peak_theta", "peak_phi")}
                    | {"split": split(r), "fixed_frame_peak_error_deg": float(error),
                       "beam_centered_peak_error_deg": r["peak_theta"],
                       "theta0_db": float(r["gain"][0, 0]), "maximum_db": float(r["gain"].max())})
    return rows


def training_samples(records):
    # Broadside azimuth commands have the same geometric steering direction.
    # Aggregate their repeats rather than giving that direction 12-fold weight.
    broadside = np.median([r["gain"] for r in records if r["pitch"] == 0 and r["azimuth"] < 360], axis=0)
    groups = [(0, 0, broadside)] + [(r["pitch"], r["azimuth"], r["gain"])
             for r in records if split(r) == "train" and r["pitch"] != 0]
    delta, theta, pitch, gains = [], [], [], []
    for p, a, g in groups:
        g = g-g.max()
        mask = (g >= -15) & (THETA > 0)
        indices = np.flatnonzero(mask)
        # Deterministic angular subsampling, bounded per steering direction.
        indices = indices[np.linspace(0, len(indices)-1, min(140, len(indices))).astype(int)]
        delta.append((UV-direction(p, a)[:2, None, None]).reshape(2, -1)[:, indices])
        theta.append(THETA.ravel()[indices])
        pitch.append(np.full(len(indices), p))
        gains.append(g.ravel()[indices])
    return np.concatenate(delta, axis=1), np.concatenate(theta), np.concatenate(pitch), np.concatenate(gains)


def optimize(rows, cols, samples, lattice, max_nfev=35, initial=None):
    x, theta, pitches, observed = samples
    orientation = 0 if rows * (np.sqrt(3)/2 if lattice == "triangular" else 1) >= cols else 90
    lower = [.3, orientation-15, 0, 0, 0]
    upper = [.9, orientation+15, .95, .95, 8]
    start = [.55, orientation, .1, .1, 2] if initial is None else initial
    result = least_squares(lambda v: relative_db(v, rows, cols, x, theta, pitches, lattice)-observed,
                           start, bounds=(lower, upper), max_nfev=max_nfev,
                           ftol=1e-7, xtol=1e-7, gtol=1e-7)
    return dict(rows=rows, cols=cols, sites=rows*cols, ports=rows*cols,
                lattice=lattice, parameters=result.x.tolist(),
                training_rmse_db=float(np.sqrt(np.mean(result.fun**2))),
                converged=bool(result.success), evaluations=int(result.nfev))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent/"results")
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    records = load_data()
    audit_rows = audit(records)
    write_csv(out/"data_audit.csv", audit_rows)
    samples = training_samples(records)
    print(f"Training samples: {len(samples[-1])}, 19 independent steering directions; holdout untouched", flush=True)
    candidates = []
    # ±10% around 500 spatial sites is a user-informed PRIOR, not a measured count.
    for rows in range(12, 51):
        for cols in range(10, 43):
            if not 450 <= rows*cols <= 550:
                continue
            candidates.append(optimize(rows, cols, samples, "triangular"))
        if rows % 5 == 0:
            print(f"Rows through {rows}: {len(candidates)} candidate lattices", flush=True)
    candidates.sort(key=lambda c: c["training_rmse_db"])
    # Resolve optimizer budgets for the leading hypotheses using TRAINING only.
    for index in range(min(12, len(candidates))):
        c = candidates[index]
        candidates[index] = optimize(c["rows"], c["cols"], samples, "triangular", 200, c["parameters"])
    candidates.sort(key=lambda c: c["training_rmse_db"])
    selected = candidates[0]
    selected |= dict(frequency_hz=FC, parameter_names=["spacing_lambda", "rotation_deg", "row_taper", "col_taper", "element_cos_power"],
                     coordinate_convention="theta from +z array normal; phi from +x toward +y",
                     element_peak_gain_dbi_assumed=0.0,
                     polarization="LHCP or RHCP selectable; not measured; one port per element per direction",
                     counting="about 500 RX elements plus about 500 TX elements; NOT 1000 RX elements",
                     site_count_prior=[450, 550], selection="training loss only; no holdout tuning")
    (out/"physical_config.json").write_text(json.dumps(selected, indent=2)+"\n")
    ranking = []
    for c in candidates:
        ranking.append({k: c[k] for k in ("rows", "cols", "sites", "ports", "training_rmse_db", "converged", "evaluations")}
                       | dict(zip(selected["parameter_names"], c["parameters"])))
    write_csv(out/"candidate_ranking.csv", ranking)
    # Rectangular control is diagnostic, never used to change the requested lattice.
    rectangular = []
    for rows, cols in [(32, 16), (16, 32), (34, 15), (15, 34)]:
        rectangular.append(optimize(rows, cols, samples, "rectangular", 100))
    rectangular.sort(key=lambda c: c["training_rmse_db"])
    control = rectangular[0]
    (out/"rectangular_control.json").write_text(json.dumps(rectangular, indent=2)+"\n")
    print("Selected:", selected, flush=True)

    results, controls, predictions = [], [], []
    for index, r in enumerate(records):
        prediction = predict(selected, r["pitch"], r["azimuth"])
        predictions.append(prediction.astype(np.float32))
        meta = {k: r[k] for k in ("filename", "pitch", "azimuth")} | {"split": split(r)}
        results.append(meta | metrics(r["gain"], prediction, r["pitch"], r["azimuth"]))
        controls.append(meta | metrics(r["gain"], predict(control, r["pitch"], r["azimuth"]), r["pitch"], r["azimuth"]))
    write_csv(out/"per_beam_metrics.csv", results)
    write_csv(out/"rectangular_control_metrics.csv", controls)
    # Store predictions only: original measurements remain the authoritative inputs.
    np.savez_compressed(out/"predicted_patterns.npz", gains=np.stack(predictions),
                        pitch=[r["pitch"] for r in records], azimuth=[r["azimuth"] for r in records])
    summary = dict(files=78, exact_duplicate_360_files=6, unique_file_commands=72,
                   geometrically_unique_directions=61, training_independent_directions=19,
                   training_sample_count=len(samples[-1]), candidate_count=len(candidates),
                   parameters_within_5pct_loss=[dict(rows=c["rows"],cols=c["cols"],sites=c["sites"])
                       for c in candidates if c["training_rmse_db"] <= selected["training_rmse_db"]*1.05],
                   runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__))
    for label in ("train", "holdout_pitch", "holdout_azimuth", "duplicate_360"):
        group = [r for r in results if r["split"] == label]
        names = [k for k in group[0] if k.endswith("_db") or k.endswith("_deg")]
        summary[label] = {"files":len(group)} | {k:float(np.nanmean([r[k] for r in group])) for k in names}
    held = [r for r in results if r["split"].startswith("holdout")]
    summary["holdout_mean_high_gain_rmse_db"] = float(np.mean([r["high_gain_rmse_db"] for r in held]))
    summary["holdout_mean_full_clip40_solid_angle_rmse_db"] = float(np.mean([r["full_clip40_solid_angle_rmse_db"] for r in held]))
    summary["physical_main_beam_gate_passed"] = summary["holdout_mean_high_gain_rmse_db"] < 1
    summary["physical_full_pattern_gate_passed"] = bool(
        summary["holdout_mean_full_clip40_solid_angle_rmse_db"] < 3 and
        np.mean([r["above_minus20_rmse_db"] for r in held]) < 3)
    summary["quality"] = "main-beam approximation; full sidelobe/null reproduction not established"
    summary["fixed_frame_peak_error_mean_deg"] = float(np.mean([r["fixed_frame_peak_error_deg"] for r in audit_rows]))
    summary["beam_centered_peak_error_mean_deg"] = float(np.mean([r["beam_centered_peak_error_deg"] for r in audit_rows]))
    (out/"summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")

    # Share one geometry and taper across all steering commands. No per-beam fitting.
    p, taper = geometry(selected)
    commands = np.array([(r["pitch"], r["azimuth"]) for r in records if r["azimuth"] < 360])
    U = direction(commands[:, 0], commands[:, 1])
    receive_weights = taper[:, None] * np.exp(2j*np.pi*FC/299792458 * (p.T@U)) / np.linalg.norm(taper)
    savemat(out/"physical_model.mat", dict(fc=FC, p_site=p, taper=taper, commands=commands,
            W_receive=receive_weights, element_cos_power=selected["parameters"][4],
            rows=selected["rows"], cols=selected["cols"], parameters=selected["parameters"],
            Nsite=p.shape[1], Nports=p.shape[1],
            cp_jones=np.array([[1,1],[1j,-1j]])/np.sqrt(2)), do_compression=True)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
