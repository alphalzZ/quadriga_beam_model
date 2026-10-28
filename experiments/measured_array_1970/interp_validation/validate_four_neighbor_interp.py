"""Leave-one-out validation of the four-neighbour transported-power interpolation.

The method under test (see time_domain/make_allfit_link_arrays.m):
  each neighbour's full power pattern is rotated by the minimal 3D rotation
  from its scan command to the target direction, then the rotated patterns
  are averaged with weights in the linear power domain.

Validation design:
  * truth is a measured state (the 78 files are on a 10 deg x 30 deg grid);
  * the held-out state beam_a is removed and reconstructed from its four
    surrounding measured neighbours beam_w/x/y/z at (p-10,a-30),(p-10,a+30),
    (p+10,a-30),(p+10,a+30), i.e. the target sits at the cell centre with
    equal weights 0.25;
  * pitch 0 and 50 have no neighbour on one side, so they cannot be held out
    with a centred cell and are reported as excluded (26 files); the 48 files
    with pitch 10..40 and azimuth 0..330 are the held-out test set. Azimuth 360
    is an exact duplicate of 0 and is not used as a separate case.

Usage (from any directory):
  python validate_four_neighbor_interp.py            # parity check + validation
"""

import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.io import loadmat
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DATA = REPO / "data" / "measured_patterns"
TABLES = REPO / "experiments/measured_array_1970/all_data_fit/results/calibration_tables.mat"
OCTAVE_REF = HERE / "results" / "octave_reference_C1.mat"
RESULTS = HERE / "results"

PITCHES = [0, 10, 20, 30, 40, 50]
AZIMUTHS = list(range(0, 360, 30))
THETA = np.arange(91, dtype=float)  # 0..90 deg, row index
PHI = np.arange(360, dtype=float)  # 0..359 deg, column index
FLOOR_DB = -80.0


# ----------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------
def load_pattern(pitch, azimuth):
    """Measured gain in dB, shape (91, 360) indexed [theta, phi]."""
    az = azimuth % 360
    f = DATA / f"pattern_Pitch_{pitch}_Azimuth_{az}.txt"
    return np.loadtxt(f, skiprows=6, delimiter="\t")[:, 1:]


def check_tables_match_files():
    """The MATLAB code reads baseline_db + residual_db; confirm it equals the txt data."""
    tables = loadmat(TABLES, variable_names=["baseline_db", "residual_db", "commands"])
    cmd = tables["commands"]
    worst = 0.0
    for k in range(cmd.shape[0]):
        p, a = int(cmd[k, 0]), int(cmd[k, 1])
        table = tables["baseline_db"][:, :, k] + tables["residual_db"][:, :, k]
        worst = max(worst, float(np.abs(table - load_pattern(p, a)).max()))
    return worst


# ----------------------------------------------------------------------------
# Method: direct port of make_allfit_link_arrays.m
# ----------------------------------------------------------------------------
def dirvec(theta_deg, phi_deg):
    t = np.deg2rad(theta_deg)
    p = np.deg2rad(phi_deg)
    return np.stack([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)], axis=0)


def align_rotation(source, target):
    """Minimal rotation taking unit vector source to target (Rodrigues form)."""
    v = np.cross(source, target)
    c = float(np.clip(np.dot(source, target), -1.0, 1.0))
    assert c > -1 + 1e-10, "Antipodal transport is undefined"
    V = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + V + V @ V / (1 + c)


def transported_power(corners, target_theta, target_phi, query_theta, query_phi):
    """Linear power at query directions from weighted, transported corner patterns.

    corners: list of (pitch, azimuth, weight, gain_db(91,360)).
    query_theta/query_phi: arrays of identical shape, degrees, canonical frame.
    """
    u = dirvec(np.ravel(query_theta), np.ravel(query_phi))
    target = dirvec(target_theta, target_phi)
    power = np.zeros(u.shape[1])
    for pitch, azimuth, weight, gain_db in corners:
        if weight == 0:
            continue
        Q = align_rotation(dirvec(pitch, azimuth), target)
        query = Q.T @ u
        tq = np.degrees(np.arccos(np.clip(query[2], -1, 1)))
        aq = np.mod(np.degrees(np.arctan2(query[1], query[0])), 360)
        tq[(tq > 90) & (tq < 90 + 1e-8)] = 90  # numerical boundary only
        pg = 10 ** (gain_db / 10)
        pg = np.concatenate([pg, pg[:, :1]], axis=1)  # close the azimuth loop
        interp = RegularGridInterpolator(
            (THETA, np.arange(361, dtype=float)), pg, bounds_error=False, fill_value=0.0
        )
        power += weight * interp(np.column_stack([tq, aq]))
    return power.reshape(np.shape(query_theta))


def four_neighbour_corners(pitch, azimuth):
    """Centred cell around (pitch, azimuth): four equal-weight neighbours."""
    out = []
    for dp in (-10, 10):
        for da in (-30, 30):
            out.append((pitch + dp, (azimuth + da) % 360, 0.25, load_pattern(pitch + dp, azimuth + da)))
    return out


def predict_grid(corners, pitch, azimuth):
    """Predicted gain (dB, floored) on the measured 91 x 360 grid."""
    TH, PH = np.meshgrid(THETA, PHI, indexing="ij")
    p = transported_power(corners, pitch, azimuth, TH, PH)
    return np.maximum(10 * np.log10(np.maximum(p, 1e-300)), FLOOR_DB)


# ----------------------------------------------------------------------------
# Baselines (no transport)
# ----------------------------------------------------------------------------
def baseline_db_average(corners):
    return sum(w * g for _, _, w, g in corners)


def baseline_linear_average(corners):
    p = sum(w * 10 ** (g / 10) for _, _, w, g in corners)
    return np.maximum(10 * np.log10(np.maximum(p, 1e-300)), FLOOR_DB)


# ----------------------------------------------------------------------------
# Parity check against the MATLAB function run in Octave (target C1)
# ----------------------------------------------------------------------------
def parity_check_against_octave():
    if not OCTAVE_REF.exists():
        return None
    from scipy.io import loadmat as lm

    ref = lm(OCTAVE_REF)
    az = ref["az_grid_deg"].ravel()
    el = ref["el_grid_deg"].ravel()
    AZ, EL = np.meshgrid(az, el)  # rows = elevation, like MATLAB meshgrid
    # MATLAB u = [cos(EL)cos(AZ); cos(EL)sin(AZ); sin(EL)], so z = sin(EL), theta = acos(z).
    u_z = np.sin(np.deg2rad(EL))
    theta_q = np.degrees(np.arccos(np.clip(u_z, -1, 1)))
    phi_q = np.mod(np.degrees(np.arctan2(np.sin(np.deg2rad(AZ)) * np.cos(np.deg2rad(EL)),
                                          np.cos(np.deg2rad(AZ)) * np.cos(np.deg2rad(EL)))), 360)
    commands = ref["commands"]
    weights = ref["weights"].ravel()
    corners = [
        (int(c[0]), int(c[1]), float(w), load_pattern(int(c[0]), int(c[1])))
        for c, w in zip(commands, weights)
    ]
    # Target C1 is the off-grid direction (45.43, 188.62).
    p_py = transported_power(corners, 45.43, 188.62, theta_q, phi_q)
    p_py = p_py.reshape(AZ.shape)
    # Points with z<0 are the unmeasured rear hemisphere: MATLAB sets them to 0.
    p_py[u_z < -1e-12] = 0
    p_ml = ref["power_grid"]
    ok = np.isfinite(p_ml)
    rel = np.max(np.abs(p_py[ok] - p_ml[ok]))
    scale = np.max(np.abs(p_ml[ok]))
    return {"max_abs_diff_linear": float(rel), "peak_linear_octave": float(scale)}


# ----------------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------------
SOLID = np.sin(np.deg2rad(THETA))[:, None] * np.ones((1, 360)) * np.deg2rad(1) ** 2


def weighted_rms(err, mask):
    w = SOLID * mask
    return float(np.sqrt(np.sum(w * err ** 2) / np.sum(w)))


def unit_of(idx):
    th, ph = np.unravel_index(idx, (91, 360))
    return dirvec(THETA[th], PHI[ph])


def metrics(pred_db, meas_db):
    err = pred_db - meas_db
    m3 = meas_db >= -3
    m10 = meas_db >= -10
    ones = np.ones_like(meas_db, dtype=bool)
    p3 = pred_db >= -3
    omega = lambda m: float(np.sum(SOLID * m))  # noqa: E731
    peak_p = int(np.argmax(pred_db))
    peak_m = int(np.argmax(meas_db))
    cosang = np.clip(np.dot(unit_of(peak_p), unit_of(peak_m)), -1, 1)
    return {
        "rms_db_hemisphere": weighted_rms(err, ones),
        "rms_db_within_m10": weighted_rms(err, m10),
        "rms_db_within_m3": weighted_rms(err, m3),
        "rms_db_m20_to_m10": weighted_rms(err, (meas_db < -10) & (meas_db >= -20)),
        "rms_db_m40_to_m20": weighted_rms(err, (meas_db < -20) & (meas_db >= -40)),
        "rms_db_below_m40": weighted_rms(err, meas_db < -40),
        "bias_db_within_m10": float(np.sum(SOLID * m10 * err) / np.sum(SOLID * m10)),
        "max_abs_db_within_m10": float(np.max(np.abs(err[m10]))),
        "peak_level_err_db": float(pred_db.max() - meas_db.max()),
        "peak_dir_err_deg": float(np.degrees(np.arccos(cosang))),
        "beam3db_solid_angle_ratio": omega(p3) / omega(m3),
        "beam3db_iou": omega(p3 & m3) / omega(p3 | m3),
    }


METHODS = ["transport_linear", "db_average_no_transport", "linear_average_no_transport"]


def evaluate_case(pitch, azimuth):
    corners = four_neighbour_corners(pitch, azimuth)
    meas = load_pattern(pitch, azimuth)
    preds = {
        "transport_linear": predict_grid(corners, pitch, azimuth),
        "db_average_no_transport": baseline_db_average(corners),
        "linear_average_no_transport": baseline_linear_average(corners),
    }
    return {name: metrics(p, meas) for name, p in preds.items()}, preds, meas, corners


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    summary = {"table_vs_txt_max_abs_db": check_tables_match_files()}
    summary["octave_parity"] = parity_check_against_octave()

    cases = [(p, a) for p in [10, 20, 30, 40] for a in AZIMUTHS]
    rows = []
    for p, a in cases:
        res, _, _, _ = evaluate_case(p, a)
        for method, m in res.items():
            rows.append({"pitch": p, "azimuth": a, "method": method, **m})

    with open(RESULTS / "per_case_metrics.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    keys = [k for k in rows[0] if k not in ("pitch", "azimuth", "method")]
    agg = {}
    for method in METHODS:
        sub = [r for r in rows if r["method"] == method]
        agg[method] = {}
        for k in keys:
            vals = np.array([r[k] for r in sub])
            agg[method][k] = {"median": float(np.median(vals)), "worst": float(vals[np.argmax(np.abs(vals))])}
    summary["n_held_out_cases"] = len(cases)
    summary["aggregate"] = agg

    # Example beam_a = Pitch 20 / Azimuth 90 with neighbours w/x/y/z
    res, preds, meas, corners = evaluate_case(20, 90)
    summary["example_beam_a_20_90"] = res
    summary["example_corners"] = [[c[0], c[1]] for c in corners]
    with open(RESULTS / "summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    make_figures(rows, preds, meas, corners)
    print(json.dumps(summary, indent=2, ensure_ascii=False)[:4000])


def make_figures(rows, preds, meas, corners):
    # Figure 1: example beam_a
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
    pred = preds["transport_linear"]
    err = pred - meas
    im0 = ax[0, 0].pcolormesh(PHI, THETA, meas, vmin=-40, vmax=0, shading="auto", cmap="viridis")
    ax[0, 0].set_title("beam_a 实测 (dB)")
    ax[0, 1].pcolormesh(PHI, THETA, pred, vmin=-40, vmax=0, shading="auto", cmap="viridis")
    ax[0, 1].set_title("四邻居插值 (dB)")
    im2 = ax[0, 2].pcolormesh(PHI, THETA, err, vmin=-10, vmax=10, shading="auto", cmap="RdBu_r")
    ax[0, 2].set_title("插值 − 实测 (dB)")
    for a in ax[0]:
        a.set_xlabel("φ (deg)")
        a.set_ylabel("θ (deg)")
        a.plot(90, 20, "w+", ms=12)
    fig.colorbar(im0, ax=ax[0, :2], shrink=0.8)
    fig.colorbar(im2, ax=ax[0, 2], shrink=0.8)

    # theta cut at phi=90
    j = 90
    ax[1, 0].plot(THETA, meas[:, j], "k-", lw=2, label="实测 beam_a")
    ax[1, 0].plot(THETA, pred[:, j], "r--", label="四邻居插值（转换+线性功率）")
    ax[1, 0].plot(THETA, preds["db_average_no_transport"][:, j], ":", color="tab:blue", label="dB 直接平均（无转换）")
    ax[1, 0].set_xlabel("θ (deg), φ=90")
    ax[1, 0].set_ylabel("dB")
    ax[1, 0].set_ylim(-60, 2)
    ax[1, 0].legend(fontsize=8)
    ax[1, 0].grid(alpha=0.3)

    # phi cut at theta=20
    i = 20
    ax[1, 1].plot(PHI, meas[i, :], "k-", lw=2, label="实测 beam_a")
    ax[1, 1].plot(PHI, pred[i, :], "r--", label="四邻居插值")
    ax[1, 1].plot(PHI, preds["db_average_no_transport"][i, :], ":", color="tab:blue", label="dB 直接平均（无转换）")
    ax[1, 1].plot(PHI, preds["linear_average_no_transport"][i, :], "-.", color="tab:green", label="线性平均（无转换）")
    ax[1, 1].set_xlabel("φ (deg), θ=20")
    ax[1, 1].set_ylabel("dB")
    ax[1, 1].set_ylim(-60, 2)
    ax[1, 1].legend(fontsize=8)
    ax[1, 1].grid(alpha=0.3)

    ax[1, 2].hist(err[meas >= -10], bins=60, color="tab:orange", alpha=0.8)
    ax[1, 2].set_xlabel("误差 (dB), 实测 ≥ −10 dB 区域")
    ax[1, 2].set_ylabel("计数")
    ax[1, 2].grid(alpha=0.3)

    fig.suptitle(
        "beam_a = Pitch 20°/Azimuth 90°，邻居 (10,60) (10,120) (30,60) (30,120)，权重各 0.25", fontsize=12
    )
    fig.tight_layout()
    fig.savefig(RESULTS / "example_beam_a_20_90.png", dpi=130)
    plt.close(fig)

    # Figure 2: leave-one-out summary
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    pitches = [10, 20, 30, 40]
    grid = np.zeros((len(pitches), len(AZIMUTHS)))
    for r in rows:
        if r["method"] != "transport_linear":
            continue
        grid[pitches.index(r["pitch"]), AZIMUTHS.index(r["azimuth"])] = r["rms_db_within_m10"]
    im = ax[0].imshow(grid, cmap="magma_r", aspect="auto")
    ax[0].set_xticks(range(len(AZIMUTHS)))
    ax[0].set_xticklabels(AZIMUTHS, rotation=90, fontsize=7)
    ax[0].set_yticks(range(len(pitches)))
    ax[0].set_yticklabels(pitches)
    ax[0].set_xlabel("Azimuth (deg)")
    ax[0].set_ylabel("Pitch (deg)")
    ax[0].set_title("留一验证 RMS 误差 (实测 ≥ −10 dB 区域, dB)")
    fig.colorbar(im, ax=ax[0])

    names = ["transport_linear", "db_average_no_transport", "linear_average_no_transport"]
    labels = ["转换+线性功率（本方法）", "dB 直接平均", "线性平均（无转换）"]
    data = [[r["rms_db_within_m10"] for r in rows if r["method"] == n] for n in names]
    ax[1].boxplot(data, tick_labels=labels, showfliers=True)
    ax[1].set_ylabel("RMS 误差 (dB), 实测 ≥ −10 dB")
    ax[1].set_title("48 个留出状态的误差分布")
    ax[1].tick_params(axis="x", labelsize=8)
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(RESULTS / "loo_summary.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
