"""Power-only identification of a CP-capable triangular phased-array hypothesis.

File Pitch/Azimuth are steering commands; matrix theta/phi are observation
angles in a fixed frame. +z is the array normal; x/y span the aperture.
TX and RX each have about 500 elements; CP choice does NOT double this count.
The closed-form AF is exactly the weighted finite lattice, not a fitted beam
shape exponent. Absolute gain, CP purity and complex phase are not measured.
"""
from pathlib import Path
import hashlib
import re

import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data/measured_patterns"
FC = 1.97e9
C0 = 299792458.0
THETA, PHI = np.meshgrid(np.arange(91.0), np.arange(360.0), indexing="ij")
UV = np.array([np.sin(np.deg2rad(THETA)) * np.cos(np.deg2rad(PHI)),
               np.sin(np.deg2rad(THETA)) * np.sin(np.deg2rad(PHI))])


def load_data():
    records = []
    for path in DATA.glob("pattern_*.txt"):
        pitch, azimuth = map(int, re.fullmatch(r"pattern_Pitch_(\d+)_Azimuth_(\d+)", path.stem).groups())
        axes = np.loadtxt(path, skiprows=5, max_rows=1, usecols=range(1, 361))
        raw = np.loadtxt(path, skiprows=6)
        assert raw.shape == (91, 361) and np.array_equal(raw[:, 0], np.arange(91))
        assert np.array_equal(axes, np.arange(360)) and np.isfinite(raw).all()
        g = raw[:, 1:]
        assert np.ptp(g[0]) < 1e-10 and g.max() <= 0.01 and g.min() >= -65
        peak = np.unravel_index(np.argmax(g), g.shape)
        records.append(dict(pitch=pitch, azimuth=azimuth, filename=path.name,
                            sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                            gain=g, peak_theta=int(peak[0]), peak_phi=int(peak[1])))
    records.sort(key=lambda r: (r["pitch"], r["azimuth"]))
    assert len(records) == 78
    assert {(r["pitch"], r["azimuth"]) for r in records} == {
        (p, a) for p in range(0, 51, 10) for a in range(0, 361, 30)}
    table = {(r["pitch"], r["azimuth"]): r for r in records}
    for pitch in range(0, 51, 10):
        assert np.array_equal(table[pitch, 0]["gain"], table[pitch, 360]["gain"])
    return records


def split(record):
    p, a = record["pitch"], record["azimuth"]
    if a == 360:
        return "duplicate_360"
    if p == 0 or (p in (10, 30, 50) and a % 60 == 0):
        return "train"
    return "holdout_pitch" if p in (20, 40) else "holdout_azimuth"


def direction(theta, phi):
    t, p = np.deg2rad(theta), np.deg2rad(phi)
    return np.array([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)])


def dirichlet(n, t):
    # Wrap first so that sinc(t) never has a zero at a nonzero integer.
    t = (t + 0.5) % 1.0 - 0.5
    return np.exp(1j * np.pi * (n - 1) * t) * n * np.sinc(n * t) / np.sinc(t)


def tapered_sum(n, t, alpha):
    # Positive raised-cosine voltage taper: w[m]=1-alpha*cos(2*pi*m/(n-1)).
    return dirichlet(n, t) - alpha / 2 * (
        dirichlet(n, t + 1 / (n - 1)) + dirichlet(n, t - 1 / (n - 1)))


def array_factor(params, rows, cols, delta_uv, lattice="triangular"):
    d, rotation, ar, ac, _ = params
    c, s = np.cos(np.deg2rad(rotation)), np.sin(np.deg2rad(rotation))
    row_space = np.sqrt(3) / 2 if lattice == "triangular" else 1.0
    tx = d * row_space * (c * delta_uv[0] + s * delta_uv[1])
    ty = d * (-s * delta_uv[0] + c * delta_uv[1])
    row_af = tapered_sum(rows, tx, ar)
    if lattice == "triangular":
        phase = np.exp(1j * np.pi * ty)
        row_af = (1 + phase) / 2 * row_af + (1 - phase) / 2 * tapered_sum(rows, tx + 0.5, ar)
    return row_af * tapered_sum(cols, ty, ac) / ((rows - ar) * (cols - ac))


def relative_db(params, rows, cols, delta_uv, theta, pitch, lattice="triangular"):
    af = array_factor(params, rows, cols, delta_uv, lattice)
    # Coherent AF and a cosine^q POWER element pattern. No empirical floor.
    pattern = 20 * np.log10(np.maximum(abs(af), 1e-15))
    element = 10 * params[4] * np.log10(
        np.maximum(np.cos(np.deg2rad(theta)), 1e-12) / np.cos(np.deg2rad(pitch)))
    return pattern + element


def predict(config, pitch, azimuth):
    delta = UV - direction(pitch, azimuth)[:2, None, None]
    if config.get("footprint") == "ellipse":
        p, taper = geometry(config)
        af = np.empty(THETA.size, dtype=complex)
        flat = delta.reshape(2, -1)
        for start in range(0, flat.shape[1], 2048):
            phase = 2j*np.pi*FC/C0 * (p[:2].T@flat[:, start:start+2048])
            af[start:start+2048] = taper@np.exp(phase)/taper.sum()
        af = af.reshape(THETA.shape)
        g = 20*np.log10(np.maximum(abs(af), 1e-15)) + 10*config["parameters"][4]*np.log10(
            np.maximum(np.cos(np.deg2rad(THETA)),1e-12)/np.cos(np.deg2rad(pitch)))
    else:
        g = relative_db(config["parameters"], config["rows"], config["cols"],
                        delta, THETA, pitch, config.get("lattice", "triangular"))
    # Normalize the WHOLE sampled pattern, as in the supplied normalized data.
    # A grating lobe is not hidden by normalizing just the commanded direction.
    return np.maximum(g - g.max(), -120.0)


def geometry(config):
    rows, cols = config["rows"], config["cols"]
    d, rotation, ar, ac, _ = config["parameters"]
    m, n = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    x = d * (np.sqrt(3) / 2 if config.get("lattice", "triangular") == "triangular" else 1) * m
    y = d * (n + (0.5 * (m % 2) if config.get("lattice", "triangular") == "triangular" else 0))
    x, y = x - x.mean(), y - y.mean()
    keep = np.ones(x.shape, dtype=bool)
    if config.get("footprint") == "ellipse":
        rx = d*np.sqrt(3)/2*(rows-1)/2
        ry = d*(cols-1)/2
        keep = (x/rx)**2+(y/ry)**2 <= 1+1e-12
        x, y = x[keep], y[keep]
        x, y = x-x.mean(), y-y.mean()
    c, s = np.cos(np.deg2rad(rotation)), np.sin(np.deg2rad(rotation))
    p = np.array([(c*x-s*y).ravel(), (s*x+c*y).ravel(), np.zeros(x.size)]) * C0 / FC
    taper = ((1-ar*np.cos(2*np.pi*m/(rows-1))) * (1-ac*np.cos(2*np.pi*n/(cols-1))))[keep].ravel()
    return p, taper


def power_interpolator(gain):
    # Interpolate POWER, including the periodic azimuth seam. Never treat dB as field.
    power = 10.0 ** (gain / 10)
    extended = np.column_stack([power, power[:, 0]])
    return RegularGridInterpolator((np.arange(91), np.arange(361)), extended,
                                   bounds_error=False, fill_value=np.nan)


def great_circle_cut(gain, pitch, azimuth, tangent="radial"):
    t, a = np.deg2rad(pitch), np.deg2rad(azimuth)
    u0 = direction(pitch, azimuth)
    e = (np.array([np.cos(t)*np.cos(a), np.cos(t)*np.sin(a), -np.sin(t)])
         if tangent == "radial" else np.array([-np.sin(a), np.cos(a), 0.0]))
    offset = np.linspace(-15, 15, 601)
    r = np.deg2rad(offset)
    u = u0[:, None]*np.cos(r) + e[:, None]*np.sin(r)
    theta = np.rad2deg(np.arccos(np.clip(u[2], -1, 1)))
    phi = np.rad2deg(np.arctan2(u[1], u[0])) % 360
    power = power_interpolator(gain)(np.column_stack([theta, phi]))
    return offset, 10*np.log10(np.maximum(power, 1e-12))


def width_3db(offset, gain):
    k = int(np.nanargmax(gain))
    level = gain[k] - 3
    left = k
    right = k
    while left > 0 and gain[left-1] >= level:
        left -= 1
    while right < len(gain)-1 and gain[right+1] >= level:
        right += 1
    if left == 0 or right == len(gain)-1:
        return float("nan")
    lo = np.interp(level, [gain[left-1], gain[left]], [offset[left-1], offset[left]])
    hi = np.interp(level, [gain[right+1], gain[right]], [offset[right+1], offset[right]])
    return float(hi-lo)


def metrics(measured, predicted, pitch, azimuth):
    observed = measured - measured.max()
    valid = THETA > 0  # theta=0 repeats the same physical point 360 times
    solid_angle = np.sin(np.deg2rad(THETA))
    result = {}
    for label, mask, floor in (
        ("high_gain", valid & (observed >= -10), -65),
        ("above_minus20", valid & (observed >= -20), -65),
        ("full_clip40", valid, -40),
        ("full_clip65", valid, -65),
    ):
        delta = np.maximum(predicted, floor) - np.maximum(observed, floor)
        result[f"{label}_rmse_db"] = float(np.sqrt(np.mean(delta[mask]**2)))
        result[f"{label}_solid_angle_rmse_db"] = float(np.sqrt(
            np.average(delta[mask]**2, weights=solid_angle[mask])))
    im = np.unravel_index(np.argmax(observed), observed.shape)
    ip = np.unravel_index(np.argmax(predicted), predicted.shape)
    result["peak_error_deg"] = float(np.rad2deg(np.arccos(np.clip(
        np.dot(direction(*im), direction(*ip)), -1, 1))))
    for tangent in ("radial", "tangential"):
        x, m = great_circle_cut(observed, pitch, azimuth, tangent)
        _, p = great_circle_cut(predicted, pitch, azimuth, tangent)
        result[f"measured_{tangent}_width_deg"] = width_3db(x, m)
        result[f"model_{tangent}_width_deg"] = width_3db(x, p)
    return result
