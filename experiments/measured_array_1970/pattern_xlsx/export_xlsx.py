"""Export calibrated sheet beams in the original 92-row x 361-column layout."""
import csv
import hashlib
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
from openpyxl import Workbook, load_workbook
from openpyxl.comments import Comment
from openpyxl.styles import Font, PatternFill
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
FLOOR_DB = -65.0


def main():
    # Verify the source-file convention before constructing matching workbooks.
    reference = HERE.parents[2] / "data/measured_patterns/pattern_Pitch_0_Azimuth_0.txt"
    lines = reference.read_text(encoding="utf-8").splitlines()
    assert lines[5].split("\t") == ["theta\\phi"] + [str(x) for x in range(360)]
    raw = np.loadtxt(reference, delimiter="\t", skiprows=6)
    assert raw.shape == (91, 361) and np.array_equal(raw[:, 0], np.arange(91))
    d = loadmat(OUT / "sheet_pattern_grid.mat", simplify_cells=True)
    assert d["power"].shape == (91, 360, 14)
    assert np.array_equal(d["theta_deg"], np.arange(91))
    assert np.array_equal(d["phi_deg"], np.arange(360))
    records = []
    files = []
    header_fill = PatternFill("solid", fgColor="E7EEF7")
    for k, name in enumerate(d["beam_name"]):
        power = d["power"][:, :, k]
        reference_peak = float(power.max())
        normalized = 10 * np.log10(np.maximum(power / reference_peak, 10 ** (FLOOR_DB / 10)))
        values = np.round(np.clip(normalized, FLOOR_DB, 0), 2)
        pitch, azi = float(d["pitch_deg"][k]), float(d["azimuth_deg"][k])
        filename = f"pattern_Pitch_{pitch:g}_Azimuth_{azi:g}.xlsx"
        metadata = {"beam": str(name), "filename": filename, "pitch_deg": pitch, "azimuth_deg": azi,
                    "fc_hz": float(d["fc_hz"]), "receive_elements": int(d["elements"]),
                    "polarization": str(d["polarization"]), "extrapolated": bool(d["extrapolated"][k]),
                    "peak_error_deg": float(d["peak_error_deg"][k]),
                    "normalization_reference_power": reference_peak,
                    "floor_db": FLOOR_DB, "floor_count": int(np.sum(normalized <= FLOOR_DB)),
                    "data_origin": "full-data fitted model; not a new measurement"}
        wb = Workbook()
        ws = wb.active
        ws.title = "Sheet1"
        ws.append(["theta\\phi"] + list(range(360)))
        for theta, row in enumerate(values):
            ws.append([theta] + row.tolist())
        ws.freeze_panes = "B2"
        ws["A1"].comment = Comment(json.dumps(metadata, ensure_ascii=False, indent=2), "QuaDRiGa export")
        ws.column_dimensions["A"].width = 12
        for cell in ws[1]:
            cell.fill = header_fill
            cell.font = Font(bold=True)
        for row in ws.iter_rows(min_row=2):
            row[0].fill = header_fill
            for cell in row[1:]:
                cell.number_format = "0.00"
        path = OUT / filename
        wb.save(path)
        # Re-open every workbook and compare every exported numeric cell.
        check = load_workbook(path, read_only=True, data_only=True)
        assert check.sheetnames == ["Sheet1"]
        sheet = check["Sheet1"]
        assert (sheet.max_row, sheet.max_column) == (92, 361)
        actual = list(sheet.values)
        assert actual[0] == tuple(["theta\\phi"] + list(range(360)))
        assert np.array_equal([row[0] for row in actual[1:]], np.arange(91))
        actual_values = np.asarray([row[1:] for row in actual[1:]], dtype=float)
        np.testing.assert_array_equal(actual_values, values)
        assert actual_values.max() == 0 and actual_values.min() >= FLOOR_DB
        assert np.isfinite(actual_values).all()
        check.close()
        metadata["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        records.append(metadata)
        files.append(path)

    manifest = {"beams": records, "theta_deg": [0, 90, 1], "phi_deg": [0, 359, 1],
                "data_shape": [91, 360], "sheet_shape": [92, 361],
                "sheet_name": "Sheet1", "normalized_peak_db": 0, "decimal_places": 2,
                "normalization": "each complete 91x360 pattern peak; before rounding",
                "validation": "all 458640 gain cells and all axis labels re-opened and compared exactly"}
    manifest_path = OUT / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    index_path = OUT / "beam_index.csv"
    with index_path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=["beam", "filename", "pitch_deg", "azimuth_deg", "extrapolated"],
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    bundle = OUT / "14_sheet_beam_patterns_xlsx.zip"
    with ZipFile(bundle, "w", ZIP_DEFLATED) as archive:
        for path in files + [index_path, manifest_path, HERE / "README.md"]:
            archive.write(path, path.name)
    with ZipFile(bundle) as archive:
        assert archive.testzip() is None
        assert len([n for n in archive.namelist() if n.endswith(".xlsx")]) == 14
    print("PASS: 14 XLSX workbooks, all 458640 gain cells checked, and ZIP bundle verified")


if __name__ == "__main__":
    main()
