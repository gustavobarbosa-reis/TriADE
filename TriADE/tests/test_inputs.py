import types

import numpy as np
import pandas as pd
import pytest

from triade.hplc import add_hplc, load_hplc, mark_indeterminate, parse_exclusions
from triade.plates import get_plate_format
from triade.runs import read_run
from triade.samples import add_masses, apply_row_starting_materials, load_samples
from triade.settings import Settings
from triade.synthetic import write_mzml

P384 = get_plate_format("384")


def test_sample_aliases(tmp_path):
    f = tmp_path / "s.csv"
    pd.DataFrame({"Position": ["a01", "B2"], "SMILES": ["CCO", "CCN"], "Code": ["x1", "x2"],
                  "Injection order": [2, 1]}).to_csv(f, index=False)
    df = load_samples(f, P384)
    assert list(df["Well"]) == ["A1", "B2"]
    assert list(df["Sample_Code"]) == ["x1", "x2"]
    assert list(df["Injection_Order"]) == [2, 1]
    assert (df["Plate"] == "1").all()


def test_sample_errors(tmp_path):
    f = tmp_path / "s.csv"
    pd.DataFrame({"Well": ["A1", "Q1"], "Product_SMILES": ["C", "C"]}).to_csv(f, index=False)
    with pytest.raises(ValueError, match="not valid"):
        load_samples(f, P384)
    pd.DataFrame({"Well": ["A1", "A01"], "Product_SMILES": ["C", "C"]}).to_csv(f, index=False)
    with pytest.raises(ValueError, match="more than once"):
        load_samples(f, P384)
    pd.DataFrame({"Well": ["A1"]}).to_csv(f, index=False)
    with pytest.raises(ValueError, match="Product_SMILES"):
        load_samples(f, P384)


def test_excel_plate_sheets_and_sm_by_row(tmp_path):
    f = tmp_path / "s.xlsx"
    with pd.ExcelWriter(f) as xw:
        pd.DataFrame({"Well": ["A1", "B1"], "Product SMILES": ["CCOc1ccccn1", "CCNc1ccccn1"]}).to_excel(
            xw, sheet_name="Plate 1", index=False)
        pd.DataFrame({"Well": ["A1"], "Product SMILES": ["CCOc1ccccn1"]}).to_excel(xw, sheet_name="Plate 2", index=False)
    df = load_samples(f, P384)
    assert sorted(zip(df["Plate"], df["Well"])) == [("1", "A1"), ("1", "B1"), ("2", "A1")]
    rows = tmp_path / "rows.csv"
    pd.DataFrame({"Row": ["A", "B"], "SMILES": ["Clc1ccccn1", "Brc1ccccn1"]}).to_csv(rows, index=False)
    df = apply_row_starting_materials(df, rows, P384)
    assert list(df["SM_SMILES"]) == ["Clc1ccccn1", "Brc1ccccn1", "Clc1ccccn1"]
    m = add_masses(df, ["[M+H]+", "[M+Na]+"], "positive")
    assert m.loc[0, "MH_Theoretical_mz"] == pytest.approx(123.0684 + 1.00728, abs=1e-3)
    assert m.loc[0, "SM_Masses"][0] == pytest.approx(113.00323, abs=1e-4)


def test_exclusions():
    ex = parse_exclusions(["2:9", "1:A3; 1:C"], P384)
    assert ("2", "A9") in ex and ("2", "P9") in ex and ("1", "A3") in ex and ("1", "C24") in ex
    assert len(ex) == 16 + 1 + 24
    with pytest.raises(ValueError):
        parse_exclusions(["9"], P384)


def test_hplc_formats(tmp_path):
    a = tmp_path / "plate 1 conv.csv"
    pd.DataFrame({"Well": ["A1", "A2", "A3"], "Conversion_%": [50, 2, np.nan]}).to_csv(a, index=False)
    b = tmp_path / "presence.csv"
    pd.DataFrame({"Plate": [2, 2], "Well": ["A1", "A2"], "Product": ["yes", "no"]}).to_csv(b, index=False)
    c = tmp_path / "instrument plate 3.csv"
    pd.DataFrame({"Position": ["A1", "A2"], "UV Ratio": [30.0, 0.0], "MS Ratio": [45.0, 0.0],
                  "Product Area": [100.0, 0.0], "SM Area": [50.0, 0.0]}).to_csv(c, index=False)
    h = load_hplc([a, b, c], P384, 5.0)
    get = lambda p, w: h[(h.Plate == p) & (h.Well == w)].iloc[0]
    assert get("1", "A1").HPLC_Conversion == 50
    assert get("2", "A1").HPLC_Conversion == 100 and get("2", "A2").HPLC_Conversion == 0
    assert get("3", "A1").HPLC_Conversion == 45
    assert get("3", "A2").HPLC_No_Signal == 1 and np.isnan(get("3", "A2").HPLC_Conversion)
    h = mark_indeterminate(h, {("1", "A1")})
    assert get("1", "A1").HPLC_Indeterminate == 1


def test_add_hplc_counts_missing_as_not_detected():
    res = pd.DataFrame({"Well": ["A1", "A2", "A3", "A4"], "Detected": [1, 0, 1, 0], "Confidence_Level": [5, 1, 4, 2]})
    hplc = pd.DataFrame({"Plate": "1", "Well": ["A1", "A2", "A3"], "HPLC_Conversion": [80, 0, np.nan],
                         "HPLC_No_Signal": [0, 0, 1], "HPLC_Indeterminate": [0, 0, 0]})
    m = add_hplc(res, hplc, "1", Settings())
    assert list(res["HPLC_Product"]) == [1, 0, 0, 0]
    assert list(res["Agreement"]) == ["both", "neither", "Echo-MS only", "neither"]
    assert m["Agreement_%"] == 75.0


def test_mzml_roundtrip_and_time_units(tmp_path):
    f = tmp_path / "x.mzML"
    spectra = [(np.array([100.0, 200.0]), np.array([1.0, 2.0 + k])) for k in range(5)]
    write_mzml(f, np.arange(5) * 0.5, spectra, "-")
    run = read_run(f)
    assert np.allclose(run.t, np.arange(5) * 0.5, atol=1e-4)
    assert run.polarities == {"-"}
    assert run.tic[3] == pytest.approx(1 + 2 + 3)
    # a file that stores the time in seconds
    text = f.read_text().replace('unitAccession="UO:0000031" unitName="minute"', 'unitAccession="UO:0000010" unitName="second"')
    for k in range(5):
        text = text.replace(f'value="{k * 0.5 / 60:.8f}"', f'value="{k * 0.5:.6f}"', 1)
    g = tmp_path / "y.mzML"
    g.write_text(text)
    assert np.allclose(read_run(g).t, np.arange(5) * 0.5, atol=1e-4)


def test_mzml_does_not_depend_on_pyteomics(tmp_path, monkeypatch):
    """A broken or incomplete pyteomics installation (e.g. without psims) must not stop the reading."""
    from triade import runs

    class Broken:
        def __init__(self, *a, **k):
            raise RuntimeError("Parsing PSI formats requires `psims`.")

    monkeypatch.setattr(runs, "HAVE_PYTEOMICS", True)
    monkeypatch.setattr(runs, "_pyteomics_mzml", types.SimpleNamespace(MzML=Broken), raising=False)
    f = tmp_path / "x.mzML"
    write_mzml(f, np.arange(3) * 0.5, [(np.array([100.0]), np.array([5.0]))] * 3, "+")
    assert len(read_run(f).t) == 3
