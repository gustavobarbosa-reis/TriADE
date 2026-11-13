import zipfile

import pytest

from triade import Inputs, Settings, run_analysis
from triade.cli import main as cli_main


def _levels_ok(res, s):
    """Higher confidence levels must hold a larger share of real products."""
    r = res.results["1"]
    share = {k: r.loc[r.Confidence_Level == k, "HPLC_Product"].mean() for k in (1, 2, 3, 4, 5)
             if (r.Confidence_Level == k).sum() >= 3}
    assert share.get(5, 1) >= 0.9
    assert share.get(1, 0) <= 0.2


def test_96_with_hplc(plate96, tmp_path):
    res = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"], hplc=[plate96["hplc"]]),
                       Settings(plate_format="96"), tmp_path / "out")
    pa = res.plates["1"]
    assert pa.order_name.startswith("row-serpentine")
    assert pa.first_slot == 6                      # 6 lead pulses in the synthetic run
    assert abs(pa.period - 2.5) < 0.005
    r = res.results["1"]
    assert len(r) == 96 and set(r.Confidence_Level) <= {1, 2, 3, 4, 5}
    m = res.summary.iloc[0]
    assert m["Agreement_%"] >= 85
    _levels_ok(res, res.settings)
    for p in ["Summary.csv", "settings.json", "All_results.xlsx", "log.txt",
              "1_detection/Plate1_detection.csv", "1_detection/Plate1_vs_HPLC-MS.jpg",
              "2_confidence_levels/Plate1_confidence_levels.jpg", "2_confidence_levels/Plate1_levels_vs_HPLC-MS.jpg",
              "3_triage/HPLC-MS_triage_savings.csv", "3_triage/HPLC-MS_worklist_scenario_4.csv",
              "4_time_charts/Plate1_time_chart_rows.jpg", "5_QC/Plate1_QC_TIC_ejections.jpg",
              "2_confidence_levels/Plate1_HPLC-MS_conversion.jpg", "5_QC/Plate1_run_QC.csv",
              "6_statistics/Levels_with_product_CI95.csv", "6_statistics/Threshold_calibration.csv",
              "6_statistics/probability_model.json", "7_spectra/Spectra_by_well.pdf"]:
        assert (res.out_dir / p).exists(), p
    assert res.zip_path and zipfile.ZipFile(res.zip_path).testzip() is None
    tri = res.triage.set_index(res.triage.Scenario.str[0])
    assert tri.loc["0", "HPLC_injections"] == 96
    assert tri.loc["4", "HPLC_injections"] == len(res.worklist)
    assert "False_present" in tri.columns


def test_without_hplc(plate96, tmp_path):
    res = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"]), Settings(plate_format="96"),
                       tmp_path / "out", figures=False, make_zip=False)
    assert not res.has_hplc
    assert "Agreement" not in res.results["1"].columns
    assert "False_present" not in res.triage.columns
    assert res.zip_path is None


def test_two_plates_ids_from_file_names(two_plates, tmp_path):
    res = run_analysis(Inputs(mzml=list(reversed(two_plates["mzml"])), samples=two_plates["samples"],
                              hplc=two_plates["hplc"], hplc_exclude=["2:1"]),
                       Settings(plate_format="96"), tmp_path / "out", figures=False, make_zip=False)
    assert sorted(res.results) == ["1", "2"]
    r2 = res.results["2"]
    assert (r2.loc[r2.Column == 1, "HPLC_Result"].str.startswith("indeterminate")).all()
    assert len(res.worklist) == res.triage.loc[4, "HPLC_injections"]


def test_1536_injection_order(plate1536_subset, tmp_path):
    d = plate1536_subset
    res = run_analysis(Inputs(mzml=[d["mzml"]], samples=d["samples"], hplc=[d["hplc"]]),
                       Settings(plate_format="1536"), tmp_path / "out", make_zip=False)
    pa = res.plates["1"]
    assert pa.order_name.startswith("injection order")
    r = res.results["1"]
    assert {"AE", "AF"} <= set(r.Row)
    assert res.summary.iloc[0]["Agreement_%"] >= 80


def test_negative_mode(negative384, tmp_path):
    d = negative384
    s = Settings(polarity="negative", ions=["[M-H]-"])
    res = run_analysis(Inputs(mzml=[d["mzml"]], samples=d["samples"], hplc=[d["hplc"]]), s, tmp_path / "out",
                       figures=False, make_zip=False)
    assert res.summary.iloc[0]["Polarities_acquired"] == "-"
    assert res.summary.iloc[0]["Agreement_%"] >= 80


def test_wrong_polarity_warns(plate96, tmp_path):
    s = Settings(plate_format="96", polarity="negative", ions=["[M-H]-"])
    res = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"]), s, tmp_path / "out",
                       figures=False, make_zip=False)
    assert any("polarity" in w for w in res.warnings)
    assert res.results["1"]["Detected"].sum() == 0


def test_extra_adducts_never_remove_detections(plate96, tmp_path):
    """Adding adducts can only add detections (and therefore false positives), never remove them."""
    base = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"]), Settings(plate_format="96"),
                        tmp_path / "a", figures=False, make_zip=False).results["1"]
    many = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"]),
                        Settings(plate_format="96", ions=["[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[2M+H]+"]),
                        tmp_path / "b", figures=False, make_zip=False).results["1"]
    assert (many["Detected"].to_numpy() >= base["Detected"].to_numpy()).all()


def test_plate_mismatch_is_reported(plate96, tmp_path):
    with pytest.raises(ValueError, match="not in the sample table"):
        run_analysis(Inputs(mzml=[plate96["mzml"], plate96["mzml"]], samples=plate96["samples"],
                            mzml_plates=["1", "7"]), Settings(plate_format="96"), tmp_path / "out",
                     figures=False, make_zip=False)


def test_cli(plate96, tmp_path, capsys):
    rc = cli_main(["--mzml", str(plate96["mzml"]), "--samples", str(plate96["samples"]), "--plate-format", "96",
                   "--hplc", str(plate96["hplc"]), "--cutoff", "300000", "--out", str(tmp_path / "o"),
                   "--no-figures", "--no-zip", "-q"])
    assert rc == 0
    assert Settings.load(tmp_path / "o" / "settings.json").cutoff_cps == 300000
    assert cli_main(["--list-ions"]) == 0
    assert "[M+Na]+" in capsys.readouterr().out
