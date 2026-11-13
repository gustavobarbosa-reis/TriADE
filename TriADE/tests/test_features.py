"""Reagent library / template input, plate maps, not-measurable wells, isotope patterns, statistics,
palettes and formats, the high-resolution mode and random access to spectra."""
import numpy as np
import openpyxl
import pandas as pd
import pytest

from triade import Inputs, Settings, run_analysis
from triade import isotopes, stats
from triade._mzml import MzML, read_spectrum_at
from triade.adducts import IONS
from triade.figures import PALETTES
from triade.plates import get_plate_format
from triade.runs import inspect_mzml, instrument_check
from triade.samples import load_samples, plate_grid, preview
from triade.synthetic import make_dataset
from triade.templates import write_template

P96 = get_plate_format("96")


# ---------------------------------------------------------------- inputs
def _template_from_library(lib: pd.DataFrame, path, with_hplc=True):
    """The synthetic plate written in the template format: reagent codes + Samples + HPLC-MS sheets."""
    sms = {s: f"SM{i + 1}" for i, s in enumerate(lib["SM_SMILES"].unique())}
    amines = {}
    rows = []
    for r in lib.itertuples(index=False):
        core = r.SM_SMILES[2:]                      # "Cl" + core
        amine = r.Product_SMILES[: -len(core)]
        code = amines.setdefault(amine, f"AM{len(amines) + 1}")
        rows.append({"Plate": r.Plate, "Well": r.Well, "Sample_Code": r.Sample_Code,
                     "Product_SMILES": r.Product_SMILES, "Reagent_1": sms[r.SM_SMILES], "Reagent_2": code})
    with pd.ExcelWriter(path) as xw:
        pd.DataFrame([{"Code": c, "SMILES": s, "Role": "starting material"} for s, c in sms.items()]
                     + [{"Code": c, "SMILES": a + "[H]", "Role": "reagent"} for a, c in amines.items()]
                     ).to_excel(xw, sheet_name="Reagents", index=False)
        pd.DataFrame(rows).to_excel(xw, sheet_name="Samples", index=False)
        if with_hplc:
            lib[["Plate", "Well", "Conversion_%"]].to_excel(xw, sheet_name="HPLC-MS", index=False)
    return path


def test_template_end_to_end(plate96, tmp_path):
    xlsx = _template_from_library(plate96["library"], tmp_path / "template.xlsx")
    df = load_samples(xlsx, P96)
    assert {"Reagent_1", "Reagent_2", "Reagent_1_Code", "Reagent_2_Code"} <= set(df.columns)
    assert df["Reagent_1"].str.startswith("Cl").all()
    assert (df["Background_Key"] == df["Reagent_1"]).all()           # the shared reagent drives the background
    res = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=xlsx), Settings(plate_format="96"),
                       tmp_path / "out", figures=False, make_zip=False)
    assert res.has_hplc                                              # taken from the HPLC-MS sheet
    assert res.summary.iloc[0]["Agreement_%"] >= 85


def test_written_template_is_readable(tmp_path):
    for fmt in ("96", "384", "1536"):
        p = write_template(tmp_path / f"t{fmt}.xlsx", fmt)
        df = load_samples(p, get_plate_format(fmt))
        assert len(df) == 6 and (df["Reagent_1_Code"] != "").all()


def test_plate_map_sheets(tmp_path):
    p = write_template(tmp_path / "t.xlsx", "96")
    wb = openpyxl.load_workbook(p)
    del wb["Samples"]
    wb.create_sheet("Samples").append(["Plate", "Well", "Sample_Code", "Product_SMILES"])
    prod, r1 = wb["Map Product P1"], wb["Map Reagent_1 P1"]
    for col, smi in ((2, "C1COCCN1c1ncccn1"), (3, "CN(C)c1ncccn1")):     # wells A1, A2
        prod.cell(row=2, column=col, value=smi)
        r1.cell(row=2, column=col, value="SM-A")
    wb.save(p)
    df = load_samples(p, P96)
    assert sorted(df["Well"]) == ["A1", "A2"]
    assert (df["Reagent_1"] == "Clc1ncccn1").all()
    assert plate_grid(df, P96, "1", "Reagent_1_Code").loc["A", 1] == "SM-A"


def test_preview_reports_problems(tmp_path):
    f = tmp_path / "s.csv"
    pd.DataFrame({"Well": ["A1", "A2", "A3"], "Product_SMILES": ["CCOc1ncccn1", "C1CC((", ""],
                  "Reagent_1": ["Clc1ncccn1", "NOT-A-CODE", ""]}).to_csv(f, index=False)
    df = load_samples(f, P96)
    table, issues = preview(df, ["[M+H]+"], "positive", {"1": (130.0, 600.0)})
    text = " | ".join(issues)
    assert "A2: the product" in text and "A2: Reagent_1 'NOT-A-CODE'" in text and "A3: no product" in text
    assert "A1: [M+H]+ outside the scanned range" in text            # CCOc1ncccn1 [M+H]+ = 125.07


# ---------------------------------------------------------------- not measurable
def test_not_measurable_goes_to_hplc(tmp_path):
    d = make_dataset(tmp_path, "96", seed=2, mz_range=(200.0, 650.0))
    res = run_analysis(Inputs(mzml=[d["mzml"]], samples=d["samples"], hplc=[d["hplc"]]), Settings(plate_format="96"),
                       tmp_path / "out", figures=False, make_zip=False)
    r = res.results["1"]
    nm = r[r["Confidence_Level"] == 0]
    assert len(nm) > 0 and (nm["MH_Theoretical_mz"] < 200.6).all() and (nm["Call"] == "NOT MEASURABLE").all()
    assert set(nm["Well"]) <= set(res.worklist["Well"])                 # always sent to HPLC-MS
    assert any("not measurable" in w for w in res.warnings)


# ---------------------------------------------------------------- isotopes
def test_isotope_patterns():
    mh = IONS["[M+H]+"]
    p = {k: rel for k, rel, _ in isotopes.pattern("C23H27N5O2", mh)}
    assert p[1] == pytest.approx(0.271, abs=0.01)
    cl = isotopes.pattern("C10H10ClN3", mh)
    k, rel, shift = isotopes.confirmation_peak(cl)
    assert k == 2 and rel == pytest.approx(0.33, abs=0.02) and shift == pytest.approx(1.997, abs=0.01)
    br2 = {k: rel for k, rel, _ in isotopes.pattern("C6H4Br2", IONS["[M+H]+"])}
    assert br2[2] > 1.9                                              # M+2 is the most intense peak
    na = isotopes.ion_composition("C5H5N", IONS["[M+Na]+"])
    assert na == {"C": 5, "H": 5, "N": 1, "Na": 1}
    z2 = isotopes.pattern("C30H40N6", IONS["[M+2H]2+"])
    assert z2[1][2] == pytest.approx(0.5017, abs=0.001)               # isotope spacing 1/z


def test_isotope_columns_and_relative_rule(plate96, tmp_path):
    a = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"]), Settings(plate_format="96"),
                     tmp_path / "a", figures=False, make_zip=False).results["1"]
    det = a[a["Confidence_Level"] >= 4]
    assert det["Isotope_Pattern_Score"].median() > 0.6
    assert det["Isotope_Confirmation_Peak"].str.startswith("M+1").all()
    b = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"]),
                     Settings(plate_format="96", m1_rule="relative"), tmp_path / "b", figures=False,
                     make_zip=False).results["1"]
    assert len(b) == 96 and b["Detected"].sum() >= 0.8 * a["Detected"].sum()


# ---------------------------------------------------------------- statistics
def test_wilson_and_mcnemar():
    assert stats.wilson(36, 37) == (86.2, 99.5)
    assert stats.wilson(0, 10)[0] == 0.0
    m = stats.mcnemar([1, 1, 1, 0, 0, 1], [1, 0, 0, 0, 0, 1])
    assert m["only_A"] == 2 and m["only_B"] == 0 and m["p_value"] == 0.5
    assert stats.auc([3, 2, 1, 0], [1, 1, 0, 0]) == 1.0


def test_statistics_outputs(two_plates, tmp_path):
    res = run_analysis(Inputs(mzml=two_plates["mzml"], samples=two_plates["samples"], hplc=two_plates["hplc"]),
                       Settings(plate_format="96"), tmp_path / "out", figures=False, make_zip=False)
    st = res.stats
    assert set(st["levels"]["Plate"]) == {"1", "2", "all"}
    assert (st["levels"]["CI95_low_%"] <= st["levels"]["With_product_%"]).all()
    cal = st["calibration"]
    assert len(cal["table"]) == 3 and "plate 1" in set(cal["table"]["Tested_on"])
    val = st["probability"]["validation"]
    assert len(val) == 3 and (val["AUC"] > 0.8).all()
    r = res.results["1"]
    assert r["P_Product"].between(0, 1).all()
    assert "Cohen_kappa_CI95" in res.summary.columns
    assert "False_present_CI95_%" in res.triage.columns
    # the saved model applied to a plate without HPLC-MS
    model = res.out_dir / "6_statistics" / "probability_model.json"
    assert model.exists()
    res2 = run_analysis(Inputs(mzml=[two_plates["mzml"][0]], samples=two_plates["samples"], mzml_plates=["1"]),
                        Settings(plate_format="96"), tmp_path / "out2", figures=False, make_zip=False)
    assert "P_Product" not in res2.results["1"].columns
    m = st["probability"]["model"].save(tmp_path / "model.json")
    res3 = run_analysis(Inputs(mzml=[two_plates["mzml"][0]], samples=two_plates["samples"], mzml_plates=["1"],
                               probability_model=m), Settings(plate_format="96"), tmp_path / "out3", figures=False,
                        make_zip=False)
    p = res3.results["1"]["P_Product"]
    assert p.notna().all() and stats.auc(p, two_plates["library"].set_index(["Plate", "Well"]).loc["1"]
                                         .loc[res3.results["1"]["Well"], "Conversion_%"].gt(5).astype(int)) > 0.8


def test_run_qc_internal_standard(plate96, tmp_path):
    res = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"]),
                       Settings(plate_format="96", internal_standard_mz=279.09), tmp_path / "o", figures=False,
                       make_zip=False)
    q = res.stats["run_qc"]["1"]
    assert len(q) == 96 and q["IS_Intensity"].gt(0).all()
    assert "IS_CV_%" in res.summary.columns and "Run_QC_Flag" in res.results["1"].columns


def test_compare_two_analyses(plate96, tmp_path):
    a = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"], hplc=[plate96["hplc"]]),
                     Settings(plate_format="96"), tmp_path / "a", figures=False, make_zip=False)
    b = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"], hplc=[plate96["hplc"]]),
                     Settings(plate_format="96", ions=["[M+H]+", "[M+Na]+"]), tmp_path / "b", figures=False,
                     make_zip=False)
    c = stats.compare(a.results, b.results)
    assert {"McNemar_p_detection", "McNemar_p_accuracy", "Detected_only_B"} <= set(c.columns)


# ---------------------------------------------------------------- figures
@pytest.mark.parametrize("palette,fmt", [(p, f) for p, f in zip(list(PALETTES) + ["custom"],
                                                                 ["jpg", "png", "tiff", "svg", "pdf", "png"])])
def test_palettes_and_formats(plate96, tmp_path, palette, fmt):
    s = Settings(plate_format="96", palette=palette, figure_format=fmt, figure_dpi=150, spectra_pdf="none",
                 custom_colors={"level5": "#ff0000", "seq": "magma"} if palette == "custom" else {})
    res = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"], hplc=[plate96["hplc"]]), s,
                       tmp_path / "o", make_zip=False)
    figs = res.figures["1"]
    assert len(figs) >= 9 and all(f.suffix == f".{fmt}" and f.stat().st_size > 1000 for _, f, _ in figs)
    assert (res.out_dir / f"2_confidence_levels/Plate1_HPLC-MS_conversion.{fmt}").exists()


def test_spectra_pdf_levels(plate96, tmp_path):
    s = Settings(plate_format="96", spectra_pdf="levels", spectra_pdf_levels=[5])
    res = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"]), s, tmp_path / "o", make_zip=False)
    pdf = res.out_dir / "7_spectra" / "Spectra_by_well.pdf"
    n5 = int((res.results["1"]["Confidence_Level"] == 5).sum())
    assert pdf.exists() and pdf.read_bytes().count(b"/Type /Page\n") + pdf.read_bytes().count(b"/Type /Page ") >= n5 > 0


# ---------------------------------------------------------------- mzML access and instrument
def test_random_access_spectrum(plate96):
    with MzML(str(plate96["mzml"])) as r:
        seq = [s for _, s in zip(range(40), r)]
    for i in (0, 17, 39):
        s = read_spectrum_at(plate96["mzml"], i)
        assert np.allclose(s["m/z array"], seq[i]["m/z array"])
        assert np.allclose(s["intensity array"], seq[i]["intensity array"])


def test_high_resolution_mode(tmp_path):
    d = make_dataset(tmp_path, "96", seed=0, high_resolution=True, impurity_rate=0.3, adduct_trap_rate=0.0)
    info = inspect_mzml(d["mzml"])
    assert info["centroid"] and info["instrument_model"] == "synthetic QTOF"
    assert instrument_check(info, "unit") and not instrument_check(info, "high")
    out = {}
    for inst in ("high", "unit"):
        res = run_analysis(Inputs(mzml=[d["mzml"]], samples=d["samples"], hplc=[d["hplc"]]),
                           Settings(plate_format="96", instrument=inst), tmp_path / inst, figures=(inst == "high"),
                           make_zip=False)
        out[inst] = res.summary.iloc[0]
    assert out["high"]["EchoMS_only"] < out["unit"]["EchoMS_only"]     # near-isobaric impurities rejected
    assert out["high"]["Cohen_kappa"] > out["unit"]["Cohen_kappa"]


def test_cli_new_options(plate96, tmp_path):
    from triade.cli import main
    assert main(["--template", str(tmp_path / "t.xlsx"), "--plate-format", "96"]) == 0
    assert (tmp_path / "t.xlsx").exists()
    rc = main(["--mzml", str(plate96["mzml"]), "--samples", str(plate96["samples"]), "--plate-format", "96",
               "--palette", "viridis", "--format", "png", "--dpi", "150", "--spectra-pdf", "none", "--m1-rule",
               "relative", "--out", str(tmp_path / "o"), "--no-zip", "-q"])
    assert rc == 0
    s = Settings.load(tmp_path / "o" / "settings.json")
    assert (s.palette, s.figure_format, s.m1_rule) == ("viridis", "png", "relative")


def test_font_scale(plate96, tmp_path):
    from triade import charts, figures
    with pytest.raises(ValueError):
        Settings(font_scale=3).validate()
    res = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"], hplc=[plate96["hplc"]]),
                       Settings(plate_format="96", font_scale=1.5, spectra_pdf="none"), tmp_path / "o", make_zip=False)
    assert figures.FONT["scale"] == 1.5 and len(res.figures["1"]) >= 9
    spec = charts.plate_map(res.results["1"], P96, figures.Style(font_scale=1.5))
    assert spec["config"]["axis"]["labelFontSize"] == 15
    figures.Style()                                   # back to normal for the other tests
    assert figures.FONT["scale"] == 1.0


# ---------------------------------------------------------------- round 3
def test_replicates_and_controls(tmp_path):
    from triade import replicates
    d = make_dataset(tmp_path, "96", seed=3, replicates=3, controls=True)
    s = pd.read_csv(d["samples"])
    assert s["Replicate"].max() == 3 and (s["Control"].fillna("") != "").sum() == 3 * 24
    res = run_analysis(Inputs(mzml=[d["mzml"]], samples=d["samples"], hplc=[d["hplc"]]),
                       Settings(plate_format="96", replicate_rule="majority"), tmp_path / "o", figures=False,
                       make_zip=False)
    pa = res.plates["1"]
    assert sorted(pa.replicates) == [1, 2, 3]
    assert len({x.first_slot for x in pa.replicates.values()}) == 1
    r = res.results["1"]
    assert len(r) == 96 and (r["Replicates"] == 3).all()
    m = res.summary.iloc[0]
    assert m["Replicate_kappa_Fleiss"] > 0.7 and m["Plate_accepted"] == "yes"
    assert m["Positive_detected_%"] == 100
    assert not set(res.worklist["Well"]) & set(r.loc[r["Control"] != "", "Well"])     # controls not injected
    # the 'all' rule detects no more wells than 'majority'
    res_all = run_analysis(Inputs(mzml=[d["mzml"]], samples=d["samples"]),
                           Settings(plate_format="96", replicate_rule="all"), tmp_path / "o2", figures=False,
                           make_zip=False)
    assert res_all.results["1"]["Detected"].sum() <= r["Detected"].sum()
    # a failing control set
    bad = r.copy()
    bad.loc[bad["Control"] == "negative", "Detected"] = 1
    assert replicates.control_check(bad, Settings())["Plate_accepted"] == "no"


def test_side_products(plate96, tmp_path):
    from triade.sideproducts import parse_rules
    rules = parse_rules(["hydrolysis: Reagent_1: -Cl +OH", "oxidation: Product: +O  # comment"])
    assert rules[0]["delta_mass"] == pytest.approx(-34.96885 + 15.99491 + 1.00783, abs=1e-3)
    assert rules[1]["key"] == "oxidation"
    with pytest.raises(ValueError):
        parse_rules(["bad rule"])
    with pytest.raises(ValueError):
        Settings(side_products=["x: Product: Cl"]).validate()
    res = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"]),
                       Settings(plate_format="96", side_products=["hydrolysis: Reagent_1: -Cl +OH"], spectra_pdf="none"),
                       tmp_path / "o", make_zip=False)
    r = res.results["1"]
    assert {"SP_hydrolysis_mz", "SP_hydrolysis_Detected"} <= set(r.columns)
    assert (res.out_dir / "8_side_products" / "Plate1_hydrolysis.jpg").exists()


def test_review_worklist_and_excel(plate96, tmp_path):
    from triade import review
    from triade.triage import hplc_worklist
    base = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"]), Settings(plate_format="96"),
                        tmp_path / "a", figures=False, make_zip=False)
    r = base.results["1"]
    in_list = base.worklist["Well"].iloc[0]
    out_list = r.loc[r["Confidence_Level"] == 5, "Well"].iloc[0]
    rv = review.set_decision(review.empty(), "1", in_list, "present", "clear")
    rv = review.set_decision(rv, "1", out_list, "hplc", "check")
    rv.to_csv(tmp_path / "rv.csv", index=False)
    res = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"], review=tmp_path / "rv.csv"),
                       Settings(plate_format="96", blank_every=5, transfer_volume_nl=50), tmp_path / "b",
                       figures=False, make_zip=False)
    wl = set(res.worklist["Well"])
    assert in_list not in wl and out_list in wl and len(wl) == len(base.worklist)
    x = res.out_dir / "3_triage" / "HPLC-MS_transfer_and_sequence_scenario_4.xlsx"
    tr = pd.read_excel(x, sheet_name="Transfer list")
    seq = pd.read_excel(x, sheet_name="Injection sequence")
    assert len(tr) == len(wl) and tr["Destination_Well"].iloc[0] == "A1" and tr["Volume_nL"].iloc[0] == 50
    assert (seq["Type"] == "Blank").sum() == (len(wl) - 1) // 5
    with pytest.raises(ValueError):
        review.set_decision(review.empty(), "1", "A1", "maybe")


def test_history_and_control_charts(tmp_path):
    from triade import history
    hist = tmp_path / "history.csv"
    for k in range(7):
        d = make_dataset(tmp_path / f"d{k}", "96", seed=10 + k, plate_id=str(k + 1))
        res = run_analysis(Inputs(mzml=[d["mzml"]], samples=d["samples"], history=hist), Settings(plate_format="96"),
                           tmp_path / f"o{k}", figures=(k == 6), make_zip=False)
    h = pd.read_csv(hist)
    assert len(h) == 7 and {"Detected_%", "Pulse_median"} <= set(h.columns)
    cc = res.control_charts
    assert len(cc) and cc["UCL"].notna().any()
    assert (res.out_dir / "5_QC" / "Control_charts.jpg").exists()
    # a clear outlier is flagged
    h2 = pd.concat([h, h.iloc[[-1]].assign(mzML="x.mzML", Plate="9", **{"Detected_%": 0.0})], ignore_index=True)
    lim = history.control_limits(h2)
    assert (lim[(lim["mzML"] == "x.mzML") & (lim["Metric"] == "Detected_%")]["Flag"] != "").all()


def test_report_and_publication_mode(two_plates, tmp_path):
    from triade import figures
    s = Settings(plate_format="96", figure_width_mm=85, show_titles=False, well_shape="circle", font_family="Arial",
                 spectra_pdf="none")
    res = run_analysis(Inputs(mzml=two_plates["mzml"], samples=two_plates["samples"], hplc=two_plates["hplc"]), s,
                       tmp_path / "o", make_zip=False)
    html = (res.out_dir / "Report.html").read_text()
    assert "Methods" in html and "data:image/png;base64" in html and "TriADE" in html
    assert (res.out_dir / "Report.pdf").stat().st_size > 10000
    caps = (res.out_dir / "Figure_captions.txt").read_text()
    assert "confidence levels" in caps
    assert (res.out_dir / "2_confidence_levels" / "Plate1_levels_and_HPLC-MS_conversion.jpg").exists()
    from PIL import Image
    w = Image.open(res.out_dir / "2_confidence_levels" / "Plate1_confidence_levels.jpg").size[0]
    assert w < 1.4 * 85 / 25.4 * 300            # about one journal column at 300 dpi (plus margins)
    figures.Style()


def test_axis_labels_from_reagent_codes(plate96, tmp_path):
    from triade.pipeline import _axis_labels
    xlsx = _template_from_library(plate96["library"], tmp_path / "t.xlsx", with_hplc=False)
    df = load_samples(xlsx, P96)
    rows, cols = _axis_labels(df, P96)
    assert rows["A"].startswith("SM") and cols[1].startswith("AM")


# ---------------------------------------------------------------- standard method
def test_standard_method_and_non_standard_marking(plate96, tmp_path):
    import json
    from triade.settings import STANDARD_CRITERIA
    assert Settings().is_standard and Settings(ions=["[M+H]+", "[M+Na]+"], palette="viridis").is_standard
    s = Settings(plate_format="96", cutoff_cps=300000, min_spectra=3)
    assert set(s.non_standard()) == {"cutoff_cps", "min_spectra"}
    assert "cutoff_cps" in STANDARD_CRITERIA and "ions" not in STANDARD_CRITERIA
    res = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"]), s, tmp_path / "o",
                       figures=False, make_zip=False)
    assert res.summary.iloc[0]["Method"] == "non-standard"
    assert any(w.startswith("NON-STANDARD PARAMETERS") for w in res.warnings)
    saved = json.loads((tmp_path / "o" / "settings.json").read_text())
    assert saved["_method"].startswith("NON-STANDARD PARAMETERS") and Settings.load(tmp_path / "o" / "settings.json") == s
    std = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"]), Settings(plate_format="96"),
                       tmp_path / "s", figures=False, make_zip=False)
    assert std.summary.iloc[0]["Method"] == "standard" and not any("NON-STANDARD" in w for w in std.warnings)


def test_three_category_map_has_no_question_mark(plate96, tmp_path):
    import matplotlib.pyplot as plt
    from triade import figures
    seen = []
    orig = figures.PlateFigure.cell

    def spy(self, well, color, text=None, *a, **k):
        seen.append(text)
        return orig(self, well, color, text, *a, **k)
    figures.PlateFigure.cell = spy
    try:
        res = run_analysis(Inputs(mzml=[plate96["mzml"]], samples=plate96["samples"]), Settings(plate_format="96"),
                           tmp_path / "o", figures=False, make_zip=False)
        figures.heatmap_call(res.results["1"], P96, tmp_path / "m", "t", "s", figures.Style())
    finally:
        figures.PlateFigure.cell = orig
        plt.close("all")
    assert "?" not in seen and (res.results["1"]["Call"] == "CHECK SPECTRUM").any()
