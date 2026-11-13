"""Command line:  triade --mzml plate1.mzML --samples samples.csv --out results [options]

Run `triade --help` for every option, `triade --list-ions` for the ion catalogue and
`triade-app` for the browser interface."""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .adducts import IONS
from .pipeline import Inputs, run_analysis
from .settings import Settings
from .triage import SCENARIOS


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="triade", description="Product detection in Echo-MS plate runs (mzML), "
                                "with confidence levels and HPLC-MS triage.")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("--list-ions", action="store_true", help="list the available ions/adducts and exit")
    p.add_argument("--template", type=Path, metavar="FILE.xlsx",
                   help="write the Excel input template (for --plate-format) and exit")
    io = p.add_argument_group("inputs")
    io.add_argument("--mzml", nargs="+", type=Path, help="one mzML file per plate")
    io.add_argument("--mzml-plates", nargs="+", help="plate ID of each mzML file (default: from the file name)")
    io.add_argument("--samples", type=Path, help="sample table (CSV/TSV/Excel): Plate, Well, Sample_Code, "
                    "Product_SMILES, SM_SMILES, Injection_Order")
    io.add_argument("--sm-by-row", type=Path, help="optional table Row | SMILES (| Plate) with one starting "
                    "material per plate row")
    io.add_argument("--hplc", nargs="+", type=Path, default=[], help="optional HPLC-MS result files")
    io.add_argument("--hplc-plates", nargs="+", help="plate ID of each HPLC-MS file")
    io.add_argument("--hplc-exclude", nargs="+", default=[], help="HPLC-MS results to treat as indeterminate: "
                    "plate:column, plate:well or plate:row (e.g. 2:9 2:A9 1:C)")
    io.add_argument("--model", type=Path, help="probability_model.json of an earlier run, applied to these plates")
    io.add_argument("--out", type=Path, default=Path("triade_results"), help="output folder")
    st = p.add_argument_group("settings (anything not given here comes from --settings or the defaults)")
    ns = p.add_argument_group("non-standard detection criteria (the standard method fixes them; changing any of "
                              "them marks every output NON-STANDARD PARAMETERS)")
    st.add_argument("--settings", type=Path, help="settings JSON (e.g. settings.json of a previous run)")
    st.add_argument("--instrument", choices=["unit", "high"],
                    help="unit = triple/single quadrupole; high = QTOF/Orbitrap (experimental, centroided data)")
    ns.add_argument("--ppm", type=float, dest="mass_tolerance_ppm", help="m/z tolerance in ppm (high resolution)")
    st.add_argument("--plate-format", choices=["96", "384", "1536"])
    st.add_argument("--polarity", choices=["positive", "negative", "both"])
    st.add_argument("--ions", nargs="+", metavar="ION", help="ions to search, e.g. '[M+H]+' '[M+Na]+'. "
                    "Each extra adduct increases false positives: select only those you expect.")
    ns.add_argument("--cutoff", type=float, dest="cutoff_cps", help="M-peak cutoff in cps (standard 200000)")
    ns.add_argument("--level5", type=float, dest="level5_cps", help="evidence for level 5 Present (cps)")
    ns.add_argument("--level3", type=float, dest="level3_cps", help="evidence for level 3 Check spectrum (cps)")
    ns.add_argument("--tolerance", type=float, dest="mass_tolerance_da", help="m/z tolerance (Da)")
    ns.add_argument("--min-spectra", type=int, dest="min_spectra")
    ns.add_argument("--m1-rule", choices=["absolute", "relative"], dest="m1_rule",
                    help="M+1 check: absolute (>= fraction of the cutoff, default) or relative (vs. the formula)")
    st.add_argument("--period", type=float, dest="ejection_period_s", help="seconds between ejections")
    st.add_argument("--order", dest="ejection_order", help="auto, table or a plate pattern (e.g. row-serpentine)")
    ns.add_argument("--background", dest="background_group", choices=["starting_material", "row", "column", "off"])
    st.add_argument("--hplc-run-min", type=float, dest="hplc_run_min")
    st.add_argument("--hplc-flow", type=float, dest="hplc_flow_ml_min")
    ns.add_argument("--hplc-min-pct", type=float, dest="hplc_product_min_pct")
    st.add_argument("--target", type=float, dest="target_false_present_pct",
                    help="threshold calibration: acceptable errors in %% of wells (default 1)")
    st.add_argument("--internal-standard", type=float, dest="internal_standard_mz", help="internal standard m/z")
    st.add_argument("--palette", choices=["default", "viridis", "cividis", "okabe-ito", "grayscale"])
    st.add_argument("--format", dest="figure_format", choices=["jpg", "png", "tiff", "svg", "pdf"])
    st.add_argument("--dpi", type=int, dest="figure_dpi")
    st.add_argument("--font-scale", type=float, dest="font_scale",
                    help="text size of the figures (1 = normal, 0.8 smaller, 1.5 larger; 0.5-2)")
    st.add_argument("--spectra-pdf", dest="spectra_pdf", choices=["none", "check", "levels", "all"],
                    help="PDF with the spectrum of each well: none, check (level 3, default), levels, all")
    st.add_argument("--spectra-levels", nargs="+", type=int, dest="spectra_pdf_levels",
                    help="levels for --spectra-pdf levels (0-5)")
    st.add_argument("--scenario", type=int, dest="triage_scenario", choices=range(len(SCENARIOS)),
                    help="triage scenario for the HPLC-MS worklist: " + "; ".join(n for n, *_ in SCENARIOS))
    p.add_argument("--no-figures", action="store_true")
    p.add_argument("--no-zip", action="store_true")
    p.add_argument("-q", "--quiet", action="store_true")
    return p


SETTING_ARGS = ["instrument", "mass_tolerance_ppm", "plate_format", "polarity", "ions", "cutoff_cps", "level5_cps",
                "level3_cps", "mass_tolerance_da", "min_spectra", "m1_rule", "ejection_period_s", "ejection_order",
                "background_group", "hplc_run_min", "hplc_flow_ml_min", "hplc_product_min_pct", "triage_scenario",
                "target_false_present_pct", "internal_standard_mz", "palette", "figure_format", "figure_dpi", "font_scale",
                "spectra_pdf", "spectra_pdf_levels"]


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    if a.list_ions:
        for ion in IONS.values():
            print(f"{ion.name:<18} {ion.polarity}  {ion.note}")
        return 0
    if a.template:
        from .templates import write_template
        print(f"Template written: {write_template(a.template, a.plate_format or '384')}")
        return 0
    if not a.mzml or not a.samples:
        build_parser().error("--mzml and --samples are required")
    logging.basicConfig(level=logging.WARNING if a.quiet else logging.INFO, format="%(message)s")
    d = Settings.load(a.settings).to_dict() if a.settings else Settings().to_dict()
    for k in SETTING_ARGS:
        v = getattr(a, k)
        if v is not None:
            d[k] = v
    if a.polarity and not a.ions and not a.settings:
        from .adducts import default_ions
        d["ions"] = default_ions(a.polarity)
    s = Settings.from_dict(d)
    try:
        res = run_analysis(Inputs(mzml=a.mzml, samples=a.samples, mzml_plates=a.mzml_plates,
                                  starting_materials_by_row=a.sm_by_row, hplc=a.hplc, hplc_plates=a.hplc_plates,
                                  hplc_exclude=a.hplc_exclude, probability_model=a.model),
                           s, a.out, figures=not a.no_figures, make_zip=not a.no_zip)
    except (ValueError, FileNotFoundError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    print(f"Detection criteria: {s.method_label()}")
    print(res.summary.to_string(index=False))
    for w in res.warnings:
        print("WARNING:", w, file=sys.stderr)
    print(f"\nResults in {res.out_dir}" + (f" and {res.zip_path}" if res.zip_path else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
