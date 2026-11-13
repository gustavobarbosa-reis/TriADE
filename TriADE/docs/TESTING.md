# Testing

## 1. Automated tests

Open a terminal in the program folder and activate its environment as in the README (Installation, steps 3–6),
then type:

```
pytest
```

Add `-v` (`pytest -v`) to see the name of each test.

`pytest` must end with `passed` and no `failed`. One test (`test_webapp_with_streamlit`) runs only when
Streamlit is installed; it is reported as `skipped` otherwise.

| File | What it checks |
|---|---|
| `tests/test_plates_adducts.py` | plate formats (96/384/1536, rows A–AF), well names, ejection patterns, ion m/z |
| `tests/test_chem_settings.py` | formulas and masses from SMILES, radical correction, settings save/load and validation |
| `tests/test_inputs.py` | sample tables (CSV/Excel, alternative column names, errors), starting materials by row, HPLC-MS formats, exclusions, mzML reading (minutes and seconds, with or without pyteomics) |
| `tests/test_pipeline.py` | complete runs on synthetic plates: ejection period and order, grading against the known truth, HPLC-MS comparison, triage and worklist, two plates, 1536 wells with a shuffled injection order, negative mode, wrong polarity, the effect of extra adducts, the command line |
| `tests/test_features.py` | Excel template and reagent library, plate-map sheets, preview, not-measurable wells, isotope patterns (Cl, Br₂, Na adduct, 2+ ions), relative M+1 rule, confidence intervals, threshold calibration, probability model (fit, save, apply), run QC with an internal standard, McNemar comparison, every palette and file format, spectra PDF, random access to spectra, the high-resolution mode, new command-line options |
| `tests/test_features.py` (round 3) | replicate ejections (3 per well, joint placement, Fleiss κ, rules), controls and plate acceptance, side-product rules and maps, manual review and the worklist, Excel transfer list and injection sequence, project history and control charts, report (HTML/PDF), publication mode, reagent labels |
| `tests/test_webapp.py` | the browser interface: loads, warns about adducts, runs the synthetic example, preview, results, spectrum viewer (well click, raw spectrum, window maximum), statistics, palettes/formats and the high-resolution mode |

The synthetic plates are written by `triade.synthetic`: made-up SNAr products (one heteroaryl chloride per
row, one amine per column), a reactivity model that sets each well's conversion (the "HPLC-MS" truth), and
an mzML with ejection pulses, product and starting-material ions with isotopes, sodium adducts, noise, and
a few traps (impurities on the product [M+H]+ or [M+Na]+ m/z). `--high-resolution` writes QTOF-like centroid
data (m/z noise 0.4 mDa).

The same tests run on GitHub for Windows, macOS and Linux and Python 3.10–3.13
(`.github/workflows/tests.yml`), with and without RDKit and Streamlit.

## 2. Testing the browser interface by hand

```
triade-app
```

1. Sidebar: plate format **96**. *Inputs* → **Load a synthetic example** → **Run the analysis**.
2. *Preview*: the reagents on the plate grid and the m/z of each ion.
3. *Results*: summary, triage table, worklist, the maps of the plate and the zip download.
4. *Spectrum viewer*: click a well on the map, then a point of the trace; switch between "as analysed" and
   "raw mzML"; check the dashed lines at the calculated m/z; download a PNG.
5. *Statistics*: confidence intervals, threshold calibration, probability model, run QC.
6. *Results* → **Keep this analysis for comparison**; add `[M+Na]+` to the ions (a red warning appears), run
   again and read the comparison table.
7. Sidebar *Figures*: choose **Viridis** and **SVG**, run again, open the zip.
8. *Spectrum viewer*: give a well the decision *Send to HPLC-MS*, save it, and check that it appears in the worklist
   of *Results*; download `review.csv`.
9. *Results*: download the Excel transfer list + sequence and the report (HTML and PDF).
10. Sidebar: type a profile name and **Save profile**; choose it again in *Saved profile*.
11. Advanced → *Replicates, controls and side products*: add `hydrolysis: Reagent_1: -Cl +OH` and run again
    (the synthetic example has P/N/B controls in columns 1, 11 and 12).
12. Advanced → *Figures*: publication mode, round wells; run again and open `Figure_captions.txt` in the zip.

## 3. Validating on your own plates

Before using the triage in routine work, run the program on at least one plate for which you also have
HPLC-MS results, ideally of the same chemistry, plate format and acquisition method:

1. Run with the defaults (`[M+H]+` or `[M-H]-` only) and the HPLC-MS results.
2. Check `5_QC/`: the shaded ejection windows must sit on the TIC pulses, and the ejection order in the log
   must be the one you used (or give `Injection_Order`).
3. Look at `2_confidence_levels/Plate*_levels_vs_HPLC-MS`: the share of wells with product should fall from
   level 5 to level 1. Look at the confidence intervals in `6_statistics/Levels_with_product_CI95.csv`.
4. Read `6_statistics/Threshold_calibration.csv`: the rows "plates except X → plate X" show whether thresholds
   chosen on some plates still work on another.
5. Open `3_triage/HPLC-MS_triage_savings.csv`: choose the scenario whose false presents and missed products
   you can accept, and set it in the settings (`triage_scenario`).
6. Read `3_triage/Discrepancies_EchoMS_vs_HPLC-MS.csv` and look at those wells in the spectrum viewer.
7. If you consider an extra adduct, compare the two analyses (McNemar).
8. Keep the `settings.json` (and `probability_model.json`) of the chosen configuration for the next plates.

Wells without an HPLC-MS result count as not detected, so make sure the HPLC-MS results cover every well.
