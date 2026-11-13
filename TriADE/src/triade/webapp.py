"""Browser interface (Streamlit). Start it with   triade-app   or   streamlit run src/triade/webapp.py

Everything runs on your own computer: the files you upload are copied to a temporary folder and are
never sent anywhere else."""
from __future__ import annotations

import io
import json
import shutil
import tempfile
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

from triade import __version__, ui
from triade import review as review_mod
from triade.adducts import IONS, PRIMARY, default_ions
from triade.figures import PALETTE_NAMES, PALETTES, Style
from triade.pipeline import Inputs, run_analysis
from triade.plates import get_plate_format
from triade.runs import inspect_mzml, instrument_check, plate_id_from_name
from triade.settings import LEVEL_LABELS, Settings
from triade.triage import SCENARIOS, hplc_worklist

POLARITIES = {"Positive": "positive", "Negative": "negative", "Both (polarity switching)": "both"}
INSTRUMENTS = {"Triple / single quadrupole — unit resolution": "unit",
               "QTOF / Orbitrap — high resolution (experimental)": "high"}
ROLE_NAMES = {"level5": "Level 5 Present", "level4": "Level 4 Likely present", "level3": "Level 3 Check spectrum",
              "level2": "Level 2 Likely absent", "level1": "Level 1 Not detected", "level0": "Not measurable",
              "detected": "Detected", "not_detected": "Not detected", "check": "Check spectrum (3 categories)",
              "both": "Both detected", "neither": "Neither", "echo_only": "Echo-MS only", "ref_only": "HPLC-MS only",
              "line": "Lines (main)", "line2": "Lines (second)", "band": "Shaded bands"}
LV = lambda k: f"{k if k else '×'} {LEVEL_LABELS[k]}"  # noqa: E731

st.set_page_config(page_title="TriADE", page_icon="🧪", layout="wide")
st.markdown(ui.CSS, unsafe_allow_html=True)
state = st.session_state
if "workdir" not in state:
    state.workdir = tempfile.mkdtemp(prefix="triade_")
if "review" not in state:
    state.review = review_mod.empty()
WORK = Path(state.workdir)


def save_upload(upload, sub: str) -> Path:
    d = WORK / sub
    d.mkdir(parents=True, exist_ok=True)
    p = d / Path(upload.name).name
    p.write_bytes(upload.getbuffer())
    return p


def load_example(fmt: str, polarity: str, high: bool) -> None:
    from triade.synthetic import make_dataset
    rows = None if fmt != "1536" else ["A", "B", "AE", "AF"]
    d = make_dataset(WORK / f"example_{fmt}_{polarity}_{high}", fmt, polarity if polarity != "both" else "positive",
                     rows=rows, injection_order=(fmt == "1536"), high_resolution=high, controls=True)
    state.example = {"mzml": [d["mzml"]], "samples": d["samples"], "hplc": [d["hplc"]]}


def example_kit(fmt: str) -> bytes:
    from triade.synthetic import make_dataset
    from triade.templates import write_template
    d = make_dataset(WORK / f"kit_{fmt}", fmt, rows=None if fmt != "1536" else ["A", "B", "AE", "AF"],
                     injection_order=(fmt == "1536"), controls=True)
    tpl = write_template(WORK / f"kit_{fmt}" / "input_template.xlsx", fmt)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for p in (d["mzml"], d["samples"], d["hplc"], tpl):
            z.write(p, Path(p).name)
    return buf.getvalue()


def figure_bytes(fig_obj, fmt: str) -> bytes:
    buf = io.BytesIO()
    fig_obj.savefig(buf, format=fmt, dpi=200, bbox_inches="tight")
    return buf.getvalue()


def vega(spec, key, select=True):
    if select:
        try:
            return st.vega_lite_chart(spec, use_container_width=True, on_select="rerun", key=key)
        except TypeError:   # older Streamlit: no selections
            pass
    st.vega_lite_chart(spec, use_container_width=True)
    return None


def picked(event, name, field):
    try:
        sel = event["selection"][name] if isinstance(event, dict) else event.selection[name]
        if isinstance(sel, list) and sel:
            return sel[0].get(field)
        if isinstance(sel, dict) and sel.get(field):
            v = sel[field]
            return v[0] if isinstance(v, list) else v
    except (KeyError, AttributeError, TypeError, IndexError):
        pass
    return None


def show_error(err: Exception) -> None:
    msg, hint = ui.friendly(err)
    st.error(f"**{msg}**\n\n*How to fix:* {hint}")


# ======================================================================================
# Sidebar
# ======================================================================================
with st.sidebar:
    st.markdown("### Method")
    profiles = ui.list_profiles()
    loaded = Settings()
    pc1, pc2 = st.columns([3, 1])
    with pc1:
        prof = st.selectbox("Saved profile", ["—"] + profiles, help="Named settings kept on this computer.")
    if prof != "—":
        try:
            loaded = Settings.from_dict(ui.load_profile(prof))
        except (ValueError, TypeError, OSError) as e:
            st.error(f"Could not read the profile: {e}")
        with pc2:
            st.write("")
            if st.button("🗑", help="Delete this profile"):
                ui.delete_profile(prof)
                st.rerun()
    up = st.file_uploader("…or load a settings.json", type=["json"], key="settings_upload")
    if up is not None:
        try:
            loaded = Settings.from_dict(json.loads(up.getvalue().decode("utf-8")))
        except (ValueError, TypeError) as e:
            st.error(f"Could not read the settings: {e}")

    st.markdown("### Essentials")
    inst_label = st.radio("Mass spectrometer", list(INSTRUMENTS), index=list(INSTRUMENTS.values()).index(loaded.instrument),
                          help="Unit resolution: 0.1 Da grid, tolerances in Da. High resolution: centroid peaks, "
                               "tolerances in ppm.")
    instrument = INSTRUMENTS[inst_label]
    if instrument == "high":
        st.info("High-resolution mode is **experimental** (tested on synthetic data). The mzML files must be "
                "centroided (msConvert filter *peakPicking*).")
    fmts = ["96", "384", "1536"]
    plate_format = st.selectbox("Plate format (wells)", fmts, index=fmts.index(loaded.plate_format))
    pol_label = st.radio("Ionisation mode", list(POLARITIES), index=list(POLARITIES.values()).index(loaded.polarity),
                         horizontal=False)
    polarity = POLARITIES[pol_label]
    allowed = [n for n, i in IONS.items() if polarity == "both" or i.polarity == {"positive": "+", "negative": "-"}[polarity]]
    default = [i for i in loaded.ions if i in allowed] or default_ions(polarity)
    ions = st.multiselect("Ions / adducts / clusters", allowed, default=default,
                          help="m/z = (n·M + adduct) / |z|. The protonated / deprotonated molecule is the most specific.")
    st.warning(
        "**Select clusters and adducts carefully.** Every extra ion is one more m/z per well where an unrelated "
        "signal — a reagent adduct, a DMSO or solvent cluster, another product, carryover — can pass the criteria. "
        "Extra ions can raise the number of **false positives sharply**, and false positives decide which wells are "
        "*not* sent to HPLC-MS. Add an ion only if you expect it, and check the choice against HPLC-MS on one plate.")
    extra = [i for i in ions if i not in PRIMARY.values()]
    if extra:
        st.error(f"{len(extra)} ion(s) beyond [M+H]+ / [M-H]- selected: {', '.join(extra)}. Wells detected only "
                 "through these ions are graded 3 'Check spectrum', but they still count as detected in the binary map.")
    std = Settings()

    st.markdown("### Advanced")
    with st.expander("Ion list explained"):
        st.dataframe(pd.DataFrame([{"Ion": i.name, "Mode": i.polarity, "Type": i.group, "Note": i.note}
                                   for n, i in IONS.items() if n in allowed]), hide_index=True)
    with st.expander("Detection criteria — standard method (fixed)", expanded=not loaded.is_standard):
        st.caption("The detection and grading criteria are fixed, so that analyses made by different people and on "
                   "different days are comparable. Only the ions above are chosen by the user.")
        crit_rows = [("Cutoff (M peak, level 4)", f"{std.cutoff_cps:,.0f} cps"),
                     ("Level 5 'Present'", f"≥ {std.level5_cps:,.0f} cps"),
                     ("Level 3 'Check spectrum'", f"≥ {std.level3_cps:,.0f} cps"),
                     ("Level 1 'Not detected': raw signal", f"≤ {std.level1_raw_max_cps:,.0f} cps"),
                     ("Failed well with raw signal above → level 3", f"{std.level2_raw_max_cps:,.0f} cps"),
                     ("M+1 confirmation", f"≥ {100 * std.m1_fraction_of_cutoff:g}% of the cutoff"),
                     ("m/z tolerance", f"± {std.mass_tolerance_da:g} Da (unit) / ± {std.mass_tolerance_ppm:g} ppm (high)"),
                     ("M→M+1 shift tolerance", f"± {std.isotope_shift_tolerance_da:g} Da"),
                     ("Same ion in at least", f"{std.min_spectra} spectra"),
                     ("Background control", f"≥ {std.background_ratio:g}× the wells sharing Reagent_1"),
                     ("Carryover → level 3", f"previous ejection ≥ {std.carryover_ratio:g}×"),
                     ("Ejection window", f"{std.window_pre_s:g} s before to {std.window_post_s:g} s after the apex"),
                     ("HPLC-MS product", f"conversion ≥ {std.hplc_product_min_pct:g}%"),
                     ("Replicate ejections", f"detected by the {std.replicate_rule} of replicates"),
                     ("Plate acceptance", f"positive controls ≥ {std.control_min_positive_pct:g}%, "
                                          f"negatives/blanks ≤ {std.control_max_negative_pct:g}% detected")]
        st.dataframe(pd.DataFrame(crit_rows, columns=["Criterion", "Standard value"]), hide_index=True)
        unlock = st.checkbox("Use non-standard criteria (requires validation)", value=not loaded.is_standard,
                             help="Only to adapt the method to another instrument. Every output is then marked "
                                  "NON-STANDARD PARAMETERS with the values changed.")
    src = loaded if unlock else std
    if unlock:
        with st.expander("Non-standard criteria", expanded=True):
            st.error("**Non-standard criteria.** The results, maps, report and settings.json will be marked "
                     "NON-STANDARD PARAMETERS. Validate the criteria against HPLC-MS before using them for triage.")
            cutoff = st.number_input("Cutoff: minimum M-peak intensity (cps) = level 4", min_value=1_000.0,
                                     value=float(src.cutoff_cps), step=10_000.0)
            level5 = st.number_input("Level 5 'Present' from (cps)", min_value=1_000.0, value=float(src.level5_cps),
                                     step=100_000.0)
            level3 = st.number_input("Level 3 'Check spectrum' from (cps)", min_value=0.0, value=float(src.level3_cps),
                                     step=10_000.0)
            level1 = st.number_input("Level 1 'Not detected': raw signal up to (cps)", min_value=0.0,
                                     value=float(src.level1_raw_max_cps), step=10_000.0)
            level2 = st.number_input("Raw signal above which a failed well is level 3 (cps)", min_value=0.0,
                                     value=float(src.level2_raw_max_cps), step=50_000.0)
            rules_m1 = ["absolute", "relative"]
            m1_rule = st.radio("M+1 confirmation", rules_m1, index=rules_m1.index(src.m1_rule),
                               format_func=lambda r: {"absolute": "≥ fraction of the cutoff",
                                                      "relative": "vs. the M+1 expected from the formula"}[r])
            m1_frac = st.number_input("M+1 intensity, fraction of the cutoff", 0.0, 1.0,
                                      float(src.m1_fraction_of_cutoff), 0.05)
            m1_rel = st.number_input("M+1 vs. expected: at least this fraction", 0.05, 1.5,
                                     float(src.m1_relative_fraction), 0.05)
            m1_floor = st.number_input("M+1 vs. expected: minimum intensity (cps)", 0.0, 1e6, float(src.m1_floor_cps),
                                       5_000.0)
            tol = st.number_input("m/z tolerance (± Da, unit resolution)", 0.05, 2.0, float(src.mass_tolerance_da), 0.05)
            ppm = st.number_input("m/z tolerance (± ppm, high resolution)", 0.5, 100.0, float(src.mass_tolerance_ppm), 0.5)
            iso_tol = st.number_input("M→M+1 shift tolerance (± Da)", 0.05, 1.0, float(src.isotope_shift_tolerance_da),
                                      0.05)
            min_spectra = int(st.number_input("Same ion in at least N spectra", 1, 10, int(src.min_spectra), 1))
            bg_opts = ["starting_material", "row", "column", "off"]
            bg_group = st.selectbox("Background control: compare with wells sharing the same", bg_opts,
                                    index=bg_opts.index(src.background_group),
                                    format_func=lambda o: {"starting_material": "shared reagent (Reagent_1)",
                                                           "row": "plate row", "column": "plate column", "off": "off"}[o])
            bg_ratio = st.number_input("… and require at least this × their median signal", 0.0, 20.0,
                                       float(src.background_ratio), 0.5)
            carry = st.number_input("Carryover: previous ejection ≥ this × the well's signal → level 3", 0.0, 50.0,
                                    float(src.carryover_ratio), 0.5)
            pre = st.number_input("Window before the apex (s)", 0.05, 5.0, float(src.window_pre_s), 0.05)
            post = st.number_input("Window after the apex (s)", 0.05, 5.0, float(src.window_post_s), 0.05)
            min_pct = st.number_input("HPLC-MS conversion that counts as product (%)", 0.0, 100.0,
                                      float(src.hplc_product_min_pct), 1.0)
            rr = ["any", "majority", "all"]
            rep_rule = st.radio("Replicate ejections: a well is detected when", rr, index=rr.index(src.replicate_rule),
                                format_func=lambda r: {"any": "any replicate detects it", "majority": "most replicates do",
                                                       "all": "all replicates do"}[r])
            cmin = st.number_input("Plate accepted if positive controls detected ≥ (%)", 0.0, 100.0,
                                   float(src.control_min_positive_pct), 5.0)
            cmax = st.number_input("… and negative controls / blanks detected ≤ (%)", 0.0, 100.0,
                                   float(src.control_max_negative_pct), 5.0)
    else:
        cutoff, level5, level3 = std.cutoff_cps, std.level5_cps, std.level3_cps
        level1, level2, m1_rule = std.level1_raw_max_cps, std.level2_raw_max_cps, std.m1_rule
        m1_frac, m1_rel, m1_floor = std.m1_fraction_of_cutoff, std.m1_relative_fraction, std.m1_floor_cps
        tol, ppm, iso_tol = std.mass_tolerance_da, std.mass_tolerance_ppm, std.isotope_shift_tolerance_da
        min_spectra, bg_group, bg_ratio, carry = std.min_spectra, std.background_group, std.background_ratio, \
            std.carryover_ratio
        pre, post, min_pct = std.window_pre_s, std.window_post_s, std.hplc_product_min_pct
        rep_rule, cmin, cmax = std.replicate_rule, std.control_min_positive_pct, std.control_max_negative_pct
    with st.expander("Ejection timing"):
        period = st.number_input("Seconds between ejections (0 = from the file name 'delay NNNN' or the data)",
                                 0.0, 30.0, float(loaded.ejection_period_s or 0.0), 0.1)
        order_opts = ["auto", "table", "row-serpentine", "row-major", "column-serpentine", "column-major"]
        order = st.selectbox("Ejection order", order_opts,
                             index=order_opts.index(loaded.ejection_order) if loaded.ejection_order in order_opts else 0,
                             help="auto: the Injection_Order column when present, otherwise the best-fitting pattern.")
    with st.expander("Side products"):
        sp_text = st.text_area("Side products to search (one per line: name: Product or Reagent_1: change)",
                               "\n".join(loaded.side_products),
                               placeholder="hydrolysis: Reagent_1: -Cl +OH\noxidation: Product: +O")
    with st.expander("HPLC-MS, triage and export"):
        run_min = st.number_input("HPLC-MS run time per injection (min)", 0.5, 120.0, float(loaded.hplc_run_min), 0.5)
        flow = st.number_input("HPLC-MS flow rate (mL/min)", 0.01, 10.0, float(loaded.hplc_flow_ml_min), 0.05)
        scen_names = [n for n, *_ in SCENARIOS]
        scenario = scen_names.index(st.selectbox("Triage scenario for the HPLC-MS worklist", scen_names,
                                                 index=int(loaded.triage_scenario)))
        tpf = st.selectbox("Destination plate of the transfer list", ["96", "384"],
                           index=["96", "384"].index(loaded.transfer_plate_format)
                           if loaded.transfer_plate_format in ("96", "384") else 0)
        tvol = st.number_input("Transfer volume (nL, 0 = not written)", 0.0, 100000.0, float(loaded.transfer_volume_nl or 0), 2.5)
        blank_every = int(st.number_input("Injection sequence: a blank every N samples (0 = none)", 0, 1000,
                                          int(loaded.blank_every), 1))
    with st.expander("Statistics and QC"):
        target = st.number_input("Threshold calibration: acceptable errors (% of wells)", 0.1, 20.0,
                                 float(loaded.target_false_present_pct), 0.1)
        is_mz = st.number_input("Internal standard m/z (0 = none)", 0.0, 5000.0, float(loaded.internal_standard_mz or 0.0),
                                0.1, format="%.4f")
        hist_path = st.text_input("Project history file (on this computer, for control charts across plates)",
                                  state.get("hist_path", ""), placeholder=str(Path.home() / "triade_history.csv"),
                                  help="A CSV that gains one row per plate at every run; leave empty to skip.")
        state.hist_path = hist_path
    with st.expander("Figures"):
        pal_keys = list(PALETTE_NAMES)
        palette = st.selectbox("Colour palette", pal_keys, index=pal_keys.index(loaded.palette),
                               format_func=lambda k: PALETTE_NAMES[k])
        custom_colors = dict(loaded.custom_colors or {})
        if palette == "custom":
            base = PALETTES["default"]
            for role, label in ROLE_NAMES.items():
                custom_colors[role] = st.color_picker(label, custom_colors.get(role, base[role]), key=f"col_{role}")
            custom_colors["seq"] = st.selectbox("Intensity / conversion scale", ["viridis", "cividis", "magma", "Blues",
                                                                                  "Greys", "YlGnBu"], index=0)
        fmts_out = ["jpg", "png", "tiff", "svg", "pdf"]
        fig_fmt = st.selectbox("File format", fmts_out, index=fmts_out.index(loaded.figure_format),
                               help="SVG and PDF are vector formats (best for papers); TIFF is often asked by journals.")
        dpi = int(st.select_slider("Resolution (dpi)", [150, 300, 600],
                                   value=loaded.figure_dpi if loaded.figure_dpi in (150, 300, 600) else 300))
        fonts = [0.8, 0.9, 1.0, 1.1, 1.25, 1.5, 1.75, 2.0]
        font_scale = float(st.select_slider("Font size", fonts, value=min(fonts, key=lambda f: abs(f - loaded.font_scale)),
                                            format_func=lambda f: f"{f:g}×"))
        fams = ["default", "Arial", "Helvetica", "Times New Roman", "DejaVu Sans"]
        font_family = st.selectbox("Font", fams, index=fams.index(loaded.font_family) if loaded.font_family in fams else 0)
        well_shape = st.radio("Wells", ["square", "circle"], index=["square", "circle"].index(loaded.well_shape),
                              horizontal=True)
        cell_values = st.checkbox("Values inside the cells (intensity and conversion maps)", value=loaded.cell_values)
        axis_labels = st.checkbox("Reagent codes along the plate axes", value=loaded.axis_reagent_labels)
        highlight = st.checkbox("Outline discrepancies, reviewed wells, controls and QC flags", value=loaded.highlight_wells)
        pub = st.checkbox("Publication mode (journal width, no titles — captions in a text file)",
                          value=loaded.figure_width_mm is not None or not loaded.show_titles)
        if pub:
            widths = {"One column (85 mm)": 85.0, "1.5 columns (120 mm)": 120.0, "Two columns (180 mm)": 180.0}
            wl = st.radio("Figure width", list(widths), index=2 if loaded.figure_width_mm in (None, 180.0) else
                          list(widths.values()).index(loaded.figure_width_mm) if loaded.figure_width_mm in widths.values() else 2)
            width_mm, show_titles = widths[wl], False
        else:
            width_mm, show_titles = None, True
        pdf_opts = ["none", "check", "levels", "all"]
        spectra_pdf = st.selectbox("PDF with the spectrum of each well", pdf_opts, index=pdf_opts.index(loaded.spectra_pdf),
                                   format_func=lambda o: {"none": "no PDF", "check": "level 3 'Check spectrum' only",
                                                          "levels": "levels chosen below", "all": "all wells (large)"}[o])
        pdf_levels = st.multiselect("Levels in the PDF", [5, 4, 3, 2, 1, 0], default=list(loaded.spectra_pdf_levels),
                                    format_func=LV) if spectra_pdf == "levels" else list(loaded.spectra_pdf_levels)
        report_on = st.checkbox("Write a report (HTML and PDF)", value=loaded.report)

settings_error = None
try:
    settings = Settings(
        instrument=instrument, plate_format=plate_format, polarity=polarity, ions=ions,
        ejection_period_s=period or None, ejection_order=order, mass_tolerance_da=tol, mass_tolerance_ppm=ppm,
        isotope_shift_tolerance_da=iso_tol, min_spectra=min_spectra, cutoff_cps=cutoff, m1_fraction_of_cutoff=m1_frac,
        m1_rule=m1_rule, m1_relative_fraction=m1_rel, m1_floor_cps=m1_floor, background_group=bg_group,
        background_ratio=bg_ratio, level5_cps=level5, level3_cps=level3, level1_raw_max_cps=level1,
        level2_raw_max_cps=level2, carryover_ratio=carry, window_pre_s=pre, window_post_s=post,
        hplc_product_min_pct=min_pct, hplc_run_min=run_min, hplc_flow_ml_min=flow, triage_scenario=scenario,
        target_false_present_pct=target, internal_standard_mz=is_mz or None, replicate_rule=rep_rule,
        control_min_positive_pct=cmin, control_max_negative_pct=cmax,
        side_products=[x for x in sp_text.splitlines() if x.strip()], transfer_plate_format=tpf,
        transfer_volume_nl=tvol or None, blank_every=blank_every, palette=palette,
        custom_colors=custom_colors if palette == "custom" else {}, figure_format=fig_fmt, figure_dpi=dpi,
        cell_values=cell_values, font_scale=font_scale, font_family=font_family, well_shape=well_shape,
        axis_reagent_labels=axis_labels, highlight_wells=highlight, figure_width_mm=width_mm, show_titles=show_titles,
        spectra_pdf=spectra_pdf, spectra_pdf_levels=pdf_levels or [3], report=report_on).validate()
except ValueError as e:
    settings, settings_error = None, str(e)

with st.sidebar:
    st.markdown("### Save these settings")
    pname = st.text_input("Profile name", placeholder="e.g. SNAr · DMSO · positive · 384")
    if st.button("Save profile", disabled=settings is None or not pname.strip()):
        try:
            ui.save_profile(pname, settings)
            st.success(f"Saved: {pname}")
        except (ValueError, OSError) as e:
            st.error(str(e))
    if settings is not None:
        st.download_button("Download settings.json", json.dumps(settings.to_dict(), indent=2), "settings.json",
                           "application/json")

plate = get_plate_format(plate_format)
res = state.get("result")

# ======================================================================================
# Header and status of the steps
# ======================================================================================
st.markdown(ui.header_html(__version__), unsafe_allow_html=True)

samples_ok = state.get("samples_ok", False)
n_issues = state.get("n_issues", 0)
tabs = st.tabs([
    ":material/upload_file: 1 · Inputs" + (" ✓" if state.get("inputs_ready") else ""),
    ":material/fact_check: 2 · Preview" + (f" ⚠ {n_issues}" if n_issues else (" ✓" if samples_ok else "")),
    ":material/grid_on: 3 · Results" + (" ✓" if res is not None else ""),
    ":material/troubleshoot: 4 · Spectrum viewer",
    ":material/query_stats: 5 · Statistics",
    ":material/help: Help"])
tab_in, tab_prev, tab_res, tab_view, tab_stat, tab_help = tabs

# --------------------------------------------------------------------------------------- inputs
with tab_in:
    c0, c00 = st.columns(2)
    with c0:
        from triade.templates import write_template
        tpl = write_template(WORK / f"input_template_{plate_format}.xlsx", plate_format)
        st.download_button("Download the input template (Excel)", tpl.read_bytes(), "triade_input_template.xlsx",
                           "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", icon=":material/download:")
        st.caption("Sheets: Instructions, Reagents (library), Samples (one row per well; Control and Replicate "
                   "optional), Map … (plate-shaped grids), HPLC-MS. The yellow rows are an example.")
    with c00:
        if st.button("Prepare a complete synthetic example (mzML + tables)", icon=":material/science:"):
            state.kit = example_kit(plate_format)
        if state.get("kit"):
            st.download_button("Download the example files (.zip)", state.kit, f"triade_example_{plate_format}.zip",
                               "application/zip", icon=":material/folder_zip:")
        st.caption("Made-up molecules and spectra, to see what each file looks like.")

    st.markdown("#### Required")
    c1, c2 = st.columns(2)
    with c1:
        mz_up = st.file_uploader("mzML files (one per plate)", type=["mzml", "mzML"], accept_multiple_files=True)
        st.caption("Convert with ProteoWizard msConvert (MS1, no MS-Numpress). 'Plate 2' and 'delay 2500' in the name "
                   "set the plate ID and the 2.5 s interval. High-resolution data: centroided.")
    with c2:
        samples_up = st.file_uploader("Sample table: the Excel template, or a CSV/TSV/Excel list",
                                      type=["csv", "tsv", "txt", "xlsx", "xls"])
        st.caption("Plate, Well, Sample_Code, Product_SMILES, Reagent_1, Reagent_2, …, Injection_Order, Control, "
                   "Replicate (the last three optional).")
    mz_plate_ids = []
    if mz_up:
        st.markdown("**Plate ID of each mzML file** (must match the Plate column of the sample table)")
        cols = st.columns(min(len(mz_up), 4))
        for k, f in enumerate(mz_up):
            with cols[k % len(cols)]:
                mz_plate_ids.append(st.text_input(f.name, value=plate_id_from_name(f.name) or str(k + 1),
                                                  key=f"mzpid_{k}_{f.name}"))

    with st.expander("Optional inputs", expanded=False):
        c3, c4 = st.columns(2)
        with c3:
            sm_up = st.file_uploader("Starting materials by plate row (Row, SMILES[, Plate])",
                                     type=["csv", "tsv", "txt", "xlsx"])
            hplc_up = st.file_uploader("HPLC-MS results (if not in the template's HPLC-MS sheet)",
                                       type=["csv", "tsv", "txt", "xlsx"], accept_multiple_files=True)
            excl = st.text_input("HPLC-MS results to treat as indeterminate", "",
                                 help="plate:column, plate:well or plate:row — e.g. 2:9, 1:A3, 1:C")
        with c4:
            model_up = st.file_uploader("Probability model (probability_model.json of an earlier run)", type=["json"])
            review_up = st.file_uploader("Review decisions (review.csv from the spectrum viewer)", type=["csv", "xlsx"])
            if review_up is not None and state.get("review_loaded") != review_up.name:
                try:
                    state.review = review_mod.load(save_upload(review_up, "review"))
                    state.review_loaded = review_up.name
                except ValueError as e:
                    show_error(e)
        hplc_plate_ids = []
        if hplc_up:
            cols = st.columns(min(len(hplc_up), 4))
            for k, f in enumerate(hplc_up):
                with cols[k % len(cols)]:
                    hplc_plate_ids.append(st.text_input(f"Plate of {f.name}", value=plate_id_from_name(f.name) or "",
                                                        key=f"hppid_{k}_{f.name}"))

    use_example = not mz_up and not samples_up and bool(state.get("example"))
    ready = (bool(mz_up) and samples_up is not None) or use_example
    state.inputs_ready = ready
    st.markdown(ui.checklist([
        (bool(mz_up) or use_example, "mzML file(s) uploaded"),
        (samples_up is not None or use_example, "sample table uploaded"),
        (settings is not None, "settings valid" + (f" — {settings_error}" if settings_error else "")),
        (state.get("n_issues", 0) == 0, "no problems in the Preview" if not state.get("n_issues") else
         f"{state.get('n_issues')} problem(s) in the Preview — check them first")]), unsafe_allow_html=True)
    b1, b2, b3 = st.columns([1.2, 1.2, 2])
    with b1:
        if st.button("Load a synthetic example", icon=":material/experiment:"):
            with st.spinner("Writing a synthetic plate…"):
                load_example(plate_format, polarity, instrument == "high")
            st.rerun()
    with b2:
        make_figs = st.checkbox("Draw the figures", value=True)
    if state.get("example") and use_example:
        st.info("Synthetic example loaded (made-up molecules and spectra); it is used because no files are uploaded.")
    if st.button("Run the analysis", type="primary", disabled=not ready or settings is None, icon=":material/play_arrow:"):
        try:
            if use_example:
                ex = state.example
                inputs = Inputs(mzml=ex["mzml"], samples=ex["samples"], hplc=ex["hplc"])
            else:
                inputs = Inputs(
                    mzml=[save_upload(f, "mzml") for f in mz_up], mzml_plates=mz_plate_ids,
                    samples=save_upload(samples_up, "samples"),
                    starting_materials_by_row=save_upload(sm_up, "samples") if sm_up else None,
                    hplc=[save_upload(f, "hplc") for f in hplc_up or []],
                    hplc_plates=[h or None for h in hplc_plate_ids] or None,
                    hplc_exclude=[excl] if excl.strip() else [],
                    probability_model=save_upload(model_up, "model") if model_up else None)
            if len(state.review):
                rp = WORK / "review" / "review.csv"
                rp.parent.mkdir(exist_ok=True)
                state.review.to_csv(rp, index=False)
                inputs.review = rp
            if state.get("hist_path", "").strip():
                inputs.history = Path(state.hist_path.strip()).expanduser()
            shutil.rmtree(WORK / "results", ignore_errors=True)
            bar = st.progress(0.0, text="Starting")
            t0 = time.time()
            state.result = run_analysis(inputs, settings, WORK / "results", figures=make_figs, make_zip=True,
                                        progress=lambda f, m: bar.progress(min(max(f, 0.0), 1.0),
                                                                           text=m + ui.eta_text(f, time.time() - t0)))
            st.success(f"Done in {time.time() - t0:.0f} s — open Results, Spectrum viewer and Statistics.")
            res = state.result
        except (ValueError, FileNotFoundError, NotImplementedError, KeyError) as e:
            show_error(e)

# --------------------------------------------------------------------------------------- preview
with tab_prev:
    from triade.samples import apply_row_starting_materials, load_samples, plate_grid, preview, reagent_columns
    src_samples = save_upload(samples_up, "samples") if samples_up is not None else \
        (state.example["samples"] if state.get("example") else None)
    if src_samples is None:
        st.info("Upload the sample table (or load the synthetic example) to see what the program understood.")
        state.samples_ok, state.n_issues = False, 0
    else:
        try:
            df = load_samples(src_samples, plate)
            if sm_up is not None:
                df = apply_row_starting_materials(df, save_upload(sm_up, "samples"), plate)
            ranges, checks = {}, []
            mz_files = [(save_upload(f, "mzml"), pid) for f, pid in zip(mz_up or [], mz_plate_ids)] if mz_up else \
                [(Path(p), "1") for p in (state.example["mzml"] if use_example else [])]
            for path, pid in mz_files:
                info = inspect_mzml(path)
                ranges[str(pid)] = (info["mz_min"], info["mz_max"])
                checks += instrument_check(info, instrument)
                st.caption(f"📄 {path.name}: plate {pid}, m/z {info['mz_min']:.1f}–{info['mz_max']:.1f}, "
                           f"{'positive' if info['polarities'] == {'+'} else 'negative' if info['polarities'] == {'-'} else 'both polarities'}"
                           f", {'centroid' if info['centroid'] else 'profile' if info['centroid'] is False else 'unknown'} data"
                           f"{', ' + info['instrument_model'] if info['instrument_model'] else ''}")
            for c in checks:
                st.warning(c)
            table, issues = preview(df, ions, polarity, ranges or None)
            state.samples_ok, state.n_issues = True, len(issues)
            k = st.columns(5)
            for col, (lab, val) in zip(k, [("Wells", df.drop_duplicates(["Plate", "Well"]).shape[0]),
                                            ("Plates", df["Plate"].nunique()),
                                            ("Reagent columns", len(reagent_columns(df))),
                                            ("Controls", int((df["Control"] != "").sum()) if "Control" in df else 0),
                                            ("Replicates", int(df["Replicate"].max()) if "Replicate" in df else 1)]):
                col.markdown(ui.card(lab, val), unsafe_allow_html=True)
            if issues:
                st.error(f"{len(issues)} problem(s) found")
                st.dataframe(pd.DataFrame({"Problem": issues}), hide_index=True)
            else:
                st.success("No problems found in the sample table.")
            pid = st.selectbox("Plate", sorted(df["Plate"].unique()), key="prev_plate")
            fields = [c for c in ["Sample_Code", "Formula"] + [f"{c}_Code" for c in reagent_columns(df)]
                      + reagent_columns(df) + ["Control", "Injection_Order"] if c in table.columns or c in df.columns]
            field_ = st.selectbox("Show on the plate", fields, key="prev_field")
            src = table if field_ in table.columns else df[df["Replicate"] == 1] if "Replicate" in df else df
            st.dataframe(plate_grid(src, plate, pid, field_))
            st.markdown("**Every well: formula and m/z of each selected ion**")
            st.dataframe(table[table["Plate"] == pid], hide_index=True)
        except ValueError as e:
            state.samples_ok = False
            show_error(e)

# --------------------------------------------------------------------------------------- results
if res is not None:
    review_mod.apply(res.results, state.review)
with tab_res:
    if res is None:
        st.info("Run the analysis first.")
    else:
        s_ = res.settings
        for w in res.warnings:
            st.warning(w)
        work = hplc_worklist(res.results, s_.triage_scenario)
        n_all = sum(len(r) for r in res.results.values())
        tri4 = res.triage.iloc[s_.triage_scenario]
        cards = [("Wells", n_all, f"{len(res.results)} plate(s)"),
                 ("Level 5 present", int(sum((r["Confidence_Level"] == 5).sum() for r in res.results.values())), ""),
                 ("To HPLC-MS", len(work), f"scenario {s_.triage_scenario}"),
                 ("Time saved", f"{(n_all - len(work)) * s_.hplc_run_min / 60:.1f} h",
                  f"{(n_all - len(work)) * s_.hplc_run_min * s_.hplc_flow_ml_min:,.0f} mL mobile phase")]
        if res.has_hplc:
            cards.append(("Agreement / κ", f"{res.summary['Agreement_%'].mean():.1f}%",
                          "κ " + ", ".join(f"{v:.2f}" for v in res.summary["Cohen_kappa"])))
        if "Plate_accepted" in res.summary.columns:
            acc = res.summary["Plate_accepted"].fillna("–")
            cards.append(("Controls", "OK" if (acc == "yes").all() else "FAILED",
                          ", ".join(f"P{p}: {a}" for p, a in zip(res.summary["Plate"], acc))))
        cols = st.columns(len(cards))
        for col, (lab, val, sub) in zip(cols, cards):
            col.markdown(ui.card(lab, val, sub), unsafe_allow_html=True)
        st.write("")
        d1, d2, d3, d4 = st.columns(4)
        if res.zip_path and Path(res.zip_path).exists():
            d1.download_button("All results (.zip)", Path(res.zip_path).read_bytes(), Path(res.zip_path).name,
                               "application/zip", type="primary", icon=":material/folder_zip:")
        for col, name, mime in ((d2, "Report.html", "text/html"), (d3, "Report.pdf", "application/pdf")):
            p = Path(res.out_dir) / name
            if p.exists():
                col.download_button(name, p.read_bytes(), name, mime, icon=":material/description:")
        xb = io.BytesIO()
        review_mod.write_excel(xb, work, s_, SCENARIOS[s_.triage_scenario][0])
        d4.download_button("HPLC-MS transfer list + sequence (Excel)", xb.getvalue(),
                           f"HPLC-MS_transfer_and_sequence_scenario_{s_.triage_scenario}.xlsx",
                           "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", icon=":material/table:")

        st.markdown("#### Plates")
        pid = st.selectbox("Plate", list(res.results), key="res_plate")
        r = res.results[pid]
        figs = res.figures.get(pid, [])
        shown = [(t, prev) for t, f, prev in figs if prev is not None and Path(prev).exists()]
        if shown:
            ftabs = st.tabs([t for t, _ in shown])
            for tab, (t, f) in zip(ftabs, shown):
                with tab:
                    st.image(str(f))
        elif figs:
            st.info(f"Figures saved as {s_.figure_format.upper()} (in the zip); this format cannot be shown in the browser.")
        with st.expander("Summary table of the plates"):
            st.dataframe(res.summary, hide_index=True)
        st.markdown("#### HPLC-MS triage")
        st.caption("Which confidence levels are trusted without HPLC-MS; not-measurable wells and controls are handled "
                   "apart. " + ("False present = trusted as present but no product in HPLC-MS; products missed = trusted "
                                "as absent but product in HPLC-MS (with 95% confidence intervals)." if res.has_hplc else
                                "Upload HPLC-MS results to see the errors each scenario makes."))
        st.dataframe(res.triage, hide_index=True)
        st.markdown(f"**Worklist — {SCENARIOS[s_.triage_scenario][0]}: {len(work)} injections** "
                    f"(manual review decisions included: {len(state.review)})")
        st.dataframe(work, hide_index=True)
        st.markdown("#### Wells")
        lv = st.multiselect("Show confidence levels", [5, 4, 3, 2, 1, 0], default=[5, 4, 3, 2, 1, 0], format_func=LV)
        cols_w = ["Well", "Sample_Code", "Control", "Confidence_Level", "Confidence_Label", "Call", "Review_Decision",
                  "P_Product", "Level_Note", "Evidence_Intensity", "Raw_Ejection_Intensity", "Ions_Detected",
                  "Replicate_Levels", "Isotope_Theoretical", "Isotope_Observed", "Isotope_Pattern_Score", "SM_Interference",
                  "Echo_Product_Ratio_%", "Run_QC_Flag", "HPLC_Conversion", "HPLC_Result", "Agreement"] + \
            [c for c in r.columns if c.startswith("SP_") and c.endswith("_Detected")]
        st.dataframe(r.loc[r["Confidence_Level"].isin(lv), [c for c in cols_w if c in r.columns]], hide_index=True)
        if not res.discrepancies.empty:
            with st.expander(f"Discrepancies with {s_.reference_label} ({len(res.discrepancies)})"):
                st.dataframe(res.discrepancies, hide_index=True)
        st.markdown("#### Compare with another analysis")
        st.caption("Keep this analysis, change a setting, run again: changed calls and McNemar's test.")
        if st.button("Keep this analysis for comparison", icon=":material/compare_arrows:"):
            state.kept = {"results": {k: v.copy() for k, v in res.results.items()},
                          "label": f"{', '.join(s_.ions)}; cutoff {s_.cutoff_cps:,.0f}"}
        kept = state.get("kept")
        if kept:
            from triade.stats import compare
            st.caption(f"A = kept ({kept['label']}); B = current ({', '.join(s_.ions)}; cutoff {s_.cutoff_cps:,.0f})")
            st.dataframe(compare(kept["results"], res.results), hide_index=True)

# --------------------------------------------------------------------------------------- spectrum viewer
with tab_view:
    if res is None:
        st.info("Run the analysis first.")
    else:
        from triade import charts, spectra
        from triade.figures import well_page
        s_ = res.settings
        style = Style.from_settings(s_)
        top = st.columns([1, 1, 2])
        with top[0]:
            vp = st.selectbox("Plate", list(res.plates), key="view_plate")
        pa_main = res.plates[vp]
        reps = sorted(pa_main.replicates) if pa_main.replicates else []
        rep = None
        if reps:
            with top[1]:
                rep = st.selectbox("Replicate", reps, key="view_rep")
        pa = pa_main.replicates[rep] if rep is not None else pa_main
        r = pa_main.results
        wells = list(r["Well"])
        if state.get("view_well") not in wells:
            state.view_well = wells[0]
        cmap_col, info_col = st.columns([3, 2])
        with cmap_col:
            st.caption("Click a well on the map, or choose it in the list.")
            ev = vega(charts.plate_map(r, plate, style, state.view_well), key=f"map_{vp}")
            w_click = picked(ev, "well", "well")
            if w_click and w_click in wells and w_click != state.get("last_click"):
                state.view_well, state.last_click = w_click, w_click
        with info_col:
            labels = {w: f"{w} — {LV(int(lv))}" for w, lv in zip(r["Well"], r["Confidence_Level"])}
            well = st.selectbox("Well", wells, index=wells.index(state.view_well), format_func=labels.get)
            state.view_well = well
            row = r.set_index("Well").loc[well]
            st.markdown(f"**{well} · {row.get('Sample_Code', '')} — level {labels[well].split(' — ')[1]}**")
            facts = {"Call": row["Call"],
                     "Evidence (cps)": f"{row['Evidence_Intensity']:,.0f}" if np.isfinite(row["Evidence_Intensity"]) else "–",
                     "Raw signal (cps)": f"{row['Raw_Ejection_Intensity']:,.0f}" if np.isfinite(row["Raw_Ejection_Intensity"]) else "–",
                     "Ions detected": row.get("Ions_Detected", "") or "none",
                     "Isotopes expected": row.get("Isotope_Theoretical", ""),
                     "Isotopes observed": row.get("Isotope_Observed", "") or "–",
                     "Pattern score (0-1)": row.get("Isotope_Pattern_Score", ""),
                     "Note": row.get("Level_Note", ""), "Reagent interference": row.get("SM_Interference", ""),
                     "Product formula": row.get("Formula", "")}
            for extra_col, lab in (("Control", "Control"), ("Replicate_Levels", "Replicate levels"),
                                   ("P_Product", "P(product)"), ("Run_QC_Flag", "Run QC")):
                if extra_col in row.index and str(row[extra_col]) not in ("", "nan"):
                    facts[lab] = row[extra_col]
            if "HPLC_Result" in row.index:
                facts["HPLC-MS"] = f"{row['HPLC_Result']}" + (f" ({row['HPLC_Conversion']:.1f}%)"
                                                             if np.isfinite(row.get("HPLC_Conversion", np.nan)) else "")
            st.dataframe(pd.DataFrame({"": list(facts), "Value": [str(v) for v in facts.values()]}), hide_index=True)

            st.markdown("**Your decision for this well**")
            cur = state.review[(state.review["Plate"] == str(vp)) & (state.review["Well"] == well)]
            opts = ["", "present", "absent", "hplc"]
            dec = st.radio("Decision", opts, index=opts.index(cur["Decision"].iloc[0]) if len(cur) else 0,
                           format_func=lambda o: "no decision" if not o else review_mod.DECISIONS[o],
                           key=f"dec_{vp}_{well}", horizontal=False)
            com = st.text_input("Comment", cur["Comment"].iloc[0] if len(cur) else "", key=f"com_{vp}_{well}")
            if st.button("Save decision", icon=":material/check:", key=f"save_{vp}_{well}"):
                state.review = review_mod.set_decision(state.review, str(vp), well, dec, com)
                st.success("Saved — the worklist in Results is updated.")
            if len(state.review):
                st.download_button(f"Download review.csv ({len(state.review)} decisions)",
                                   state.review.to_csv(index=False), "review.csv", "text/csv",
                                   icon=":material/download:")

        tr = spectra.well_trace(pa, well, s_)
        if not tr:
            st.warning("No ejection was assigned to this well.")
        else:
            st.markdown("**Ejection trace** — click a point to show that spectrum (the band is the window analysed)")
            key_scan = f"scan_{vp}_{rep}_{well}"
            if key_scan not in state:
                state[key_scan] = tr["best_scan"] if tr["best_scan"] is not None else int(tr["scan"][len(tr["scan"]) // 2])
            ev = vega(charts.trace(tr, style, state[key_scan]), key=f"trace_{vp}_{rep}_{well}")
            sc_click = picked(ev, "scan", "scan")
            if sc_click is not None:
                state[key_scan] = int(sc_click)
            scans = [int(x) for x in tr["scan"]]
            state[key_scan] = st.select_slider(
                "Spectrum (scan)", scans, value=state[key_scan] if state[key_scan] in scans else scans[0],
                format_func=lambda k: f"{pa.run.t[k]:.2f} s" + (" · in window" if k in tr["window_scans"] else ""))
            mode = st.radio("Show", ["one spectrum (as analysed)", "one spectrum (raw mzML)", "maximum of the window",
                                     "mean of the window"], horizontal=True)
            scan = state[key_scan]
            if mode.startswith("one spectrum (raw"):
                try:
                    mz, it = spectra.spectrum(pa, scan, "raw")
                except Exception as e:  # noqa: BLE001 - file moved or unreadable
                    st.error(f"Could not read the raw spectrum: {e}")
                    mz, it = spectra.spectrum(pa, scan)
            elif mode.startswith("one"):
                mz, it = spectra.spectrum(pa, scan)
            else:
                mz, it = spectra.mean_spectrum(pa, tr["window_scans"], "max" if mode.startswith("max") else "mean")
            marks = spectra.markers(pa, well, s_)
            lo, hi = spectra.view_range(marks)
            wide = st.checkbox("Whole scanned range", value=False)
            view = (float(pa.run.mz_min), float(pa.run.mz_max)) if wide else (lo, hi)
            m_obs = {}
            for m in marks:
                if m["kind"] == "ion":
                    near = np.abs(mz - m["mz"]) <= pa.run.tol(m["mz"], s_.mass_tolerance_da)
                    m_obs[m["ion"]] = float(it[near].max()) if near.any() else 0.0
            st.markdown("**Spectrum** — dashed: calculated m/z of each ion analysed; bars: theoretical isotope peaks; "
                        "dotted: reagent ions; horizontal: cutoff. Drag to zoom, double-click to reset.")
            vega(charts.spectrum(mz, it, marks, style, s_.cutoff_cps, view, m_obs), key=f"spec_{vp}_{rep}_{well}",
                 select=False)
            page = well_page(pa, well, s_, style, scan=scan,
                             mode="raw" if mode.startswith("one spectrum (raw") else "analysed")
            e1, e2, e3 = st.columns(3)
            e1.download_button("PNG", figure_bytes(page, "png"), f"Plate{vp}_{well}_spectrum.png", "image/png",
                               icon=":material/image:")
            e2.download_button("SVG", figure_bytes(page, "svg"), f"Plate{vp}_{well}_spectrum.svg", "image/svg+xml",
                               icon=":material/image:")
            e3.download_button("Spectrum (CSV)", pd.DataFrame({"mz": mz, "intensity": it}).to_csv(index=False),
                               f"Plate{vp}_{well}_scan{scan}.csv", "text/csv", icon=":material/table:")
            import matplotlib.pyplot as plt
            plt.close(page)

# --------------------------------------------------------------------------------------- statistics
with tab_stat:
    if res is None:
        st.info("Run the analysis first.")
    else:
        stt = res.stats
        lt = stt.get("levels")
        if lt is not None and len(lt):
            st.markdown("#### Share of wells with product at each level (95% confidence intervals)")
            st.dataframe(lt, hide_index=True)
        cal = stt.get("calibration") or {}
        if cal:
            st.markdown("#### Threshold calibration")
            st.caption(f"Lowest 'present' threshold and highest 'absent' threshold that keep the errors ≤ "
                       f"{cal['target_pct']:g}% of the wells; with several plates, fitted on the others and tested on each.")
            st.dataframe(cal["table"], hide_index=True)
            a, b, c = st.columns(3)
            a.markdown(ui.card("Suggested level 5", f"{cal['suggested_level5_cps']:,.0f}"
                               if np.isfinite(cal["suggested_level5_cps"]) else "–", "cps"), unsafe_allow_html=True)
            b.markdown(ui.card("Suggested 'absent' limit", f"{cal['suggested_absent_raw_cps']:,.0f}"
                               if np.isfinite(cal["suggested_absent_raw_cps"]) else "–", "cps"), unsafe_allow_html=True)
            c.markdown(ui.card("AUC evidence / raw", f"{cal['auc_evidence']} / {cal['auc_raw']}"), unsafe_allow_html=True)
        for t, f, prev in res.figures.get("all", []):
            if prev is not None and Path(prev).exists():
                st.markdown(f"#### {t}")
                st.image(str(prev))
        prob = stt.get("probability") or {}
        if prob.get("model") is not None:
            st.markdown("#### Probability of product per well")
            st.caption("Logistic regression on the Echo-MS measurements; 'plates except X → plate X' rows are the "
                       "honest check (the model never saw that plate). AUC 1 = perfect, 0.5 = chance; Brier: lower is "
                       "better. P_Product is in the results.")
            if len(prob["validation"]):
                st.dataframe(prob["validation"], hide_index=True)
            mpath = Path(res.out_dir) / "6_statistics" / "probability_model.json"
            if mpath.exists():
                st.download_button("Download the model (apply it to new plates)", mpath.read_bytes(),
                                   "probability_model.json", "application/json", icon=":material/download:")
        if len(res.control_charts):
            st.markdown("#### Control charts across plates")
            flagged = res.control_charts[res.control_charts["Flag"].fillna("") != ""]
            st.caption(f"{len(res.history)} plate run(s) in the project history; {len(flagged)} flagged point(s).")
            if len(flagged):
                st.dataframe(flagged[["Run", "Run_date", "mzML", "Plate", "Label", "Value", "Flag"]], hide_index=True)
        elif state.get("hist_path"):
            st.caption("Control charts appear when the project history has at least 6 plate runs.")
        st.markdown("#### Run QC")
        for qp, q in (stt.get("run_qc") or {}).items():
            if q is None or q.empty:
                continue
            flagged = q[q["Run_QC_Flag"] != ""]
            st.markdown(f"Plate {qp}: {len(flagged)} flagged ejection(s)")
            if len(flagged):
                st.dataframe(flagged, hide_index=True)

# --------------------------------------------------------------------------------------- help
with tab_help:
    st.markdown(f"""
### Files the program needs
| File | Required | What it holds |
|---|---|---|
| mzML | yes, one per plate | the Echo-MS run (msConvert, MS1, no MS-Numpress; high resolution: centroided) |
| Sample table | yes | the Excel template or a list: Plate, Well, Sample_Code, Product_SMILES, Reagent_1, Reagent_2, …, Injection_Order, Control, Replicate |
| HPLC-MS results | no | Plate, Well, Conversion_% — in the template's HPLC-MS sheet or separate files |
| Review decisions | no | review.csv saved from the spectrum viewer |
| Probability model | no | probability_model.json from a run with HPLC-MS |

**Reagents.** Write each reagent once in the Reagents sheet (Code, SMILES) and use the code elsewhere. *Reagent_1* should
be the reagent shared by a series of wells: the background control compares each well with the wells sharing it.
**Controls**: mark wells as positive, negative or blank in a Control column; the plate is accepted when the controls
behave. **Replicates**: number repeated ejections of a well in a Replicate column, with the Injection_Order of every
ejection.

### m/z range
A product whose ion falls outside the scanned range is **× Not measurable** and always sent to HPLC-MS. Below
~200 Da solvent ions are dense (more false positives); small molecules have a weak M+1 peak.

### Standard method
The detection and grading criteria are fixed (sidebar → *Detection criteria — standard method*): only the ions are
chosen, so analyses are comparable. Non-standard criteria can be unlocked to adapt the method to another instrument;
every output is then marked **NON-STANDARD PARAMETERS**.

### Levels
| Level | Meaning | Rule (defaults) |
|---|---|---|
| 5 | {LEVEL_LABELS[5]} | validated evidence ≥ 1,000,000 cps |
| 4 | {LEVEL_LABELS[4]} | validated evidence ≥ cutoff (200,000 cps) |
| 3 | {LEVEL_LABELS[3]} | weaker evidence, or a reagent ion on the product ion, carryover, only a secondary adduct |
| 2 | {LEVEL_LABELS[2]} | no validated ion, raw signal above level 1 |
| 1 | {LEVEL_LABELS[1]} | raw signal ≤ 100,000 cps |
| × | {LEVEL_LABELS[0]} | cannot be measured in this run: sent to HPLC-MS |

### Isotope pattern
Shown for every well (expected, observed, score) and drawn in the spectrum viewer; it does not change the level.

### Review, export and report
In the spectrum viewer, record a decision for any well (confirmed present, confirmed absent, send to HPLC-MS); the
worklist follows it. *Results* gives the worklist as an Excel transfer list and injection sequence (generic columns,
easy to map to any instrument), and the analysis report (HTML/PDF).

### Privacy
The analysis runs on this computer; uploads are copied to `{WORK}`. Profiles are kept in `{ui.profile_dir()}`.
""")

st.markdown(ui.footer_html(__version__), unsafe_allow_html=True)
