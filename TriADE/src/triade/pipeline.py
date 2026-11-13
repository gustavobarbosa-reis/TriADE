"""The whole analysis in one call: read the runs, assign ejections to wells, search the product ions,
grade each well, compare with HPLC-MS (optional) and write every table and figure.

    from triade import Inputs, Settings, run_analysis
    res = run_analysis(Inputs(mzml=["plate1.mzML"], samples="samples.csv"), Settings(), "results")
"""
from __future__ import annotations

import logging
import math
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import __version__
from . import figures as fig
from . import history as history_mod
from . import replicates, sideproducts, stats
from . import review as review_mod
from .adducts import PRIMARY, ions_for
from .detection import PlateAnalysis, analyse_plate, first_slot_for_table
from .hplc import add_hplc, load_hplc, mark_indeterminate, parse_exclusions
from .plates import get_plate_format
from .runs import plate_id_from_name, read_run
from .samples import add_masses, apply_row_starting_materials, load_samples
from .settings import LEVEL_LABELS, Settings
from .templates import template_hplc
from .triage import SCENARIOS, discrepancies, hplc_worklist, triage_table

log = logging.getLogger("triade")


@dataclass
class Inputs:
    """Input files. Plate IDs link each mzML (and HPLC-MS file) to the 'Plate' of the sample table;
    when not given they are read from the file name ('... Plate 2 ...') or numbered 1, 2, ..."""
    mzml: list
    samples: str | Path
    mzml_plates: list | None = None
    starting_materials_by_row: str | Path | None = None   # optional table Row | SMILES (| Plate)
    hplc: list = field(default_factory=list)
    hplc_plates: list | None = None
    hplc_exclude: list = field(default_factory=list)       # e.g. ["2:9"] = plate 2, column 9 indeterminate
    probability_model: str | Path | None = None            # a saved model (probability_model.json) to apply
    review: str | Path | None = None                       # manual review decisions (review.csv)
    history: str | Path | None = None                      # project history CSV (control charts across plates)


@dataclass
class Results:
    out_dir: Path
    settings: Settings
    plates: dict                 # plate id -> PlateAnalysis
    results: dict                # plate id -> result frame
    summary: pd.DataFrame
    triage: pd.DataFrame
    worklist: pd.DataFrame
    discrepancies: pd.DataFrame
    warnings: list
    files: list
    zip_path: Path | None = None
    figures: dict = field(default_factory=dict)       # plate id (or "all") -> [(title, file, preview png)]
    stats: dict = field(default_factory=dict)
    review: pd.DataFrame = field(default_factory=lambda: review_mod.empty())
    history: pd.DataFrame = field(default_factory=pd.DataFrame)
    control_charts: pd.DataFrame = field(default_factory=pd.DataFrame)
    axis_labels: dict = field(default_factory=dict)   # plate id -> (row labels, column labels) for the maps

    @property
    def has_hplc(self) -> bool:
        return all("HPLC_Product" in r.columns for r in self.results.values()) and bool(self.results)

    @property
    def all_results(self) -> pd.DataFrame:
        return pd.concat(self.results.values(), ignore_index=True) if self.results else pd.DataFrame()


def _criteria_text(s: Settings, acquired: set) -> str:
    sel = [i for i in ions_for(s.ions, s.polarity) if i.polarity in acquired]
    ions = ", ".join(i.name for i in sel) or "none"
    txt = (f"Ions {ions}  ·  ±{s.mass_tolerance_da:g} Da  ·  M→M+1 shift 1/|z| ± {s.isotope_shift_tolerance_da:g} Da  ·  "
           f"same ion in ≥ {s.min_spectra} spectra")
    if s.background_group != "off" and s.background_ratio > 0:
        txt += f"  ·  ≥ {s.background_ratio:g}x {s.background_group.replace('_', ' ')} background"
    if not s.is_standard:
        txt += "  ·  NON-STANDARD PARAMETERS"
    return txt


def _plate_ids(paths, given, fallback_prefix=""):
    ids = []
    for k, p in enumerate(paths):
        g = given[k] if given and k < len(given) and given[k] not in (None, "") else None
        ids.append(str(g) if g is not None else (plate_id_from_name(Path(p).name) or f"{fallback_prefix}{k + 1}"))
    return ids


def _axis_labels(df: pd.DataFrame, plate) -> tuple:
    """Reagent codes for the map axes: a row label when every well of the row has the same Reagent_1 code,
    a column label when every well of the column has the same Reagent_2 code."""
    rows, cols = {}, {}
    rc = next((c for c in ("Reagent_1_Code",) if c in df.columns), None)
    cc = next((c for c in ("Reagent_2_Code",) if c in df.columns), None)
    pos = df["Well"].map(plate.split)
    if rc:
        for r, g in df.groupby(pos.map(lambda x: x[0])):
            v = set(g[rc]) - {""}
            if len(v) == 1:
                rows[r] = v.pop()
    if cc:
        for c, g in df.groupby(pos.map(lambda x: x[1])):
            v = set(g[cc]) - {""}
            if len(v) == 1:
                cols[int(c)] = v.pop()
    return (rows or None, cols or None)


def _report(res, out: Path) -> list:
    """Report figures are drawn as PNG whatever the chosen format."""
    from dataclasses import replace
    from . import report
    rd = out.parent / f".{out.name}_report"
    shutil.rmtree(rd, ignore_errors=True)
    rd.mkdir()
    st = fig.Style.from_settings(res.settings)
    st = replace(st, fmt="png", dpi=110, width_mm=None, titles=True)
    plate = get_plate_format(res.settings.plate_format)
    figs = {}
    try:
        for pid, r in res.results.items():
            labs = res.axis_labels.get(pid) if res.settings.axis_reagent_labels else None
            items = [("Confidence levels", fig.heatmap_levels(r, plate, rd / f"lv{pid}", f"Plate {pid} — confidence levels",
                                                              "", st, labels=labs))]
            if "HPLC_Conversion" in r.columns:
                items.append(("Echo-MS levels and HPLC-MS conversion", fig.composite_levels_conversion(
                    r, plate, rd / f"cmp{pid}", f"Plate {pid}", res.settings.reference_label, st, labels=labs)))
                items.append(("Levels vs. HPLC-MS", fig.levels_vs_reference_chart(r, rd / f"lvr{pid}", pid,
                                                                                   res.settings, st)))
            figs[pid] = items
        return report.write(res, figs, out)
    finally:
        fig.Style.from_settings(res.settings)        # restore the chosen style for later figures
        shutil.rmtree(rd, ignore_errors=True)


def _write_csv(df, path, files):
    df.to_csv(path, index=False, encoding="utf-8-sig")
    files.append(path)


def _preview(path: Path, preview_dir: Path) -> Path | None:
    """A PNG copy for the interface when the chosen format cannot be shown in a browser."""
    if path.suffix.lower() in (".jpg", ".png"):
        return path
    try:
        from PIL import Image
        if path.suffix.lower() == ".tiff":
            preview_dir.mkdir(exist_ok=True)
            p = preview_dir / (path.stem + ".png")
            Image.open(path).convert("RGB").save(p)
            return p
    except Exception:  # pragma: no cover
        return None
    return None


def run_analysis(inputs: Inputs, settings: Settings | None = None, out_dir="triade_results",
                 figures: bool = True, make_zip: bool = True, progress=None) -> Results:
    """Runs everything and writes the results to out_dir. progress(fraction, message) is called on the way."""
    s = (settings or Settings()).validate()
    plate = get_plate_format(s.plate_format)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    say = progress or (lambda f, m: None)
    files: list = []
    warnings: list = []
    figs: dict = {}
    stat_out: dict = {}

    handler = logging.FileHandler(out / "log.txt", mode="w", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(levelname)s  %(message)s"))
    log.addHandler(handler)
    t_start = time.time()
    try:
        log.info("TriADE %s", __version__)
        if not s.is_standard:
            msg = (f"{s.method_label()}: the detection criteria differ from the standard method, so these results "
                   "are not directly comparable with standard analyses. Validate them against HPLC-MS.")
            warnings.append(msg); log.warning(msg)
        else:
            log.info("Detection criteria: standard method.")
        sel = ions_for(s.ions, s.polarity)
        log.info("Instrument: %s; plate format %s wells; polarity %s; ions %s",
                 "high resolution (experimental)" if s.instrument == "high" else "unit resolution",
                 plate.name, s.polarity, ", ".join(i.name for i in sel))
        if len(sel) > 1:
            log.info("More than one ion per product: every extra adduct raises the chance of false positives.")

        # ---------- samples ----------
        say(0.02, "Reading the sample table")
        samples = load_samples(inputs.samples, plate)
        if inputs.starting_materials_by_row:
            samples = apply_row_starting_materials(samples, inputs.starting_materials_by_row, plate)
        samples = add_masses(samples, s.ions, s.polarity)

        mz_ids = _plate_ids(inputs.mzml, inputs.mzml_plates)
        if len(set(mz_ids)) != len(mz_ids):
            raise ValueError(f"two mzML files have the same plate ID ({', '.join(mz_ids)}); set the plate of each file")
        known = sorted(samples["Plate"].unique())
        if len(mz_ids) == 1 and mz_ids[0] not in known and len(known) == 1:
            log.info("One mzML file and one plate in the sample table: plate %s used for %s", known[0],
                     Path(inputs.mzml[0]).name)
            mz_ids = [known[0]]
        missing = [p for p in mz_ids if p not in known]
        if missing:
            raise ValueError(f"plate(s) {', '.join(missing)} of the mzML files are not in the sample table "
                             f"(plates there: {', '.join(known)}). Set the plate of each mzML file.")

        # ---------- runs and detection ----------
        plates: dict[str, PlateAnalysis] = {}
        metrics = {pid: {} for pid in mz_ids}
        axis_labels: dict = {}
        rules = sideproducts.parse_rules(s.side_products)
        echo_minutes = 0.0
        n = len(inputs.mzml)
        for k, (path, pid) in enumerate(zip(inputs.mzml, mz_ids)):
            say(0.05 + 0.55 * k / n, f"Plate {pid}: reading {Path(path).name}")
            run = read_run(path, s.mz_bin_da, s)
            warnings += [f"Plate {pid}: {w}" for w in getattr(run, "warnings", [])]
            echo_minutes += run.duration_min
            df = samples[samples["Plate"] == pid].reset_index(drop=True)
            log.info("Plate %s: %d wells", pid, len(df))
            if not (set(run.polarities) & {i.polarity for i in sel}):
                msg = (f"Plate {pid}: the run has only {'positive' if '+' in run.polarities else 'negative'} scans "
                       "but none of the selected ions has that polarity: nothing can be detected.")
                warnings.append(msg); log.warning(msg)
            say(0.05 + 0.55 * (k + 0.4) / n, f"Plate {pid}: searching the product ions")
            reps = sorted(df["Replicate"].unique()) if "Replicate" in df.columns else [1]
            if len(reps) > 1:
                df = df.assign(_Order_Base=df["Injection_Order"].min())
                first = first_slot_for_table(run, df, s)
                pa = analyse_plate(run, df[df["Replicate"] == reps[0]].reset_index(drop=True), plate, pid, s,
                                   first_slot=first)
                pa.replicates = {reps[0]: pa}
                for r in reps[1:]:
                    pa.replicates[r] = analyse_plate(run, df[df["Replicate"] == r].reset_index(drop=True), plate,
                                                     pid, s, first_slot=pa.first_slot)
                frames = {r: x.results for r, x in pa.replicates.items()}
                pa.results = replicates.combine(frames, s)
                metrics[pid].update(replicates.replicate_summary(frames, pa.results))
                log.info("  %d replicate ejections per well combined (%s rule); Fleiss kappa %s", len(reps),
                         s.replicate_rule, metrics[pid].get("Replicate_kappa_Fleiss"))
            else:
                pa = analyse_plate(run, df, plate, pid, s)
            if rules:
                sp = sideproducts.search(pa, df.drop_duplicates("Well"), rules, s)
                pa.results = pa.results.merge(sp, on="Well", how="left")
                for rl in rules:
                    metrics[pid][f"SideProduct_{rl['key']}_wells"] = int(pa.results[f"SP_{rl['key']}_Detected"].sum())
            axis_labels[pid] = _axis_labels(df, plate)
            warnings += [f"Plate {pid}: {w}" for w in pa.warnings]
            n_nm = int((pa.results["Confidence_Level"] == 0).sum())
            if n_nm:
                msg = (f"Plate {pid}: {n_nm} well(s) not measurable (product ion outside the scanned range "
                       f"m/z {run.mz_min:.0f}-{run.mz_max:.0f}, no structure, polarity not acquired or no ejection): "
                       "they are always sent to HPLC-MS.")
                warnings.append(msg); log.warning(msg)
            plates[pid] = pa
        results = {pid: pa.results for pid, pa in plates.items()}

        # ---------- review and controls ----------
        review_df = review_mod.load(inputs.review) if inputs.review else review_mod.empty()
        review_mod.apply(results, review_df)

        # ---------- HPLC-MS (optional) ----------
        ref_names = {}
        hplc_paths, hplc_ids = list(inputs.hplc), inputs.hplc_plates
        if not hplc_paths:
            sheet = template_hplc(inputs.samples)
            if sheet is not None:
                tmp = out / "_hplc_from_template.csv"
                sheet.to_csv(tmp, index=False)
                hplc_paths, hplc_ids = [tmp], None
                log.info("HPLC-MS results taken from the HPLC-MS sheet of %s", Path(inputs.samples).name)
        if hplc_paths:
            say(0.62, "Comparing with HPLC-MS")
            hp_ids = _plate_ids(hplc_paths, hplc_ids)
            if len(hplc_paths) == 1 and len(plates) == 1 and hp_ids[0] not in plates:
                hp_ids = [next(iter(plates))]
            hplc = load_hplc(hplc_paths, plate, s.hplc_product_min_pct, hp_ids)
            hplc = mark_indeterminate(hplc, parse_exclusions(inputs.hplc_exclude, plate))
            for path, pid in zip(hplc_paths, hp_ids):
                ref_names.setdefault(str(pid), []).append(Path(path).name)
            for pid, res in results.items():
                if hplc.empty or pid not in set(hplc["Plate"]):
                    msg = f"Plate {pid}: no {s.reference_label} results for this plate."
                    warnings.append(msg); log.warning(msg)
                    continue
                metrics[pid].update(add_hplc(res, hplc, pid, s))
                y, p = res["HPLC_Product"], res["Detected"]
                lo, hi = stats.bootstrap_kappa(y, p)
                metrics[pid]["Cohen_kappa_CI95"] = f"{lo}-{hi}"
                log.info("Plate %s vs %s: %s%% agreement, kappa %s (95%% CI %s-%s)", pid, s.reference_label,
                         metrics[pid]["Agreement_%"], metrics[pid]["Cohen_kappa"], lo, hi)
            (out / "_hplc_from_template.csv").unlink(missing_ok=True)

        for pid, res in results.items():
            cc = replicates.control_check(res, s)
            if cc:
                metrics[pid].update(cc)
                if cc["Plate_accepted"] == "no":
                    msg = f"Plate {pid} NOT accepted by its controls: {cc['Acceptance_note']}."
                    warnings.append(msg); log.warning(msg)

        # ---------- statistics ----------
        say(0.68, "Statistics")
        model = stats.ProbabilityModel.load(inputs.probability_model) if inputs.probability_model else None
        prob = stats.probability(results, model)
        stat_out["probability"] = prob
        stat_out["levels"] = stats.level_table(results)
        stat_out["calibration"] = stats.calibrate(results, s)
        stat_out["run_qc"] = {pid: stats.run_qc(pa, s) for pid, pa in plates.items()}
        for pid, q in stat_out["run_qc"].items():
            if not q.empty:
                results[pid]["Run_QC_Flag"] = results[pid]["Well"].map(dict(zip(q["Well"], q["Run_QC_Flag"]))).fillna("")
                d = stats.drift(q)
                if d:
                    metrics[pid]["Pulse_drift_spearman"] = d["spearman_vs_time"]
                    metrics[pid]["Pulse_last_vs_first_quarter_%"] = d["last_quarter_vs_first_%"]
                if "IS_Intensity" in q and q["IS_Intensity"].notna().any():
                    di = stats.drift(q, "IS_Intensity")
                    metrics[pid]["IS_CV_%"] = round(100 * q["IS_Intensity"].std() / q["IS_Intensity"].mean(), 1)
                    metrics[pid]["IS_drift_spearman"] = di.get("spearman_vs_time")
                nflag = int((q["Run_QC_Flag"] != "").sum())
                if nflag:
                    msg = f"Plate {pid}: {nflag} ejection(s) flagged by the run QC (weak pulse or internal standard)."
                    warnings.append(msg); log.warning(msg)

        # ---------- summary, triage, worklist ----------
        say(0.72, "Triage and tables")
        summary_rows = []
        for pid, res in results.items():
            pa = plates[pid]
            cnt = res["Confidence_Level"].value_counts()
            row = {"Plate": pid, "mzML": pa.run.path.name, "Method": "standard" if s.is_standard else "non-standard",
                   "Wells": len(res),
                   "Run_min": round(pa.run.duration_min, 2), "Scan_range_mz": f"{pa.run.mz_min:.1f}-{pa.run.mz_max:.1f}",
                   "Ejection_period_s": round(pa.period, 4),
                   "Ejection_order": pa.order_name, "First_ejection_s": round(float(pa.grid[pa.first_slot]), 2),
                   "Polarities_acquired": "".join(sorted(pa.run.polarities)),
                   "Ions": ", ".join(i.name for i in sel), "Cutoff_cps": s.cutoff_cps,
                   "Detected": int(res["Detected"].sum()),
                   "Not_processed": int((res["Classification"] == "NOT PROCESSED").sum())}
            for lv in (5, 4, 3, 2, 1, 0):
                row[f"Level_{lv}_{LEVEL_LABELS[lv].replace(' ', '_')}"] = int(cnt.get(lv, 0))
            row.update(metrics.get(pid, {}))
            summary_rows.append(row)
        summary = pd.DataFrame(summary_rows)
        tri = stats.triage_ci(triage_table(results, echo_minutes, s))
        work = hplc_worklist(results, s.triage_scenario)
        disc = discrepancies(results)

        # ---------- files ----------
        say(0.76, "Writing tables")
        d_det, d_lv, d_tri = out / "1_detection", out / "2_confidence_levels", out / "3_triage"
        d_time, d_qc, d_st = out / "4_time_charts", out / "5_QC", out / "6_statistics"
        for d in (d_det, d_lv, d_tri, d_time, d_qc, d_st):
            d.mkdir(exist_ok=True)
        lv_cols = ["Plate", "Well", "Sample_Code", "Product_SMILES", "Formula", "Monoisotopic_Mass", "Confidence_Level",
                   "Confidence_Label", "Call", "Level_Note", "P_Product", "Evidence_Intensity", "Raw_Ejection_Intensity",
                   "Previous_Ejection_Intensity", "Ions_Detected", "SM_Interference", "Echo_Product_Ratio_%",
                   "Isotope_Ion", "Isotope_Theoretical", "Isotope_Observed", "Isotope_Pattern_Score",
                   "Run_QC_Flag", "HPLC_Conversion", "HPLC_Result", "HPLC_Product", "Agreement"]
        for pid, res in results.items():
            _write_csv(res, d_det / f"Plate{pid}_detection.csv", files)
            _write_csv(res[[c for c in lv_cols if c in res.columns]], d_lv / f"Plate{pid}_confidence_levels.csv", files)
            _write_csv(plates[pid].order_scores, d_qc / f"Plate{pid}_ejection_order_hypotheses.csv", files)
            q = stat_out["run_qc"].get(pid)
            if q is not None and not q.empty:
                _write_csv(q, d_qc / f"Plate{pid}_run_QC.csv", files)
        _write_csv(summary, out / "Summary.csv", files)
        _write_csv(tri, d_tri / "HPLC-MS_triage_savings.csv", files)
        _write_csv(work, d_tri / f"HPLC-MS_worklist_scenario_{s.triage_scenario}.csv", files)
        files.append(review_mod.write_excel(d_tri / f"HPLC-MS_transfer_and_sequence_scenario_{s.triage_scenario}.xlsx",
                                            work, s, SCENARIOS[s.triage_scenario][0]))
        if len(review_df):
            _write_csv(review_df, d_tri / "Review_decisions.csv", files)
        if not disc.empty:
            _write_csv(disc, d_tri / f"Discrepancies_EchoMS_vs_{s.reference_label}.csv", files)
        if not stat_out["levels"].empty:
            _write_csv(stat_out["levels"], d_st / "Levels_with_product_CI95.csv", files)
        cal = stat_out["calibration"]
        if cal:
            _write_csv(cal["table"], d_st / "Threshold_calibration.csv", files)
            _write_csv(cal["present_curve"], d_st / "Present_threshold_curve.csv", files)
            _write_csv(cal["absent_curve"], d_st / "Absent_threshold_curve.csv", files)
        if not prob["validation"].empty:
            _write_csv(prob["validation"], d_st / "Probability_model_validation.csv", files)
        if prob.get("model") is not None and inputs.probability_model is None:
            files.append(prob["model"].save(d_st / "probability_model.json"))
        s.save(out / "settings.json"); files.append(out / "settings.json")
        if warnings:
            (out / "warnings.txt").write_text("\n".join(warnings) + "\n", encoding="utf-8")
            files.append(out / "warnings.txt")
        xlsx = out / "All_results.xlsx"
        with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
            summary.to_excel(xw, sheet_name="Summary", index=False)
            for pid, res in results.items():
                res.to_excel(xw, sheet_name=f"Plate {pid}"[:31], index=False)
            tri.to_excel(xw, sheet_name="Triage", index=False)
            work.to_excel(xw, sheet_name="HPLC-MS worklist", index=False)
            if not disc.empty:
                disc.to_excel(xw, sheet_name="Discrepancies", index=False)
            if not stat_out["levels"].empty:
                stat_out["levels"].to_excel(xw, sheet_name="Levels CI95", index=False)
            if cal:
                cal["table"].to_excel(xw, sheet_name="Threshold calibration", index=False)
            pd.DataFrame(list(s.to_dict().items()), columns=["Setting", "Value"]).astype(str).to_excel(
                xw, sheet_name="Settings", index=False)
        files.append(xlsx)

        # ---------- project history and control charts ----------
        hist, charts_df = pd.DataFrame(), pd.DataFrame()
        tmp_res = Results(out, s, plates, results, summary, tri, work, disc, warnings, files, None, {}, stat_out)
        if inputs.history:
            hist = history_mod.update(inputs.history, history_mod.rows_for_run(tmp_res))
            charts_df = history_mod.control_limits(hist)
            _write_csv(hist, d_qc / "Project_history.csv", files)
            if not charts_df.empty:
                _write_csv(charts_df, d_qc / "Control_chart_limits.csv", files)
                cur = charts_df[charts_df["mzML"].isin(summary["mzML"]) & (charts_df["Flag"] != "")]
                for r in cur.itertuples(index=False):
                    msg = f"Plate {r.Plate}: control chart — {r.Label} = {r.Value:.3g} ({r.Flag})."
                    warnings.append(msg); log.warning(msg)

        # ---------- figures ----------
        if figures:
            style = fig.Style.from_settings(s)
            fig.CAPTIONS.clear()
            prev_dir = out.parent / f".{out.name}_preview"
            shutil.rmtree(prev_dir, ignore_errors=True)

            def add(key, title, path):
                if path is not None:
                    figs.setdefault(key, []).append((title, Path(path), _preview(Path(path), prev_dir)))

            allv = pd.concat([r["Log10_Max_Validated_Intensity"] for r in results.values()]).dropna()
            vmin = math.log10(s.cutoff_cps)
            vmax = max(math.ceil(allv.max() * 10) / 10 if len(allv) else vmin + 1, vmin + 0.5)
            ref = s.reference_label
            for k, (pid, res) in enumerate(results.items()):
                say(0.78 + 0.14 * k / len(results), f"Plate {pid}: figures")
                pa = plates[pid]
                labs = axis_labels.get(pid) if s.axis_reagent_labels else None
                crit = _criteria_text(s, pa.run.polarities)
                cut = f"cutoff {s.cutoff_cps:,.0f} cps"
                n_proc = int((res["Confidence_Level"] > 0).sum())
                nd = int(res["Detected"].sum())
                stem = d_det / f"Plate{pid}"
                add(pid, "Confidence levels", fig.heatmap_levels(
                    res, plate, d_lv / f"Plate{pid}_confidence_levels", f"Plate {pid} — confidence levels (Echo-MS)",
                    crit + f"\n5: evidence ≥ {s.level5_cps:,.0f} cps · 4: ≥ {s.cutoff_cps:,.0f} · "
                           f"3: ≥ {s.level3_cps:,.0f} or ambiguous · 2: raw signal > {s.level1_raw_max_cps:,.0f} · "
                           f"1: raw ≤ {s.level1_raw_max_cps:,.0f} cps · ×: not measurable (sent to {ref})", style, labels=labs))
                add(pid, "Detected / check / not detected", fig.heatmap_call(
                    res, plate, f"{stem}_3_categories", f"Plate {pid} — detected / check spectrum / not detected ({cut})",
                    crit + "\nCheck spectrum: confidence level 3 (weak or ambiguous evidence, reagent ion on the "
                           "product ion, possible carryover, or only a secondary adduct passed).", style, labels=labs))
                add(pid, "log10 intensity", fig.heatmap_log10(
                    res, plate, f"{stem}_log10_intensity",
                    f"Plate {pid} — log10 intensity of validated product ions ({cut})",
                    crit + f"\nHighest validated M-peak intensity per well (≥ 1 spectrum). "
                           f"{int(res['Log10_Max_Validated_Intensity'].notna().sum())}/{len(res)} wells with signal; "
                           f"{nd} meet the ≥{s.min_spectra}-spectra rule.", vmin, vmax, style, labels=labs))
                add(pid, "Binary detection", fig.heatmap_binary(
                    res, plate, f"{stem}_binary",
                    f"Plate {pid} — product detection ({cut}): {nd}/{n_proc} detected "
                    f"({100 * nd / n_proc if n_proc else 0:.1f}%)",
                    crit + (f"\nM peak ≥ {s.cutoff_cps:,.0f} cps; M+1 peak ≥ {s.cutoff_cps * s.m1_fraction_of_cutoff:,.0f} cps."
                            if s.m1_rule == "absolute" else
                            f"\nM peak ≥ {s.cutoff_cps:,.0f} cps; M+1 ≥ {s.m1_relative_fraction:g} × the M+1 expected "
                            f"from the formula and ≥ {s.m1_floor_cps:,.0f} cps."), style, labels=labs))
                if "Agreement" in res.columns:
                    m = metrics[pid]
                    add(pid, f"Echo-MS vs. {ref}", fig.heatmap_agreement(
                        res, plate, f"{stem}_vs_{ref}",
                        f"Plate {pid} — Echo-MS vs. {ref} ({cut}): {m['Agreement_%']}% agreement, "
                        f"Cohen's kappa {m['Cohen_kappa']} (95% CI {m.get('Cohen_kappa_CI95', '')})",
                        crit + f"\n{ref}: {', '.join(ref_names.get(pid, []))} ({m['Reference_with_product']} of "
                               f"{m['Reference_wells']} wells with product, conversion ≥ {s.hplc_product_min_pct:g}%; "
                               f"wells without an {ref} result counted as not detected).", ref, style, labels=labs))
                    add(pid, f"{ref} conversion", fig.heatmap_conversion(
                        res, plate, d_lv / f"Plate{pid}_{ref}_conversion", f"Plate {pid} — {ref} conversion (%)",
                        f"{', '.join(ref_names.get(pid, []))}; product counted from {s.hplc_product_min_pct:g}% "
                        "conversion", ref, style, labels=labs))
                    lt = stat_out["levels"]
                    sub = lt[lt["Plate"] == pid] if len(lt) else lt
                    ci = (dict(zip(sub["Level"].astype(int), zip(sub["CI95_low_%"], sub["CI95_high_%"])))
                          if len(sub) else None)
                    add(pid, "Echo-MS levels and HPLC-MS conversion", fig.composite_levels_conversion(
                        res, plate, d_lv / f"Plate{pid}_levels_and_{ref}_conversion",
                        f"Plate {pid} — Echo-MS confidence levels and {ref} conversion", ref, style, labels=labs))
                    add(pid, f"Levels vs. {ref}", fig.levels_vs_reference_chart(
                        res, d_lv / f"Plate{pid}_levels_vs_{ref}", pid, s, style, ci))
                add(pid, "Run time chart", fig.time_chart(pa, plate, d_time / f"Plate{pid}_time_chart_rows", style))
                add(pid, "QC: ejections on the TIC", fig.qc_tic(pa, d_qc / f"Plate{pid}_QC_TIC_ejections", s, style))
                add(pid, "QC: mean ejection profile", fig.qc_profile(pa, d_qc / f"Plate{pid}_QC_mean_ejection_profile",
                                                                     s, style))
                q = stat_out["run_qc"].get(pid)
                if q is not None and not q.empty:
                    add(pid, "QC: pulse height over the run", fig.run_qc_chart(q, pid, d_qc / f"Plate{pid}_run_QC",
                                                                               style))
            if cal:
                add("all", "Threshold calibration", fig.calibration_chart(cal, d_st / "Threshold_calibration", s, style))
            if rules:
                d_sp2 = out / "8_side_products"
                d_sp2.mkdir(exist_ok=True)
                for pid, res in results.items():
                    for rl in rules:
                        add(pid, f"Side product: {rl['name']}", fig.heatmap_side_product(
                            res, plate, d_sp2 / f"Plate{pid}_{rl['key']}", rl, style,
                            labels=axis_labels.get(pid) if s.axis_reagent_labels else None))
            if not charts_df.empty:
                add("all", "Control charts across plates", fig.control_chart_figure(charts_df, d_qc / "Control_charts",
                                                                                    style))
            if s.spectra_pdf != "none":
                say(0.93, "Spectra PDF")
                levels = {"check": [3], "levels": list(s.spectra_pdf_levels), "all": None}[s.spectra_pdf]
                d_sp = out / "7_spectra"
                d_sp.mkdir(exist_ok=True)
                pdf = fig.spectra_pdf(plates, d_sp / "Spectra_by_well.pdf", s, style, levels)
                if pdf:
                    files.append(pdf)
            files += sorted(p for p in out.rglob("*") if p.suffix in (".jpg", ".png", ".tiff", ".svg", ".pdf"))

        if figures and fig.CAPTIONS:
            (out / "Figure_captions.txt").write_text(
                "\n\n".join(f"{Path(stem).name}\n{title}\n{sub}".strip() for stem, title, sub in fig.CAPTIONS),
                encoding="utf-8")
            files.append(out / "Figure_captions.txt")
        result = Results(out, s, plates, results, summary, tri, work, disc, warnings, files, None, figs, stat_out,
                         review_df, hist, charts_df, axis_labels)
        if figures and s.report:
            say(0.96, "Report")
            files += _report(result, out)

        zip_path = None
        if make_zip:
            say(0.99, "Packing the results")
            log.info("Done in %.0f s.", time.time() - t_start)
            handler.flush()
            zip_path = Path(shutil.make_archive(str(out.parent / out.name), "zip", root_dir=out))
        say(1.0, "Done")
        result.zip_path = zip_path
        return result
    finally:
        log.removeHandler(handler)
        handler.close()


def scenario_names() -> list:
    return [name for name, *_ in SCENARIOS]


def primary_ion_names() -> list:
    return list(PRIMARY.values())


__all__ = ["Inputs", "Results", "run_analysis", "scenario_names", "primary_ion_names"]
