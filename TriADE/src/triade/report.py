"""Report of an analysis: Report.html (self-contained, opens in any browser) and Report.pdf, with a methods
paragraph written from the settings, the plate summaries, controls, triage, warnings and the main figures."""
from __future__ import annotations

import base64
import html
import io
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import __version__
from .adducts import ions_for
from .triage import SCENARIOS


def methods_text(res) -> str:
    s = res.settings
    pa = next(iter(res.plates.values()))
    ions = ", ".join(i.name for i in ions_for(s.ions, s.polarity))
    tol = f"±{s.mass_tolerance_ppm:g} ppm" if s.instrument == "high" else f"±{s.mass_tolerance_da:g} Da"
    m1 = (f"a distinct M+1 peak ≥ {s.m1_fraction_of_cutoff * 100:g}% of the cutoff" if s.m1_rule == "absolute" else
          f"an M+1 peak ≥ {s.m1_relative_fraction:g} × the M+1 expected from the formula (and ≥ {s.m1_floor_cps:,.0f} cps)")
    bg = ("" if s.background_group == "off" else
          f", and the ejection signal was ≥ {s.background_ratio:g}× the median of the same m/z in the wells sharing the "
          f"same {'reagent' if s.background_group == 'starting_material' else s.background_group}")
    txt = (f"Acoustic-ejection MS runs ({len(res.plates)} plate(s), {s.plate_format}-well format, "
           f"{'high-resolution' if s.instrument == 'high' else 'unit-resolution'} data, {s.polarity} mode) were analysed "
           f"with TriADE {__version__}. The interval between ejections was fitted to the total ion current "
           f"({pa.period:.3f} s) and each ejection was assigned to its well ({pa.order_name}). For each well the "
           f"product ion(s) {ions} were searched within {tol} of the calculated m/z in every spectrum from "
           f"{s.window_pre_s:g} s before to {s.window_post_s:g} s after the ejection apex. An ion was considered "
           f"validated when the M peak was ≥ {s.cutoff_cps:,.0f} cps with {m1} in at least {s.min_spectra} spectra{bg}. "
           f"Wells were graded in five confidence levels: 5 Present (evidence ≥ {s.level5_cps:,.0f} cps), 4 Likely "
           f"present (≥ {s.cutoff_cps:,.0f} cps), 3 Check spectrum (≥ {s.level3_cps:,.0f} cps, or a reagent ion on the "
           f"product ion, carryover ≥ {s.carryover_ratio:g}× from the previous ejection, or only a secondary adduct), "
           f"2 Likely absent (raw signal > {s.level1_raw_max_cps:,.0f} cps) and 1 Not detected; wells that could not be "
           f"measured were sent to HPLC-MS.")
    txt += (" The standard detection criteria of the program were used." if s.is_standard else
            f" Non-standard detection criteria were used: {s.method_label()[len('NON-STANDARD PARAMETERS ('):-1]}.")
    if res.has_hplc:
        txt += (f" HPLC-MS conversion ≥ {s.hplc_product_min_pct:g}% defined the presence of the product; wells without "
                "an HPLC-MS result were counted as not detected.")
    if any("Replicates" in r.columns for r in res.results.values()):
        txt += f" Replicate ejections were combined with the '{s.replicate_rule}' rule."
    return txt


def _png(path) -> str:
    return base64.b64encode(Path(path).read_bytes()).decode()


def _table_html(df: pd.DataFrame, max_rows: int = 60) -> str:
    return df.head(max_rows).to_html(index=False, border=0, classes="t", na_rep="–", float_format=lambda v: f"{v:,.4g}")


def write(res, figs: dict, out_dir: Path) -> list[Path]:
    """figs: plate id -> [(title, png path)] rendered for the report."""
    s = res.settings
    out_dir = Path(out_dir)
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    keep = [c for c in ["Plate", "mzML", "Wells", "Scan_range_mz", "Detected", "Level_5_Present", "Level_4_Likely_present",
                        "Level_3_Check_spectrum", "Level_2_Likely_absent", "Level_1_Not_detected", "Level_0_Not_measurable",
                        "Agreement_%", "Cohen_kappa", "Cohen_kappa_CI95", "Plate_accepted", "Acceptance_note",
                        "Replicate_kappa_Fleiss"] if c in res.summary.columns]
    summ = res.summary[keep]
    tri_cols = [c for c in ["Scenario", "HPLC_injections", "HPLC_injections_%_of_wells", "HPLC_time_h", "Time_saved_h",
                            "Net_time_saved_h_minus_EchoMS", "Mobile_phase_saved_mL", "False_present", "Products_missed"]
                if c in res.triage.columns]
    tri = res.triage[tri_cols]
    scen = SCENARIOS[s.triage_scenario][0]
    n_rev = sum(int((r.get("Review_Decision", pd.Series(dtype=str)).fillna("") != "").sum()) for r in res.results.values())

    # ---- HTML
    parts = [f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>TriADE report</title>
<style>body{{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;color:#2b2b29;max-width:1100px;margin:32px auto;
padding:0 20px;line-height:1.5}}h1{{font-size:26px;margin-bottom:0}}h2{{font-size:19px;margin-top:34px;border-bottom:1px solid
#e4e3df;padding-bottom:4px}}.muted{{color:#6b6a66}}table.t{{border-collapse:collapse;font-size:12.5px;margin:8px 0}}
table.t th{{background:#eef3fa;text-align:left;padding:5px 8px}}table.t td{{padding:4px 8px;border-top:1px solid #eeeeec}}
img{{max-width:100%;margin:10px 0;border:1px solid #eeeeec}}.warn{{background:#fff7e6;border-left:4px solid #e8b04b;
padding:6px 10px;margin:4px 0;font-size:13px}}</style></head><body>
<h1>TriADE report</h1><p class="muted">{now} · TriADE {__version__}</p>
<h2>Methods</h2><p>{html.escape(methods_text(res))}</p>
<h2>Plates</h2>{_table_html(summ)}"""]
    if res.warnings:
        parts.append("<h2>Warnings</h2>" + "".join(f'<div class="warn">{html.escape(w)}</div>' for w in res.warnings))
    parts.append(f"<h2>HPLC-MS triage</h2><p>Worklist scenario: <b>{html.escape(scen)}</b> — "
                 f"{len(res.worklist)} wells to inject" + (f"; {n_rev} wells reviewed by hand" if n_rev else "") +
                 f".</p>{_table_html(tri)}")
    for pid, items in figs.items():
        parts.append(f"<h2>Plate {html.escape(str(pid))}</h2>")
        for title, p in items:
            parts.append(f"<p class='muted'>{html.escape(title)}</p><img src='data:image/png;base64,{_png(p)}'>")
    parts.append("<p class='muted' style='margin-top:40px'>Generated by TriADE (MIT license). Settings in "
                 "settings.json; full tables in All_results.xlsx.</p></body></html>")
    h = out_dir / "Report.html"
    h.write_text("".join(parts), encoding="utf-8")

    # ---- PDF
    from matplotlib.backends.backend_pdf import PdfPages
    import textwrap
    pdf = out_dir / "Report.pdf"
    with PdfPages(pdf) as pp:
        fig = plt.figure(figsize=(8.27, 11.69))
        renderer = fig.canvas.get_renderer()

        def put(y, text, **kw):
            """Writes text at y (figure fraction, top-aligned) and returns the y just below it."""
            tx = fig.text(0.08, y, text, va="top", **kw)
            bb = tx.get_window_extent(renderer=renderer).transformed(fig.transFigure.inverted())
            return bb.y0 - 0.012

        y = put(0.96, "TriADE report", fontsize=18, fontweight="bold")
        y = put(y + 0.004, f"{now} · TriADE {__version__}", fontsize=9, color="#6b6a66") - 0.01
        y = put(y, "Methods", fontsize=12, fontweight="bold")
        y = put(y, textwrap.fill(methods_text(res), 110), fontsize=8.2, linespacing=1.45) - 0.01
        y = put(y, "HPLC-MS triage", fontsize=12, fontweight="bold")
        y = put(y, textwrap.fill(f"Worklist: {scen} — {len(res.worklist)} wells to inject"
                                 + (f"; {n_rev} wells reviewed by hand." if n_rev else "."), 110), fontsize=8.5)
        cells = tri.astype(str).values.tolist()
        for row in cells:
            row[0] = textwrap.fill(row[0], 26)
        heads = [textwrap.fill(c.replace("_", " "), 11) for c in tri.columns]
        h = 0.03 * (len(cells) + 1.6)
        ax = fig.add_axes([0.08, y - h, 0.84, h]); ax.axis("off")
        tb = ax.table(cellText=cells, colLabels=heads, loc="upper left", cellLoc="left",
                      colWidths=[0.26] + [0.74 / (len(heads) - 1)] * (len(heads) - 1))
        tb.auto_set_font_size(False); tb.set_fontsize(6.3); tb.scale(1, 2.3)
        if res.warnings:
            fig.text(0.08, 0.04, textwrap.fill("Warnings: " + " | ".join(res.warnings), 140)[:1500], fontsize=6.2,
                     color="#8a2b2a", va="bottom")
        pp.savefig(fig); plt.close(fig)
        fig = plt.figure(figsize=(11.69, 8.27))
        ax = fig.add_axes([0.04, 0.1, 0.92, 0.8]); ax.axis("off")
        ax.set_title("Plates", loc="left", fontsize=13, fontweight="bold")
        tb = ax.table(cellText=[[textwrap.fill(v, 22) for v in row] for row in summ.astype(str).values.tolist()],
                      colLabels=[textwrap.fill(c.replace("_", " "), 10) for c in summ.columns],
                      loc="upper left", cellLoc="left")
        tb.auto_set_font_size(False); tb.set_fontsize(6.5); tb.scale(1, 3.0)
        pp.savefig(fig); plt.close(fig)
        for pid, items in figs.items():
            for title, p in items:
                img = plt.imread(p)
                hgt, wid = img.shape[:2]
                fig = plt.figure(figsize=(11.69, 8.27))
                ax = fig.add_axes([0.03, 0.03, 0.94, 0.9]); ax.imshow(img); ax.axis("off")
                fig.text(0.03, 0.96, f"Plate {pid} — {title}", fontsize=11, fontweight="bold")
                pp.savefig(fig, dpi=150); plt.close(fig)
    return [h, pdf]
