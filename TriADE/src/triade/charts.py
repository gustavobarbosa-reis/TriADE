"""Interactive charts for the browser interface, as Vega-Lite specifications (drawn by Streamlit's
st.vega_lite_chart; no extra library). Each function returns a plain dict."""
from __future__ import annotations

import numpy as np

from . import spectra
from .settings import LEVEL_LABELS

def _config(style) -> dict:
    k = float(getattr(style, "font_scale", 1.0))
    return {"view": {"stroke": None},
            "axis": {"labelColor": "#6b6a66", "titleColor": "#2b2b29", "gridColor": "#eeeeec",
                     "labelFontSize": 10 * k, "titleFontSize": 11 * k},
            "legend": {"labelFontSize": 10 * k, "titleFontSize": 11 * k},
            "text": {"fontSize": 11 * k}}


def _num(v):
    return None if v is None or (isinstance(v, float) and not np.isfinite(v)) else float(v)


def plate_map(res, plate, style, selected: str | None = None) -> dict:
    """Clickable plate map coloured by confidence level (selection name 'well')."""
    vals = []
    for r in res.itertuples(index=False):
        row, col = plate.split(r.Well)
        lv = int(r.Confidence_Level)
        vals.append({"well": r.Well, "row": row, "col": col, "level": lv,
                     "label": f"{lv if lv else '×'} {LEVEL_LABELS[lv]}", "code": str(getattr(r, "Sample_Code", "")),
                     "sel": r.Well == selected})
    levels = [5, 4, 3, 2, 1, 0]
    size = 22 if plate.n_cols <= 12 else 16 if plate.n_cols <= 24 else 10
    return {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "data": {"values": vals}, "config": _config(style),
        "width": size * plate.n_cols, "height": size * plate.n_rows,
        "params": [{"name": "well", "select": {"type": "point", "fields": ["well"], "on": "click"}}],
        "mark": {"type": "rect", "stroke": "#000000", "strokeWidth": 0.4},
        "encoding": {
            "x": {"field": "col", "type": "ordinal", "sort": list(range(1, plate.n_cols + 1)),
                  "axis": {"orient": "top", "title": None, "labelAngle": 0}},
            "y": {"field": "row", "type": "ordinal", "sort": list(plate.rows), "axis": {"title": None}},
            "color": {"field": "label", "type": "nominal", "title": "Level",
                      "scale": {"domain": [f"{k if k else '×'} {LEVEL_LABELS[k]}" for k in levels],
                                "range": [style.level(k) for k in levels]}},
            "strokeWidth": {"condition": {"test": "datum.sel", "value": 3}, "value": 0.4},
            "tooltip": [{"field": "well", "title": "Well"}, {"field": "code", "title": "Sample"},
                        {"field": "label", "title": "Level"}]}}


def trace(tr: dict, style, selected_scan: int | None = None) -> dict:
    """Ejection trace: TIC and the product-ion signals around the well; each scan is a clickable point
    (selection name 'scan')."""
    rows = []
    for i, (t, sc) in enumerate(zip(tr["t"], tr["scan"])):
        rows.append({"t": float(t), "scan": int(sc), "series": "TIC", "value": float(tr["tic"][i])})
        for name, y in tr["ions"].items():
            rows.append({"t": float(t), "scan": int(sc), "series": f"{name} (product)", "value": float(y[i])})
    series = ["TIC"] + [f"{n} (product)" for n in tr["ions"]]
    colors = ["#b0afaa", style["line"], style["line2"], "#6b6a66", "#2b2b29"][:len(series)]
    labels = [{"t": t, "label": w} for t, w in tr["labels"]]
    sel_t = None
    if selected_scan is not None:
        sel_t = next((r["t"] for r in rows if r["scan"] == selected_scan), None)
    layers = [
        {"data": {"values": [{"a": tr["window"][0], "b": tr["window"][1]}]},
         "mark": {"type": "rect", "color": style["band"], "opacity": 0.8},
         "encoding": {"x": {"field": "a", "type": "quantitative"}, "x2": {"field": "b"}}},
        {"mark": {"type": "line", "strokeWidth": 1.4},
         "encoding": {"x": {"field": "t", "type": "quantitative", "title": "time (s)", "scale": {"zero": False}},
                      "y": {"field": "value", "type": "quantitative", "title": "intensity (cps)",
                            "scale": {"type": "symlog", "constant": 10000}},
                      "color": {"field": "series", "type": "nominal", "title": None,
                                "scale": {"domain": series, "range": colors}}}},
        {"params": [{"name": "scan", "select": {"type": "point", "fields": ["scan"], "on": "click", "nearest": True}}],
         "mark": {"type": "point", "filled": True, "size": 40},
         "encoding": {"x": {"field": "t", "type": "quantitative"}, "y": {"field": "value", "type": "quantitative"},
                      "color": {"field": "series", "type": "nominal"},
                      "tooltip": [{"field": "series"}, {"field": "t", "format": ".2f", "title": "time (s)"},
                                  {"field": "scan"}, {"field": "value", "format": ",.0f", "title": "cps"}]}},
        {"data": {"values": labels}, "mark": {"type": "text", "dy": 4, "color": "#2b2b29", "baseline": "top"},
         "encoding": {"x": {"field": "t", "type": "quantitative"}, "y": {"value": 0}, "text": {"field": "label"}}},
    ]
    if sel_t is not None:
        layers.append({"data": {"values": [{"t": sel_t}]}, "mark": {"type": "rule", "strokeDash": [4, 3], "color": "#2b2b29"},
                       "encoding": {"x": {"field": "t", "type": "quantitative"}}})
    return {"$schema": "https://vega.github.io/schema/vega-lite/v5.json", "data": {"values": rows},
            "config": _config(style), "height": 230, "width": "container", "layer": layers}


def spectrum(mz, it, marks, style, cutoff: float, view=None, m_obs: dict | None = None) -> dict:
    """Stick spectrum with the calculated m/z of the product ions (dashed), the theoretical isotope peaks
    (bars scaled to the observed M), the reagent ions (dotted) and the cutoff. Drag to zoom, double-click
    to reset."""
    lo, hi = view or spectra.view_range(marks)
    sel = (mz >= lo - 50) & (mz <= hi + 50)
    peaks = [{"mz": float(a), "i": float(b)} for a, b in zip(mz[sel], it[sel])]
    ions = [{"mz": m["mz"], "label": f"{m['label']} {m['mz']:.3f}"} for m in marks if m["kind"] == "ion"]
    reag = [{"mz": m["mz"], "label": m["label"]} for m in marks if m["kind"] == "reagent" and lo - 50 <= m["mz"] <= hi + 50]
    iso = [{"mz": m["mz"], "i": m["rel"] * (m_obs or {}).get(m["ion"], 0.0), "label": m["label"],
            "pct": 100 * m["rel"]} for m in marks if m["kind"] == "isotope" and (m_obs or {}).get(m["ion"], 0) > 0]
    x = {"field": "mz", "type": "quantitative", "title": "m/z", "scale": {"domain": [lo, hi]}}
    layers = [
        {"data": {"values": peaks}, "params": [{"name": "zoom", "select": "interval", "bind": "scales"}],
         "mark": {"type": "rule", "color": style["line"], "strokeWidth": 1.2},
         "encoding": {"x": x, "y": {"field": "i", "type": "quantitative", "title": "intensity (cps)"},
                      "y2": {"datum": 0},
                      "tooltip": [{"field": "mz", "format": ".4f", "title": "m/z"},
                                  {"field": "i", "format": ",.0f", "title": "cps"}]}},
        {"data": {"values": ions}, "mark": {"type": "rule", "strokeDash": [6, 4], "color": style["line2"],
                                            "strokeWidth": 1.4},
         "encoding": {"x": {"field": "mz", "type": "quantitative"},
                      "tooltip": [{"field": "label", "title": "calculated"}]}},
        {"data": {"values": ions}, "mark": {"type": "text", "angle": 270, "align": "left", "dx": 4, "dy": -4,
                                            "color": style["line2"], "fontSize": 11 * float(getattr(style, "font_scale", 1.0))},
         "encoding": {"x": {"field": "mz", "type": "quantitative"}, "y": {"value": 6}, "text": {"field": "label"}}},
        {"data": {"values": [{"y": cutoff}]}, "mark": {"type": "rule", "strokeDash": [2, 2], "color": "#2b2b29"},
         "encoding": {"y": {"field": "y", "type": "quantitative"}, "tooltip": [{"field": "y", "title": "cutoff"}]}},
    ]
    if reag:
        layers.append({"data": {"values": reag}, "mark": {"type": "rule", "strokeDash": [1, 3], "color": "#6b6a66"},
                       "encoding": {"x": {"field": "mz", "type": "quantitative"},
                                    "tooltip": [{"field": "label", "title": "reagent ion"},
                                                {"field": "mz", "format": ".3f"}]}})
    if iso:
        layers.append({"data": {"values": iso}, "mark": {"type": "tick", "orient": "horizontal", "thickness": 3, "size": 18,
                                                         "color": style["line2"]},
                       "encoding": {"x": {"field": "mz", "type": "quantitative"},
                                    "y": {"field": "i", "type": "quantitative"},
                                    "tooltip": [{"field": "label", "title": "theoretical"},
                                                {"field": "pct", "format": ".1f", "title": "% of M"}]}})
    return {"$schema": "https://vega.github.io/schema/vega-lite/v5.json", "config": _config(style),
            "height": 320, "width": "container", "layer": layers}
