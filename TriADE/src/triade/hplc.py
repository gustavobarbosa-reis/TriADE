"""Optional HPLC-MS (LC-MS) results, used to check the Echo-MS calls.

Accepted files (CSV/TSV/Excel), one per plate or one with a Plate column:
  * a simple table: Well | Conversion_% (or Product_%), optionally Plate and Indeterminate (1/0);
  * a table with Well | Product (1/0, yes/no) when only presence is known;
  * an instrument export with a Position/Well column and product/(product + reagent) "Ratio" columns in %
    (the largest ratio is used) and, optionally, peak "Area" columns: a well with every area at zero had
    no HPLC-MS signal at all.
Wells without an HPLC-MS result (no signal, marked indeterminate, or missing from the file) count as
"not detected" in every comparison.
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .plates import PlateFormat
from .samples import ALIASES, key, read_table

log = logging.getLogger("triade")

CONV_KEYS = ("conversion", "conversion%", "conversionpct", "conversionpercent", "product%", "productpct",
             "productpercent", "yield", "yield%", "conversao")
PRODUCT_KEYS = ("product", "productpresent", "present", "detected", "hplcproduct", "lcmsproduct")
IND_KEYS = ("indeterminate", "indeterminado", "excluded", "exclude")


def _plate_from_name(name: str):
    m = re.search(r"plate\s*[_-]?\s*([A-Za-z0-9]+)", name, re.I)
    return m.group(1) if m else None


def _as_bool(v):
    if pd.isna(v):
        return np.nan
    s = str(v).strip().lower()
    if s in ("1", "1.0", "yes", "y", "true", "sim", "present", "detected"):
        return 1.0
    if s in ("0", "0.0", "no", "n", "false", "nao", "não", "absent", "not detected"):
        return 0.0
    return np.nan


def load_hplc(paths, plate: PlateFormat, present_pct: float, plate_ids=None) -> pd.DataFrame:
    """Returns Plate, Well, HPLC_Conversion (%), HPLC_No_Signal, HPLC_Indeterminate (+ UV/MS ratios when present).
    plate_ids: optional list with the plate of each file (when the file has no Plate column)."""
    frames = []
    for n, path in enumerate(paths):
        path = Path(path)
        t = read_table(path)
        t.columns = [str(c).strip() for c in t.columns]
        keys = {key(c): c for c in t.columns}
        well = next((keys[k] for k in ALIASES["Well"] if k in keys), None)
        if well is None:
            raise ValueError(f"{path.name}: the well column (Well or Position) is missing")
        pc = next((keys[k] for k in ALIASES["Plate"] if k in keys), None)
        given = plate_ids[n] if plate_ids and n < len(plate_ids) and plate_ids[n] else None
        if pc is not None:
            plates = t[pc].map(lambda v: str(int(v)) if isinstance(v, (int, float)) and float(v).is_integer() else str(v))
        else:
            pid = given or _plate_from_name(path.name) or "1"
            plates = pd.Series([str(pid)] * len(t))
        out = pd.DataFrame({"Plate": plates.values, "Well": t[well].map(plate.parse).values})
        conv_col = next((keys[k] for k in CONV_KEYS if k in keys), None)
        prod_col = next((keys[k] for k in PRODUCT_KEYS if k in keys), None)
        ratios = [c for c in t.columns if "ratio" in c.lower()]
        if conv_col is not None:
            out["HPLC_Conversion"] = pd.to_numeric(t[conv_col], errors="coerce").values
        elif ratios:
            for c in ratios:
                k = key(c)
                name = "HPLC_UV_ratio" if "uv" in k else "HPLC_MS_ratio" if "ms" in k else f"HPLC_{k}"
                out[name] = pd.to_numeric(t[c], errors="coerce").values
            out["HPLC_Conversion"] = out[[c for c in out.columns if c.startswith("HPLC_") and "ratio" in c]].max(axis=1)
        elif prod_col is not None:
            # presence only: 100 = product, 0 = no product
            out["HPLC_Conversion"] = t[prod_col].map(_as_bool).values * 100.0
        else:
            raise ValueError(f"{path.name}: no Conversion_%, Product or 'Ratio' column found")
        areas = [c for c in t.columns if "area" in c.lower() and "%" not in c and c not in ratios]
        if areas:
            a = t[areas].apply(pd.to_numeric, errors="coerce").fillna(0)
            out["HPLC_No_Signal"] = (a == 0).all(axis=1).astype(int).values
        else:
            out["HPLC_No_Signal"] = 0
        ind_col = next((keys[k] for k in IND_KEYS if k in keys), None)
        out["HPLC_Indeterminate"] = (t[ind_col].map(_as_bool).fillna(0).astype(int).values if ind_col else 0)
        out.loc[(out["HPLC_No_Signal"] == 1) | (out["HPLC_Indeterminate"] == 1), "HPLC_Conversion"] = np.nan
        out = out.dropna(subset=["Well"])
        frames.append(out)
        log.info("HPLC-MS %s: plate(s) %s, %d wells, %d with no signal", path.name,
                 sorted(out["Plate"].unique()), len(out), int(out["HPLC_No_Signal"].sum()))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def parse_exclusions(items, plate: PlateFormat) -> set:
    """'2:9' = plate 2, column 9; '2:A9' = plate 2, well A9; '1:C' = plate 1, row C (one per item or
    separated by commas/semicolons/new lines)."""
    ex = set()
    for it in items or []:
        for tok in re.split(r"[,;\n]", str(it)):
            tok = tok.split("#")[0].strip()
            if not tok:
                continue
            if ":" not in tok:
                raise ValueError(f"exclusion {tok!r}: use plate:column, plate:well or plate:row (e.g. 2:9, 2:A9, 1:C)")
            pl, what = (x.strip() for x in tok.split(":", 1))
            what = what.upper()
            if what.isdigit():
                ex |= {(pl, f"{r}{int(what)}") for r in plate.rows}
            elif what in plate.rows:
                ex |= {(pl, f"{what}{c}") for c in range(1, plate.n_cols + 1)}
            else:
                w = plate.parse(what)
                if w is None:
                    raise ValueError(f"exclusion {tok!r}: {what} is not a column, row or well of a {plate.name}-well plate")
                ex.add((pl, w))
    return ex


def mark_indeterminate(hplc: pd.DataFrame, excl: set) -> pd.DataFrame:
    if hplc.empty or not excl:
        return hplc
    h = hplc.copy()
    m = [(p, w) in excl for p, w in zip(h["Plate"], h["Well"])]
    h.loc[m, "HPLC_Indeterminate"] = 1
    h.loc[m, "HPLC_Conversion"] = np.nan
    log.info("HPLC-MS: %d wells marked as indeterminate (counted as not detected).", int(sum(m)))
    return h


def add_hplc(res: pd.DataFrame, hplc: pd.DataFrame, plate_id: str, s) -> dict:
    """Adds the HPLC-MS columns, Reference and Agreement to the plate result; returns the metrics."""
    lc = hplc[hplc["Plate"] == plate_id].drop(columns="Plate").drop_duplicates("Well").set_index("Well")
    if lc.empty:
        return {}
    for c in lc.columns:
        res[c] = res["Well"].map(lc[c])
    conv = res["HPLC_Conversion"]
    res["HPLC_Product"] = (conv.fillna(0) >= s.hplc_product_min_pct).astype(int)
    ind = res["HPLC_Indeterminate"].fillna(0) == 1
    nos = res["HPLC_No_Signal"].fillna(0) == 1
    res["HPLC_Result"] = np.select(
        [ind, nos, conv.isna(), conv >= s.hplc_product_min_pct],
        ["indeterminate (counted as not detected)", "no signal (counted as not detected)",
         "no data (counted as not detected)", "product"], "no product")
    ref = s.reference_label
    yy, pp = res["HPLC_Product"].to_numpy(), res["Detected"].to_numpy()
    res["Reference"] = yy
    res["Agreement"] = np.where(yy == pp, np.where(pp == 1, "both", "neither"),
                                np.where(pp == 1, "Echo-MS only", f"{ref} only"))
    samples = (res["Control"].fillna("") == "").to_numpy() if "Control" in res.columns else np.ones(len(res), bool)
    y, p = yy[samples], pp[samples]          # controls are not part of the comparison
    ind, nos, conv = ind[samples], nos[samples], conv[samples]
    po = float((y == p).mean())
    pe = y.mean() * p.mean() + (1 - y.mean()) * (1 - p.mean())
    d = {"Reference_wells": len(y), "Reference_with_product": int(y.sum()),
         "Agreement_%": round(100 * po, 1), "Cohen_kappa": round((po - pe) / (1 - pe), 3) if pe < 1 else np.nan,
         "Both_detected": int(((y == 1) & (p == 1)).sum()), "Neither": int(((y == 0) & (p == 0)).sum()),
         "EchoMS_only": int(((y == 0) & (p == 1)).sum()), f"{ref}_only": int(((y == 1) & (p == 0)).sum()),
         f"{ref}_no_signal": int(nos.sum()), f"{ref}_indeterminate": int(ind.sum())}
    rs = res[samples]
    for k in (5, 4, 3, 2, 1):
        sel = rs["Confidence_Level"] == k
        d[f"Level_{k}_with_product_%"] = round(100 * rs.loc[sel, "HPLC_Product"].mean(), 1) if sel.any() else np.nan
    d["Spearman_level_vs_conversion"] = round(rs["Confidence_Level"].corr(conv.fillna(0), method="spearman"), 3)
    if "Echo_Product_Ratio_%" in rs.columns and rs["Echo_Product_Ratio_%"].notna().sum() > 2:
        m = rs["Echo_Product_Ratio_%"].notna()
        d["Spearman_EchoRatio_vs_conversion"] = round(
            rs.loc[m, "Echo_Product_Ratio_%"].corr(conv[m].fillna(0), method="spearman"), 3)
    return d
