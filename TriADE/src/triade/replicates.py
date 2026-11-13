"""Replicate ejections of the same well, and plate controls.

Replicates: each replicate is analysed on its own; the well is then called detected when 'any', the
'majority' or 'all' replicates detect it (replicate_rule), and its level is the median of the replicate
levels (the lower one for an even number), moved to 3 'Check spectrum' when it contradicts the call.

Controls: wells marked positive (product known to be there), negative (no product) or blank. A plate is
accepted when enough positive controls are detected and few enough negatives/blanks are.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .settings import LEVEL_LABELS


def _call(lv: int) -> str:
    return "DETECTED" if lv >= 4 else "CHECK SPECTRUM" if lv == 3 else "NOT MEASURABLE" if lv == 0 else "NOT DETECTED"


def combine(frames: dict, s) -> pd.DataFrame:
    """frames: replicate number -> result frame of that replicate. Returns one row per well."""
    reps = sorted(frames)
    allr = pd.concat([frames[r] for r in reps], ignore_index=True)
    rows = []
    for well, g in allr.groupby("Well", sort=False):
        g = g.sort_values("Replicate")
        meas = g[g["Confidence_Level"] > 0]
        n = len(meas)
        n_det = int(meas["Detected"].sum())
        if s.replicate_rule == "any":
            det = n_det >= 1
        elif s.replicate_rule == "all":
            det = n > 0 and n_det == n
        else:
            det = n_det > n / 2
        if n:
            levels = sorted(meas["Confidence_Level"].astype(int))
            lv = levels[(n - 1) // 2]
        else:
            lv = 0
        notes = []
        if det and lv <= 2:
            lv, notes = 3, ["replicates disagree (detected by the replicate rule, low median level)"]
        elif not det and lv >= 4:
            lv, notes = 3, ["replicates disagree (median level high, not detected by the replicate rule)"]
        cand = g[g["Confidence_Level"] == lv] if (g["Confidence_Level"] == lv).any() else g
        rep_row = cand.sort_values("Evidence_Intensity", ascending=False, na_position="last").iloc[0].copy()
        ev = meas["Evidence_Intensity"].astype(float)
        rep_row["Confidence_Level"] = lv
        rep_row["Confidence_Label"] = LEVEL_LABELS[lv]
        rep_row["Call"] = _call(lv)
        rep_row["Detected"] = int(det)
        rep_row["Level_Note"] = "; ".join([x for x in [str(rep_row.get("Level_Note", "") or "")] + notes if x])
        rep_row["Replicates"] = len(g)
        rep_row["Replicate_Levels"] = ";".join(str(int(v)) for v in g["Confidence_Level"])
        rep_row["Replicate_Detected"] = ";".join(str(int(v)) for v in g["Detected"])
        rep_row["Replicate_Evidence_CV_%"] = (round(100 * ev.std() / ev.mean(), 1)
                                              if n > 1 and ev.mean() > 0 else np.nan)
        rep_row["Replicate"] = "all"
        rows.append(rep_row)
    out = pd.DataFrame(rows).reset_index(drop=True)
    return out


def fleiss_kappa(frames: dict) -> float:
    """Agreement of the binary calls between replicates (Fleiss' kappa, wells with all replicates)."""
    reps = sorted(frames)
    m = pd.concat([frames[r].set_index("Well")["Detected"].rename(r) for r in reps], axis=1).dropna()
    if len(reps) < 2 or m.empty:
        return float("nan")
    k = len(reps)
    yes = m.sum(axis=1).to_numpy(float)
    no = k - yes
    p_i = (yes * (yes - 1) + no * (no - 1)) / (k * (k - 1))
    p_bar = p_i.mean()
    pj = np.array([yes.sum(), no.sum()]) / (len(m) * k)
    pe = float((pj ** 2).sum())
    return round(float((p_bar - pe) / (1 - pe)), 3) if pe < 1 else float("nan")


def replicate_summary(frames: dict, combined: pd.DataFrame) -> dict:
    if len(frames) < 2:
        return {}
    cv = combined["Replicate_Evidence_CV_%"]
    agree = combined["Replicate_Detected"].map(lambda v: len(set(v.split(";"))) == 1)
    return {"Replicates": len(frames), "Replicate_kappa_Fleiss": fleiss_kappa(frames),
            "Replicate_call_agreement_%": round(100 * agree.mean(), 1),
            "Replicate_evidence_CV_median_%": round(float(cv.median()), 1) if cv.notna().any() else np.nan}


# ------------------------------------------------------------------------------------------
# Controls
# ------------------------------------------------------------------------------------------
def control_check(res: pd.DataFrame, s) -> dict:
    """Detection of the control wells and the plate acceptance."""
    c = res[res["Control"].fillna("") != ""] if "Control" in res.columns else res.iloc[:0]
    if c.empty:
        return {}
    out = {}
    pos = c[c["Control"] == "positive"]
    neg = c[c["Control"].isin(["negative", "blank"])]
    if len(pos):
        out["Positive_controls"] = len(pos)
        out["Positive_detected_%"] = round(100 * pos["Detected"].mean(), 1)
    if len(neg):
        out["Negative_blank_controls"] = len(neg)
        out["Negative_blank_detected_%"] = round(100 * neg["Detected"].mean(), 1)
        out["Negative_blank_median_raw_cps"] = float(neg["Raw_Ejection_Intensity"].median())
    reasons = []
    if len(pos) and out["Positive_detected_%"] < s.control_min_positive_pct:
        reasons.append(f"{out['Positive_detected_%']:g}% of positive controls detected (< {s.control_min_positive_pct:g}%)")
    if len(neg) and out["Negative_blank_detected_%"] > s.control_max_negative_pct:
        reasons.append(f"{out['Negative_blank_detected_%']:g}% of negative controls/blanks detected "
                       f"(> {s.control_max_negative_pct:g}%)")
    out["Plate_accepted"] = "no" if reasons else "yes"
    out["Acceptance_note"] = "; ".join(reasons)
    return out
