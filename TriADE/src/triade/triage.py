"""HPLC-MS triage: which wells still need HPLC-MS under each scenario, how much time and mobile phase
that saves and, when HPLC-MS results are available, which errors each scenario makes."""
from __future__ import annotations

import pandas as pd

SCENARIOS = [
    ("0: everything to HPLC-MS (no triage)", set(), set()),
    ("1: Echo-MS only rules out level 1", set(), {1}),
    ("2: Echo-MS only rules out levels 1-2", set(), {1, 2}),
    ("3: level 5 = present; level 1 = absent", {5}, {1}),
    ("4: level 5 = present; levels 1-2 = absent", {5}, {1, 2}),
    ("5: levels 4-5 = present; levels 1-2 = absent", {4, 5}, {1, 2}),
]


def samples_only(r: pd.DataFrame) -> pd.DataFrame:
    """Control wells (positive, negative, blank) are not samples: they are left out of the triage."""
    return r[r["Control"].fillna("") == ""] if "Control" in r.columns else r


def triage_table(results: dict, echo_minutes: float, s) -> pd.DataFrame:
    """results: plate id -> result frame (with Confidence_Level and, optionally, HPLC_Product)."""
    results = {pid: samples_only(r) for pid, r in results.items()}
    total = sum(len(r) for r in results.values())
    has_ref = all("HPLC_Product" in r.columns for r in results.values())
    rows = []
    for name, pres, aus in SCENARIOS:
        row = {"Scenario": name}
        runs = fp = fn = fn30 = 0
        for pid, r in results.items():
            lv = r["Confidence_Level"]
            sp, sa = lv.isin(pres), lv.isin(aus)
            k = int((~(sp | sa)).sum())
            row[f"Plate_{pid}_HPLC_injections"] = k
            runs += k
            if has_ref:
                y = r["HPLC_Product"] == 1
                a, b = int((sp & ~y).sum()), int((sa & y).sum())
                row[f"Plate_{pid}_false_present"], row[f"Plate_{pid}_products_missed"] = a, b
                fp += a; fn += b
                fn30 += int((sa & y & (r["HPLC_Conversion"].fillna(0) >= 30)).sum())
        saved = total - runs
        row.update({
            "HPLC_injections": runs, "HPLC_injections_%_of_wells": round(100 * runs / total, 1) if total else 0,
            "HPLC_time_h": round(runs * s.hplc_run_min / 60, 1),
            "Mobile_phase_used_mL": round(runs * s.hplc_run_min * s.hplc_flow_ml_min),
            "Injections_saved": saved, "Saved_%": round(100 * saved / total, 1) if total else 0,
            "Time_saved_h": round(saved * s.hplc_run_min / 60, 1),
            "Net_time_saved_h_minus_EchoMS": round((saved * s.hplc_run_min - (echo_minutes if saved else 0)) / 60, 1),
            "Mobile_phase_saved_mL": round(saved * s.hplc_run_min * s.hplc_flow_ml_min)})
        if has_ref:
            row.update({"False_present": fp, "False_present_%_of_wells": round(100 * fp / total, 1),
                        "Products_missed": fn, "Products_missed_%_of_wells": round(100 * fn / total, 1),
                        "Products_missed_conversion_ge_30%": fn30})
        rows.append(row)
    t = pd.DataFrame(rows)
    t.insert(1, "Assumptions",
             f"{total} wells; {s.hplc_run_min:g} min per HPLC-MS injection at {s.hplc_flow_ml_min:g} mL/min; "
             f"Echo-MS time {echo_minutes:.1f} min in total; no blanks, equilibration or standards"
             + ("; wells without an HPLC-MS result count as not detected" if has_ref else ""))
    return t


def discrepancies(results: dict) -> pd.DataFrame:
    keep = ["Plate", "Well", "Sample_Code", "Confidence_Level", "Confidence_Label", "Level_Note", "Evidence_Intensity",
            "Raw_Ejection_Intensity", "Previous_Ejection_Intensity", "Detected", "HPLC_Result", "HPLC_Conversion",
            "Ions_Detected", "Echo_Product_Ratio_%", "SM_Interference", "Reason"]
    parts = []
    for r in results.values():
        if "HPLC_Product" not in r.columns:
            continue
        y = r["HPLC_Product"] == 1
        d = r.copy()
        d["Discrepancy"] = ""
        d.loc[(d["Confidence_Level"] >= 4) & ~y, "Discrepancy"] = "false present (level 4-5, no product in HPLC-MS)"
        d.loc[(d["Confidence_Level"] <= 2) & y, "Discrepancy"] = "product missed (level 1-2, product in HPLC-MS)"
        parts.append(d[d["Discrepancy"] != ""][["Discrepancy"] + [c for c in keep if c in d.columns]])
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts).sort_values(["Discrepancy", "Plate", "Confidence_Level"], ascending=[True, True, False])


def hplc_worklist(results: dict, scenario: int) -> pd.DataFrame:
    """The wells that still need HPLC-MS under a scenario, in plate order (a ready-to-use injection list).
    Manual review decisions take precedence: 'present'/'absent' leave the list, 'hplc' joins it."""
    name, pres, aus = SCENARIOS[scenario]
    parts = []
    cols = ["Plate", "Well", "Sample_Code", "Product_SMILES", "Confidence_Level", "Confidence_Label", "Level_Note"]
    for r in results.values():
        r = samples_only(r)
        need = ~(r["Confidence_Level"].isin(pres) | r["Confidence_Level"].isin(aus))
        if "Review_Decision" in r.columns:
            dec = r["Review_Decision"].fillna("")
            need = (need & ~dec.isin(["present", "absent"])) | (dec == "hplc")
        sub = r.loc[need, cols + (["Review_Decision"] if "Review_Decision" in r.columns else [])]
        parts.append(sub)
    out = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=cols)
    out.insert(0, "HPLC_Injection", range(1, len(out) + 1))
    return out
