"""Manual review of wells (from the spectrum viewer) and the Excel files for moving the wells that still need
HPLC-MS: a transfer list (source well -> destination plate/well) and an injection sequence."""
from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd

from .plates import get_plate_format

DECISIONS = {"present": "Confirmed present (no HPLC-MS)", "absent": "Confirmed absent (no HPLC-MS)",
             "hplc": "Send to HPLC-MS"}
COLUMNS = ["Plate", "Well", "Decision", "Comment", "Reviewed_at"]


def empty() -> pd.DataFrame:
    return pd.DataFrame(columns=COLUMNS)


def load(path) -> pd.DataFrame:
    t = pd.read_csv(path, dtype=str).fillna("") if str(path).lower().endswith(".csv") else \
        pd.read_excel(path, dtype=str).fillna("")
    missing = [c for c in ("Plate", "Well", "Decision") if c not in t.columns]
    if missing:
        raise ValueError(f"{Path(path).name}: review file without the column(s) {', '.join(missing)}")
    bad = sorted(set(t["Decision"]) - set(DECISIONS) - {""})
    if bad:
        raise ValueError(f"{Path(path).name}: unknown decision(s) {', '.join(bad)} (use present, absent or hplc)")
    for c in COLUMNS:
        if c not in t.columns:
            t[c] = ""
    return t[COLUMNS][t["Decision"] != ""].drop_duplicates(["Plate", "Well"], keep="last").reset_index(drop=True)


def set_decision(review: pd.DataFrame, plate: str, well: str, decision: str, comment: str = "") -> pd.DataFrame:
    """Adds or replaces the decision of a well ('' removes it)."""
    r = review[~((review["Plate"] == str(plate)) & (review["Well"] == well))]
    if decision:
        if decision not in DECISIONS:
            raise ValueError(f"unknown decision {decision!r}")
        r = pd.concat([r, pd.DataFrame([{"Plate": str(plate), "Well": well, "Decision": decision, "Comment": comment,
                                         "Reviewed_at": datetime.now().strftime("%Y-%m-%d %H:%M")}])],
                      ignore_index=True)
    return r.reset_index(drop=True)


def apply(results: dict, review: pd.DataFrame | None) -> None:
    """Adds Review_Decision / Review_Comment to the result frames (in place)."""
    rv = review if review is not None else empty()
    for pid, r in results.items():
        sub = rv[rv["Plate"] == str(pid)].set_index("Well")
        r["Review_Decision"] = r["Well"].map(sub["Decision"]).fillna("") if len(sub) else ""
        r["Review_Comment"] = r["Well"].map(sub["Comment"]).fillna("") if len(sub) else ""


# ------------------------------------------------------------------------------------------
# Excel export of the HPLC-MS worklist
# ------------------------------------------------------------------------------------------
def transfer_list(worklist: pd.DataFrame, s) -> pd.DataFrame:
    """Each well still needing HPLC-MS -> a destination plate and well (filled row by row)."""
    dest = get_plate_format(s.transfer_plate_format)
    wells = dest.wells()
    rows = []
    for i, r in enumerate(worklist.itertuples(index=False)):
        rows.append({"Transfer": i + 1, "Source_Plate": r.Plate, "Source_Well": r.Well, "Sample_Code": r.Sample_Code,
                     "Destination_Plate": f"HPLC-{i // len(wells) + 1}", "Destination_Well": wells[i % len(wells)],
                     "Volume_nL": s.transfer_volume_nl if s.transfer_volume_nl else ""})
    return pd.DataFrame(rows, columns=["Transfer", "Source_Plate", "Source_Well", "Sample_Code", "Destination_Plate",
                                       "Destination_Well", "Volume_nL"])


def injection_sequence(transfer: pd.DataFrame, s) -> pd.DataFrame:
    """Injection order on the HPLC-MS, with a blank every N samples if asked."""
    rows, n = [], 0
    for r in transfer.itertuples(index=False):
        if s.blank_every and n and n % int(s.blank_every) == 0:
            rows.append({"Type": "Blank", "Sample_Name": "BLANK", "Sample_Code": "", "Source_Plate": "",
                         "Source_Well": "", "Vial_Plate": "", "Vial_Position": ""})
        rows.append({"Type": "Sample", "Sample_Name": r.Sample_Code or f"P{r.Source_Plate}-{r.Source_Well}",
                     "Sample_Code": r.Sample_Code, "Source_Plate": r.Source_Plate, "Source_Well": r.Source_Well,
                     "Vial_Plate": r.Destination_Plate, "Vial_Position": r.Destination_Well})
        n += 1
    seq = pd.DataFrame(rows, columns=["Type", "Sample_Name", "Sample_Code", "Source_Plate", "Source_Well",
                                      "Vial_Plate", "Vial_Position"])
    seq.insert(0, "Injection", range(1, len(seq) + 1))
    seq["Run_min"] = s.hplc_run_min
    return seq


def write_excel(path, worklist: pd.DataFrame, s, scenario_name: str) -> Path:
    from openpyxl.styles import Font, PatternFill
    from openpyxl.utils import get_column_letter
    tr = transfer_list(worklist, s)
    seq = injection_sequence(tr, s)
    notes = pd.DataFrame({"Item": ["Scenario", "Wells to inject", "Destination plate format", "Blank every",
                                   "HPLC-MS time (h)", "Mobile phase (mL)", "Columns"],
                          "Value": [scenario_name, len(worklist), f"{s.transfer_plate_format}-well",
                                    s.blank_every or "none", round(len(worklist) * s.hplc_run_min / 60, 1),
                                    round(len(worklist) * s.hplc_run_min * s.hplc_flow_ml_min),
                                    "Generic layout: copy or map the columns to the import format of your "
                                    "liquid handler and HPLC-MS software."]})
    path = path if hasattr(path, "write") else Path(path)       # a file path or an in-memory buffer
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for name, df in (("Transfer list", tr), ("Injection sequence", seq), ("Wells", worklist), ("Notes", notes)):
            df.to_excel(xw, sheet_name=name, index=False)
            ws = xw.sheets[name]
            for j, col in enumerate(df.columns, 1):
                c = ws.cell(row=1, column=j)
                c.font, c.fill = Font(bold=True), PatternFill("solid", fgColor="D9E6F7")
                width = max([len(str(col))] + [len(str(v)) for v in df[col].head(200)]) + 2
                ws.column_dimensions[get_column_letter(j)].width = min(max(width, 8), 70)
            ws.freeze_panes = "A2"
    return path
