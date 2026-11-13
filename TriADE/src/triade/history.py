"""Project history: one row of quality metrics per plate run, kept in a CSV file across analyses, and
control charts (Shewhart and CUSUM) that show when the method drifts away from its usual behaviour."""
from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

METRICS = {   # column -> (label, direction that is a problem: "low", "high" or "both")
    "Pulse_median": ("Ejection pulse height (median, × baseline)", "low"),
    "Pulse_drift_spearman": ("Pulse drift over the run (Spearman ρ)", "both"),
    "IS_CV_%": ("Internal standard CV (%)", "high"),
    "Detected_%": ("Wells detected (%)", "both"),
    "Level5_%": ("Wells at level 5 (%)", "both"),
    "Positive_detected_%": ("Positive controls detected (%)", "low"),
    "Negative_blank_detected_%": ("Negative controls / blanks detected (%)", "high"),
    "Agreement_%": ("Agreement with HPLC-MS (%)", "low"),
    "Replicate_kappa_Fleiss": ("Replicate agreement (Fleiss κ)", "low"),
}


def run_date(path) -> str:
    """Acquisition date from the mzML (run startTimeStamp), else the file date."""
    try:
        with open(path, "rb") as f:
            head = f.read(20000).decode("utf-8", "ignore")
        m = re.search(r'startTimeStamp="([^"]+)"', head)
        if m:
            return m.group(1)[:19].replace("T", " ")
        return datetime.fromtimestamp(Path(path).stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")
    except OSError:
        return ""


def rows_for_run(res) -> pd.DataFrame:
    """Metrics of each plate of a Results object."""
    rows = []
    for _, srow in res.summary.iterrows():
        pid = srow["Plate"]
        r = res.results[pid]
        q = res.stats.get("run_qc", {}).get(pid)
        n = max(int((r["Confidence_Level"] > 0).sum()), 1)
        row = {"Run_date": run_date(res.plates[pid].run.path), "Analysed_at": datetime.now().strftime("%Y-%m-%d %H:%M"),
               "mzML": srow["mzML"], "Plate": pid, "Wells": int(srow["Wells"]),
               "Detected_%": round(100 * r["Detected"].sum() / n, 1),
               "Level5_%": round(100 * (r["Confidence_Level"] == 5).sum() / n, 1),
               "Pulse_median": round(float(q["Pulse_Height"].median()), 3) if q is not None and len(q) else np.nan}
        for c in ("Pulse_drift_spearman", "IS_CV_%", "Positive_detected_%", "Negative_blank_detected_%",
                  "Agreement_%", "Replicate_kappa_Fleiss"):
            row[c] = srow[c] if c in srow.index else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def update(path, new_rows: pd.DataFrame) -> pd.DataFrame:
    """Appends the new runs to the history CSV (a run analysed again replaces its earlier row)."""
    path = Path(path)
    old = pd.read_csv(path, dtype={"Plate": str}) if path.exists() and path.stat().st_size else pd.DataFrame()
    h = pd.concat([old, new_rows], ignore_index=True)
    h = h.drop_duplicates(["mzML", "Plate"], keep="last")
    h = h.sort_values(["Run_date", "mzML", "Plate"]).reset_index(drop=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    h.to_csv(path, index=False)
    return h


def control_limits(h: pd.DataFrame, min_baseline: int = 5) -> pd.DataFrame:
    """Per metric and run: centre and ±3 SD limits from the earlier runs (at least min_baseline of them),
    and a standardised CUSUM (k = 0.5, h = 4) to catch small persistent shifts."""
    out = []
    for col, (label, side) in METRICS.items():
        if col not in h.columns or h[col].notna().sum() < 2:
            continue
        cpos = cneg = 0.0
        for i, v in enumerate(h[col]):
            prev = h[col].iloc[:i].dropna()
            row = {"Run": i + 1, "Run_date": h["Run_date"].iloc[i], "mzML": h["mzML"].iloc[i],
                   "Plate": h["Plate"].iloc[i], "Metric": col, "Label": label, "Value": v}
            if len(prev) >= min_baseline and np.isfinite(v):
                mu, sd = prev.mean(), prev.std(ddof=1)
                sd = sd if sd > 0 else max(abs(mu) * 0.01, 1e-9)
                z = (v - mu) / sd
                cpos, cneg = max(0.0, cpos + z - 0.5), min(0.0, cneg + z + 0.5)
                flag = []
                if (z > 3 and side in ("high", "both")) or (z < -3 and side in ("low", "both")):
                    flag.append("outside ±3 SD")
                if (cpos > 4 and side in ("high", "both")) or (cneg < -4 and side in ("low", "both")):
                    flag.append("CUSUM shift")
                row.update({"Centre": mu, "LCL": mu - 3 * sd, "UCL": mu + 3 * sd, "z": round(z, 2),
                            "Flag": "; ".join(flag)})
            else:
                row.update({"Centre": np.nan, "LCL": np.nan, "UCL": np.nan, "z": np.nan,
                            "Flag": "" if np.isfinite(v) else "",})
            out.append(row)
    return pd.DataFrame(out)
