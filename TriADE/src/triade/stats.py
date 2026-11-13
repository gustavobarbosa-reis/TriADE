"""Statistics for the triage:

* Wilson 95% confidence intervals for the share of wells with product at each level and for the triage
  error rates; a bootstrap interval for Cohen's kappa;
* threshold calibration: the lowest 'present' threshold and the highest 'absent' threshold that keep the
  errors under a target, checked plate by plate when there is more than one plate;
* a per-well probability of product (logistic regression on the Echo-MS measurements), validated by
  leaving one plate out, that can be saved and applied to new plates;
* run QC: ejection pulse height over the run (spray drops) and an optional internal standard;
* McNemar's test to compare two analyses of the same wells.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

Z95 = 1.959964


# ------------------------------------------------------------------------------------------
# Confidence intervals
# ------------------------------------------------------------------------------------------
def wilson(k: int, n: int, z: float = Z95) -> tuple[float, float]:
    """Wilson score interval of a proportion k/n, in %."""
    if n <= 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (round(100 * max(0.0, centre - half), 1), round(100 * min(1.0, centre + half), 1))


def kappa(y, p) -> float:
    y, p = np.asarray(y), np.asarray(p)
    po = float((y == p).mean())
    pe = y.mean() * p.mean() + (1 - y.mean()) * (1 - p.mean())
    return (po - pe) / (1 - pe) if pe < 1 else float("nan")


def bootstrap_kappa(y, p, n_boot: int = 2000, seed: int = 0) -> tuple[float, float]:
    y, p = np.asarray(y), np.asarray(p)
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        i = rng.integers(0, len(y), len(y))
        k = kappa(y[i], p[i])
        if np.isfinite(k):
            vals.append(k)
    if not vals:
        return (float("nan"), float("nan"))
    return (round(float(np.percentile(vals, 2.5)), 3), round(float(np.percentile(vals, 97.5)), 3))


def _samples(results: dict) -> dict:
    from .triage import samples_only
    return {pid: samples_only(r) for pid, r in results.items()}


def level_table(results: dict) -> pd.DataFrame:
    """Per plate (and all plates): wells and share with HPLC-MS product at each level, with 95% CI."""
    results = _samples(results)
    rows = []
    frames = [(pid, r) for pid, r in results.items() if "HPLC_Product" in r.columns]
    if not frames:
        return pd.DataFrame()
    if len(frames) > 1:
        frames.append(("all", pd.concat([r for _, r in frames])))
    for pid, r in frames:
        for k in (5, 4, 3, 2, 1, 0):
            sel = r["Confidence_Level"] == k
            n = int(sel.sum())
            if not n:
                continue
            w = int(r.loc[sel, "HPLC_Product"].sum())
            lo, hi = wilson(w, n)
            rows.append({"Plate": pid, "Level": k, "Wells": n, "With_product": w,
                         "With_product_%": round(100 * w / n, 1), "CI95_low_%": lo, "CI95_high_%": hi})
    return pd.DataFrame(rows)


def triage_ci(tri: pd.DataFrame) -> pd.DataFrame:
    """Adds 95% CIs to the false-present and missed-product rates of the triage table."""
    if "False_present" not in tri.columns:
        return tri
    t = tri.copy()
    total = int(round(t.loc[0, "HPLC_injections"] * 100 / t.loc[0, "HPLC_injections_%_of_wells"]))
    for col in ("False_present", "Products_missed"):
        ci = [wilson(int(v), total) for v in t[col]]
        t[f"{col}_CI95_%"] = [f"{lo:g}-{hi:g}" for lo, hi in ci]
    return t


# ------------------------------------------------------------------------------------------
# Threshold calibration
# ------------------------------------------------------------------------------------------
def _present_curve(ev, y, total):
    order = np.unique(ev[ev > 0])[::-1]
    out = []
    for T in order:
        called = ev >= T
        fp = int((called & (y == 0)).sum())
        out.append({"threshold_cps": float(T), "called_present": int(called.sum()), "false_present": fp,
                    "false_present_%_of_wells": round(100 * fp / total, 2)})
    return pd.DataFrame(out)


def _absent_curve(raw, ev, y, level3, total):
    cand = np.unique(raw)
    out = []
    for R in cand:
        called = (raw <= R) & (ev < level3)
        miss = int((called & (y == 1)).sum())
        out.append({"threshold_cps": float(R), "called_absent": int(called.sum()), "products_missed": miss,
                    "products_missed_%_of_wells": round(100 * miss / total, 2)})
    return pd.DataFrame(out)


def calibrate(results: dict, s) -> dict:
    """Suggested thresholds for 'present' (level 5, on the evidence) and 'absent' (levels 1-2, on the raw
    signal) that keep false presents / missed products <= target (% of wells). With several plates the
    threshold found on the other plates is also tested on each plate (out-of-plate check)."""
    frames = {pid: r for pid, r in _samples(results).items() if "HPLC_Product" in r.columns}
    if not frames:
        return {}
    target = float(s.target_false_present_pct)

    def fit(r):
        r = r[r["Confidence_Level"] > 0]
        ev = r["Evidence_Intensity"].fillna(0).to_numpy()
        raw = r["Raw_Ejection_Intensity"].fillna(0).to_numpy()
        y = r["HPLC_Product"].to_numpy()
        total = len(r)
        pc = _present_curve(ev, y, total)
        ok = pc[pc["false_present_%_of_wells"] <= target]
        t_present = float(ok["threshold_cps"].min()) if len(ok) else float("nan")
        ac = _absent_curve(raw, ev, y, s.level3_cps, total)
        ok = ac[ac["products_missed_%_of_wells"] <= target]
        t_absent = float(ok["threshold_cps"].max()) if len(ok) else float("nan")
        return t_present, t_absent, pc, ac

    def apply(r, tp, ta):
        r = r[r["Confidence_Level"] > 0]
        ev = r["Evidence_Intensity"].fillna(0); raw = r["Raw_Ejection_Intensity"].fillna(0)
        y = r["HPLC_Product"]
        pres = ev >= tp if np.isfinite(tp) else ev < -1
        aus = (raw <= ta) & (ev < s.level3_cps) if np.isfinite(ta) else raw < -1
        n = len(r)
        return {"called_present": int(pres.sum()), "false_present": int((pres & (y == 0)).sum()),
                "called_absent": int(aus.sum()), "products_missed": int((aus & (y == 1)).sum()),
                "to_HPLC_MS": int(n - pres.sum() - aus.sum()), "wells": n}

    all_r = pd.concat(frames.values())
    tp, ta, pc, ac = fit(all_r)
    rows = [{"Fitted_on": "all plates", "Tested_on": "all plates (same data)",
             "Present_threshold_cps": tp, "Absent_threshold_cps": ta, **apply(all_r, tp, ta)}]
    if len(frames) > 1:
        for pid, r in frames.items():
            others = pd.concat([x for q, x in frames.items() if q != pid])
            tpo, tao, *_ = fit(others)
            rows.append({"Fitted_on": f"plates except {pid}", "Tested_on": f"plate {pid}",
                         "Present_threshold_cps": tpo, "Absent_threshold_cps": tao, **apply(r, tpo, tao)})
    table = pd.DataFrame(rows)
    for c in ("false_present", "products_missed"):
        table[f"{c}_%"] = (100 * table[c] / table["wells"]).round(2)
    current = apply(all_r, s.level5_cps, s.level2_raw_max_cps)   # the present settings: level 5 / levels 1-2
    table["Note"] = ""
    if np.isfinite(tp) and tp < s.cutoff_cps:
        table.loc[0, "Note"] = (f"the present threshold is below the cutoff ({s.cutoff_cps:,.0f}); level 5 cannot be "
                                "lower than the cutoff, so the suggestion for level 5 is the cutoff itself")
    return {"table": table, "present_curve": pc, "absent_curve": ac, "target_pct": target,
            "suggested_level5_cps": max(tp, s.cutoff_cps) if np.isfinite(tp) else tp, "suggested_absent_raw_cps": ta, "current": current,
            "auc_evidence": auc(all_r["Evidence_Intensity"].fillna(0), all_r["HPLC_Product"]),
            "auc_raw": auc(all_r["Raw_Ejection_Intensity"].fillna(0), all_r["HPLC_Product"])}


def auc(score, y) -> float:
    """Area under the ROC curve (probability that a well with product scores higher than one without)."""
    s, y = np.asarray(score, float), np.asarray(y)
    pos, neg = s[y == 1], s[y == 0]
    if not len(pos) or not len(neg):
        return float("nan")
    ranks = pd.Series(np.concatenate([pos, neg])).rank().to_numpy()
    return round(float((ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))), 3)


# ------------------------------------------------------------------------------------------
# Probability of product per well
# ------------------------------------------------------------------------------------------
FEATURES = ["log10_evidence", "log10_raw", "echo_ratio", "echo_ratio_missing", "log10_prev_ratio",
            "sm_interference", "isotope_score", "isotope_missing", "valid_spectra"]


def features(r: pd.DataFrame) -> pd.DataFrame:
    ev = r["Evidence_Intensity"].fillna(0).clip(lower=0)
    raw = r["Raw_Ejection_Intensity"].fillna(0).clip(lower=0)
    prev = r["Previous_Ejection_Intensity"].fillna(0).clip(lower=0) if "Previous_Ejection_Intensity" in r else 0 * ev
    ratio = r["Echo_Product_Ratio_%"] if "Echo_Product_Ratio_%" in r else pd.Series(np.nan, index=r.index)
    iso = r["Isotope_Pattern_Score"] if "Isotope_Pattern_Score" in r else pd.Series(np.nan, index=r.index)
    vs = [c for c in r.columns if c.endswith("_Valid_Spectra")]
    valid = r[vs].fillna(0).max(axis=1) if vs else 0 * ev
    return pd.DataFrame({
        "log10_evidence": np.log10(ev + 1e3), "log10_raw": np.log10(raw + 1e3),
        "echo_ratio": ratio.fillna(0) / 100, "echo_ratio_missing": ratio.isna().astype(float),
        "log10_prev_ratio": np.log10((prev + 1e3) / (raw + 1e3)),
        "sm_interference": r["SM_Interference"].fillna("").astype(str).str.len().gt(0).astype(float)
        if "SM_Interference" in r else 0.0,
        "isotope_score": iso.fillna(0), "isotope_missing": iso.isna().astype(float),
        "valid_spectra": np.minimum(valid, 10) / 10}, index=r.index)


def _fit_logistic(X, y, l2=1.0, iters=50):
    Xb = np.column_stack([np.ones(len(X)), X])
    w = np.zeros(Xb.shape[1])
    reg = l2 * np.eye(len(w)); reg[0, 0] = 0
    for _ in range(iters):
        p = 1 / (1 + np.exp(-np.clip(Xb @ w, -30, 30)))
        g = Xb.T @ (p - y) + reg @ w
        H = Xb.T @ (Xb * (p * (1 - p))[:, None]) + reg
        step = np.linalg.solve(H, g)
        w -= step
        if np.abs(step).max() < 1e-8:
            break
    return w


class ProbabilityModel:
    def __init__(self, mean, std, weights, info=None):
        self.mean, self.std, self.weights = np.asarray(mean), np.asarray(std), np.asarray(weights)
        self.info = info or {}

    @classmethod
    def fit(cls, r: pd.DataFrame, l2: float = 1.0, info=None):
        X = features(r)[FEATURES].to_numpy(float)
        mean, std = X.mean(axis=0), X.std(axis=0)
        std[std == 0] = 1.0
        w = _fit_logistic((X - mean) / std, r["HPLC_Product"].to_numpy(float), l2)
        return cls(mean, std, w, info)

    def predict(self, r: pd.DataFrame) -> np.ndarray:
        X = (features(r)[FEATURES].to_numpy(float) - self.mean) / self.std
        z = self.weights[0] + X @ self.weights[1:]
        p = 1 / (1 + np.exp(-np.clip(z, -30, 30)))
        p[(r["Confidence_Level"] == 0).to_numpy()] = np.nan
        return p

    def save(self, path) -> Path:
        d = {"features": FEATURES, "mean": self.mean.tolist(), "std": self.std.tolist(),
             "weights": self.weights.tolist(), "info": self.info}
        Path(path).write_text(json.dumps(d, indent=2), encoding="utf-8")
        return Path(path)

    @classmethod
    def load(cls, path):
        d = json.loads(Path(path).read_text(encoding="utf-8"))
        if d.get("features") != FEATURES:
            raise ValueError("this probability model was made by another version of the program (different features)")
        return cls(d["mean"], d["std"], d["weights"], d.get("info"))


def brier(p, y) -> float:
    p, y = np.asarray(p, float), np.asarray(y, float)
    m = np.isfinite(p)
    return round(float(np.mean((p[m] - y[m]) ** 2)), 4) if m.any() else float("nan")


def probability(results: dict, model: ProbabilityModel | None = None) -> dict:
    """Fits the model on the plates with HPLC-MS (leave-one-plate-out check when there are >= 2) and adds
    P_Product to every plate; with `model`, applies a saved model instead of fitting."""
    out = {"validation": pd.DataFrame(), "model": model}
    train = {pid: r[r["Confidence_Level"] > 0] for pid, r in _samples(results).items() if "HPLC_Product" in r.columns}
    if model is None:
        if not train:
            return out
        all_t = pd.concat(train.values())
        if all_t["HPLC_Product"].nunique() < 2 or len(all_t) < 40:
            out["note"] = "too few wells or only one class in HPLC-MS: no probability model"
            return out
        rows = []
        if len(train) > 1:
            for pid, r in train.items():
                others = pd.concat([x for q, x in train.items() if q != pid])
                if others["HPLC_Product"].nunique() < 2:
                    continue
                p = ProbabilityModel.fit(others).predict(r)
                rows.append({"Fitted_on": f"plates except {pid}", "Tested_on": f"plate {pid}", "Wells": len(r),
                             "AUC": auc(np.nan_to_num(p), r["HPLC_Product"]), "Brier": brier(p, r["HPLC_Product"])})
        model = ProbabilityModel.fit(all_t, info={"plates": sorted(train), "wells": len(all_t)})
        p_all = model.predict(all_t)
        rows.append({"Fitted_on": "all plates", "Tested_on": "all plates (same data)", "Wells": len(all_t),
                     "AUC": auc(np.nan_to_num(p_all), all_t["HPLC_Product"]), "Brier": brier(p_all, all_t["HPLC_Product"])})
        out["validation"], out["model"] = pd.DataFrame(rows), model
    for r in results.values():
        r["P_Product"] = np.round(model.predict(r), 3)
    return out


# ------------------------------------------------------------------------------------------
# Run QC
# ------------------------------------------------------------------------------------------
def run_qc(pa, s) -> pd.DataFrame:
    """Per well: TIC pulse height of its ejection (relative to the plate median) and, when given, the
    internal-standard signal. Wells in a stretch where the spray dropped are flagged (< 30% of the
    running median of the pulse heights)."""
    from .runs import ion_window_max, window_scans
    run, grid = pa.run, pa.grid
    rows = []
    is_x = None
    if s.internal_standard_mz and run.mz_min <= s.internal_standard_mz <= run.mz_max:
        is_x = ion_window_max(run, float(s.internal_standard_mz), s.mass_tolerance_da)
    for w, sl in sorted(pa.slot_of.items(), key=lambda kv: kv[1]):
        if sl >= len(grid):
            continue
        a, b = window_scans(run.t, grid[sl], -s.window_pre_s, s.window_post_s)
        rows.append({"Well": w, "Ejection_Time_s": round(float(grid[sl]), 2),
                     "Pulse_Height": float(pa.rel_tic[a:b].max()) if b > a else np.nan,
                     "IS_Intensity": float(is_x[a:b].max()) if is_x is not None and b > a else np.nan})
    q = pd.DataFrame(rows)
    if q.empty:
        return q
    med = q["Pulse_Height"].median()
    q["Pulse_vs_Median"] = (q["Pulse_Height"] / med).round(2) if med > 0 else np.nan
    run_med = q["Pulse_Height"].rolling(25, center=True, min_periods=5).median()
    q["Run_QC_Flag"] = np.where(q["Pulse_Height"] < 0.3 * run_med, "weak ejection pulse", "")
    if is_x is not None:
        im = q["IS_Intensity"].median()
        q["IS_vs_Median"] = (q["IS_Intensity"] / im).round(2) if im > 0 else np.nan
        low = q["IS_Intensity"] < 0.5 * im
        q.loc[low, "Run_QC_Flag"] = (q.loc[low, "Run_QC_Flag"] + "; internal standard < 50% of median").str.strip("; ")
    return q


def drift(q: pd.DataFrame, col: str = "Pulse_Height") -> dict:
    """Spearman correlation of a QC signal with ejection time (a steady fall = drift)."""
    if q.empty or q[col].notna().sum() < 5:
        return {}
    rho = q[col].corr(q["Ejection_Time_s"], method="spearman")
    first, last = q[col].iloc[: len(q) // 4].median(), q[col].iloc[-len(q) // 4:].median()
    return {"spearman_vs_time": round(float(rho), 3),
            "last_quarter_vs_first_%": round(100 * last / first, 1) if first else float("nan")}


# ------------------------------------------------------------------------------------------
# Comparing two analyses
# ------------------------------------------------------------------------------------------
def mcnemar(a, b) -> dict:
    """Exact McNemar test on paired binary outcomes a, b (e.g. detected with setting A and with B)."""
    from scipy.stats import binomtest
    a, b = np.asarray(a).astype(bool), np.asarray(b).astype(bool)
    n01, n10 = int((~a & b).sum()), int((a & ~b).sum())
    p = binomtest(min(n01, n10), n01 + n10, 0.5).pvalue if n01 + n10 else 1.0
    return {"only_A": n10, "only_B": n01, "p_value": round(float(p), 4)}


def compare(res_a: dict, res_b: dict, label_a="A", label_b="B") -> pd.DataFrame:
    """Two analyses of the same wells: changed calls, and (with HPLC-MS) whether the accuracy differs."""
    rows = []
    for pid in sorted(set(res_a) & set(res_b)):
        x = res_a[pid].merge(res_b[pid], on="Well", suffixes=("_A", "_B"))
        row = {"Plate": pid, "Wells": len(x),
               "Level_changed": int((x["Confidence_Level_A"] != x["Confidence_Level_B"]).sum())}
        m = mcnemar(x["Detected_A"], x["Detected_B"])
        row.update({f"Detected_only_{label_a}": m["only_A"], f"Detected_only_{label_b}": m["only_B"],
                    "McNemar_p_detection": m["p_value"]})
        if "HPLC_Product_A" in x.columns:
            y = x["HPLC_Product_A"]
            ca, cb = x["Detected_A"] == y, x["Detected_B"] == y
            m2 = mcnemar(ca, cb)
            row.update({f"Correct_{label_a}_%": round(100 * ca.mean(), 1), f"Correct_{label_b}_%": round(100 * cb.mean(), 1),
                        f"Right_only_{label_a}": m2["only_A"], f"Right_only_{label_b}": m2["only_B"],
                        "McNemar_p_accuracy": m2["p_value"]})
        rows.append(row)
    return pd.DataFrame(rows)
