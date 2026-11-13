"""Per-well detection: ejection order, ion search with the isotope check, background control,
confidence levels."""
from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .adducts import IONS, PRIMARY, SM_INTERFERENCE_IONS, SM_MAIN_ION, Ion, ions_for
from .plates import PlateFormat
from . import isotopes
from .runs import Run, fit_ejection_grid, ion_window_max, period_from_name, window_scans
from .settings import LEVEL_LABELS, Settings

log = logging.getLogger("triade")


@dataclass
class PlateAnalysis:
    plate_id: str
    run: Run
    results: pd.DataFrame
    period: float
    grid: np.ndarray
    rel_tic: np.ndarray
    order_name: str
    first_slot: int
    slot_of: dict
    order_scores: pd.DataFrame
    blank: tuple
    warnings: list = field(default_factory=list)
    replicates: dict = field(default_factory=dict)   # replicate number -> PlateAnalysis (replicate ejections)


# ------------------------------------------------------------------------------------------
# Ejection order
# ------------------------------------------------------------------------------------------
def blank_window(P: float, s: Settings, warnings: list):
    """Blank window (start, end) relative to the apex: ends where the next well's window starts and has
    the width of the ejection window, without starting before blank_start_s."""
    end = P - s.window_pre_s
    width = s.window_pre_s + s.window_post_s
    start = max(s.blank_start_s, end - width)
    if end - start < 0.5 * width:
        msg = (f"A {P:.2f} s period leaves only {max(end - start, 0):.2f} s of blank window between ejections: "
               "the false-positive estimate from the blank window is weak.")
        log.warning("  " + msg); warnings.append(msg)
    return start, max(end, start)


def _primary_mz(df: pd.DataFrame, pol: str) -> np.ndarray:
    ion = IONS[PRIMARY[pol]]
    return np.array([ion.mz(M) if np.isfinite(M) else np.nan for M in df["Monoisotopic_Mass"]])


def detect_order(run: Run, df: pd.DataFrame, grid: np.ndarray, plate: PlateFormat, s: Settings, warnings: list,
                 first_slot: int | None = None):
    """Finds which pulse belongs to which well. Uses the injection order of the sample table when given;
    otherwise tests the usual ejection patterns. For each candidate (order, first pulse) it counts the
    wells whose protonated/deprotonated product is >= align_threshold in >= 2 spectra of their ejection."""
    n_slots = len(grid)
    bounds = np.array([window_scans(run.t, g, -s.window_pre_s, s.window_post_s) for g in grid])
    wells = list(df["Well"])
    hit = np.zeros((len(wells), n_slots), dtype=bool)
    for pol in ("+", "-"):
        is_pol = run.pol == pol
        if not is_pol.any():
            continue
        for k, mz in enumerate(_primary_mz(df, pol)):
            if not np.isfinite(mz):
                continue
            above = (ion_window_max(run, mz, s.mass_tolerance_da) >= s.align_threshold_cps) & is_pol
            cs = np.concatenate([[0], np.cumsum(above)])
            hit[k] |= (cs[bounds[:, 1]] - cs[bounds[:, 0]]) >= 2

    has_table_order = df["Injection_Order"].notna().any()
    positions = {}
    if s.ejection_order == "table" or (s.ejection_order == "auto" and has_table_order):
        if not has_table_order:
            raise ValueError("ejection_order is 'table' but the sample table has no Injection_Order for this plate")
        sub = df[df["Injection_Order"].notna()].sort_values("Injection_Order")
        base = df["_Order_Base"].iloc[0] if "_Order_Base" in df.columns else sub["Injection_Order"].min()
        positions = {w: int(o - base) for w, o in zip(sub["Well"], sub["Injection_Order"])}
        orders = {"injection order from the sample table": list(sub["Well"])}
    else:
        # the Echo skips wells that are not in the survey: keep only the wells of the sample table
        in_table = set(wells)
        orders = {n: [w for w in o if w in in_table] for n, o in plate.ejection_orders().items()}
        if s.ejection_order != "auto":
            pick = [n for n in orders if n.startswith(s.ejection_order)]
            if not pick:
                raise ValueError(f"unknown ejection order {s.ejection_order!r}; options: auto, table, "
                                 + ", ".join(n.split(" (")[0] for n in orders))
            orders = {pick[0]: orders[pick[0]]}
    widx = {w: k for k, w in enumerate(wells)}
    rows = []
    for name, order in orders.items():
        pos = positions if name.startswith("injection order") else {w: i for i, w in enumerate(order)}
        span = max(pos.values()) + 1 if pos else 0
        if span > n_slots:
            msg = (f"The run has {n_slots} ejection pulses but the order '{name}' needs {span}: "
                   "check the plate format, the ejection interval and the injection order.")
            warnings.append(msg); log.warning("  " + msg)
            continue
        ok = [w for w in order if w in widx]
        ki = np.array([widx[w] for w in ok], dtype=int)
        rel_slots = np.array([pos[w] for w in ok], dtype=int)
        for first in range(0, n_slots - span + 1):
            rows.append({"order": name, "first_slot": first, "first_ejection_s": round(float(grid[first]), 2),
                         "wells_explained": int(hit[ki, first + rel_slots].sum()) if len(ki) else 0})
    if not rows:
        raise ValueError("No ejection order fits the number of pulses found in the run.")
    scores = pd.DataFrame(rows).sort_values("wells_explained", ascending=False).reset_index(drop=True)
    best = scores.iloc[0]
    order_name, first = best["order"], int(best["first_slot"])
    if s.first_slot is not None:
        first = int(s.first_slot)
    if first_slot is not None:          # replicate ejections: the first pulse found for replicate 1
        first = int(first_slot)
    second = scores.iloc[1]["wells_explained"] if len(scores) > 1 else 0
    log.info("  Ejection order: %s; first well at pulse %d (t = %.2f s). Explains %d wells (second-best: %d).",
             order_name, first, grid[first], best["wells_explained"], second)
    if best["wells_explained"] == 0:
        msg = ("No ejection-order hypothesis explains the data: the well assignment is arbitrary. "
               "Give the injection order in the sample table.")
        warnings.append(msg); log.warning("  " + msg)
    elif best["wells_explained"] < 1.1 * second:
        msg = "Small margin between ejection-order hypotheses; check the assignment against the sample list."
        warnings.append(msg); log.warning("  " + msg)
    order = orders[order_name]
    pos = positions if order_name.startswith("injection order") else {w: i for i, w in enumerate(order)}
    slot_of = {w: first + pos[w] for w in order}
    return order_name, first, slot_of, scores


def first_slot_for_table(run: Run, df: pd.DataFrame, s: Settings) -> int:
    """First ejection pulse for a table order that repeats wells (replicate ejections): every row is scored
    at its own position (Injection_Order - first order), so the replicates are placed together."""
    hint = s.ejection_period_s or period_from_name(run.path.name)
    _, grid, _ = fit_ejection_grid(run, hint)
    bounds = np.array([window_scans(run.t, g, -s.window_pre_s, s.window_post_s) for g in grid])
    rows = df[df["Injection_Order"].notna()]
    pos = (rows["Injection_Order"] - rows["Injection_Order"].min()).astype(int).to_numpy()
    hit = np.zeros((len(rows), len(grid)), dtype=bool)
    for pol in ("+", "-"):
        is_pol = run.pol == pol
        if not is_pol.any():
            continue
        for k, mz in enumerate(_primary_mz(rows, pol)):
            if not np.isfinite(mz):
                continue
            above = (ion_window_max(run, mz, s.mass_tolerance_da) >= s.align_threshold_cps) & is_pol
            cs = np.concatenate([[0], np.cumsum(above)])
            hit[k] |= (cs[bounds[:, 1]] - cs[bounds[:, 0]]) >= 2
    span = int(pos.max()) + 1
    scores = [int(hit[np.arange(len(rows)), first + pos].sum()) for first in range(0, len(grid) - span + 1)]
    if not scores:
        raise ValueError(f"The run has {len(grid)} ejection pulses but the injection order needs {span}.")
    return int(np.argmax(scores)) if s.first_slot is None else int(s.first_slot)


# ------------------------------------------------------------------------------------------
# Ion search in one spectrum
# ------------------------------------------------------------------------------------------
def _local_maxima(x: np.ndarray, lo: int, hi: int) -> np.ndarray:
    idx = np.arange(max(lo, 0), min(hi, len(x) - 1) + 1)
    if idx.size == 0:
        return idx
    left = x[np.maximum(idx - 1, 0)]
    right = x[np.minimum(idx + 1, len(x) - 1)]
    v = x[idx]
    return idx[(v > 0) & (v >= left) & (v >= right)]


def _best_peak(x: np.ndarray, lo: int, hi: int, target_bin: float):
    cand = _local_maxima(x, lo, hi)
    if cand.size == 0:
        return None
    v = x[cand]
    top = cand[v == v.max()]
    return int(top[np.argmin(np.abs(top - target_bin))])


def ion_peaks(x: np.ndarray, bin0: int, bin_da: float, target_mz: float, spacing: float, s: Settings):
    """M peak (the highest within ±tolerance of the theoretical m/z) and the highest M+1 that is a distinct
    peak (a valley between them) at M + spacing ± tolerance. Returns bin indices (m, c); None when absent."""
    tol = s.mass_tolerance_da
    shift_tol = s.isotope_shift_tolerance_da * spacing
    tb = target_mz / bin_da - bin0
    lo = math.ceil(tb - tol / bin_da - 1e-6)
    hi = math.floor(tb + tol / bin_da + 1e-6)
    m = _best_peak(x, lo, hi, tb)
    if m is None:
        return None, None
    t1 = tb + spacing / bin_da
    lo1 = max(math.ceil(t1 - tol / bin_da - 1e-6), math.ceil(m + (spacing - shift_tol) / bin_da - 1e-6))
    hi1 = min(math.floor(t1 + tol / bin_da + 1e-6), math.floor(m + (spacing + shift_tol) / bin_da + 1e-6))
    cand = _local_maxima(x, lo1, hi1)
    cand = [c for c in cand[cand > m] if x[m + 1:c].min() < x[c]]
    if not cand:
        return m, None
    return m, max(cand, key=lambda b: (x[b], -abs(b - (m + spacing / bin_da))))


def evaluate_scan(x, run: Run, target_mz: float, spacing: float, s: Settings, m1_rel: float | None = None):
    """State 0 = no M peak >= cutoff; 1 = M ok but no valid M+1; 2 = validated.
    score = the highest cutoff at which this spectrum would be validated (0 without M+1).
    M+1 rule 'absolute' (default): M+1 >= m1_fraction_of_cutoff x cutoff. Rule 'relative': M+1 >=
    m1_relative_fraction x the M+1 expected from the formula (m1_rel x M) and >= m1_floor_cps."""
    m, c = ion_peaks(x, run.bin0, run.bin_da, target_mz, spacing, s)
    if m is None:
        return 0, None, 0.0
    if s.m1_rule == "relative" and m1_rel:
        m1_ok = c is not None and x[c] >= max(s.m1_relative_fraction * m1_rel * x[m], s.m1_floor_cps)
        score = float(x[m]) if m1_ok else 0.0
    else:
        frac = s.m1_fraction_of_cutoff
        m1_need = x[c] / frac if c is not None and frac > 0 else x[m]
        score = float(min(x[m], m1_need)) if c is not None else 0.0
        m1_ok = c is not None and x[c] >= s.cutoff_cps * frac
    if x[m] < s.cutoff_cps:
        return 0, None, score
    data = {"I_M": float(x[m]), "mz_M": (m + run.bin0) * run.bin_da, "bin_M": int(m)}
    if not m1_ok:
        return 1, data, score
    data.update({"I_M1": float(x[c]), "mz_M1": (c + run.bin0) * run.bin_da, "shift": (c - m) * run.bin_da})
    return 2, data, score


C13 = 1.0033548


def evaluate_scan_hr(run: Run, k: int, target_mz: float, spacing: float, s: Settings, m1_rel: float | None = None):
    """High-resolution version of evaluate_scan, on the centroid peaks of scan k: M = the most intense peak
    within ± ppm of the theoretical m/z; M+1 = the most intense peak at M + 1.00335/|z| ± max(ppm, 3 mDa)
    (13C/15N fine structure). Same states, scores and M+1 rules as in the unit-resolution mode."""
    mz, it = run.peaks(k)
    tol = run.tol(target_mz, s.mass_tolerance_da)
    i0, i1 = np.searchsorted(mz, target_mz - tol), np.searchsorted(mz, target_mz + tol, "right")
    if i1 <= i0:
        return 0, None, 0.0
    j = i0 + int(np.argmax(it[i0:i1]))
    I_M, mz_M = float(it[j]), float(mz[j])
    t1 = mz_M + C13 * spacing
    tol1 = max(run.tol(t1, s.mass_tolerance_da), 0.003)
    j0, j1 = np.searchsorted(mz, t1 - tol1), np.searchsorted(mz, t1 + tol1, "right")
    I_1 = float(it[j0:j1].max()) if j1 > j0 else None
    mz_1 = float(mz[j0 + int(np.argmax(it[j0:j1]))]) if j1 > j0 else None
    if s.m1_rule == "relative" and m1_rel:
        m1_ok = I_1 is not None and I_1 >= max(s.m1_relative_fraction * m1_rel * I_M, s.m1_floor_cps)
        score = I_M if m1_ok else 0.0
    else:
        frac = s.m1_fraction_of_cutoff
        score = float(min(I_M, I_1 / frac if frac > 0 else I_M)) if I_1 is not None else 0.0
        m1_ok = I_1 is not None and I_1 >= s.cutoff_cps * frac
    if I_M < s.cutoff_cps:
        return 0, None, score
    data = {"I_M": I_M, "mz_M": mz_M}
    if not m1_ok:
        return 1, data, score
    data.update({"I_M1": I_1, "mz_M1": mz_1, "shift": mz_1 - mz_M})
    return 2, data, score


def ion_status(run: Run, mz: float, ion: Ion, s: Settings) -> str:
    if not np.isfinite(mz):
        return "NO STRUCTURE"
    if not np.any(run.pol == ion.polarity):
        return "NOT ACQUIRED (no " + ("negative" if ion.polarity == "-" else "positive") + "-mode scans)"
    if mz - s.mass_tolerance_da < run.mz_min - 1e-6 or mz + ion.isotope_spacing + s.mass_tolerance_da > run.mz_max + 1e-6:
        return f"OUT OF SCAN RANGE ({run.mz_min:.0f}-{run.mz_max:.0f})"
    return ""


def analyse_window(run: Run, a: int, b: int, mz: float, ion: Ion, s: Settings, m1_rel: float | None = None):
    """(validated spectra, spectra with M >= cutoff, best spectrum, evidence), where evidence = the highest
    cutoff at which >= min_spectra spectra would still be validated."""
    n_valid = n_m = 0
    best = None
    scores = []
    for k in range(a, b):
        if run.pol[k] != ion.polarity:
            continue
        if run.sparse:
            st, d, score = evaluate_scan_hr(run, k, mz, ion.isotope_spacing, s, m1_rel)
        else:
            st, d, score = evaluate_scan(run.dense[k], run, mz, ion.isotope_spacing, s, m1_rel)
        scores.append(score)
        n_m += st >= 1
        if st == 2:
            n_valid += 1
            if best is None or d["I_M"] > best["I_M"]:
                best = dict(d, scan=k)
    k = max(s.min_spectra, 1)
    evidence = sorted(scores)[-k] if len(scores) >= k else 0.0
    return n_valid, n_m, best, evidence


def isotope_info(run: Run, a: int, b: int, ion, pat, mz: float, s: Settings) -> dict:
    """Observed vs theoretical isotope pattern of one ion in the ejection window: the median ratio of each
    isotope peak to M over the spectra whose M peak is >= level3_cps."""
    out = {"Isotope_Ion": ion.name if ion else "", "Isotope_Theoretical": "", "Isotope_Observed": "",
           "Isotope_Confirmation_Peak": "", "Isotope_Pattern_Score": np.nan, "Isotope_Spectra_Used": 0}
    if not pat or not np.isfinite(mz) or b <= a:
        return out
    k_conf, rel_conf, _ = isotopes.confirmation_peak(pat)
    out["Isotope_Theoretical"] = isotopes.describe(pat)
    out["Isotope_Confirmation_Peak"] = f"M+{k_conf} ({100 * rel_conf:.0f}% of M)"
    tol = s.mass_tolerance_da
    ratios = []
    if run.sparse:
        for k in range(a, b):
            if run.pol[k] != ion.polarity:
                continue
            pm, pi = run.peaks(k)
            tl = run.tol(mz, tol)
            i0, i1 = np.searchsorted(pm, mz - tl), np.searchsorted(pm, mz + tl, "right")
            if i1 <= i0:
                continue
            j = i0 + int(np.argmax(pi[i0:i1]))
            if pi[j] < s.level3_cps:
                continue
            r = {}
            for kk, _, shift in pat:
                if kk == 0:
                    continue
                c = pm[j] + shift
                w = max(run.tol(c, tol), 0.003)
                j0, j1 = np.searchsorted(pm, c - w), np.searchsorted(pm, c + w, "right")
                r[kk] = float(pi[j0:j1].max()) / float(pi[j]) if j1 > j0 else 0.0
            ratios.append(r)
    else:
        tb = mz / run.bin_da - run.bin0
        lo, hi = math.ceil(tb - tol / run.bin_da - 1e-6), math.floor(tb + tol / run.bin_da + 1e-6)
        for k in range(a, b):
            if run.pol[k] != ion.polarity:
                continue
            x = run.dense[k]
            m = _best_peak(x, lo, hi, tb)
            if m is not None and x[m] >= s.level3_cps:
                ratios.append(isotopes.observed_ratios(x, m, pat, run.bin_da))
    if ratios:
        med = {kk: float(np.median([r.get(kk, 0.0) for r in ratios])) for kk in ratios[0]}
        out["Isotope_Observed"] = isotopes.describe_obs(med, pat)
        out["Isotope_Pattern_Score"] = round(isotopes.pattern_score(med, pat), 2)
        out["Isotope_Spectra_Used"] = len(ratios)
    return out


def confidence_level(evidence: float, raw: float, s: Settings) -> int:
    if evidence >= s.level5_cps:
        return 5
    if evidence >= s.level4_cps:
        return 4
    if evidence >= s.level3_cps or raw > s.level2_raw_max_cps:
        return 3
    return 2 if raw > s.level1_raw_max_cps else 1


# ------------------------------------------------------------------------------------------
# Background groups
# ------------------------------------------------------------------------------------------
def background_groups(df: pd.DataFrame, plate: PlateFormat, s: Settings) -> dict:
    """well -> label of its background group (wells sharing the starting material, row or column)."""
    if s.background_group == "off":
        return {}
    out = {}
    keys = df["Background_Key"] if "Background_Key" in df.columns else df["SM_SMILES"]
    for w, sm in zip(df["Well"], keys):
        r, c = plate.split(w)
        if s.background_group == "starting_material":
            # the shared reagent (Reagent_1 as a rule); without one, fall back to the plate row
            out[w] = f"SM:{sm}" if sm else f"row:{r}"
        elif s.background_group == "row":
            out[w] = f"row:{r}"
        else:
            out[w] = f"col:{c}"
    return out


# ------------------------------------------------------------------------------------------
# Plate
# ------------------------------------------------------------------------------------------
def analyse_plate(run: Run, df: pd.DataFrame, plate: PlateFormat, plate_id: str, s: Settings,
                  first_slot: int | None = None) -> PlateAnalysis:
    warnings: list = []
    hint = s.ejection_period_s or period_from_name(run.path.name)
    P, grid, rel = fit_ejection_grid(run, hint)
    log.info("  Ejection period: %.4f s (%d possible pulses in the run)", P, len(grid))
    order_name, first, slot_of, scores = detect_order(run, df, grid, plate, s, warnings, first_slot)
    blank = blank_window(P, s, warnings)
    rank = {w: sl - first + 1 for w, sl in slot_of.items()}     # position in the ejection sequence
    acquired = run.polarities
    sel = ions_for(s.ions, s.polarity)
    primaries = [IONS[PRIMARY[p]] for p in ("+", "-") if p in acquired]
    groups = background_groups(df, plate, s)
    ion_mz = {ion.name: dict(zip(df["Well"], df[f"{ion.key}_Theoretical_mz"])) for ion in sel}
    masses = dict(zip(df["Well"], df["Monoisotopic_Mass"]))
    win = {w: window_scans(run.t, grid[sl], -s.window_pre_s, s.window_post_s)
           for w, sl in slot_of.items() if sl < len(grid)}
    tol = s.mass_tolerance_da

    def group_background(well, target, pol):
        g = groups.get(well)
        if g is None:
            return np.nan
        x = ion_window_max(run, target, tol)
        vals = []
        for v, gv in groups.items():
            if v == well or gv != g or v not in win:
                continue
            near = any(abs(ion_mz[i.name][v] - target) <= s.background_min_distance_da
                       for i in sel if i.polarity == pol and np.isfinite(ion_mz[i.name][v]))
            if near:
                continue
            a, b = win[v]
            if b > a:
                vals.append(float(x[a:b].max()))
        return float(np.median(vals)) if vals else np.nan

    out = []
    for r in df.itertuples(index=False):
        w = r.Well
        row_label, col = plate.split(w)
        slot = slot_of.get(w)
        apex = grid[slot] if slot is not None and slot < len(grid) else np.nan
        row = {"Plate": plate_id, "Well": w, "Row": row_label, "Column": col,
               "Ejection_Order": rank.get(w), "Ejection_Time_s": round(float(apex), 2) if np.isfinite(apex) else np.nan,
               "Replicate": getattr(r, "Replicate", 1), "Control": getattr(r, "Control", ""),
               "Sample_Code": r.Sample_Code, "Product_SMILES": r.Product_SMILES, "SMILES_Note": r.SMILES_Note,
               "Formula": r.Formula, "Monoisotopic_Mass": round(r.Monoisotopic_Mass, 4) if np.isfinite(r.Monoisotopic_Mass) else np.nan}
        a, b = window_scans(run.t, apex, -s.window_pre_s, s.window_post_s) if np.isfinite(apex) else (0, 0)
        ba, bb = window_scans(run.t, apex, *blank) if np.isfinite(apex) else (0, 0)
        row["Scans_In_Window"] = int(b - a)
        M = r.Monoisotopic_Mass

        # starting-material interferences on the selected product ions (acquired polarities only)
        hits, primary_hit = [], False
        for smM in r.SM_Masses:
            for pol in acquired:
                for name, d in SM_INTERFERENCE_IONS[pol].items():
                    for ion in sel:
                        if ion.polarity != pol or not np.isfinite(M):
                            continue
                        pmz = ion.mz(M)
                        sm_tol = s.sm_interference_da if not run.sparse else 2 * run.tol(pmz, 0.0)
                        if abs(smM + d - pmz) <= sm_tol:
                            hits.append(f"{name} {smM + d:.2f} is {smM + d - pmz:+.2f} Da from product {ion.name}")
                            primary_hit |= ion.name in PRIMARY.values()
        row["SM_SMILES"], row["SM_Formula"] = r.SM_SMILES, r.SM_Formula
        row["SM_Interference"] = "; ".join(hits)

        detected_ions, blank_ions, reasons, best_I = [], [], [], 0.0
        pats = {}
        evidence_I = evidence_primary = 0.0
        searchable = False
        for ion in sel:
            k = ion.key
            mz = ion_mz[ion.name][w]
            row[f"{k}_Theoretical_mz"] = round(mz, 4) if np.isfinite(mz) else np.nan
            status = ion_status(run, mz, ion, s)
            if status or b <= a:
                row.update({f"{k}_Status": status or "NO EJECTION WINDOW", f"{k}_Valid_Spectra": np.nan,
                            f"{k}_Spectra_M_Above_Cutoff": np.nan, f"{k}_Max_Intensity": np.nan,
                            f"{k}_Observed_mz": np.nan, f"{k}_M1_Observed_mz": np.nan, f"{k}_M1_Intensity": np.nan,
                            f"{k}_Observed_Shift": np.nan, f"{k}_Raw_Ejection_Intensity": np.nan,
                            f"{k}_Ejection_vs_Blank_Ratio": np.nan, f"{k}_Background": np.nan,
                            f"{k}_Ejection_vs_Background": np.nan})
                continue
            pat = isotopes.pattern(r.Formula, ion) if r.Formula else []
            m1_rel = next((rel for k_, rel, _ in pat if k_ == 1), None)
            pats[ion.name] = pat
            n_valid, n_m, best, evidence = analyse_window(run, a, b, mz, ion, s, m1_rel)
            xmax = ion_window_max(run, mz, tol)
            ej_raw = float(xmax[a:b].max())
            bl_raw = float(xmax[ba:bb].max()) if bb > ba else 0.0
            bg = group_background(w, mz, ion.polarity) if s.background_group != "off" else np.nan
            bg_ratio = round(ej_raw / max(bg if np.isfinite(bg) else 0.0, s.background_floor_cps), 2)
            bg_ok = s.background_group == "off" or s.background_ratio <= 0 or bg_ratio >= s.background_ratio
            det = n_valid >= s.min_spectra and bg_ok
            row.update({
                f"{k}_Status": "DETECTED" if det else "NOT DETECTED",
                f"{k}_Valid_Spectra": n_valid, f"{k}_Spectra_M_Above_Cutoff": n_m,
                f"{k}_Max_Intensity": best["I_M"] if best else 0.0,
                f"{k}_Observed_mz": round(best["mz_M"], 1) if best else np.nan,
                f"{k}_M1_Observed_mz": round(best["mz_M1"], 1) if best else np.nan,
                f"{k}_M1_Intensity": best["I_M1"] if best else np.nan,
                f"{k}_Observed_Shift": round(best["shift"], 2) if best else np.nan,
                f"{k}_Raw_Ejection_Intensity": ej_raw,
                f"{k}_Ejection_vs_Blank_Ratio": round(ej_raw / max(bl_raw, s.background_floor_cps), 2),
                f"{k}_Background": bg, f"{k}_Ejection_vs_Background": bg_ratio,
            })
            searchable = True
            if best:
                best_I = max(best_I, best["I_M"])
            if bg_ok:
                evidence_I = max(evidence_I, evidence)
                if ion.name in PRIMARY.values():
                    evidence_primary = max(evidence_primary, evidence)
            if det:
                detected_ions.append(ion.name)
            elif n_valid >= s.min_spectra:
                reasons.append(f"{ion.name}: valid in {n_valid} spectra but only {bg_ratio}x the background of wells "
                               f"with the same {s.background_group.replace('_', ' ')} (< {s.background_ratio:g}x)")
            elif n_valid == 1:
                reasons.append(f"{ion.name}: valid in 1 spectrum only")
            elif n_m > 0:
                reasons.append(f"{ion.name}: M peak without a valid M+1")
            nv_blank, *_ = analyse_window(run, ba, bb, mz, ion, s, m1_rel)
            if nv_blank >= s.min_spectra:
                blank_ions.append(ion.name)

        # raw signal of the protonated / deprotonated product (whatever ions were selected)
        raw_I = 0.0
        prim_mz = {}
        for ion in primaries:
            mz = ion.mz(M) if np.isfinite(M) else np.nan
            prim_mz[ion.name] = mz
            if not ion_status(run, mz, ion, s) and b > a:
                raw_I = max(raw_I, float(ion_window_max(run, mz, tol)[a:b].max()))
        detected = bool(detected_ions)
        row["Ions_Detected"] = "; ".join(detected_ions)
        row["Max_Validated_Intensity"] = best_I if best_I > 0 else np.nan
        row["Log10_Max_Validated_Intensity"] = round(math.log10(best_I), 3) if best_I > 0 else np.nan
        row["Detected"] = int(detected)
        if not np.isfinite(M):
            row["Classification"], row["Reason"] = "NOT PROCESSED", r.SMILES_Note
        else:
            row["Classification"] = "DETECTED" if detected else "NOT DETECTED"
            row["Reason"] = "" if detected else ("; ".join(reasons) if reasons else "no ion signal above the cutoff")
        row["Blank_Window_Also_Passes"] = "; ".join(blank_ions)

        # isotope pattern of the main ion (information only; it does not change the level)
        info_ion = next((i for i in sel if i.name in PRIMARY.values() and i.name in pats), None) or \
            next((i for i in sel if i.name in pats), None)
        row.update(isotope_info(run, a, b, info_ion, pats.get(info_ion.name) if info_ion else None,
                                ion_mz[info_ion.name][w] if info_ion else np.nan, s))

        # confidence level and downgrades
        notes = []
        lv = confidence_level(evidence_I, raw_I, s) if searchable and np.isfinite(M) else None
        if lv is None:
            lv = 0
            notes.append("not measurable (no structure, product m/z outside the scanned range, polarity not acquired "
                         "or no ejection window): send to HPLC-MS")
        prev_I = np.nan
        if slot is not None and slot >= 1 and prim_mz:
            pa, pb = window_scans(run.t, grid[slot - 1], -s.window_pre_s, s.window_post_s)
            vals = [float(ion_window_max(run, mz, tol)[pa:pb].max()) for mz in prim_mz.values()
                    if np.isfinite(mz) and pb > pa]
            prev_I = max(vals) if vals else np.nan
        row["Previous_Ejection_Intensity"] = prev_I
        if lv >= 4 and np.isfinite(prev_I) and s.carryover_ratio > 0 and prev_I >= s.carryover_ratio * max(raw_I, 1.0):
            notes.append(f"possible carryover: previous ejection has {prev_I / max(raw_I, 1.0):.1f}x this signal "
                         "at the product m/z")
            lv = 3
        if lv > 3 and primary_hit:
            notes.append("starting-material ion overlaps the protonated/deprotonated product")
            lv = 3
        primary_selected = any(i.name in PRIMARY.values() for i in sel if i.polarity in acquired)
        if lv >= 4 and primary_selected and evidence_primary < s.level4_cps:
            only = [i for i in detected_ions if i not in PRIMARY.values()]
            notes.append(f"only {', '.join(only) or 'secondary ions'} pass the criteria")
            lv = 3
        row["Evidence_Intensity"] = evidence_I if searchable else np.nan
        row["Raw_Ejection_Intensity"] = raw_I if searchable else np.nan
        row["Confidence_Level"] = lv
        row["Confidence_Label"] = LEVEL_LABELS[lv]
        row["Level_Note"] = "; ".join(notes)
        row["Call"] = ("DETECTED" if lv >= 4 else "CHECK SPECTRUM" if lv == 3 else
                       "NOT MEASURABLE" if lv == 0 else "NOT DETECTED")

        # starting-material signal and the product / (product + starting material) ratio
        if r.SM_Masses:
            pol = "+" if "+" in acquired else "-"
            name, d = SM_MAIN_ION[pol]
            sm_raw, ok = 0.0, False
            for smM in r.SM_Masses:
                smz = smM + d
                st = ion_status(run, smz, IONS[PRIMARY[pol]], s)
                if not st and b > a:
                    sm_raw += float(ion_window_max(run, smz, tol)[a:b].max()); ok = True
            row["SM_Ion"] = name
            row["SM_Raw_Ejection_Intensity"] = sm_raw if ok else np.nan
            row["Echo_Product_Ratio_%"] = (round(100 * raw_I / (raw_I + sm_raw), 1)
                                           if ok and searchable and (raw_I + sm_raw) > 0 else np.nan)
        else:
            row["SM_Ion"], row["SM_Raw_Ejection_Intensity"], row["Echo_Product_Ratio_%"] = "", np.nan, np.nan
        out.append(row)

    res = pd.DataFrame(out)
    res["_r"] = res["Row"].map(plate.rows.index)
    res = res.sort_values(["_r", "Column"]).drop(columns="_r").reset_index(drop=True)
    return PlateAnalysis(plate_id, run, res, P, grid, rel, order_name, first, slot_of, scores, blank, warnings)
