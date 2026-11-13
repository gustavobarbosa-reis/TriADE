"""Data for the spectrum viewer and the spectra PDF: the ejection trace around a well, any spectrum of the
run (as analysed or as stored in the mzML) and the m/z markers (product ions, isotopes, reagent ions)."""
from __future__ import annotations

import numpy as np

from . import _mzml, isotopes
from .adducts import SM_INTERFERENCE_IONS, ions_for
from .runs import ion_window_max, window_scans


def well_trace(pa, well: str, s, neighbours: int = 2) -> dict:
    """Times, TIC and the extracted signal of each selected ion around the ejection of a well."""
    run, res = pa.run, pa.results.set_index("Well")
    slot = pa.slot_of.get(well)
    if slot is None or slot >= len(pa.grid):
        return {}
    apex = float(pa.grid[slot])
    span = (neighbours + 0.5) * pa.period
    i0, i1 = np.searchsorted(run.t, apex - span), np.searchsorted(run.t, apex + span)
    a, b = window_scans(run.t, apex, -s.window_pre_s, s.window_post_s)
    traces = {}
    for ion in ions_for(s.ions, s.polarity):
        if ion.polarity not in run.polarities:
            continue
        col = f"{ion.key}_Theoretical_mz"
        mz = res.at[well, col] if col in res.columns else np.nan
        if np.isfinite(mz) and run.mz_min <= mz <= run.mz_max:
            y = ion_window_max(run, mz, s.mass_tolerance_da)[i0:i1].astype(float)
            y[run.pol[i0:i1] != ion.polarity] = 0.0          # polarity switching: only scans of this polarity
            traces[ion.name] = y
    by_slot = {sl: w for w, sl in pa.slot_of.items()}
    labels = [(float(pa.grid[sl]), by_slot[sl]) for sl in range(max(slot - neighbours, 0),
                                                                 min(slot + neighbours + 1, len(pa.grid)))
              if sl in by_slot]
    in_window = list(range(a, b))
    best = max(in_window, key=lambda k: max((traces[n][k - i0] for n in traces), default=run.tic[k]),
               default=None) if in_window else None
    return {"t": run.t[i0:i1], "scan": np.arange(i0, i1), "tic": run.tic[i0:i1].astype(float), "ions": traces,
            "window": (apex - s.window_pre_s, apex + s.window_post_s),
            "blank": (apex + pa.blank[0], apex + pa.blank[1]), "apex": apex, "labels": labels,
            "window_scans": in_window, "best_scan": best}


def spectrum(pa, scan: int, mode: str = "analysed"):
    """(m/z, intensity) of one scan of the run: 'analysed' = the program's 0.1 Da grid (what decided the
    level); 'raw' = the spectrum as stored in the mzML file."""
    run = pa.run
    if mode == "raw":
        sp = _mzml.read_spectrum_at(run.path, int(run.file_index[scan]))
        return np.asarray(sp["m/z array"], float), np.asarray(sp["intensity array"], float)
    if run.sparse:
        mz, it = run.peaks(scan)
        return mz.astype(float), it.astype(float)
    x = run.dense[scan]
    nz = np.nonzero(x)[0]
    return (nz + run.bin0) * run.bin_da, x[nz].astype(float)


def mean_spectrum(pa, scans, how: str = "max"):
    """Maximum or mean of several scans (on the analysis grid; peaks merged at 0.1 mDa in high resolution)."""
    run = pa.run
    scans = list(scans)
    if run.sparse:
        import pandas as pd
        parts = [run.peaks(k) for k in scans]
        if not parts:
            return np.array([]), np.array([])
        d = pd.DataFrame({"mz": np.round(np.concatenate([p[0] for p in parts]), 4),
                          "it": np.concatenate([p[1] for p in parts]).astype(float)})
        g = d.groupby("mz")["it"]
        v = g.max() if how == "max" else g.sum() / len(scans)
        return v.index.to_numpy(float), v.to_numpy(float)
    block = run.dense[scans]
    x = block.max(axis=0) if how == "max" else block.mean(axis=0)
    nz = np.nonzero(x)[0]
    return (nz + run.bin0) * run.bin_da, x[nz].astype(float)


def markers(pa, well: str, s) -> list[dict]:
    """Theoretical m/z to draw on the spectrum: each selected product ion, its isotope peaks (with the
    expected intensity relative to M) and the reagent ions checked as interferences."""
    row = pa.results.set_index("Well").loc[well]
    out = []
    formula = row.get("Formula", "")
    M = row.get("Monoisotopic_Mass", np.nan)
    pols = pa.run.polarities
    for ion in ions_for(s.ions, s.polarity):
        if ion.polarity not in pols or not np.isfinite(M):
            continue
        mz = ion.mz(M)
        out.append({"mz": mz, "label": ion.name, "kind": "ion", "rel": 1.0, "ion": ion.name})
        if formula:
            for k, rel, shift in isotopes.pattern(formula, ion):
                if k > 0 and rel >= 0.01:
                    out.append({"mz": mz + shift, "label": f"{ion.name} M+{k}", "kind": "isotope", "rel": rel,
                                "ion": ion.name})
    sm_masses = row.get("SM_Masses", None)
    if sm_masses is None:
        from .chem import molecule, split_smiles_list
        sm_masses = [m.mass for m in (molecule(x) for x in split_smiles_list(row.get("SM_SMILES", ""))) if m]
    for i, smM in enumerate(sm_masses or [], 1):
        for pol in pols:
            for name, d in SM_INTERFERENCE_IONS[pol].items():
                out.append({"mz": smM + d, "label": f"reagent {i} {name.replace('SM', 'R')}", "kind": "reagent",
                            "rel": np.nan, "ion": ""})
    return out


def view_range(marks: list[dict], margin: float = 4.0) -> tuple[float, float]:
    prod = [m["mz"] for m in marks if m["kind"] in ("ion", "isotope")]
    if not prod:
        return (0.0, 0.0)
    return (min(prod) - margin, max(prod) + margin)
