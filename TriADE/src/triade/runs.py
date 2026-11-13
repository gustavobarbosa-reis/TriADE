"""Reading an Echo-MS run (mzML) into a dense scan × m/z grid, and finding the ejection pulses."""
from __future__ import annotations

import logging
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

try:  # pragma: no cover - depends on the environment
    from pyteomics import mzml as _pyteomics_mzml
    HAVE_PYTEOMICS = True
except ImportError:  # pragma: no cover
    HAVE_PYTEOMICS = False

from . import _mzml

log = logging.getLogger("triade")


@dataclass
class Run:
    path: Path
    t: np.ndarray            # time of each scan (s)
    pol: np.ndarray          # '+' or '-' per scan
    dense: np.ndarray        # (n_scans, n_bins) maximum intensity on an m/z grid
    bin0: int                # index of the first bin (m/z = bin * bin_da)
    bin_da: float
    mz_min: float
    mz_max: float
    tic: np.ndarray = field(default=None)
    file_index: np.ndarray = field(default=None)   # position of each scan among the spectra of the file
    # high-resolution mode: every centroid peak kept (no grid)
    sparse: bool = False
    ppm: float | None = None
    scan_ptr: np.ndarray = field(default=None)    # peaks of scan k = flat[scan_ptr[k]:scan_ptr[k + 1]]
    flat_mz: np.ndarray = field(default=None)
    flat_it: np.ndarray = field(default=None)
    flat_scan: np.ndarray = field(default=None)
    mz_order: np.ndarray = field(default=None)    # flat indices sorted by m/z
    mz_sorted: np.ndarray = field(default=None)
    warnings: list = field(default_factory=list)

    def peaks(self, k: int):
        """Centroid peaks (m/z, intensity) of scan k in the high-resolution mode."""
        a, b = self.scan_ptr[k], self.scan_ptr[k + 1]
        return self.flat_mz[a:b], self.flat_it[a:b]

    def tol(self, mz: float, tol_da: float) -> float:
        """Tolerance in Da at this m/z: ppm in the high-resolution mode, otherwise tol_da."""
        return mz * self.ppm * 1e-6 if self.sparse else tol_da
    xic: dict = field(default_factory=dict)   # cache of ion_window_max per bin range

    @property
    def polarities(self) -> set:
        return set(np.unique(self.pol))

    @property
    def duration_min(self) -> float:
        return float(self.t[-1] - self.t[0]) / 60.0


def fill_missing_times(t: np.ndarray) -> np.ndarray:
    """Scans without 'scan start time' get a time interpolated from their neighbours."""
    miss = ~np.isfinite(t)
    if miss.all():
        raise ValueError("No scan has a 'scan start time'; the ejection pulses cannot be located.")
    if miss.any():
        idx = np.arange(len(t)); ok = ~miss
        dt = np.median(np.diff(t[ok])) if ok.sum() > 1 else 0.0
        t = np.interp(idx, idx[ok], t[ok])
        first, last = idx[ok][0], idx[ok][-1]
        t[:first] = t[first] - dt * (first - idx[:first])
        t[last + 1:] = t[last] + dt * (idx[last + 1:] - last)
        log.warning("  %d scans without 'scan start time': time interpolated.", int(miss.sum()))
    if np.any(np.diff(t) < 0):
        log.warning("  Scan times out of order: the per-well windows may be wrong.")
    return t


def _collect(reader, path: Path):
    """Time, polarity, m/z and intensity arrays of every MS1 spectrum."""
    ts, pols, mzs, its, idx = [], [], [], [], []
    with reader(str(path)) as r:
        for n, s in enumerate(r):
            if int(s.get("ms level", 1)) != 1:
                continue
            idx.append(n)
            scan = s["scanList"]["scan"][0]
            ts.append(_seconds(scan.get("scan start time")))
            pols.append("-" if "negative scan" in s else "+")
            mzs.append(np.asarray(s["m/z array"], float))
            its.append(np.asarray(s["intensity array"], float))
    return ts, pols, mzs, its, idx


def _read_spectra(path: Path):
    """The built-in reader is used first (it needs only lxml). pyteomics is tried only for files it
    cannot read (MS-Numpress compression)."""
    try:
        return _collect(_mzml.MzML, path)
    except NotImplementedError as e:
        if not HAVE_PYTEOMICS:
            raise ValueError(f"{path.name}: {e}") from None
        try:
            return _collect(_pyteomics_mzml.MzML, path)
        except Exception as e2:  # pragma: no cover - depends on the pyteomics installation
            raise ValueError(f"{path.name}: {e} (pyteomics could not read it either: {e2})") from None


def _seconds(value) -> float:
    """'scan start time' in seconds (mzML stores minutes or seconds; pyteomics keeps the unit)."""
    if value is None:
        return np.nan
    unit = str(getattr(value, "unit_info", "minute") or "minute").lower()
    return float(value) * (1.0 if unit.startswith("s") else 60.0)


def read_run(path, bin_da: float = 0.1, settings=None) -> Run:
    path = Path(path)
    log.info("Reading %s", path.name)
    ts, pols, mzs, its, idx = _read_spectra(path)
    if not ts:
        raise ValueError(f"{path.name}: no MS1 spectra found")
    ts = fill_missing_times(np.array(ts))
    if settings is not None and getattr(settings, "instrument", "unit") == "high":
        return _sparse_run(path, ts, pols, mzs, its, idx, settings)
    allmz = np.concatenate([m for m in mzs if m.size]) if any(m.size for m in mzs) else np.array([0.0])
    b_all = np.rint(allmz / bin_da).astype(np.int64)
    bin0, bin1 = int(b_all.min()), int(b_all.max())
    dense = np.zeros((len(ts), bin1 - bin0 + 1), dtype=np.float32)
    for i, (mz, it) in enumerate(zip(mzs, its)):
        if mz.size:
            b = np.rint(mz / bin_da).astype(np.int64) - bin0
            np.maximum.at(dense[i], b, it.astype(np.float32))
    run = Run(path, np.array(ts), np.array(pols), dense, bin0, bin_da, float(allmz.min()), float(allmz.max()))
    run.tic = dense.sum(axis=1)
    run.file_index = np.array(idx)
    pos = int((run.pol == "+").sum())
    log.info("  %d MS1 scans (%d positive, %d negative), m/z %.1f–%.1f, %.1f min, median Δt %.3f s",
             len(ts), pos, len(ts) - pos, run.mz_min, run.mz_max, run.t[-1] / 60, np.median(np.diff(run.t)))
    return run


def _sparse_run(path, ts, pols, mzs, its, idx, s) -> Run:
    """High-resolution mode: keeps the centroid peaks of every scan, indexed by m/z for fast extraction."""
    sizes = np.array([m.size for m in mzs])
    ptr = np.concatenate([[0], np.cumsum(sizes)])
    flat_mz = np.concatenate([np.sort(m) for m in mzs]) if sizes.sum() else np.array([0.0])
    flat_it = np.concatenate([it[np.argsort(m)] for m, it in zip(mzs, its)]).astype(np.float32) \
        if sizes.sum() else np.array([0.0], np.float32)
    flat_scan = np.repeat(np.arange(len(mzs)), sizes)
    order = np.argsort(flat_mz, kind="stable")
    run = Run(path, np.array(ts), np.array(pols), None, 0, 0.0, float(flat_mz.min()), float(flat_mz.max()))
    run.sparse, run.ppm = True, float(s.mass_tolerance_ppm)
    run.scan_ptr, run.flat_mz, run.flat_it, run.flat_scan = ptr, flat_mz, flat_it, flat_scan
    run.mz_order, run.mz_sorted = order, flat_mz[order]
    run.tic = np.array([it.sum() for it in its], dtype=np.float32)
    run.file_index = np.array(idx)
    med = np.median(sizes) if sizes.size else 0
    log.info("  %d MS1 scans, high-resolution mode (%d peaks per scan, median), m/z %.4f-%.4f, %.1f min",
             len(ts), med, run.mz_min, run.mz_max, run.t[-1] / 60)
    if med > 20000:
        run.warnings.append(f"{path.name}: {med:.0f} points per spectrum: the data look like profile data. Centroid "
                            "them in msConvert (filter peakPicking) for the high-resolution mode.")
    return run


def period_from_name(name: str) -> float | None:
    """'... delay 2500 ...' -> 2.5 s between ejections."""
    m = re.search(r"delay\s*(\d+)", name, re.I)
    return int(m.group(1)) / 1000.0 if m else None


def plate_id_from_name(name: str) -> str | None:
    m = re.search(r"plate\s*[_-]?\s*([A-Za-z0-9]+)", name, re.I)
    return m.group(1) if m else None


def fit_ejection_grid(run: Run, period_hint: float | None):
    """Fits a periodic comb (period P, phase) to the normalised TIC. Returns (P, apex times, rel. TIC)."""
    from scipy.ndimage import median_filter
    t = run.t
    n = max(3, int(round(5.0 / np.median(np.diff(t)))))
    rel = run.tic / np.maximum(median_filter(run.tic, size=n, mode="nearest"), 1.0) - 1.0

    def comb(periods):
        return np.array([np.abs(np.sum(rel * np.exp(2j * np.pi * t / P))) for P in periods])

    # free search: the largest period with a high score is the fundamental (P/2, P/3 also score)
    wide = np.arange(0.5, 5.0, 0.0005)
    s_wide = comb(wide)
    P = wide[s_wide >= 0.8 * s_wide.max()].max()
    if period_hint:
        near = np.arange(period_hint * 0.97, period_hint * 1.03, 0.0002)
        s_near = comb(near)
        if s_near.max() >= 0.8 * s_wide.max():
            P = near[np.argmax(s_near)]
        else:
            log.warning("  The period %.3f s (from the file name or the settings) does not match the TIC; "
                        "using the period found in the data, %.4f s.", period_hint, P)
    phi = np.angle(np.sum(rel * np.exp(2j * np.pi * t / P)))
    t0 = (phi / (2 * np.pi) * P) % P
    grid = np.arange(t0, t[-1] + 1e-9, P)
    return float(P), grid, rel


def window_scans(t: np.ndarray, center: float, start: float, end: float):
    a = np.searchsorted(t, center + start, side="left")
    b = np.searchsorted(t, center + end, side="right")
    return int(a), int(b)


def ion_window_max(run: Run, target_mz: float, tol: float) -> np.ndarray:
    """Maximum intensity per scan in target ± tol (all scans). In the high-resolution mode the tolerance is
    the ppm tolerance of the run."""
    if run.sparse:
        tol = run.tol(target_mz, tol)
        key = ("hr", round(target_mz, 5), round(tol, 6))
        if key not in run.xic:
            i0 = np.searchsorted(run.mz_sorted, target_mz - tol, "left")
            i1 = np.searchsorted(run.mz_sorted, target_mz + tol, "right")
            out = np.zeros(len(run.t), dtype=np.float32)
            sel = run.mz_order[i0:i1]
            np.maximum.at(out, run.flat_scan[sel], run.flat_it[sel])
            run.xic[key] = out
        return run.xic[key]
    lo = math.ceil((target_mz - tol) / run.bin_da - 1e-6) - run.bin0
    hi = math.floor((target_mz + tol) / run.bin_da + 1e-6) - run.bin0
    lo, hi = max(lo, 0), min(hi, run.dense.shape[1] - 1)
    if (lo, hi) not in run.xic:
        run.xic[(lo, hi)] = (np.zeros(len(run.t), dtype=np.float32) if hi < lo
                             else run.dense[:, lo:hi + 1].max(axis=1))
    return run.xic[(lo, hi)]


def inspect_mzml(path, n_spectra: int = 40) -> dict:
    """Quick look at the first spectra of a file (no full read): scanned m/z range, polarities, centroid or
    profile data, instrument model and the typical spacing between m/z points. Used for the preview and to
    check the instrument type chosen in the settings."""
    from lxml import etree
    path = Path(path)
    info = {"file": path.name, "mz_min": np.inf, "mz_max": -np.inf, "polarities": set(), "centroid": None,
            "instrument_model": "", "median_spacing_da": np.nan, "spectra_checked": 0}
    spacings = []
    for _, el in etree.iterparse(str(path), events=("end",), tag=("{*}cvParam", "{*}spectrum"), huge_tree=True):
        tag = etree.QName(el).localname
        if tag == "cvParam":
            parent = el.getparent()
            pname = etree.QName(parent).localname if parent is not None else ""
            acc = el.get("accession", "")
            if pname == "instrumentConfiguration" and el.get("name") and not info["instrument_model"] \
                    and acc not in ("MS:1000031",):
                info["instrument_model"] = el.get("name")
            continue
        s = _mzml._spectrum_dict(el)
        el.clear()
        if int(s.get("ms level", 1) or 1) != 1:
            continue
        mz = np.asarray(s["m/z array"], float)
        info["polarities"].add("-" if "negative scan" in s else "+")
        if "centroid spectrum" in s:
            info["centroid"] = True
        elif "profile spectrum" in s:
            info["centroid"] = False
        for k in ("lowest observed m/z", "highest observed m/z"):
            if k in s and s[k] is not None:
                v = float(s[k])
                info["mz_min"], info["mz_max"] = min(info["mz_min"], v), max(info["mz_max"], v)
        if mz.size:
            info["mz_min"], info["mz_max"] = min(info["mz_min"], mz.min()), max(info["mz_max"], mz.max())
            if mz.size > 3:
                d = np.diff(np.sort(mz))
                d = d[d > 0]
                if d.size:
                    spacings.append(np.median(d))
        info["spectra_checked"] += 1
        if info["spectra_checked"] >= n_spectra:
            break
    if spacings:
        info["median_spacing_da"] = float(np.median(spacings))
    return info


def instrument_check(info: dict, instrument: str) -> list[str]:
    """Warnings when the file does not look like the instrument type chosen in the settings."""
    out = []
    model = (info.get("instrument_model") or "").lower()
    hr_words = ("tof", "orbitrap", "exactive", "q exactive", "synapt", "xevo g2", "impact", "timstof", "zenotof",
                "tripletof", "x500", "fticr", "lumos", "exploris", "astral")
    looks_hr = any(w in model for w in hr_words)
    looks_unit = any(w in model for w in ("triple quad", "qtrap", "tsq", "quattro", "6500", "7500", "5500", "4500",
                                          "xevo tq", "single quad", "lcms-80", "lcms-20"))
    sp = info.get("median_spacing_da", np.nan)
    if instrument == "unit" and (looks_hr or (info.get("centroid") and np.isfinite(sp) and sp < 0.02)):
        out.append(f"{info['file']}: the data look like high resolution ({info.get('instrument_model') or 'fine m/z spacing'}"
                   "): consider the QTOF / Orbitrap instrument setting.")
    if instrument == "high":
        if looks_unit:
            out.append(f"{info['file']}: the instrument ({info['instrument_model']}) looks like a unit-resolution "
                       "quadrupole: consider the triple-quadrupole setting.")
        if info.get("centroid") is False:
            out.append(f"{info['file']}: profile (not centroided) data. Convert again with msConvert and the filter "
                       "'peakPicking' (vendor centroiding) for the high-resolution mode.")
    return out
