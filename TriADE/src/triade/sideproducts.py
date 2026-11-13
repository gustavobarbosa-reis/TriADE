"""Search for expected side products in every well.

A rule is written  "name: target: change", for example
    hydrolysis: Reagent_1: -Cl +OH        (the starting material with Cl replaced by OH)
    bis-adduct: Product: +C4H8NO -H       (a second amine on the product)
    oxidation: Product: +O
target is Product or Reagent_1, Reagent_2, ...; change is a list of formulas added (+) or removed (-).
The side product is searched as the protonated / deprotonated molecule with the same spectrum criteria as
the product (M and M+1 in >= min_spectra spectra), and must be clearly above the median signal of the same
m/z in the other wells of the plate.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from ._smiles import MONO
from .adducts import IONS, PRIMARY
from .chem import molecule
from .isotopes import composition
from .runs import ion_window_max, window_scans

_RULE = re.compile(r"^\s*([^:]+?)\s*:\s*(Product|Reagent_\d+)\s*:\s*(.+?)\s*$", re.I)
_TERM = re.compile(r"([+-])\s*([A-Z][A-Za-z0-9]*)")


def parse_rules(rules) -> list[dict]:
    out = []
    for raw in rules or []:
        for line in str(raw).splitlines():
            line = line.split("#")[0].strip()
            if not line:
                continue
            m = _RULE.match(line)
            if not m:
                raise ValueError(f"side-product rule {line!r}: use 'name: Product|Reagent_1: -Cl +OH'")
            name, target, change = m.groups()
            target = "Product" if target.lower() == "product" else "Reagent_" + target.split("_")[1]
            terms = _TERM.findall(change)
            if not terms or re.sub(r"[+-]\s*[A-Z][A-Za-z0-9]*|\s", "", change):
                raise ValueError(f"side-product rule {line!r}: the change must be like '-Cl +OH' or '+O'")
            delta = {}
            for sign, formula in terms:
                for el, n in composition(formula).items():
                    if el not in MONO:
                        raise ValueError(f"side-product rule {line!r}: unknown element {el}")
                    delta[el] = delta.get(el, 0) + (n if sign == "+" else -n)
            key = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_") or f"SP{len(out) + 1}"
            out.append({"name": name, "key": key, "target": target, "change": change,
                        "delta_mass": sum(MONO[el] * n for el, n in delta.items()), "delta": delta})
    keys = [r["key"] for r in out]
    if len(set(keys)) != len(keys):
        raise ValueError("two side-product rules have the same name")
    return out


def _target_mass(row, target) -> float:
    if target == "Product":
        return row.get("Monoisotopic_Mass", np.nan)
    m = molecule(row.get(target, ""))
    return m.mass if m else np.nan


def search(pa, samples: pd.DataFrame, rules: list[dict], s) -> pd.DataFrame:
    """One row per well: m/z, evidence, valid spectra and call for each side product."""
    from .detection import analyse_window
    run, res = pa.run, pa.results
    info = samples.set_index("Well")
    ions = [IONS[PRIMARY[p]] for p in ("+", "-") if p in run.polarities]
    windows = {w: window_scans(run.t, pa.grid[sl], -s.window_pre_s, s.window_post_s)
               for w, sl in pa.slot_of.items() if sl < len(pa.grid)}
    out = pd.DataFrame({"Well": res["Well"]})
    for rule in rules:
        mzs, evs, nvs, dets = [], [], [], []
        for w in res["Well"]:
            row = info.loc[w] if w in info.index else {}
            M = _target_mass(row, rule["target"]) + rule["delta_mass"] if len(row) else np.nan
            best = (np.nan, 0.0, 0, 0)
            if np.isfinite(M) and w in windows:
                a, b = windows[w]
                for ion in ions:
                    mz = ion.mz(M)
                    if not (run.mz_min + s.mass_tolerance_da <= mz <= run.mz_max - 1.5):
                        continue
                    nv, _, _, ev = analyse_window(run, a, b, mz, ion, s)
                    x = ion_window_max(run, mz, s.mass_tolerance_da)
                    plate_bg = np.median([x[aa:bb].max() for aa, bb in windows.values() if bb > aa])
                    own = float(x[a:b].max()) if b > a else 0.0
                    ok = nv >= s.min_spectra and own >= s.background_ratio * max(plate_bg, s.background_floor_cps)
                    if ev >= best[1] or np.isnan(best[0]):
                        best = (round(mz, 4), float(ev), int(nv), int(ok))
            mzs.append(best[0]); evs.append(best[1] if np.isfinite(best[0]) else np.nan)
            nvs.append(best[2]); dets.append(best[3])
        k = rule["key"]
        out[f"SP_{k}_mz"], out[f"SP_{k}_Evidence"] = mzs, evs
        out[f"SP_{k}_Valid_Spectra"], out[f"SP_{k}_Detected"] = nvs, dets
    return out
