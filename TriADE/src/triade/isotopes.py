"""Theoretical isotope patterns of the product ions, and how well a spectrum matches them.

Used as information (pattern score, confirmation peak, overlay in the spectrum viewer) and, with the
'relative' M+1 rule, to set the M+1 intensity expected from the formula. Peaks are grouped by nominal
offset (M, M+1, M+2, ...), each at the abundance-weighted mean mass shift: this is what a unit-resolution
instrument sees, and a good approximation at QTOF resolution (13C/15N fine structure is not resolved).
"""
from __future__ import annotations

import re
from functools import lru_cache

import numpy as np

# element: [(nominal offset, exact mass shift (Da), abundance)]  (IUPAC representative abundances)
ISOTOPES = {
    "H": [(0, 0.0, 0.999885), (1, 1.006277, 0.000115)],
    "B": [(0, 0.0, 0.801), (-1, -0.996370, 0.199)],
    "C": [(0, 0.0, 0.9893), (1, 1.003355, 0.0107)],
    "N": [(0, 0.0, 0.99636), (1, 0.997035, 0.00364)],
    "O": [(0, 0.0, 0.99757), (1, 1.004217, 0.00038), (2, 2.004245, 0.00205)],
    "F": [(0, 0.0, 1.0)], "Na": [(0, 0.0, 1.0)], "P": [(0, 0.0, 1.0)], "I": [(0, 0.0, 1.0)],
    "Si": [(0, 0.0, 0.92223), (1, 0.999568, 0.04685), (2, 1.996844, 0.03092)],
    "S": [(0, 0.0, 0.9499), (1, 0.999388, 0.0075), (2, 1.995796, 0.0425), (4, 3.995010, 0.0001)],
    "Cl": [(0, 0.0, 0.7576), (2, 1.997050, 0.2424)],
    "Br": [(0, 0.0, 0.5069), (2, 1.997953, 0.4931)],
    "K": [(0, 0.0, 0.932581), (2, 1.998119, 0.067302)],
}
# atoms added to (+) or removed from (-) n·M by each ion
ION_ATOMS = {
    "[M+H]+": {"H": 1}, "[M+Na]+": {"Na": 1}, "[M+K]+": {"K": 1}, "[M+NH4]+": {"N": 1, "H": 4},
    "[M+H-H2O]+": {"H": -1, "O": -1}, "[M+MeOH+H]+": {"C": 1, "H": 5, "O": 1},
    "[M+ACN+H]+": {"C": 2, "H": 4, "N": 1}, "[M+DMSO+H]+": {"C": 2, "H": 7, "O": 1, "S": 1},
    "[2M+H]+": {"H": 1}, "[2M+Na]+": {"Na": 1}, "[M+2H]2+": {"H": 2},
    "[M-H]-": {"H": -1}, "[M+Cl]-": {"Cl": 1}, "[M+HCOO]-": {"C": 1, "H": 1, "O": 2},
    "[M+CH3COO]-": {"C": 2, "H": 3, "O": 2}, "[M-H-H2O]-": {"H": -3, "O": -1}, "[2M-H]-": {"H": -1},
}


def composition(formula: str) -> dict:
    out = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula or ""):
        out[el] = out.get(el, 0) + int(n or 1)
    return out


def ion_composition(formula: str, ion) -> dict:
    c = {el: n * ion.n_m for el, n in composition(formula).items()}
    for el, n in ION_ATOMS.get(ion.name, {}).items():
        c[el] = c.get(el, 0) + n
    return {el: n for el, n in c.items() if n > 0}


@lru_cache(maxsize=4096)
def _pattern(items: tuple, kmax: int) -> tuple:
    dist = {0: (1.0, 0.0)}           # offset -> (probability, probability-weighted shift)
    for el, n in items:
        iso = ISOTOPES.get(el, [(0, 0.0, 1.0)])
        for _ in range(n):
            new = {}
            for k, (p, ws) in dist.items():
                for off, d, a in iso:
                    kk = k + off
                    if kk < 0 or kk > kmax:
                        continue
                    P, W = new.get(kk, (0.0, 0.0))
                    new[kk] = (P + p * a, W + (ws + p * d) * a)
            dist = new
    p0 = dist[0][0]
    return tuple((k, p / p0, (w / p if p else 0.0)) for k, (p, w) in sorted(dist.items()))


def pattern(formula: str, ion, kmax: int = 4) -> list[tuple[int, float, float]]:
    """[(offset k, intensity relative to M, m/z shift from M)] of the ion."""
    comp = ion_composition(formula, ion)
    z = abs(ion.charge)
    return [(k, rel, shift / z) for k, rel, shift in _pattern(tuple(sorted(comp.items())), kmax)]


def confirmation_peak(pat) -> tuple[int, float, float]:
    """The most intense isotope peak after M (M+1 for CHNO, M+2 with Cl/Br)."""
    return max((p for p in pat if p[0] > 0), key=lambda p: p[1], default=(1, 0.0, 1.00336))


def observed_ratios(x: np.ndarray, m_bin: int, pat, bin_da: float, window_da: float = 0.3) -> dict:
    """Intensity of each isotope peak relative to the M peak in one spectrum (dense grid)."""
    out = {}
    if x[m_bin] <= 0:
        return out
    for k, _, shift in pat:
        if k == 0:
            continue
        c = m_bin + shift / bin_da
        lo, hi = int(np.ceil(c - window_da / bin_da)), int(np.floor(c + window_da / bin_da))
        seg = x[max(lo, 0):min(hi, len(x) - 1) + 1]
        out[k] = float(seg.max()) / float(x[m_bin]) if seg.size else 0.0
    return out


def pattern_score(obs: dict, pat, min_rel: float = 0.05) -> float:
    """0-1: mean agreement min(obs/theo, theo/obs) over the isotope peaks expected at >= min_rel of M."""
    vals = []
    for k, rel, _ in pat:
        if k == 0 or rel < min_rel or k not in obs:
            continue
        o = obs[k]
        vals.append(0.0 if o <= 0 else min(o / rel, rel / o))
    return float(np.mean(vals)) if vals else float("nan")


def describe(pat, min_rel: float = 0.02) -> str:
    return "; ".join(f"M+{k} {100 * rel:.0f}%" for k, rel, _ in pat if k > 0 and rel >= min_rel)


def describe_obs(obs: dict, pat, min_rel: float = 0.02) -> str:
    return "; ".join(f"M+{k} {100 * obs[k]:.0f}%" for k, rel, _ in pat if k > 0 and rel >= min_rel and k in obs)
