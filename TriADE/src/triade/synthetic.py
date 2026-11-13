"""Synthetic Echo-MS plates for examples and tests (no real data).

A made-up SNAr library: each plate row has one heteroaryl chloride, each column one amine. The product
of row r and column c is the amine attached where the chlorine was. A simple reactivity model sets each
well's conversion, which is written as the "HPLC-MS" truth. The mzML holds one centroided MS1 scan
every `scan_dt_s` with:
  * a solvent/matrix pulse at every ejection (so the ejection grid can be fitted);
  * the product [M+H]+ (or [M-H]-) with its 13C M+1 and a smaller [M+Na]+, scaled by conversion;
  * the starting material [M+H]+ with 37Cl and M+1, scaled by (1 - conversion);
  * random chemical noise;
  * a few traps in wells without product: an impurity on the product [M+H]+ m/z, or one on the product
    [M+Na]+ m/z (found only when [M+Na]+ is selected).

    python -m triade.synthetic --plate 96 --out example_data
"""
from __future__ import annotations

import argparse
import base64
import re
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

from .adducts import IONS
from .chem import molecule
from .plates import get_plate_format

# heteroaryl chlorides written as "Cl" + core; the product is amine + core (ring bond labels 8/9 avoid
# clashes with the amine's ring labels)
CORES = [
    "c8ncc(cn8)C(=O)OCC", "c8nc(C)cc(n8)C", "c8ccc(cn8)C#N", "c8nc9ccccc9o8", "c8nc9ccccc9s8",
    "c8ccnc(n8)C(F)(F)F", "c8ccc(nn8)C", "c8ncc(F)c(n8)OC", "c8nccc9ccccc89", "c8ccc(cc8C#N)C(=O)OC",
    "c8nc(cs8)C(=O)OCC", "c8ccc(cn8)C(=O)N", "c8ncnc9[nH]ccc89", "c8nc(ccc8C(=O)OC)C", "c8ccc9ncccc9n8",
    "c8nc(nc(n8)OC)OC", "c8cc(ncn8)C", "c8ccc(nc8)C(=O)OC", "c8ccc(c(n8)C)C#N", "c8nc9cccnc9[nH]8",
    "c8ncc(cn8)Br", "c8ccc(cc8[N+](=O)[O-])C", "c8nccc(n8)N9CCOCC9", "c8ccc(nc8)S(C)(=O)=O",
    "c8nc(C)nc(c8)C", "c8ccnc9ccccc89", "c8ccc(cn8)OC", "c8nc(CC)cc(n8)CC", "c8ccc9OCOc9c8", "c8ncc(C)cn8",
    "c8cc(C)nc(C)n8", "c8ncc(Cl)cn8",
]
AMINES = [   # written so that the amine nitrogen is the last atom (it bonds to the core)
    "C1COCCN1", "C1CCN(CC1)", "C1CCCN1", "CN(C)", "CN1CCN(CC1)", "c1ccccc1CN", "C1CC1N", "OCCN",
    "OC1CCN(CC1)", "c1ccccc1CN(C)", "CC(C)CN", "C1CSCCN1", "C1CCCCC1N", "c1ccccc1CCN", "C1CC2(OCCO2)CCN1",
    "C1CCOC1CN", "CC(C)(C)OC(=O)N1CCN(CC1)", "o1cccc1CN", "CN(C)CCN", "CCOC(=O)C1CCN(CC1)", "CC(C)N",
    "COCCN", "OC1CCN(C1)", "c1cnccc1CN",
] * 2

MATRIX_MZ = {"+": [(279.09, 4e5), (391.28, 2.5e5), (413.27, 1.5e5)], "-": [(255.23, 3e5), (283.26, 2e5)]}
SODIUM_FRACTION = 0.15


def _n_carbon(formula: str) -> int:
    m = re.search(r"C(\d*)(?![a-z])", formula)
    return int(m.group(1) or 1) if m else 0


def library(plate_format="96", seed=0, rows=None, cols=None) -> pd.DataFrame:
    """Sample table (Plate, Well, Sample_Code, Product_SMILES, SM_SMILES) + Conversion_% (the truth)."""
    plate = get_plate_format(plate_format)
    rng = np.random.default_rng(seed)
    rows = list(rows or plate.rows)
    cols = list(cols or range(1, plate.n_cols + 1))
    row_react = {r: rng.uniform(-1.5, 2.5) for r in rows}
    col_react = {c: rng.uniform(-2.0, 2.0) for c in cols}
    out = []
    for i, r in enumerate(plate.rows):
        if r not in rows:
            continue
        core = CORES[i % len(CORES)]
        for c in cols:
            amine = AMINES[(c - 1) % len(AMINES)]
            z = row_react[r] + col_react[c] + rng.normal(0, 0.8)
            conv = 100 / (1 + np.exp(-2.0 * z)) if rng.random() > 0.08 else 0.0
            conv = 0.0 if conv < 3 else round(conv, 1)
            out.append({"Plate": "1", "Well": f"{r}{c}", "Sample_Code": f"SYN-{r}{c:02d}",
                        "Product_SMILES": amine + core, "SM_SMILES": "Cl" + core, "Conversion_%": conv})
    return pd.DataFrame(out)


def _mzml_header(n: int, instrument_model: str = "") -> str:
    return ('<?xml version="1.0" encoding="utf-8"?>\n'
            '<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0" id="synthetic">\n'
            '<cvList count="2"><cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" '
            'URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>'
            '<cv id="UO" fullName="Unit Ontology" URI="https://raw.githubusercontent.com/bio-ontology-research-group/'
            'unit-ontology/master/unit.obo"/></cvList>\n'
            '<fileDescription><fileContent><cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>'
            '</fileContent></fileDescription>\n'
            '<softwareList count="1"><software id="triade_synthetic" version="1"><cvParam cvRef="MS" '
            'accession="MS:1000799" name="custom unreleased software tool" value="TriADE synthetic"/>'
            '</software></softwareList>\n'
            '<instrumentConfigurationList count="1"><instrumentConfiguration id="IC1"><cvParam cvRef="MS" '
            'accession="MS:1000031" name="instrument model" value=""/>'
            + (f'<cvParam cvRef="MS" accession="MS:1000000" name="{instrument_model}" value=""/>' if instrument_model else '')
            + '</instrumentConfiguration>'
            '</instrumentConfigurationList>\n'
            '<dataProcessingList count="1"><dataProcessing id="DP1"><processingMethod order="1" '
            'softwareRef="triade_synthetic"><cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML" '
            'value=""/></processingMethod></dataProcessing></dataProcessingList>\n'
            '<run id="run1" defaultInstrumentConfigurationRef="IC1">\n'
            f'<spectrumList count="{n}" defaultDataProcessingRef="DP1">\n')


def _array(values: np.ndarray, accession: str, name: str, unit: str) -> str:
    data = base64.b64encode(zlib.compress(values.astype("<f8").tobytes())).decode()
    return (f'<binaryDataArray encodedLength="{len(data)}">'
            '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>'
            '<cvParam cvRef="MS" accession="MS:1000574" name="zlib compression" value=""/>'
            f'<cvParam cvRef="MS" accession="{accession}" name="{name}" value="" {unit}/>'
            f'<binary>{data}</binary></binaryDataArray>')


def write_mzml(path, times_s, spectra, polarity="+", instrument_model="") -> None:
    """spectra: list of (mz array, intensity array), one per scan."""
    pol = ('<cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value=""/>' if polarity == "+" else
           '<cvParam cvRef="MS" accession="MS:1000129" name="negative scan" value=""/>')
    with open(path, "w", encoding="utf-8") as f:
        f.write(_mzml_header(len(spectra), instrument_model))
        for i, (t, (mz, it)) in enumerate(zip(times_s, spectra)):
            f.write(f'<spectrum index="{i}" id="scan={i + 1}" defaultArrayLength="{len(mz)}">'
                    '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>'
                    '<cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>'
                    f'{pol}<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>'
                    '<scanList count="1"><cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>'
                    f'<scan><cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{t / 60:.8f}" '
                    'unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/></scan></scanList>'
                    '<binaryDataArrayList count="2">'
                    + _array(mz, "MS:1000514", "m/z array",
                             'unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"')
                    + _array(it, "MS:1000515", "intensity array",
                             'unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"')
                    + '</binaryDataArrayList></spectrum>\n')
        f.write('</spectrumList>\n</run>\n</mzML>\n')


def simulate_run(lib: pd.DataFrame, order: list, path, polarity="+", period_s=2.5, scan_dt_s=0.17,
                 lead_pulses=6, tail_pulses=6, mz_range=(100.0, 650.0), seed=0, product_scale=6e6,
                 sm_scale=3e6, carryover=0.0, impurity_rate=0.05, adduct_trap_rate=0.12, mz_noise_da=0.02,
                 instrument_model=""):
    """Writes an mzML of the plate ejected in `order` (list of wells); returns the apex times (s)."""
    rng = np.random.default_rng(seed)
    proton = 1.007276
    info = lib.set_index("Well")
    n_pulses = lead_pulses + len(order) + tail_pulses
    t_end = (n_pulses + 1) * period_s
    times = np.arange(0.0, t_end, scan_dt_s) + rng.uniform(0, 0.01)
    apex = [(w, (lead_pulses + k + 0.5) * period_s) for k, w in enumerate(order)]   # a well may appear twice
    sigma = 0.16
    # ions of every well: list of (mz, amplitude)
    ions = {}
    sign = 1 if polarity == "+" else -1
    for w in dict.fromkeys(order):
        r = info.loc[w]
        prod, sm = molecule(r["Product_SMILES"]), molecule(r["SM_SMILES"])
        conv = float(r["Conversion_%"]) / 100
        rf = rng.lognormal(0, 0.5)              # response factor of this product
        lst = []
        if prod is not None and conv > 0:
            a = product_scale * conv * rf
            mz = prod.mass + sign * proton
            c13 = 0.0107 * _n_carbon(prod.formula)
            lst += [(mz, a), (mz + 1.003355, a * c13), (mz + 2.00671, a * c13 * c13 / 2)]
            if polarity == "+":
                na = prod.mass + IONS["[M+Na]+"].delta
                lst += [(na, a * SODIUM_FRACTION), (na + 1.003355, a * SODIUM_FRACTION * c13)]
        if sm is not None and conv < 1:
            a = sm_scale * (1 - conv) * rng.lognormal(0, 0.3)
            mz = sm.mass + sign * proton
            c13 = 0.0107 * _n_carbon(sm.formula)
            lst += [(mz, a), (mz + 1.00336, a * c13), (mz + 1.99705, a * 0.32)]
        # traps in wells without product: an impurity isobaric with the product ion (a false positive
        # whatever the ions selected) or one on the product [M+Na]+ (a false positive only if [M+Na]+ is
        # selected) -- the reason to choose adducts carefully
        is_control = str(r.get("Control", "") or "") in ("negative", "blank")
        if prod is not None and conv == 0 and polarity == "+" and not is_control:
            u = rng.random()
            if u < impurity_rate:
                mz, a = prod.mass + proton + rng.normal(0, 0.15), product_scale * rng.uniform(0.1, 0.5)
                lst += [(mz, a), (mz + 1.00336, a * 0.15)]
            elif u < impurity_rate + adduct_trap_rate:
                mz, a = prod.mass + IONS["[M+Na]+"].delta + rng.normal(0, 0.15), product_scale * rng.uniform(0.1, 0.4)
                lst += [(mz, a), (mz + 1.00336, a * 0.15)]
        ions[w] = [(m, i) for m, i in lst if mz_range[0] <= m <= mz_range[1]]
    pol = "+" if polarity == "+" else "-"
    matrix = MATRIX_MZ[pol]
    pulses = (np.arange(n_pulses) + 0.5) * period_s
    spectra = []
    for t in times:
        mzs, its = [], []
        # matrix ions at every pulse (also before/after the plate)
        k = int(np.argmin(np.abs(pulses - t)))
        g = np.exp(-0.5 * ((t - pulses[k]) / sigma) ** 2)
        for m, a in matrix:
            mzs.append(m + rng.normal(0, mz_noise_da / 2)); its.append(a * (0.15 + g) * rng.lognormal(0, 0.1))
        # wells ejected near this time (this one and a fading tail of the previous one)
        for w, tw in apex:
            dt = t - tw
            if -4 * sigma < dt < 4 * sigma:
                f = np.exp(-0.5 * (dt / sigma) ** 2)
            elif carryover and period_s - 4 * sigma < dt < period_s + 4 * sigma:
                f = carryover * np.exp(-0.5 * ((dt - period_s) / sigma) ** 2)
            else:
                continue
            for m, a in ions[w]:
                v = a * f * rng.lognormal(0, 0.15)
                if v > 500:
                    mzs.append(m + rng.normal(0, mz_noise_da)); its.append(v)
        n_noise = 60
        mzs += list(rng.uniform(*mz_range, n_noise)); its += list(rng.lognormal(np.log(8e3), 0.7, n_noise))
        o = np.argsort(mzs)
        spectra.append((np.asarray(mzs)[o], np.asarray(its)[o]))
    write_mzml(path, times, spectra, pol, instrument_model)
    return apex


def make_dataset(out_dir, plate_format="96", polarity="positive", seed=0, period_s=2.5, rows=None, cols=None,
                 injection_order=False, order=None, plate_id="1", high_resolution=False, replicates=1,
                 controls=False, **kw) -> dict:
    """Writes <out>/samples.csv, <out>/Plate <id> delay <ms> synthetic.mzML and <out>/hplc_plate<id>.csv.
    Returns the paths. injection_order=True writes the Injection_Order column (a shuffled order when
    order is None), as for a cherry-picked plate."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    plate = get_plate_format(plate_format)
    lib = library(plate_format, seed, rows, cols)
    lib["Plate"] = plate_id
    wells = list(lib["Well"])
    if order is None:
        if injection_order:
            order = list(np.random.default_rng(seed + 1).permutation(wells))
        else:   # row-serpentine over the wells in the table
            order = [w for w in next(iter(plate.ejection_orders().values())) if w in set(wells)]
    if controls:   # first well of each row = positive control (product added), last = blank, second-last = negative
        lib["Control"] = ""
        rows_used = sorted({plate.split(w)[0] for w in wells}, key=plate.rows.index)
        cols_used = sorted({plate.split(w)[1] for w in wells})
        for r in rows_used:
            lib.loc[lib["Well"] == f"{r}{cols_used[0]}", ["Control", "Conversion_%"]] = ["positive", 100.0]
            lib.loc[lib["Well"] == f"{r}{cols_used[-1]}", ["Control", "Conversion_%"]] = ["blank", 0.0]
            lib.loc[lib["Well"] == f"{r}{cols_used[-2]}", ["Control", "Conversion_%"]] = ["negative", 0.0]
    if replicates > 1:
        order = list(order) * int(replicates)
        injection_order = True
    pol = "+" if polarity == "positive" else "-"
    mz_path = out / f"Plate {plate_id} delay {int(round(period_s * 1000))} synthetic {plate.name}-well.mzML"
    if high_resolution:
        kw.setdefault("mz_noise_da", 0.0004)
        kw.setdefault("instrument_model", "synthetic QTOF")
        mz_path = mz_path.with_name(mz_path.stem + " QTOF.mzML")
    simulate_run(lib, order, mz_path, pol, period_s=period_s, seed=seed, **kw)
    samples = lib.drop(columns="Conversion_%")
    if replicates > 1:
        n = len(order) // int(replicates)
        samples = pd.concat([samples.assign(Replicate=r + 1, Injection_Order=samples["Well"].map(
            {w: r * n + k + 1 for k, w in enumerate(order[:n])})) for r in range(int(replicates))], ignore_index=True)
    elif injection_order:
        pos = {w: k + 1 for k, w in enumerate(order)}
        samples["Injection_Order"] = samples["Well"].map(pos)
    s_path = out / f"samples_plate{plate_id}.csv"
    samples.to_csv(s_path, index=False)
    h_path = out / f"hplc_plate{plate_id}.csv"
    lib[["Plate", "Well", "Conversion_%"]].to_csv(h_path, index=False)
    return {"mzml": mz_path, "samples": s_path, "hplc": h_path, "order": order, "library": lib}


def main(argv=None):
    p = argparse.ArgumentParser(description="Write a synthetic Echo-MS plate (mzML + sample table + HPLC-MS truth).")
    p.add_argument("--plate", default="96", choices=["96", "384", "1536"])
    p.add_argument("--polarity", default="positive", choices=["positive", "negative"])
    p.add_argument("--period", type=float, default=2.5)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--injection-order", action="store_true", help="shuffled ejection order given in the table")
    p.add_argument("--high-resolution", action="store_true", help="QTOF-like centroid data (m/z noise 0.4 mDa)")
    p.add_argument("--out", default="example_data")
    a = p.parse_args(argv)
    d = make_dataset(a.out, a.plate, a.polarity, a.seed, a.period, injection_order=a.injection_order,
                     high_resolution=a.high_resolution)
    for k in ("mzml", "samples", "hplc"):
        print(f"{k:8s} {d[k]}")


if __name__ == "__main__":
    main()
