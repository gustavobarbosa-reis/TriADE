"""Sample table: which product (and reagents) is in each well, and the injection order.

Accepted inputs (column names are matched ignoring case, spaces and accents):

1. The Excel template (`triade --template` / "Download the input template" in the interface):
     * Reagents  – library: Code | SMILES | Name | Role. Write each reagent once and use its code elsewhere.
     * Samples   – one row per well: Plate | Well | Sample_Code | Product_SMILES | Reagent_1 | Reagent_2 | ...
                   | Injection_Order. Reagent cells hold a code from the Reagents sheet or a SMILES.
     * Map ...   – optional plate-map sheets (a grid shaped like the plate) named
                   "Map <field> P<plate>", e.g. "Map Product P1", "Map Reagent_1 P1", "Map Sample_Code P1".
     * HPLC-MS   – optional results: Plate | Well | Conversion_%.
2. One table (CSV/TSV/Excel) with a row per well, as the Samples sheet. The older SM_SMILES column
   (several SMILES separated by ';') is still accepted and split into Reagent_1, Reagent_2, ...
3. An Excel workbook with one sheet per plate ("Plate 1", "Plate 2", ...).

Reagent_1 should be the reagent shared by a series of wells (e.g. the starting material of a plate row):
it defines the background control. Every reagent is checked as a possible interference.
Reagents can also come from a separate table with one SMILES per plate row (Row | SMILES, optional Plate).
"""
from __future__ import annotations

import logging
import re
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd

from .adducts import ions_for
from .chem import molecule, split_smiles_list
from .plates import PlateFormat

log = logging.getLogger("triade")

ALIASES = {
    "Plate": ("plate", "plateid", "placa", "platename"),
    "Well": ("well", "wells", "position", "positions", "posicao", "posicoes", "poco", "pocos", "wellposition"),
    "Sample_Code": ("samplecode", "sampleid", "code", "codes", "codigo", "codigos", "id", "sample", "compoundid"),
    "Product_SMILES": ("productsmiles", "smiles", "smile", "product", "productsmile", "targetsmiles", "produto"),
    "SM_SMILES": ("smsmiles", "startingmaterialsmiles", "startingmaterials", "startingmaterial", "reagentsmiles",
                  "reactantsmiles", "materialdepartida", "reagents"),
    "Injection_Order": ("injectionorder", "ejectionorder", "order", "injection", "sequence", "ordem"),
    "Sheet_MW": ("totalmolweight", "molweight", "mw", "molecularweight", "massamolecular", "pesomolecular"),
    "Control": ("control", "controltype", "welltype", "wellrole", "controle"),
    "Replicate": ("replicate", "rep", "replica", "replicata", "replicateno"),
}
CONTROL_VALUES = {"positive": "positive", "pos": "positive", "+": "positive", "positivecontrol": "positive",
                  "negative": "negative", "neg": "negative", "-": "negative", "negativecontrol": "negative",
                  "blank": "blank", "branco": "blank", "empty": "blank", "positivo": "positive",
                  "negativo": "negative"}
_REAGENT_COL = re.compile(r"^(?:reagent|reactant|sm|startingmaterial|materialdepartida|reagente)(\d+)(?:smiles|code)?$")
TEMPLATE_SHEETS = {"samples", "reagents"}


def key(text) -> str:
    """Column name without accents, spaces, capitals or punctuation ("Posição" -> "posicao")."""
    t = unicodedata.normalize("NFKD", str(text)).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]", "", t.lower())


def read_table(path, sheet=0) -> pd.DataFrame:
    """CSV/TSV (comma, semicolon or tab, with or without BOM) or Excel."""
    path = Path(path)
    if path.suffix.lower() in (".xlsx", ".xlsm", ".xls"):
        return pd.read_excel(path, sheet_name=sheet)
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        head = f.readline()
    sep = max([",", ";", "\t"], key=head.count) if any(c in head for c in ",;\t") else ","
    return pd.read_csv(path, sep=sep, encoding="utf-8-sig")


def reagent_columns(df: pd.DataFrame) -> list[str]:
    cols = [c for c in df.columns if re.fullmatch(r"Reagent_\d+", str(c))]
    return sorted(cols, key=lambda c: int(c.split("_")[1]))


def _rename(df: pd.DataFrame) -> pd.DataFrame:
    cols = {}
    for c in df.columns:
        k = key(c)
        m = _REAGENT_COL.match(k)
        if m:
            cols[c] = f"Reagent_{int(m.group(1))}"
            continue
        for name, keys in ALIASES.items():
            if k in keys and name not in cols.values():
                cols[c] = name
                break
    # "Injection" is the sample code in some instrument exports: only take it as the order if numeric
    out = df.rename(columns=cols)
    if "Injection_Order" in out.columns and not pd.api.types.is_numeric_dtype(out["Injection_Order"]):
        out = out.rename(columns={"Injection_Order": "Sample_Code" if "Sample_Code" not in out.columns else "_injection"})
    return out


def _header_row(raw: pd.DataFrame) -> int:
    for i in range(min(15, len(raw))):
        keys = {key(v) for v in raw.iloc[i].tolist()}
        if keys & set(ALIASES["Well"]) and keys & set(ALIASES["Product_SMILES"]):
            return i
    raise ValueError("No header with a well column (Well/Position) and a SMILES column was found in the first 15 rows.")


def _sheet_plate(name: str) -> str | None:
    m = re.search(r"(?:plate|placa)\s*[_-]?\s*([A-Za-z0-9]+)\s*$", str(name).strip(), re.I) or \
        re.search(r"(?:plate|placa)\s*[_-]?\s*(\d+)", str(name), re.I)
    return m.group(1) if m else None


# ------------------------------------------------------------------------------------------
# Template workbook
# ------------------------------------------------------------------------------------------
def read_library(path) -> dict:
    """Code -> SMILES from the Reagents sheet (empty when there is none)."""
    path = Path(path)
    if path.suffix.lower() not in (".xlsx", ".xlsm", ".xls"):
        return {}
    sheets = {key(s): s for s in pd.ExcelFile(path).sheet_names}
    if "reagents" not in sheets:
        return {}
    t = pd.read_excel(path, sheet_name=sheets["reagents"])
    keys = {key(c): c for c in t.columns}
    cc = next((keys[k] for k in ("code", "reagentcode", "id", "codigo") if k in keys), None)
    sc = next((keys[k] for k in ("smiles", "reagentsmiles") if k in keys), None)
    if cc is None or sc is None:
        raise ValueError(f"{path.name}, sheet Reagents: the columns Code and SMILES are needed")
    lib = {}
    for c, s in zip(t[cc], t[sc]):
        if pd.notna(c) and pd.notna(s) and str(c).strip():
            code = str(c).strip()
            if code in lib and lib[code] != str(s).strip():
                raise ValueError(f"{path.name}, sheet Reagents: the code {code!r} is used for two different SMILES")
            lib[code] = str(s).strip()
    return lib


_MAP_SHEET = re.compile(r"^map\s+(.+?)\s+(?:p|plate|placa)\s*([A-Za-z0-9]+)$", re.I)


def _map_field(text: str) -> str | None:
    k = key(text)
    m = _REAGENT_COL.match(k)
    if m:
        return f"Reagent_{int(m.group(1))}"
    for name in ("Product_SMILES", "Sample_Code", "Injection_Order"):
        if k in ALIASES[name] or k == key(name):
            return name
    return None


def read_plate_map(raw: pd.DataFrame, plate: PlateFormat) -> dict:
    """Grid sheet -> {well: value}. The grid has the column numbers in one row and the row letters in
    the first column below it; anything above (title, notes) is ignored."""
    for i in range(min(20, len(raw))):
        nums = [str(v).strip() for v in raw.iloc[i].tolist()]
        cols = {j: int(float(v)) for j, v in enumerate(nums) if re.fullmatch(r"\d+(\.0)?", v)}
        if len(cols) >= 2 and sorted(cols.values())[:2] == [1, 2]:
            out = {}
            for r in range(i + 1, len(raw)):
                lab = str(raw.iat[r, 0]).strip().upper()
                if lab not in plate.rows:
                    continue
                for j, c in cols.items():
                    v = raw.iat[r, j]
                    if pd.notna(v) and str(v).strip() and c <= plate.n_cols:
                        out[f"{lab}{c}"] = v.strip() if isinstance(v, str) else v
            return out
    raise ValueError("no row with the column numbers 1, 2, 3, ... was found")


def _load_template(path: Path, plate: PlateFormat) -> pd.DataFrame:
    xl = pd.ExcelFile(path)
    frames = []
    names = {key(s): s for s in xl.sheet_names}
    if "samples" in names:
        t = _rename(pd.read_excel(path, sheet_name=names["samples"]).dropna(axis=1, how="all"))
        t = t.dropna(how="all")
        if len(t):
            frames.append(t)
    maps = {}
    for sheet in xl.sheet_names:
        m = _MAP_SHEET.match(str(sheet).strip())
        if not m:
            continue
        field_name = _map_field(m.group(1))
        if field_name is None:
            raise ValueError(f"{path.name}, sheet {sheet!r}: unknown field {m.group(1)!r} (use Product, Sample_Code, "
                             "Reagent_1, Reagent_2, ... or Injection_Order)")
        try:
            grid = read_plate_map(pd.read_excel(path, sheet_name=sheet, header=None), plate)
        except ValueError as e:
            raise ValueError(f"{path.name}, sheet {sheet!r}: {e}") from None
        for w, v in grid.items():
            maps.setdefault((m.group(2), w), {})[field_name] = v
    if maps:
        mt = pd.DataFrame([{"Plate": p, "Well": w, **vals} for (p, w), vals in maps.items()])
        mt["Well"] = mt["Well"].map(plate.parse)
        if frames:
            base = frames[0]
            base["Plate"] = base["Plate"].map(_plate_str) if "Plate" in base.columns else "1"
            base["Well"] = base["Well"].map(plate.parse)
            merged = base.merge(mt, on=["Plate", "Well"], how="outer", suffixes=("", "_map"))
            for c in [c for c in merged.columns if c.endswith("_map")]:
                orig = c[:-4]
                merged[orig] = merged[orig].where(merged[orig].notna() & (merged[orig].astype(str).str.strip() != ""),
                                                  merged[c])
                merged = merged.drop(columns=c)
            frames = [merged]
        else:
            frames = [mt]
    if not frames:
        raise ValueError(f"{path.name}: the Samples sheet is empty and there are no 'Map ...' sheets")
    return pd.concat(frames, ignore_index=True)


def _plate_str(v) -> str:
    if isinstance(v, (int, float, np.integer, np.floating)) and not pd.isna(v) and float(v).is_integer():
        return str(int(v))
    return str(v).strip()


# ------------------------------------------------------------------------------------------
# Loading
# ------------------------------------------------------------------------------------------
def load_samples(path, plate: PlateFormat, default_plate: str = "1") -> pd.DataFrame:
    """Reads the sample table into a standard frame (one row per well)."""
    path = Path(path)
    library = {}
    if path.suffix.lower() in (".xlsx", ".xlsm", ".xls"):
        sheets = pd.ExcelFile(path).sheet_names
        if TEMPLATE_SHEETS & {key(s) for s in sheets} or any(_MAP_SHEET.match(str(s).strip()) for s in sheets):
            library = read_library(path)
            df = _load_template(path, plate)
        else:
            frames = []
            plate_sheets = [(s, _sheet_plate(s)) for s in sheets if _sheet_plate(s)]
            use = plate_sheets if plate_sheets else [(sheets[0], None)]
            for sheet, pid in use:
                raw = pd.read_excel(path, sheet_name=sheet, header=None)
                h = _header_row(raw)
                d = _rename(pd.read_excel(path, sheet_name=sheet, header=h).dropna(axis=1, how="all"))
                if pid is not None and "Plate" not in d.columns:
                    d["Plate"] = pid
                frames.append(d)
            df = pd.concat(frames, ignore_index=True)
    else:
        df = _rename(read_table(path))
    return standardise(df, plate, default_plate, source=path.name, library=library)


def _resolve(value, library: dict) -> tuple[str, str]:
    """(SMILES, code) of a reagent/product cell: a library code or a SMILES."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "", ""
    v = str(value).strip()
    if not v or v.lower() in ("nan", "none"):
        return "", ""
    if v in library:
        return library[v], v
    up = {k.upper(): k for k in library}
    if v.upper() in up:
        k = up[v.upper()]
        return library[k], k
    return v, ""


def standardise(df: pd.DataFrame, plate: PlateFormat, default_plate: str = "1", source: str = "sample table",
                library: dict | None = None) -> pd.DataFrame:
    library = library or {}
    for need in ("Well", "Product_SMILES"):
        if need not in df.columns:
            raise ValueError(f"{source}: the column {need} is missing (accepted names: {', '.join(ALIASES[need])})")
    df = df.dropna(subset=["Well"]).copy()
    if "Plate" not in df.columns:
        df["Plate"] = default_plate
    df["Plate"] = df["Plate"].map(lambda v: default_plate if pd.isna(v) else _plate_str(v))
    wells = df["Well"].map(plate.parse)
    bad = df.loc[wells.isna(), "Well"].astype(str).tolist()
    if bad:
        raise ValueError(f"{source}: {len(bad)} well name(s) are not valid for a {plate.name}-well plate: "
                         + ", ".join(bad[:10]) + (" ..." if len(bad) > 10 else ""))
    df["Well"] = wells
    if "Replicate" in df.columns:
        df["Replicate"] = pd.to_numeric(df["Replicate"], errors="coerce").fillna(1).astype(int)
    else:
        df["Replicate"] = 1
    dup = df[df.duplicated(["Plate", "Well", "Replicate"], keep=False)]
    if not dup.empty:
        raise ValueError(f"{source}: wells listed more than once: "
                         + ", ".join(f"plate {p} {w}" for p, w in dup[["Plate", "Well"]].drop_duplicates().values[:10])
                         + ". For replicate ejections of the same well, number them in a Replicate column (1, 2, 3).")
    if df["Replicate"].max() > 1 and "Injection_Order" not in df.columns:
        raise ValueError(f"{source}: replicate ejections need the Injection_Order column (the order of every ejection)")
    if "Sample_Code" not in df.columns:
        df["Sample_Code"] = ""
    df["Sample_Code"] = df["Sample_Code"].map(lambda v: "" if pd.isna(v) else _plate_str(v))
    df["Product_SMILES"] = [_resolve(v, library)[0] for v in df["Product_SMILES"]]

    # reagents: Reagent_k columns, and/or the older SM_SMILES column split on ';'
    reag = reagent_columns(df)
    if "SM_SMILES" in df.columns:
        extra = [split_smiles_list(v) for v in df["SM_SMILES"]]
        start = len(reag)
        for k in range(max((len(e) for e in extra), default=0)):
            col = f"Reagent_{start + k + 1}"
            df[col] = [e[k] if k < len(e) else "" for e in extra]
        df = df.drop(columns="SM_SMILES")
        reag = reagent_columns(df)
    for c in reag:
        resolved = [_resolve(v, library) for v in df[c]]
        df[c] = [r[0] for r in resolved]
        df[f"{c}_Code"] = [r[1] for r in resolved]
    df["SM_SMILES"] = ["; ".join(s for s in row if s) for row in df[reag].itertuples(index=False)] if reag else ""
    if "Injection_Order" in df.columns:
        df["Injection_Order"] = pd.to_numeric(df["Injection_Order"], errors="coerce")
    else:
        df["Injection_Order"] = np.nan
    if "Control" in df.columns:
        bad_c = [v for v in df["Control"].dropna().unique() if str(v).strip() and key(v) not in CONTROL_VALUES]
        if bad_c:
            raise ValueError(f"{source}: Control must be positive, negative or blank (found {', '.join(map(str, bad_c[:5]))})")
        df["Control"] = df["Control"].map(lambda v: CONTROL_VALUES.get(key(v), "") if pd.notna(v) else "")
    else:
        df["Control"] = ""
    keep = ["Plate", "Well", "Replicate", "Sample_Code", "Product_SMILES", "SM_SMILES", "Injection_Order", "Control"] \
        + reag + [f"{c}_Code" for c in reag]
    if "Sheet_MW" in df.columns:
        keep.append("Sheet_MW")
    out = df[keep].reset_index(drop=True)
    out["Background_Key"] = background_keys(out)
    log.info("Sample table %s: %d wells on %d plate(s); %d reagent column(s)%s", source, len(out),
             out["Plate"].nunique(), len(reag), f"; {len(library)} reagents in the library" if library else "")
    return out


def background_keys(df: pd.DataFrame) -> list:
    """The reagent that defines each well's background group: per plate, the first reagent column whose
    values are shared by series of wells (median group size >= 3). Empty when none is shared."""
    keys = [""] * len(df)
    reag = reagent_columns(df)
    for pid, sub in df.groupby("Plate"):
        chosen = None
        for c in reag:
            v = sub[c][sub[c] != ""]
            if len(v) and v.value_counts().median() >= 3:
                chosen = c
                break
        if chosen:
            for i in sub.index:
                keys[i] = sub.at[i, chosen]
    return keys


def apply_row_starting_materials(df: pd.DataFrame, path, plate: PlateFormat) -> pd.DataFrame:
    """Fills the reagents from a table with one starting material per plate row (Row, SMILES, optional
    Plate); wells that already have reagents keep theirs."""
    t = _rename(read_table(path))
    keys = {key(c): c for c in t.columns}
    rc = next((keys[k] for k in ("row", "linha", "fileira", "platerow") if k in keys), None)
    sc = "Product_SMILES" if "Product_SMILES" in t.columns else ("SM_SMILES" if "SM_SMILES" in t.columns else None)
    if rc is None or sc is None:
        raise ValueError(f"{Path(path).name}: expected the columns Row and SMILES")
    t = pd.DataFrame({"Row": t[rc].astype(str).str.strip().str.upper(), "SMILES": t[sc].astype(str).str.strip(),
                      "Plate": t["Plate"].astype(str) if "Plate" in t.columns else ""})
    out = df.copy()
    if "Reagent_1" not in out.columns:
        out["Reagent_1"], out["Reagent_1_Code"] = "", ""
    row_of = out["Well"].map(lambda w: plate.split(w)[0])
    for i in out.index:
        if out.at[i, "SM_SMILES"]:
            continue
        cand = t[(t["Row"] == row_of[i]) & t["Plate"].isin(["", out.at[i, "Plate"]])]
        if not cand.empty:
            smis = split_smiles_list(cand.sort_values("Plate").iloc[-1]["SMILES"])
            out.at[i, "Reagent_1"] = smis[0] if smis else ""
            out.at[i, "SM_SMILES"] = "; ".join(smis)
    out["Background_Key"] = background_keys(out)
    log.info("Starting materials by row from %s: %d wells filled", Path(path).name, int((out["SM_SMILES"] != "").sum()))
    return out


def add_masses(df: pd.DataFrame, ion_names, polarity: str) -> pd.DataFrame:
    """Formula, monoisotopic mass and theoretical m/z of each selected ion; reagent masses."""
    df = df.copy()
    forms, mono, avg, notes = [], [], [], []
    for smi in df["Product_SMILES"]:
        m = molecule(smi)
        if m is None:
            forms.append(""); mono.append(np.nan); avg.append(np.nan)
            notes.append("missing SMILES" if not smi or smi.lower() in ("nan", "none") else "invalid SMILES")
        else:
            forms.append(m.formula); mono.append(m.mass); avg.append(m.average_mass); notes.append(m.note)
    df["Formula"], df["Monoisotopic_Mass"], df["SMILES_Note"] = forms, mono, notes
    for ion in ions_for(ion_names, polarity):
        df[f"{ion.key}_Theoretical_mz"] = [ion.mz(M) if np.isfinite(M) else np.nan for M in df["Monoisotopic_Mass"]]
    n_fix = int(df["SMILES_Note"].str.startswith("radical").sum())
    if n_fix:
        log.warning("%d SMILES with a radical atom corrected (H added).", n_fix)
    if "Sheet_MW" in df.columns:   # the spreadsheet mass checks the SMILES (warning only)
        mw = pd.to_numeric(df["Sheet_MW"], errors="coerce")
        bad = df.loc[(mw - np.array(avg, dtype=float)).abs() > 0.5, "Well"].tolist()
        if bad:
            log.warning("Molecular weight in the table differs by more than 0.5 from the SMILES in %d well(s): %s",
                        len(bad), ", ".join(bad[:20]) + (" ..." if len(bad) > 20 else ""))
    sm_masses, sm_forms = [], []
    for text in df["SM_SMILES"]:
        ms = [molecule(s) for s in split_smiles_list(text)]
        sm_masses.append([m.mass for m in ms if m is not None])
        sm_forms.append("; ".join(m.formula for m in ms if m is not None))
    df["SM_Masses"], df["SM_Formula"] = sm_masses, sm_forms
    return df


def preview(df: pd.DataFrame, ion_names, polarity: str, scan_ranges: dict | None = None):
    """What the program understood from the inputs, before the analysis.
    Returns (table with one row per well, list of problems). scan_ranges: plate -> (mz_min, mz_max)."""
    x = add_masses(df, ion_names, polarity)
    ions = ions_for(ion_names, polarity)
    issues = []
    reag = reagent_columns(x)
    for r in x.itertuples(index=False):
        where = f"plate {r.Plate} {r.Well}"
        if not r.Product_SMILES:
            issues.append(f"{where}: no product SMILES")
        elif not np.isfinite(r.Monoisotopic_Mass):
            issues.append(f"{where}: the product {r.Product_SMILES!r} is not a valid SMILES nor a code of the "
                          "Reagents sheet")
        for c in reag:
            smi = getattr(r, c)
            if smi and molecule(smi) is None:
                issues.append(f"{where}: {c} {smi!r} is not a valid SMILES nor a code of the Reagents sheet")
        if scan_ranges and r.Plate in scan_ranges and np.isfinite(r.Monoisotopic_Mass):
            lo, hi = scan_ranges[r.Plate]
            out = [i.name for i in ions if not (lo <= getattr(r, f"{i.key}_Theoretical_mz") <= hi)]
            if out:
                issues.append(f"{where}: {', '.join(out)} outside the scanned range m/z {lo:.0f}-{hi:.0f} "
                              "(not measurable: the well will be sent to HPLC-MS)")
    cols = ["Plate", "Well", "Sample_Code", "Product_SMILES", "Formula", "Monoisotopic_Mass"] \
        + [f"{i.key}_Theoretical_mz" for i in ions] + [c for c in reag] + ["SM_Formula", "Background_Key", "Injection_Order"]
    table = x[[c for c in cols if c in x.columns]].copy()
    for c in table.columns:
        if c.endswith("_mz") or c == "Monoisotopic_Mass":
            table[c] = table[c].round(4)
    return table, issues


def plate_grid(df: pd.DataFrame, plate: PlateFormat, plate_id: str, column: str) -> pd.DataFrame:
    """A column of the sample table laid out like the plate (rows x columns), for the preview."""
    sub = df[df["Plate"] == plate_id]
    grid = pd.DataFrame("", index=list(plate.rows), columns=list(range(1, plate.n_cols + 1)))
    for w, v in zip(sub["Well"], sub[column]):
        r, c = plate.split(w)
        grid.at[r, c] = "" if pd.isna(v) else str(v)
    used_rows = [r for r in grid.index if (grid.loc[r] != "").any()]
    used_cols = [c for c in grid.columns if (grid[c] != "").any()]
    return grid.loc[used_rows, used_cols] if used_rows else grid.iloc[:0, :0]
