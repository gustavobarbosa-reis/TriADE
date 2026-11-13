"""The Excel input template: instructions, reagent library, sample list, plate maps and HPLC-MS sheet,
filled with a small made-up example."""
from __future__ import annotations

from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

from .plates import get_plate_format

HEAD = PatternFill("solid", fgColor="D9E6F7")
EXAMPLE = PatternFill("solid", fgColor="FFF7D6")
THIN = Side(style="thin", color="B0AFAA")
BOX = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)

REAGENTS = [
    ("SM-A", "Clc1ncccn1", "2-chloropyrimidine", "starting material"),
    ("SM-B", "Clc1ccc(cn1)C#N", "6-chloronicotinonitrile", "starting material"),
    ("AM-1", "C1COCCN1", "morpholine", "reagent"),
    ("AM-2", "CNC", "dimethylamine", "reagent"),
    ("AM-3", "C1CCNCC1", "piperidine", "reagent"),
]
SAMPLES = [   # plate, well, code, product, reagent 1, reagent 2, injection order
    ("1", "A1", "S-0001", "C1COCCN1c1ncccn1", "SM-A", "AM-1", 1),
    ("1", "A2", "S-0002", "CN(C)c1ncccn1", "SM-A", "AM-2", 2),
    ("1", "A3", "S-0003", "C1CCN(CC1)c1ncccn1", "SM-A", "AM-3", 3),
    ("1", "B3", "S-0004", "C1COCCN1c1ccc(cn1)C#N", "SM-B", "AM-1", 4),
    ("1", "B2", "S-0005", "CN(C)c1ccc(cn1)C#N", "SM-B", "AM-2", 5),
    ("1", "B1", "S-0006", "C1CCN(CC1)c1ccc(cn1)C#N", "SM-B", "AM-3", 6),
]
HPLC = [("1", "A1", 95.2), ("1", "A2", 40.0), ("1", "A3", 0.0), ("1", "B1", 12.5), ("1", "B2", None), ("1", "B3", 88.1)]

INSTRUCTIONS = [
    ("TriADE — input template", "title"),
    ("Fill this workbook and upload it as the 'Sample table'. The yellow rows are a made-up example: replace "
     "or delete them. Column names must stay as they are.", ""),
    ("", ""),
    ("Sheet 'Reagents' (optional, recommended)", "head"),
    ("A library: each reagent once, with a short code (Code), its SMILES, a name and its role. In the other "
     "sheets you can then write the code instead of the SMILES.", ""),
    ("", ""),
    ("Sheet 'Samples' (required, unless you use the Map sheets)", "head"),
    ("One row per well.", ""),
    ("  Plate — plate ID (1, 2, ...). It must match the plate of each mzML file. Needed when there is more than one plate.", ""),
    ("  Well — A1, B12, AF48 ... (validated for the plate format).", ""),
    ("  Sample_Code — your identifier of the sample (optional).", ""),
    ("  Product_SMILES — the expected product (SMILES, or a code of the Reagents sheet). Required.", ""),
    ("  Reagent_1, Reagent_2, ... — the reagents of the well (code or SMILES), one per column; add columns "
     "Reagent_3, Reagent_4 if needed.", ""),
    ("     Reagent_1 should be the reagent shared by a series of wells (e.g. the starting material of a plate "
     "row): the background control compares each well with the wells that share it.", ""),
    ("     Every reagent is checked as a possible interference on the product ions (e.g. its Na+ or DMSO "
     "adduct falling on the product m/z).", ""),
    ("  Injection_Order — order in which the wells were ejected (1, 2, 3 ...). Optional: leave it empty for "
     "standard patterns (row by row, serpentine, column by column); give it for cherry-picked or custom orders.", ""),
    ("", ""),
    ("Sheets 'Map <field> P<plate>' (optional alternative to the Samples list)", "head"),
    ("A grid shaped like the plate: column numbers in the first row, row letters in the first column. "
     "<field> is Product, Sample_Code, Reagent_1, Reagent_2 ... and <plate> the plate ID, e.g. "
     "'Map Reagent_1 P1'. Copy a Map sheet for each field and plate you need. Values given in the Samples "
     "sheet take precedence.", ""),
    ("", ""),
    ("Sheet 'HPLC-MS' (optional)", "head"),
    ("Plate | Well | Conversion_% (0-100). An empty conversion = no HPLC-MS result; such wells count as not "
     "detected in the comparison. HPLC-MS files can also be uploaded separately.", ""),
    ("", ""),
    ("mzML files", "head"),
    ("One per plate. Convert the instrument files with ProteoWizard msConvert, MS1 spectra, without MS-Numpress "
     "compression. A file name containing 'Plate 2' and 'delay 2500' sets the plate ID and the 2.5 s interval "
     "between ejections. High-resolution data (QTOF/Orbitrap) must be centroided (msConvert filter "
     "'peakPicking').", ""),
]


def _style_header(ws, ncol):
    for j in range(1, ncol + 1):
        c = ws.cell(row=1, column=j)
        c.font, c.fill, c.border = Font(bold=True), HEAD, BOX


def write_template(path, plate_format="384") -> Path:
    plate = get_plate_format(plate_format)
    wb = Workbook()
    ws = wb.active
    ws.title = "Instructions"
    for i, (text, kind) in enumerate(INSTRUCTIONS, 1):
        c = ws.cell(row=i, column=1, value=text)
        c.alignment = Alignment(wrap_text=True, vertical="top")
        if kind == "title":
            c.font = Font(bold=True, size=14)
        elif kind == "head":
            c.font = Font(bold=True, size=11, color="184F95")
    ws.column_dimensions["A"].width = 130

    ws = wb.create_sheet("Reagents")
    ws.append(["Code", "SMILES", "Name", "Role"])
    for r in REAGENTS:
        ws.append(list(r))
    _style_header(ws, 4)
    for row in ws.iter_rows(min_row=2, max_row=1 + len(REAGENTS)):
        for c in row:
            c.fill = EXAMPLE
    for col, w in zip("ABCD", (12, 40, 28, 20)):
        ws.column_dimensions[col].width = w

    ws = wb.create_sheet("Samples")
    ws.append(["Plate", "Well", "Sample_Code", "Product_SMILES", "Reagent_1", "Reagent_2", "Injection_Order"])
    for r in SAMPLES:
        ws.append(list(r))
    _style_header(ws, 7)
    for row in ws.iter_rows(min_row=2, max_row=1 + len(SAMPLES)):
        for c in row:
            c.fill = EXAMPLE
    for col, w in zip("ABCDEFG", (8, 8, 14, 40, 14, 14, 16)):
        ws.column_dimensions[col].width = w
    ws.freeze_panes = "A2"

    for field in ("Product", "Reagent_1", "Reagent_2"):
        ws = wb.create_sheet(f"Map {field} P1")
        ws.cell(row=1, column=1, value="").border = BOX
        for c in range(1, plate.n_cols + 1):
            cell = ws.cell(row=1, column=c + 1, value=c)
            cell.font, cell.fill, cell.border = Font(bold=True), HEAD, BOX
            ws.column_dimensions[get_column_letter(c + 1)].width = 14 if plate.n_cols <= 24 else 9
        for i, r in enumerate(plate.rows, 2):
            cell = ws.cell(row=i, column=1, value=r)
            cell.font, cell.fill, cell.border = Font(bold=True), HEAD, BOX
        ws.freeze_panes = "B2"

    ws = wb.create_sheet("HPLC-MS")
    ws.append(["Plate", "Well", "Conversion_%"])
    for r in HPLC:
        ws.append(list(r))
    _style_header(ws, 3)
    for row in ws.iter_rows(min_row=2, max_row=1 + len(HPLC)):
        for c in row:
            c.fill = EXAMPLE
    path = Path(path)
    wb.save(path)
    return path


def template_hplc(path):
    """The HPLC-MS sheet of a template workbook as a frame (None when absent or empty)."""
    import pandas as pd
    from .samples import key
    path = Path(path)
    if path.suffix.lower() not in (".xlsx", ".xlsm", ".xls"):
        return None
    names = {key(s): s for s in pd.ExcelFile(path).sheet_names}
    sheet = names.get("hplcms") or names.get("hplc") or names.get("lcms")
    if sheet is None:
        return None
    t = pd.read_excel(path, sheet_name=sheet).dropna(how="all")
    return t if len(t) else None
