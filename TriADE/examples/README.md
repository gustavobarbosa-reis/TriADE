# Examples

**`input_template.xlsx`** is the recommended input: one workbook with an Instructions sheet, a Reagents
library (each reagent once, with a code), the Samples list (one row per well, reagents by code), optional
plate-shaped Map sheets and an optional HPLC-MS sheet. The yellow rows are a made-up example. Make one for
another plate format with `triade --template my_template.xlsx --plate-format 1536` or with *Download the input
template* in the interface.

Simpler CSV templates (open them in Excel or any spreadsheet program):

| File | Use |
|---|---|
| `sample_table_template.csv` | the sample table: one row per well |
| `starting_materials_by_row_template.csv` | optional: one starting material per plate row, when the sample table has no `SM_SMILES` |
| `hplc_results_template.csv` | optional HPLC-MS results (an empty conversion = no result, counted as not detected) |

`Injection_Order` is optional: leave the column out when the plate was ejected in a standard pattern
(row by row, serpentine or not, or column by column); the program finds the pattern. Give it for
cherry-picked plates or custom orders (the template above is a serpentine over two rows).

## A complete synthetic plate

No real data are distributed with the program. To try it, write a synthetic plate (made-up molecules and
spectra, with a known "HPLC-MS" conversion for every well):

```bash
triade-synthetic --plate 384 --seed 11 --out example_data
triade --mzml "example_data/Plate 1 delay 2500 synthetic 384-well.mzML" \
        --samples example_data/samples_plate1.csv --hplc example_data/hplc_plate1.csv --out example_results
```

or press *Load a synthetic example* in the browser interface. Options: `--plate 96|384|1536`,
`--polarity positive|negative`, `--injection-order` (shuffled order given in the table), `--period` (s).
