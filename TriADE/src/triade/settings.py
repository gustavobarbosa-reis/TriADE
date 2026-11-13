"""All analysis settings in one place. Every field can be changed from the web interface, the
command line (--settings file.json) or Python, and is saved with the results (settings.json)."""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path

from .adducts import IONS, default_ions
from .plates import PLATE_FORMATS


@dataclass
class Settings:
    # --- instrument, plate and acquisition -------------------------------------------------
    instrument: str = "unit"             # "unit" (triple/single quadrupole) or "high" (QTOF/Orbitrap, experimental)
    plate_format: str = "384"            # "96", "384" or "1536"
    polarity: str = "positive"           # "positive", "negative" or "both"
    ions: list = field(default_factory=lambda: default_ions("positive"))   # ions searched for each product
    ejection_period_s: float | None = None   # time between ejections; None = from 'delay NNNN' in the file name or the data
    ejection_order: str = "auto"         # "auto", "table" (Injection_Order column) or one of the plate patterns
    first_slot: int | None = None        # force the pulse of the first well (0 = first pulse found)

    # --- detection criteria ----------------------------------------------------------------
    mass_tolerance_da: float = 0.5       # ± around the theoretical m/z
    isotope_shift_tolerance_da: float = 0.5   # ± on the observed M -> M+1 shift (shift = 1/|z| Da)
    min_spectra: int = 2                 # the same validated ion in at least this many spectra
    cutoff_cps: float = 200_000          # minimum M-peak intensity (also the level-4 threshold)
    m1_fraction_of_cutoff: float = 0.5   # the M+1 peak needs >= this fraction of the cutoff ('absolute' rule)
    m1_rule: str = "absolute"            # "absolute" (as above) or "relative" (M+1 vs the M+1 expected from the formula)
    m1_relative_fraction: float = 0.7    # 'relative' rule: M+1 >= this x (expected M+1/M) x M ...
    m1_floor_cps: float = 20_000         # ... and >= this intensity
    mass_tolerance_ppm: float = 10.0     # high-resolution mode: ± ppm around each theoretical m/z

    # --- background control -------------------------------------------------------------
    # Compares each well with the same m/z in wells that share its starting material (or row, column).
    background_group: str = "starting_material"   # "starting_material", "row", "column" or "off"
    background_ratio: float = 2.0        # ejection signal >= this x the group background
    background_floor_cps: float = 20_000
    background_min_distance_da: float = 1.5   # wells whose own product ions are this close are left out

    # --- confidence levels ----------------------------------------------------------------
    level5_cps: float = 1_000_000        # evidence >= this -> 5 Present
    level3_cps: float = 80_000           # evidence >= this (below the cutoff) -> 3 Check spectrum
    level2_raw_max_cps: float = 500_000  # raw [M+H]+/[M-H]- above this without passing -> 3
    level1_raw_max_cps: float = 100_000  # raw at or below this -> 1 Not detected; between -> 2 Likely absent
    carryover_ratio: float = 5.0         # previous ejection >= this x the well's signal -> level 3
    sm_interference_da: float = 1.0      # starting-material ion this close to a product ion -> warning / level 3

    # --- ejection windows (s from the pulse apex) ---------------------------------------------
    window_pre_s: float = 0.45
    window_post_s: float = 0.60
    blank_start_s: float = 1.00
    align_threshold_cps: float = 100_000
    mz_bin_da: float = 0.1

    # --- HPLC-MS comparison and triage ---------------------------------------------------------
    hplc_product_min_pct: float = 5.0    # HPLC-MS conversion from which a product counts as present
    hplc_run_min: float = 5.0            # one HPLC-MS injection (min)
    hplc_flow_ml_min: float = 0.6        # HPLC-MS flow rate (mL/min)
    reference_label: str = "HPLC-MS"
    triage_scenario: int = 4             # scenario used for the HPLC-MS worklist (see triage.SCENARIOS)
    target_false_present_pct: float = 1.0   # threshold calibration: acceptable false presents (% of wells)

    # --- replicates, controls, side products ------------------------------------------------------
    replicate_rule: str = "majority"     # replicate ejections: a well is detected if "any", "majority" or "all" detect it
    control_min_positive_pct: float = 100.0   # plate accepted if >= this % of positive controls are detected ...
    control_max_negative_pct: float = 0.0     # ... and <= this % of negative controls / blanks are detected
    side_products: list = field(default_factory=list)   # rules "name: target: change", e.g. "hydrolysis: Reagent_1: -Cl +OH"

    # --- HPLC-MS export ------------------------------------------------------------------------
    transfer_plate_format: str = "96"    # destination plate of the transfer list
    transfer_volume_nl: float | None = None   # volume per transfer (empty = not written)
    blank_every: int = 0                 # injection sequence: a blank every N samples (0 = none)
    internal_standard_mz: float | None = None   # optional internal standard m/z for the run QC

    # --- figures ----------------------------------------------------------------------------------
    palette: str = "default"             # default, viridis, cividis, okabe-ito, grayscale, custom
    custom_colors: dict = field(default_factory=dict)   # palette 'custom': role -> colour (see figures.ROLES)
    figure_format: str = "jpg"           # jpg, png, tiff, svg, pdf
    figure_dpi: int = 300
    cell_values: bool = True             # write the values inside the cells of the intensity / conversion maps
    font_scale: float = 1.0              # text size of every figure (1 = normal; 0.8 smaller, 1.5 larger)
    spectra_pdf: str = "check"           # spectra of each well in a PDF: "none", "check" (level 3), "levels", "all"
    figure_width_mm: float | None = None  # publication mode: figure width (85 = one journal column, 180 = two)
    font_family: str = "default"         # default, Arial, Helvetica, Times New Roman, DejaVu Sans
    show_titles: bool = True             # False: no title/subtitle in the figures (captions written to a text file)
    well_shape: str = "square"           # plate maps: "square" or "circle"
    axis_reagent_labels: bool = True     # plate maps: reagent codes next to the row letters / column numbers
    highlight_wells: bool = True         # plate maps: outline discrepancies, reviewed wells, controls and QC flags
    report: bool = True                  # write Report.html and Report.pdf
    spectra_pdf_levels: list = field(default_factory=lambda: [3])   # levels for spectra_pdf = "levels"

    # ---------------------------------------------------------------------------------------
    def validate(self) -> "Settings":
        if str(self.plate_format) not in PLATE_FORMATS:
            raise ValueError(f"plate_format must be one of {', '.join(PLATE_FORMATS)}")
        self.plate_format = str(self.plate_format)
        if self.polarity not in ("positive", "negative", "both"):
            raise ValueError("polarity must be 'positive', 'negative' or 'both'")
        bad = [i for i in self.ions if i not in IONS]
        if bad:
            raise ValueError(f"unknown ion(s): {', '.join(bad)}")
        if not self.ions:
            raise ValueError("select at least one ion")
        if self.background_group not in ("starting_material", "row", "column", "off"):
            raise ValueError("background_group must be 'starting_material', 'row', 'column' or 'off'")
        if not (self.level1_raw_max_cps <= self.level2_raw_max_cps):
            raise ValueError("level1_raw_max_cps must be <= level2_raw_max_cps")
        if not (self.level3_cps <= self.cutoff_cps <= self.level5_cps):
            raise ValueError("thresholds must satisfy level3_cps <= cutoff_cps <= level5_cps")
        if self.min_spectra < 1:
            raise ValueError("min_spectra must be >= 1")
        if self.instrument not in ("unit", "high"):
            raise ValueError("instrument must be 'unit' or 'high'")
        if self.m1_rule not in ("absolute", "relative"):
            raise ValueError("m1_rule must be 'absolute' or 'relative'")
        if self.palette not in ("default", "viridis", "cividis", "okabe-ito", "grayscale", "custom"):
            raise ValueError("palette must be default, viridis, cividis, okabe-ito, grayscale or custom")
        if str(self.figure_format).lower().lstrip(".") not in ("jpg", "jpeg", "png", "tiff", "tif", "svg", "pdf"):
            raise ValueError("figure_format must be jpg, png, tiff, svg or pdf")
        self.figure_format = {"jpeg": "jpg", "tif": "tiff"}.get(str(self.figure_format).lower().lstrip("."),
                                                                 str(self.figure_format).lower().lstrip("."))
        if not 72 <= int(self.figure_dpi) <= 1200:
            raise ValueError("figure_dpi must be between 72 and 1200")
        if self.replicate_rule not in ("any", "majority", "all"):
            raise ValueError("replicate_rule must be any, majority or all")
        if self.well_shape not in ("square", "circle"):
            raise ValueError("well_shape must be square or circle")
        if str(self.transfer_plate_format) not in PLATE_FORMATS:
            raise ValueError(f"transfer_plate_format must be one of {', '.join(PLATE_FORMATS)}")
        self.transfer_plate_format = str(self.transfer_plate_format)
        if self.figure_width_mm is not None and not 40 <= float(self.figure_width_mm) <= 400:
            raise ValueError("figure_width_mm must be between 40 and 400 (or empty)")
        if not 0 <= int(self.blank_every) <= 1000:
            raise ValueError("blank_every must be >= 0")
        from .sideproducts import parse_rules
        parse_rules(self.side_products)          # raises ValueError on a malformed rule
        if not 0.5 <= float(self.font_scale) <= 2.0:
            raise ValueError("font_scale must be between 0.5 and 2.0")
        if self.spectra_pdf not in ("none", "check", "levels", "all"):
            raise ValueError("spectra_pdf must be none, check, levels or all")
        if self.mass_tolerance_ppm <= 0:
            raise ValueError("mass_tolerance_ppm must be > 0")
        if not 0 <= int(self.triage_scenario) <= 5:
            raise ValueError("triage_scenario must be between 0 and 5")
        return self

    @property
    def level4_cps(self) -> float:
        return self.cutoff_cps

    def to_dict(self) -> dict:
        return asdict(self)

    # --- standard method -------------------------------------------------------------------------------
    def non_standard(self) -> dict:
        """Detection criteria that differ from the standard method: name -> (value used, standard value)."""
        std = Settings()
        out = {}
        for name in STANDARD_CRITERIA:
            v, ref = getattr(self, name), getattr(std, name)
            same = (abs(float(v) - float(ref)) <= 1e-9 * max(1.0, abs(float(ref)))
                    if isinstance(ref, (int, float)) and not isinstance(ref, bool) and v is not None else v == ref)
            if not same:
                out[name] = (v, ref)
        return out

    @property
    def is_standard(self) -> bool:
        return not self.non_standard()

    def method_label(self) -> str:
        ns = self.non_standard()
        if not ns:
            return "standard method"
        return "NON-STANDARD PARAMETERS (" + ", ".join(f"{k} = {v[0]:g}" if isinstance(v[0], (int, float))
                                                       else f"{k} = {v[0]}" for k, v in ns.items()) + ")"

    def standard_criteria(self) -> dict:
        """The criteria of the standard method, as fixed values."""
        std = Settings()
        return {name: getattr(std, name) for name in STANDARD_CRITERIA}

    def save(self, path: Path) -> None:
        d = self.to_dict()
        d["_method"] = self.method_label()            # informative only; ignored when the file is read back
        Path(path).write_text(json.dumps(d, indent=2), encoding="utf-8")

    @classmethod
    def from_dict(cls, d: dict) -> "Settings":
        d = {k: v for k, v in d.items() if not str(k).startswith("_")}
        names = {f.name for f in fields(cls)}
        unknown = set(d) - names
        if unknown:
            raise ValueError(f"unknown setting(s): {', '.join(sorted(unknown))}")
        return cls(**d).validate()

    @classmethod
    def load(cls, path: Path) -> "Settings":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


# The detection and grading criteria of the standard method. They are fixed in the interface (only the ions
# are chosen by the user); changing any of them is possible but every output is then marked
# "NON-STANDARD PARAMETERS" with the values that differ.
STANDARD_CRITERIA = (
    "mass_tolerance_da", "mass_tolerance_ppm", "isotope_shift_tolerance_da", "min_spectra", "cutoff_cps",
    "m1_fraction_of_cutoff", "m1_rule", "m1_relative_fraction", "m1_floor_cps", "background_group",
    "background_ratio", "background_floor_cps", "background_min_distance_da", "level5_cps", "level3_cps",
    "level2_raw_max_cps", "level1_raw_max_cps", "carryover_ratio", "sm_interference_da", "window_pre_s",
    "window_post_s", "blank_start_s", "align_threshold_cps", "mz_bin_da", "hplc_product_min_pct",
    "replicate_rule", "control_min_positive_pct", "control_max_negative_pct",
)

LEVEL_LABELS = {5: "Present", 4: "Likely present", 3: "Check spectrum", 2: "Likely absent", 1: "Not detected",
                0: "Not measurable"}
