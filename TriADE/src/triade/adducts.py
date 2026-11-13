"""Catalogue of ions (adducts and clusters) that can be searched for each product.

Each ion has the form  m/z = (n·M + delta) / |z|,  where M is the neutral monoisotopic mass.

Choosing ions carefully matters: every extra ion is one more m/z per well where an unrelated
signal (a starting-material adduct, a solvent cluster, another product, carryover) can pass the
criteria. On the plates used to develop this tool, adding [M+Na]+ produced more false positives
than true detections. Select only ions you expect for your chemistry, carrier solvent and
instrument, and check the result against HPLC-MS on at least one plate.
"""
from __future__ import annotations

from dataclasses import dataclass

ELECTRON = 0.000548579909
PROTON = 1.007276467
H2O = 18.010564684
NA = 22.989769282
K = 38.963706487
NH3 = 17.026549101
CL = 34.968852682
MEOH = 32.026214748
ACN = 41.026549101
DMSO = 78.013935745
FORMIC = 46.005479308
ACETIC = 60.021129372


@dataclass(frozen=True)
class Ion:
    name: str            # e.g. "[M+H]+"
    key: str             # short column prefix, e.g. "MH"
    n_m: int             # number of M in the ion (2 for dimers)
    delta: float         # mass added to n·M (Da), electrons included
    charge: int          # signed charge
    group: str           # "primary", "metal", "solvent", "dimer", "multiply charged", "loss"
    note: str            # when to use it / why it is risky

    @property
    def polarity(self) -> str:
        return "+" if self.charge > 0 else "-"

    def mz(self, mass: float) -> float:
        return (self.n_m * mass + self.delta) / abs(self.charge)

    @property
    def isotope_spacing(self) -> float:
        """Distance between M and M+1 on the m/z axis (1/|z| Da)."""
        return 1.0 / abs(self.charge)


_IONS = [
    # positive mode
    Ion("[M+H]+", "MH", 1, PROTON, 1, "primary",
        "Protonated molecule. The default and most specific choice in positive mode."),
    Ion("[M+Na]+", "MNa", 1, NA - ELECTRON, 1, "metal",
        "Sodium adduct; common in acoustic ejection (no column to remove Na+). Raises false positives: "
        "sodium adducts of starting materials and other species fall on many product masses."),
    Ion("[M+K]+", "MK", 1, K - ELECTRON, 1, "metal",
        "Potassium adduct; usually weak. Same risk as [M+Na]+."),
    Ion("[M+NH4]+", "MNH4", 1, NH3 + PROTON, 1, "metal",
        "Ammonium adduct; only expected with ammonium salts in the sample or carrier solvent."),
    Ion("[M+H-H2O]+", "MHmH2O", 1, PROTON - H2O, 1, "loss",
        "In-source water loss (alcohols, carboxylic acids). Low specificity."),
    Ion("[M+MeOH+H]+", "MMeOHH", 1, MEOH + PROTON, 1, "solvent",
        "Methanol cluster; carrier-solvent dependent."),
    Ion("[M+ACN+H]+", "MACNH", 1, ACN + PROTON, 1, "solvent",
        "Acetonitrile cluster; only with acetonitrile in the carrier solvent."),
    Ion("[M+DMSO+H]+", "MDMSOH", 1, DMSO + PROTON, 1, "solvent",
        "DMSO cluster; samples in DMSO form it, and DMSO clusters of the starting materials can land "
        "exactly on product masses."),
    Ion("[2M+H]+", "2MH", 2, PROTON, 1, "dimer",
        "Protonated dimer; appears at high concentration. Rarely needed for detection."),
    Ion("[2M+Na]+", "2MNa", 2, NA - ELECTRON, 1, "dimer",
        "Sodiated dimer; appears at high concentration."),
    Ion("[M+2H]2+", "M2H", 1, 2 * PROTON, 2, "multiply charged",
        "Doubly protonated (two basic sites, larger molecules). Its isotope spacing is 0.5 Da, which "
        "unit-resolution quadrupoles resolve poorly."),
    # negative mode
    Ion("[M-H]-", "MmH", 1, -PROTON, -1, "primary",
        "Deprotonated molecule. The default and most specific choice in negative mode."),
    Ion("[M+Cl]-", "MCl", 1, CL + ELECTRON, -1, "metal",
        "Chloride adduct; expected with chlorinated solvents or chloride salts."),
    Ion("[M+HCOO]-", "MFA", 1, FORMIC - PROTON, -1, "solvent",
        "Formate adduct; common when the carrier solvent contains formic acid."),
    Ion("[M+CH3COO]-", "MAc", 1, ACETIC - PROTON, -1, "solvent",
        "Acetate adduct; only with acetic acid or acetate buffers."),
    Ion("[M-H-H2O]-", "MmHmH2O", 1, -PROTON - H2O, -1, "loss",
        "In-source water loss. Low specificity."),
    Ion("[2M-H]-", "2MmH", 2, -PROTON, -1, "dimer",
        "Deprotonated dimer; appears at high concentration."),
]

IONS: dict[str, Ion] = {i.name: i for i in _IONS}
PRIMARY = {"+": "[M+H]+", "-": "[M-H]-"}


def ions_for(names, polarity: str = "both") -> list[Ion]:
    """The selected ions, in catalogue order, restricted to a polarity ('positive', 'negative', 'both')."""
    pol = {"positive": {"+"}, "negative": {"-"}, "both": {"+", "-"}}[polarity]
    unknown = [n for n in names if n not in IONS]
    if unknown:
        raise ValueError(f"unknown ion(s): {', '.join(unknown)}; available: {', '.join(IONS)}")
    return [i for i in _IONS if i.name in set(names) and i.polarity in pol]


def default_ions(polarity: str) -> list[str]:
    return {"positive": ["[M+H]+"], "negative": ["[M-H]-"], "both": ["[M+H]+", "[M-H]-"]}[polarity]


# Ions of the starting materials checked as possible interferences on the product ions,
# whatever the ions selected for the products (they form anyway).
SM_INTERFERENCE_IONS = {
    "+": {"[SM+H]+": PROTON, "[SM+H]+ 37Cl": PROTON + 1.99705, "[SM+Na]+": NA - ELECTRON,
          "[SM+K]+": K - ELECTRON, "[SM+NH4]+": 18.033823, "[SM+MeOH+H]+": MEOH + PROTON,
          "[SM+DMSO+H]+": DMSO + PROTON},
    "-": {"[SM-H]-": -PROTON, "[SM-H]- 37Cl": -PROTON + 1.99705, "[SM+Cl]-": CL + ELECTRON,
          "[SM+HCOO]-": FORMIC - PROTON, "[SM+CH3COO]-": ACETIC - PROTON},
}
SM_MAIN_ION = {"+": ("[SM+H]+", PROTON), "-": ("[SM-H]-", -PROTON)}
