"""Formula and masses from SMILES (RDKit when installed, otherwise the built-in reader)."""
from __future__ import annotations

from dataclasses import dataclass

try:  # pragma: no cover - depends on the environment
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    HAVE_RDKIT = True
except ImportError:  # pragma: no cover
    HAVE_RDKIT = False

from . import _smiles


@dataclass
class Molecule:
    smiles: str
    formula: str
    mass: float            # monoisotopic, neutral
    average_mass: float
    note: str = ""         # e.g. "radical atom corrected"


def _missing(smiles) -> bool:
    return smiles is None or str(smiles).strip().lower() in ("", "nan", "none", "na")


def molecule(smiles) -> Molecule | None:
    """Reads a SMILES. Atoms written as radicals (e.g. '[C]' in an -OCHF2 exported from some drawing
    programs) get the missing hydrogens, with a note. Returns None for a missing or invalid SMILES."""
    if _missing(smiles):
        return None
    smi = str(smiles).strip()
    fixed = []
    if HAVE_RDKIT:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        for a in m.GetAtoms():
            if a.GetNumRadicalElectrons():
                fixed.append(f"{a.GetSymbol()}{a.GetIdx()}")
                a.SetNumRadicalElectrons(0)
                a.SetNoImplicit(False)
        if fixed:
            Chem.SanitizeMol(m)
        mol = Molecule(smi, rdMolDescriptors.CalcMolFormula(m), Descriptors.ExactMolWt(m), Descriptors.MolWt(m))
    else:
        m = _smiles.parse(smi)
        if m is None:
            return None
        for a in m.GetAtoms():
            if a.GetNumRadicalElectrons():
                fixed.append(f"{a.GetSymbol()}{a.GetIdx()}")
                a.SetNumRadicalElectrons(0)
                a.SetNoImplicit(False)
        mol = Molecule(smi, _smiles.formula(m), _smiles.exact_mass(m), _smiles.average_mass(m))
    if fixed:
        mol.note = "radical atom in SMILES corrected by adding H (" + ", ".join(fixed) + ")"
    return mol


def split_smiles_list(text) -> list[str]:
    """Several SMILES in one cell, separated by ';' or '|' (a '.' inside a SMILES is kept: it is a salt
    or a mixture written as one structure)."""
    if _missing(text):
        return []
    return [s.strip() for s in str(text).replace("|", ";").split(";") if s.strip()]
