import json

import pytest

from triade.chem import molecule, split_smiles_list
from triade.settings import Settings


def test_caffeine():
    m = molecule("Cn1cnc2c1c(=O)n(C)c(=O)n2C")
    assert m.formula == "C8H10N4O2"
    assert m.mass == pytest.approx(194.08038, abs=1e-4)


def test_chlorine_and_invalid():
    m = molecule("Clc1ccccn1")
    assert m.formula == "C5H4ClN"
    assert m.mass == pytest.approx(113.00323, abs=1e-4)
    assert molecule("not a smiles((") is None
    assert molecule("") is None


@pytest.mark.parametrize("smi,formula", [
    ("c1ccncc1", "C5H5N"), ("c1cc[nH]c1", "C4H5N"), ("Cn1ccnc1", "C4H6N2"), ("c1ccoc1", "C4H4O"),
    ("c1ccsc1", "C4H4S"), ("O=c1cc[nH]cc1", "C5H5NO"), ("Cc1ccc(cc1)[N+](=O)[O-]", "C7H7NO2"),
    ("c1ccc2[nH]ccc2c1", "C8H7N"), ("C1COCCN1c1ncccn1", "C8H11N3O"), ("CC(=O)N[C@@H](C)C(=O)O", "C5H9NO3")])
def test_formulas(smi, formula):
    assert molecule(smi).formula == formula


def test_radical_is_fixed():
    m = molecule("F[C](F)Oc1ccccc1")
    assert m is not None and m.note.startswith("radical")
    assert m.formula == "C7H6F2O"


def test_split_list():
    assert split_smiles_list("CCO; CCN | c1ccccc1") == ["CCO", "CCN", "c1ccccc1"]
    assert split_smiles_list(None) == []


def test_settings_roundtrip(tmp_path):
    s = Settings(plate_format="1536", polarity="negative", ions=["[M-H]-", "[M+Cl]-"], cutoff_cps=150000)
    s.validate().save(tmp_path / "s.json")
    s2 = Settings.load(tmp_path / "s.json")
    assert s2 == s
    assert json.loads((tmp_path / "s.json").read_text())["ions"] == ["[M-H]-", "[M+Cl]-"]


@pytest.mark.parametrize("bad", [dict(plate_format="100"), dict(polarity="x"), dict(ions=[]),
                                 dict(ions=["[M+Xe]+"]), dict(level3_cps=5e5), dict(min_spectra=0),
                                 dict(background_group="plate"), dict(triage_scenario=9)])
def test_settings_validation(bad):
    with pytest.raises(ValueError):
        Settings(**bad).validate()


def test_unknown_setting():
    with pytest.raises(ValueError):
        Settings.from_dict({"cutof": 1})
