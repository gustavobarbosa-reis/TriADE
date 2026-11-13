"""Shared synthetic data sets (generated once per test session; no real data is used)."""
import pandas as pd
import pytest

from triade.synthetic import make_dataset


@pytest.fixture(scope="session")
def plate96(tmp_path_factory):
    return make_dataset(tmp_path_factory.mktemp("p96"), "96", seed=0)


@pytest.fixture(scope="session")
def two_plates(tmp_path_factory):
    d = tmp_path_factory.mktemp("two")
    a = make_dataset(d, "96", seed=1, plate_id="1")
    b = make_dataset(d, "96", seed=2, plate_id="2")
    samples = d / "samples_all.csv"
    pd.concat([pd.read_csv(a["samples"]), pd.read_csv(b["samples"])]).to_csv(samples, index=False)
    return {"mzml": [a["mzml"], b["mzml"]], "samples": samples, "hplc": [a["hplc"], b["hplc"]],
            "library": pd.concat([a["library"], b["library"]])}


@pytest.fixture(scope="session")
def plate1536_subset(tmp_path_factory):
    """Four rows of a 1536-well plate (incl. AE and AF) ejected in a shuffled order given in the table."""
    return make_dataset(tmp_path_factory.mktemp("p1536"), "1536", seed=3, period_s=2.0,
                        rows=["A", "P", "AE", "AF"], injection_order=True)


@pytest.fixture(scope="session")
def negative384(tmp_path_factory):
    return make_dataset(tmp_path_factory.mktemp("neg"), "384", polarity="negative", seed=5, rows=list("ABCD"))
