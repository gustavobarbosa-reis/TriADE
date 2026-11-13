import pytest

from triade.adducts import IONS, ions_for
from triade.plates import get_plate_format

PROTON = 1.007276


@pytest.mark.parametrize("fmt,rows,cols,last", [("96", 8, 12, "H12"), ("384", 16, 24, "P24"), ("1536", 32, 48, "AF48")])
def test_formats(fmt, rows, cols, last):
    p = get_plate_format(fmt)
    w = p.wells()
    assert len(w) == rows * cols and w[-1] == last
    for order in p.ejection_orders().values():
        assert sorted(order) == sorted(w)


def test_parse_wells():
    p = get_plate_format("1536")
    assert p.parse("a01") == "A1"
    assert p.parse(" AF 48 ") == "AF48"
    assert p.parse("AG1") is None and p.parse("A49") is None
    assert get_plate_format("384").parse("Q1") is None
    assert p.position("AE3") == (30, 2)
    with pytest.raises(ValueError):
        get_plate_format("100")


def test_ion_mz():
    M = 300.0
    assert IONS["[M+H]+"].mz(M) == pytest.approx(M + PROTON, abs=1e-6)
    assert IONS["[M-H]-"].mz(M) == pytest.approx(M - PROTON, abs=1e-6)
    assert IONS["[M+Na]+"].mz(M) == pytest.approx(M + 22.989218, abs=1e-4)
    assert IONS["[2M+H]+"].mz(M) == pytest.approx(2 * M + PROTON, abs=1e-6)
    assert IONS["[M+2H]2+"].mz(M) == pytest.approx((M + 2 * PROTON) / 2, abs=1e-6)
    assert IONS["[M+2H]2+"].isotope_spacing == pytest.approx(0.5)


def test_ions_for_polarity():
    names = ["[M+H]+", "[M+Na]+", "[M-H]-"]
    assert [i.name for i in ions_for(names, "positive")] == ["[M+H]+", "[M+Na]+"]
    assert [i.name for i in ions_for(names, "negative")] == ["[M-H]-"]
    with pytest.raises(ValueError):
        ions_for(["[M+Xe]+"], "positive")
