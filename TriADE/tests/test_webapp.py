"""Smoke tests of the browser interface.

* test_webapp_with_stub runs the script with a minimal stand-in for Streamlit (always available), pressing
  "Load a synthetic example" and then "Run the analysis".
* test_webapp_with_streamlit does the same with Streamlit's own test runner when Streamlit is installed.
"""
import runpy
import sys
import types
from pathlib import Path

import pytest

APP = Path(__file__).resolve().parents[1] / "src" / "triade" / "webapp.py"


class _State(dict):
    __getattr__ = dict.get

    def __setattr__(self, k, v):
        self[k] = v


class _Ctx:
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def progress(self, *a, **k):
        return self

    def __getattr__(self, name):          # column.metric(...), column.download_button(...) ...
        return lambda *a, **k: None


def make_stub(pressed=(), choices=None, state=None):
    choices = choices or {}
    log = {"error": [], "warning": [], "dataframe": 0, "image": 0}
    st = types.ModuleType("streamlit")
    st.session_state = state if state is not None else _State()
    st.sidebar = _Ctx()

    def choose(label, options, index=0, **k):
        return choices.get(label, list(options)[index])

    def multiselect(label, options, default=None, **k):
        return choices.get(label, list(default or []))

    def number_input(label, min_value=None, max_value=None, value=None, step=None, **k):
        return choices.get(label, value)

    st.selectbox = st.radio = choose

    def select_slider(label, options, value=None, **k):
        return choices.get(label, value if value is not None else list(options)[0])
    st.select_slider = select_slider
    st.color_picker = lambda label, value=None, **k: choices.get(label, value)
    log["vega"] = []

    def vega_lite_chart(spec, **k):
        import json
        json.dumps(spec)                     # the spec must be plain JSON
        log["vega"].append(spec)
        return choices.get("_vega_event")
    st.vega_lite_chart = vega_lite_chart
    st.metric = lambda *a, **k: None
    st.multiselect = multiselect
    st.number_input = number_input
    st.text_input = lambda label, value="", **k: choices.get(label, value)
    st.text_area = lambda label, value="", **k: choices.get(label, value)
    st.rerun = lambda *a, **k: None
    st.checkbox = lambda label, value=False, **k: choices.get(label, value)
    st.file_uploader = lambda *a, **k: None
    st.button = lambda label, **k: label in pressed
    st.download_button = lambda *a, **k: False
    st.columns = lambda spec, **k: [_Ctx() for _ in range(spec if isinstance(spec, int) else len(spec))]
    st.tabs = lambda labels: [_Ctx() for _ in labels]
    st.expander = st.spinner = lambda *a, **k: _Ctx()
    st.progress = lambda *a, **k: _Ctx()
    st.error = lambda msg, **k: log["error"].append(str(msg))
    st.warning = lambda msg, **k: log["warning"].append(str(msg))

    def count(name):
        def f(*a, **k):
            log[name] += 1
        return f
    st.dataframe = count("dataframe")
    st.image = count("image")
    for name in ("set_page_config", "title", "header", "subheader", "markdown", "caption", "info", "success",
                 "divider", "write", "json"):
        setattr(st, name, lambda *a, **k: None)
    return st, log


def _run(monkeypatch, **kw):
    import os
    monkeypatch.setenv("HOME", os.environ.get("TRIADE_TEST_HOME", "/tmp/triade_test_home"))
    st, log = make_stub(**kw)
    monkeypatch.setitem(sys.modules, "streamlit", st)
    runpy.run_path(str(APP), run_name="__main__")
    return st, log


def test_webapp_with_stub(monkeypatch):
    choices = {"Plate format (wells)": "96"}
    st, log = _run(monkeypatch, choices=choices)                       # first page load
    assert not log["error"]
    assert any("Select clusters and adducts carefully" in w for w in log["warning"])
    state = st.session_state
    _run(monkeypatch, choices=choices, state=state, pressed={"Load a synthetic example"})
    assert state.example
    st, log = _run(monkeypatch, choices=choices, state=state, pressed={"Run the analysis"})
    assert not log["error"], log["error"]
    res = state.result
    assert res is not None and len(res.results["1"]) == 96 and res.has_hplc
    st, log = _run(monkeypatch, choices=choices, state=state)          # results, viewer, statistics
    assert log["image"] >= 6 and log["dataframe"] >= 4
    assert len(log["vega"]) >= 3                                        # plate map, trace, spectrum
    # click a well on the map (event as returned by Streamlit)
    well = state.result.results["1"]["Well"].iloc[5]
    ev = {"selection": {"well": [{"well": well}], "scan": []}}
    st, log = _run(monkeypatch, choices={**choices, "_vega_event": ev}, state=state)
    assert state.view_well == well and not log["error"]
    for mode in ("one spectrum (raw mzML)", "maximum of the window"):
        st, log = _run(monkeypatch, choices={**choices, "Show": mode}, state=state)
        assert not log["error"], log["error"]


def test_webapp_palettes_and_formats(monkeypatch):
    state = _State()
    choices = {"Plate format (wells)": "96", "Colour palette": "viridis", "File format": "svg",
               "PDF with the spectrum of each well": "all"}
    _run(monkeypatch, choices=choices, state=state, pressed={"Load a synthetic example"})
    st, log = _run(monkeypatch, choices=choices, state=state, pressed={"Run the analysis"})
    assert not log["error"], log["error"]
    res = state.result
    assert all(str(f).endswith(".svg") for _, f, _ in res.figures["1"])
    assert (res.out_dir / "7_spectra" / "Spectra_by_well.pdf").exists()


def test_webapp_high_resolution(monkeypatch):
    state = _State()
    choices = {"Plate format (wells)": "96",
               "Mass spectrometer": "QTOF / Orbitrap — high resolution (experimental)"}
    _run(monkeypatch, choices=choices, state=state, pressed={"Load a synthetic example"})
    st, log = _run(monkeypatch, choices=choices, state=state, pressed={"Run the analysis"})
    assert not log["error"], log["error"]
    assert state.result.settings.instrument == "high" and state.result.plates["1"].run.sparse


def test_webapp_extra_adduct_warning(monkeypatch):
    _, log = _run(monkeypatch, choices={"Ions / adducts / clusters": ["[M+H]+", "[M+Na]+"]})
    assert any("[M+Na]+" in e for e in log["error"])


def test_webapp_with_streamlit():
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(str(APP), default_timeout=300)
    at.run()
    assert not at.exception
    # widgets are found by their label, so the test does not depend on their order on the page
    next(sb for sb in at.selectbox if sb.label == "Plate format (wells)").set_value("96").run()
    next(b for b in at.button if b.label == "Load a synthetic example").click().run()
    assert not at.exception
    next(b for b in at.button if b.label == "Run the analysis").click().run()
    assert not at.exception, at.exception
    assert not [e.value for e in at.error], [e.value for e in at.error]
    assert at.session_state["result"] is not None
    at.run()                                   # results, spectrum viewer and statistics drawn again
    assert not at.exception, at.exception


def test_webapp_review_profiles_and_export(monkeypatch, tmp_path):
    from triade import ui
    monkeypatch.setenv("TRIADE_TEST_HOME", str(tmp_path))
    state = _State()
    choices = {"Plate format (wells)": "96", "Profile name": "my method"}
    _run(monkeypatch, choices=choices, state=state, pressed={"Load a synthetic example"})
    st, log = _run(monkeypatch, choices=choices, state=state, pressed={"Run the analysis", "Save profile"})
    assert not log["error"], log["error"]
    monkeypatch.setenv("HOME", str(tmp_path))
    assert "my method" in ui.list_profiles()
    res = state.result
    well = res.results["1"]["Well"].iloc[3]
    state.view_well = well
    st, log = _run(monkeypatch, choices={**choices, "Decision": "hplc", "Comment": "looks odd"}, state=state,
                   pressed={"Save decision"})
    assert not log["error"], log["error"]
    assert list(state.review["Well"]) == [well] and state.review["Decision"].iloc[0] == "hplc"
    # the profile is offered and loads
    st, log = _run(monkeypatch, choices={**choices, "Saved profile": "my method"}, state=state)
    assert not log["error"], log["error"]
    assert (res.out_dir / "Report.html").exists() and (res.out_dir / "Report.pdf").exists()


def test_friendly_errors():
    from triade import ui
    msg, hint = ui.friendly(ValueError("plate(s) 3 of the mzML files are not in the sample table"))
    assert "Plate column" in hint
    assert ui.friendly(ValueError("something else"))[1].startswith("Check the inputs")


def test_webapp_standard_criteria_locked(monkeypatch):
    state = _State()
    choices = {"Plate format (wells)": "96"}
    _run(monkeypatch, choices=choices, state=state, pressed={"Load a synthetic example"})
    # without the unlock, values typed in the (hidden) fields cannot change the criteria
    st, log = _run(monkeypatch, choices={**choices, "Cutoff: minimum M-peak intensity (cps) = level 4": 50000.0},
                   state=state, pressed={"Run the analysis"})
    assert not log["error"], log["error"]
    assert state.result.settings.is_standard
    # with the unlock the change is applied and flagged
    ch = {**choices, "Use non-standard criteria (requires validation)": True,
          "Cutoff: minimum M-peak intensity (cps) = level 4": 300000.0}
    st, log = _run(monkeypatch, choices=ch, state=state, pressed={"Run the analysis"})
    assert state.result.settings.cutoff_cps == 300000 and not state.result.settings.is_standard
    assert any("Non-standard criteria" in e for e in log["error"])          # the red notice in the sidebar
