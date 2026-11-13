"""Minimal mzML reader used when pyteomics is not installed. Yields, for each <spectrum>, a dict with
the keys used by triade.runs: 'ms level', 'positive scan'/'negative scan', 'scanList' ->
{'scan': [{'scan start time': min}]}, 'm/z array', 'intensity array'. Supports 32/64-bit floats and
integers, with or without zlib compression (not MS-Numpress)."""
import base64
import zlib
import numpy as np
from lxml import etree

_NS_TAG = lambda el: etree.QName(el).localname


def _decode(bda):
    dtype, comp, kind = np.float64, False, None
    text = None
    for ch in bda:
        tag = _NS_TAG(ch)
        if tag == "cvParam":
            name = ch.get("name")
            if name == "32-bit float":
                dtype = np.float32
            elif name == "64-bit float":
                dtype = np.float64
            elif name == "32-bit integer":
                dtype = np.int32
            elif name == "64-bit integer":
                dtype = np.int64
            elif name == "zlib compression":
                comp = True
            elif name and name.startswith("MS-Numpress"):
                raise NotImplementedError("MS-Numpress compression is not supported; convert the file again without Numpress.")
            elif name in ("m/z array", "intensity array", "time array"):
                kind = name
        elif tag == "binary":
            text = ch.text or ""
    raw = base64.b64decode(text) if text else b""
    if comp and raw:
        raw = zlib.decompress(raw)
    return kind, np.frombuffer(raw, dtype=dtype).astype(np.float64)


def _time_in_minutes(p):
    v = float(p.get("value"))
    unit = (p.get("unitName") or "").lower()
    return v / 60.0 if unit.startswith("second") else v


def _spectrum_dict(el) -> dict:
    s = {"id": el.get("id"), "index": int(el.get("index", -1))}
    scan = {}
    for ch in el:
        tag = _NS_TAG(ch)
        if tag == "cvParam":
            s[ch.get("name")] = ch.get("value")
        elif tag == "scanList":
            for sc in ch.iter("{*}scan"):
                for p in sc.iter("{*}cvParam"):
                    if p.get("name") == "scan start time":
                        scan["scan start time"] = _time_in_minutes(p)
                break
        elif tag == "binaryDataArrayList":
            for bda in ch:
                kind, arr = _decode(bda)
                if kind:
                    s[kind] = arr
    s["scanList"] = {"scan": [scan]}
    s.setdefault("m/z array", np.array([]))
    s.setdefault("intensity array", np.array([]))
    return s


class MzML:
    def __init__(self, path):
        self.path = path

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def __iter__(self):
        for _, el in etree.iterparse(self.path, events=("end",), tag="{*}spectrum", huge_tree=True):
            yield _spectrum_dict(el)
            el.clear()
            while el.getprevious() is not None:
                del el.getparent()[0]


_OFFSETS: dict = {}


def read_spectrum_at(path, index: int) -> dict:
    """One spectrum by its position in the file (0-based), without parsing the rest of the file. The
    positions of the <spectrum> elements are found once per file and cached."""
    import mmap
    import os
    import re
    path = str(path)
    key = (path, os.path.getmtime(path))
    with open(path, "rb") as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
        if key not in _OFFSETS:
            _OFFSETS[key] = [m.start() for m in re.finditer(rb"<spectrum[\s>]", mm)]
        offs = _OFFSETS[key]
        if not 0 <= index < len(offs):
            raise IndexError(f"spectrum {index} not in file ({len(offs)} spectra)")
        start = offs[index]
        end = mm.find(b"</spectrum>", start) + len(b"</spectrum>")
        frag = mm[start:end]
    return _spectrum_dict(etree.fromstring(frag, parser=etree.XMLParser(huge_tree=True, recover=True)))
