"""Plate maps and charts. Every figure takes its colours from a palette (default, viridis, cividis,
Okabe-Ito, grayscale or custom) and is saved in the chosen format (jpg, png, tiff, svg, pdf)."""
from __future__ import annotations

import textwrap
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import colormaps  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, to_rgb  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

from .plates import PlateFormat  # noqa: E402
from .settings import LEVEL_LABELS  # noqa: E402

INK = "#2b2b29"
MUTED = "#6b6a66"
GRID = "#e4e3df"
CELL_EDGE = "#000000"      # cell border colour
CELL_LW = 0.35             # cell border width (pt)

# colour roles used by the figures; a custom palette may override any of them
ROLES = ["level5", "level4", "level3", "level2", "level1", "level0", "detected", "not_detected", "check",
         "both", "neither", "echo_only", "ref_only", "seq", "line", "line2", "band", "empty"]

PALETTES = {
    "default": {
        "level5": "#184f95", "level4": "#86b6ef", "level3": "#e4e3df", "level2": "#eb9992", "level1": "#8a2b2a",
        "level0": "#ffffff", "detected": "#184f95", "not_detected": "#e4e3df", "check": "#e8b04b",
        "both": "#184f95", "neither": "#e4e3df", "echo_only": "#eb6834", "ref_only": "#9b59d0",
        "seq": ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"],
        "line": "#184f95", "line2": "#eb6834", "band": "#cde2fb", "empty": "#ffffff"},
    "viridis": {
        "level5": "#fde725", "level4": "#5ec962", "level3": "#21918c", "level2": "#3b528b", "level1": "#440154",
        "level0": "#ffffff", "detected": "#fde725", "not_detected": "#440154", "check": "#21918c",
        "both": "#fde725", "neither": "#440154", "echo_only": "#5ec962", "ref_only": "#3b528b",
        "seq": "viridis", "line": "#31688e", "line2": "#35b779", "band": "#e3f1ec", "empty": "#ffffff"},
    "cividis": {
        "level5": "#fee838", "level4": "#bcaf6f", "level3": "#7d7c78", "level2": "#434e6c", "level1": "#00224e",
        "level0": "#ffffff", "detected": "#fee838", "not_detected": "#00224e", "check": "#7d7c78",
        "both": "#fee838", "neither": "#00224e", "echo_only": "#bcaf6f", "ref_only": "#434e6c",
        "seq": "cividis", "line": "#00224e", "line2": "#bcaf6f", "band": "#ecebe4", "empty": "#ffffff"},
    "okabe-ito": {
        "level5": "#0072b2", "level4": "#56b4e9", "level3": "#f0e442", "level2": "#e69f00", "level1": "#d55e00",
        "level0": "#ffffff", "detected": "#0072b2", "not_detected": "#e0e0e0", "check": "#f0e442",
        "both": "#0072b2", "neither": "#e0e0e0", "echo_only": "#d55e00", "ref_only": "#cc79a7",
        "seq": ["#e8f3fa", "#9fd1f0", "#56b4e9", "#2a8fcb", "#0072b2", "#004a75"],
        "line": "#0072b2", "line2": "#d55e00", "band": "#e8f3fa", "empty": "#ffffff"},
    "grayscale": {
        "level5": "#1a1a1a", "level4": "#555555", "level3": "#8c8c8c", "level2": "#c4c4c4", "level1": "#f0f0f0",
        "level0": "#ffffff", "detected": "#1a1a1a", "not_detected": "#f0f0f0", "check": "#8c8c8c",
        "both": "#1a1a1a", "neither": "#f0f0f0", "echo_only": "#6e6e6e", "ref_only": "#b5b5b5",
        "seq": "Greys", "line": "#333333", "line2": "#8c8c8c", "band": "#eeeeee", "empty": "#ffffff"},
}
PALETTE_NAMES = {"default": "Default (blue / red)", "viridis": "Viridis", "cividis": "Cividis",
                 "okabe-ito": "Okabe-Ito (colour-blind safe)", "grayscale": "Grayscale (print)", "custom": "Custom"}


FONT = {"scale": 1.0, "base": 1.0, "width_in": None, "titles": True}   # state of the figure being drawn
CAPTIONS: list = []        # (file stem, title, subtitle) of figures drawn without titles (publication mode)


def _fs(size: float) -> float:
    return size * FONT["scale"]


def _grow(w: float, h: float) -> tuple:
    """Chart size. Publication mode: the chosen width (text shrinks with it, down to ~55 %). Otherwise a
    little larger when the text is enlarged, so labels do not collide."""
    if FONT["width_in"]:
        f = FONT["width_in"] / w
        FONT["scale"] = FONT["base"] * min(1.0, max(0.55, f))
        return (FONT["width_in"], h * f)
    FONT["scale"] = FONT["base"]
    g = max(1.0, 1 + 0.6 * (FONT["base"] - 1))
    return (w * g, h * g)


def _wrap(text: str, width: int) -> str:
    """Wraps a title to the space it has at the current text size."""
    return "\n".join(textwrap.wrap(text, max(20, int(width / max(FONT["scale"], 1e-6))), break_long_words=False))


def _tight(fig, **kw):
    """tight_layout without the warning matplotlib gives when a very narrow figure (publication mode) cannot
    fit every label: the figure is still saved, with the default spacing."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout(**kw)


def _suptitle(fig, text, stem=None, **kw):
    if FONT["titles"]:
        fig.suptitle(text, **kw)
    else:
        CAPTIONS.append((str(stem or ""), text.replace("\n", " "), ""))


@dataclass
class Style:
    colors: dict = field(default_factory=lambda: dict(PALETTES["default"]))
    fmt: str = "jpg"
    dpi: int = 300
    cell_values: bool = True
    font_scale: float = 1.0
    width_mm: float | None = None        # publication mode: figure width
    font_family: str = "default"
    titles: bool = True
    well_shape: str = "square"
    highlight: bool = True

    def __post_init__(self):
        FONT.update(scale=float(self.font_scale), base=float(self.font_scale),
                    width_in=float(self.width_mm) / 25.4 if self.width_mm else None, titles=bool(self.titles))
        plt.rcParams["font.size"] = 10 * float(self.font_scale)
        fam = {"default": None, "Arial": "Arial", "Helvetica": "Helvetica", "Times New Roman": "Times New Roman",
               "DejaVu Sans": "DejaVu Sans"}.get(self.font_family)
        plt.rcParams["font.family"] = "sans-serif" if fam in (None, "Arial", "Helvetica", "DejaVu Sans") else "serif"
        if fam:
            key = "font.sans-serif" if plt.rcParams["font.family"] == ["sans-serif"] else "font.serif"
            plt.rcParams[key] = [fam] + [f for f in plt.rcParams[key] if f != fam]

    @classmethod
    def from_settings(cls, s) -> "Style":
        base = "default" if s.palette == "custom" else s.palette
        colors = dict(PALETTES[base])
        if s.palette == "custom":
            colors.update({k: v for k, v in (s.custom_colors or {}).items() if k in ROLES and v})
        return cls(colors, s.figure_format, int(s.figure_dpi), bool(s.cell_values), float(s.font_scale),
                   s.figure_width_mm, s.font_family, bool(s.show_titles), s.well_shape, bool(s.highlight_wells))

    def __getitem__(self, role):
        return self.colors[role]

    @property
    def cmap(self):
        seq = self.colors["seq"]
        if isinstance(seq, str):
            return colormaps[seq] if seq != "Greys" else LinearSegmentedColormap.from_list(
                "greys", ["#f0f0f0", "#1a1a1a"])
        return LinearSegmentedColormap.from_list("seq", list(seq))

    def level(self, k: int) -> str:
        return self.colors[f"level{int(k)}"]

    def activate(self):
        """Makes this style's font size the current one (when several styles are in use)."""
        self.__post_init__()
        return self

    def save(self, fig, stem, tight=True) -> Path:
        """Saves to stem + extension of the chosen format."""
        path = Path(f"{stem}.{self.fmt}")
        kw = {"dpi": self.dpi, "facecolor": "white"}
        if tight:
            kw.update(bbox_inches="tight", pad_inches=0.25)
        if self.fmt == "jpg":
            kw["pil_kwargs"] = {"quality": 92}
        elif self.fmt == "tiff":
            kw["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(path, format={"jpg": "jpeg"}.get(self.fmt, self.fmt), **kw)
        plt.close(fig)
        return path


def text_color(face) -> str:
    """Black or white text, whichever reads better on the cell colour."""
    r, g, b = to_rgb(face)
    lum = 0.2126 * r ** 2.2 + 0.7152 * g ** 2.2 + 0.0722 * b ** 2.2
    return INK if lum > 0.25 else "white"


def _swatch(color, label):
    return Patch(facecolor=color, edgecolor=CELL_EDGE, lw=CELL_LW, label=label)


class PlateFigure:
    """A plate map: one square (or circle) per well with a thin black border. With ax given, it draws into an
    existing axes (composite figures) and finish_axes() completes it without saving."""

    def __init__(self, plate: PlateFormat, style: Style, colorbar: bool = False, ax=None, fig=None):
        self.plate, self.style = plate, style
        big = plate.n_cols > 24
        width = min(18.0, max(9.5, 0.5 * plate.n_cols + 1.5)) + (1.0 if colorbar else 0)
        height = width * plate.n_rows / plate.n_cols * 0.86 + 2.4
        if FONT["width_in"]:
            f = FONT["width_in"] / width
            FONT["scale"] = FONT["base"] * min(1.0, max(0.55, f))
            width, height = FONT["width_in"], height * f
        else:
            FONT["scale"] = FONT["base"]
        if ax is None:
            self.fig = plt.figure(figsize=(width, height))
            self.ax = self.fig.add_axes([0.05, 0.08, 0.9, 0.72])
            self.ax.set_anchor("N")               # the plate sits right under the title
        else:
            self.fig, self.ax = fig, ax
        self.glyphs = plate.n_cols <= 24          # text inside the cells only up to 384 wells
        self.font = (6 if big else 9) * FONT["scale"]
        self.drawn = set()
        self.marks = {}

    def _patch(self, i, j, face, edge=CELL_EDGE, lw=CELL_LW, ls="-", z=1):
        if self.style.well_shape == "circle":
            from matplotlib.patches import Circle
            return Circle((j + 0.5, i + 0.5), 0.43, facecolor=face, edgecolor=edge, lw=lw, ls=ls, zorder=z)
        return Rectangle((j, i), 1, 1, facecolor=face, edgecolor=edge, lw=lw, ls=ls, zorder=z)

    def cell(self, well, color, text=None, text_color_=None, small=False):
        i, j = self.plate.position(well)
        self.drawn.add(well)
        self.ax.add_patch(self._patch(i, j, color))
        if text and self.glyphs:
            self.ax.text(j + 0.5, i + 0.53, text, ha="center", va="center",
                         fontsize=(self.font - 1.5 * FONT["scale"]) if small else (self.font + 0.5 * FONT["scale"]),
                         fontweight="normal" if small else "bold", color=text_color_ or text_color(color), zorder=3)

    def highlights(self, res, discrepancies: bool = True):
        """Outlines: discrepancy with HPLC-MS (thick, unless discrepancies=False), reviewed (dashed), control
        (corner letter), QC flag (dot)."""
        if not self.style.highlight:
            return
        get = lambda c: res[c].fillna("").astype(str) if c in res.columns else pd.Series("", index=res.index)
        disc = get("Agreement").str.endswith(" only") & discrepancies
        for w, d, rv, c, q in zip(res["Well"], disc, get("Review_Decision"), get("Control"), get("Run_QC_Flag")):
            i, j = self.plate.position(w)
            if d:
                self.ax.add_patch(self._patch(i, j, "none", edge=INK, lw=2.2, z=4))
                self.marks["disc"] = True
            if rv:
                if self.style.well_shape == "square":
                    self.ax.add_patch(Rectangle((j + 0.09, i + 0.09), 0.82, 0.82, facecolor="none",
                                                edgecolor=self.style["line2"], lw=1.3, ls=(0, (2, 1.5)), zorder=4))
                else:
                    from matplotlib.patches import Circle
                    self.ax.add_patch(Circle((j + 0.5, i + 0.5), 0.36, facecolor="none", edgecolor=self.style["line2"],
                                             lw=1.3, ls=(0, (2, 1.5)), zorder=4))
                self.marks["review"] = True
            if c and self.glyphs:
                self.ax.text(j + 0.1, i + 0.12, {"positive": "P", "negative": "N", "blank": "B"}[c], ha="left",
                             va="top", fontsize=self.font * 0.6, color=INK, fontweight="bold", zorder=5,
                             bbox={"boxstyle": "square,pad=0.08", "fc": "white", "ec": "none", "alpha": 0.85})
                self.marks["control"] = True
            if q:
                self.ax.plot(j + 0.85, i + 0.15, marker="o", ms=2.6 * FONT["scale"], color=self.style["line2"],
                             mec="white", mew=0.4, zorder=5)
                self.marks["qc"] = True

    def mark_handles(self):
        from matplotlib.lines import Line2D
        h = []
        if self.marks.get("disc"):
            h.append(Patch(facecolor="none", edgecolor=INK, lw=2.2, label="differs from HPLC-MS"))
        if self.marks.get("review"):
            h.append(Patch(facecolor="none", edgecolor=self.style["line2"], lw=1.3, ls="--", label="reviewed by hand"))
        if self.marks.get("control"):
            h.append(Patch(facecolor="white", edgecolor="#b0afaa", lw=0.5, label="P / N / B control"))
        if self.marks.get("qc"):
            h.append(Line2D([], [], marker="o", ls="", color=self.style["line2"], label="run QC flag"))
        return h

    def finish_axes(self, handles=None, ncol=None, row_labels=None, col_labels=None):
        ax, p = self.ax, self.plate
        empty = [w for w in p.wells() if w not in self.drawn]
        for w in empty:          # wells that are not in the sample table
            i, j = p.position(w)
            ax.add_patch(self._patch(i, j, self.style["empty"], edge="#d6d5d0"))
        handles = list(handles or [])
        if empty and handles:
            handles.append(Patch(facecolor=self.style["empty"], edgecolor="#b0afaa", lw=CELL_LW,
                                 label=f"not in the sample table ({len(empty)})"))
        if handles:
            handles += self.mark_handles()
        ax.set_xlim(0, p.n_cols); ax.set_ylim(p.n_rows, 0)
        ax.set_xticks(np.arange(p.n_cols) + 0.5, [str(c) for c in range(1, p.n_cols + 1)], fontsize=self.font,
                      color=INK)
        ylab = [f"{r} · {row_labels[r]}" if row_labels and row_labels.get(r) else r for r in p.rows]
        ax.set_yticks(np.arange(p.n_rows) + 0.5, ylab, fontsize=self.font, color=INK)
        ax.xaxis.tick_top(); ax.tick_params(length=0)
        if col_labels:
            sec = ax.secondary_xaxis("bottom")
            sec.set_xticks(np.arange(p.n_cols) + 0.5, [col_labels.get(c, "") for c in range(1, p.n_cols + 1)],
                           fontsize=self.font * 0.85, color=MUTED, rotation=90)
            sec.tick_params(length=0)
            for sp in sec.spines.values():
                sp.set_visible(False)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_aspect("equal")
        if handles:
            ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.01 - (0.09 if col_labels else 0)),
                      ncol=ncol or (min(len(handles), 5) if p.n_cols > 12 else min(3, len(handles))),
                      frameon=False, fontsize=_fs(10))

    def finish(self, stem, title, subtitle, handles=None, ncol=None, row_labels=None, col_labels=None) -> Path:
        self.finish_axes(handles, ncol, row_labels, col_labels)
        if FONT["titles"]:
            self.fig.text(0.05, 0.965, title, ha="left", va="top", fontsize=_fs(14), fontweight="bold", color=INK)
            self.fig.text(0.05, 0.92, subtitle, ha="left", va="top", fontsize=_fs(9.5), color=MUTED, linespacing=1.5)
        else:
            CAPTIONS.append((str(stem), title, subtitle))
        return self.style.save(self.fig, stem)


def _colorbar(f: PlateFigure, cmap, norm, label):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    cax = make_axes_locatable(f.ax).append_axes("right", size="1.6%", pad=0.15)
    cb = f.fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax)
    cb.set_label(label, color=INK, fontsize=_fs(10))
    cb.outline.set_visible(False); cb.ax.tick_params(colors=MUTED, labelsize=_fs(9))


# ------------------------------------------------------------------------------------------
# Plate maps
# ------------------------------------------------------------------------------------------
def _lab(labels):
    return {"row_labels": labels[0], "col_labels": labels[1]} if labels else {}


def heatmap_log10(res, plate, stem, title, subtitle, vmin, vmax, style: Style, labels=None) -> Path:
    cmap, norm = style.cmap, plt.Normalize(vmin, vmax)
    f = PlateFigure(plate, style, colorbar=True)
    for w, v in zip(res["Well"], res["Log10_Max_Validated_Intensity"]):
        if np.isfinite(v):
            face = cmap(norm(v))
            f.cell(w, face, f"{v:.1f}" if style.cell_values else None, small=True)
        else:
            f.cell(w, style["not_detected"])
    _colorbar(f, cmap, norm, "log10 intensity (cps)")
    return f.finish(stem, title, subtitle, **_lab(labels), handles=[_swatch(style["not_detected"], "not detected (no validated signal)")])


def heatmap_binary(res, plate, stem, title, subtitle, style: Style, labels=None) -> Path:
    f = PlateFigure(plate, style)
    for w, d in zip(res["Well"], res["Detected"]):
        f.cell(w, style["detected"] if d == 1 else style["not_detected"])
    n_det = int(res["Detected"].sum())
    return f.finish(stem, title, subtitle, **_lab(labels), handles=[_swatch(style["detected"], f"Detected ({n_det})"),
                                            _swatch(style["not_detected"], f"Not detected ({len(res) - n_det})")])


def heatmap_levels(res, plate, stem, title, subtitle, style: Style, labels=None) -> Path:
    f = PlateFigure(plate, style)
    for w, lv in zip(res["Well"], res["Confidence_Level"]):
        lv = int(lv)
        f.cell(w, style.level(lv), "×" if lv == 0 else str(lv))
    f.highlights(res)
    cnt = res["Confidence_Level"].value_counts()
    ks = [5, 4, 3, 2, 1] + ([0] if cnt.get(0, 0) else [])
    return f.finish(stem, title, subtitle, **_lab(labels), handles=[_swatch(style.level(k), f"{k if k else '×'} {LEVEL_LABELS[k]} "
                                                                     f"({int(cnt.get(k, 0))})") for k in ks])


def heatmap_call(res, plate, stem, title, subtitle, style: Style, labels=None) -> Path:
    colors = {"DETECTED": style["detected"], "CHECK SPECTRUM": style["check"], "NOT DETECTED": style["not_detected"],
              "NOT MEASURABLE": style.level(0)}
    lab = {"DETECTED": "Detected", "CHECK SPECTRUM": "Check spectrum", "NOT DETECTED": "Not detected",
           "NOT MEASURABLE": "Not measurable (×)"}
    f = PlateFigure(plate, style)
    for w, c in zip(res["Well"], res["Call"]):        # the colour alone marks 'Check spectrum'
        f.cell(w, colors[c], "×" if c == "NOT MEASURABLE" else None)
    f.highlights(res, discrepancies=False)
    cnt = res["Call"].value_counts()
    keys = [k for k in colors if k != "NOT MEASURABLE" or cnt.get(k, 0)]
    return f.finish(stem, title, subtitle, **_lab(labels), handles=[_swatch(colors[k], f"{lab[k]}: {int(cnt.get(k, 0))}") for k in keys])


def heatmap_agreement(res, plate, stem, title, subtitle, ref_label, style: Style, labels=None) -> Path:
    classes = {"both": ("Both detected", style["both"], ""), "neither": ("Neither", style["neither"], ""),
               "Echo-MS only": ("Echo-MS only", style["echo_only"], "+"),
               f"{ref_label} only": (f"{ref_label} only", style["ref_only"], "−")}
    f = PlateFigure(plate, style)
    for w, a in zip(res["Well"], res["Agreement"]):
        _, color, glyph = classes[a]
        f.cell(w, color, glyph or None)
    cnt = res["Agreement"].value_counts()
    return f.finish(stem, title, subtitle, **_lab(labels), handles=[_swatch(c, f"{lbl}{' (' + g + ')' if g else ''}: {int(cnt.get(k, 0))}")
                                            for k, (lbl, c, g) in classes.items()])


def heatmap_conversion(res, plate, stem, title, subtitle, ref_label, style: Style, labels=None) -> Path:
    """The HPLC-MS conversion of each well (0-100 %), with the value written in the cell."""
    cmap, norm = style.cmap, plt.Normalize(0, 100)
    f = PlateFigure(plate, style, colorbar=True)
    for w, v, ns, ind in zip(res["Well"], res["HPLC_Conversion"], res.get("HPLC_No_Signal", 0 * res["Detected"]),
                             res.get("HPLC_Indeterminate", 0 * res["Detected"])):
        if np.isfinite(v):
            face = cmap(norm(v))
            f.cell(w, face, f"{v:.0f}" if style.cell_values else None, small=True)
        else:
            f.cell(w, style["empty"], "n/a" if style.cell_values else None, text_color_=MUTED, small=True)
    _colorbar(f, cmap, norm, f"{ref_label} conversion (%)")
    n_na = int(res["HPLC_Conversion"].isna().sum())
    return f.finish(stem, title, subtitle, **_lab(labels), handles=[_swatch(style["empty"], f"no {ref_label} result (no signal, "
                                                                    f"indeterminate or missing): {n_na}")] if n_na else [])


def composite_levels_conversion(res, plate, stem, title, ref_label, style: Style, labels=None) -> Path:
    """Panel a: Echo-MS confidence levels; panel b: HPLC-MS conversion — same plate, same layout."""
    stacked = plate.n_cols > 24
    one = min(18.0, max(9.5, 0.5 * plate.n_cols + 1.5))
    w, h = (one + 1.2, 2 * (one * plate.n_rows / plate.n_cols * 0.86 + 1.6)) if stacked else \
        (2 * one + 1.6, one * plate.n_rows / plate.n_cols * 0.86 + 2.6)
    size = _grow(w, h)
    fig = plt.figure(figsize=size)
    if stacked:
        ax1, ax2 = fig.add_axes([0.06, 0.55, 0.86, 0.38]), fig.add_axes([0.06, 0.08, 0.86, 0.38])
    else:
        ax1, ax2 = fig.add_axes([0.04, 0.12, 0.43, 0.72]), fig.add_axes([0.53, 0.12, 0.43, 0.72])
    scale = FONT["scale"]
    f1 = PlateFigure(plate, style, ax=ax1, fig=fig)
    for wl, lv in zip(res["Well"], res["Confidence_Level"]):
        lv = int(lv)
        f1.cell(wl, style.level(lv), "×" if lv == 0 else str(lv))
    f1.highlights(res)
    cnt = res["Confidence_Level"].value_counts()
    ks = [5, 4, 3, 2, 1] + ([0] if cnt.get(0, 0) else [])
    f1.finish_axes([_swatch(style.level(k), f"{k if k else '×'} {LEVEL_LABELS[k]} ({int(cnt.get(k, 0))})") for k in ks],
                   ncol=3, **_lab(labels))
    FONT["scale"] = scale
    f2 = PlateFigure(plate, style, ax=ax2, fig=fig)
    cmap, norm = style.cmap, plt.Normalize(0, 100)
    for wl, v in zip(res["Well"], res["HPLC_Conversion"]):
        if np.isfinite(v):
            f2.cell(wl, cmap(norm(v)), f"{v:.0f}" if style.cell_values else None, small=True)
        else:
            f2.cell(wl, style["empty"], "n/a" if style.cell_values else None, text_color_=MUTED, small=True)
    f2.finish_axes([_swatch(style["empty"], f"no {ref_label} result")] if res["HPLC_Conversion"].isna().any() else [],
                   **_lab(labels))
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    cax = make_axes_locatable(ax2).append_axes("right", size="2.5%", pad=0.1)
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax)
    cb.set_label(f"{ref_label} conversion (%)", fontsize=_fs(9)); cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=_fs(8), colors=MUTED)
    for ax, letter, sub in ((ax1, "a", "Echo-MS confidence level"), (ax2, "b", f"{ref_label} conversion (%)")):
        ax.text(-0.06, 1.12, letter, transform=ax.transAxes, fontsize=_fs(16), fontweight="bold", color=INK,
                va="bottom")
        ax.text(0.0, 1.12, sub, transform=ax.transAxes, fontsize=_fs(10.5), color=INK, va="bottom")
    if FONT["titles"]:
        fig.text(0.04, 0.985, title, ha="left", va="top", fontsize=_fs(13), fontweight="bold", color=INK)
    else:
        CAPTIONS.append((str(stem), title, "a: Echo-MS confidence level; b: " + f"{ref_label} conversion (%)"))
    return style.save(fig, stem)


# ------------------------------------------------------------------------------------------
# Charts
# ------------------------------------------------------------------------------------------
def _clean(ax):
    ax.tick_params(colors=MUTED)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def levels_vs_reference_chart(res, stem, plate_id, s, style: Style, ci: dict | None = None) -> Path:
    ref = s.reference_label
    x = res.copy()
    x["HPLC_Conversion"] = x["HPLC_Conversion"].fillna(0)
    levels = ([0] if (x["Confidence_Level"] == 0).any() else []) + [1, 2, 3, 4, 5]
    pct = [100 * x.loc[x["Confidence_Level"] == k, "HPLC_Product"].mean() if (x["Confidence_Level"] == k).any() else 0.0
           for k in levels]
    ns = [int((x["Confidence_Level"] == k).sum()) for k in levels]
    names = [f"{k if k else '×'}\n" + "\n".join(textwrap.wrap(LEVEL_LABELS[k], 12 if FONT["scale"] <= 1.1 else 8, break_long_words=False))
             for k in levels]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=_grow(13, 5.4), gridspec_kw={"wspace": 0.28})
    xs = range(len(levels))
    a1.bar(xs, pct, width=0.62, color=[style.level(k) for k in levels], edgecolor="#b0afaa", lw=0.5)
    for i, (k, p, n) in enumerate(zip(levels, pct, ns)):
        top = p
        if ci and k in ci and n:
            lo, hi = ci[k]
            a1.plot([i, i], [lo, hi], color=INK, lw=1.2)
            a1.plot([i - 0.08, i + 0.08], [lo, lo], color=INK, lw=1.2)
            a1.plot([i - 0.08, i + 0.08], [hi, hi], color=INK, lw=1.2)
            top = max(p, hi)
        a1.text(i, top + 2, f"{p:.0f}%\n(n = {n})", ha="center", va="bottom", fontsize=_fs(9.5), color=INK)
    a1.set_ylim(0, 122); a1.set_yticks(range(0, 101, 20))
    a1.set_ylabel(f"wells with product in {ref} (%)\nconversion ≥ {s.hplc_product_min_pct:g}%", color=INK)
    a1.set_title(_wrap(f"Share of wells with the product in {ref}" + (" (bars = 95% CI)" if ci else ""), 55),
                 loc="left", fontsize=_fs(11), color=INK)
    rng = np.random.default_rng(0)
    for i, k in enumerate(levels):
        y = x.loc[x["Confidence_Level"] == k, "HPLC_Conversion"].to_numpy()
        a2.scatter(i + rng.uniform(-0.22, 0.22, len(y)), y, s=14, color=style["line"], alpha=0.45, lw=0)
        if len(y):
            a2.plot([i - 0.3, i + 0.3], [np.median(y)] * 2, color=INK, lw=2)
    a2.set_ylim(-3, 103); a2.set_ylabel(f"{ref} conversion (%)", color=INK)
    a2.set_title(_wrap(f"{ref} conversion of each well (line = median)", 55), loc="left", fontsize=_fs(11),
                 color=INK)
    for ax in (a1, a2):
        ax.set_xticks(list(xs), names, fontsize=_fs(9.5), color=INK)
        ax.set_xlabel("confidence level (Echo-MS)", color=INK)
        _clean(ax)
        ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    fig.subplots_adjust(top=0.84 - 0.08 * max(0.0, FONT["scale"] - 1))
    i5, i1 = levels.index(5), levels.index(1)
    _suptitle(fig, _wrap(f"Plate {plate_id} — confidence levels vs. {ref} ({len(x)} wells): level 5 = {pct[i5]:.0f}% "
                       f"with product, level 1 = {pct[i1]:.0f}%", 120), x=0.06, y=0.99, ha="left", fontsize=_fs(13),
                 fontweight="bold", color=INK, stem=stem)
    n_ns = int(res["HPLC_No_Signal"].fillna(0).sum()) if "HPLC_No_Signal" in res else 0
    n_ind = int(res["HPLC_Indeterminate"].fillna(0).sum()) if "HPLC_Indeterminate" in res else 0
    if n_ns or n_ind:
        txt = f"{n_ns} wells with no {ref} signal" + (f" and {n_ind} marked as indeterminate" if n_ind else "")
        fig.text(0.06, 0.915, txt + " are counted as not detected (no product).", ha="left", va="bottom",
                 fontsize=_fs(10), color=MUTED)
    return style.save(fig, stem)


def qc_tic(pa, stem, s, style: Style) -> Path:
    run, grid, slot_of = pa.run, pa.grid, pa.slot_of
    first_slot, last_slot = min(slot_of.values()), max(slot_of.values())
    last_slot = min(last_slot, len(grid) - 1)
    fig, axes = plt.subplots(3, 1, figsize=_grow(14, 10))
    ax = axes[0]
    ax.plot(run.t, run.tic, lw=0.4, color=style["line"])
    ax.axvspan(grid[first_slot] - s.window_pre_s, grid[last_slot] + s.window_post_s, color=style["band"],
               alpha=0.5, lw=0)
    ax.set_yscale("log"); ax.set_xlabel("time (s)"); ax.set_ylabel("TIC (cps)")
    ax.set_title(f"Plate {pa.plate_id} — full TIC; shaded = {len(slot_of)} assigned ejections "
                 f"(period {pa.period:.4f} s)", loc="left", fontsize=_fs(11))
    by_slot = {sl: w for w, sl in slot_of.items()}
    for ax, (t_a, t_b) in zip(axes[1:], [(grid[first_slot] - 12, grid[first_slot] + 50),
                                        (grid[last_slot] - 50, grid[last_slot] + 12)]):
        m = (run.t >= t_a) & (run.t <= t_b)
        if not m.any():
            continue
        ax.plot(run.t[m], run.tic[m], ".-", lw=0.7, ms=2.5, color=style["line"])
        for sl, g in enumerate(grid):
            if t_a <= g <= t_b:
                own = sl in by_slot
                ax.axvspan(g - s.window_pre_s, g + s.window_post_s, color=style["line"] if own else "#b0afaa",
                           alpha=0.15 if own else 0.08, lw=0)
                if own:
                    ax.text(g, run.tic[m].max() * 1.02, by_slot[sl], ha="center", va="bottom", fontsize=_fs(7), color=INK)
        ax.set_xlim(t_a, t_b); ax.set_xlabel("time (s)"); ax.set_ylabel("TIC (cps)")
        ax.set_ylim(top=run.tic[m].max() * 1.12)
    axes[1].set_title(f"Start of run — ejection order: {pa.order_name}; band = spectra used for each well",
                      loc="left", fontsize=_fs(11))
    axes[2].set_title("End of run", loc="left", fontsize=_fs(11))
    _tight(fig)
    return style.save(fig, stem, tight=False)


def qc_profile(pa, stem, s, style: Style) -> Path:
    """Mean ejection profile: TIC / baseline and the normalised XIC of the expected primary ion, aligned on
    the pulse apex, with the ejection and blank windows."""
    from .runs import ion_window_max
    run, grid = pa.run, pa.grid
    half = min(1.25, pa.period)
    tt = np.arange(-half, half, 0.02)
    res = pa.results
    mz_col = next((c for c in ("MH_Theoretical_mz", "MmH_Theoretical_mz") if c in res.columns), None)
    target = dict(zip(res["Well"], res[mz_col])) if mz_col else {}
    tic_prof, xic_prof = [], []
    for w, sl in pa.slot_of.items():
        if sl >= len(grid):
            continue
        g = grid[sl]
        tic_prof.append(np.interp(g + tt, run.t, pa.rel_tic + 1.0))
        mz = target.get(w, np.nan)
        if np.isfinite(mz):
            x = ion_window_max(run, mz, s.mass_tolerance_da)
            local = x[(run.t > g - half) & (run.t < g + half)]
            if local.size and local.max() > 0:
                xic_prof.append(np.interp(g + tt, run.t, x / local.max()))
    fig, ax = plt.subplots(figsize=_grow(8, 4.5))
    ax.axvspan(-s.window_pre_s, s.window_post_s, color=style["band"], alpha=0.7, lw=0,
               label="ejection window (spectra analysed)")
    ax.axvspan(pa.blank[0], min(pa.blank[1], half), color=GRID, alpha=0.7, lw=0,
               label="blank window (false-positive check)")
    if tic_prof:
        ax.plot(tt, np.mean(tic_prof, axis=0), color=MUTED, lw=2, label=f"TIC / baseline (mean of {len(tic_prof)})")
    if xic_prof:
        ax.plot(tt, 1 + np.mean(xic_prof, axis=0), color=style["line"], lw=2,
                label="expected primary-ion XIC (normalised, +1)")
    ax.set_xlabel("time from ejection apex (s)"); ax.set_ylabel("relative signal")
    ax.set_title(f"Plate {pa.plate_id} — mean ejection profile", loc="left", fontsize=_fs(11))
    ax.legend(frameon=False, fontsize=_fs(8), loc="upper right")
    _tight(fig)
    return style.save(fig, stem, tight=False)


def time_chart(pa, plate, stem, style: Style) -> Path:
    """TIC of the whole run with the time span of each plate row in alternating bands."""
    run, grid, P = pa.run, pa.grid, pa.period
    slots = {}
    for w, sl in pa.slot_of.items():
        if sl < len(grid):
            slots.setdefault(plate.split(w)[0], []).append((sl, w))
    tmin = run.t / 60.0
    fig, ax = plt.subplots(figsize=_grow(17, 5.6))
    ymax = float(run.tic.max()) * 1.16
    many = plate.n_rows > 16
    for i, r in enumerate(plate.rows):
        ss = sorted(slots.get(r, []))
        if not ss:
            continue
        a = (grid[ss[0][0]] - P / 2) / 60.0
        b = (grid[ss[-1][0]] + P / 2) / 60.0
        if i % 2 == 0:
            ax.axvspan(a, b, color=style["band"], lw=0, zorder=0)
        ax.axvline(a, color="#b0afaa", lw=0.5, zorder=1); ax.axvline(b, color="#b0afaa", lw=0.5, zorder=1)
        mid = (a + b) / 2
        ax.text(mid, ymax * 0.985, r, ha="center", va="top", fontsize=_fs(8) if many else 12, fontweight="bold", color=INK)
        if not many:
            ax.text(mid, ymax * 0.925, f"{ss[0][1]}→{ss[-1][1]}", ha="center", va="top", fontsize=_fs(6.5), color=MUTED)
    ax.plot(tmin, run.tic, lw=0.45, color=style["line"], zorder=2)
    ax.set_xlim(tmin[0], tmin[-1]); ax.set_ylim(0, ymax)
    ax.set_xlabel("time (min)"); ax.set_ylabel("TIC (million cps)")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v / 1e6:g}"))
    _clean(ax)
    ax.set_title(f"Plate {pa.plate_id} — total ion current over the run; bands = time window of each plate row "
                 f"(one ejection every {P:.3f} s)", loc="left", fontsize=_fs(11), color=INK)
    _tight(fig)
    return style.save(fig, stem, tight=False)


# ------------------------------------------------------------------------------------------
# One well: trace and spectrum (spectra PDF and viewer export)
# ------------------------------------------------------------------------------------------
def well_page(pa, well, s, style: Style, scan=None, mode="analysed", fig=None):
    """Figure with the ejection trace (top) and one spectrum with the theoretical m/z (bottom)."""
    from . import spectra
    tr = spectra.well_trace(pa, well, s)
    row = pa.results.set_index("Well").loc[well]
    fig = fig or plt.figure(figsize=_grow(11.7, 8.3))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.5], hspace=0.38, top=0.80, bottom=0.08, left=0.08, right=0.97)
    a1, a2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    lv = int(row["Confidence_Level"])
    head = (f"Plate {pa.plate_id} · {well} · {row.get('Sample_Code', '')} — level {lv if lv else '×'} "
            f"{LEVEL_LABELS[lv]}")
    sub = [f"Evidence {row.get('Evidence_Intensity', np.nan):,.0f} cps · raw {row.get('Raw_Ejection_Intensity', np.nan):,.0f} cps"
           f" · ions detected: {row.get('Ions_Detected', '') or 'none'}"]
    if isinstance(row.get("Isotope_Theoretical"), str) and row.get("Isotope_Theoretical"):
        sub.append(f"Isotopes {row.get('Isotope_Ion', '')}: expected {row['Isotope_Theoretical']}; observed "
                   f"{row.get('Isotope_Observed', '') or '–'} (score {row.get('Isotope_Pattern_Score', np.nan):.2f})")
    if row.get("Level_Note"):
        sub.append(str(row["Level_Note"]))
    if "HPLC_Result" in row.index:
        sub.append(f"{s.reference_label}: {row['HPLC_Result']}"
                   + (f" ({row['HPLC_Conversion']:.1f}%)" if np.isfinite(row.get("HPLC_Conversion", np.nan)) else ""))
    fig.text(0.08, 0.975, head, fontsize=_fs(13), fontweight="bold", color=INK, va="top")
    fig.text(0.08, 0.945, "\n".join(textwrap.shorten(x, 170) for x in sub), fontsize=_fs(8.5), color=MUTED, va="top",
             linespacing=1.45)
    if not tr:
        a1.text(0.5, 0.5, "no ejection assigned to this well", ha="center", transform=a1.transAxes)
        return fig
    a1.axvspan(*tr["window"], color=style["band"], alpha=0.8, lw=0, label="ejection window")
    a1.plot(tr["t"], tr["tic"], color="#b0afaa", lw=1, label="TIC")
    colors = [style["line"], style["line2"], "#6b6a66", "#2b2b29"]
    for (name, y), c in zip(tr["ions"].items(), colors):
        a1.plot(tr["t"], y, color=c, lw=1.4, marker=".", ms=3, label=f"{name} (product)")
    scan = tr["best_scan"] if scan is None else scan
    if scan is not None:
        a1.axvline(pa.run.t[scan], color=INK, lw=0.8, ls="--")
    for t0, w in tr["labels"]:
        a1.text(t0, 1.01, w, transform=a1.get_xaxis_transform(), ha="center", va="bottom", fontsize=_fs(8),
                color=INK if w == well else MUTED, fontweight="bold" if w == well else "normal")
    a1.set_yscale("symlog", linthresh=1e4); a1.set_ylim(bottom=0)
    a1.set_xlabel("time (s)"); a1.set_ylabel("intensity (cps)")
    a1.legend(frameon=False, fontsize=_fs(7.5), loc="upper right", ncol=2)
    _clean(a1)
    if scan is None:
        return fig
    mz, it = spectra.spectrum(pa, scan, mode)
    marks = spectra.markers(pa, well, s)
    lo, hi = spectra.view_range(marks)
    sel = (mz >= lo) & (mz <= hi)
    a2.vlines(mz[sel], 0, it[sel], color=style["line"], lw=1.0)
    ymax = max(float(it[sel].max()) if sel.any() else 0.0, s.cutoff_cps) * 1.25
    m_obs = {}
    for m in marks:
        if not lo <= m["mz"] <= hi:
            continue
        if m["kind"] == "ion":
            near = sel & (np.abs(mz - m["mz"]) <= pa.run.tol(m["mz"], s.mass_tolerance_da))
            m_obs[m["ion"]] = float(it[near].max()) if near.any() else 0.0
            a2.axvline(m["mz"], color=style["line2"], lw=1, ls="--")
            a2.text(m["mz"], ymax * 0.97, f"{m['label']}\n{m['mz']:.2f}", rotation=0, ha="center", va="top",
                    fontsize=_fs(7.5), color=style["line2"])
        elif m["kind"] == "reagent":
            a2.axvline(m["mz"], color=MUTED, lw=0.7, ls=":")
            a2.text(m["mz"], ymax * 0.80, m["label"], rotation=90, ha="right", va="top", fontsize=_fs(6), color=MUTED)
    for m in marks:     # theoretical isotope peaks, scaled to the observed M of the same ion
        if m["kind"] == "isotope" and lo <= m["mz"] <= hi and m_obs.get(m["ion"], 0) > 0:
            a2.plot(m["mz"], m["rel"] * m_obs[m["ion"]], marker="_", ms=14, mew=2, color=style["line2"])
    a2.axhline(s.cutoff_cps, color=INK, lw=0.7, ls="-.")
    a2.text(hi, s.cutoff_cps, " cutoff", va="bottom", ha="right", fontsize=_fs(7), color=INK)
    a2.set_xlim(lo, hi); a2.set_ylim(0, ymax)
    a2.set_xlabel("m/z"); a2.set_ylabel("intensity (cps)")
    grid_txt = "as analysed, centroid peaks" if pa.run.sparse else "as analysed, 0.1 Da grid"
    a2.set_title(f"Spectrum at {pa.run.t[scan]:.2f} s ({grid_txt if mode == 'analysed' else 'raw mzML'}); "
                 "dashed = calculated m/z, bars = theoretical isotopes, dotted = reagent ions", loc="left",
                 fontsize=_fs(9), color=INK)
    _clean(a2)
    return fig


def spectra_pdf(plates: dict, path, s, style: Style, levels=None) -> Path | None:
    """One page per well (trace + spectrum) for the wells whose level is in `levels`."""
    from matplotlib.backends.backend_pdf import PdfPages
    pages = [(pid, w) for pid, pa in plates.items() for w, lv in zip(pa.results["Well"], pa.results["Confidence_Level"])
             if levels is None or int(lv) in levels]
    if not pages:
        return None
    with PdfPages(path) as pdf:
        for pid, w in pages:
            fig = well_page(plates[pid], w, s, style)
            pdf.savefig(fig)
            plt.close(fig)
    return Path(path)


def calibration_chart(cal: dict, stem, s, style: Style) -> Path:
    """Errors vs threshold: false presents for the 'present' threshold, missed products for the 'absent' one."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=_grow(13, 4.8), gridspec_kw={"wspace": 0.25})
    pc, ac, tgt = cal["present_curve"], cal["absent_curve"], cal["target_pct"]
    if len(pc):
        a1.step(pc["threshold_cps"], pc["false_present_%_of_wells"], where="post", color=style["line"], lw=1.6)
        a1.axvline(s.level5_cps, color=MUTED, ls="--", lw=1, label=f"current level 5: {s.level5_cps:,.0f}")
        if np.isfinite(cal["suggested_level5_cps"]):
            a1.axvline(cal["suggested_level5_cps"], color=style["line2"], lw=1.4,
                       label=f"suggested: {cal['suggested_level5_cps']:,.0f}")
    a1.axhline(tgt, color=INK, ls=":", lw=1, label=f"target {tgt:g}% of wells")
    a1.set_xscale("log"); a1.set_xlabel("'present' threshold on the evidence (cps)")
    a1.set_ylabel("false presents (% of wells)")
    a1.set_title(_wrap("Trusting wells above a threshold as present", 50), loc="left", fontsize=_fs(11))
    if len(ac):
        a2.step(ac["threshold_cps"].clip(lower=1e3), ac["products_missed_%_of_wells"], where="post",
                color=style["line"], lw=1.6)
        a2.axvline(s.level2_raw_max_cps, color=MUTED, ls="--", lw=1, label=f"current levels 1-2: ≤ {s.level2_raw_max_cps:,.0f}")
        if np.isfinite(cal["suggested_absent_raw_cps"]):
            a2.axvline(max(cal["suggested_absent_raw_cps"], 1e3), color=style["line2"], lw=1.4,
                       label=f"suggested: ≤ {cal['suggested_absent_raw_cps']:,.0f}")
    a2.axhline(tgt, color=INK, ls=":", lw=1, label=f"target {tgt:g}% of wells")
    a2.set_xscale("log"); a2.set_xlabel("'absent' threshold on the raw signal (cps)")
    a2.set_ylabel("products missed (% of wells)")
    a2.set_title(_wrap("Trusting wells below a threshold as absent", 50), loc="left", fontsize=_fs(11))
    for ax in (a1, a2):
        _clean(ax); ax.legend(frameon=False, fontsize=_fs(8.5)); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
    _suptitle(fig, _wrap(f"Threshold calibration against {s.reference_label} (AUC evidence {cal['auc_evidence']}, "
                       f"raw signal {cal['auc_raw']})", 110), x=0.06, ha="left", fontsize=_fs(13), fontweight="bold",
                 color=INK, stem=stem)
    fig.subplots_adjust(top=0.84 - 0.08 * max(0.0, FONT["scale"] - 1))
    return style.save(fig, stem)


def run_qc_chart(q, plate_id, stem, style: Style) -> Path:
    has_is = "IS_vs_Median" in q.columns and q["IS_vs_Median"].notna().any()
    fig, axes = plt.subplots(2 if has_is else 1, 1, figsize=_grow(14, 6.5 if has_is else 4), squeeze=False)
    ax = axes[0, 0]
    ax.plot(q["Ejection_Time_s"] / 60, q["Pulse_vs_Median"], ".", ms=3, color=style["line"])
    ax.plot(q["Ejection_Time_s"] / 60, q["Pulse_vs_Median"].rolling(25, center=True, min_periods=5).median(),
            color=style["line2"], lw=1.6, label="running median (25 ejections)")
    flag = q["Run_QC_Flag"] != ""
    ax.plot(q.loc[flag, "Ejection_Time_s"] / 60, q.loc[flag, "Pulse_vs_Median"], "x", color=INK, label="flagged")
    ax.set_ylabel("ejection pulse height\n(× plate median)"); ax.legend(frameon=False, fontsize=_fs(8))
    ax.set_title(f"Plate {plate_id} — ejection pulse height over the run", loc="left", fontsize=_fs(11))
    if has_is:
        a2 = axes[1, 0]
        a2.plot(q["Ejection_Time_s"] / 60, q["IS_vs_Median"], ".", ms=3, color=style["line"])
        a2.axhline(0.5, color=INK, ls=":", lw=1)
        a2.set_ylabel("internal standard\n(× median)"); a2.set_title("Internal standard", loc="left", fontsize=_fs(11))
    for a in axes[:, 0]:
        _clean(a); a.set_xlabel("time (min)"); a.grid(color=GRID, lw=0.6)
    _tight(fig)
    return style.save(fig, stem, tight=False)


def heatmap_side_product(res, plate, stem, rule, style: Style, labels=None) -> Path:
    k = rule["key"]
    f = PlateFigure(plate, style)
    det = res[f"SP_{k}_Detected"].fillna(0).astype(int)
    for w, d in zip(res["Well"], det):
        f.cell(w, style["detected"] if d else style["not_detected"])
    return f.finish(stem, f"Plate {res['Plate'].iloc[0]} — side product: {rule['name']} ({int(det.sum())} wells)",
                    f"{rule['target']} {rule['change']}  ·  searched as [M+H]+ / [M-H]- with the product criteria "
                    "and ≥ the plate background", **_lab(labels),
                    handles=[_swatch(style["detected"], f"detected ({int(det.sum())})"),
                             _swatch(style["not_detected"], f"not detected ({int(len(det) - det.sum())})")])


def control_chart_figure(cc, stem, style: Style) -> Path:
    """One small panel per metric: values run after run, centre line and ±3 SD limits from earlier runs."""
    metrics = [m for m in cc["Metric"].unique() if cc.loc[cc["Metric"] == m, "Value"].notna().sum() >= 2]
    n = max(len(metrics), 1)
    ncol = 2 if n > 1 else 1
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=_grow(13, 2.6 * nrow + 0.6), squeeze=False)
    for ax, m in zip(axes.flat, metrics):
        d = cc[cc["Metric"] == m]
        x = d["Run"].to_numpy()
        ax.plot(x, d["Value"], "-", color="#b0afaa", lw=0.8)
        flag = d["Flag"].fillna("") != ""
        ax.plot(x[~flag.to_numpy()], d.loc[~flag, "Value"], "o", ms=4, color=style["line"])
        ax.plot(x[flag.to_numpy()], d.loc[flag, "Value"], "o", ms=6, color=style["line2"], label="flagged")
        ax.plot(x, d["Centre"], color=INK, lw=0.9)
        ax.plot(x, d["UCL"], color=MUTED, lw=0.8, ls="--")
        ax.plot(x, d["LCL"], color=MUTED, lw=0.8, ls="--")
        ax.set_title(_wrap(d["Label"].iloc[0], 60), loc="left", fontsize=_fs(10))
        ax.set_xlabel("run (plate)", fontsize=_fs(9))
        _clean(ax); ax.grid(color=GRID, lw=0.5)
    for ax in list(axes.flat)[len(metrics):]:
        ax.axis("off")
    _suptitle(fig, "Control charts across plates (centre and ±3 SD from the earlier runs; flagged = outside the "
                   "limits or CUSUM shift)", stem=stem, x=0.04, ha="left", fontsize=_fs(12), fontweight="bold")
    _tight(fig, rect=(0, 0, 1, 0.95))
    return style.save(fig, stem, tight=False)
