"""Look of the browser interface: logo, page styles, summary cards, method profiles and friendly error
messages. Kept apart from webapp.py so the page code stays readable."""
from __future__ import annotations

import json
import re
from pathlib import Path

LOGO_SVG = """<svg xmlns="http://www.w3.org/2000/svg" width="46" height="46" viewBox="0 0 46 46" aria-hidden="true">
<rect x="2" y="7" width="42" height="32" rx="6" fill="#eef3fa" stroke="#1f5fae" stroke-width="2"/>
<g fill="#9ec5f4">{dots}</g>
<polyline points="6,33 14,33 17,30 20,33 23,33 26,13 29,33 33,33 35,28 37,33 40,33" fill="none" stroke="#1f5fae"
 stroke-width="2.2" stroke-linejoin="round" stroke-linecap="round"/></svg>""".replace(
    "{dots}", "".join(f'<circle cx="{9 + 7 * c}" cy="{14 + 6 * r}" r="1.9"/>' for r in range(3) for c in range(5)))

CSS = """<style>
header[data-testid="stHeader"]{background:transparent}
.block-container{padding-top:3.6rem;max-width:1400px}
.ems-head{display:flex;align-items:center;gap:14px;margin-bottom:.2rem}
.ems-head h1{font-size:1.75rem;margin:0;padding:0;font-weight:700;letter-spacing:-.01em}
.ems-head p{margin:0;color:#5f6b7a;font-size:.92rem}
.ems-card{background:#f6f8fb;border:1px solid #e3e8ef;border-radius:10px;padding:12px 14px;height:100%}
.ems-card .v{font-size:1.45rem;font-weight:700;color:#1f2933;line-height:1.2}
.ems-card .l{font-size:.78rem;color:#5f6b7a;text-transform:uppercase;letter-spacing:.04em}
.ems-card .s{font-size:.8rem;color:#5f6b7a}
.ems-check{border-left:4px solid #1f5fae;background:#f6f8fb;padding:8px 14px;border-radius:6px;margin:6px 0 12px}
.ems-check li{margin:2px 0}
.ems-foot{color:#8a94a3;font-size:.8rem;margin-top:2.5rem;border-top:1px solid #e3e8ef;padding-top:.6rem}
</style>"""

THEME_ARGS = ["--theme.base", "light", "--theme.primaryColor", "#1f5fae", "--theme.backgroundColor", "#ffffff",
              "--theme.secondaryBackgroundColor", "#f3f6fa", "--theme.textColor", "#1f2933",
              "--theme.font", "sans serif"]


def header_html(version: str) -> str:
    return (f'<div class="ems-head">{LOGO_SVG}<div><h1>TriADE</h1>'
            f'<p>Triage for acoustic droplet ejection MS · product detection in plate runs · confidence levels · HPLC-MS triage · '
            f'v{version}</p></div></div>')


def footer_html(version: str) -> str:
    return (f'<div class="ems-foot">TriADE {version} · MIT license · runs on this computer, nothing is '
            'uploaded to the internet · if you use it in published work, please cite the repository (CITATION.cff)</div>')


def card(label: str, value, sub: str = "") -> str:
    return (f'<div class="ems-card"><div class="l">{label}</div><div class="v">{value}</div>'
            f'<div class="s">{sub}</div></div>')


def checklist(items: list[tuple[bool, str]]) -> str:
    lis = "".join(f"<li>{'✅' if ok else '⬜'} {text}</li>" for ok, text in items)
    return f'<div class="ems-check"><b>Before running</b><ul style="margin:4px 0 0 0;padding-left:18px">{lis}</ul></div>'


# ------------------------------------------------------------------------------------------
# Method profiles (named settings kept on this computer)
# ------------------------------------------------------------------------------------------
def profile_dir() -> Path:
    d = Path.home() / ".triade" / "profiles"
    d.mkdir(parents=True, exist_ok=True)
    return d


def list_profiles() -> list[str]:
    return sorted(p.stem for p in profile_dir().glob("*.json"))


def _safe(name: str) -> str:
    n = re.sub(r"[^\w .,+()-]", "_", name).strip()
    if not n:
        raise ValueError("give the profile a name")
    return n[:80]


def save_profile(name: str, settings) -> Path:
    p = profile_dir() / f"{_safe(name)}.json"
    settings.save(p)
    return p


def load_profile(name: str) -> dict:
    return json.loads((profile_dir() / f"{_safe(name)}.json").read_text(encoding="utf-8"))


def delete_profile(name: str) -> None:
    (profile_dir() / f"{_safe(name)}.json").unlink(missing_ok=True)


# ------------------------------------------------------------------------------------------
# Friendly errors
# ------------------------------------------------------------------------------------------
HINTS = [
    (r"not in the sample table", "Check the Plate column of the sample table and the plate ID typed next to each mzML "
                                 "file in *Inputs*: they must be the same (e.g. 1 and 1)."),
    (r"well name\(s\) are not valid", "The wells do not exist in the chosen plate format. Check *Plate format* in "
                                      "the sidebar (96, 384 or 1536) and the Well column."),
    (r"listed more than once", "A well appears twice. For replicate ejections add a Replicate column (1, 2, 3) and "
                               "the Injection_Order of every ejection."),
    (r"column (Well|Product_SMILES) is missing", "Use the input template (Inputs → Download the input template) or "
                                                 "rename the columns as in the template."),
    (r"No ejection order fits|ejection pulses but the order", "The number of ejections in the run does not match the "
     "wells: check the plate format and the interval between ejections, or give Injection_Order."),
    (r"MS-Numpress", "Convert the file again with msConvert without Numpress compression (zlib is fine)."),
    (r"no MS1 spectra", "The mzML has no MS1 spectra: convert again keeping MS1 (msConvert filter msLevel 1)."),
    (r"same plate ID", "Two mzML files got the same plate: change the plate ID next to one of them."),
    (r"thresholds must satisfy", "Level 3 ≤ cutoff ≤ level 5: adjust the custom thresholds."),
    (r"side-product rule", "Write each rule as  name: Product (or Reagent_1): -Cl +OH  — one per line."),
    (r"Control must be", "In the Control column use positive, negative or blank (or leave it empty)."),
    (r"review file", "The review file needs the columns Plate, Well and Decision (present, absent or hplc)."),
    (r"replicate ejections need", "Add the Injection_Order column with the order of every ejection."),
]


def friendly(err: Exception) -> tuple[str, str]:
    msg = str(err)
    for pattern, hint in HINTS:
        if re.search(pattern, msg, re.I):
            return msg, hint
    return msg, "Check the inputs in the *Preview* tab; the log (log.txt in the results) has more detail."


def eta_text(fraction: float, elapsed: float) -> str:
    if fraction <= 0.03 or elapsed < 2:
        return ""
    left = elapsed * (1 - fraction) / fraction
    return f" · about {left:.0f} s left" if left < 90 else f" · about {left / 60:.1f} min left"
