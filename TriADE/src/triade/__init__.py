"""TriADE (triage for acoustic droplet ejection mass spectrometry): product detection in acoustic-ejection mass spectrometry (Echo-MS) plate runs,
with confidence levels and HPLC-MS triage."""
__version__ = "1.4.0"

from .adducts import IONS, PRIMARY, default_ions  # noqa: E402,F401
from .plates import PLATE_FORMATS, get_plate_format  # noqa: E402,F401
from .settings import LEVEL_LABELS, Settings  # noqa: E402,F401
from .pipeline import Inputs, Results, run_analysis  # noqa: E402,F401
