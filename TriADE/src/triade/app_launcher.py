"""`triade-app`: opens the browser interface (Streamlit) on this computer."""
import sys
from pathlib import Path


def main() -> None:
    try:
        from streamlit.web import cli as stcli
    except ImportError:
        sys.exit('The browser interface needs Streamlit: pip install "triade[app]"  (or: pip install streamlit)')
    app = Path(__file__).with_name("webapp.py")
    from triade.ui import THEME_ARGS
    sys.argv = ["streamlit", "run", str(app), "--server.maxUploadSize", "4000",
                "--browser.gatherUsageStats", "false", "--server.address", "localhost"] + THEME_ARGS + sys.argv[1:]
    sys.exit(stcli.main())


if __name__ == "__main__":
    main()
