#!/bin/bash
# TriADE: double-click (macOS) or run ./start_mac_linux.command (Linux) to open the browser
# interface. The first time it creates a Python environment (.venv) and installs the libraries
# (needs internet, takes a few minutes).
cd "$(dirname "$0")" || exit 1
PY=""
for v in python3.12 python3.13 python3.11 python3.10 python3; do
    if command -v "$v" >/dev/null 2>&1 && "$v" -c 'import sys; sys.exit(0 if (3,10) <= sys.version_info[:2] <= (3,13) else 1)'; then
        PY="$v"; break
    fi
done
if [ -z "$PY" ]; then
    echo "Python 3.10 to 3.13 was not found. Install Python 3.12 from https://www.python.org/downloads/ and try again."
    read -r -p "Press Enter to close."; exit 1
fi
if [ ! -x ".venv/bin/python" ]; then
    echo "First run: installing the libraries. This takes a few minutes..."
    "$PY" -m venv .venv && .venv/bin/python -m pip install --quiet --upgrade pip \
        && .venv/bin/python -m pip install --quiet -e ".[app]" || {
        echo "Installation failed. Check the internet connection, delete the .venv folder and try again."
        read -r -p "Press Enter to close."; exit 1; }
fi
echo "Opening the interface in your browser. Press Ctrl+C (or close this window) to stop."
exec .venv/bin/triade-app
