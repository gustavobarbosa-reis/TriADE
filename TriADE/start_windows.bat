@echo off
rem TriADE: double-click to open the browser interface.
rem The first time it creates a Python environment for the program (.venv folder) and installs the
rem libraries (needs internet, takes a few minutes). After that it only starts.
chcp 65001 >nul
cd /d "%~dp0"

set "PY="
for %%V in (3.12 3.13 3.11 3.10) do (
    if not defined PY (
        py -%%V -c "import sys" >nul 2>&1 && set "PY=py -%%V"
    )
)
if not defined PY (
    python -c "import sys; sys.exit(0 if (3, 10) <= sys.version_info[:2] <= (3, 13) else 1)" >nul 2>&1 && set "PY=python"
)
if not defined PY (
    echo Python 3.10 to 3.13 was not found.
    echo Install Python 3.12 from https://www.python.org/downloads/ ^(tick "Add python.exe to PATH"^) and try again.
    pause
    exit /b 1
)

if not exist ".venv\Scripts\python.exe" (
    echo First run: installing the libraries. This takes a few minutes...
    %PY% -m venv .venv || goto :failed
    ".venv\Scripts\python.exe" -m pip install --quiet --upgrade pip
    ".venv\Scripts\python.exe" -m pip install --quiet -e ".[app]" || goto :failed
)
echo Opening the interface in your browser. Close this window to stop the program.
".venv\Scripts\triade-app.exe"
exit /b 0

:failed
echo Installation failed. Check the internet connection, delete the .venv folder and try again.
pause
exit /b 1
