#!/bin/bash
cd "$(dirname "$0")"

echo "============================================"
echo "  Cluster-Scout - Beta Launcher"
echo "============================================"
echo

if ! command -v uv &> /dev/null; then
    echo "[1/3] uv not found - installing it now, please wait..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="$HOME/.local/bin:$PATH"
    if ! command -v uv &> /dev/null; then
        echo
        echo "ERROR: Could not find or install uv automatically."
        echo "Please install it manually from https://astral.sh/uv and re-run this script."
        read -n 1 -s -r -p "Press any key to exit..."
        exit 1
    fi
    echo "  uv installed successfully."
else
    echo "[1/3] uv found."
fi

echo
echo "[2/3] Installing/updating dependencies - first run may take a few minutes..."
if ! uv sync; then
    echo
    echo "ERROR: Dependency installation failed. See the output above for details."
    read -n 1 -s -r -p "Press any key to exit..."
    exit 1
fi

echo
echo "[3/3] Launching Cluster-Scout..."
echo
uv run app.py
status=$?
if [ $status -ne 0 ]; then
    echo
    echo "Cluster-Scout closed with an error (see above)."
    echo "If you need help, send a screenshot of this window to Matt."
    read -n 1 -s -r -p "Press any key to exit..."
fi
