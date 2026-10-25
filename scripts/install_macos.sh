#!/bin/bash
# Native macOS setup with Homebrew: OpenFOAM.app (openfoam.com builds for macOS,
# https://github.com/gerlero/openfoam-app), FreeCAD, ParaView and a Python venv.
# Apple silicon gets the latest OpenFOAM; Intel Macs are limited to v2506.
# Docker (scripts/docker_run.sh) is the tested alternative on any Mac.
set -euo pipefail
cd "$(dirname "$0")/.."
case "$PWD" in *" "*) echo "ERROR: path contains spaces ($PWD); OpenFOAM cannot handle that. Clone elsewhere."; exit 1;; esac
command -v brew >/dev/null || { echo "Install Homebrew first: https://brew.sh"; exit 1; }

if [ "$(uname -m)" = "arm64" ]; then
    brew install gerlero/openfoam/openfoam
else
    brew install gerlero/openfoam/openfoam@2506
fi
brew install --cask freecad paraview
xcode-select -p >/dev/null 2>&1 || xcode-select --install   # compiler for the fvOption

python3 -m venv .venv
.venv/bin/pip install -r requirements.txt

cat <<'MSG'

Done. To work with the pipeline:
    openfoam                      # starts a shell with the OpenFOAM environment
    source .venv/bin/activate
    ./openfoam/Allwmake           # compile the PCM fvOption (once)
    make test && make smoke
FreeCAD and ParaView are found in /Applications automatically; override with
PCMPIPE_FREECADCMD / PCMPIPE_PVPYTHON if they live elsewhere.
MSG
