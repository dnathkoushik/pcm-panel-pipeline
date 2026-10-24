#!/bin/bash
# One-time setup inside WSL2 Ubuntu 24.04 (or native Ubuntu 22.04/24.04).
# Installs OpenFOAM (openfoam.com) v2406, FreeCAD, ParaView and the Python deps,
# then compiles the PCM enthalpy-porosity fvOption library.
set -euo pipefail
OF_VERSION=${OF_VERSION:-2406}

sudo apt-get update
sudo apt-get install -y curl ca-certificates
curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash
sudo apt-get update
sudo apt-get install -y "openfoam${OF_VERSION}-default" \
    freecad paraview \
    python3-numpy python3-matplotlib python3-yaml

BASHRC_LINE="source /usr/lib/openfoam/openfoam${OF_VERSION}/etc/bashrc"
grep -qxF "$BASHRC_LINE" ~/.bashrc || echo "$BASHRC_LINE" >> ~/.bashrc
# shellcheck disable=SC1090
source "/usr/lib/openfoam/openfoam${OF_VERSION}/etc/bashrc" || true

cd "$(dirname "$0")/.."
case "$PWD" in *" "*) echo "ERROR: path contains spaces ($PWD); copy the project to ~/ first"; exit 1;; esac
./openfoam/Allwmake
echo
echo "Checks:"
command -v chtMultiRegionFoam snappyHexMesh freecadcmd pvpython
ls -l "$FOAM_USER_LIBBIN/libpcmEnthalpyPorosity.so"
