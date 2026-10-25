#!/bin/bash
# Native setup for Ubuntu 22.04/24.04 and Debian 12 (also inside WSL2).
# Installs OpenFOAM (openfoam.com), FreeCAD, ParaView and the Python packages,
# then compiles the PCM enthalpy-porosity fvOption library.
# Other distributions: use Docker (scripts/docker_run.sh) instead.
set -euo pipefail
OF_VERSION=${OF_VERSION:-2406}

cd "$(dirname "$0")/.."
case "$PWD" in *" "*) echo "ERROR: path contains spaces ($PWD); OpenFOAM cannot handle that. Clone elsewhere."; exit 1;; esac

sudo apt-get update
sudo apt-get install -y curl ca-certificates gnupg software-properties-common
# Ubuntu 24.04 dropped FreeCAD from its archive: use the official FreeCAD PPA there
if [ -z "$(apt-cache policy freecad 2>/dev/null | sed -n 's/^ *Candidate: //p' | grep -v none)" ] \
   && grep -qi ubuntu /etc/os-release; then
    sudo add-apt-repository -y ppa:freecad-maintainers/freecad-stable
fi
curl -fsSL https://dl.openfoam.com/add-debian-repo.sh | sudo bash
sudo apt-get update
sudo apt-get install -y "openfoam${OF_VERSION}-default" \
    freecad paraview xvfb xauth \
    python3-numpy python3-matplotlib python3-yaml

BASHRC_LINE="source /usr/lib/openfoam/openfoam${OF_VERSION}/etc/bashrc"
grep -qxF "$BASHRC_LINE" ~/.bashrc || echo "$BASHRC_LINE" >> ~/.bashrc
set +u
# shellcheck disable=SC1090
source "/usr/lib/openfoam/openfoam${OF_VERSION}/etc/bashrc"
set -u

./openfoam/Allwmake
echo
echo "Installed:"
command -v chtMultiRegionFoam snappyHexMesh freecadcmd pvpython
ls -l "$FOAM_USER_LIBBIN/libpcmEnthalpyPorosity.so"
echo "Open a new terminal (or 'source ~/.bashrc'), then: make test && make smoke"
