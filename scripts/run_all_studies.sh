#!/bin/bash
# Runs every study end to end. Expect days of CPU time at the medium mesh;
# start with `smoke` and `validation`. Usage: scripts/run_all_studies.sh [NPROCS]
set -euo pipefail
NP=${1:-8}
cd "$(dirname "$0")/.."
[ -n "${WM_PROJECT_DIR:-}" ] || source /usr/lib/openfoam/openfoam2406/etc/bashrc
for study in smoke validation shapes delta_T shell_thickness protrusion_ratio gap volume; do
    echo "################ $study"
    python3 -m pcmpipe all --study "$study" --np "$NP" --np-mesh "$NP" --keep-going
done
