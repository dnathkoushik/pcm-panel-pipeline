#!/bin/bash
# CAD, mesh and solve one Phase 2 case written by `python3 -m roofpipe pcm`.
#
#   bash scripts/run_pcm_case.sh <run_dir> [n_procs]      (from roof/)
#
#   1. CAD     FreeCAD (cad/build_roof_pcm_freecad.py) if freecadcmd exists, else
#              python3 -m roofpipe stl  ->  <run_dir>/cad/pcm_body.stl
#   2. mesh    blockMesh, surfaceFeatureExtract, snappyHexMesh, splitMeshRegions
#              -cellZones (regions concrete | pcm), checkMesh per region
#   3. mass    python3 -m roofpipe finalize: PCM density x V_exact / V_mesh
#   4. solve   scripts/run_slab_case.sh (decomposePar, chtMultiRegionFoam, reconstructPar)
#   5. prune   field folders before the analysed (last) day, except 0
#
# Idempotent like run_slab_case.sh: what exists is not redone, the solver continues.
# MESH_ONLY=1 stops after the mesh and the PCM mass (to inspect a geometry).
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
roof="$(dirname "$here")"
# shellcheck source=../slurm/env.sh
source "$roof/slurm/env.sh"

run_dir=$(cd "${1:?usage: run_pcm_case.sh <run_dir> [n_procs]}" && pwd)
np=${2:-1}
case_dir="$run_dir/case"
name=$(basename "$run_dir")

step() {                                    # step <log file> <command ...>
    local log=$1
    shift
    echo "[$(date '+%F %T')] $name: $*  > $log"
    if ! foam "$@" > "$log" 2>&1; then
        echo "FAILED: $*. Last lines of $(pwd)/$log:" >&2
        tail -n 40 "$log" >&2
        return 1
    fi
}

# ---- 1. CAD
if [ ! -s "$run_dir/cad/pcm_body.stl" ]; then
    cd "$run_dir/cad"
    if foam bash -c 'command -v freecadcmd' >/dev/null 2>&1; then
        echo "[$(date '+%F %T')] $name: FreeCAD  > cad/log.freecad"
        if ! ROOF_PCM_JSON="$run_dir/cad/body.json" foam freecadcmd \
                "$roof/cad/build_roof_pcm_freecad.py" > log.freecad 2>&1 \
           || [ ! -s pcm_body.stl ]; then
            echo "FAILED: FreeCAD, see $run_dir/cad/log.freecad" >&2
            tail -n 20 log.freecad >&2
            exit 1
        fi
    else
        echo "[$(date '+%F %T')] $name: no freecadcmd, writing the STL with roofpipe/stl.py"
        (cd "$roof" && python3 -m roofpipe stl "$run_dir")
    fi
fi

# ---- 2. mesh (the solver's controlDict is swapped out while meshing)
cd "$case_dir"
if [ -f system/controlDict.solve ]; then       # an interrupted mesh run: restore first
    mv system/controlDict.solve system/controlDict
fi
if [ ! -s constant/concrete/polyMesh/owner ] && [ ! -s constant/concrete/polyMesh/owner.gz ]; then
    mkdir -p constant/triSurface
    cp "$run_dir/cad/pcm_body.stl" constant/triSurface/pcm_body.stl
    cp system/controlDict system/controlDict.solve
    cp system/controlDict.mesh system/controlDict
    step log.blockMesh blockMesh
    step log.surfaceFeatureExtract surfaceFeatureExtract
    step log.snappyHexMesh snappyHexMesh -overwrite
    step log.splitMeshRegions splitMeshRegions -cellZones -overwrite
    for r in concrete pcm; do
        for f in cellLevel pointLevel level0Edge surfaceIndex refinementHistory; do
            rm -f "constant/$r/polyMesh/$f"
        done
    done
    rm -rf constant/polyMesh 0/cellToRegion 0/*/cellToRegion
    mv system/controlDict.solve system/controlDict
    for r in concrete pcm; do
        step "log.checkMesh.$r" checkMesh -region "$r" -constant
        if grep -q "Failed [0-9]* mesh checks" "log.checkMesh.$r"; then
            echo "warning: checkMesh reports failed checks for region $r (see log.checkMesh.$r)" >&2
        fi
    done
fi

# ---- 3. PCM mass from the meshed volume (once)
if ! grep -q '"meshed"' "$run_dir/case_meta.json"; then
    (cd "$roof" && python3 -m roofpipe finalize "$run_dir")
fi

if [ -n "${MESH_ONLY:-}" ]; then
    echo "[$(date '+%F %T')] $name: meshed (MESH_ONLY set, not solving)"
    exit 0
fi

# ---- 4. solve
bash "$roof/scripts/run_slab_case.sh" "$case_dir" "$np"

# ---- 5. keep only the analysed day's fields (and 0)
end=$(sed -n 's/^endTime[[:space:]]*\([0-9.eE+-]*\);.*/\1/p' system/controlDict)
python3 - "$case_dir" "$end" <<'EOF'
import os, shutil, sys
case, end = sys.argv[1], float(sys.argv[2])
for d in os.listdir(case):
    try:
        t = float(d)
    except ValueError:
        continue
    if 0.0 < t < end - 86400.0 - 1e-6:
        shutil.rmtree(os.path.join(case, d))
EOF
echo "[$(date '+%F %T')] $name: done"
