#!/bin/bash
# Mesh and solve one generated roof case with chtMultiRegionFoam.
#
#   bash scripts/run_slab_case.sh <case_dir> [n_procs]
#
# Idempotent: called again (e.g. by a resubmitted job after a time-limit kill) it
# skips the mesh and the decomposition that exist and continues from the latest
# written time (controlDict: startFrom latestTime). Needs the environment of
# slurm/env.sh, which it sources if the caller has not.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../slurm/env.sh
source "$here/../slurm/env.sh"

case_dir=${1:?usage: run_slab_case.sh <case_dir> [n_procs]}
np=${2:-1}
cd "$case_dir"
read -r -a launcher <<< "$MPI_LAUNCHER"

step() {                                    # step <log file> <command ...>
    local log=$1
    shift
    echo "[$(date '+%F %T')] $*  > $log"
    if ! foam "$@" > "$log" 2>&1; then
        echo "FAILED: $*. Last lines of $case_dir/$log:" >&2
        tail -n 40 "$log" >&2
        return 1
    fi
}

keep_old_log() {                            # a restart keeps the previous solver log
    if [ -f "$1" ]; then mv "$1" "$1.$(date +%Y%m%d-%H%M%S)"; fi
}

# ---- mesh: one block, copied to the region `concrete`
if [ ! -s constant/concrete/polyMesh/owner ] && [ ! -s constant/concrete/polyMesh/owner.gz ]; then
    step log.blockMesh blockMesh
    rm -rf constant/concrete/polyMesh
    mkdir -p constant/concrete
    cp -r constant/polyMesh constant/concrete/polyMesh
    step log.checkMesh.concrete checkMesh -region concrete
fi

# ---- solve
keep_old_log log.chtMultiRegionFoam
if [ "$np" -gt 1 ]; then
    dict_np=$(sed -n 's/^numberOfSubdomains[[:space:]]*\([0-9][0-9]*\);.*/\1/p' system/decomposeParDict)
    if [ "$dict_np" != "$np" ]; then
        echo "system/decomposeParDict is written for $dict_np ranks, the job has $np:" \
             "regenerate the case with --np $np, or submit with --ntasks=$dict_np" >&2
        exit 2
    fi
    if [ -d processor0 ]; then
        n_dec=$(find . -maxdepth 1 -type d -name 'processor*' | wc -l)
        if [ "$n_dec" -ne "$np" ]; then
            echo "the case is decomposed into $n_dec parts, the job has $np ranks" >&2
            exit 2
        fi
        echo "processor directories found: continuing the decomposed run"
    else
        step log.decomposePar decomposePar -allRegions -force -latestTime
    fi
    step log.chtMultiRegionFoam "${launcher[@]}" "$np" chtMultiRegionFoam -parallel
    step log.reconstructPar reconstructPar -allRegions -newTimes
    rm -rf processor*
else
    step log.chtMultiRegionFoam chtMultiRegionFoam
fi

if ! grep -q '^End' log.chtMultiRegionFoam; then
    echo "chtMultiRegionFoam did not print End: see $case_dir/log.chtMultiRegionFoam" >&2
    exit 1
fi
echo "[$(date '+%F %T')] solver finished: $case_dir"
