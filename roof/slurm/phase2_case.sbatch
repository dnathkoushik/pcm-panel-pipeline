#!/bin/bash
#SBATCH --job-name=roof-p2
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#
# PARAM Shakti: set the partition, and the account if your allocation needs one:
##SBATCH --partition=<partition>
##SBATCH --account=<account>
#
# Phase 2: ONE PCM case (CASE_ID, e.g. pyr_N_h070): case set-up, CAD, mesh, solver,
# gates. About 216 000 cells; 15 min on 4 cores of a workstation, so 2 h is ample.
# Submit every case with slurm/phase2_submit_all.sh; one case by hand:
#   mkdir -p logs && sbatch --export=ALL,CASE_ID=pyr_N_h070 slurm/phase2_case.sbatch
# A resubmitted job continues the case. Delete <runs>/<CASE_ID> to start over.
set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"
if [ ! -d roofpipe ] || [ ! -f config/roof.yaml ]; then
    echo "Submit from the roof/ directory of the repository (the job started in $PWD)" >&2
    exit 2
fi
: "${CASE_ID:?set CASE_ID, e.g. sbatch --export=ALL,CASE_ID=pyr_N_h070 ...}"
# shellcheck source=env.sh
source slurm/env.sh
roof_check_env
roof_check_lib

NP=${SLURM_NTASKS:-1}
RUN_DIR=${ROOF_RUNS_DIR:-$PWD/runs}/$CASE_ID

echo "== Phase 2 case $CASE_ID: $NP MPI rank(s), $(hostname), $(date '+%F %T')"
if [ -d "$RUN_DIR/case" ]; then
    echo "== $RUN_DIR/case exists: continuing it"
else
    set_args=()
    for kv in ${ROOF_SET:-}; do
        set_args+=(--set "$kv")
    done
    python3 -m roofpipe pcm --case "$CASE_ID" --np "$NP" ${set_args[@]+"${set_args[@]}"}
fi

bash scripts/run_pcm_case.sh "$RUN_DIR" "$NP"

echo "== Phase 2 gates"
python3 -m roofpipe check "$RUN_DIR"
echo "== done $(date '+%F %T'): $RUN_DIR/results"
