#!/bin/bash
#SBATCH --job-name=roof-p1-baseline
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
#
# PARAM Shakti: set the partition, and the account if your allocation needs one
# (`sinfo -s` lists the partitions; without --partition the default one is used):
##SBATCH --partition=<partition>
##SBATCH --account=<account>
#
# Phase 1: the pure concrete roof slab (the paper's Ref-SU) over run.n_days design
# days, then the Phase 1 gates. The case is small (54 000 cells, 8640 steps): a few
# minutes on 8 cores.
#
#   cd pcm-panel-pipeline/roof
#   sbatch slurm/phase1_baseline.sbatch
#
# Variants: CASE_NAME (run folder) and ROOF_SET (space-separated KEY=VALUE overrides
# of config/roof.yaml), e.g. the paper's insulated test rig:
#   sbatch --export=ALL,CASE_NAME=baseline_specimen,ROOF_SET="indoor.mode=specimen" \
#          slurm/phase1_baseline.sbatch
#
# A resubmitted job continues an existing case. Delete <runs>/<name>/case to start over.
set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"
if [ ! -d roofpipe ] || [ ! -f config/roof.yaml ]; then
    echo "Submit from the roof/ directory of the repository (the job started in $PWD)" >&2
    exit 2
fi
# shellcheck source=env.sh
source slurm/env.sh
roof_check_env

NP=${SLURM_NTASKS:-1}
NAME=${CASE_NAME:-baseline_concrete}
RUN_DIR=${ROOF_RUNS_DIR:-$PWD/runs}/$NAME

echo "== Phase 1 baseline: $NAME, $NP MPI rank(s), $(hostname), $(date '+%F %T')"
if [ -d "$RUN_DIR/case" ]; then
    echo "== $RUN_DIR/case exists: continuing it"
else
    set_args=()
    for kv in ${ROOF_SET:-}; do
        set_args+=(--set "$kv")
    done
    python3 -m roofpipe baseline --name "$NAME" --np "$NP" ${set_args[@]+"${set_args[@]}"}
fi

bash scripts/run_slab_case.sh "$RUN_DIR/case" "$NP"

echo "== Phase 1 gates"
python3 -m roofpipe check "$RUN_DIR"
echo "== done $(date '+%F %T'): $RUN_DIR/results"
