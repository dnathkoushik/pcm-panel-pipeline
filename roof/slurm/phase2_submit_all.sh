#!/bin/bash
# Phase 2: submit every PCM case (3 shapes x 2 orientations x 5 heights = 30) as its
# own SLURM job. Run on the login node, from roof/:
#
#   bash slurm/phase2_submit_all.sh                 # all cases not yet passed
#   bash slurm/phase2_submit_all.sh pyr_ hem_N      # only ids containing these texts
#   DRY=1 bash slurm/phase2_submit_all.sh           # print the sbatch lines only
#
# It first compiles the PCM library (libpcmEnthalpyPorosity.so) once, so the 30 jobs
# do not race to build it. A case whose results/summary.json says "passed": true is
# skipped; delete <runs>/<id> to run it again.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=env.sh
source slurm/env.sh
roof_check_env
roof_check_lib --build

runs=${ROOF_RUNS_DIR:-$PWD/runs}
mkdir -p logs
mapfile -t ids < <(python3 -m roofpipe cases | awk 'NF > 3 {print $1}')
n=0
for id in "${ids[@]}"; do
    if [ "$#" -gt 0 ]; then
        hit=0
        for pat in "$@"; do [[ $id == *"$pat"* ]] && hit=1; done
        [ "$hit" -eq 1 ] || continue
    fi
    if grep -q '"passed": true' "$runs/$id/results/summary.json" 2>/dev/null; then
        echo "skip $id (passed)"
        continue
    fi
    cmd=(sbatch --job-name="roof-$id" --export="ALL,CASE_ID=$id" slurm/phase2_case.sbatch)
    if [ -n "${DRY:-}" ]; then
        echo "${cmd[*]}"
    else
        "${cmd[@]}"
    fi
    n=$((n + 1))
done
echo "$n job(s) $([ -n "${DRY:-}" ] && echo "to submit" || echo "submitted"); logs in roof/logs/"
