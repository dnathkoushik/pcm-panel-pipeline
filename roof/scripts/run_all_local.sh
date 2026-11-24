#!/bin/bash
# Phase 2 on a workstation, without SLURM: every PCM case, JOBS at a time with NP
# MPI ranks each (default 3 x 4 = 12 cores), inside the repository's Docker image.
#
#   bash roof/scripts/run_all_local.sh                 (from the repository root)
#   JOBS=2 NP=6 bash roof/scripts/run_all_local.sh pyr_
#
# Cases that already passed are skipped; each case logs to roof/logs/<id>.log.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
JOBS=${JOBS:-3}
NP=${NP:-4}
mkdir -p roof/logs
scripts/docker_run.sh bash -c 'cd roof && python3 -m roofpipe cases' | awk 'NF > 3 {print $1}' > roof/logs/.cases
todo=()
while read -r id; do
    if [ "$#" -gt 0 ]; then
        hit=0
        for pat in "$@"; do [[ $id == *"$pat"* ]] && hit=1; done
        [ "$hit" -eq 1 ] || continue
    fi
    if grep -q '"passed": true' "roof/runs/$id/results/summary.json" 2>/dev/null; then
        echo "skip $id (passed)"
        continue
    fi
    todo+=("$id")
done < roof/logs/.cases
echo "${#todo[@]} case(s), $JOBS at a time, $NP ranks each"
printf '%s\n' "${todo[@]}" | xargs -P "$JOBS" -I{} bash -c '
    id={}
    scripts/docker_run.sh bash -c "cd roof && { [ -d runs/$id/case ] || python3 -m roofpipe pcm --case $id --np '"$NP"'; } && bash scripts/run_pcm_case.sh runs/$id '"$NP"' && python3 -m roofpipe check $id" > roof/logs/$id.log 2>&1 \
        && echo "$(date +%T) done   $id: $(tail -1 roof/logs/$id.log)" \
        || echo "$(date +%T) FAILED $id: see roof/logs/$id.log"'
