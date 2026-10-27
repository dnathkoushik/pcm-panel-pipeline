#!/bin/bash
# Long runs on machines with a job time limit (GitHub-hosted runners: 6 h per job).
# Runs CAD, mesh and the solver for the matching case(s) of one study in screening
# mode, inside the container, and stops the solver cleanly (stopAt writeNow) once
# BUDGET_MIN minutes have passed. If runs/<case>/case already exists (restart files
# from an earlier run), only the solver runs, continuing from the last written time
# (controlDict: startFrom latestTime).
#
#   scripts/ci_long_run.sh STUDY CASE NP BUDGET_MIN [--set KEY=VALUE ...]
#   scripts/ci_long_run.sh shapes pyr_ 4 280 --set run.write_interval_min=50
#
# Results, 3D views and plots are separate calls: pcmpipe results | post | plots.
set -euo pipefail
study=$1 case=$2 np=$3 budget=$4
shift 4
cd "$(dirname "$0")/.."
pcm=(scripts/docker_run.sh python3 -m pcmpipe)
opts=(--study "$study" --mode screening --case "$case" --np "$np" --np-mesh "$np" "$@")

dicts() { ls runs/*"$case"*/case/system/controlDict 2>/dev/null || true; }

# Rewrite in place (same inode), so the solver sees the change however it watches files.
set_stop() {
    local text
    text=$(sed "s/^stopAt .*/stopAt          $2;/" "$1")
    printf '%s\n' "$text" > "$1"
}

if [ -n "$(dicts)" ]; then
    echo "Restart files found: continuing from the last written time"
    for f in $(dicts); do
        set_stop "$f" endTime
        log="$(dirname "$(dirname "$f")")/log.chtMultiRegionFoam"
        if [ -f "$log" ]; then mv "$log" "$log.$(date +%s)"; fi
    done
    steps=(run)
else
    steps=(cad mesh setup run)
fi

deadline=$(( $(date +%s) + 60 * budget ))
(
    while [ "$(date +%s)" -lt "$deadline" ]; do sleep 20; done
    while true; do                       # keep checking: the case may be set up late
        for f in $(dicts); do
            if ! grep -q '^stopAt .*writeNow' "$f"; then
                set_stop "$f" writeNow
                echo "Wall-clock budget of $budget min used: stopAt writeNow in $f"
            fi
        done
        sleep 20
    done
) &
watchdog=$!
trap 'kill "$watchdog" 2>/dev/null || true' EXIT

status=0
for s in "${steps[@]}"; do
    "${pcm[@]}" "$s" "${opts[@]}" || { status=$?; break; }
done

# After a solver failure the written times are still decomposed: reconstruct them
# so that results and post-processing can use what was computed.
if [ "$status" -ne 0 ] && ls -d runs/*"$case"*/case/processor0 >/dev/null 2>&1; then
    echo "Solver failed: reconstructing the written times"
    scripts/docker_run.sh python3 -c "
import glob, os, subprocess
for d in glob.glob('runs/*$case*/case'):
    if os.path.isdir(os.path.join(d, 'processor0')):
        subprocess.run(['reconstructPar', '-allRegions', '-newTimes'], cwd=d)
" || true
fi
exit "$status"
