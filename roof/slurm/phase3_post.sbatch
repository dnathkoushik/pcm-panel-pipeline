#!/bin/bash
#SBATCH --job-name=roof-p3-post
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#
# PARAM Shakti: set the partition, and the account if your allocation needs one:
##SBATCH --partition=<partition>
##SBATCH --account=<account>
#
# Phase 3: after the Phase 2 jobs. For every finished case:
#   pvpython post/pv_extract.py   mid-plane slices (CSV), T contours with the melting
#                                 front (PNG), surface metrics from the fields
#   python3 post/collect.py       all cases against the baseline: summary_all.csv,
#                                 timeseries_*.csv and the paper-style figures
#                                 -> roof/results/phase3/
#   mkdir -p logs && sbatch slurm/phase3_post.sbatch
# PURGE=1 also deletes each case's field folders after extraction (keeps 0 and the
# latest): a finished case is then a few MB.
# pvpython: PVPYTHON (default: pvpython, from a ParaView module you add to
# PYTHON_MODULE in slurm/env.sh). Without a display post/pvpython_headless.sh starts a
# private Xvfb; an offscreen ParaView build (EGL/OSMesa) needs none.
set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"
if [ ! -d roofpipe ] || [ ! -f config/roof.yaml ]; then
    echo "Submit from the roof/ directory of the repository (the job started in $PWD)" >&2
    exit 2
fi
# shellcheck source=env.sh
source slurm/env.sh
runs=${ROOF_RUNS_DIR:-$PWD/runs}
read -r -a pv <<< "${PVPYTHON:-pvpython}"
purge=()
[ -n "${PURGE:-}" ] && purge=(--purge-fields)

n=0
for d in "$runs"/baseline_concrete "$runs"/*_[NI]_h[0-9][0-9][0-9]; do
    [ -f "$d/results/summary.json" ] || continue
    if [ -f "$d/results/pv/surface_metrics.csv" ]; then
        echo "skip $(basename "$d") (extracted)"
        continue
    fi
    echo "== $(basename "$d")"
    PVPYTHON="${pv[*]}" foam bash post/pvpython_headless.sh post/pv_extract.py "$d" \
        ${purge[@]+"${purge[@]}"} \
        || echo "warning: pv_extract failed for $d" >&2
    n=$((n + 1))
done
echo "== pvpython on $n case(s); collecting"
python3 post/collect.py
