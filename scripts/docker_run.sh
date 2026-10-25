#!/bin/bash
# Run any command inside the pipeline container with this repository mounted at
# /work. Builds the image on first use. Linux, macOS and Windows (Git Bash/WSL).
#
#   scripts/docker_run.sh                                  # interactive shell
#   scripts/docker_run.sh python3 -m pcmpipe all --study smoke --np 2
#
set -euo pipefail
IMAGE=${PCM_IMAGE:-pcm-panel-pipeline}
cd "$(dirname "$0")/.."
REPO="$(pwd)"

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
    echo "Building $IMAGE (first run, ~10 min) ..."
    docker build --platform linux/amd64 -t "$IMAGE" .
fi

# HOME inside the repo (git-ignored) keeps FreeCAD settings and the compiled
# fvOption library between runs, whatever uid the container runs as.
args=(--rm --platform linux/amd64 -v "$REPO":/work -w /work -e HOME=/work/.docker-home)
[ -t 0 ] && [ -t 1 ] && args+=(-it)
# On Linux, run as the host user so files in runs/ are owned by you.
# Docker Desktop (macOS/Windows) maps file ownership itself.
if [ "$(uname -s)" = "Linux" ]; then
    args+=(--user "$(id -u):$(id -g)")
fi

exec docker run "${args[@]}" "$IMAGE" "${@:-bash}"
