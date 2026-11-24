# shellcheck shell=bash
# Environment for the roof-embedded PCM jobs on PARAM Shakti (or any SLURM cluster).
# Sourced by slurm/*.sbatch and scripts/run_slab_case.sh. Edit the settings block
# once, or export the variables before `sbatch` (the job inherits them).
#
# OpenFOAM must be the openfoam.com (ESI) release, v2306 or newer: v2406 is what the
# repository's CI and Docker image use. Foundation releases (openfoam.org, version
# numbers without a "v") are NOT supported. Give ONE of these; the first one set wins:
#   OPENFOAM_BASHRC  an installed etc/bashrc, e.g. $HOME/OpenFOAM/OpenFOAM-v2406/etc/bashrc
#   OPENFOAM_MODULE  environment module(s), space-separated and loaded in order, e.g.
#                    "gcc/9.3 openmpi/4.1 openfoam/v2406"      (`module avail openfoam`)
#   PCM_SIF          an Apptainer/Singularity image of the repository's Docker image
#                    (single node, MPI inside the container); roof/README.md shows how
#                    to build it. Leave all three empty if OpenFOAM is already loaded.

# --------------------------------------------------------------- settings
: "${OPENFOAM_BASHRC:=}"
: "${OPENFOAM_MODULE:=}"
: "${PCM_SIF:=}"
: "${PCM_SIF_BASHRC:=/usr/lib/openfoam/openfoam2406/etc/bashrc}"   # path inside the image
: "${PYTHON_MODULE:=}"      # module(s) with python3 and PyYAML, e.g. anaconda3 (empty: system python3)
: "${MPI_LAUNCHER:=}"       # the rank count is appended. Empty: "mpirun -np"; e.g. "srun -n"
# export ROOF_RUNS_DIR=/scratch/$USER/roof_runs    # uncomment to write the runs on /scratch
# ---------------------------------------------------------------------------

if [ -n "${ROOF_ENV_LOADED:-}" ]; then
    return 0
fi
MPI_LAUNCHER=${MPI_LAUNCHER:-mpirun -np}
export OPENFOAM_BASHRC OPENFOAM_MODULE PCM_SIF PCM_SIF_BASHRC PYTHON_MODULE MPI_LAUNCHER

# Run "$@" without -e/-u: OpenFOAM's bashrc and module systems are not written for them.
_roof_relaxed() {
    local _old=$- _status=0
    set +eu
    "$@" || _status=$?
    case $_old in *e*) set -e ;; esac
    case $_old in *u*) set -u ;; esac
    return $_status
}

# Source a file with NO positional arguments: OpenFOAM's bashrc reads its arguments
# as name=value settings and as further files to source.
_roof_source() {
    local _file=$1
    set --
    # shellcheck disable=SC1090
    source "$_file"
}

_roof_module() {
    local _names=("$@") f
    if ! type module >/dev/null 2>&1; then
        for f in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh /usr/share/lmod/lmod/init/bash; do
            if [ -f "$f" ]; then
                _roof_source "$f"
                break
            fi
        done
    fi
    module load "${_names[@]}"
}

if [ -n "$PYTHON_MODULE" ]; then
    # shellcheck disable=SC2086  # several modules, space-separated, are loaded in order
    _roof_relaxed _roof_module $PYTHON_MODULE
fi

if [ -n "$OPENFOAM_BASHRC" ]; then
    _roof_relaxed _roof_source "$OPENFOAM_BASHRC"
    foam() { "$@"; }
elif [ -n "$OPENFOAM_MODULE" ]; then
    # shellcheck disable=SC2086
    _roof_relaxed _roof_module $OPENFOAM_MODULE
    # some modules only define the bashrc location
    if [ -z "${WM_PROJECT_DIR:-}" ]; then
        for _roof_f in "${FOAM_BASH:-}" "${FOAM_BASHRC:-}"; do
            if [ -n "$_roof_f" ] && [ -f "$_roof_f" ]; then
                _roof_relaxed _roof_source "$_roof_f"
                break
            fi
        done
    fi
    foam() { "$@"; }
elif [ -n "$PCM_SIF" ]; then
    ROOF_APPTAINER=$(command -v apptainer || command -v singularity || true)
    if [ -z "$ROOF_APPTAINER" ]; then
        echo "slurm/env.sh: PCM_SIF is set, but neither apptainer nor singularity is on PATH" >&2
        return 1
    fi
    export ROOF_APPTAINER
    # Inside the image: source its OpenFOAM (with the command hidden from it: the
    # bashrc reads its arguments as settings and files to source), and hide the job's
    # SLURM/PMI variables so that the container's own mpirun starts the ranks on this
    # node (the image has no srun).
    foam() {
        # shellcheck disable=SC2016
        "$ROOF_APPTAINER" exec ${ROOF_BIND:+--bind "$ROOF_BIND"} "$PCM_SIF" bash -c '
            set +u
            bashrc=$0
            cmd=("$@")
            set --
            source "$bashrc" >/dev/null 2>&1
            set -- "${cmd[@]}"
            for v in $(compgen -e); do
                case $v in SLURM_*|PMI_*|PMIX_*) unset "$v" ;; esac
            done
            exec "$@"' "$PCM_SIF_BASHRC" "$@"
    }
else
    foam() { "$@"; }        # OpenFOAM already in the environment (roof_check_env tells)
fi
export -f foam

# Fails (non-zero) unless chtMultiRegionFoam of an openfoam.com release and
# python3 with PyYAML are available.
roof_check_env() {
    local v
    v=$(foam bash -c 'command -v chtMultiRegionFoam >/dev/null 2>&1 && echo "${WM_PROJECT_VERSION:-unknown}"') || true
    if [ -z "$v" ]; then
        echo "OpenFOAM (chtMultiRegionFoam) not found: set OPENFOAM_BASHRC, OPENFOAM_MODULE" \
             "or PCM_SIF in slurm/env.sh" >&2
        return 1
    fi
    case $v in
        v[0-9][0-9][0-9][0-9]*)
            if [ "${v:1:4}" -lt 2306 ]; then
                echo "warning: OpenFOAM $v is older than v2306, the oldest release tested" >&2
            fi ;;
        *)
            echo "OpenFOAM '$v' is not an openfoam.com (ESI) release (vYYMM): not supported" >&2
            return 1 ;;
    esac
    echo "OpenFOAM $v: $(foam bash -c 'command -v chtMultiRegionFoam')"
    if ! python3 -c 'import yaml' >/dev/null 2>&1; then
        echo "python3 cannot import yaml: pip install --user pyyaml, or set PYTHON_MODULE" >&2
        return 1
    fi
    echo "Python: $(python3 --version 2>&1), PyYAML found; MPI launcher: $MPI_LAUNCHER <n>"
}

export ROOF_ENV_LOADED=1
