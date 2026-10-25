#!/bin/bash
# Container entry point: give the current uid a home and a passwd entry, load
# OpenFOAM, build the PCM fvOption library if needed, then run the command.
export HOME="${HOME:-/tmp/home}"
mkdir -p "$HOME" 2>/dev/null || { export HOME=/tmp/home; mkdir -p "$HOME"; }
if ! whoami >/dev/null 2>&1; then
    echo "pcm:x:$(id -u):$(id -g):pcm:${HOME}:/bin/bash" >> /etc/passwd
fi
export USER="$(whoami)"

# shellcheck disable=SC1090
source "/usr/lib/openfoam/openfoam${OF_VERSION}/etc/bashrc"

if [ -x /work/openfoam/Allwmake ] \
   && [ ! -f "$FOAM_USER_LIBBIN/libpcmEnthalpyPorosity.so" ]; then
    echo "[pcm-entrypoint] compiling libpcmEnthalpyPorosity.so ..."
    ( cd /work/openfoam && ./Allwmake > "$HOME/log.wmake" 2>&1 ) \
        || { echo "[pcm-entrypoint] wmake failed:"; tail -40 "$HOME/log.wmake"; exit 1; }
fi
exec "$@"
