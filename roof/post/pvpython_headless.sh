#!/bin/bash
# pvpython without a display: if DISPLAY is unset and Xvfb exists, start a private
# Xvfb for this call (xvfb-run can hang waiting for Xvfb's ready signal inside
# containers). An offscreen ParaView build (EGL/OSMesa) needs none of this: set
# PVPYTHON to it and it is used as is.
#
#   bash roof/post/pvpython_headless.sh roof/post/pv_extract.py roof/runs/<case>
set -euo pipefail
read -r -a pv <<< "${PVPYTHON:-pvpython}"
if [ -n "${DISPLAY:-}" ] || ! command -v Xvfb >/dev/null 2>&1; then
    exec "${pv[@]}" "$@"
fi
for n in $(seq 90 140); do                       # a free display number
    [ -e "/tmp/.X11-unix/X$n" ] || [ -e "/tmp/.X$n-lock" ] || break
done
Xvfb ":$n" -screen 0 1920x1080x24 -nolisten tcp >/dev/null 2>&1 &
xpid=$!
trap 'kill "$xpid" 2>/dev/null || true' EXIT
for _ in $(seq 100); do
    [ -e "/tmp/.X11-unix/X$n" ] && break
    sleep 0.1
done
DISPLAY=":$n" "${pv[@]}" "$@"
