# Complete toolchain for the pipeline: OpenFOAM (openfoam.com), FreeCAD, ParaView
# (headless via Xvfb) and the Python packages. Works on Linux, macOS (Intel and
# Apple Silicon, via amd64 emulation) and Windows (Docker Desktop).
#
#   docker build --platform linux/amd64 -t pcm-panel-pipeline .
#   scripts/docker_run.sh python3 -m pcmpipe all --study smoke --np 2
#
FROM ubuntu:24.04

ARG OF_VERSION=2406
ENV OF_VERSION=${OF_VERSION} \
    DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    MPLBACKEND=Agg
# pvpython needs an X display for screenshots; Xvfb provides a virtual one
ENV PCMPIPE_PVPYTHON="xvfb-run -a -s '-screen 0 1920x1080x24' pvpython"
# Open MPI inside containers: no CMA single-copy, allow more ranks than cores
ENV OMPI_MCA_btl_vader_single_copy_mechanism=none \
    OMPI_MCA_rmaps_base_oversubscribe=1

# FreeCAD is not in the Ubuntu 24.04 archive: use the official FreeCAD PPA
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl gnupg software-properties-common \
 && add-apt-repository -y ppa:freecad-maintainers/freecad-stable \
 && curl -fsSL https://dl.openfoam.com/add-debian-repo.sh | bash \
 && apt-get update \
 && apt-get install -y --no-install-recommends \
        "openfoam${OF_VERSION}-default" \
        freecad paraview xvfb xauth \
        python3 python3-numpy python3-matplotlib python3-yaml \
        make \
 && rm -rf /var/lib/apt/lists/*

# Let any host uid (docker run --user) get a passwd entry; Open MPI and
# OpenFOAM look the user up.
RUN chmod 666 /etc/passwd
COPY docker/entrypoint.sh /usr/local/bin/pcm-entrypoint
RUN chmod 755 /usr/local/bin/pcm-entrypoint

USER ubuntu
ENV HOME=/home/ubuntu USER=ubuntu
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/pcm-entrypoint"]
CMD ["bash"]
