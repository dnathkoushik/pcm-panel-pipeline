# Shortcuts for Linux / macOS. Native targets need a loaded OpenFOAM environment;
# docker-* targets need only Docker.  `make help` lists them.
PYTHON ?= python3
NP ?= 4
STUDY ?= validation
# MODE: paper (resolved shell, Sec. 3) | screening (thin-shell resistance, coarse mesh)
MODE ?= paper
# SMOKE: smoke (paper set-up) | smoke_roof (n-octadecane, Kharagpur roof loads)
SMOKE ?= smoke
IMAGE ?= pcm-panel-pipeline

.PHONY: help test lib smoke study plots docker-build docker-shell docker-test docker-smoke docker-study clean

help:  ## list targets
	@grep -E '^[a-z-]+:.*## ' $(MAKEFILE_LIST) | sed -E 's/:.*## /	/' | expand -t 16

test:  ## Python tests: geometry, latent-heat scheme, offline pipeline (no OpenFOAM)
	$(PYTHON) -m unittest discover -s tests -v

lib:  ## compile the PCM fvOption library (OpenFOAM env required)
	./openfoam/Allwmake

smoke: lib  ## short plumbing run: CAD, mesh, solve, post, plots (SMOKE=smoke_roof: roof loads)
	$(PYTHON) -m pcmpipe all --study $(SMOKE) --mode $(MODE) --np $(NP) --np-mesh $(NP)

study: lib  ## full study, e.g. make study STUDY=volume MODE=screening NP=8
	$(PYTHON) -m pcmpipe all --study $(STUDY) --mode $(MODE) --np $(NP) --np-mesh $(NP) --keep-going

plots:  ## re-plot a study from existing results
	$(PYTHON) -m pcmpipe plots --study $(STUDY) --mode $(MODE)

docker-build:  ## build the container (OpenFOAM + FreeCAD + ParaView)
	docker build --platform linux/amd64 -t $(IMAGE) .

docker-shell:  ## interactive shell in the container, repo mounted at /work
	scripts/docker_run.sh

docker-test:  ## Python tests inside the container
	scripts/docker_run.sh make test

docker-smoke:  ## smoke run inside the container
	scripts/docker_run.sh make smoke NP=$(NP) MODE=$(MODE) SMOKE=$(SMOKE)

docker-study:  ## full study inside the container, e.g. make docker-study STUDY=validation NP=8
	scripts/docker_run.sh make study STUDY=$(STUDY) NP=$(NP) MODE=$(MODE)

clean:  ## remove all results and the container home folder
	rm -rf runs .docker-home
