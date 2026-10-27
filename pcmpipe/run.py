"""Step 4a -- execution: decompose, run chtMultiRegionFoam, reconstruct."""
from __future__ import annotations

import glob
import os
import shutil

from . import case_setup, config, foam, mesh


def build_library(dry: bool = False):
    """Compile libpcmEnthalpyPorosity.so once per OpenFOAM installation."""
    lib = os.path.join(os.environ.get("FOAM_USER_LIBBIN", ""), "libpcmEnthalpyPorosity.so")
    if os.path.exists(lib) and not dry:
        return
    foam.run(["./Allwmake"], os.path.join(config.ROOT, "openfoam"), log="wmake", dry=dry)


def run_case(case: dict, n_procs: int | None = None, dry: bool = False):
    cdir = os.path.join(case["dir"], "case")
    n = n_procs or case["run"]["n_procs"]
    if not dry:
        foam.check_openfoam()
        build_library()
    if n > 1:
        case_setup.write_decomposition(cdir, n, case["run"]["decomposition"], mesh.regions(case))
        # -latestTime: a fresh case has only 0/; a restarted one needs only the time it
        # continues from (controlDict: startFrom latestTime)
        foam.run(["decomposePar", "-allRegions", "-force", "-latestTime"], cdir, dry=dry)
        foam.run([*foam.tool(case, "mpirun"), "-np", str(n), "chtMultiRegionFoam", "-parallel"],
                 cdir, dry=dry)
        foam.run(["reconstructPar", "-allRegions", "-newTimes"], cdir, dry=dry)
        if not dry:
            for d in glob.glob(os.path.join(cdir, "processor*")):
                shutil.rmtree(d)
    else:
        foam.run(["chtMultiRegionFoam"], cdir, dry=dry)
