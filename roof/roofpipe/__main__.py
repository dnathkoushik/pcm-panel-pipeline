"""python3 -m roofpipe <command>   (run from the roof/ directory)

Phase 1
  baseline   write the pure concrete slab case, its boundary tables and the 1D
             reference             [--name NAME] [--np N] [--set KEY=VALUE ...] [--force]
  check      gates on a run: complete, periodic, 1D reference (pure slab), energy
                                   RUN (name under the runs dir, or a path)
Phase 2
  cases      list the PCM cases (shape x orientation x height) and their volumes
  pcm        write PCM case(s)     --case ID [--case ID ...] | --all   [--np N] [--force]
  stl        PCM body surface without FreeCAD                     RUN
  finalize   after meshing: PCM density so the meshed mass is exact   RUN
"""
from __future__ import annotations

import argparse
import sys

from . import config


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="python3 -m roofpipe",
                                description="Roof-embedded PCM study: OpenFOAM case set-up")
    p.add_argument("--config", default=None, help="configuration file (default: config/roof.yaml)")
    sub = p.add_subparsers(dest="cmd", required=True)

    def with_set(sp):
        sp.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                        help="override a configuration value, e.g. --set indoor.mode=specimen")

    b = sub.add_parser("baseline", help="Phase 1: write the pure concrete slab case")
    b.add_argument("--name", default="baseline_concrete", help="run folder name")
    b.add_argument("--np", type=int, default=None, help="MPI ranks (default: run.n_procs)")
    b.add_argument("--force", action="store_true", help="rewrite an existing case")
    with_set(b)

    c = sub.add_parser("check", help="gates and daily metrics of a run")
    c.add_argument("run", help="run name (under the runs dir) or path")
    c.add_argument("--no-plot", action="store_true", help="skip the quick-look PNG")

    ls = sub.add_parser("cases", help="Phase 2: list the PCM cases")
    with_set(ls)

    pc = sub.add_parser("pcm", help="Phase 2: write PCM case(s)")
    pc.add_argument("--case", action="append", default=[], metavar="ID",
                    help="case id from `cases`, e.g. pyr_N_h070 (repeatable)")
    pc.add_argument("--all", action="store_true", help="every case of pcm_study")
    pc.add_argument("--np", type=int, default=None, help="MPI ranks (default: run.n_procs)")
    pc.add_argument("--force", action="store_true", help="rewrite existing cases")
    with_set(pc)

    st = sub.add_parser("stl", help="Phase 2: PCM body STL without FreeCAD")
    st.add_argument("run")
    fi = sub.add_parser("finalize", help="Phase 2: PCM density from the meshed volume")
    fi.add_argument("run")

    a = p.parse_args(argv)
    if a.cmd == "baseline":
        from . import case
        cfg = config.load(a.config, a.set)
        n = a.np or int(cfg["run"]["n_procs"])
        if n < 1:
            raise SystemExit("--np must be >= 1")
        case.build_baseline(cfg, a.name, n, a.force)
        return 0
    if a.cmd == "check":
        from . import check
        return check.main(config.run_dir(a.run), plot=not a.no_plot)
    if a.cmd in ("cases", "pcm"):
        from . import geometry
        cfg = config.load(a.config, a.set)
        bodies = geometry.study_cases(cfg)
        if a.cmd == "cases":
            for bd in bodies:
                print(f"{bd.case_id}  {geometry.SHAPES[bd.shape]:<11} {bd.orientation}  "
                      f"h {bd.height * 1e3:5.1f} mm  V {bd.volume * 1e6:7.1f} mL  "
                      f"z {min(bd.z_base, bd.z_apex) * 1e3:5.1f}..{max(bd.z_base, bd.z_apex) * 1e3:5.1f} mm  "
                      f"side margin {bd.margin_xy * 1e3:5.1f} mm")
            print(f"{len(bodies)} cases")
            return 0
        from . import pcm_case
        want = {bd.case_id: bd for bd in bodies}
        ids = list(want) if a.all else a.case
        unknown = [i for i in ids if i not in want]
        if unknown or not ids:
            raise SystemExit(f"unknown case id(s) {unknown}; `python3 -m roofpipe cases` lists them"
                             if unknown else "give --case ID or --all")
        n = a.np or int(cfg["run"]["n_procs"])
        for i in ids:
            pcm_case.build(cfg, want[i], n, a.force)
        return 0
    if a.cmd == "stl":
        from . import stl
        stl.write(config.run_dir(a.run))
        return 0
    if a.cmd == "finalize":
        from . import pcm_case
        pcm_case.finalize(config.run_dir(a.run))
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
