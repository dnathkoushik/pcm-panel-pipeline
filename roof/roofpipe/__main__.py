"""python3 -m roofpipe <command>   (run from the roof/ directory)

  baseline   Phase 1: write the pure concrete slab case, its boundary tables and the
             1D reference           [--name NAME] [--np N] [--set KEY=VALUE ...] [--force]
  check      Phase 1 gates on a finished (or partly finished) run: periodicity, 1D
             reference, energy balance              RUN (name under the runs dir, or path)
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

    b = sub.add_parser("baseline", help="Phase 1: write the pure concrete slab case")
    b.add_argument("--name", default="baseline_concrete", help="run folder name")
    b.add_argument("--np", type=int, default=None, help="MPI ranks (default: run.n_procs)")
    b.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                   help="override a configuration value, e.g. --set indoor.mode=specimen")
    b.add_argument("--force", action="store_true", help="rewrite an existing case")

    c = sub.add_parser("check", help="Phase 1 gates on a run")
    c.add_argument("run", help="run name (under the runs dir) or path")
    c.add_argument("--no-plot", action="store_true", help="skip the quick-look PNG")

    a = p.parse_args(argv)
    if a.cmd == "baseline":
        from . import case
        cfg = config.load(a.config, a.set)
        n = a.np or int(cfg["run"]["n_procs"])
        if n < 1:
            raise SystemExit("--np must be >= 1")
        case.build_baseline(cfg, a.name, n, a.force)
        return 0
    from . import check
    return check.main(config.run_dir(a.run), plot=not a.no_plot)


if __name__ == "__main__":
    sys.exit(main())
