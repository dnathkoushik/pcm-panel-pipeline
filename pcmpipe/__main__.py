"""Command line:  python -m pcmpipe <step> --study <name> [options]

Steps: list | cad | mesh | setup | run | results | post | plots | all
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import traceback

from . import cad, case_setup, config, mesh, results, run

POST_DIR = os.path.join(config.ROOT, "post")


def pv_post(case: dict, dry: bool):
    exe = case["tools"]["pvpython"]
    times = [str(t) for t in case["post"]["contour_times_min"]]
    for script, extra in (("pv_contours.py", ["--times", *times]), ("pv_extract_transient.py", [])):
        cmd = [exe, os.path.join(POST_DIR, script), case["dir"], *extra]
        print("  " + " ".join(cmd))
        if not dry:
            subprocess.run(cmd, check=True)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="pcmpipe", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["list", "cad", "mesh", "setup", "run", "results",
                                     "post", "plots", "all"])
    ap.add_argument("--study", default="validation", help="study name from config/studies.yaml")
    ap.add_argument("--case", default=None, help="only cases whose id contains this text")
    ap.add_argument("--np", type=int, default=None, help="MPI ranks (overrides run.n_procs)")
    ap.add_argument("--np-mesh", type=int, default=1, help="MPI ranks for snappyHexMesh")
    ap.add_argument("--dry-run", action="store_true", help="write files, print commands only")
    ap.add_argument("--keep-going", action="store_true", help="continue after a failed case")
    a = ap.parse_args(argv)

    if a.step == "list":
        for s in config.study_names():
            print(f"{s}:")
            for c in config.expand_study(s):
                print(f"    {c['id']}")
        return 0

    cases = [c for c in config.expand_study(a.study) if not a.case or a.case in c["id"]]
    if not cases:
        sys.exit("no matching cases")

    if a.step == "plots":
        return subprocess.call([sys.executable, os.path.join(POST_DIR, "plot_results.py"),
                                "--study", a.study])

    steps = ["cad", "mesh", "setup", "run", "results", "post"] if a.step == "all" else [a.step]
    failed = []
    for c in cases:
        print(f"=== {c['id']} ({a.study}) ===")
        config.save_case(c)
        try:
            for s in steps:
                print(f"- {s}")
                if s == "cad":
                    cad.build(c, dry=a.dry_run)
                elif s == "mesh":
                    mesh.build(c, n_procs=a.np_mesh, dry=a.dry_run)
                elif s == "setup":
                    case_setup.build(c)
                elif s == "run":
                    run.run_case(c, n_procs=a.np, dry=a.dry_run)
                elif s == "results" and not a.dry_run:
                    m = results.metrics(c["dir"], c["run"]["melt_threshold_100"])
                    print("  " + json.dumps({k: m.get(k) for k in (
                        "t80_min", "t100_min", "average_melting_rate_g_min",
                        "normalised_rate_g_min_L")}))
                elif s == "post":
                    pv_post(c, a.dry_run)
        except Exception as e:   # noqa: BLE001
            failed.append(c["id"])
            traceback.print_exc()
            if not a.keep_going:
                return 1
    if a.step == "all" and not a.dry_run:
        subprocess.call([sys.executable, os.path.join(POST_DIR, "plot_results.py"),
                         "--study", a.study])
    if failed:
        print(f"failed cases: {failed}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
