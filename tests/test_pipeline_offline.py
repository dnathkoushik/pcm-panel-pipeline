"""Offline end-to-end test of the Python side of the pipeline (no OpenFOAM needed).

Writes every case (dictionaries + 0/constant/system) for several studies into a
temporary runs folder, fakes the function-object output a chtMultiRegionFoam run
would leave behind, then runs `results` and the plotting script and checks the
outputs. Runs on Linux, macOS and Windows.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T100_MIN = {"pyr": 292.0, "tet": 353.0, "hem": 320.0}     # paper values for pyr/tet


def fake_function_objects(case_dir, t100_min):
    """postProcessing files with the same layout/headers OpenFOAM writes."""
    T100 = 60.0 * t100_min
    t = np.arange(0.0, 1.1 * T100, 30.0)
    s = np.clip(t / T100, 0.0, 1.0)
    g = np.clip(1.25 * s - 0.25 * s**3, 0.0, 1.0)
    Tm = 278.15 + 9.0 * (1.0 - np.exp(-t / 600.0)) + 6.0 * s
    base = os.path.join(case_dir, "case", "postProcessing", "pcm")
    os.makedirs(os.path.join(base, "pcmAverages", "0"))
    with open(os.path.join(base, "pcmAverages", "0", "volFieldValue.dat"), "w") as f:
        f.write("# Region type : all\n# Time\tvolAverage(liquidFraction)\tvolAverage(T)\n")
        for row in zip(t, g, Tm):
            f.write("%g\t%.6f\t%.4f\n" % row)
    os.makedirs(os.path.join(base, "pcmProbes", "0"))
    with open(os.path.join(base, "pcmProbes", "0", "T"), "w") as f:
        f.write("# Probe 0 (0 0 0)\n# Probe 1 (0 0 0)\n# Probe 2 (0 0 0)\n#  Time  0  1  2\n")
        for a, c in zip(t, Tm):
            f.write("%g %.3f %.3f %.3f\n" % (a, c + 1.5, c + 0.5, c - 0.5))


class OfflinePipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.env = dict(os.environ, PCMPIPE_RUNS_DIR=cls.tmp.name, MPLBACKEND="Agg")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def pcmpipe(self, *args):
        r = subprocess.run([sys.executable, "-m", "pcmpipe", *args], cwd=ROOT, env=self.env,
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, msg=f"pcmpipe {' '.join(args)}\n{r.stdout}\n{r.stderr}")
        return r.stdout

    def test_1_dry_run_every_study(self):
        out = self.pcmpipe("list")
        studies = [line[:-1] for line in out.splitlines() if line and not line.startswith(" ")]
        self.assertIn("validation", studies)
        for study in studies:
            self.pcmpipe("mesh", "--study", study, "--dry-run")
            self.pcmpipe("setup", "--study", study)
        case = os.path.join(self.tmp.name, "pyr_V250_r50_t1_g0_dT20", "case")
        for rel in ("system/snappyHexMeshDict", "system/controlDict", "constant/pcm/fvOptions",
                    "constant/regionProperties", "0/air/T", "0/pcm/T", "0/shell/T"):
            self.assertTrue(os.path.isfile(os.path.join(case, rel)), rel)

    def test_1b_screening_mode(self):
        self.pcmpipe("mesh", "--study", "shapes", "--mode", "screening", "--dry-run")
        self.pcmpipe("setup", "--study", "shapes", "--mode", "screening")
        case = os.path.join(self.tmp.name, "hem_V250_r50_t1_g0_dT20_screen", "case")
        with open(os.path.join(case, "constant", "regionProperties")) as f:
            self.assertIn("solid ()", f.read())
        self.assertFalse(os.path.exists(os.path.join(case, "0", "shell")))
        with open(os.path.join(case, "0", "pcm", "T")) as f:
            T = f.read()
        self.assertIn("pcm_to_air", T)
        self.assertIn("thicknessLayers (0.001)", T)
        with open(os.path.join(case, "system", "snappyHexMeshDict")) as f:
            snappy = f.read()
        self.assertNotIn("panel_outer", snappy)
        self.assertIn("pcm_cavity_air", snappy)

    def test_2_results_and_plots(self):
        for study in ("validation", "volume"):
            self.pcmpipe("setup", "--study", study)
        for d in os.listdir(self.tmp.name):
            if "_r50_t1_g0_dT20" in d and d.endswith("dT20"):
                fake_function_objects(os.path.join(self.tmp.name, d), T100_MIN[d[:3]] * (
                    float(d.split("_V")[1].split("_")[0]) / 250.0) ** 0.6)
        self.pcmpipe("results", "--study", "validation")
        with open(os.path.join(self.tmp.name, "pyr_V250_r50_t1_g0_dT20", "metrics.json")) as f:
            m = json.load(f)
        self.assertAlmostEqual(m["t100_min"], 292.0, delta=1.0)
        self.assertGreater(m["t100_min"], m["t80_min"])

        for study in ("validation", "volume"):
            r = subprocess.run([sys.executable, os.path.join(ROOT, "post", "plot_results.py"),
                                "--study", study], cwd=ROOT, env=self.env, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, msg=r.stdout + r.stderr)
        figs = os.path.join(self.tmp.name, "_figures")
        for rel in ("validation/liquid_fraction_and_rate.png", "validation/probe_temperatures.png",
                    "volume/melting_time_bars.png", "volume/normalised_rate_paper.png",
                    "volume/normalised_rate_physical.png", "volume/summary.csv"):
            self.assertTrue(os.path.isfile(os.path.join(figs, rel)), rel)


if __name__ == "__main__":
    unittest.main(verbosity=2)
