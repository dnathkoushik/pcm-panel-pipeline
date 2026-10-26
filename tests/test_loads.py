"""Roof thermal loads (pcmpipe/loads.py) and material / loads / --set configuration."""
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pcmpipe import config, loads  # noqa: E402

SPEC = config._load("loads.yaml")["kharagpur_roof"]


class SolarAndWeather(unittest.TestCase):
    def test_kharagpur_sun_in_may(self):
        rise, sset, noon = loads.solar_times(135, 22.33, 87.31, 5.5)
        self.assertAlmostEqual(noon, 11.65, delta=0.1)          # IST is set for 82.5 E
        self.assertAlmostEqual(sset - rise, 13.2, delta=0.2)    # day length mid-May, 22 N
        self.assertAlmostEqual(rise, 5.0, delta=0.15)

    def test_synthetic_day(self):
        h = np.arange(0.0, 24.0, 0.01)
        T, I = loads.weather(SPEC, h)
        w = SPEC["weather"]
        self.assertAlmostEqual(T.max(), w["T_max_C"], places=2)
        self.assertAlmostEqual(T.min(), w["T_min_C"], places=2)
        self.assertAlmostEqual(h[T.argmax()], w["t_Tmax_h"], delta=0.02)
        self.assertAlmostEqual(I.sum() * 0.01 / 1000.0, w["daily_global_horizontal_kWh_m2"],
                               delta=0.01)
        self.assertEqual(I[h < 4.5].max(), 0.0)

    def test_epw_reader(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "x.epw")
            with open(path, "w") as f:
                f.write("HEADER\n" * 8)
                for day, peak in ((1, 35.0), (2, 39.0)):
                    for hr in range(1, 25):
                        T = peak - 6 + 6 * np.sin(np.pi * (hr - 6) / 12) if 6 <= hr <= 18 else peak - 8
                        row = [2001, 5, day, hr, 60, "x", f"{T:.1f}", 20, 60, 100000, 0, 0, 400,
                               max(0.0, 900 * np.sin(np.pi * (hr - 6) / 12)) if 6 < hr < 18 else 0]
                        f.write(",".join(str(v) for v in row) + "\n")
            day, hrs, T, I = loads.read_epw_day(path, 5)
            self.assertEqual(day, 2)
            self.assertEqual(hrs[0], 0.5)
            self.assertAlmostEqual(T.max(), 39.0, places=1)


class Slab(unittest.TestCase):
    LAYERS = [{"thickness": 0.125, "k": 1.4, "rho": 2100.0, "cp": 840.0}]

    def test_steady_state_matches_resistances(self):
        h = np.arange(0.0, 24.0, 0.5)
        Ts, q = loads.slab_periodic(h, np.full_like(h, 50.0), np.full_like(h, 30.0),
                                    self.LAYERS, 17.0, 6.13)
        U = 1.0 / (1 / 17.0 + 0.125 / 1.4 + 1 / 6.13)
        np.testing.assert_allclose(q, U * 20.0, rtol=1e-3)
        np.testing.assert_allclose(Ts, 30.0 + U * 20.0 / 6.13, atol=1e-2)

    def test_mass_delays_and_damps_the_ceiling(self):
        p = loads.profile(SPEC)
        clock = np.mod(p["clock_h"], 24.0)
        lag = clock[p["T_ceiling_C"].argmax()] - clock[p["T_solair_C"].argmax()]
        self.assertTrue(2.0 < lag < 6.0, lag)
        swing_in = p["T_ceiling_C"].max() - p["T_ceiling_C"].min()
        swing_out = p["T_solair_C"].max() - p["T_solair_C"].min()
        self.assertLess(swing_in, 0.6 * swing_out)
        self.assertTrue(np.all(p["T_ceiling_C"] < p["T_solair_C"].max()))
        self.assertEqual(p["time_s"][0], 0.0)
        self.assertAlmostEqual(p["time_s"][-1], 86400.0)


class Configuration(unittest.TestCase):
    def test_octadecane_roof_cases(self):
        cases = config.expand_study("octadecane_roof")
        self.assertEqual([c["id"] for c in cases], ["pyr_V250_r50_t1_g0_kgp_c18",
                                                    "tet_V250_r50_t1_g0_kgp_c18",
                                                    "hem_V250_r50_t1_g0_kgp_c18"])
        c = cases[0]
        pcm, ic = c["physics"]["pcm"], c["physics"]["initial_conditions"]
        self.assertEqual((pcm["T_solidus"], pcm["T_liquidus"]), (300.63, 301.63))
        self.assertEqual(pcm["latent_heat"], 236980.0)
        self.assertEqual(pcm["A_mush"], 1.0e6)                     # paper value kept
        self.assertEqual(c["physics"]["solver"]["nOuterCorrectors"], 3)
        self.assertAlmostEqual(ic["T_pcm"], 297.15)
        self.assertAlmostEqual(ic["T_air"] - 273.15, loads.room_temperature_C(SPEC, [6.0])[0])
        screen = config.expand_study("octadecane_roof", "screening")[0]
        self.assertEqual(screen["id"], "pyr_V250_r50_t1_g0_kgp_c18_screen")

    def test_paper_cases_unchanged(self):
        c = config.expand_study("validation")[0]
        self.assertEqual(c["id"], "pyr_V250_r50_t1_g0_dT20")
        self.assertEqual(c["physics"]["pcm"]["latent_heat"], 115000.0)
        self.assertEqual(c["physics"]["initial_conditions"]["T_air"], 298.15)
        self.assertIsNone(c["loads_spec"])

    def test_set_overrides_and_wall_layers(self):
        sets = config.parse_sets(["mesh.wall_layer_passes={air: 1, pcm: 2}", "run.n_procs=4",
                                  "pcm_material=n_octadecane"])
        c = config.expand_study("shapes", "screening", sets)[0]
        self.assertEqual(c["mesh"]["active"]["wall_layer_passes"], {"air": 1, "pcm": 2})
        self.assertEqual(c["run"]["n_procs"], 4)
        self.assertTrue(c["id"].endswith("_c18_screen"))
        for name in ("smoke", "screen", "coarse", "medium", "fine"):
            c = config.base_case()
            self.assertEqual(c["mesh"]["presets"][name]["wall_layer_passes"], {"air": 0, "pcm": 0})
        with self.assertRaises(ValueError):
            config.parse_sets(["no_equals_sign"])

    def test_foam_table(self):
        s = loads.foam_table([0, 600], [300.0, 301.25])
        self.assertEqual(s, "table ((0 300.000) (600 301.250))")


if __name__ == "__main__":
    unittest.main(verbosity=2)
