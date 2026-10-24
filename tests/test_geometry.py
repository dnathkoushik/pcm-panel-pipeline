"""Checks the geometry solver against numbers published in Yang et al. (2023)."""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pcmpipe.geometry import build_panel, polygon_area  # noqa: E402

BASE = dict(shape="pyramid", volume_ml=250, protrusion_ratio=0.5, shell_thickness_mm=1.0,
            gap_mm=0.0, n_x=2, n_y=2, pyramid_edge_angle_deg=43.0,
            tetrahedron_edge_angle_deg=54.7356, hemisphere_fill=0.95)


def panel(**kw):
    d = dict(BASE)
    d.update(kw)
    return build_panel(d)


class PaperGeometry(unittest.TestCase):
    def test_envelope_volume(self):
        for shape in ("pyramid", "tetrahedron", "hemisphere"):
            for V in (250, 500, 750):
                g = panel(shape=shape, volume_ml=V)
                self.assertAlmostEqual(g.envelope_volume_mm3 / 1000, V, places=6)

    def test_effective_area_table3(self):
        # Table 3: 232 cm2 (pyramid), 282 cm2 (tetrahedron)
        self.assertAlmostEqual(panel(shape="pyramid").exposed_area_mm2 / 100, 232, delta=232 * 0.03)
        self.assertAlmostEqual(panel(shape="tetrahedron").exposed_area_mm2 / 100, 282, delta=282 * 0.01)

    def test_tetra_footprint_fig4(self):
        # Fig. 4b: 12.8 cm x 6.4 cm parallelogram
        g = panel(shape="tetrahedron")
        xs = [p[0] for p in g.footprint]
        a = math.dist(g.footprint[0], g.footprint[3])
        self.assertAlmostEqual(a / 10, 6.4, delta=0.05)
        self.assertAlmostEqual(math.dist(g.footprint[0], g.footprint[1]) / 10, 12.8, delta=0.1)
        self.assertGreater(max(xs) - min(xs), 0)

    def test_edge_angles_sec45(self):
        # Sec. 4.5: pyramid 29.2/43.0/52.6 deg, tetrahedron 40.3/54.7/63.2 deg
        for r, pyr, tet in ((0.3, 29.2, 40.3), (0.5, 43.0, 54.7), (0.7, 52.6, 63.2)):
            self.assertAlmostEqual(panel(protrusion_ratio=r).edge_angle_deg, pyr, delta=0.1)
            self.assertAlmostEqual(panel(shape="tetrahedron", protrusion_ratio=r).edge_angle_deg,
                                   tet, delta=0.1)

    def test_area_per_volume_sec47(self):
        # Sec. 4.7, tetrahedron: 1.13, 0.895, 0.782 cm2/mL (self-similar scaling)
        for V, ref in ((250, 1.13), (500, 0.895), (750, 0.782)):
            g = panel(shape="tetrahedron", volume_ml=V)
            self.assertAlmostEqual(g.exposed_area_mm2 / 100 / V, ref, delta=0.01)

    def test_ratio_keeps_footprint(self):
        a = polygon_area(panel(protrusion_ratio=0.3).footprint)
        b = polygon_area(panel(protrusion_ratio=0.7).footprint)
        self.assertAlmostEqual(a, b, places=6)

    def test_gap_preserves_volume_and_grows_area(self):
        for shape in ("pyramid", "tetrahedron", "hemisphere"):
            g0, g5 = panel(shape=shape), panel(shape=shape, gap_mm=5)
            self.assertAlmostEqual(g0.envelope_volume_mm3, g5.envelope_volume_mm3, places=3)
            self.assertGreater(g5.exposed_area_mm2, g0.exposed_area_mm2)

    def test_probes_inside_protrusion(self):
        g = panel()
        zs = [p[2] for p in g.probes.values()]
        self.assertTrue(all(g.z_tip * 1e-3 < z <= -g.slab_thickness * 1e-3 + 1e-12 for z in zs))


if __name__ == "__main__":
    unittest.main(verbosity=2)
