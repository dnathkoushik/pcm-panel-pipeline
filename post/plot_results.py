"""Comparative plots in the style of Yang et al. (2023), Figs. 8, 10, 11-13, 15, 17, 18.

    python post/plot_results.py --study delta_T [--rate-convention physical|paper]

Reads runs/<case>/metrics.json and series.csv (written by `pcmpipe results`) and
writes figures to runs/_figures/<study>/.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pcmpipe import config  # noqa: E402
from pcmpipe import results as res  # noqa: E402

SHAPES = ["pyramid", "tetrahedron", "hemisphere"]
COLOR = {"pyramid": "#1f3fbf", "tetrahedron": "#e01010", "hemisphere": "#1a9641"}
LINESTYLE = {"pyramid": "-", "tetrahedron": "--", "hemisphere": ":"}
PARAM_LABEL = {
    "delta_T": ("Initial temperature difference (K)", "{:g}"),
    "panel.shell_thickness_mm": ("Shell thickness (mm)", "{:g}"),
    "panel.protrusion_ratio": ("Proportion (%)", "{:.0f}"),
    "panel.gap_mm": ("Gap (mm)", "{:g}"),
    "panel.volume_ml": ("Module volume (mL)", "{:g}"),
    "mesh.preset": ("Mesh preset", "{}"),
}
REF = json.load(open(os.path.join(config.ROOT, "validation", "paper_reference.json")))

plt.rcParams.update({"font.size": 11, "axes.linewidth": 1.0, "legend.frameon": False,
                     "xtick.direction": "in", "ytick.direction": "in"})


def load(study):
    spec = config._load("studies.yaml")[study]
    params = [k for k in spec.get("vary", {}) if k != "panel.shape"]
    rows = []
    for c in config.expand_study(study):
        mpath = os.path.join(c["dir"], "metrics.json")
        if not os.path.exists(mpath):
            if os.path.exists(os.path.join(c["dir"], "case_meta.json")):
                try:
                    res.metrics(c["dir"], c["run"]["melt_threshold_100"])
                except FileNotFoundError:
                    continue
            else:
                continue
        m = json.load(open(mpath))
        with open(os.path.join(c["dir"], "series.csv")) as f:
            data = np.array([[float(v) for v in r] for r in list(csv.reader(f))[1:]])
        m["series"] = data
        m["param"] = {p: config.get_dotted(c, p) for p in params}
        if params and params[0] == "panel.protrusion_ratio":
            m["param"][params[0]] *= 100
        m["probes_data"] = res.probe_series(os.path.join(c["dir"], "case"))
        rows.append(m)
    return params, rows


def savefig(fig, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}.png")


def fig_melting_bars(rows, param, outdir):
    """Figs. 11/12/13/15/17/18 (left): 80 % solid bars, 100 % hatched bars."""
    label, fmt = PARAM_LABEL.get(param, (param, "{}"))
    vals = sorted({r["param"][param] for r in rows}, key=lambda v: (str(type(v)), v))
    shapes = [s for s in SHAPES if any(r["shape"] == s for r in rows)]
    nb = 2 * len(shapes)
    w = 0.8 / nb
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for k, (frac, hatch) in enumerate((("t80_min", None), ("t100_min", "\\\\\\\\"))):
        for j, s in enumerate(shapes):
            ys = []
            for v in vals:
                hit = [r for r in rows if r["shape"] == s and r["param"][param] == v]
                ys.append(hit[0].get(frac) or np.nan if hit else np.nan)
            pos = np.arange(len(vals)) - 0.4 + w * (k * len(shapes) + j + 0.5)
            ax.bar(pos, ys, w, color="white" if hatch else COLOR[s], edgecolor=COLOR[s],
                   hatch=hatch, label=f"{s.capitalize()}, {'80' if k == 0 else '100'}%")
    ax.set_xticks(range(len(vals)), [fmt.format(v) for v in vals])
    ax.set_xlabel(label)
    ax.set_ylabel("Melting time (min)")
    ax.legend(fontsize=8, loc="upper right")
    savefig(fig, outdir, "melting_time_bars")


def _styles(rows, param):
    """Colour by shape when one value per shape, else by parameter value."""
    vals = sorted({r["param"].get(param) for r in rows}, key=lambda v: (str(type(v)), v)) if param else [None]
    cmap = plt.get_cmap("viridis", max(2, len(vals)))
    fmt = PARAM_LABEL.get(param, (param, "{}"))[1]
    out = []
    for r in rows:
        v = r["param"].get(param) if param else None
        color = COLOR[r["shape"]] if len(vals) == 1 else cmap(vals.index(v))
        label = r["shape"].capitalize() + (f", {fmt.format(v)}" if param else "")
        out.append((r, color, label))
    return out


def fig_liquid_fraction(rows, param, outdir):
    """Fig. 11 (right): liquid fraction curves for each shape / parameter value."""
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for r, color, label in _styles(rows, param):
        ax.plot(r["series"][:, 0], r["series"][:, 1], LINESTYLE[r["shape"]], color=color, label=label)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Liquid fraction")
    ax.set_ylim(0, 1.05)
    ax.set_xlim(left=0)
    ax.legend(fontsize=8, loc="lower right")
    savefig(fig, outdir, "liquid_fraction_curves")


def fig_mean_temperature(rows, param, outdir):
    """Figs. 13/15/17 (right): volume-averaged PCM temperature."""
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for r, color, label in _styles(rows, param):
        ax.plot(r["series"][:, 0], r["series"][:, 2] - 273.15, LINESTYLE[r["shape"]],
                color=color, label=label)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_xlim(left=0)
    ax.legend(fontsize=8, loc="lower right")
    savefig(fig, outdir, "mean_pcm_temperature")


def fig_rate_twin(rows, outdir):
    """Fig. 10: liquid fraction (left) and transient melting rate (right)."""
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax2 = ax.twinx()
    for i, r in enumerate(rows):
        t, g, _, rate = r["series"].T
        c = COLOR[r["shape"]]
        ax.plot(t, g, "-", color=c, label=r["shape"].capitalize())
        ax2.plot(t, rate, "-.", color=c)
        for key, y0 in (("t80_min", 0.83), ("t100_min", 1.03)):
            if r.get(key):
                ax.axvline(r[key], color=c, lw=0.8, ls=":")
                ax.text(r[key], y0 + 0.045 * i, f"t = {r[key]:.0f} min", color=c,
                        ha="right", fontsize=8)
    ax.axhline(0.8, color="k", lw=0.6, ls=":")
    ax.axhline(1.0, color="k", lw=0.6, ls=":")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Liquid fraction")
    ax2.set_ylabel("Melting rate (g min$^{-1}$)")
    ax.set_ylim(0, 1.03 + 0.045 * len(rows) + 0.03)
    ax2.set_ylim(bottom=0)
    ax.set_xlim(left=0)
    ax.legend(loc="center right")
    savefig(fig, outdir, "liquid_fraction_and_rate")


def fig_normalised_rate(rows, param, outdir, convention):
    """Fig. 18 (right): normalised melting rate vs. module volume."""
    key = "normalised_rate_paper_g_min_L" if convention == "paper" else "normalised_rate_g_min_L"
    vals = sorted({r["param"][param] for r in rows})
    shapes = [s for s in SHAPES if any(r["shape"] == s for r in rows)]
    w = 0.8 / len(shapes)
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for j, s in enumerate(shapes):
        ys = []
        for v in vals:
            hit = [r for r in rows if r["shape"] == s and r["param"][param] == v]
            ys.append(hit[0].get(key, np.nan) if hit else np.nan)
        pos = np.arange(len(vals)) - 0.4 + w * (j + 0.5)
        ax.bar(pos, ys, w, color="white", edgecolor=COLOR[s], hatch="////", label=s.capitalize())
        ax.plot(pos, ys, "--", color=COLOR[s], marker="o", ms=3)
        for x, y in zip(pos, ys):
            if np.isfinite(y):
                ax.text(x, y / 2, f"{y:.3f}", color=COLOR[s], ha="center", fontsize=7, rotation=90)
    for s in shapes:   # paper values as markers
        ref = REF["normalised_rate_g_min_L"].get(s)
        if ref and convention == "paper":
            ax.scatter([vals.index(v) for v in (250, 500, 750) if v in vals],
                       [ref[str(v)] for v in (250, 500, 750) if v in vals],
                       marker="x", color=COLOR[s], zorder=5, label=f"{s.capitalize()} (paper)")
    ax.set_xticks(range(len(vals)), [f"{v:g}" for v in vals])
    ax.set_xlabel("Module volume (mL)")
    ax.set_ylabel("Normalised melting rate (g min$^{-1}$ L$^{-1}$)")
    ax.set_title(f"convention: {convention}", fontsize=9)
    ax.legend(fontsize=8)
    savefig(fig, outdir, f"normalised_rate_{convention}")


def fig_probes(rows, outdir):
    """Fig. 8: thermocouple temperatures (CFD) with the paper's melting times."""
    shapes = [r for r in rows if r["probes_data"].size]
    if not shapes:
        return
    fig, axes = plt.subplots(1, len(shapes), figsize=(6.0 * len(shapes), 4.5), squeeze=False)
    for ax, r in zip(axes[0], shapes):
        d = r["probes_data"]
        for i, (name, col) in enumerate((("top", "red"), ("mid", "blue"), ("bot", "green"))):
            if d.shape[1] > i + 1:
                ax.plot(d[:, 0] / 60, d[:, i + 1] - 273.15, color=col, label=f"{name}_CFD")
        ref = REF["melting_time_min"].get(r["shape"])
        if ref:
            for k, ls in (("t80", ":"), ("t100", "--")):
                ax.axvline(ref[k], color="gray", ls=ls, lw=0.8, label=f"paper {k} = {ref[k]} min")
        ax.set_title(r["shape"].capitalize())
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Temperature (Â°C)")
        ax.set_ylim(0, 30)
        ax.legend(fontsize=8, loc="lower right")
    savefig(fig, outdir, "probe_temperatures")


def write_table(rows, params, outdir):
    os.makedirs(outdir, exist_ok=True)
    keys = ["id", "shape", *params, "pcm_mass_g", "exposed_area_cm2", "t80_min", "t100_min",
            "average_melting_rate_g_min", "normalised_rate_g_min_L",
            "normalised_rate_paper_g_min_L", "average_melting_rate_per_area_g_min_cm2"]
    with open(os.path.join(outdir, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for r in rows:
            w.writerow([r["param"].get(k, r.get(k)) for k in keys])
    print(f"  wrote summary.csv ({len(rows)} cases)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True)
    ap.add_argument("--rate-convention", choices=["physical", "paper", "both"], default="both")
    a = ap.parse_args()
    params, rows = load(a.study)
    if not rows:
        sys.exit(f"no finished cases with results for study '{a.study}'")
    outdir = os.path.join(config.RUNS_DIR, "_figures", a.study)
    param = params[0] if params else None
    write_table(rows, params, outdir)
    fig_liquid_fraction(rows, param, outdir)
    fig_mean_temperature(rows, param, outdir)
    if param:
        fig_melting_bars(rows, param, outdir)
    else:
        fig_rate_twin(rows, outdir)
        fig_probes(rows, outdir)
    if param == "panel.volume_ml":
        for conv in (["physical", "paper"] if a.rate_convention == "both" else [a.rate_convention]):
            fig_normalised_rate(rows, param, outdir, conv)


if __name__ == "__main__":
    main()
