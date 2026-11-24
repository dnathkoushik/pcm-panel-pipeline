"""Phase 3: every case against the pure concrete baseline, as tables and the figures of
Abass & Muthulingam (2024): temperature vs time (Figs. 7-8), heat flux vs time (Fig. 9),
maximum heat-gain reduction (Eq. 8, Fig. 12) and thermal loads (Eq. 9, Fig. 13).

    python3 post/collect.py [--baseline baseline_concrete]        (from roof/)

Reads each run's results/ (written by `python3 -m roofpipe check`; the time series come
from the in-solver function objects, every 300 s). Writes roof/results/phase3/:
    summary_all.csv                  one row per case: loads, peaks, lag, HGR_max, PCM state
    timeseries_<quantity>.csv        analysed day, one column per case:
        T_roof_C, T_ceiling_C        roof and ceiling surface temperatures
        q_roof_in_W_m2               outdoor heat flux, into the slab
        q_ceiling_gain_W_m2          indoor heat flux, into the room
        load_cumulative_Wh_m2        cumulative cooling load (time integral of the gain > 0)
        liquid_fraction              PCM cases
    fig1_temperature_h070.png ... fig6_ceiling_T_all.png
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOF = os.path.dirname(HERE)
sys.path.insert(0, ROOF)
from roofpipe import config, geometry  # noqa: E402

# validated categorical slots 1-3 (dataviz reference palette; all-pairs pass, light mode)
SHAPE_COLOR = {"tet": "#2a78d6", "pyr": "#eb6834", "hem": "#1baf7a"}
SHAPE_NAME = {"tet": "Tetrahedron", "pyr": "Pyramid", "hem": "Hemisphere"}
ORIENT_NAME = {"N": "Normal (base on the ceiling side)", "I": "Inverted (base on the roof side)"}
HEIGHT_RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]   # blue 250..650
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
REF_COLOR = INK


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    return {h: [float(r[i]) for r in rows[1:]] for i, h in enumerate(rows[0])}


def load_run(runs, name):
    res = os.path.join(runs, name, "results")
    try:
        with open(os.path.join(res, "summary.json"), encoding="utf-8") as f:
            summ = json.load(f)
        ts = read_csv(os.path.join(res, "timeseries_lastday.csv"))
        with open(os.path.join(runs, name, "case_meta.json"), encoding="utf-8") as f:
            meta = json.load(f)
    except FileNotFoundError:
        return None
    return {"summary": summ, "ts": ts, "meta": meta}


def cumulative(t, q):
    out, acc = [0.0], 0.0
    for i in range(1, len(t)):
        acc += 0.5 * (max(q[i], 0.0) + max(q[i - 1], 0.0)) * (t[i] - t[i - 1]) / 3600.0
        out.append(acc)
    return out


def write_wide(path, hours, columns):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["hour"] + list(columns))
        for i, h in enumerate(hours):
            w.writerow([f"{h:.4f}"] + [f"{columns[c][i]:.6g}" if i < len(columns[c]) else ""
                                       for c in columns])


def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "axes.titlecolor": INK,
        "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "font.size": 10, "axes.titlesize": 11, "legend.frameon": False,
        "lines.linewidth": 1.8, "lines.solid_capstyle": "round",
    })
    return plt


def clock_axis(ax):
    ax.set_xlim(0, 24)
    ax.set_xticks(range(0, 25, 3), [f"{h:02d}:00" for h in range(0, 25, 3)])
    ax.set_xlabel("Clock time on the analysed day (IST)")


def label_end(ax, x, y, text, color):
    ax.annotate(text, (x[-1], y[-1]), xytext=(4, 0), textcoords="offset points",
                va="center", fontsize=8, color=INK2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", default=None,
                    help="pure-slab run to compare with (default: the one with the PCM study's "
                         "time step: baseline_concrete or baseline_dt<dt>)")
    ap.add_argument("--config", default=None)
    a = ap.parse_args()
    cfg = config.load(a.config)
    runs = config.runs_dir()
    out = os.path.join(ROOF, "results", "phase3")
    os.makedirs(out, exist_ok=True)

    if a.baseline is None:
        dt = float(cfg["pcm_study"]["run"]["deltaT"])
        for name in ("baseline_concrete", f"baseline_dt{dt:g}"):
            r = load_run(runs, name)
            if r is not None and abs(r["meta"]["time"]["deltaT_s"] - dt) < 1e-9:
                a.baseline = name
                break
        else:
            raise SystemExit(f"no pure-slab run with the PCM study's step ({dt:g} s): "
                             f"python3 -m roofpipe baseline --name baseline_dt{dt:g} "
                             f"--set run.deltaT={dt:g}, run it, then check it")
    base = load_run(runs, a.baseline)
    if base is None:
        raise SystemExit(f"no results for the baseline {a.baseline} in {runs}")
    bodies = geometry.study_cases(cfg)
    cases = {}
    for b in bodies:
        r = load_run(runs, b.case_id)
        if r is not None:
            r["body"] = b
            cases[b.case_id] = r
    print(f"baseline + {len(cases)} of {len(bodies)} PCM cases have results")

    # ---------------------------------------------------------------- tables
    bm = base["summary"]["metrics"]
    b_peak = bm["q_ceiling_gain_W_m2"]["max"]
    b_cool = bm["daily_Wh_m2"]["cooling_load"]
    b_heat = bm["daily_Wh_m2"]["heating_load"]
    b_tmax = bm["T_ceiling_C"]["max"]
    rows = [{"case_id": a.baseline, "shape": "none", "orientation": "", "height_mm": "",
             "volume_mL": 0.0, "pcm_mass_g": 0.0, "passed": base["summary"]["passed"],
             "peak_gain_W_m2": b_peak, "peak_gain_at_h": bm["q_ceiling_gain_W_m2"]["at_h"],
             "HGR_max_pct": 0.0, "cooling_load_Wh_m2": b_cool, "heating_load_Wh_m2": b_heat,
             "total_load_Wh_m2": b_cool + b_heat, "cooling_load_reduction_pct": 0.0,
             "ceiling_T_max_C": b_tmax, "ceiling_T_max_reduction_K": 0.0,
             "ceiling_T_min_C": bm["T_ceiling_C"]["min"], "time_lag_h": bm["time_lag_h"],
             "decrement_factor": bm["decrement_factor"], "lf_min": "", "lf_max": "",
             "latent_cycled_Wh_m2": ""}]
    for cid, r in cases.items():
        m, b = r["summary"]["metrics"], r["body"]
        pk = m["q_ceiling_gain_W_m2"]["max"]
        cool = m["daily_Wh_m2"]["cooling_load"]
        heat = m["daily_Wh_m2"]["heating_load"]
        ms = r["meta"].get("meshed", {})
        rows.append({
            "case_id": cid, "shape": b.shape, "orientation": b.orientation,
            "height_mm": round(b.height * 1e3), "volume_mL": b.volume * 1e6,
            "pcm_mass_g": ms.get("pcm_mass_kg", float("nan")) * 1e3,
            "passed": r["summary"]["passed"],
            "peak_gain_W_m2": pk, "peak_gain_at_h": m["q_ceiling_gain_W_m2"]["at_h"],
            "HGR_max_pct": 100.0 * (b_peak - pk) / b_peak,          # paper Eq. 8
            "cooling_load_Wh_m2": cool, "heating_load_Wh_m2": heat, "total_load_Wh_m2": cool + heat,
            "cooling_load_reduction_pct": 100.0 * (b_cool - cool) / b_cool,
            "ceiling_T_max_C": m["T_ceiling_C"]["max"],
            "ceiling_T_max_reduction_K": b_tmax - m["T_ceiling_C"]["max"],
            "ceiling_T_min_C": m["T_ceiling_C"]["min"], "time_lag_h": m["time_lag_h"],
            "decrement_factor": m["decrement_factor"],
            "lf_min": m["pcm"]["liquid_fraction"]["min"], "lf_max": m["pcm"]["liquid_fraction"]["max"],
            "latent_cycled_Wh_m2": m["pcm"]["latent_cycled_Wh_m2"],
        })
    with open(os.path.join(out, "summary_all.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})

    hours = base["ts"]["hour"]
    allruns = {a.baseline: base, **cases}
    for q in ("T_roof_C", "T_ceiling_C", "q_roof_in_W_m2", "q_ceiling_gain_W_m2"):
        write_wide(os.path.join(out, f"timeseries_{q}.csv"), hours,
                   {n: r["ts"][q] for n, r in allruns.items()})
    write_wide(os.path.join(out, "timeseries_load_cumulative_Wh_m2.csv"), hours,
               {n: cumulative(r["ts"]["time_s"], r["ts"]["q_ceiling_gain_W_m2"])
                for n, r in allruns.items()})
    if cases:
        write_wide(os.path.join(out, "timeseries_liquid_fraction.csv"), hours,
                   {n: r["ts"]["liquid_fraction"] for n, r in cases.items()})

    # --------------------------------------------------------------- figures
    plt = style()
    bt = base["ts"]
    ref_h = 70                                          # the paper's PCM mass (618 g)

    def pick(shape, orient, h):
        return cases.get(geometry.case_id(shape, orient, h * 1e-3))

    for fig_name, quantity, ylab, title in (
            ("fig1_temperature_h070.png", "T", "Temperature (°C)",
             "Temperature variation, PCM height 70 mm (0.72 L, 618 g)"),
            ("fig2_heatflux_h070.png", "q", "Heat flux into the room (W/m²)",
             "Indoor heat flux, PCM height 70 mm (0.72 L, 618 g)")):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
        for ax, o in zip(axes, ("N", "I")):
            if quantity == "T":
                ax.plot(bt["hour"], bt["T_out_C"], color=MUTED, lw=1.2, ls="--", label="T_amb (outdoor air)")
                ax.plot(bt["hour"], bt["T_roof_C"], color=INK2, lw=1.2, label="T_ext (roof surface, no PCM)")
                ax.plot(bt["hour"], bt["T_ceiling_C"], color=REF_COLOR, lw=2.0, label="T_int,Ref (no PCM)")
                key = "T_ceiling_C"
            else:
                ax.plot(bt["hour"], bt["q_ceiling_gain_W_m2"], color=REF_COLOR, lw=2.0, label="q_Ref (no PCM)")
                ax.axhline(0.0, color=MUTED, lw=0.8)
                key = "q_ceiling_gain_W_m2"
            for sh in ("tet", "pyr", "hem"):
                r = pick(sh, o, ref_h)
                if r is None:
                    continue
                ax.plot(r["ts"]["hour"], r["ts"][key], color=SHAPE_COLOR[sh], lw=2.0,
                        label=f"{'T_int,Exp' if quantity == 'T' else 'q_Exp'} {SHAPE_NAME[sh].lower()}")
            ax.set_title(ORIENT_NAME[o])
            clock_axis(ax)
        axes[0].set_ylabel(ylab)
        axes[1].legend(fontsize=8, loc="upper left")
        fig.suptitle(title, fontsize=12)
        fig.tight_layout()
        fig.savefig(os.path.join(out, fig_name), dpi=160)
        plt.close(fig)

    if cases:
        # Fig 3: thermal loads (paper Fig. 13). Top: daily cooling load per case against the
        # plain slab; bottom: the load each case saves (the difference is small next to the load)
        hs = [int(h) for h in cfg["pcm_study"]["heights_mm"]]
        fig, axes = plt.subplots(2, 2, figsize=(12, 7.6), sharex=True, sharey="row")
        width = 0.26
        for col, o in enumerate(("N", "I")):
            for k, sh in enumerate(("tet", "pyr", "hem")):
                xs, ys = [], []
                for i, h in enumerate(hs):
                    r = pick(sh, o, h)
                    if r is not None:
                        xs.append(i + (k - 1) * width)
                        ys.append(r["summary"]["metrics"]["daily_Wh_m2"]["cooling_load"])
                for row, vals in ((0, ys), (1, [b_cool - y for y in ys])):
                    axes[row, col].bar(xs, vals, width=width - 0.03, color=SHAPE_COLOR[sh],
                                       label=SHAPE_NAME[sh] if (row, col) == (0, 0) else None,
                                       edgecolor=SURFACE, linewidth=1.0)
                for x, y in zip(xs, ys):
                    axes[1, col].annotate(f"{b_cool - y:.0f}", (x, b_cool - y), xytext=(0, 2),
                                          textcoords="offset points", ha="center", fontsize=7,
                                          color=INK2)
            top_ax = axes[0, col]
            top_ax.axhline(b_cool, color=REF_COLOR, lw=1.5, ls="--")
            top_ax.annotate(f"no PCM: {b_cool:.0f} Wh/m²", (len(hs) - 0.5, b_cool), xytext=(0, 4),
                            textcoords="offset points", ha="right", fontsize=8, color=INK2)
            top_ax.set_ylim(0, b_cool * 1.12)
            top_ax.set_title(ORIENT_NAME[o])
            vols = [geometry.body(cfg, "pyr", o, h * 1e-3).volume * 1e6 for h in hs]
            axes[1, col].set_xticks(range(len(hs)), [f"{h} mm\n{v:.0f} mL" for h, v in zip(hs, vols)])
            axes[1, col].set_xlabel("PCM height and volume (equal for every shape)")
            for row in (0, 1):
                axes[row, col].grid(axis="x", visible=False)
        axes[0, 0].set_ylabel("Daily cooling load (Wh/m²)")
        axes[1, 0].set_ylabel("Load saved against no PCM (Wh/m²)")
        heat_note = ("" if all(r["summary"]["metrics"]["daily_Wh_m2"]["heating_load"] < 0.5
                               for r in [base, *cases.values()])
                     else "  (heating loads in summary_all.csv)")
        fig.suptitle("Thermal load through the roof, room held at "
                     f"{cfg['indoor']['T_set_C']:g} °C{heat_note}", fontsize=12)
        fig.legend(loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=3, fontsize=9)
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        fig.savefig(os.path.join(out, "fig3_thermal_load.png"), dpi=160)
        plt.close(fig)

        # Fig 4: HGR_max vs volume (paper Eq. 8, Fig. 12)
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
        for ax, o in zip(axes, ("N", "I")):
            ends = []
            for sh in ("tet", "pyr", "hem"):
                pts = [(r["volume_mL"], r["HGR_max_pct"]) for r in rows[1:]
                       if r["shape"] == sh and r["orientation"] == o]
                pts.sort()
                if pts:
                    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=SHAPE_COLOR[sh],
                            marker="o", ms=6, label=SHAPE_NAME[sh])
                    ends.append((pts[-1][1], pts[-1][0], SHAPE_NAME[sh]))
            # one end label per cluster: lines that end within 0.5 % share a label
            ends.sort()
            groups = []
            for y, x, name in ends:
                if groups and y - groups[-1][0][-1] < 0.5:
                    groups[-1][0].append(y)
                    groups[-1][2].append(name)
                else:
                    groups.append(([y], x, [name]))
            for ys_, x, names in groups:
                label_end(ax, [x], [sum(ys_) / len(ys_)], " / ".join(names), INK2)
            lo, hi = ax.get_xlim()
            ax.set_xlim(lo, hi + 0.25 * (hi - lo))          # room for the end labels
            ax.axhline(0.0, color=MUTED, lw=0.8)
            ax.set_title(ORIENT_NAME[o])
            ax.set_xlabel("PCM volume per 300 × 300 mm cell (mL)")
        axes[0].set_ylabel("Maximum heat-gain reduction HGR_max (%)")
        axes[1].legend(fontsize=8, loc="upper left")
        fig.suptitle("Peak heat gain through the ceiling, reduction against the plain slab (Eq. 8)",
                     fontsize=12)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "fig4_hgr_vs_volume.png"), dpi=160)
        plt.close(fig)

        # Fig 5: liquid fraction over the day, h = 70 mm
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.0), sharey=True)
        for ax, o in zip(axes, ("N", "I")):
            for sh in ("tet", "pyr", "hem"):
                r = pick(sh, o, ref_h)
                if r is not None:
                    ax.plot(r["ts"]["hour"], r["ts"]["liquid_fraction"], color=SHAPE_COLOR[sh],
                            lw=2.0, label=SHAPE_NAME[sh])
            ax.set_ylim(0, 1)
            ax.set_title(ORIENT_NAME[o])
            clock_axis(ax)
        axes[0].set_ylabel("PCM liquid fraction (–)")
        axes[1].legend(fontsize=8, loc="upper left")
        fig.suptitle("OM37 melting and solidification over the day, PCM height 70 mm", fontsize=12)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "fig5_liquid_fraction_h070.png"), dpi=160)
        plt.close(fig)

        # Fig 6: ceiling temperature, every case (small multiples)
        fig, axes = plt.subplots(2, 3, figsize=(14, 7.2), sharex=True, sharey=True)
        for i, o in enumerate(("N", "I")):
            for j, sh in enumerate(("tet", "pyr", "hem")):
                ax = axes[i][j]
                ax.plot(bt["hour"], bt["T_ceiling_C"], color=REF_COLOR, lw=1.6, ls="--",
                        label="no PCM")
                for k, h in enumerate(hs):
                    r = pick(sh, o, h)
                    if r is not None:
                        ax.plot(r["ts"]["hour"], r["ts"]["T_ceiling_C"], color=HEIGHT_RAMP[k],
                                lw=1.6, label=f"h {h} mm")
                ax.set_title(f"{SHAPE_NAME[sh]}, {'normal' if o == 'N' else 'inverted'}", fontsize=10)
                clock_axis(ax)
                if i == 0:
                    ax.set_xlabel("")
            axes[i][0].set_ylabel("Ceiling surface T (°C)")
        axes[0][2].legend(fontsize=8, loc="upper left")
        fig.suptitle("Ceiling surface temperature, all 30 cases (darker = more PCM)", fontsize=12)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "fig6_ceiling_T_all.png"), dpi=160)
        plt.close(fig)

        # Fig 7: the melting front, normal and inverted pyramid at 70 mm, from the
        # pv_extract contours (skipped until post/pv_extract.py has run on both cases)
        times = [("12p0h", "12:00"), ("15p0h", "15:00"), ("18p0h", "18:00"), ("21p0h", "21:00")]
        pair = [(geometry.case_id("pyr", o, ref_h * 1e-3), name) for o, name in
                (("N", "Normal"), ("I", "Inverted"))]
        pngs = {(cid, tk): os.path.join(runs, cid, "results", "pv", f"contour_T_{tk}.png")
                for cid, _ in pair for tk, _ in times}
        if all(os.path.isfile(p) for p in pngs.values()):
            from matplotlib import image as mpimg
            crop = (slice(92, 668), slice(130, 1270))          # the slab in pv_extract's 1400 x 760
            fig = plt.figure(figsize=(15, 4.6))
            gs = fig.add_gridspec(2, len(times) + 1, width_ratios=[1] * len(times) + [0.32],
                                  wspace=0.04, hspace=0.06)
            for i, (cid, name) in enumerate(pair):
                for j, (tk, label) in enumerate(times):
                    ax = fig.add_subplot(gs[i, j])
                    ax.imshow(mpimg.imread(pngs[(cid, tk)])[crop])
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.grid(False)
                    for sp in ax.spines.values():
                        sp.set_visible(False)
                    if i == 0:
                        ax.set_title(label, fontsize=11)
                    if j == 0:
                        ax.set_ylabel(name, fontsize=11, color=INK)
            cb = fig.add_subplot(gs[:, -1])
            cb.imshow(mpimg.imread(pngs[(pair[0][0], times[0][0])])[488:756, 1318:1400])
            cb.axis("off")
            fig.suptitle(f"Melting front, pyramid at {ref_h} mm: vertical mid-plane, roof on top; "
                         "white: liquid fraction 0.5, cyan: PCM outline", fontsize=12)
            fig.savefig(os.path.join(out, "fig7_melting_front_pyr_h070.png"), dpi=160,
                        bbox_inches="tight")
            plt.close(fig)

    print(f"-> {out}")
    for r in rows:
        if r["shape"] != "none":
            continue
        print(f"  baseline: peak gain {r['peak_gain_W_m2']:.1f} W/m2 at {r['peak_gain_at_h']:.2f} h, "
              f"cooling load {r['cooling_load_Wh_m2']:.0f} Wh/m2")
    best = sorted(rows[1:], key=lambda r: -r["HGR_max_pct"])[:5]
    for r in best:
        print(f"  {r['case_id']}: HGR_max {r['HGR_max_pct']:+.1f} %, cooling load "
              f"{r['cooling_load_Wh_m2']:.0f} Wh/m2 ({r['cooling_load_reduction_pct']:+.1f} %), "
              f"liquid fraction {r['lf_min']:.2f}..{r['lf_max']:.2f}")


if __name__ == "__main__":
    main()
