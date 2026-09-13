#!/usr/bin/env python3
"""Rebuild clinical Figure 6 with collision-free stacked sensitivity forests."""
from __future__ import annotations

import importlib.util
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
CONTENT = SCRIPT_DIR.parents[3]
PORTABLE_FIGURE_ROOT = CONTENT / "paper_upgrade_v1/13_biotransformation_resolved_manuscript_v1/figures_submission_detail_optimized_v1"
SOURCE_DIR = PORTABLE_FIGURE_ROOT / "source/figure6"
SOURCE_SCRIPT = PORTABLE_FIGURE_ROOT / "scripts/make_figure5_redesign_v1.py"
OUT = Path(os.environ.get("FMH_FIGURE_OUT", str(SCRIPT_DIR.parent / "reproduced/Figure6"))).resolve()
QA = OUT / "qa"
STEM = "Figure6_clinical_evidence_FINAL"

spec = importlib.util.spec_from_file_location("figure6_source", SOURCE_SCRIPT)
source = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(source)
source.SOURCE = SOURCE_DIR

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica", "sans-serif"],
    "font.size": 7.2,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def plot_sensitivity_pair(ax_plot, ax_tab, sens, outcome, title, xlim, xlab, show_xlabel=True):
    scenarios = [
        ("primary_digitised_endpoint", "Digitised endpoint"),
        ("ancova_replacement", "ANCOVA replacement"),
        ("complete_exclusion", "TRIAL028 excluded"),
    ]
    y = np.array([2.0, 1.0, 0.0])
    lines = []
    ax_plot.axvline(0, color=source.C["faint"], ls="--", lw=0.85, zorder=1)
    for yi, (scenario_id, _) in zip(y, scenarios):
        row = sens.loc[(sens.scenario == scenario_id) & (sens.outcome == outcome)].iloc[0]
        ax_plot.errorbar(
            row.estimate,
            yi,
            xerr=[[row.estimate - row.ci_lb], [row.ci_ub - row.estimate]],
            fmt="D",
            color=source.C["teal"],
            ecolor=source.C["teal"],
            elinewidth=1.15,
            capsize=2.0,
            ms=5.2,
            markeredgecolor=source.C["ink"],
            markeredgewidth=0.4,
            clip_on=True,
            zorder=3,
        )
        lines.append(f"k={int(row.k)}\n{row.estimate:.2f} [{row.ci_lb:.2f}, {row.ci_ub:.2f}]")

    ax_plot.set_yticks(y, [label for _, label in scenarios], fontproperties=source.FP, fontsize=5.9)
    ax_plot.tick_params(axis="y", pad=2.5)
    ax_plot.set_xlim(*xlim)
    ax_plot.set_ylim(-0.45, 2.45)
    ax_plot.set_xlabel(xlab if show_xlabel else "", fontproperties=source.FP, fontsize=6.2, labelpad=1)
    ax_plot.set_title(title, loc="left", fontsize=6.8, fontweight="bold", color=source.C["teal"] if outcome == "HbA1c" else source.C["purple"], pad=2, fontproperties=source.FP)
    source.style(ax_plot, "x")

    ax_tab.set_xlim(0, 1)
    ax_tab.set_ylim(ax_plot.get_ylim())
    ax_tab.axis("off")
    ax_tab.text(0.02, 2.35, "k; MD [95% CI]", fontsize=5.4, color=source.C["mute"], fontproperties=source.FP, va="bottom")
    for yi, line in zip(y, lines):
        ax_tab.text(0.02, yi, line, va="center", ha="left", fontsize=5.3, color=source.C["ink"], fontproperties=source.FP, linespacing=1.06)


def export(fig, axes, panel_ids, source_files):
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    source.require_matplotlib_panel_alignment(
        fig,
        axes=axes,
        panel_ids=panel_ids,
        row_groups=[{"id": "primary-forests", "panels": ["c", "d"]}],
        column_groups=[],
        exemptions=[{
            "panels": ["c", "d"],
            "checks": ["panel-width", "column"],
            "reason": "Primary forests use plot and numerical-table splits with outcome-specific x-axis widths.",
        }],
        json_out=QA / f"{STEM}_alignment.json",
        overlay_svg=QA / f"{STEM}_alignment_overlay.svg",
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        require_panel_labels=True,
        strict=True,
    )
    outputs = {
        "pdf": OUT / f"{STEM}.pdf",
        "svg": OUT / f"{STEM}.svg",
        "tiff": OUT / f"{STEM}.tiff",
        "preview": OUT / f"{STEM}_preview.png",
    }
    # Keep the physical canvas fixed at the intended submission width.  A tight
    # bounding box would enlarge the page to include external labels and then
    # silently shrink type when the journal places the figure at 183 mm.
    fig.savefig(outputs["pdf"])
    fig.savefig(outputs["svg"])
    fig.savefig(outputs["tiff"], dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(outputs["preview"], dpi=320)
    plt.close(fig)
    (OUT / f"{STEM}_manifest.json").write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "layout_change": "Only the RoB 2 completion-status note was updated after signed SY/XL consensus; values and layout were not changed.",
        "source_files": [str(path) for path in source_files],
        "outputs": {key: str(path) for key, path in outputs.items()},
    }, indent=2), encoding="utf-8")


def main():
    screen = pd.read_csv(SOURCE_DIR / "Figure5_screening_flow.tsv", sep="\t")
    maturity = pd.read_csv(SOURCE_DIR / "Figure5_maturity47.tsv", sep="\t")
    meta = source.load_meta()
    sensitivity = pd.read_csv(SOURCE_DIR / "Figure5_TRIAL028_sensitivity.tsv", sep="\t")
    rob = pd.read_csv(SOURCE_DIR / "Figure5_TRIAL028_RoB2.tsv", sep="\t")

    fig = plt.figure(figsize=(7.205, 9.8))
    outer = fig.add_gridspec(3, 1, height_ratios=[1.05, 1.40, 1.37], left=0.142, right=0.975, top=0.955, bottom=0.055, hspace=0.15)

    top = outer[0].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.55)
    ax_a = fig.add_subplot(top[0, 0])
    b_col = top[0, 1].subgridspec(2, 1, height_ratios=[4.6, 1.25], hspace=0.25)
    ax_b = fig.add_subplot(b_col[0, 0])
    ax_b_note = fig.add_subplot(b_col[1, 0])

    mid = outer[1].subgridspec(2, 1, height_ratios=[12.0, 1.15], hspace=0.26)
    mid_plots = mid[0].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.46)
    c_split = mid_plots[0, 0].subgridspec(1, 2, width_ratios=[2.4, 1.45], wspace=0.10)
    d_split = mid_plots[0, 1].subgridspec(1, 2, width_ratios=[2.4, 1.45], wspace=0.10)
    ax_c = fig.add_subplot(c_split[0, 0]); ax_c_tab = fig.add_subplot(c_split[0, 1])
    ax_d = fig.add_subplot(d_split[0, 0]); ax_d_tab = fig.add_subplot(d_split[0, 1])
    ax_cd_leg = fig.add_subplot(mid[1])

    bottom = outer[2].subgridspec(1, 2, width_ratios=[1.48, 1.0], wspace=0.40)
    e_col = bottom[0, 0].subgridspec(4, 1, height_ratios=[0.42, 2.15, 2.15, 0.62], hspace=0.52)
    ax_e_title = fig.add_subplot(e_col[0]); ax_e_title.axis("off")
    h_split = e_col[1].subgridspec(1, 2, width_ratios=[2.15, 1.35], wspace=0.11)
    f_split = e_col[2].subgridspec(1, 2, width_ratios=[2.15, 1.35], wspace=0.11)
    ax_e_h = fig.add_subplot(h_split[0, 0]); ax_e_ht = fig.add_subplot(h_split[0, 1])
    ax_e_f = fig.add_subplot(f_split[0, 0]); ax_e_ft = fig.add_subplot(f_split[0, 1])
    ax_e_note = fig.add_subplot(e_col[3]); ax_e_note.axis("off")

    f_col = bottom[0, 1].subgridspec(3, 1, height_ratios=[0.42, 4.63, 0.62], hspace=0.22)
    ax_f_pad = fig.add_subplot(f_col[0]); ax_f_pad.axis("off")
    ax_f = fig.add_subplot(f_col[1])
    ax_f_note = fig.add_subplot(f_col[2])

    for ax, label in [(ax_a, "a"), (ax_c, "c"), (ax_d, "d"), (ax_f, "f")]:
        source.panel_label(ax, label)
    ax_b.text(-0.10, 1.02, "b", transform=ax_b.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left", color=source.C["ink"], fontproperties=source.FP)
    ax_e_title.text(0.0, 0.35, "e   Saffron sensitivity to TRIAL028 handling", fontsize=8.4, fontweight="bold", color=source.C["ink"], fontproperties=source.FP, va="center")

    source.plot_a(ax_a, screen)
    # Start exclusion connectors just outside the main boxes.  At the fixed
    # physical width, the longest count line otherwise touches a connector.
    for connector in ax_a.lines:
        xdata = np.asarray(connector.get_xdata(), dtype=float).copy()
        if xdata.size == 2:
            xdata[0] += 0.025
            connector.set_xdata(xdata)
    source.plot_b(ax_b, maturity, ax_b_note)
    # Keep the maturity caveat inside its own note row on the fixed-width page.
    ax_b_note.clear(); ax_b_note.axis("off")
    ax_b_note.text(
        0.0, 0.58,
        "Formula signal is not attributable single-material efficacy.\nZero ≠ ineffective.",
        transform=ax_b_note.transAxes, fontsize=5.4, color=source.C["mute"],
        fontproperties=source.FP, va="center", ha="left", linespacing=1.15,
    )
    source.plot_forest(ax_c, ax_c_tab, meta, "HbA1c", "Primary HbA1c estimates", "Mean difference (pp)", (-3.1, 2.4))
    source.plot_forest(ax_d, ax_d_tab, meta, "FPG", "Primary FPG estimates", "Mean difference (mg/dL)", (-72, 72))
    source.plot_cd_legend(ax_cd_leg)
    # The upper mini-forest omits a redundant x-axis title so it cannot collide
    # with the lower FPG section at the final 183-mm submission width.
    plot_sensitivity_pair(ax_e_h, ax_e_ht, sensitivity, "HbA1c", "HbA1c", (-0.70, 0.55), "HbA1c MD (pp)", show_xlabel=False)
    plot_sensitivity_pair(ax_e_f, ax_e_ft, sensitivity, "FPG", "FPG", (-22, 16), "FPG MD (mg/dL)")
    ax_e_note.text(0.0, 0.62, "Row 1 = primary saffron pools in c/d. All panel-e CIs cross zero.", fontsize=5.5, color=source.C["mute"], fontproperties=source.FP, va="center")
    source.plot_f(ax_f, rob, ax_f_note)
    # Replace the source legend with a compact, fixed-width key; the original
    # three-column legend extended beyond the canvas after physical sizing.
    ax_f_note.clear(); ax_f_note.axis("off")
    rob_key = [
        (0.02, source.C["low"], "+  Low risk"),
        (0.34, source.C["some"], "!  Some concerns"),
        (0.77, source.C["high"], "−  High risk"),
    ]
    for x, colour, label in rob_key:
        ax_f_note.add_patch(mpl.patches.Rectangle((x, 0.56), 0.035, 0.18, transform=ax_f_note.transAxes, facecolor=colour, edgecolor="none"))
        ax_f_note.text(x + 0.047, 0.65, label, transform=ax_f_note.transAxes, fontsize=5.3, color=source.C["ink"], fontproperties=source.FP, va="center", ha="left")
    ax_f_note.text(0.0, 0.18, "Bias risk only (not efficacy).\nIndependent duplicate assessment completed.", transform=ax_f_note.transAxes, fontsize=5.1, color=source.C["mute"], fontproperties=source.FP, va="center", ha="left", linespacing=1.12)

    fig.canvas.draw()
    pos_c = ax_c.get_position(); pos_b = ax_b.get_position(); pos_bn = ax_b_note.get_position()
    ax_a.set_position([pos_c.x0, pos_bn.y0, pos_b.x0 - 0.055 - pos_c.x0, (pos_b.y0 + pos_b.height) - pos_bn.y0])
    fig.suptitle("Clinical evidence maturity and TRIAL028 sensitivity", x=pos_c.x0, y=0.988, ha="left", fontsize=10.8, fontweight="bold", color=source.C["ink"], fontproperties=source.FP)

    export(fig, [ax_c, ax_d], ["c", "d"], [
        SOURCE_DIR / "Figure5_screening_flow.tsv",
        SOURCE_DIR / "Figure5_maturity47.tsv",
        SOURCE_DIR / "Figure5_primary_meta.tsv",
        SOURCE_DIR / "Figure5_TRIAL028_sensitivity.tsv",
        SOURCE_DIR / "Figure5_TRIAL028_RoB2.tsv",
    ])


if __name__ == "__main__":
    main()
