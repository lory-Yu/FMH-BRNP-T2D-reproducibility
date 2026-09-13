#!/usr/bin/env python3
"""Rebuild Figure 3 with panel-d row labels contained inside panel d."""
from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_ROOT = SCRIPT_DIR.parent
OUT = Path(os.environ.get(
    "FMH_FIGURE_OUT",
    str(FIGURE_ROOT / "reproduced/Figure4"),
)).resolve()
ORIGINAL_SCRIPTS = SCRIPT_DIR
os.environ["FMH_FIGURE_OUT"] = str(OUT)
sys.path.insert(0, str(ORIGINAL_SCRIPTS))

import make_figure3_evidence_completion as base  # noqa: E402
import redesign_common as common  # noqa: E402

C = common.C
FP_SANS = common.FP_SANS
STEM = os.environ.get("FMH_FIGURE_STEM", "Figure3_three_cohort_butyrate_cd_boundary_corrected_v1")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def plot_d_contained(ax_label, ax_plot, ax_tab, ax_foot, effects, pooled, met, expanded) -> None:
    """Draw panel d with a dedicated, clipped label column inside its grid cell."""
    primary = base.primary_rows(effects, "PWY-5676")
    mcard = primary[primary.cohort == "MetaCardis"].iloc[0]
    m_no = met[
        (met.pathway_id == "PWY-5676")
        & (met.model == "disease_no_metformin_full_DEU_FRA")
        & (met.metric == "positive_abundance")
    ].iloc[0]
    m_adj = met[
        (met.pathway_id == "PWY-5676")
        & (met.model == "disease_metformin_adjusted_full_DEU_FRA")
        & (met.metric == "positive_abundance")
    ].iloc[0]
    pool_no = base.pooled_row(pooled, "PWY-5676", "no_metformin_sensitivity")

    rows = [
        ("MetaCardis primary", float(mcard.transformed_estimate), float(mcard.transformed_ci95_low), float(mcard.transformed_ci95_high), float(mcard.p_value), C["warm"]),
        ("No metformin", float(m_no.transformed_estimate), float(m_no.transformed_CI_lower), float(m_no.transformed_CI_upper), float(m_no.p_value), C["teal"]),
        ("Metformin-adj.", float(m_adj.transformed_estimate), float(m_adj.transformed_CI_lower), float(m_adj.transformed_CI_upper), float(m_adj.p_value), C["teal"]),
        ("3-cohort non-met.", float(pool_no.abundance_ratio), float(pool_no.ratio_CI_lower), float(pool_no.ratio_CI_upper), float(pool_no.p_value), C["blue"]),
    ]
    y = np.arange(len(rows))[::-1].astype(float)

    ax_label.set_xlim(0, 1)
    ax_label.set_ylim(-0.55, 3.55)
    ax_label.axis("off")
    for yi, row in zip(y, rows):
        ax_label.text(
            0.98,
            yi,
            row[0],
            ha="right",
            va="center",
            fontsize=6.2,
            color=C["mute"],
            fontproperties=FP_SANS,
            clip_on=True,
        )

    for yi, (_, est, lo, hi, _, col) in zip(y, rows):
        eb = ax_plot.errorbar(
            est,
            yi,
            xerr=[[est - lo], [hi - est]],
            fmt="o",
            ms=5.2,
            color=col,
            ecolor=col,
            capsize=2.0,
            elinewidth=1.15,
            zorder=3,
            clip_on=True,
        )
        for artist in (eb[0], *eb[1], *eb[2]):
            artist.set_clip_on(True)

    base._style_log_ratio_axis(ax_plot)
    ax_plot.set_ylim(-0.55, 3.55)
    ax_plot.set_yticks([])
    ax_plot.set_xlabel("PWY-5676 abundance ratio", fontproperties=FP_SANS, fontsize=7.2, labelpad=1)
    ax_plot.set_title("Metformin-related sensitivity", loc="left", fontweight="bold", color=C["ink"], pad=6, fontsize=8.8, fontproperties=FP_SANS)
    common.style(ax_plot, "x")

    ax_tab.set_ylim(-0.55, 3.55)
    common.draw_stats_cards(
        ax_tab,
        y,
        [f"{est:.3f} [{lo:.3f}–{hi:.3f}]\nP={p:.3g}" for _, est, lo, hi, p, _ in rows],
        header="Estimate; P",
        header_y=3.45,
        face_cycle=["#F3F0EA", "#EEF3F6", "#F3F0EA", "#E8EEF2"],
        half_h=0.38,
        fontsize=5.8,
        linespacing=1.08,
    )

    ax_foot.set_xlim(0, 1)
    ax_foot.set_ylim(0, 1)
    ax_foot.axis("off")
    n_t = int(expanded.pathways_tested.iloc[0])
    n_s = int(expanded.bh_fdr_significant.iloc[0])
    ax_foot.add_patch(patches.Rectangle((0.02, 0.08), 0.96, 0.70, facecolor="#F7F1E6", edgecolor="none", clip_on=False))
    ax_foot.text(0.50, 0.55, f"{n_s}/{n_t} SCFA pathways BH-FDR < 0.05", ha="center", va="center", fontsize=5.9, fontweight="bold", color=C["ink"], fontproperties=FP_SANS)
    ax_foot.text(0.50, 0.28, "Subset sensitivity; not a causal metformin effect.", ha="center", va="center", fontsize=5.6, color=C["mute"], fontproperties=FP_SANS)


def plot_forest_contained(ax_label, ax_plot, ax_tab, effects, pooled, pathway, title, color) -> None:
    """Forest panel with cohort labels confined to a dedicated internal column."""
    df = base.primary_rows(effects, pathway)
    names = ["Qin 2012", "Karlsson 2013", "MetaCardis 2020"]
    keys = ["Qin", "Karlsson", "MetaCardis"]
    y = np.array([3.0, 2.0, 1.0])
    lines = []

    ax_plot.axhspan(-0.35, 0.45, color=C["band"], zorder=0, lw=0)
    for yi, key in zip(y, keys):
        row = df[df["cohort"] == key].iloc[0]
        est = float(row.transformed_estimate)
        lo = float(row.transformed_ci95_low)
        hi = float(row.transformed_ci95_high)
        ax_plot.errorbar(est, yi, xerr=[[est - lo], [hi - est]], fmt="o", ms=5.4, color=color, ecolor=color, elinewidth=1.25, capsize=2.4, zorder=3)
        lines.append(f"n={int(row.n)}\n{est:.3f} [{lo:.3f}–{hi:.3f}]")

    summary = base.pooled_row(pooled, pathway)
    est = float(summary.abundance_ratio)
    lo = float(summary.ratio_CI_lower)
    hi = float(summary.ratio_CI_upper)
    ax_plot.errorbar(est, 0.05, xerr=[[est - lo], [hi - est]], fmt="D", ms=6.2, color=C["ink"], ecolor=C["ink"], elinewidth=1.35, capsize=2.4, zorder=4)
    lines.append(f"REML mKH; max(1,q)\n{est:.3f} [{lo:.3f}–{hi:.3f}]\nP={float(summary.p_value):.3f}; q={float(summary.q_BH_within_family_test):.3f}\nI2={float(summary.I2):.1f}%")

    ax_label.set_xlim(0, 1)
    ax_label.set_ylim(-0.75, 3.65)
    ax_label.axis("off")
    for yi, label in zip([3.0, 2.0, 1.0, 0.05], names + ["Pooled"]):
        ax_label.text(0.98, yi, label, ha="right", va="center", fontsize=6.4, color=C["mute"], fontproperties=FP_SANS, clip_on=True)

    base._style_log_ratio_axis(ax_plot)
    ax_plot.set_ylim(-0.75, 3.65)
    ax_plot.set_yticks([])
    ax_plot.set_xlabel("Abundance ratio (T2D/control)", fontproperties=FP_SANS, fontsize=7.2, labelpad=1)
    ax_plot.set_title(title, loc="left", fontweight="bold", color=C["ink"], pad=6, fontsize=8.8, fontproperties=FP_SANS)
    common.style(ax_plot, "x")

    ax_tab.set_ylim(-0.75, 3.65)
    common.draw_stats_cards(ax_tab, [3.0, 2.0, 1.0, 0.05], lines, header="Estimate [95% CI]", header_y=3.55, face_cycle=["#F3F0EA", "#EEF3F6", "#F3F0EA", "#E8EEF2"], half_h=0.40, fontsize=5.5, linespacing=1.05)


def export(fig, axes, panel_ids, sources) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    qa = OUT / "qa"
    qa.mkdir(exist_ok=True)
    fig.canvas.draw()
    common.require_matplotlib_panel_alignment(
        fig,
        axes=axes,
        panel_ids=panel_ids,
        row_groups=[{"id": "top", "panels": ["a", "b"]}, {"id": "bottom", "panels": ["c", "d"]}],
        column_groups=[],
        exemptions=[{
            "panels": ["a", "b", "c", "d"],
            "checks": ["panel-width", "column"],
            "reason": "Each corner contains internal label, plot, table and footer columns of intentionally unequal width.",
        }],
        json_out=qa / f"{STEM}_alignment.json",
        overlay_svg=qa / f"{STEM}_alignment_overlay.svg",
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
    fig.savefig(outputs["pdf"], bbox_inches="tight")
    fig.savefig(outputs["svg"], bbox_inches="tight")
    fig.savefig(outputs["tiff"], dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(outputs["preview"], dpi=320, bbox_inches="tight")
    plt.close(fig)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "change": "Panel-d row labels moved from outboard y tick labels to a dedicated internal label column; all numeric inputs unchanged.",
        "sources": [str(p) for p in sources],
        "outputs": {k: str(v) for k, v in outputs.items()},
    }
    (OUT / f"{STEM}_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> None:
    files = base.freeze()
    context = pd.read_csv(files["context"], sep="\t")
    effects = pd.read_csv(files["effects"], sep="\t")
    pooled = pd.read_csv(files["pooled"], sep="\t")
    met = pd.read_csv(files["met"], sep="\t")
    expanded = pd.read_csv(files["expanded"], sep="\t")

    fig = plt.figure(figsize=(7.205, 7.5))
    outer = fig.add_gridspec(2, 2, left=0.12, right=0.985, top=0.915, bottom=0.055, wspace=0.30, hspace=0.26, width_ratios=[1.0, 1.40], height_ratios=[1.0, 1.0])

    def corner(spec, width_ratios=(1.25, 1.05)):
        nest = spec.subgridspec(2, 1, height_ratios=[5.4, 0.72], hspace=0.18)
        top = nest[0].subgridspec(1, 2, width_ratios=width_ratios, wspace=0.14)
        return fig.add_subplot(top[0, 0]), fig.add_subplot(top[0, 1]), fig.add_subplot(nest[1])

    ax_a, ax_a_n, ax_a_foot = corner(outer[0, 0], width_ratios=(3.4, 0.75))
    b_nest = outer[0, 1].subgridspec(2, 1, height_ratios=[5.4, 0.72], hspace=0.18)
    b_top = b_nest[0].subgridspec(1, 3, width_ratios=(1.05, 1.40, 1.08), wspace=0.10)
    ax_b_label = fig.add_subplot(b_top[0, 0])
    ax_b_plot = fig.add_subplot(b_top[0, 1])
    ax_b_tab = fig.add_subplot(b_top[0, 2])
    ax_b_foot = fig.add_subplot(b_nest[1])
    ax_c_plot, ax_c_tab, ax_c_foot = corner(outer[1, 0])

    d_nest = outer[1, 1].subgridspec(2, 1, height_ratios=[5.4, 0.72], hspace=0.18)
    d_top = d_nest[0].subgridspec(1, 3, width_ratios=(1.12, 1.38, 1.05), wspace=0.10)
    ax_d_label = fig.add_subplot(d_top[0, 0])
    ax_d_plot = fig.add_subplot(d_top[0, 1])
    ax_d_tab = fig.add_subplot(d_top[0, 2])
    ax_d_foot = fig.add_subplot(d_nest[1])

    common.panel_label(ax_a, "a")
    common.panel_label(ax_b_label, "b")
    common.panel_label(ax_c_plot, "c")
    common.panel_label(ax_d_label, "d")

    base.plot_a(ax_a, ax_a_n, context)
    ax_a_foot.set_xlim(0, 1); ax_a_foot.set_ylim(0, 1); ax_a_foot.axis("off")
    ax_a_foot.plot([0.18, 0.26], [0.62, 0.62], color=C["blue"], lw=5.5, solid_capstyle="butt")
    ax_a_foot.text(0.28, 0.62, "Control", va="center", fontsize=6.2, fontproperties=FP_SANS)
    ax_a_foot.plot([0.52, 0.60], [0.62, 0.62], color=C["warm"], lw=5.5, solid_capstyle="butt")
    ax_a_foot.text(0.62, 0.62, "T2D", va="center", fontsize=6.2, fontproperties=FP_SANS)
    ax_a_foot.text(0.50, 0.28, "Total eligible 1,243.", ha="center", va="center", fontsize=5.5, color=C["mute"], fontproperties=FP_SANS)

    plot_forest_contained(ax_b_label, ax_b_plot, ax_b_tab, effects, pooled, "PWY-5676", "PWY-5676", C["warm"])
    ax_b_foot.set_xlim(0, 1); ax_b_foot.set_ylim(0, 1); ax_b_foot.axis("off")
    ax_b_foot.text(0.50, 0.45, "Shared log abundance-ratio scale (b–d).", ha="center", va="center", fontsize=5.8, color=C["mute"], fontproperties=FP_SANS)

    base.plot_forest_split(ax_c_plot, ax_c_tab, effects, pooled, "CENTFERM-PWY", "CENTFERM-PWY", C["teal"])
    ax_c_foot.set_xlim(0, 1); ax_c_foot.set_ylim(0, 1); ax_c_foot.axis("off")
    ax_c_foot.text(0.50, 0.45, "Primary demographic-adjusted model.", ha="center", va="center", fontsize=5.8, color=C["mute"], fontproperties=FP_SANS)

    plot_d_contained(ax_d_label, ax_d_plot, ax_d_tab, ax_d_foot, effects, pooled, met, expanded)
    fig.suptitle("Cross-cohort heterogeneity in butyrate-related pathway associations", x=0.09, y=0.965, ha="left", fontsize=10.5, fontweight="bold", color=C["ink"], fontproperties=FP_SANS)

    export(fig, [ax_a, ax_b_label, ax_c_plot, ax_d_label], list("abcd"), list(files.values()))


if __name__ == "__main__":
    main()
