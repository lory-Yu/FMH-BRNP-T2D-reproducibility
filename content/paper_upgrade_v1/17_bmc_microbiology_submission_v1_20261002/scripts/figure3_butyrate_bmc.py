#!/usr/bin/env python3
"""Build BMC Microbiology Figure 3 from frozen three-cohort source tables.

The displayed medication analyses are sensitivity analyses for ascertainment and
subset composition.  They are deliberately not labelled as causal metformin
effects or formal interactions.
"""
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

PACKAGE = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("FMH_FIGURE_OUT", str(PACKAGE / "figures"))).resolve()
SOURCE = PACKAGE / "source_data" / "figure3"
ORIGINAL_SCRIPTS = Path(__file__).resolve().parent
os.environ["FMH_FIGURE_OUT"] = str(OUT)
sys.path.insert(0, str(ORIGINAL_SCRIPTS))

import make_figure3_evidence_completion as base  # noqa: E402
import redesign_common as common  # noqa: E402

C = common.C
FP_SANS = common.FP_SANS
STEM = os.environ.get("FMH_FIGURE_STEM", "Figure3_butyrate_pathways_BMC")

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
        ("No recorded\nmetformin", float(m_no.transformed_estimate), float(m_no.transformed_CI_lower), float(m_no.transformed_CI_upper), float(m_no.p_value), C["teal"]),
        ("Metformin-adjusted", float(m_adj.transformed_estimate), float(m_adj.transformed_CI_lower), float(m_adj.transformed_CI_upper), float(m_adj.p_value), C["teal"]),
        ("3-cohort\nnon-users", float(pool_no.abundance_ratio), float(pool_no.ratio_CI_lower), float(pool_no.ratio_CI_upper), float(pool_no.p_value), C["blue"]),
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
            fontsize=5.8,
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
    ax_plot.set_title("Medication-ascertainment\nsensitivity", loc="left", fontweight="bold", color=C["ink"], pad=4, fontsize=8.3, fontproperties=FP_SANS)
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
        fontsize=5.3,
        linespacing=1.08,
    )

    ax_foot.set_xlim(0, 1)
    ax_foot.set_ylim(0, 1)
    ax_foot.axis("off")
    n_t = int(expanded.pathways_tested.iloc[0])
    n_s = int(expanded.bh_fdr_significant.iloc[0])
    ax_foot.add_patch(patches.Rectangle((0.02, 0.04), 0.96, 0.72, facecolor="#F7F1E6", edgecolor="none", clip_on=False))
    ax_foot.text(0.50, 0.55, f"{n_s}/{n_t} SCFA pathways BH-FDR < 0.05", ha="center", va="center", fontsize=5.7, fontweight="bold", color=C["ink"], fontproperties=FP_SANS)
    ax_foot.text(0.50, 0.22, "Sensitivity to medication ascertainment; not a causal drug effect.", ha="center", va="center", fontsize=5.3, color=C["mute"], fontproperties=FP_SANS)


def plot_forest_contained(ax_label, ax_plot, ax_tab, effects, pooled, pathway, title, color) -> None:
    """Forest panel with cohort labels confined to a dedicated internal column."""
    df = base.primary_rows(effects, pathway)
    names = ["Qin 2012", "Karlsson 2013", "MetaCardis"]
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
    lines.append(f"REML + modified KH\n{est:.3f} [{lo:.3f}–{hi:.3f}]\nP={float(summary.p_value):.3f}; q={float(summary.q_BH_within_family_test):.3f}\nI²={float(summary.I2):.1f}%")

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
        "png_600dpi": OUT / f"{STEM}.png",
        "preview": OUT / f"{STEM}_preview.png",
    }
    fig.savefig(outputs["pdf"], facecolor="white")
    fig.savefig(outputs["svg"], facecolor="white")
    fig.savefig(outputs["tiff"], dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(outputs["png_600dpi"], dpi=600, facecolor="white")
    fig.savefig(outputs["preview"], dpi=320, facecolor="white")
    plt.close(fig)
    legend = (
        "Figure 3. Cross-cohort heterogeneity in butyrate-related pathway associations. "
        "(a) Eligible T2D and control participants in the three source cohorts; model-specific sample sizes were lower when covariates or positive-abundance values were unavailable. "
        "(b,c) Cohort-specific demographic-adjusted abundance ratios and 95% confidence intervals for PWY-5676 and CENTFERM-PWY among detected complete cases. Circles denote cohort estimates and diamonds denote random-effects summaries fitted by REML with modified Knapp-Hartung inference; q values are Benjamini-Hochberg adjusted within the two-pathway family. The lower note strips report the complementary detection models: PWY-5676 was estimable in MetaCardis (OR 2.009, 95% CI 1.211–3.334) but unstable in Qin and invariant in Karlsson; CENTFERM-PWY was estimable in Qin (OR 0.490, 95% CI 0.163–1.474) and MetaCardis (OR 1.203, 95% CI 0.737–1.964) but unstable in Karlsson. "
        "(d) PWY-5676 sensitivity analyses based on the MetaCardis primary estimate, exclusion of participants with recorded metformin use, covariate adjustment for metformin, and a three-cohort analysis restricted to recorded non-users. The medication panels assess ascertainment and subset sensitivity, not a causal drug effect or treatment-by-pathway interaction. None of 25 expanded short-chain-fatty-acid pathways passed BH-FDR < 0.05."
    )
    (OUT / f"{STEM}_legend.txt").write_text(legend + "\n", encoding="utf-8")
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "change": "BMC Microbiology layout; detection-model results surfaced in panels b and c; MetaCardis display label standardised; panel-d labels contained; medication wording limited to ascertainment sensitivity; all numeric inputs unchanged.",
        "sources": [str(p) for p in sources],
        "outputs": {k: str(v) for k, v in outputs.items()},
        "final_size_mm": [170.0, 173.99],
        "legend": legend,
    }
    (OUT / f"{STEM}_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> None:
    files = {
        "context": SOURCE / "Figure3_cohort_context.tsv",
        "effects": SOURCE / "Figure3_cohort_effects.tsv",
        "pooled": SOURCE / "three_cohort_REML_KH_rebuilt.tsv",
        "met": SOURCE / "Figure3_metformin_sensitivity.tsv",
        "expanded": SOURCE / "Figure3_expanded_SCFA_audit.tsv",
    }
    missing = [str(p) for p in files.values() if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing Figure 3 source tables: " + ", ".join(missing))
    context = pd.read_csv(files["context"], sep="\t")
    # Display-only standardisation. The frozen source table remains unchanged.
    context["cohort"] = context["cohort"].replace({"MetaCardis 2020": "MetaCardis"})
    effects = pd.read_csv(files["effects"], sep="\t")
    pooled = pd.read_csv(files["pooled"], sep="\t")
    met = pd.read_csv(files["met"], sep="\t")
    expanded = pd.read_csv(files["expanded"], sep="\t")

    fig = plt.figure(figsize=(6.6929, 6.85))
    outer = fig.add_gridspec(2, 2, left=0.15, right=0.985, top=0.895, bottom=0.055, wspace=0.30, hspace=0.26, width_ratios=[1.0, 1.40], height_ratios=[1.0, 1.0])

    def corner(spec, width_ratios=(1.25, 1.05)):
        nest = spec.subgridspec(2, 1, height_ratios=[5.15, 1.0], hspace=0.25)
        top = nest[0].subgridspec(1, 2, width_ratios=width_ratios, wspace=0.14)
        return fig.add_subplot(top[0, 0]), fig.add_subplot(top[0, 1]), fig.add_subplot(nest[1])

    ax_a, ax_a_n, ax_a_foot = corner(outer[0, 0], width_ratios=(3.4, 0.75))
    b_nest = outer[0, 1].subgridspec(2, 1, height_ratios=[5.15, 1.0], hspace=0.25)
    b_top = b_nest[0].subgridspec(1, 3, width_ratios=(1.05, 1.40, 1.08), wspace=0.10)
    ax_b_label = fig.add_subplot(b_top[0, 0])
    ax_b_plot = fig.add_subplot(b_top[0, 1])
    ax_b_tab = fig.add_subplot(b_top[0, 2])
    ax_b_foot = fig.add_subplot(b_nest[1])
    ax_c_plot, ax_c_tab, ax_c_foot = corner(outer[1, 0])

    d_nest = outer[1, 1].subgridspec(2, 1, height_ratios=[5.15, 1.0], hspace=0.25)
    d_top = d_nest[0].subgridspec(1, 3, width_ratios=(1.34, 1.24, 1.16), wspace=0.10)
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
    ax_a_foot.plot([0.18, 0.26], [0.44, 0.44], color=C["blue"], lw=5.5, solid_capstyle="butt")
    ax_a_foot.text(0.28, 0.44, "Control", va="center", fontsize=6.2, fontproperties=FP_SANS)
    ax_a_foot.plot([0.52, 0.60], [0.44, 0.44], color=C["warm"], lw=5.5, solid_capstyle="butt")
    ax_a_foot.text(0.62, 0.44, "T2D", va="center", fontsize=6.2, fontproperties=FP_SANS)
    ax_a_foot.text(0.50, 0.09, "Total eligible 1,243.", ha="center", va="center", fontsize=5.5, color=C["mute"], fontproperties=FP_SANS)

    plot_forest_contained(ax_b_label, ax_b_plot, ax_b_tab, effects, pooled, "PWY-5676", "PWY-5676", C["warm"])
    ax_b_foot.set_xlim(0, 1); ax_b_foot.set_ylim(0, 1); ax_b_foot.axis("off")
    ax_b_foot.add_patch(patches.Rectangle((0.01, 0.03), 0.98, 0.88, facecolor="#F3F6F7", edgecolor="none"))
    ax_b_foot.text(
        0.04,
        0.73,
        "Detection OR (95% CI): MetaCardis 2.01 (1.21–3.33);",
        ha="left",
        va="center",
        fontsize=5.05,
        color=C["ink"],
        fontproperties=FP_SANS,
    )
    ax_b_foot.text(
        0.04,
        0.45,
        "Qin unstable; Karlsson invariant.",
        ha="left",
        va="center",
        fontsize=5.05,
        color=C["mute"],
        fontproperties=FP_SANS,
    )
    ax_b_foot.text(
        0.04,
        0.22,
        "Forest: positive abundance, detected complete cases.",
        ha="left",
        va="center",
        fontsize=5.05,
        color=C["mute"],
        fontproperties=FP_SANS,
    )

    base.plot_forest_split(ax_c_plot, ax_c_tab, effects, pooled, "CENTFERM-PWY", "CENTFERM-PWY", C["teal"])
    for txt in ax_c_tab.texts:
        label = txt.get_text().replace("REML mKH; max(1,q)", "REML + modified KH").replace("I2=", "I²=")
        txt.set_text(label)
    ax_c_foot.set_xlim(0, 1); ax_c_foot.set_ylim(0, 1); ax_c_foot.axis("off")
    ax_c_foot.add_patch(patches.Rectangle((0.01, 0.03), 0.98, 0.88, facecolor="#F3F6F7", edgecolor="none"))
    ax_c_foot.text(
        0.04,
        0.73,
        "Detection OR (95% CI): Qin 0.49 (0.16–1.47);",
        ha="left",
        va="center",
        fontsize=5.05,
        color=C["ink"],
        fontproperties=FP_SANS,
    )
    ax_c_foot.text(
        0.04,
        0.45,
        "MetaCardis 1.20 (0.74–1.96); Karlsson unstable.",
        ha="left",
        va="center",
        fontsize=5.05,
        color=C["mute"],
        fontproperties=FP_SANS,
    )
    ax_c_foot.text(
        0.04,
        0.22,
        "Forest: positive abundance, detected complete cases.",
        ha="left",
        va="center",
        fontsize=5.05,
        color=C["mute"],
        fontproperties=FP_SANS,
    )

    plot_d_contained(ax_d_label, ax_d_plot, ax_d_tab, ax_d_foot, effects, pooled, met, expanded)
    fig.suptitle("Cross-cohort heterogeneity in butyrate-related pathway associations", x=0.09, y=0.982, ha="left", fontsize=10.2, fontweight="bold", color=C["ink"], fontproperties=FP_SANS)

    export(fig, [ax_a, ax_b_label, ax_c_plot, ax_d_label], list("abcd"), list(files.values()))


if __name__ == "__main__":
    main()
