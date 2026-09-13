#!/usr/bin/env python3
"""Nature redesign Figure 3: 2x2 corner layout with separated plot/table columns."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches
from matplotlib.ticker import FixedFormatter, FixedLocator, NullFormatter

sys.path.insert(0, str(Path(__file__).resolve().parent))
from redesign_common import (  # noqa: E402
    C,
    FP_SANS,
    OUT,
    ROOT,
    SOURCE,
    draw_stats_cards,
    ensure_dirs,
    panel_label,
    save_figure,
    style,
)

STEM = "Figure3_three_cohort_butyrate_verified_v2"
PORTABLE_SOURCE = SOURCE / "figure4"
X_LIM = (0.55, 2.45)


def freeze() -> dict[str, Path]:
    ensure_dirs()
    paths = {
        "effects": PORTABLE_SOURCE / "Figure3_cohort_effects.tsv",
        "pooled": PORTABLE_SOURCE / "Figure3_pooled_effects.tsv",
        "met": PORTABLE_SOURCE / "Figure3_metformin_sensitivity.tsv",
        "context": PORTABLE_SOURCE / "Figure3_cohort_context.tsv",
        "expanded": PORTABLE_SOURCE / "Figure3_expanded_SCFA_audit.tsv",
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing portable Figure 4 inputs: " + ", ".join(missing))
    return paths


def primary_rows(effects: pd.DataFrame, pathway: str) -> pd.DataFrame:
    return effects[
        (effects["pathway_id"] == pathway)
        & (effects["variant"] == "primary_demographic_adjusted")
        & (effects["metric"] == "positive_abundance")
    ].copy()


def pooled_row(pooled: pd.DataFrame, pathway: str, family: str = "primary_demographic_adjusted") -> pd.Series:
    # Portable rebuild uses explicit method names.  Normalise the historical
    # no-metformin family label and derive the two-pathway BH family in place.
    if "variant" in pooled.columns:
        variant = "known_no_metformin" if family == "no_metformin_sensitivity" else family
        rows = pooled[(pooled.pathway_id == pathway) & (pooled.variant == variant) & (pooled.test == "modified_KH")].copy()
        if rows.empty:
            raise KeyError(f"No modified_KH row for {variant}/{pathway}")
        fam = pooled[(pooled.variant == variant) & (pooled.test == "modified_KH")].copy()
        p = fam.p_value.astype(float)
        order = np.argsort(p.to_numpy())
        ranked = p.to_numpy()[order]
        q = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
        q_full = np.empty_like(q); q_full[order] = np.clip(q, 0, 1)
        q_map = dict(zip(fam.pathway_id, q_full, strict=True))
        row = rows.iloc[0].copy()
        row["ratio_CI_lower"] = row["ratio_ci95_low"]
        row["ratio_CI_upper"] = row["ratio_ci95_high"]
        row["q_BH_within_family_test"] = q_map[pathway]
        return row
    return pooled[(pooled.pathway_id == pathway) & (pooled.meta_family == family) & (pooled.test == "adhoc")].iloc[0]


def plot_a(ax_plot, ax_n, context: pd.DataFrame) -> None:
    df = context.iloc[::-1]
    y = np.arange(len(df), dtype=float)
    # Stretch three bars across the same vertical span as forest panels
    y = np.array([2.6, 1.5, 0.4])
    ax_plot.barh(y, df["control"].to_numpy(), color=C["blue"], height=0.72, label="Control", zorder=2)
    ax_plot.barh(y, df["t2d"].to_numpy(), left=df["control"].to_numpy(), color=C["warm"], height=0.72, label="T2D", zorder=2)
    ax_plot.set_yticks(y, df["cohort"].tolist(), fontproperties=FP_SANS, fontsize=7.2)
    ax_plot.set_ylim(-0.4, 3.5)
    ax_plot.set_xlim(0, 820)
    ax_plot.set_xlabel("Eligible participants", fontproperties=FP_SANS, fontsize=7.2, labelpad=1)
    ax_plot.set_title("Cohort composition", loc="left", fontweight="bold", color=C["ink"], pad=6, fontsize=8.8, fontproperties=FP_SANS)
    style(ax_plot, "x")

    ax_n.set_ylim(ax_plot.get_ylim())
    draw_stats_cards(
        ax_n,
        y,
        [f"n={int(row.eligible)}" for _, row in df.iterrows()],
        header="n",
        header_y=3.35,
        face_cycle=[C["band"], "#EEF3F6", "#F3F0EA"],
        half_h=0.30,
        fontsize=6.4,
    )


def _style_log_ratio_axis(ax) -> None:
    ax.set_xscale("log")
    ax.set_xlim(*X_LIM)
    ax.set_xticks([0.6, 0.8, 1.0, 1.5, 2.0])
    ax.xaxis.set_major_formatter(FixedFormatter(["0.6", "0.8", "1", "1.5", "2"]))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.axvline(1.0, color=C["faint"], ls="--", lw=0.85, zorder=1)


def plot_forest_split(ax_plot, ax_tab, effects, pooled, pathway, title, color) -> None:
    """Forest points on ax_plot; statistics only on ax_tab (no overlap)."""
    df = primary_rows(effects, pathway)
    names = ["Qin 2012", "Karlsson 2013", "MetaCardis 2020"]
    keys = ["Qin", "Karlsson", "MetaCardis"]
    y = np.array([3.0, 2.0, 1.0])

    ax_plot.axhspan(-0.35, 0.45, color=C["band"], zorder=0, lw=0)
    lines = []
    for yi, key, name in zip(y, keys, names):
        row = df[df["cohort"] == key].iloc[0]
        est = float(row.transformed_estimate)
        lo = float(row.transformed_ci95_low)
        hi = float(row.transformed_ci95_high)
        ax_plot.errorbar(
            est,
            yi,
            xerr=[[est - lo], [hi - est]],
            fmt="o",
            ms=5.4,
            color=color,
            ecolor=color,
            elinewidth=1.25,
            capsize=2.4,
            zorder=3,
        )
        lines.append(f"n={int(row.n)}\n{est:.3f} [{lo:.3f}–{hi:.3f}]")

    s = pooled_row(pooled, pathway)
    est = float(s.abundance_ratio)
    lo = float(s.ratio_CI_lower)
    hi = float(s.ratio_CI_upper)
    ax_plot.errorbar(
        est,
        0.05,
        xerr=[[est - lo], [hi - est]],
        fmt="D",
        ms=6.2,  # ~1.15x study points
        color=C["ink"],
        ecolor=C["ink"],
        elinewidth=1.35,
        capsize=2.4,
        zorder=4,
    )
    lines.append(
        f"REML mKH; max(1,q)\n{est:.3f} [{lo:.3f}-{hi:.3f}]\n"
        f"P={float(s.p_value):.3f}; q={float(s.q_BH_within_family_test):.3f}\n"
        f"I2={float(s.I2):.1f}%"
    )

    _style_log_ratio_axis(ax_plot)
    ax_plot.set_ylim(-0.75, 3.65)
    ax_plot.set_yticks(list(y) + [0.05], names + ["Pooled"], fontproperties=FP_SANS, fontsize=7.0)
    ax_plot.set_xlabel("Abundance ratio (T2D/control)", fontproperties=FP_SANS, fontsize=7.2, labelpad=1)
    ax_plot.set_title(title, loc="left", fontweight="bold", color=C["ink"], pad=6, fontsize=8.8, fontproperties=FP_SANS)
    style(ax_plot, "x")

    y_pos = list(y) + [0.05]
    ax_tab.set_ylim(-0.75, 3.65)
    draw_stats_cards(
        ax_tab,
        y_pos,
        lines,
        header="Estimate [95% CI]",
        header_y=3.55,
        face_cycle=["#F3F0EA", "#EEF3F6", "#F3F0EA", "#E8EEF2"],
        half_h=0.40,
        fontsize=5.5,
        linespacing=1.05,
    )


def plot_d_split(ax_plot, ax_tab, ax_foot, effects, pooled, met, expanded) -> None:
    primary = primary_rows(effects, "PWY-5676")
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
    pool_no = pooled_row(pooled, "PWY-5676", "no_metformin_sensitivity")

    rows = [
        ("MetaCardis primary", float(mcard.transformed_estimate), float(mcard.transformed_ci95_low), float(mcard.transformed_ci95_high), float(mcard.p_value), C["warm"]),
        ("No metformin", float(m_no.transformed_estimate), float(m_no.transformed_CI_lower), float(m_no.transformed_CI_upper), float(m_no.p_value), C["teal"]),
        ("Metformin-adj.", float(m_adj.transformed_estimate), float(m_adj.transformed_CI_lower), float(m_adj.transformed_CI_upper), float(m_adj.p_value), C["teal"]),
        ("3-cohort non-met.", float(pool_no.abundance_ratio), float(pool_no.ratio_CI_lower), float(pool_no.ratio_CI_upper), float(pool_no.p_value), C["blue"]),
    ]
    y = np.arange(len(rows))[::-1].astype(float)

    for yi, (lab, est, lo, hi, p, col) in zip(y, rows):
        eb = ax_plot.errorbar(est, yi, xerr=[[est - lo], [hi - est]], fmt="o", ms=5.2, color=col, ecolor=col, capsize=2.0, elinewidth=1.15, zorder=3, clip_on=True)
        for artist in (eb[0], *eb[1], *eb[2]):
            artist.set_clip_on(True)

    _style_log_ratio_axis(ax_plot)
    ax_plot.set_ylim(-0.55, 3.55)
    ax_plot.set_yticks(y, [r[0] for r in rows], fontproperties=FP_SANS, fontsize=6.4)
    ax_plot.tick_params(axis="y", pad=5.0)
    ax_plot.set_xlabel("PWY-5676 abundance ratio", fontproperties=FP_SANS, fontsize=7.2, labelpad=1)
    ax_plot.set_title("Metformin-related sensitivity", loc="left", fontweight="bold", color=C["ink"], pad=6, fontsize=8.8, fontproperties=FP_SANS)
    style(ax_plot, "x")

    ax_tab.set_ylim(-0.55, 3.55)
    draw_stats_cards(
        ax_tab,
        y,
        [f"{est:.3f} [{lo:.3f}–{hi:.3f}]\nP={p:.3g}" for lab, est, lo, hi, p, col in rows],
        header="Estimate; P",
        header_y=3.45,
        face_cycle=["#F3F0EA", "#EEF3F6", "#F3F0EA", "#E8EEF2"],
        half_h=0.38,
        fontsize=5.8,
        linespacing=1.08,
    )

    # Footer strip under BOTH plot and table (0/25 badge + sensitivity note)
    ax_foot.set_xlim(0, 1)
    ax_foot.set_ylim(0, 1)
    ax_foot.axis("off")
    n_t = int(expanded.pathways_tested.iloc[0])
    n_s = int(expanded.bh_fdr_significant.iloc[0])
    ax_foot.add_patch(
        patches.Rectangle(
            (0.02, 0.08),
            0.96,
            0.70,
            facecolor="#F7F1E6",
            edgecolor="none",
            clip_on=False,
        )
    )
    ax_foot.text(
        0.50,
        0.55,
        f"{n_s}/{n_t} SCFA pathways BH-FDR < 0.05",
        ha="center",
        va="center",
        fontsize=5.9,
        fontweight="bold",
        color=C["ink"],
        fontproperties=FP_SANS,
    )
    ax_foot.text(
        0.50,
        0.28,
        "Subset sensitivity; not a causal metformin effect.",
        ha="center",
        va="center",
        fontsize=5.6,
        color=C["mute"],
        fontproperties=FP_SANS,
    )


def main() -> None:
    files = freeze()
    context = pd.read_csv(files["context"], sep="\t")
    effects = pd.read_csv(files["effects"], sep="\t")
    pooled = pd.read_csv(files["pooled"], sep="\t")
    met = pd.read_csv(files["met"], sep="\t")
    expanded = pd.read_csv(files["expanded"], sep="\t")

    # Compact mid gutter: kill excess white under 3a/3b; keep equal corner heights
    fig = plt.figure(figsize=(7.205, 7.5))
    outer = fig.add_gridspec(
        2,
        2,
        left=0.12,
        right=0.985,
        top=0.915,
        bottom=0.055,
        wspace=0.34,
        hspace=0.26,
        width_ratios=[1.0, 1.40],
        height_ratios=[1.0, 1.0],
    )

    def corner(spec, width_ratios=(1.25, 1.05)):
        nest = spec.subgridspec(2, 1, height_ratios=[5.4, 0.72], hspace=0.18)
        top = nest[0].subgridspec(1, 2, width_ratios=width_ratios, wspace=0.14)
        return fig.add_subplot(top[0, 0]), fig.add_subplot(top[0, 1]), fig.add_subplot(nest[1])

    ax_a, ax_a_n, ax_a_foot = corner(outer[0, 0], width_ratios=(3.4, 0.75))
    ax_b_plot, ax_b_tab, ax_b_foot = corner(outer[0, 1])
    ax_c_plot, ax_c_tab, ax_c_foot = corner(outer[1, 0])
    ax_d_plot, ax_d_tab, ax_d_foot = corner(outer[1, 1])

    panel_label(ax_a, "a")
    panel_label(ax_b_plot, "b")
    panel_label(ax_c_plot, "c")
    panel_label(ax_d_plot, "d")

    plot_a(ax_a, ax_a_n, context)
    ax_a_foot.set_xlim(0, 1)
    ax_a_foot.set_ylim(0, 1)
    ax_a_foot.axis("off")
    ax_a_foot.plot([0.18, 0.26], [0.62, 0.62], color=C["blue"], lw=5.5, solid_capstyle="butt")
    ax_a_foot.text(0.28, 0.62, "Control", va="center", fontsize=6.2, fontproperties=FP_SANS)
    ax_a_foot.plot([0.52, 0.60], [0.62, 0.62], color=C["warm"], lw=5.5, solid_capstyle="butt")
    ax_a_foot.text(0.62, 0.62, "T2D", va="center", fontsize=6.2, fontproperties=FP_SANS)
    ax_a_foot.text(0.50, 0.28, "Total eligible 1,243.", ha="center", va="center", fontsize=5.5, color=C["mute"], fontproperties=FP_SANS)

    plot_forest_split(ax_b_plot, ax_b_tab, effects, pooled, "PWY-5676", "PWY-5676", C["warm"])
    ax_b_foot.set_xlim(0, 1)
    ax_b_foot.set_ylim(0, 1)
    ax_b_foot.axis("off")
    ax_b_foot.text(0.50, 0.45, "Shared log abundance-ratio scale (b-d).", ha="center", va="center", fontsize=5.8, color=C["mute"], fontproperties=FP_SANS)

    plot_forest_split(ax_c_plot, ax_c_tab, effects, pooled, "CENTFERM-PWY", "CENTFERM-PWY", C["teal"])
    ax_c_foot.set_xlim(0, 1)
    ax_c_foot.set_ylim(0, 1)
    ax_c_foot.axis("off")
    ax_c_foot.text(0.50, 0.45, "Primary demographic-adjusted model.", ha="center", va="center", fontsize=5.8, color=C["mute"], fontproperties=FP_SANS)

    plot_d_split(ax_d_plot, ax_d_tab, ax_d_foot, effects, pooled, met, expanded)

    fig.suptitle(
        "Cross-cohort heterogeneity in butyrate-related pathway associations",
        x=0.09,
        y=0.965,
        ha="left",
        fontsize=10.5,
        fontweight="bold",
        color=C["ink"],
        fontproperties=FP_SANS,
    )

    pdf = save_figure(
        fig,
        STEM,
        axes=[ax_a, ax_b_plot, ax_c_plot, ax_d_plot],
        panel_ids=list("abcd"),
        row_groups=[{"id": "top", "panels": ["a", "b"]}, {"id": "bottom", "panels": ["c", "d"]}],
        column_groups=[],
        sources=list(files.values()),
        caption="Figure 3. Butyrate-related pathway associations were pathway- and cohort-dependent.",
        boundary="Functional potential among pathway-positive samples; not faecal butyrate concentration or flux.",
        notes=[
            "Polish v3: colored Estimate/n cards; reduced white under 3a; compact footers; equal corner heights.",
        ],
        exemptions=[
            {
                "panels": ["a", "b", "c", "d"],
                "checks": ["panel-width", "column"],
                "reason": "Corners include n/table/footer strips so plot widths differ by design.",
            }
        ],
    )
    print(pdf)


if __name__ == "__main__":
    main()
