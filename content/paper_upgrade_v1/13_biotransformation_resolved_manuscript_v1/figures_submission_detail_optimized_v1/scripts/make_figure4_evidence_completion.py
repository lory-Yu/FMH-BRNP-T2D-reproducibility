#!/usr/bin/env python3
"""Nature redesign Figure 4: Wei 2025 public-data reanalysis."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches

sys.path.insert(0, str(Path(__file__).resolve().parent))
from redesign_common import (  # noqa: E402
    C,
    FP_SANS,
    OUT,
    ROOT,
    SOURCE,
    draw_stats_cards,
    ensure_dirs,
    note_card,
    panel_label,
    save_figure,
    style,
)

STEM = "Figure4_Wei2025_same_source_verified_v2"
PORTABLE_SOURCE = SOURCE / "figure5"


def freeze() -> dict[str, Path]:
    ensure_dirs()
    mapping = {
        "alpha": "Figure3_alpha_by_donor.tsv",
        "stats": "Figure3_alpha_stats.tsv",
        "pcoa": "Figure3_pcoa.tsv",
        "cand": "Figure3_candidate_genera.tsv",
        "global": "Figure3_global_genera.tsv",
    }
    out = {key: PORTABLE_SOURCE / name.replace("Figure3_", "Figure4_") for key, name in mapping.items()}
    out["qc"] = PORTABLE_SOURCE / "Figure4_qc_summary.tsv"
    out["beta"] = PORTABLE_SOURCE / "Figure4_beta_stats.tsv"
    missing = [str(path) for path in out.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing portable Figure 5 inputs: " + ", ".join(missing))
    return out


def plot_a(ax, qc: pd.DataFrame) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Sample linkage and sequencing QC", loc="left", fontweight="bold", color=C["ink"], pad=3, fontsize=8.8, fontproperties=FP_SANS)
    colors = [C["blue"], C["teal"], C["amber"], C["warm"]]
    n = len(qc)
    start = 0.02
    gap = 0.070
    w = (0.96 - (n - 1) * gap) / n
    y0, h = 0.22, 0.58  # fill more of the top-row height
    for i, (_, row) in enumerate(qc.iterrows()):
        x = start + i * (w + gap)
        ax.add_patch(
            patches.FancyBboxPatch(
                (x, y0),
                w,
                h,
                boxstyle="round,pad=0.008,rounding_size=0.01",
                facecolor=colors[i] + "20",
                edgecolor=colors[i],
                lw=0.85,
                zorder=2,
            )
        )
        ax.text(x + w / 2, y0 + h * 0.62, row["value"], ha="center", va="center", fontsize=8.6, fontweight="bold", color=C["ink"], fontproperties=FP_SANS, zorder=3)
        ax.text(x + w / 2, y0 + h * 0.28, row["stage"], ha="center", va="center", fontsize=6.4, color=C["mute"], fontproperties=FP_SANS, zorder=3)
        if i < n - 1:
            pad = 0.012
            ax.annotate(
                "",
                xy=(x + w + gap - pad, y0 + h / 2),
                xytext=(x + w + pad, y0 + h / 2),
                arrowprops={"arrowstyle": "-|>", "lw": 0.9, "color": C["faint"], "shrinkA": 0, "shrinkB": 0},
                zorder=1,
            )
    ax.text(0.5, 0.08, "Public-data reprocessing of source-defined converter labels — not an enzyme assay.", ha="center", fontsize=6.1, color=C["mute"], fontproperties=FP_SANS)


def beeswarm(ax, alpha, stats):
    sub = alpha[alpha.rb1_group.isin(["High", "Low"])].copy()
    groups = ["Low", "High"]
    colors = [C["low"], C["high"]]
    for xi, (g, col) in enumerate(zip(groups, colors)):
        vals = sub.loc[sub.rb1_group == g, "Observed_ASV"].to_numpy(float)
        assert len(vals) == 10
        order = np.argsort(vals)
        xpos = np.zeros(len(vals))
        for k, idx in enumerate(order):
            xpos[idx] = ((k % 2) * 2 - 1) * (0.03 + 0.015 * ((k // 2) % 3))
        ax.scatter(xi + xpos, vals, s=36, color=col, edgecolors="white", lw=0.55, zorder=3, alpha=0.95)
        med = float(np.median(vals))
        q1, q3 = np.percentile(vals, [25, 75])
        ax.plot([xi - 0.18, xi + 0.18], [med, med], color=C["ink"], lw=1.6, zorder=4)
        ax.plot([xi, xi], [q1, q3], color=C["ink"], lw=1.05, zorder=4, solid_capstyle="round")
        ax.plot([xi - 0.08, xi + 0.08], [q1, q1], color=C["ink"], lw=1.0, zorder=4)
        ax.plot([xi - 0.08, xi + 0.08], [q3, q3], color=C["ink"], lw=1.0, zorder=4)
    row = stats[stats.metric == "Observed_ASV"].iloc[0]
    ax.set_xticks([0, 1], ["Low (n=10)", "High (n=10)"], fontproperties=FP_SANS, fontsize=7.0)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel("Observed ASV richness", fontproperties=FP_SANS, fontsize=7.6)
    ax.set_title("Observed richness by converter group", loc="left", fontweight="bold", color=C["ink"], pad=6, fontsize=8.8, fontproperties=FP_SANS)
    style(ax, "y")
    return row


def plot_pcoa(ax, pcoa):
    for g, col, mk in [("High", C["high"], "o"), ("Low", C["low"], "s")]:
        df = pcoa[pcoa.rb1_group == g]
        assert len(df) == 10
        ax.scatter(df.PCoA1, df.PCoA2, s=42, c=col, marker=mk, edgecolors="white", lw=0.55, label=f"{g} (n=10)", zorder=3)
        ax.scatter([df.PCoA1.mean()], [df.PCoA2.mean()], s=70, c=col, marker="X", edgecolors=C["ink"], lw=0.5, zorder=4)
    ax.axhline(0, color=C["grid"], lw=0.6)
    ax.axvline(0, color=C["grid"], lw=0.6)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("PCoA1 (Bray-Curtis)", fontproperties=FP_SANS, fontsize=7.2, labelpad=1)
    ax.set_ylabel("PCoA2", fontproperties=FP_SANS, fontsize=7.6)
    ax.set_title("Bray-Curtis ordination (High/Low n=10)", loc="left", fontweight="bold", color=C["ink"], pad=6, fontsize=8.6, fontproperties=FP_SANS)
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=C["high"], markeredgecolor="white", label="High"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=C["low"], markeredgecolor="white", label="Low"),
        Line2D([0], [0], marker="X", color="none", markerfacecolor=C["faint"], markeredgecolor=C["ink"], label="Centroid"),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper right", fontsize=5.8, handletextpad=0.25, borderaxespad=0.2)
    style(ax, None)



def plot_cand(ax_plot, ax_q, cand):
    df = cand[cand.detected.astype(str).str.lower().eq("true")].copy()
    short = {
        "Anaerostipes": "Anaerostipes",
        "norank_f_[Eubacterium] coprostanoligenes group": "E. coprostanoligenes",
        "Coprococcus": "Coprococcus",
        "Barnesiella": "Barnesiella",
        "norank_f_Oscillospiraceae": "Oscillospiraceae",
        "Family XIII AD3011 group": "Family XIII AD3011",
    }
    df["label"] = df["observed_genus"].map(short).fillna(df["observed_genus"])
    df = df.sort_values("clr_difference_high_vs_low")
    y = np.arange(len(df))
    est = df.clr_difference_high_vs_low.to_numpy(float)
    se = df.clr_se.to_numpy(float)
    lo, hi = est - 1.96 * se, est + 1.96 * se
    ax_plot.errorbar(est, y, xerr=[est - lo, hi - est], fmt="o", color=C["teal"], ecolor=C["teal"], ms=5.0, elinewidth=1.1, capsize=2.0)
    ax_plot.axvline(0, color=C["faint"], ls="--", lw=0.8)
    ax_plot.set_yticks(y, df["label"].tolist(), fontproperties=FP_SANS, fontsize=7.0, fontstyle="italic")
    ax_plot.set_xlabel("Adjusted CLR coefficient (High - Low)", fontproperties=FP_SANS, fontsize=7.2, labelpad=1)
    ax_plot.set_title("Candidate genus associations", loc="left", fontweight="bold", color=C["ink"], pad=6, fontsize=8.8, fontproperties=FP_SANS)
    style(ax_plot, "x")

    ax_q.set_ylim(ax_plot.get_ylim())
    draw_stats_cards(
        ax_q,
        y,
        [f"q={row.candidate_set_adjusted_q:.3f}" for _, row in df.iterrows()],
        header="q",
        header_y=float(y.max()) + 0.65,
        face_cycle=["#EEF3F6", "#F3F0EA"],
        half_h=0.36,
        fontsize=6.0,
    )


def plot_qrank(ax, global_df, cand):
    df = global_df.sort_values("adjusted_q").reset_index(drop=True).copy()
    df["rank"] = np.arange(1, len(df) + 1)
    df["neglog10q"] = -np.log10(df["adjusted_q"].clip(lower=1e-12))
    is_cand = df["genus"].isin(["Anaerostipes", "Coprococcus", "Barnesiella"]) | df["genus"].astype(str).str.contains(
        "Oscillospiraceae|coprostanoligenes|AD3011", case=False, na=False
    )
    ax.scatter(df.loc[~is_cand, "rank"], df.loc[~is_cand, "neglog10q"], s=12, c="#BDBDBD", edgecolors="none", zorder=2)
    ax.scatter(df.loc[is_cand, "rank"], df.loc[is_cand, "neglog10q"], s=28, c=C["teal"], edgecolors="white", lw=0.4, zorder=3)
    thr = -np.log10(0.05)
    ax.axhline(thr, color=C["red"], ls="--", lw=0.85)
    xmax = float(df["rank"].max())
    ax.text(
        0.98,
        thr - 0.18,
        "q=0.05",
        transform=ax.get_yaxis_transform(),
        color=C["red"],
        fontsize=6.2,
        va="top",
        ha="right",
        fontproperties=FP_SANS,
        clip_on=False,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.8, "alpha": 0.96},
    )
    n_pass = int((df.adjusted_q < 0.05).sum())
    ax.set_xlabel("Genera ranked by adjusted q", fontproperties=FP_SANS, fontsize=7.2, labelpad=1)
    ax.set_ylabel("-log10(BH q)", fontproperties=FP_SANS, fontsize=7.6)
    ax.set_title("Multiplicity across tested genera", loc="left", fontweight="bold", color=C["ink"], pad=6, fontsize=8.8, fontproperties=FP_SANS)
    from matplotlib.lines import Line2D

    ax.set_xlim(-2, xmax + 4)
    style(ax, "y")
    return n_pass, len(df)


def main() -> None:
    files = freeze()
    alpha = pd.read_csv(files["alpha"], sep="\t")
    stats = pd.read_csv(files["stats"], sep="\t")
    pcoa = pd.read_csv(files["pcoa"], sep="\t")
    cand = pd.read_csv(files["cand"], sep="\t")
    global_df = pd.read_csv(files["global"], sep="\t")
    qc = pd.read_csv(files["qc"], sep="\t")

    fig = plt.figure(figsize=(7.205, 8.2))
    gs = fig.add_gridspec(
        3,
        2,
        height_ratios=[0.58, 1.35, 1.30],
        left=0.12,
        right=0.92,
        top=0.925,
        bottom=0.07,
        wspace=0.42,
        hspace=0.36,
    )
    ax_a = fig.add_subplot(gs[0, :])

    b_box = gs[1, 0].subgridspec(2, 1, height_ratios=[5.0, 1.05], hspace=0.35)
    ax_b = fig.add_subplot(b_box[0])
    ax_b_foot = fig.add_subplot(b_box[1])
    c_box = gs[1, 1].subgridspec(2, 1, height_ratios=[5.0, 1.05], hspace=0.35)
    ax_c = fig.add_subplot(c_box[0])
    ax_c_foot = fig.add_subplot(c_box[1])

    d_box = gs[2, 0].subgridspec(2, 1, height_ratios=[5.0, 1.05], hspace=0.35)
    d_inner = d_box[0].subgridspec(1, 2, width_ratios=[3.2, 0.7], wspace=0.08)
    ax_d = fig.add_subplot(d_inner[0, 0])
    ax_d_q = fig.add_subplot(d_inner[0, 1])
    ax_d_foot = fig.add_subplot(d_box[1])
    e_box = gs[2, 1].subgridspec(2, 1, height_ratios=[5.0, 1.05], hspace=0.35)
    ax_e = fig.add_subplot(e_box[0])
    ax_e_foot = fig.add_subplot(e_box[1])

    axes = [ax_a, ax_b, ax_c, ax_d, ax_e]
    for ax, lab in zip(axes, "abcde"):
        panel_label(ax, lab)

    plot_a(ax_a, qc)
    row = beeswarm(ax_b, alpha, stats)
    plot_pcoa(ax_c, pcoa)
    plot_cand(ax_d, ax_d_q, cand)
    n_pass, n_tot = plot_qrank(ax_e, global_df, cand)

    for ax in (ax_b_foot, ax_c_foot, ax_d_foot, ax_e_foot):
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax_d_foot.add_patch(__import__("matplotlib.patches", fromlist=["FancyBboxPatch"]).FancyBboxPatch((0.02, 0.18), 0.96, 0.64, boxstyle="round,pad=0.02,rounding_size=0.04", facecolor="#F3F6F7", edgecolor="none"))
    ax_d_foot.text(0.50, 0.45, "Candidate-set BH; 0/6 passed.\nTaxonomic association only.", ha="center", va="center", fontsize=5.3, color=C["mute"], fontproperties=FP_SANS, linespacing=1.15)
    ax_b_foot.add_patch(__import__("matplotlib.patches", fromlist=["FancyBboxPatch"]).FancyBboxPatch((0.02, 0.15), 0.96, 0.70, boxstyle="round,pad=0.02,rounding_size=0.04", facecolor="#F3F6F7", edgecolor="none"))
    ax_b_foot.text(0.50, 0.45, f"Median {row.high_median:.0f} vs {row.low_median:.0f}\nWilcoxon P={row.wilcoxon_p:.3f} · Adj. P={row.adjusted_p:.3f}; q={row.adjusted_q:.4f}", ha="center", va="center", fontsize=5.4, fontproperties=FP_SANS, linespacing=1.15)
    ax_c_foot.add_patch(__import__("matplotlib.patches", fromlist=["FancyBboxPatch"]).FancyBboxPatch((0.02, 0.12), 0.96, 0.76, boxstyle="round,pad=0.02,rounding_size=0.04", facecolor="#F3F6F7", edgecolor="none"))
    ax_c_foot.text(0.50, 0.45, "PERMANOVA R2=0.066; P=0.131\nANOSIM R=0.058; P=0.160 · PERMDISP P=0.092", ha="center", va="center", fontsize=5.3, fontproperties=FP_SANS, linespacing=1.15)
    ax_e_foot.add_patch(__import__("matplotlib.patches", fromlist=["FancyBboxPatch"]).FancyBboxPatch((0.02, 0.18), 0.96, 0.64, boxstyle="round,pad=0.02,rounding_size=0.04", facecolor="#F3F6F7", edgecolor="none"))
    ax_e_foot.text(0.50, 0.45, f"Teal = source candidates\nn={n_tot}; q<0.05: {n_pass}/{n_tot}", ha="center", va="center", fontsize=5.3, fontproperties=FP_SANS, linespacing=1.15)

    fig.suptitle(
        "Reanalysis of microbiome associations in source-defined Rb1 converter groups",
        x=0.11,
        y=0.968,
        ha="left",
        fontsize=10.5,
        fontweight="bold",
        color=C["ink"],
        fontproperties=FP_SANS,
    )

    pdf = save_figure(
        fig,
        STEM,
        axes=axes,
        panel_ids=list("abcde"),
        row_groups=[{"id": "mid", "panels": ["b", "c"]}, {"id": "bot", "panels": ["d", "e"]}],
        column_groups=[],
        sources=list(files.values()),
        caption="Figure 4. Public-data reanalysis of the Wei 2025 Rb1 converter groups.",
        boundary="16S composition associated with source-defined labels; not beta-glucosidase activity or Rb1 to Compound K kinetics.",
        notes=[
            "Polish v3: compact 4a; equal mid/bot heights; outside legends; colored q cards; stats note cards.",
            "Exactly 10 High and 10 Low points in panel b.",
        ],
        exemptions=[
            {
                "panels": ["a", "c", "d"],
                "checks": ["column", "panel-width", "vertical-gutter"],
                "reason": "Panel a full-width; c equal-aspect PCoA; d includes q text strip.",
            }
        ],
    )
    print(pdf)


if __name__ == "__main__":
    main()
