#!/usr/bin/env python3
"""Build Figure 2 from frozen metabolism-aware v2 and AGORA2 outputs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_ROOT = SCRIPT_DIR.parent
CONTENT = FIGURE_ROOT.parents[2]
V2 = CONTENT / "paper_upgrade_v1/metabolism_aware_v2_20260907_013732"
AGORA_SUMMARY = FIGURE_ROOT / "source/figure2/12_analysis_summary.csv"
sys.path.insert(0, str(SCRIPT_DIR))
from audit_panel_alignment import require_matplotlib_panel_alignment  # noqa: E402


COLORS = {
    "ink": "#263238",
    "muted": "#66747D",
    "grid": "#D9E0E3",
    "direct": "#5C7C8A",
    "transform": "#2A8C82",
    "flux": "#B47A45",
    "clinical": "#9B6A88",
    "light_teal": "#DCEDEA",
    "light_blue": "#E5EDF0",
    "light_orange": "#F3E8DC",
    "grey": "#BCC6CB",
}

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 6.4,
        "axes.titlesize": 7.2,
        "axes.labelsize": 6.5,
        "xtick.labelsize": 5.8,
        "ytick.labelsize": 5.8,
        "legend.fontsize": 5.6,
        "axes.linewidth": 0.65,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    }
)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.16,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="top",
        ha="left",
        color=COLORS["ink"],
    )


def box(ax: plt.Axes, x: float, y: float, w: float, h: float, text: str, fc: str) -> None:
    wrapped = {
        "104 materials · 30,800 herb–CID pairs": "104 materials\n30,800 herb–CID pairs",
        "95 terminal tuples → 1,694 material–route paths": "95 terminal tuples\n1,694 material–route paths",
        "Primary subset: 572 paths · 35 transformations · 85 materials": "Primary subset: 572 paths\n35 transformations · 85 materials",
    }
    rect = mpl.patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.025",
        linewidth=0,
        edgecolor="none",
        facecolor=fc,
    )
    ax.add_patch(rect)
    ax.text(
        x + w / 2,
        y + h / 2,
        wrapped.get(text, text),
        ha="center",
        va="center",
        color=COLORS["ink"],
        linespacing=1.08,
        fontsize=6.0,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rebuild Figure 2 from frozen derived inputs.")
    parser.add_argument("--output", type=Path, default=FIGURE_ROOT / "reproduced/Figure2")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    two_axis = pd.read_csv(V2 / "02_rankings/KG_metabolism_aware_v2_two_axis.tsv", sep="\t")
    alpha = pd.read_csv(V2 / "03_validation/alpha_sensitivity_summary.tsv", sep="\t")
    scenario = pd.read_csv(V2 / "03_validation/evidence_scenario_summary.tsv", sep="\t")
    null = pd.read_csv(V2 / "03_validation/topology_null_1000.tsv", sep="\t")
    null_summary = pd.read_csv(V2 / "03_validation/topology_null_summary.tsv", sep="\t")
    agora_summary = pd.read_csv(AGORA_SUMMARY)

    assert len(two_axis) == 104
    assert set(alpha["alpha"]) == {0.0, 0.25, 0.5, 0.75, 1.0}
    assert len(null) == 1000
    counts = dict(zip(agora_summary["metric"], agora_summary["value"]))
    assert int(counts["source_chain_rows"]) == 8063
    assert int(counts["taxon_routes"]) == 300
    assert int(counts["formal_flux_routes"]) == 26
    assert int(counts["t2d_formal_flux_routes"]) == 9

    # 180 × 132 mm at the final publication size.
    fig, axes = plt.subplots(2, 3, figsize=(7.0866, 5.1969), constrained_layout=False)
    # Manual geometry keeps the six data rectangles exactly aligned. The
    # generous gutters accommodate the categorical labels in panels b and e.
    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.10, top=0.94, wspace=0.62, hspace=0.47)

    # a, implementation and evidence gate
    ax = axes[0, 0]
    panel_label(ax, "a")
    ax.set_title("CID-connected transformation axis", loc="left", fontweight="bold", pad=5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    box(ax, 0.08, 0.70, 0.84, 0.18, "104 materials · 30,800 herb–CID pairs", COLORS["light_blue"])
    box(ax, 0.08, 0.43, 0.84, 0.18, "95 terminal tuples → 1,694 material–route paths", COLORS["light_teal"])
    box(ax, 0.08, 0.16, 0.84, 0.18, "Primary subset: 572 paths · 35 transformations · 85 materials", COLORS["light_orange"])
    ax.text(0.5, 0.655, "↓", ha="center", va="center", fontsize=7, color=COLORS["muted"])
    ax.text(0.5, 0.385, "↓", ha="center", va="center", fontsize=7, color=COLORS["muted"])
    ax.text(0.5, 0.04, "Conflicting host paths retained at zero weight", ha="center", va="center", fontsize=5.4, color=COLORS["muted"])

    # b, AGORA2 feasibility gate
    ax = axes[0, 1]
    panel_label(ax, "b")
    ax.set_title("Conditional AGORA2 feasibility", loc="left", fontweight="bold", pad=5)
    labels = ["Unique taxon routes", "≥1 F1 strain", "T2D-core F1 routes"]
    vals = [300, 26, 9]
    y = np.arange(3)
    ax.barh(y, vals, color=[COLORS["grey"], COLORS["flux"], COLORS["transform"]], height=0.55)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 360)
    ax.set_xlabel("Route count")
    ax.grid(axis="x", color=COLORS["grid"], lw=0.5)
    ax.set_axisbelow(True)
    for yi, value in zip(y, vals):
        ax.text(value + 7, yi, str(value), va="center", ha="left", fontweight="bold", color=COLORS["ink"])

    # c, hero two-axis comparison
    ax = axes[0, 2]
    panel_label(ax, "c")
    ax.set_title("Two axes diverged (ρ = −0.032)", loc="left", fontweight="bold", pad=5)
    tier20 = two_axis["tiered_research_priority_rank"] <= 20
    ax.scatter(two_axis.loc[~tier20, "direct_adjusted_rank"], two_axis.loc[~tier20, "primary_transform_rank"], s=9, color=COLORS["grey"], alpha=0.7, linewidth=0)
    ax.scatter(two_axis.loc[tier20, "direct_adjusted_rank"], two_axis.loc[tier20, "primary_transform_rank"], s=17, color=COLORS["transform"], alpha=0.9, linewidth=0, label="Tiered v2 top 20")
    ax.plot([1, 104], [1, 104], color=COLORS["grid"], lw=0.8, ls="--", zorder=0)
    ax.set_xlim(0, 105)
    ax.set_ylim(105, 0)
    ax.set_xlabel("Degree-adjusted direct rank")
    ax.set_ylabel("Transformation-evidence rank")
    ax.grid(color=COLORS["grid"], lw=0.4)
    ax.set_axisbelow(True)

    # d, alpha sensitivity
    ax = axes[1, 0]
    panel_label(ax, "d")
    ax.set_title("Transformation weight changed the top 20", loc="left", fontweight="bold", pad=5)
    ax.plot(alpha["alpha"], alpha["top20_overlap_with_direct_adjusted"], marker="o", ms=4, lw=1.4, color=COLORS["transform"])
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(0, 22.5)
    ax.set_xticks(alpha["alpha"])
    ax.set_xlabel("Transformation contribution, α")
    ax.set_ylabel("Top-20 overlap with direct comparator")
    ax.grid(axis="y", color=COLORS["grid"], lw=0.5)

    # e, evidence scenario sensitivity
    ax = axes[1, 1]
    panel_label(ax, "e")
    ax.set_title("Priorities depended on evidence rules", loc="left", fontweight="bold", pad=5)
    display = ["Strict human + LOTUS", "Any herb occurrence", "+ mixed human/mouse", "All microbial", "Inclusive non-conflict"]
    values = scenario["top20_overlap_vs_primary"].astype(int).to_numpy()
    ypos = np.arange(len(display))
    ax.barh(ypos, values, color=[COLORS["transform"], "#5B9C93", "#82AEA7", "#A8C3BE", COLORS["grey"]], height=0.58)
    ax.set_yticks(ypos, display)
    ax.invert_yaxis()
    ax.set_xlim(0, 22)
    ax.set_xlabel("Top-20 overlap with strict primary set")
    ax.grid(axis="x", color=COLORS["grid"], lw=0.5)

    # f, topology null
    ax = axes[1, 2]
    panel_label(ax, "f")
    ax.set_title("Observed relation matched null", loc="left", fontweight="bold", pad=5)
    rho = null["spearman_direct_vs_transform"].to_numpy()
    observed = float(null_summary.loc[null_summary["metric"] == "spearman_direct_vs_transform", "observed"].iloc[0])
    null_mean = float(null_summary.loc[null_summary["metric"] == "spearman_direct_vs_transform", "null_mean"].iloc[0])
    ax.hist(rho, bins=22, color=COLORS["light_blue"], edgecolor=COLORS["direct"], linewidth=0.45)
    ax.axvline(observed, color=COLORS["transform"], lw=1.5)
    ax.axvline(null_mean, color=COLORS["muted"], lw=1.0, ls="--")
    ax.set_xlabel("Spearman ρ under permutation")
    ax.set_ylabel("Permutations")
    ax.text(0.98, 0.95, "observed −0.032\nnull mean −0.009", transform=ax.transAxes, ha="right", va="top", fontsize=5.4, color=COLORS["ink"])

    for ax in axes.flat:
        ax.tick_params(length=2.5, width=0.55, color=COLORS["muted"])

    stem = out / "Figure2_metabolism_aware_v2"
    alignment = require_matplotlib_panel_alignment(
        fig,
        axes=list(axes.flat),
        panel_ids=list("abcdef"),
        json_out=str(stem) + ".alignment.json",
        overlay_svg=str(stem) + ".alignment.svg",
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        require_panel_labels=True,
        strict=True,
    )

    fig.savefig(str(stem) + ".pdf")
    fig.savefig(str(stem) + ".svg")
    fig.savefig(str(stem) + ".tiff", dpi=600)
    fig.savefig(str(stem) + "_preview.png", dpi=300)
    plt.close(fig)

    source_rows = []
    for _, row in alpha.iterrows():
        source_rows.append({"panel": "d", "metric": f"alpha_{row['alpha']}", "value": row["top20_overlap_with_direct_adjusted"], "source": str(V2 / "03_validation/alpha_sensitivity_summary.tsv")})
    for _, row in scenario.iterrows():
        source_rows.append({"panel": "e", "metric": row["scenario"], "value": row["top20_overlap_vs_primary"], "source": str(V2 / "03_validation/evidence_scenario_summary.tsv")})
    for metric in ["source_chain_rows", "taxon_routes", "formal_flux_routes", "t2d_formal_flux_routes"]:
        source_rows.append({"panel": "b", "metric": metric, "value": counts[metric], "source": str(AGORA_SUMMARY.relative_to(CONTENT))})
    source_rows.extend(
        [
            {"panel": "a", "metric": "materials", "value": 104, "source": str(V2 / "04_reports/METABOLISM_AWARE_V2_REPORT.md")},
            {"panel": "a", "metric": "terminal_routes", "value": 95, "source": str(V2 / "01_routes/transformation_routes_audited_v2.tsv")},
            {"panel": "a", "metric": "material_route_paths", "value": 1694, "source": str(V2 / "01_routes/herb_transformation_paths_v2.tsv")},
            {"panel": "a", "metric": "primary_paths", "value": 572, "source": str(V2 / "01_routes/herb_transformation_paths_v2.tsv")},
            {"panel": "c", "metric": "direct_transform_spearman", "value": observed, "source": str(V2 / "03_validation/topology_null_summary.tsv")},
            {"panel": "f", "metric": "null_iterations", "value": len(null), "source": str(V2 / "03_validation/topology_null_1000.tsv")},
        ]
    )
    for row in source_rows:
        source_path = Path(row["source"])
        if source_path.is_absolute() and source_path.is_relative_to(CONTENT):
            row["source"] = str(source_path.relative_to(CONTENT))
    pd.DataFrame(source_rows).to_csv(out / "Figure2_source_data_map.tsv", sep="\t", index=False)
    (out / "Figure2_build_manifest.json").write_text(
        json.dumps(
            {
                "backend": "Python/matplotlib",
                "final_width_mm": 180,
                "final_height_mm": 132,
                "alignment_verdict": alignment.get("verdict"),
                "no_clinical_tuning": True,
                "interpretation": "Hypothesis prioritisation and conditional model feasibility only",
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
