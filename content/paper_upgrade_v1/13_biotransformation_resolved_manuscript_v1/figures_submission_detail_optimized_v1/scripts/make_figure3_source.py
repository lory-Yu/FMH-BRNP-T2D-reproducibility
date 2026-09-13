#!/usr/bin/env python3
"""Create a submission-grade figure from frozen BRNP evidence-audit outputs."""

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
AUDIT = CONTENT / "paper_upgrade_v1/kg_evidence_sensitivity_20260908_110558"
ENZYME = CONTENT / "paper_upgrade_v1/evidence_completion_human_gut_20260908_005251"
OUT = FIGURE_ROOT / "reproduced/Figure3"
sys.path.insert(0, str(SCRIPT_DIR))
from audit_panel_alignment import require_matplotlib_panel_alignment  # noqa: E402


COLORS = {
    "ink": "#263238",
    "muted": "#68757D",
    "grid": "#DCE3E6",
    "teal": "#2F8F83",
    "teal_dark": "#246F68",
    "teal_light": "#DDEDEA",
    "blue": "#5D7F8F",
    "blue_light": "#E4EDF1",
    "amber": "#B77842",
    "amber_light": "#F4E7DA",
    "plum": "#866A8B",
    "grey": "#BBC6CB",
    "grey_light": "#F1F4F5",
    "white": "#FFFFFF",
}

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 6.3,
        "axes.titlesize": 7.3,
        "axes.labelsize": 6.4,
        "xtick.labelsize": 5.6,
        "ytick.labelsize": 5.6,
        "axes.linewidth": 0.65,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    }
)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.annotate(
        label,
        xy=(0, 1),
        xycoords="axes fraction",
        xytext=(-18, 10),
        textcoords="offset points",
        fontsize=8.2,
        fontweight="bold",
        va="top",
        ha="left",
        color=COLORS["ink"],
        annotation_clip=False,
    )


def rounded_box(ax: plt.Axes, xy, width, height, face, edge="none", lw=0.6) -> None:
    ax.add_patch(
        mpl.patches.FancyBboxPatch(
            xy,
            width,
            height,
            boxstyle="round,pad=0.012,rounding_size=0.022",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
        )
    )


def load_and_check():
    scenario = pd.read_csv(AUDIT / "02_ranking/02_ranking_scenario_summary.tsv", sep="\t")
    shifts = pd.read_csv(AUDIT / "02_ranking/03_herb_rank_shifts.tsv", sep="\t")
    coverage = pd.read_csv(AUDIT / "03_agora2/01_full_model_coverage_scan.tsv", sep="\t")
    bglx = pd.read_csv(ENZYME / "cleaned/11_BglX_exact_sequence_record.tsv", sep="\t")
    conflict = pd.read_csv(ENZYME / "cleaned/10_BbBgl_coordinate_primer_conflict.tsv", sep="\t")
    qa = json.loads((AUDIT / "logs/final_qa.json").read_text(encoding="utf-8"))

    assert qa["unique_cids_audited"] == 10911
    assert qa["pubchem_problem_cids"] == 288
    assert qa["pubchem_resolved"] == 258
    assert qa["unresolved_after_pubchem"] == 30
    assert len(coverage) == 7302 and coverage["read_ok"].all()
    assert int(coverage["ec_generic_21"].sum()) == 4747
    assert int(coverage["target_positive"].sum()) == 0
    assert len(bglx) == 1 and int(bglx.iloc[0]["protein_length_aa"]) == 787
    assert set(conflict["matched_protein_id"]) == {"AHJ23176.1", "AHJ23177.1"}
    return scenario, shifts, coverage, bglx, conflict, qa


def draw_panel_a(ax: plt.Axes, qa: dict) -> None:
    panel_label(ax, "a")
    ax.set_title("Chemical-identity adjudication", loc="left", fontweight="bold", pad=7)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    rounded_box(ax, (0.04, 0.70), 0.92, 0.18, COLORS["blue_light"])
    ax.text(0.50, 0.81, f"{qa['unique_cids_audited']:,} unique CIDs", ha="center", va="center", fontsize=7.2, fontweight="bold", color=COLORS["ink"])
    ax.text(0.50, 0.745, "parent compounds and route products", ha="center", va="center", fontsize=5.4, color=COLORS["muted"])

    ax.annotate("", xy=(0.50, 0.63), xytext=(0.50, 0.70), arrowprops={"arrowstyle": "-|>", "lw": 0.7, "color": COLORS["muted"]})
    rounded_box(ax, (0.20, 0.48), 0.60, 0.15, COLORS["amber_light"])
    ax.text(0.50, 0.565, f"{qa['pubchem_problem_cids']:,} flagged CIDs", ha="center", va="center", fontsize=6.8, fontweight="bold", color=COLORS["ink"])
    ax.text(0.50, 0.515, "missing, unparsable or internally discordant", ha="center", va="center", fontsize=5.2, color=COLORS["muted"])

    ax.plot([0.50, 0.50, 0.28, 0.28], [0.48, 0.43, 0.43, 0.37], color=COLORS["muted"], lw=0.7)
    ax.plot([0.50, 0.50, 0.72, 0.72], [0.48, 0.43, 0.43, 0.37], color=COLORS["muted"], lw=0.7)
    rounded_box(ax, (0.08, 0.21), 0.40, 0.16, COLORS["teal_light"])
    rounded_box(ax, (0.52, 0.21), 0.40, 0.16, COLORS["grey_light"], edge=COLORS["grey"])
    ax.text(0.28, 0.305, f"{qa['pubchem_resolved']:,} resolved", ha="center", va="center", fontsize=6.8, fontweight="bold", color=COLORS["teal_dark"])
    ax.text(0.28, 0.253, "by PubChem record", ha="center", va="center", fontsize=5.3, color=COLORS["muted"])
    ax.text(0.72, 0.305, f"{qa['unresolved_after_pubchem']:,} unresolved", ha="center", va="center", fontsize=6.8, fontweight="bold", color=COLORS["ink"])
    ax.text(0.72, 0.253, "retained as audit flags", ha="center", va="center", fontsize=5.3, color=COLORS["muted"])
    ax.text(0.50, 0.08, "10,881 CIDs strictly usable after adjudication", ha="center", va="center", fontsize=6.0, fontweight="bold", color=COLORS["ink"])


def draw_panel_b(ax: plt.Axes, scenario: pd.DataFrame, shifts: pd.DataFrame) -> None:
    panel_label(ax, "b")
    ax.set_title("Entity rules changed rank stability", loc="left", fontweight="bold", pad=7)

    order = ["local_consensus_strict", "pubchem_resolved_strict"]
    labels = ["Local conflict-free", "PubChem-adjudicated"]
    for i, key in enumerate(order):
        sub = shifts.loc[shifts["scenario"] == key].copy()
        vals = sub["rank_delta_vs_original"].abs().to_numpy()
        # Deterministic offsets reveal coincident observations without changing data values.
        jitter = ((np.arange(len(vals)) % 13) - 6) * 0.0055
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=8, color=COLORS["grey"], alpha=0.55, linewidth=0, zorder=2)
        top_idx = np.argsort(vals)[-5:]
        ax.scatter(np.full(5, i) + jitter[top_idx], vals[top_idx], s=15, color=COLORS["amber"] if i == 0 else COLORS["teal"], linewidth=0, zorder=3)
        med = float(np.median(vals))
        ax.plot([i - 0.14, i + 0.14], [med, med], lw=1.5, color=COLORS["ink"], zorder=4)

        row = scenario.loc[scenario["scenario"] == key].iloc[0]
        note = f"Top20 {int(row['top20_overlap_vs_original'])}/20\nρ={row['spearman_rank_vs_original']:.3f}"
        ax.text(
            i,
            42.0,
            note,
            ha="center",
            va="top",
            fontsize=5.5,
            color=COLORS["ink"],
            linespacing=1.25,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.6},
            zorder=5,
        )

    ax.set_xticks([0, 1], labels)
    ax.set_xlim(-0.45, 1.45)
    ax.set_ylim(-1, 43)
    ax.set_ylabel("Absolute rank shift from primary BRNP")
    ax.grid(axis="y", color=COLORS["grid"], lw=0.5)
    ax.set_axisbelow(True)


def draw_panel_c(ax: plt.Axes, bglx: pd.DataFrame, conflict: pd.DataFrame) -> None:
    panel_label(ax, "c")
    ax.set_title("Rb1 route: one sequence-anchored step", loc="left", fontweight="bold", pad=7)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    xs = [0.08, 0.35, 0.62, 0.87]
    names = ["Rb1", "Rd", "F2", "Compound K"]
    faces = [COLORS["blue_light"], COLORS["teal_light"], COLORS["grey_light"], COLORS["grey_light"]]
    widths = [0.13, 0.13, 0.13, 0.23]
    for x, name, face, width in zip(xs, names, faces, widths):
        ax.add_patch(mpl.patches.Ellipse((x, 0.58), width=width, height=0.13, facecolor=face, edgecolor=COLORS["muted"], lw=0.6))
        ax.text(x, 0.58, name, ha="center", va="center", fontsize=5.8, fontweight="bold", color=COLORS["ink"])

    ax.annotate("", xy=(0.285, 0.58), xytext=(0.145, 0.58), arrowprops={"arrowstyle": "-|>", "lw": 1.4, "color": COLORS["teal_dark"]})
    ax.plot([0.415, 0.555], [0.58, 0.58], color=COLORS["amber"], lw=1.1, ls=(0, (3, 2)))
    ax.plot([0.685, 0.745], [0.58, 0.58], color=COLORS["amber"], lw=1.1, ls=(0, (3, 2)))

    protein = bglx.iloc[0]["protein_accession"]
    ax.text(0.215, 0.76, f"BglX · {protein}", ha="center", va="center", fontsize=5.7, fontweight="bold", color=COLORS["teal_dark"])
    ax.text(0.215, 0.69, "recombinant enzyme + HPLC", ha="center", va="center", fontsize=5.1, color=COLORS["muted"])
    ax.text(0.62, 0.75, "BbBgl activity reported", ha="center", va="center", fontsize=5.7, fontweight="bold", color=COLORS["amber"])
    ax.text(0.62, 0.69, "gene identity unresolved", ha="center", va="center", fontsize=5.1, color=COLORS["muted"])
    ax.text(0.62, 0.43, "Coordinates → AHJ23177.1\nPrimers → AHJ23176.1", ha="center", va="center", fontsize=5.2, color=COLORS["ink"], linespacing=1.25)

    rounded_box(ax, (0.14, 0.12), 0.72, 0.15, COLORS["grey_light"], edge=COLORS["grey"])
    ax.text(0.50, 0.205, "Current BRNP scoring paths", ha="center", va="center", fontsize=5.3, color=COLORS["muted"])
    ax.text(0.50, 0.155, "0 Rb1 / Rd / F2 / Compound K routes", ha="center", va="center", fontsize=6.1, fontweight="bold", color=COLORS["ink"])


def draw_panel_d(ax: plt.Axes, coverage: pd.DataFrame) -> None:
    panel_label(ax, "d")
    ax.set_title("AGORA2 route coverage (7,302 models)", loc="left", fontweight="bold", pad=7)
    labels = ["Models scanned", "Generic EC 3.2.1.21", "Target compounds", "Route-specific ECs"]
    values = [len(coverage), int(coverage["ec_generic_21"].sum()), int(coverage["any_target_entity_identifier_or_name"].sum()), int(coverage["any_route_specific_ec"].sum())]
    y = np.arange(len(labels))
    cols = [COLORS["blue"], COLORS["plum"], COLORS["amber"], COLORS["amber"]]
    for yi, value, color in zip(y, values, cols):
        ax.hlines(yi, 0, value, color=COLORS["grid"], lw=2.0, zorder=1)
        ax.scatter(value, yi, s=34 if value else 30, facecolor=color if value else COLORS["white"], edgecolor=color, lw=1.1, zorder=3)
        if value > 6500:
            label_x, label_ha = value - 180, "right"
        else:
            label_x, label_ha = value + (220 if value else 130), "left"
        if yi == len(labels) - 1:
            label_y, label_va = yi - 0.11, "bottom"
        else:
            label_y, label_va = yi + 0.11, "top"
        ax.text(
            label_x,
            label_y,
            f"{value:,}",
            ha=label_ha,
            va=label_va,
            fontsize=6.1,
            fontweight="bold",
            color=COLORS["ink"],
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.5},
            zorder=5,
        )
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(-100, 8000)
    ax.set_xlabel("Number of AGORA2 models")
    ax.set_axisbelow(True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rebuild Figure 3 from frozen evidence-audit outputs.")
    parser.add_argument("--output", type=Path, default=OUT)
    return parser.parse_args()


def main() -> None:
    global OUT
    OUT = parse_args().output.resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    scenario, shifts, coverage, bglx, conflict, qa = load_and_check()

    fig, axes = plt.subplots(2, 2, figsize=(7.0866, 5.1969), constrained_layout=False)
    fig.subplots_adjust(left=0.095, right=0.985, bottom=0.095, top=0.94, wspace=0.34, hspace=0.46)

    draw_panel_a(axes[0, 0], qa)
    draw_panel_b(axes[0, 1], scenario, shifts)
    draw_panel_c(axes[1, 0], bglx, conflict)
    draw_panel_d(axes[1, 1], coverage)

    for ax in axes.flat:
        ax.tick_params(length=2.5, width=0.55, color=COLORS["muted"])

    stem = OUT / "Figure_BRNP_entity_Rb1_evidence_audit"
    alignment = require_matplotlib_panel_alignment(
        fig,
        axes=list(axes.flat),
        panel_ids=list("abcd"),
        json_out=str(stem) + ".alignment.json",
        overlay_svg=str(stem) + ".alignment.svg",
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        require_panel_labels=True,
        strict=True,
    )

    fig.savefig(str(stem) + ".pdf")
    fig.savefig(str(stem) + ".svg")
    fig.savefig(str(stem) + ".png", dpi=600)
    fig.savefig(str(stem) + ".tiff", dpi=600)
    fig.savefig(str(stem) + "_preview.png", dpi=300)
    plt.close(fig)

    source_rows = [
        ("a", "unique_cids", qa["unique_cids_audited"], "kg_evidence_sensitivity", "01_chemistry/10_cid_consensus_final.tsv"),
        ("a", "flagged_cids", qa["pubchem_problem_cids"], "kg_evidence_sensitivity", "logs/final_qa.json"),
        ("a", "resolved_cids", qa["pubchem_resolved"], "kg_evidence_sensitivity", "logs/final_qa.json"),
        ("a", "unresolved_cids", qa["unresolved_after_pubchem"], "kg_evidence_sensitivity", "logs/final_qa.json"),
        ("a", "strictly_usable_cids", qa["unique_cids_audited"] - qa["unresolved_after_pubchem"], "kg_evidence_sensitivity", "logs/final_qa.json"),
        ("c", "BglX_nucleotide_accession", bglx.iloc[0]["nucleotide_accession"], "human_gut_evidence", "cleaned/11_BglX_exact_sequence_record.tsv"),
        ("c", "BglX_protein_accession", bglx.iloc[0]["protein_accession"], "human_gut_evidence", "cleaned/11_BglX_exact_sequence_record.tsv"),
        ("c", "BglX_supported_step", bglx.iloc[0]["directly_supported_step"], "human_gut_evidence", "cleaned/11_BglX_exact_sequence_record.tsv"),
        ("c", "BbBgl_coordinate_defined_protein", conflict.loc[conflict["definition"].str.startswith("coordinates"), "matched_protein_id"].iloc[0], "human_gut_evidence", "cleaned/10_BbBgl_coordinate_primer_conflict.tsv"),
        ("c", "BbBgl_primer_defined_protein", conflict.loc[conflict["definition"].str.startswith("exact primer"), "matched_protein_id"].iloc[0], "human_gut_evidence", "cleaned/10_BbBgl_coordinate_primer_conflict.tsv"),
        ("c", "current_scored_ginsenoside_routes", 0, "kg_evidence_sensitivity", "logs/ranking_sensitivity_summary.json"),
        ("d", "agora_models_scanned", len(coverage), "kg_evidence_sensitivity", "03_agora2/01_full_model_coverage_scan.tsv"),
        ("d", "generic_ec_3_2_1_21_models", int(coverage["ec_generic_21"].sum()), "kg_evidence_sensitivity", "03_agora2/01_full_model_coverage_scan.tsv"),
        ("d", "target_entity_models", int(coverage["any_target_entity_identifier_or_name"].sum()), "kg_evidence_sensitivity", "03_agora2/01_full_model_coverage_scan.tsv"),
        ("d", "route_specific_ec_models", int(coverage["any_route_specific_ec"].sum()), "kg_evidence_sensitivity", "03_agora2/01_full_model_coverage_scan.tsv"),
    ]
    for key in ["local_consensus_strict", "pubchem_resolved_strict"]:
        row = scenario.loc[scenario["scenario"] == key].iloc[0]
        source_rows.extend(
            [
                ("b", f"{key}_paths", int(row["n_paths_used"]), "kg_evidence_sensitivity", "02_ranking/02_ranking_scenario_summary.tsv"),
                ("b", f"{key}_routes", int(row["n_routes_used"]), "kg_evidence_sensitivity", "02_ranking/02_ranking_scenario_summary.tsv"),
                ("b", f"{key}_top20_overlap", int(row["top20_overlap_vs_original"]), "kg_evidence_sensitivity", "02_ranking/02_ranking_scenario_summary.tsv"),
                ("b", f"{key}_rho", float(row["spearman_rank_vs_original"]), "kg_evidence_sensitivity", "02_ranking/02_ranking_scenario_summary.tsv"),
            ]
        )
    source_df = pd.DataFrame(source_rows, columns=["panel", "metric", "value", "source_package", "source_relative_path"])
    package_roots = {"kg_evidence_sensitivity": AUDIT, "human_gut_evidence": ENZYME}
    source_df["source_repository_path"] = [
        str((package_roots[p] / rel).relative_to(CONTENT))
        for p, rel in zip(source_df["source_package"], source_df["source_relative_path"])
    ]
    source_df.to_csv(OUT / "Figure_BRNP_entity_Rb1_source_data.tsv", sep="\t", index=False)

    (OUT / "Figure_BRNP_entity_Rb1_build_manifest.json").write_text(
        json.dumps(
            {
                "backend": "Python/matplotlib",
                "final_width_mm": 180,
                "final_height_mm": 132,
                "alignment_verdict": alignment.get("verdict"),
                "visual_jitter": "deterministic index offsets; data values unchanged",
                "data_rows_excluded": 0,
                "causal_claim": False,
                "interpretation": "Evidence curation, sequence anchoring and model coverage; not enzyme or clinical validation",
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
