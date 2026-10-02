#!/usr/bin/env python3
"""Figure 1: evidence gates across microbial transformation reconstruction.

The figure deliberately separates three operations with different counting units:
record deduplication, taxon-to-strain model expansion, and evidence filtering.
Areas of boxes are decorative and do not encode magnitude.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT / "source_data" / "figure1" / "12_analysis_summary.csv"
OUT = ROOT / "figures"
QA = ROOT / "qa"
STEM = "Figure1_evidence_gating_microbiome"

MM = 1 / 25.4
WIDTH_MM = 170
WIDTH_IN = WIDTH_MM * MM
HEIGHT_IN = 137 * MM

INK = "#26343A"
MUTED = "#607078"
GRID = "#DDE4E5"
BG = "#FAFBFA"
TEAL = "#4E7E7B"
TEAL_LIGHT = "#DCE9E7"
BLUE = "#58759A"
BLUE_LIGHT = "#E2E9F1"
GOLD = "#B68142"
GOLD_LIGHT = "#F2E7D6"
CORAL = "#A95F50"
CORAL_LIGHT = "#F2DFDA"


def read_metrics(path: Path) -> dict[str, int]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["metric"]: int(float(row["value"])) for row in csv.DictReader(handle)}


def rounded_box(ax, xy, width, height, facecolor, edgecolor, radius=0.025, lw=0.9):
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        transform=ax.transAxes,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=lw,
        clip_on=False,
    )
    ax.add_patch(patch)
    return patch


def arrow(ax, start, end, color=MUTED, lw=1.0, style="-|>"):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            transform=ax.transAxes,
            arrowstyle=style,
            mutation_scale=9,
            linewidth=lw,
            color=color,
            shrinkA=1.5,
            shrinkB=1.5,
            clip_on=False,
        )
    )


def label_panel(ax, letter: str, title: str, subtitle: str | None = None):
    ax.text(-0.012, 1.035, letter, transform=ax.transAxes, fontsize=10.5, fontweight="bold", color=INK, va="top")
    ax.text(0.035, 1.035, title, transform=ax.transAxes, fontsize=9.2, fontweight="bold", color=INK, va="top")
    if subtitle:
        ax.text(0.035, 0.925, subtitle, transform=ax.transAxes, fontsize=6.7, color=MUTED, va="top")


def panel_a(ax, m):
    ax.set_axis_off()
    label_panel(
        ax,
        "a",
        "Database records contract after route-level deduplication",
        "Counting unit changes from composed source records to unique taxon–substrate–product routes.",
    )

    rounded_box(ax, (0.08, 0.25), 0.30, 0.43, BLUE_LIGHT, BLUE)
    ax.text(0.23, 0.575, f"{m['source_chain_rows']:,}", ha="center", va="center", transform=ax.transAxes,
            fontsize=20, fontweight="bold", color=BLUE)
    ax.text(0.23, 0.355, "database-composed\nsource records", ha="center", va="center", transform=ax.transAxes,
            fontsize=7.5, color=INK, linespacing=1.2)

    arrow(ax, (0.41, 0.47), (0.60, 0.47), color=MUTED, lw=1.3)
    ax.text(0.505, 0.57, "group identical", ha="center", va="center", transform=ax.transAxes,
            fontsize=6.5, color=MUTED)
    ax.text(0.505, 0.38, "taxon + substrate + product", ha="center", va="center", transform=ax.transAxes,
            fontsize=6.3, color=MUTED)

    rounded_box(ax, (0.63, 0.25), 0.29, 0.43, TEAL_LIGHT, TEAL)
    ax.text(0.775, 0.575, f"{m['taxon_routes']:,}", ha="center", va="center", transform=ax.transAxes,
            fontsize=20, fontweight="bold", color=TEAL)
    ax.text(0.775, 0.355, "unique microbial\ntransformation routes", ha="center", va="center", transform=ax.transAxes,
            fontsize=7.5, color=INK, linespacing=1.2)

    ax.text(0.50, 0.11, "Deduplication is a contraction; box area is not proportional to count.",
            ha="center", va="center", transform=ax.transAxes, fontsize=6.4, color=MUTED)


def panel_b(ax, m):
    ax.set_axis_off()
    label_panel(
        ax,
        "b",
        "Taxonomic mapping expands routes into strain-model evaluations",
        "This is a coverage expansion, not a screening funnel; each stage uses a different unit.",
    )

    cards = [
        (0.025, 0.28, f"{m['taxa']:,}", "source taxon\nstrings", BLUE_LIGHT, BLUE),
        (0.275, 0.28, f"{m['mapped_taxa']:,}", "mapped taxa", TEAL_LIGHT, TEAL),
        (0.525, 0.28, f"{m['strain_models']:,}", "AGORA2 strain\nmodels", GOLD_LIGHT, GOLD),
        (0.775, 0.28, f"{m['loaded_strain_route_rows']:,}", "loaded strain–route\ncases", CORAL_LIGHT, CORAL),
    ]
    for x, y, number, label, fill, edge in cards:
        rounded_box(ax, (x, y), 0.19, 0.39, fill, edge, radius=0.02)
        ax.text(x + 0.095, y + 0.255, number, ha="center", va="center", transform=ax.transAxes,
                fontsize=14.5, fontweight="bold", color=edge)
        ax.text(x + 0.095, y + 0.115, label, ha="center", va="center", transform=ax.transAxes,
                fontsize=6.7, color=INK, linespacing=1.15)
    arrow(ax, (0.218, 0.475), (0.273, 0.475))
    arrow(ax, (0.468, 0.475), (0.523, 0.475))
    arrow(ax, (0.718, 0.475), (0.773, 0.475))
    ax.text(0.245, 0.73, "taxonomy\nresolution", ha="center", transform=ax.transAxes, fontsize=5.9, color=MUTED)
    ax.text(0.495, 0.73, "model\nlookup", ha="center", transform=ax.transAxes, fontsize=5.9, color=MUTED)
    ax.text(0.745, 0.73, "route × strain\nenumeration", ha="center", transform=ax.transAxes, fontsize=5.9, color=MUTED)
    ax.text(
        0.50,
        0.12,
        f"{m['mapped_taxa']}/{m['taxa']} taxa mapped; {m['no_model_placeholder_rows']} no-model route placeholders retained "
        f"({m['matrix_rows']:,} matrix rows total).",
        ha="center",
        transform=ax.transAxes,
        fontsize=6.1,
        color=MUTED,
    )


def panel_c(ax, m):
    ax.set_axis_off()
    label_panel(
        ax,
        "c",
        "Chemical and flux criteria leave a small model-supported subset",
        "Counts remain strain–route combinations until the final collapse to unique routes.",
    )

    stages = [
        ("Loaded strain–route\ncases", m["loaded_strain_route_rows"], CORAL),
        ("Both structures\nverified", m["both_structures_verified"], GOLD),
        ("Formal F1 flux\ncombinations", m["formal_flux_rows"], TEAL),
    ]
    max_value = m["loaded_strain_route_rows"]
    y_positions = [0.68, 0.46, 0.24]
    x0, max_width = 0.19, 0.49
    for (label, value, color), y in zip(stages, y_positions):
        width = max_width * value / max_value
        ax.text(0.02, y + 0.045, label, transform=ax.transAxes, fontsize=6.7, color=INK,
                va="center", ha="left", linespacing=1.05)
        rounded_box(ax, (x0, y), width, 0.095, color, color, radius=0.012, lw=0)
        ax.text(x0 + width + 0.012, y + 0.0475, f"{value:,}", transform=ax.transAxes,
                fontsize=7.2, fontweight="bold", color=color, va="center", ha="left")

    arrow(ax, (0.69, 0.29), (0.765, 0.29), color=MUTED, lw=1.0)
    ax.text(0.675, 0.53, "collapse to unique routes", transform=ax.transAxes, fontsize=5.7,
            color=MUTED, ha="center", va="center")

    rounded_box(ax, (0.775, 0.49), 0.19, 0.23, TEAL_LIGHT, TEAL, radius=0.02)
    ax.text(0.815, 0.605, f"{m['formal_flux_routes']}", transform=ax.transAxes, ha="center", va="center",
            fontsize=13, fontweight="bold", color=TEAL)
    ax.text(0.895, 0.605, "routes with ≥1\nformal F1 result", transform=ax.transAxes, ha="center", va="center",
            fontsize=5.9, color=INK, linespacing=1.1)

    rounded_box(ax, (0.775, 0.18), 0.19, 0.23, GOLD_LIGHT, GOLD, radius=0.02)
    ax.text(0.815, 0.295, f"{m['t2d_formal_flux_routes']}", transform=ax.transAxes, ha="center", va="center",
            fontsize=13, fontweight="bold", color=GOLD)
    ax.text(0.895, 0.295, "routes linked to\nprespecified T2D core", transform=ax.transAxes, ha="center", va="center",
            fontsize=5.8, color=INK, linespacing=1.1)
    arrow(ax, (0.87, 0.485), (0.87, 0.415), color=MUTED, lw=0.9)

    ax.text(
        0.02,
        0.105,
        f"One E. coli tryptophan→indole route contributed {m['dominant_f1_route_rows']:,}/{m['formal_flux_rows']:,} "
        "F1-positive strain cases (93.3%); route counts are the primary summary.",
        transform=ax.transAxes,
        fontsize=5.6,
        color=INK,
        ha="left",
        va="bottom",
        fontstyle="italic",
    )
    ax.text(
        0.02,
        0.02,
        "F1 required eligible taxonomy, verified structures, native exchange reactions and increased product export under curated anaerobic model conditions.",
        transform=ax.transAxes,
        fontsize=5.9,
        color=MUTED,
        ha="left",
        va="bottom",
    )


def write_layout_manifest(fig, axes, path: Path):
    fig.canvas.draw()
    width_pt, height_pt = [float(v) * 72 for v in fig.get_size_inches()]
    panels = []
    for letter, ax in zip("abc", axes):
        p = ax.get_position()
        panels.append(
            {
                "id": letter,
                "bbox_pt": [p.x0 * width_pt, p.y0 * height_pt, p.x1 * width_pt, p.y1 * height_pt],
                "grid_id": "figure1-grid",
                "row_start": ord(letter) - ord("a"),
                "row_stop": ord(letter) - ord("a") + 1,
                "col_start": 0,
                "col_stop": 1,
                "panel_label": letter,
                "panel_label_anchor_pt": [(p.x0 - 0.012 * p.width) * width_pt, (p.y1 + 0.035 * p.height) * height_pt],
            }
        )
    manifest = {
        "schema_version": 1,
        "backend": "python-matplotlib",
        "figure": {"width_pt": width_pt, "height_pt": height_pt},
        "panels": panels,
        "column_groups": [{"id": "single-column", "panels": ["a", "b", "c"]}],
        "exemptions": [],
    }
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def require_matplotlib_panel_alignment(fig, axes, path: Path):
    """Render-time gate: all stacked panels must share left and right bounds."""
    fig.canvas.draw()
    positions = [ax.get_position() for ax in axes]
    tolerance = 1.5 / (float(fig.get_size_inches()[0]) * 72)
    lefts = [p.x0 for p in positions]
    rights = [p.x1 for p in positions]
    if max(lefts) - min(lefts) > tolerance or max(rights) - min(rights) > tolerance:
        raise RuntimeError("Stacked panel bounds fail the 1.5-point alignment gate")
    write_layout_manifest(fig, axes, path)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    m = read_metrics(SOURCE)

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 7,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )
    fig = plt.figure(figsize=(WIDTH_IN, HEIGHT_IN), facecolor="white")
    grid = fig.add_gridspec(3, 1, left=0.055, right=0.985, top=0.965, bottom=0.04, hspace=0.30)
    axes = [fig.add_subplot(grid[i, 0]) for i in range(3)]
    for ax in axes:
        ax.set_facecolor(BG)
    panel_a(axes[0], m)
    panel_b(axes[1], m)
    panel_c(axes[2], m)

    require_matplotlib_panel_alignment(fig, axes, QA / f"{STEM}_layout.json")
    fig.savefig(OUT / f"{STEM}.pdf", bbox_inches=None, pad_inches=0)
    fig.savefig(OUT / f"{STEM}.svg", bbox_inches=None, pad_inches=0)
    fig.savefig(OUT / f"{STEM}.png", dpi=600, bbox_inches=None, pad_inches=0)
    fig.savefig(
        OUT / f"{STEM}.tiff",
        dpi=600,
        bbox_inches=None,
        pad_inches=0,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


if __name__ == "__main__":
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig_fmh_bmc")
    main()
