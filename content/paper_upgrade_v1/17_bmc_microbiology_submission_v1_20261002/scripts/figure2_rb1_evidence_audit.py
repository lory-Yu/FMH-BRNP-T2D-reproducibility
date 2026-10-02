#!/usr/bin/env python3
"""Figure 2: public-record and model-coverage audit of the Rb1 conversion chain."""

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
SOURCE = ROOT / "source_data" / "figure2"
OUT = ROOT / "figures"
QA = ROOT / "qa"
STEM = "Figure2_Rb1_evidence_and_model_coverage"

MM = 1 / 25.4
WIDTH_MM = 170
WIDTH_IN = WIDTH_MM * MM
HEIGHT_IN = 137 * MM

INK = "#27353B"
MUTED = "#63727A"
GRID = "#DCE3E4"
BG = "#FAFBFA"
TEAL = "#4D7F7A"
TEAL_LIGHT = "#DCEAE7"
BLUE = "#587698"
BLUE_LIGHT = "#E1E9F1"
GOLD = "#B27D3E"
GOLD_LIGHT = "#F3E7D5"
RED = "#A55C52"
RED_LIGHT = "#F2DEDA"
GREY = "#AEB8BB"
GREY_LIGHT = "#EEF1F1"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def rounded_box(ax, xy, width, height, facecolor, edgecolor, radius=0.022, lw=0.9):
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


def arrow(ax, start, end, color=MUTED, lw=1.0, linestyle="solid", style="-|>"):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            transform=ax.transAxes,
            arrowstyle=style,
            mutation_scale=9,
            linewidth=lw,
            linestyle=linestyle,
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
        ax.text(0.035, 0.88, subtitle, transform=ax.transAxes, fontsize=6.7, color=MUTED, va="top")


def compound(ax, center, label, fill, edge):
    x, y = center
    rounded_box(ax, (x - 0.065, y - 0.11), 0.13, 0.22, fill, edge, radius=0.018)
    ax.text(x, y, label, transform=ax.transAxes, ha="center", va="center", fontsize=8.5,
            fontweight="bold", color=edge)


def panel_a(ax, bglx):
    ax.set_axis_off()
    label_panel(
        ax,
        "a",
        "BglX directly supports one biochemical step",
        "Exact public sequence identity and recombinant-enzyme HPLC support do not establish the complete Rb1-to-CK chain.",
    )
    y = 0.48
    xs = [0.12, 0.38, 0.64, 0.88]
    compound(ax, (xs[0], y), "Rb1", BLUE_LIGHT, BLUE)
    compound(ax, (xs[1], y), "Rd", TEAL_LIGHT, TEAL)
    compound(ax, (xs[2], y), "F2", GREY_LIGHT, GREY)
    compound(ax, (xs[3], y), "CK", GREY_LIGHT, GREY)

    arrow(ax, (xs[0] + 0.075, y), (xs[1] - 0.075, y), color=TEAL, lw=2.0)
    arrow(ax, (xs[1] + 0.075, y), (xs[2] - 0.075, y), color=GREY, lw=1.2, linestyle=(0, (3, 2)))
    arrow(ax, (xs[2] + 0.075, y), (xs[3] - 0.075, y), color=GREY, lw=1.2, linestyle=(0, (3, 2)))

    ax.text(0.25, 0.70, "directly supported", transform=ax.transAxes, ha="center", fontsize=6.7,
            fontweight="bold", color=TEAL)
    ax.text(0.25, 0.26, "BglX · HQ875470.1 / ADY62498.1\nBifidobacterium longum H-1 · purified enzyme · HPLC",
            transform=ax.transAxes, ha="center", va="top", fontsize=6.2, color=INK, linespacing=1.2)
    ax.text(0.64, 0.70, "not established by this enzyme record", transform=ax.transAxes, ha="center",
            fontsize=6.4, color=MUTED)
    ax.text(0.88, 0.26, "No native-strain knockout or\ncomplementation in this record", transform=ax.transAxes,
            ha="center", va="top", fontsize=6.1, color=MUTED)


def gene_arrow(ax, start, end, y, height, color, alpha=1.0):
    # Coordinates decrease in the transcriptional direction (minus strand).
    ax.add_patch(
        FancyArrowPatch(
            (end, y),
            (start, y),
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=height,
            color=color,
            alpha=alpha,
            transform=ax.transData,
            capstyle="round",
        )
    )


def panel_b(ax, rows):
    label_panel(
        ax,
        "b",
        "Reported coordinates and exact primers identify adjacent coding sequences",
        "The public record does not resolve a single BbBgl gene identity.",
    )
    ax.set_xlim(1_364_300, 1_368_450)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="x", labelsize=6.1, colors=MUTED, length=2.5, pad=2)
    ax.ticklabel_format(style="plain", axis="x")
    ax.set_xlabel("CP006715.1 genomic coordinate (bp)", fontsize=6.4, color=MUTED, labelpad=4)
    ax.axhline(0.50, color=GRID, lw=1.0, zorder=0)

    primer = next(row for row in rows if row["both_primers_exact"] == "True")
    coordinate = next(row for row in rows if row["both_primers_exact"] == "False")
    p_start, p_end = 1_364_534, 1_366_786
    c_start, c_end = 1_366_947, 1_368_227
    gene_arrow(ax, p_start, p_end, 0.46, 10, TEAL)
    gene_arrow(ax, c_start, c_end, 0.54, 10, RED)

    ax.text((p_start + p_end) / 2, 0.69, f"{primer['matched_locus_tag']} · {primer['matched_protein_id']}",
            ha="center", va="center", fontsize=6.5, fontweight="bold", color=TEAL)
    ax.text((p_start + p_end) / 2, 0.61, "β-glucosidase; primers 0 + 0 mismatches",
            ha="center", va="center", fontsize=5.9, color=TEAL)
    ax.text((c_start + c_end) / 2, 0.33, f"{coordinate['matched_locus_tag']} · {coordinate['matched_protein_id']}",
            ha="center", va="center", fontsize=6.5, fontweight="bold", color=RED)
    ax.text((c_start + c_end) / 2, 0.24, "hypothetical protein; primers 10 + 11 mismatches",
            ha="center", va="center", fontsize=5.9, color=RED)
    ax.text(1_366_860, 0.79, "adjacent CDSs", ha="center", fontsize=6.2, color=MUTED)
    ax.annotate("", xy=(1_366_785, 0.70), xytext=(1_366_945, 0.70),
                arrowprops={"arrowstyle": "|-|", "color": MUTED, "lw": 0.8})
    ax.text(1_368_360, 0.08, "minus strand ←", ha="right", fontsize=6.0, color=MUTED)


def panel_c(ax, scan):
    ax.set_axis_off()
    label_panel(
        ax,
        "c",
        "AGORA2 contains generic hydrolase annotation but no explicit Rb1-chain representation",
        "Complete local-release scan: 7,302/7,302 SBML models read successfully.",
    )
    total = int(scan["xml_files_scanned"])
    generic = int(scan["models_with_generic_EC_3_2_1_21_token"])
    remainder = total - generic

    rounded_box(ax, (0.035, 0.23), 0.19, 0.46, BLUE_LIGHT, BLUE)
    ax.text(0.13, 0.54, f"{total:,}", transform=ax.transAxes, ha="center", va="center",
            fontsize=18, fontweight="bold", color=BLUE)
    ax.text(0.13, 0.36, "models scanned\n0 read failures", transform=ax.transAxes, ha="center", va="center",
            fontsize=6.7, color=INK, linespacing=1.2)

    rounded_box(ax, (0.275, 0.44), 0.19, 0.25, RED_LIGHT, RED)
    ax.text(0.315, 0.565, "0", transform=ax.transAxes, ha="center", va="center",
            fontsize=13, fontweight="bold", color=RED)
    ax.text(0.395, 0.565, "matches to prespecified\nRb1/Rd/F2/CK identifiers", transform=ax.transAxes,
            ha="center", va="center", fontsize=5.7, color=INK, linespacing=1.05)

    rounded_box(ax, (0.275, 0.14), 0.19, 0.23, RED_LIGHT, RED)
    ax.text(0.315, 0.255, "0", transform=ax.transAxes, ha="center", va="center",
            fontsize=13, fontweight="bold", color=RED)
    ax.text(0.395, 0.255, "specific EC\n3.2.1.191–195", transform=ax.transAxes,
            ha="center", va="center", fontsize=5.7, color=INK, linespacing=1.05)

    ax.text(0.53, 0.69, "Generic β-glucosidase annotation (EC 3.2.1.21)", transform=ax.transAxes,
            fontsize=6.8, fontweight="bold", color=INK, va="center")
    bar_x, bar_y, bar_w, bar_h = 0.53, 0.48, 0.42, 0.13
    rounded_box(ax, (bar_x, bar_y), bar_w, bar_h, GREY_LIGHT, GREY, radius=0.012, lw=0.7)
    teal_w = bar_w * generic / total
    rounded_box(ax, (bar_x, bar_y), teal_w, bar_h, TEAL, TEAL, radius=0.012, lw=0)
    ax.text(bar_x + teal_w / 2, bar_y + bar_h / 2, f"{generic:,}", transform=ax.transAxes,
            ha="center", va="center", fontsize=7.0, fontweight="bold", color="white")
    ax.text(bar_x + teal_w + (bar_w - teal_w) / 2, bar_y + bar_h / 2, f"{remainder:,}", transform=ax.transAxes,
            ha="center", va="center", fontsize=6.6, color=MUTED)
    ax.text(bar_x, 0.41, f"{generic/total:.1%} of scanned models", transform=ax.transAxes,
            fontsize=6.2, color=TEAL, ha="left")

    rounded_box(ax, (0.53, 0.14), 0.42, 0.18, GOLD_LIGHT, GOLD, radius=0.018)
    ax.text(0.74, 0.23, "A model-coverage gap is not evidence that the\nrepresented strains lack real ginsenoside metabolism.",
            transform=ax.transAxes, ha="center", va="center", fontsize=6.3, color=INK, linespacing=1.2)


def write_layout_manifest(fig, axes, path: Path):
    fig.canvas.draw()
    width_pt, height_pt = [float(v) * 72 for v in fig.get_size_inches()]
    panels = []
    for idx, (letter, ax) in enumerate(zip("abc", axes)):
        p = ax.get_position()
        panels.append(
            {
                "id": letter,
                "bbox_pt": [p.x0 * width_pt, p.y0 * height_pt, p.x1 * width_pt, p.y1 * height_pt],
                "grid_id": "figure2-grid",
                "row_start": idx,
                "row_stop": idx + 1,
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
    bglx = read_tsv(SOURCE / "11_BglX_exact_sequence_record.tsv")[0]
    conflict = read_tsv(SOURCE / "10_BbBgl_coordinate_primer_conflict.tsv")
    scan = json.loads((SOURCE / "agora2_full_scan_summary.json").read_text(encoding="utf-8"))

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
    grid = fig.add_gridspec(3, 1, left=0.06, right=0.985, top=0.965, bottom=0.055, hspace=0.34)
    axes = [fig.add_subplot(grid[i, 0]) for i in range(3)]
    for ax in axes:
        ax.set_facecolor(BG)
    panel_a(axes[0], bglx)
    panel_b(axes[1], conflict)
    panel_c(axes[2], scan)
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
