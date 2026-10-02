#!/usr/bin/env python3
"""Build BMC Microbiology Figure 4 from frozen Wei source tables.

High and Low are source-defined converter groups.  This same-source public-data
reanalysis is not labelled as an independent validation, enzyme assay, or
measurement of Rb1-to-Compound-K flux.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


PACKAGE = Path(__file__).resolve().parents[1]
BASE_DIR = Path(__file__).resolve().parent
BASE_SCRIPT = BASE_DIR / "make_figure4_evidence_completion.py"
OUT = Path(os.environ.get("FMH_FIGURE_OUT", str(PACKAGE / "figures"))).resolve()
SOURCE = PACKAGE / "source_data" / "figure4"

sys.path.insert(0, str(BASE_DIR))
spec = importlib.util.spec_from_file_location("wei_figure_base", BASE_SCRIPT)
if spec is None or spec.loader is None:
    raise RuntimeError("Could not load the frozen Wei figure builder")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
import redesign_common  # noqa: E402


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial", "sans-serif"],
        "font.size": 7.5,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)


def save_figure_local(
    fig,
    stem,
    *,
    axes,
    panel_ids,
    row_groups,
    column_groups,
    sources,
    caption,
    boundary,
    notes=None,
    exemptions=None,
):
    """Use the established renderer but keep this revision's export contract explicit."""
    caption = (
        "Figure 4. Microbiome features in source-defined Rb1 converter groups. "
        "(a) Linkage of 50 donor labels to public reads, DADA2 quality control, retained ASVs and the primary High-versus-Low contrast; the 30 Middle donors were not included in panels b-e. "
        "(b) Individual-donor observed-ASV richness for ten Low and ten High converters; horizontal and vertical black marks show the median and interquartile range. The footer reports the Wilcoxon test and age- and sex-adjusted model with BH correction across four alpha-diversity metrics. "
        "(c) Bray-Curtis principal-coordinate analysis with PERMANOVA, ANOSIM and PERMDISP results. "
        "(d) Age- and sex-adjusted centred-log-ratio coefficients and 95% confidence intervals for six detected source-defined candidate taxonomic labels; q values use candidate-family BH correction. "
        "(e) Global BH q values for all 115 tested taxonomic features. Teal diamonds identify source-defined candidates rather than statistical significance; 0/115 passed global BH q < 0.05. High and Low are labels from the source study, and this is a same-source reanalysis rather than independent validation."
    )
    boundary = (
        "The panels describe 16S taxonomic composition associated with source-defined converter labels; "
        "they do not measure beta-glucosidase activity or donor-level Rb1-to-Compound-K kinetics."
    )
    OUT.mkdir(parents=True, exist_ok=True)
    qa_dir = OUT / "qa"
    qa_dir.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    redesign_common.require_matplotlib_panel_alignment(
        fig,
        axes=axes,
        panel_ids=panel_ids,
        row_groups=row_groups,
        column_groups=column_groups,
        exemptions=exemptions or [],
        json_out=qa_dir / f"{stem}_alignment.json",
        overlay_svg=qa_dir / f"{stem}_alignment_overlay.svg",
        tolerance_pt=2.0,
        gutter_tolerance_pt=2.0,
        require_panel_labels=True,
        strict=True,
    )
    paths = {
        "pdf": OUT / f"{stem}.pdf",
        "svg": OUT / f"{stem}.svg",
        "tiff": OUT / f"{stem}.tiff",
        "png_600dpi": OUT / f"{stem}.png",
        "preview": OUT / f"{stem}_preview.png",
    }
    fig.savefig(paths["pdf"], facecolor="white")
    fig.savefig(paths["svg"], facecolor="white")
    fig.savefig(paths["tiff"], dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(paths["png_600dpi"], dpi=600, facecolor="white")
    fig.savefig(paths["preview"], dpi=320, facecolor="white")
    plt.close(fig)

    (OUT / f"{stem}_legend.txt").write_text(caption + "\n\nBoundary: " + boundary + "\n", encoding="utf-8")
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stem": stem,
        "caption": caption,
        "boundary": boundary,
        "notes": notes or [],
        "sources": [{"path": str(p), "sha256": redesign_common.sha256(p)} for p in sources if p.exists()],
        "outputs": [{"path": str(p), "sha256": redesign_common.sha256(p)} for p in paths.values()],
        "final_size_mm": [170.0, 193.55],
        "panel_e_change": "Display semantics only; all 115 global BH q values retained unchanged.",
    }
    (OUT / f"{stem}_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return paths["pdf"]


def plot_qrank_clarified(ax, global_df, cand):
    """Preserve global BH ordering and clarify that colour denotes prior candidacy."""
    df = global_df.sort_values(["adjusted_q", "adjusted_p", "genus"], kind="mergesort").reset_index(drop=True).copy()
    df["display_rank"] = np.arange(1, len(df) + 1)
    df["neglog10q"] = -np.log10(df["adjusted_q"].clip(lower=1e-12))
    is_cand = df["genus"].isin(["Anaerostipes", "Coprococcus", "Barnesiella"]) | df["genus"].astype(str).str.contains(
        "Oscillospiraceae|coprostanoligenes|AD3011", case=False, na=False
    )
    assert int(is_cand.sum()) == 6

    ax.scatter(
        df.loc[~is_cand, "display_rank"],
        df.loc[~is_cand, "neglog10q"],
        s=12,
        c="#BDBDBD",
        edgecolors="none",
        zorder=2,
    )
    ax.scatter(
        df.loc[is_cand, "display_rank"],
        df.loc[is_cand, "neglog10q"],
        s=31,
        marker="D",
        facecolors=base.C["teal"],
        edgecolors="white",
        linewidth=0.5,
        zorder=3,
    )

    threshold = -np.log10(0.05)
    ax.axhline(threshold, color=base.C["red"], ls="--", lw=0.85)
    ax.text(
        0.98,
        threshold - 0.18,
        "Global BH q = 0.05",
        transform=ax.get_yaxis_transform(),
        color=base.C["red"],
        fontsize=6.0,
        va="top",
        ha="right",
        fontproperties=base.FP_SANS,
        clip_on=False,
    )

    outlier = df.loc[df["genus"].astype(str).str.contains("Oscillospiraceae", case=False, na=False)].iloc[0]
    ax.annotate(
        f"Oscillospiraceae\nq={outlier.adjusted_q:.3f}",
        xy=(outlier.display_rank, outlier.neglog10q),
        xytext=(82, 0.20),
        textcoords="data",
        ha="center",
        va="bottom",
        fontsize=5.5,
        color=base.C["teal"],
        fontproperties=base.FP_SANS,
    )

    n_pass = int((df["adjusted_q"] < 0.05).sum())
    ax.set_xlabel("Global BH-q rank (lowest to highest)", fontproperties=base.FP_SANS, fontsize=6.8, labelpad=1)
    ax.set_ylabel("-log10(global BH q)", fontproperties=base.FP_SANS, fontsize=7.4)
    ax.set_title("Global multiplicity across 115 features", loc="left", fontweight="bold", color=base.C["ink"], pad=6, fontsize=8.5, fontproperties=base.FP_SANS)
    ax.set_xlim(-2, float(df["display_rank"].max()) + 4)
    base.style(ax, "y")
    return n_pass, len(df)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    redesign_common.OUT = OUT
    redesign_common.QA = OUT / "qa"
    redesign_common.SCRIPTS = OUT / "scripts"
    base.OUT = OUT
    base.STEM = "Figure4_Wei_converter_groups_BMC"
    base.plot_qrank = plot_qrank_clarified
    base.save_figure = save_figure_local

    def freeze_local():
        files = {
            "alpha": SOURCE / "Figure4_alpha_by_donor.tsv",
            "stats": SOURCE / "Figure4_alpha_stats.tsv",
            "pcoa": SOURCE / "Figure4_pcoa.tsv",
            "cand": SOURCE / "Figure4_candidate_genera.tsv",
            "global": SOURCE / "Figure4_global_genera.tsv",
            "qc": SOURCE / "Figure4_qc_summary.tsv",
            "beta": SOURCE / "Figure4_beta_stats.tsv",
        }
        missing = [str(p) for p in files.values() if not p.exists()]
        if missing:
            raise FileNotFoundError("Missing Figure 4 source tables: " + ", ".join(missing))
        return files

    base.freeze = freeze_local

    original_main = base.main

    # The base function writes the explanatory foot card. Replace only its text
    # by intercepting Axes.text so the numeric panel and all source data remain unchanged.
    original_text = base.plt.Axes.text
    original_figure = base.plt.figure
    original_suptitle = base.plt.Figure.suptitle
    original_add_gridspec = base.plt.Figure.add_gridspec
    original_set_title = base.plt.Axes.set_title

    def clarified_text(self, x, y, s, *args, **kwargs):
        if s == "q":
            s = ""
        if isinstance(s, str) and s.startswith("Teal = source candidates"):
            s = "Teal diamonds = source-defined candidates, not significance\n6/115 highlighted; 0/115 passed global BH q<0.05"
            kwargs["fontsize"] = 5.1
        if isinstance(s, str) and s.startswith("PERMANOVA R2="):
            s = s.replace("PERMANOVA R2=", "PERMANOVA R²=")
        return original_text(self, x, y, s, *args, **kwargs)

    def bmc_figure(*args, **kwargs):
        if "figsize" in kwargs:
            kwargs["figsize"] = (6.6929, 7.62)  # 170 mm final width
        return original_figure(*args, **kwargs)

    def bmc_suptitle(self, t, *args, **kwargs):
        if isinstance(t, str) and "Rb1 converter groups" in t:
            t = "Microbiome features in source-defined Rb1 converter groups"
        kwargs["y"] = 0.985
        kwargs["fontsize"] = 10.2
        return original_suptitle(self, t, *args, **kwargs)

    def bmc_add_gridspec(self, *args, **kwargs):
        if kwargs.get("left") == 0.12:
            kwargs["left"] = 0.18
        if kwargs.get("right") == 0.92:
            kwargs["right"] = 0.98
        if kwargs.get("top") == 0.925:
            kwargs["top"] = 0.90
        return original_add_gridspec(self, *args, **kwargs)

    def bmc_set_title(self, label, *args, **kwargs):
        replacements = {
            "Observed richness by converter group": "Observed richness\nby converter group",
            "Bray-Curtis ordination (High/Low n=10)": "Bray-Curtis ordination\n(High/Low n=10)",
            "Candidate genus associations": "Candidate taxonomic\nassociations",
            "Candidate taxonomic associations": "Candidate taxonomic\nassociations",
            "Global multiplicity across 115 features": "Global multiplicity\nacross 115 features",
        }
        label = replacements.get(label, label)
        return original_set_title(self, label, *args, **kwargs)

    base.plt.Axes.text = clarified_text
    base.plt.figure = bmc_figure
    base.plt.Figure.suptitle = bmc_suptitle
    base.plt.Figure.add_gridspec = bmc_add_gridspec
    base.plt.Axes.set_title = bmc_set_title
    try:
        original_main()
    finally:
        base.plt.Axes.text = original_text
        base.plt.figure = original_figure
        base.plt.Figure.suptitle = original_suptitle
        base.plt.Figure.add_gridspec = original_add_gridspec
        base.plt.Axes.set_title = original_set_title


if __name__ == "__main__":
    main()
