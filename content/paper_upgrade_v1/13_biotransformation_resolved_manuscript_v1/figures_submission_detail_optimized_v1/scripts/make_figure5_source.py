#!/usr/bin/env python3
"""Create a non-destructive Wei figure revision with an unambiguous q-rank panel."""

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


SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_ROOT = SCRIPT_DIR.parent
BASE_DIR = SCRIPT_DIR
BASE_SCRIPT = BASE_DIR / "make_figure4_evidence_completion.py"
OUT = Path(os.environ.get("FMH_FIGURE_OUT", str(FIGURE_ROOT / "reproduced/Figure5"))).resolve()

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
        "png": OUT / f"{stem}_preview.png",
    }
    fig.savefig(paths["pdf"], bbox_inches="tight")
    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["tiff"], dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(paths["png"], dpi=320, bbox_inches="tight")
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
    ax.set_title("Global multiplicity across 115 genera", loc="left", fontweight="bold", color=base.C["ink"], pad=6, fontsize=8.5, fontproperties=base.FP_SANS)
    ax.set_xlim(-2, float(df["display_rank"].max()) + 4)
    base.style(ax, "y")
    return n_pass, len(df)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    redesign_common.OUT = OUT
    redesign_common.QA = OUT / "qa"
    redesign_common.SCRIPTS = OUT / "scripts"
    base.OUT = OUT
    base.STEM = "Figure5_Wei2025_qrank_clarified_v1"
    base.plot_qrank = plot_qrank_clarified
    base.save_figure = save_figure_local

    original_main = base.main

    # The base function writes the explanatory foot card. Replace only its text
    # by intercepting Axes.text so the numeric panel and all source data remain unchanged.
    original_text = base.plt.Axes.text

    def clarified_text(self, x, y, s, *args, **kwargs):
        if isinstance(s, str) and s.startswith("Teal = source candidates"):
            s = "Teal diamonds = source-defined candidates, not significance\n6/115 highlighted; 0/115 passed global BH q<0.05"
            kwargs["fontsize"] = 5.1
        return original_text(self, x, y, s, *args, **kwargs)

    base.plt.Axes.text = clarified_text
    try:
        original_main()
    finally:
        base.plt.Axes.text = original_text


if __name__ == "__main__":
    main()
