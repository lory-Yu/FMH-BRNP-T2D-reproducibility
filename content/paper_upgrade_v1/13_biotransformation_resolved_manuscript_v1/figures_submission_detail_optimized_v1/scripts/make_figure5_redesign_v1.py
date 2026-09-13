#!/usr/bin/env python3
"""Nature redesign Figure 5: clinical evidence maturity and TRIAL028 sensitivity.

Archetype: clinical quantitative grid (PRISMA-style flow + attribution bars
+ material-level forests + sensitivity forests + RoB 2 traffic lights).

Data: frozen source_data_v8 / TRIAL028 final adjudication tables.
Backend: Python (matplotlib) exclusively.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_ROOT = SCRIPT_DIR.parent
SOURCE = FIGURE_ROOT / "source/figure6"
OUT = Path(os.environ.get("FMH_FIGURE_OUT", str(FIGURE_ROOT / "reproduced/Figure6"))).resolve()
QA = OUT / "qa"
PANELS = OUT / "panels"
STEM = "Figure5_clinical_evidence_TRIAL028_redesign_v1"

sys.path.insert(0, str(SCRIPT_DIR))
from audit_panel_alignment import require_matplotlib_panel_alignment  # noqa: E402

# Reuse English herb map from redesign_common when available
sys.path.insert(0, str(SCRIPT_DIR))
try:
    from redesign_common import HERB_EN, herb_en  # noqa: E402
except Exception:  # pragma: no cover
    HERB_EN = {}

    def herb_en(name: str) -> str:
        return str(name)

SANS = Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")
FP = FontProperties(fname=str(SANS)) if SANS.exists() else FontProperties()

C = {
    "ink": "#1C1C1C",
    "mute": "#5E5E5E",
    "faint": "#9A9A9A",
    "rule": "#D0D0D0",
    "grid": "#ECECEC",
    "band": "#F7F7F7",
    "teal": "#3A7573",
    "blue": "#2F5D8A",
    "warm": "#B5673D",
    "amber": "#A67C2D",
    "purple": "#6B5B8A",
    "pool": "#2F5D8A",
    "explore": "#A67C2D",
    "single": "#5E5E5E",
    "low": "#3F7A5A",
    "some": "#C4A035",
    "high": "#A94442",
    "box": "#F3F6F7",
    "excl": "#F7F1E8",
}

MATURITY_EN = {
    "单味+复方均有": "Single-herb + formula",
    "仅单味": "Single-herb only",
    "仅含药复方信号": "Formula-only signal",
    "两层皆无": "No signal in either layer",
}
MATURITY_ORDER = [
    "Single-herb + formula",
    "Single-herb only",
    "Formula-only signal",
    "No signal in either layer",
]
MATURITY_COL = {
    "Single-herb + formula": C["teal"],
    "Single-herb only": C["blue"],
    "Formula-only signal": C["purple"],
    "No signal in either layer": C["faint"],
}

# Manuscript display order for forests
HERB_ORDER = ["西红花", "人参", "姜黄", "桑叶", "黄芪", "灵芝"]

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 7.2,
        "axes.titlesize": 8.6,
        "axes.labelsize": 7.4,
        "xtick.labelsize": 6.8,
        "ytick.labelsize": 6.8,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.unicode_minus": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    }
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def panel_label(ax, lab: str) -> None:
    off = ScaledTranslation(-12 / 72, 3 / 72, ax.figure.dpi_scale_trans)
    ax.text(
        0,
        1.0,
        lab,
        transform=ax.transAxes + off,
        fontsize=11,
        fontweight="bold",
        ha="right",
        va="bottom",
        color=C["ink"],
        clip_on=False,
        fontproperties=FP,
    )


def style(ax, grid="x") -> None:
    ax.tick_params(width=0.5, length=2.2, colors=C["mute"])
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#8A8A8A")
    if grid:
        ax.grid(axis=grid, color=C["grid"], lw=0.55, zorder=0)
        ax.set_axisbelow(True)


def load_meta() -> pd.DataFrame:
    meta = pd.read_csv(SOURCE / "Figure5_primary_meta.tsv", sep="\t")
    pooled = meta.loc[meta["test"].eq("modified_HK")].copy()
    single = meta.loc[meta["test"].eq("single_study")].drop_duplicates(["herb", "outcome"])
    out = pd.concat([pooled, single], ignore_index=True)
    out["herb_en"] = out["herb"].map(herb_en)
    out["order"] = out["herb"].map({h: i for i, h in enumerate(HERB_ORDER)})
    return out.sort_values(["outcome", "order", "herb_en"])


def write_manifest(paths: dict[str, Path]) -> Path:
    rows = [
        {
            "panel": "5a",
            "title": "Verifiable screening path",
            "chart": "PRISMA-style screening audit flow",
            "source": str(paths["screen"]),
            "sha256": sha256(paths["screen"]),
            "unit_note": "records / reports / independent trials distinguished",
            "status": "verified_from_source_data_v8",
        },
        {
            "panel": "5b",
            "title": "Evidence attribution across 47 materials",
            "chart": "horizontal attribution bar chart",
            "source": str(paths["mat"]),
            "sha256": sha256(paths["mat"]),
            "unit_note": "mutually exclusive maturity classes; n/47",
            "status": "verified_from_source_data_v8",
        },
        {
            "panel": "5c",
            "title": "Primary HbA1c estimates",
            "chart": "material-level MD forest (percentage points)",
            "source": str(paths["meta"]),
            "sha256": sha256(paths["meta"]),
            "unit_note": "modified Knapp-Hartung; k=2 exploratory; k=1 not pooled",
            "status": "verified_from_source_data_v8",
        },
        {
            "panel": "5d",
            "title": "Primary FPG estimates",
            "chart": "material-level MD forest (mg/dL)",
            "source": str(paths["meta"]),
            "sha256": sha256(paths["meta"]),
            "unit_note": "FPG harmonised to mg/dL; independent x-axis from 5c",
            "status": "verified_from_source_data_v8",
        },
        {
            "panel": "5e",
            "title": "Saffron sensitivity to TRIAL028 handling",
            "chart": "three-scenario sensitivity forests",
            "source": str(paths["sens"]),
            "sha256": sha256(paths["sens"]),
            "unit_note": "digitised / ANCOVA / excluded; never both endpoint+ANCOVA",
            "status": "verified_from_source_data_v8",
        },
        {
            "panel": "5f",
            "title": "Outcome-specific RoB 2 for TRIAL028",
            "chart": "robvis-style traffic-light matrix",
            "source": str(paths["rob"]),
            "sha256": sha256(paths["rob"]),
            "unit_note": "pending independent duplicate review (status field)",
            "status": "FINAL_PROJECT_ADJUDICATION_PENDING_INDEPENDENT_DUPLICATE_REVIEW",
        },
    ]
    man = SOURCE / "Figure5_panel_manifest.tsv"
    pd.DataFrame(rows).to_csv(man, sep="\t", index=False)
    return man


def plot_a(ax: plt.Axes, screen: pd.DataFrame) -> None:
    """PRISMA-style screening flow with visible inter-box arrows."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Verifiable screening path", loc="left", fontweight="bold", color=C["ink"], pad=4, fontsize=8.6, fontproperties=FP)

    by_id = {r.stage_id: r for r in screen.itertuples()}

    def n(sid: str) -> int:
        return int(by_id[sid].count)

    # Shorter boxes + larger gaps so arrowheads remain visible
    main = [
        (0.93, f"Multi-database archive\nn={n('U01'):,} records"),
        (0.78, f"Diabetes–herb filter\nn={n('U02'):,} → {n('U02d'):,} after dedup"),
        (0.63, f"Auditable prefilter entry\nn={n('S01'):,} records"),
        (0.48, f"Title/abstract hits\nn={n('S02'):,} rows → {n('S03'):,} records"),
        (0.33, f"Detailed screening\nn={n('S04b'):,} records"),
        (0.18, f"Eligible full texts\nn={n('S05')} papers"),
        (0.04, f"Primary synthesis\nn={n('S06')} trials · {n('S07')} endpoint rows"),
    ]
    excl = [
        (0.78, f"Dedup removed\nn={n('U02') - n('U02d'):,}"),
        (0.48, f"Pre-exclusion\nn={n('S04a')}"),
        (0.18, "Not whole-herb /\nmissing variance"),
    ]

    mx, mw, mh = 0.24, 0.40, 0.062
    ex, ew, eh = 0.66, 0.20, 0.058
    gap = 0.018  # clear space beyond box edge for arrow shaft

    def box(x, y, w, h, text, face, edge):
        ax.add_patch(
            patches.FancyBboxPatch(
                (x - w / 2, y - h / 2),
                w,
                h,
                boxstyle="round,pad=0.005,rounding_size=0.01",
                facecolor=face,
                edgecolor=edge,
                lw=0.8,
                zorder=2,
                clip_on=True,
            )
        )
        ax.text(x, y, text, ha="center", va="center", fontsize=5.4, color=C["ink"], fontproperties=FP, linespacing=1.08, zorder=3, clip_on=True)

    for y, text in main:
        box(mx, y, mw, mh, text, C["box"], C["teal"])
    for i in range(len(main) - 1):
        y0, y1 = main[i][0], main[i + 1][0]
        ax.annotate(
            "",
            xy=(mx, y1 + mh / 2 + gap),
            xytext=(mx, y0 - mh / 2 - gap),
            arrowprops={"arrowstyle": "-|>", "lw": 0.9, "color": C["mute"], "mutation_scale": 11},
            zorder=1,
            clip_on=False,
        )
    for y, text in excl:
        box(ex, y, ew, eh, text, C["excl"], C["amber"])
        ax.plot([mx + mw / 2, ex - ew / 2], [y, y], color=C["faint"], lw=0.65, clip_on=True)


def plot_b(ax: plt.Axes, mat: pd.DataFrame, ax_note: plt.Axes | None = None) -> None:
    counts = mat["evidence_maturity"].map(MATURITY_EN).value_counts()
    labels = [lab for lab in MATURITY_ORDER if lab in counts.index]
    vals = [int(counts[lab]) for lab in labels]
    y = np.arange(len(labels))[::-1]
    colors = [MATURITY_COL[lab] for lab in labels]
    ax.barh(y, vals, color=colors, height=0.62, edgecolor="white", lw=0.4, zorder=2)
    for yi, v in zip(y, vals):
        ax.text(v + 0.5, yi, f"{v}/47", va="center", ha="left", fontsize=6.5, color=C["ink"], fontproperties=FP, zorder=4)
    ax.set_yticks(y, labels, fontproperties=FP, fontsize=6.6)
    ax.tick_params(axis="y", pad=2.5)
    ax.set_xlim(0, max(vals) * 1.35)
    ax.set_xlabel("Number of materials", fontproperties=FP, fontsize=7.0, labelpad=2)
    ax.set_title("Evidence attribution across 47 materials", loc="left", fontweight="bold", color=C["ink"], pad=4, fontsize=8.6, fontproperties=FP)
    style(ax, None)
    ax.tick_params(axis="x", length=2.2, width=0.5, colors=C["mute"])

    note = "Formula signal is not attributable single-material efficacy. Zero ≠ ineffective."
    if ax_note is not None:
        ax_note.axis("off")
        ax_note.text(0.0, 0.55, note, transform=ax_note.transAxes, fontsize=5.6, color=C["mute"], fontproperties=FP, va="center", ha="left")
    else:
        ax.text(0.0, -0.22, note, transform=ax.transAxes, fontsize=5.6, color=C["mute"], fontproperties=FP, clip_on=False)


def _forest_ylab(herb_en: str, k: int) -> str:
    short = {"Mulberry leaf": "Mulberry"}.get(herb_en, herb_en)
    return f"{short} (k={k})"


def plot_forest(
    ax_plot,
    ax_tab,
    meta: pd.DataFrame,
    outcome: str,
    title: str,
    xlab: str,
    xlim: tuple[float, float],
) -> None:
    df = meta.loc[meta["outcome"].eq(outcome)].copy().sort_values("order")
    y = np.arange(len(df))[::-1].astype(float)
    ax_plot.axvline(0, color=C["faint"], ls="--", lw=0.85, zorder=1)
    lines = []
    for yi, row in zip(y, df.itertuples()):
        if row.k >= 3:
            mk, col, ms, tag = "D", C["pool"], 6.0, "Pooled"
        elif row.k == 2:
            mk, col, ms, tag = "s", C["explore"], 5.4, "Exploratory"
        else:
            mk, col, ms, tag = "o", C["single"], 5.2, "Not pooled"
        eb = ax_plot.errorbar(
            row.estimate,
            yi,
            xerr=[[row.estimate - row.ci_lb], [row.ci_ub - row.estimate]],
            fmt=mk,
            color=col,
            ecolor=col,
            elinewidth=1.1,
            capsize=1.8,
            ms=ms,
            zorder=3,
            clip_on=True,
        )
        for artist in (eb[0], *eb[1], *eb[2]):
            artist.set_clip_on(True)
        lines.append(f"{row.estimate:.2f} [{row.ci_lb:.2f}, {row.ci_ub:.2f}]\n{tag}")

    ax_plot.set_yticks(y, [_forest_ylab(r.herb_en, int(r.k)) for r in df.itertuples()], fontproperties=FP, fontsize=6.4)
    ax_plot.tick_params(axis="y", pad=5.0)
    ax_plot.set_xlim(*xlim)
    ax_plot.set_xlabel(xlab, fontproperties=FP, fontsize=6.8, labelpad=2)
    ax_plot.set_title(title, loc="left", fontweight="bold", color=C["ink"], pad=4, fontsize=8.6, fontproperties=FP)
    style(ax_plot, "x")

    ax_tab.set_xlim(0, 1)
    ax_tab.set_ylim(ax_plot.get_ylim())
    ax_tab.axis("off")
    ax_tab.text(0.02, float(y.max()) + 0.45, "MD [95% CI]", fontsize=6.2, color=C["mute"], fontproperties=FP, va="bottom")
    for yi, line in zip(y, lines):
        ax_tab.text(0.02, yi, line, va="center", ha="left", fontsize=5.8, color=C["ink"], fontproperties=FP, linespacing=1.08)


def plot_cd_legend(ax: plt.Axes) -> None:
    """Shared marker key centered under panels c and d."""
    ax.axis("off")
    handles = [
        Line2D([0], [0], marker="D", color=C["pool"], ls="none", ms=5.5, label="Pooled (k≥3)"),
        Line2D([0], [0], marker="s", color=C["explore"], ls="none", ms=5.5, label="Exploratory (k=2)"),
        Line2D([0], [0], marker="o", color=C["single"], ls="none", ms=5.5, label="Single trial (k=1)"),
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        loc="center",
        ncol=3,
        prop=FP,
        fontsize=6.6,
        columnspacing=1.6,
        handletextpad=0.35,
        title="c–d marker key",
        title_fontsize=6.2,
    )


def plot_e(ax_h, ax_ht, ax_f, ax_ft, sens: pd.DataFrame, ax_note: plt.Axes | None = None, ax_title: plt.Axes | None = None) -> None:
    """Sensitivity forests; optional title/note axes keep e block height aligned with f."""
    scenarios = [
        ("primary_digitised_endpoint", "Digitised endpoint"),
        ("ancova_replacement", "ANCOVA replacement"),
        ("complete_exclusion", "TRIAL028 excluded"),
    ]
    y = np.array([2.0, 1.0, 0.0])
    shared = [lab for _, lab in scenarios]

    for ax, ax_tab, outcome, xlim, xlab, head, show_ylab in [
        (ax_h, ax_ht, "HbA1c", (-0.70, 0.55), "HbA1c MD (pp)", "HbA1c", True),
        (ax_f, ax_ft, "FPG", (-22, 16), "FPG MD (mg/dL)", "FPG", False),
    ]:
        ax.axvline(0, color=C["faint"], ls="--", lw=0.85, zorder=1)
        lines = []
        for yi, (sid, _) in zip(y, scenarios):
            row = sens.loc[(sens.scenario == sid) & (sens.outcome == outcome)].iloc[0]
            eb = ax.errorbar(
                row.estimate,
                yi,
                xerr=[[row.estimate - row.ci_lb], [row.ci_ub - row.estimate]],
                fmt="D",
                color=C["teal"],
                ecolor=C["teal"],
                elinewidth=1.15,
                capsize=2.0,
                ms=5.4,
                zorder=3,
                clip_on=True,
                markeredgecolor=C["ink"],
                markeredgewidth=0.4,
            )
            for artist in (eb[0], *eb[1], *eb[2]):
                artist.set_clip_on(True)
            lines.append(f"k={int(row.k)}  {row.estimate:.2f} [{row.ci_lb:.2f}, {row.ci_ub:.2f}]")
        ax.set_yticks(y)
        if show_ylab:
            ax.set_yticklabels(shared, fontproperties=FP, fontsize=6.2)
            ax.tick_params(axis="y", pad=2.5)
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis="y", length=0)
        ax.set_xlim(*xlim)
        ax.set_ylim(-0.4, 2.4)
        ax.set_xlabel(xlab, fontproperties=FP, fontsize=6.5, labelpad=2)
        style(ax, "x")
        ax_tab.set_xlim(0, 1)
        ax_tab.set_ylim(ax.get_ylim())
        ax_tab.axis("off")
        ax_tab.text(0.0, 2.35, "MD [95% CI]", fontsize=5.5, color=C["mute"], fontproperties=FP, va="bottom")
        for yi, line in zip(y, lines):
            ax_tab.text(0.0, yi, line, va="center", ha="left", fontsize=5.6, color=C["ink"], fontproperties=FP)

    if ax_title is not None:
        ax_title.axis("off")
        ax_title.text(0.0, 0.35, "e   Saffron sensitivity to TRIAL028 handling", transform=ax_title.transAxes, fontsize=8.4, fontweight="bold", color=C["ink"], fontproperties=FP, va="center")
    ax_h.text(0.0, 0.98, "HbA1c", transform=ax_h.transAxes, fontsize=6.4, fontweight="bold", color=C["teal"], fontproperties=FP, va="top", clip_on=False)
    ax_f.text(0.0, 0.98, "FPG", transform=ax_f.transAxes, fontsize=6.4, fontweight="bold", color=C["purple"], fontproperties=FP, va="top", clip_on=False)

    note = "Row 1 = primary saffron pools in c/d. All CIs in panel e cross zero."
    if ax_note is not None:
        ax_note.axis("off")
        ax_note.text(0.0, 0.6, note, transform=ax_note.transAxes, fontsize=5.6, color=C["mute"], fontproperties=FP, va="center", ha="left")
    else:
        ax_h.text(0.0, -0.32, note, transform=ax_h.transAxes, fontsize=5.5, color=C["mute"], fontproperties=FP, clip_on=False)


def plot_f(ax: plt.Axes, rob: pd.DataFrame, ax_note: plt.Axes | None = None) -> None:
    domains = [
        ("D1_randomization", "D1"),
        ("D2_deviations", "D2"),
        ("D3_missing_data", "D3"),
        ("D4_measurement", "D4"),
        ("D5_selection", "D5"),
        ("overall", "Overall"),
    ]
    risk_col = {"Low risk": C["low"], "Some concerns": C["some"], "High risk": C["high"]}
    # Avoid "?": use "!" for Some concerns (colour + symbol still redundant)
    risk_sym = {"Low risk": "+", "Some concerns": "!", "High risk": "−"}
    order = ["HbA1c", "FPG"]
    df = rob.set_index("outcome").loc[order]
    ax.set_xlim(-0.5, len(domains) - 0.5)
    ax.set_ylim(-0.55, len(order) - 0.45)
    ax.set_xticks(range(len(domains)))
    ax.set_xticklabels([d[1] for d in domains], fontproperties=FP, fontsize=6.6)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(order)), order, fontproperties=FP, fontsize=7.0)
    ax.invert_yaxis()
    cell = 0.70
    for i, outcome in enumerate(order):
        for j, (col, _) in enumerate(domains):
            val = str(df.loc[outcome, col])
            ax.add_patch(
                patches.FancyBboxPatch(
                    (j - cell / 2, i - cell / 2),
                    cell,
                    cell,
                    boxstyle="round,pad=0.02,rounding_size=0.05",
                    facecolor=risk_col.get(val, C["faint"]),
                    edgecolor="white",
                    lw=0.8,
                )
            )
            ax.text(j, i, risk_sym.get(val, "·"), ha="center", va="center", fontsize=10, fontweight="bold", color="white", fontproperties=FP)
    ax.set_title("Outcome-specific RoB 2 for TRIAL028", loc="left", fontweight="bold", color=C["ink"], pad=8, fontsize=8.6, fontproperties=FP)
    ax.tick_params(length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    handles = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor=risk_col[k], markersize=7.5, label=f"{risk_sym[k]}  {k}")
        for k in ["Low risk", "Some concerns", "High risk"]
    ]
    status = str(df.iloc[0]["status"]) if "status" in df.columns else ""
    note = "Bias risk only (not efficacy)."
    if "PENDING" in status:
        note += " Pending independent duplicate review."
    if ax_note is not None:
        ax_note.axis("off")
        ax_note.legend(handles=handles, frameon=False, loc="upper left", ncol=3, prop=FP, fontsize=5.5, columnspacing=0.8)
        ax_note.text(0.0, 0.15, note, transform=ax_note.transAxes, fontsize=5.2, color=C["mute"], fontproperties=FP, va="bottom")
    else:
        ax.legend(handles=handles, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3, prop=FP, fontsize=5.6)
        ax.text(0.0, -0.42, note, transform=ax.transAxes, fontsize=5.2, color=C["mute"], fontproperties=FP, clip_on=False)


def save_all(fig, axes, panel_ids, sources) -> Path:
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig,
        axes=axes,
        panel_ids=panel_ids,
        row_groups=[
            {"id": "mid", "panels": ["c", "d"]},
        ],
        column_groups=[],
        exemptions=[
            {
                "panels": ["c", "d"],
                "checks": ["panel-width", "column"],
                "reason": "Forests use plot|table splits with unequal plot widths.",
            },
        ],
        json_out=QA / f"{STEM}_alignment.json",
        overlay_svg=QA / f"{STEM}_alignment_overlay.svg",
        tolerance_pt=2.0,
        gutter_tolerance_pt=2.0,
        require_panel_labels=True,
        strict=True,
    )
    paths = {
        "pdf": OUT / f"{STEM}.pdf",
        "svg": OUT / f"{STEM}.svg",
        "tiff": OUT / f"{STEM}.tiff",
        "png": OUT / f"{STEM}_preview.png",
    }
    fig.savefig(paths["pdf"], bbox_inches="tight")
    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["tiff"], dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(paths["png"], dpi=320, bbox_inches="tight")
    plt.close(fig)
    (OUT / f"{STEM}_legend.txt").write_text(
        "Figure 5. Clinical evidence maturity and TRIAL028 sensitivity.\n"
        "(a) Verifiable screening path. (b) Evidence attribution among 47 materials. "
        "(c,d) Primary HbA1c and FPG estimates; shared marker key under c–d. "
        "(e) Saffron TRIAL028 sensitivity; all CIs in panel e cross zero. "
        "(f) Outcome-specific RoB 2 for TRIAL028 (! = Some concerns). "
        "Negative MD favours intervention.\n",
        encoding="utf-8",
    )
    (OUT / f"{STEM}_manifest.json").write_text(
        json.dumps(
            {
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "stem": STEM,
                "sources": [{"path": str(p), "sha256": sha256(p)} for p in sources if Path(p).exists()],
                "outputs": [{"path": str(p), "sha256": sha256(p)} for p in paths.values()],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return paths["pdf"]


def main() -> None:
    for d in (OUT, SOURCE, QA, PANELS):
        d.mkdir(parents=True, exist_ok=True)

    screen = pd.read_csv(SOURCE / "Figure5_screening_flow.tsv", sep="\t")
    mat = pd.read_csv(SOURCE / "Figure5_maturity47.tsv", sep="\t")
    meta = load_meta()
    sens = pd.read_csv(SOURCE / "Figure5_TRIAL028_sensitivity.tsv", sep="\t")
    rob = pd.read_csv(SOURCE / "Figure5_TRIAL028_RoB2.tsv", sep="\t")
    man = write_manifest(
        {
            "screen": SOURCE / "Figure5_screening_flow.tsv",
            "mat": SOURCE / "Figure5_maturity47.tsv",
            "meta": SOURCE / "Figure5_primary_meta.tsv",
            "sens": SOURCE / "Figure5_TRIAL028_sensitivity.tsv",
            "rob": SOURCE / "Figure5_TRIAL028_RoB2.tsv",
        }
    )
    meta.to_csv(SOURCE / "Figure5c_d_primary_meta_display.tsv", sep="\t", index=False)

    # 180 mm wide; 4 logical bands: a|b(+note), c|d(+legend), e(+note)|f
    fig = plt.figure(figsize=(7.205, 11.0))
    outer = fig.add_gridspec(
        3,
        1,
        height_ratios=[1.10, 1.55, 1.20],
        left=0.11,
        right=0.985,
        top=0.955,
        bottom=0.045,
        hspace=0.16,
    )

    # Row 1: a | (b plot + b note) — a height matches b+note
    top = outer[0].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.55)
    ax_a = fig.add_subplot(top[0, 0])
    b_col = top[0, 1].subgridspec(2, 1, height_ratios=[4.6, 1.25], hspace=0.25)
    ax_b = fig.add_subplot(b_col[0, 0])
    ax_b_note = fig.add_subplot(b_col[1, 0])

    # Row 2: c|d forests + centered shared legend under both
    mid = outer[1].subgridspec(2, 1, height_ratios=[12.0, 1.15], hspace=0.10)
    mid_plots = mid[0].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.46)
    c_split = mid_plots[0, 0].subgridspec(1, 2, width_ratios=[2.4, 1.45], wspace=0.10)
    d_split = mid_plots[0, 1].subgridspec(1, 2, width_ratios=[2.4, 1.45], wspace=0.10)
    ax_c = fig.add_subplot(c_split[0, 0])
    ax_c_tab = fig.add_subplot(c_split[0, 1])
    ax_d = fig.add_subplot(d_split[0, 0])
    ax_d_tab = fig.add_subplot(d_split[0, 1])
    ax_cd_leg = fig.add_subplot(mid[1])

    # Row 3: mirrored columns (title/body/note) so e and f share total height
    bot = outer[2].subgridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.38)
    e_col = bot[0, 0].subgridspec(3, 1, height_ratios=[0.72, 4.2, 1.20], hspace=0.22)
    f_col = bot[0, 1].subgridspec(3, 1, height_ratios=[0.72, 4.2, 1.20], hspace=0.22)
    ax_e_title = fig.add_subplot(e_col[0])
    e_body = e_col[1].subgridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.28)
    eh = e_body[0, 0].subgridspec(1, 2, width_ratios=[2.55, 1.25], wspace=0.08)
    ef = e_body[0, 1].subgridspec(1, 2, width_ratios=[2.55, 1.25], wspace=0.08)
    ax_e_h = fig.add_subplot(eh[0, 0])
    ax_e_ht = fig.add_subplot(eh[0, 1])
    ax_e_f = fig.add_subplot(ef[0, 0])
    ax_e_ft = fig.add_subplot(ef[0, 1])
    ax_e_note = fig.add_subplot(e_col[2])
    ax_f_pad = fig.add_subplot(f_col[0]); ax_f_pad.axis("off")
    ax_f = fig.add_subplot(f_col[1])
    ax_f_note = fig.add_subplot(f_col[2])

    for ax, lab in [(ax_a, "a"), (ax_b, "b"), (ax_c, "c"), (ax_d, "d"), (ax_f, "f")]:
        panel_label(ax, lab)

    plot_a(ax_a, screen)
    plot_b(ax_b, mat, ax_b_note)
    plot_forest(ax_c, ax_c_tab, meta, "HbA1c", "Primary HbA1c estimates", "Mean difference (pp)", (-3.1, 2.4))
    plot_forest(ax_d, ax_d_tab, meta, "FPG", "Primary FPG estimates", "Mean difference (mg/dL)", (-72, 72))
    plot_cd_legend(ax_cd_leg)
    plot_e(ax_e_h, ax_e_ht, ax_e_f, ax_e_ft, sens, ax_e_note, ax_e_title)
    plot_f(ax_f, rob, ax_f_note)

    # Align a height to b+note; a left to c
    fig.canvas.draw()
    pos_c = ax_c.get_position()
    pos_b = ax_b.get_position()
    pos_bn = ax_b_note.get_position()
    b_block_y0 = pos_bn.y0
    b_block_y1 = pos_b.y0 + pos_b.height
    ax_a.set_position([pos_c.x0, b_block_y0, pos_b.x0 - 0.055 - pos_c.x0, b_block_y1 - b_block_y0])

    fig.suptitle(
        "Clinical evidence maturity and TRIAL028 sensitivity",
        x=pos_c.x0,
        y=0.988,
        ha="left",
        fontsize=10.8,
        fontweight="bold",
        color=C["ink"],
        fontproperties=FP,
    )

    pdf = save_all(
        fig,
        axes=[ax_c, ax_d],
        panel_ids=list("cd"),
        sources=[
            SOURCE / "Figure5_screening_flow.tsv",
            SOURCE / "Figure5_maturity47.tsv",
            SOURCE / "Figure5_primary_meta.tsv",
            SOURCE / "Figure5_TRIAL028_sensitivity.tsv",
            SOURCE / "Figure5_TRIAL028_RoB2.tsv",
            man,
        ],
    )
    print(pdf)


if __name__ == "__main__":
    main()
