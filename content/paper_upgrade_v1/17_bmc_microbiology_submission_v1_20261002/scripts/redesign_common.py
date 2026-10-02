#!/usr/bin/env python3
"""Shared Nature-style helpers for Figure 2–4 redesign v1."""
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
from matplotlib import font_manager, transforms
from matplotlib.font_manager import FontProperties

PACKAGE = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("FMH_FIGURE_OUT", str(PACKAGE / "figures"))).resolve()
SOURCE = PACKAGE / "source_data"
QA = OUT / "qa"
SCRIPTS = OUT / "scripts"

VENDOR = Path(__file__).resolve().parents[1] / "vendor"
sys.path.insert(0, str(VENDOR))
from audit_panel_alignment import require_matplotlib_panel_alignment  # noqa: E402

CJK_PATH = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
SANS_PATH = Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")
if CJK_PATH.exists():
    font_manager.fontManager.addfont(str(CJK_PATH))
if SANS_PATH.exists():
    font_manager.fontManager.addfont(str(SANS_PATH))

FP_CJK = FontProperties(fname=str(CJK_PATH)) if CJK_PATH.exists() else FontProperties()
FP_SANS = FontProperties(fname=str(SANS_PATH)) if SANS_PATH.exists() else FontProperties()

# Nature clinical / genomics restrained palette (2026 observations)
C = {
    "ink": "#1C1C1C",
    "mute": "#5E5E5E",
    "faint": "#9A9A9A",
    "rule": "#D0D0D0",
    "grid": "#ECECEC",
    "band": "#F5F5F5",
    "blue": "#2F5D8A",
    "blue_soft": "#A8C1D6",
    "teal": "#3A7573",
    "teal_soft": "#B7D0CE",
    "warm": "#B5673D",
    "warm_soft": "#E5C4AE",
    "amber": "#A67C2D",
    "green": "#3F7A5A",
    "red": "#A94442",
    "high": "#B5673D",
    "low": "#2F5D8A",
    "orig": "#A94442",
    "adj": "#3A7573",
}

# Short English display names for figure axes (common / pharmacopoeia usage).
HERB_EN = {
    "人参": "Ginseng",
    "丁香": "Clove",
    "白果": "Ginkgo seed",
    "黄芪": "Astragalus",
    "蜂蜜": "Honey",
    "山茱萸": "Cornus",
    "灵芝": "Ganoderma",
    "枸杞子": "Goji berry",
    "大枣": "Jujube",
    "砂仁": "Amomum",
    "西红花": "Saffron",
    "乌梅": "Mume",
    "金银花": "Honeysuckle",
    "桑椹": "Mulberry",
    "党参": "Codonopsis",
    "西洋参": "American ginseng",
    "玳玳花": "Bitter-orange flower",
    "茯苓": "Poria",
    "沙棘": "Sea buckthorn",
    "肉豆蔻": "Nutmeg",
    "桑叶": "Mulberry leaf",
    "菊花": "Chrysanthemum",
    "荷叶": "Lotus leaf",
    "木瓜": "Chaenomeles",
    "姜黄": "Turmeric",
}


def herb_en(name: str) -> str:
    key = str(name).strip()
    if key in HERB_EN:
        return HERB_EN[key]
    if all(ord(ch) < 128 for ch in key):
        return key
    return f"Herb-{abs(hash(key)) % 10000:04d}"


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial", "sans-serif"],
        "font.size": 7.5,
        "axes.titlesize": 8.8,
        "axes.labelsize": 7.6,
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "legend.fontsize": 6.8,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.unicode_minus": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    }
)

# Type scale used by figure scripts
TYPE = {
    "figure_title": 10.5,
    "panel_title": 8.8,
    "axis": 7.6,
    "tick": 7.0,
    "anno": 6.6,
    "note": 6.4,
}


def ensure_dirs() -> None:
    for d in (OUT, SOURCE, QA, SCRIPTS):
        d.mkdir(parents=True, exist_ok=True)


def panel_label(ax, label: str) -> None:
    offset = transforms.ScaledTranslation(-14 / 72, 4 / 72, ax.figure.dpi_scale_trans)
    ax.text(
        0,
        1.0,
        label,
        transform=ax.transAxes + offset,
        fontsize=11,
        fontweight="bold",
        va="bottom",
        ha="right",
        color=C["ink"],
        clip_on=False,
        fontproperties=FP_SANS,
    )


def style(ax, grid: str | None = "y") -> None:
    ax.tick_params(width=0.5, length=2.2, colors=C["mute"])
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#8A8A8A")
    if grid:
        ax.grid(axis=grid, color=C["grid"], lw=0.55, zorder=0)
        ax.set_axisbelow(True)


def draw_stats_cards(
    ax,
    y_positions,
    texts,
    *,
    header: str | None = None,
    header_y: float | None = None,
    face_cycle=None,
    edge=None,
    x0: float = 0.02,
    width: float = 0.96,
    half_h: float = 0.38,
    fontsize: float = 6.0,
    header_fontsize: float = 6.4,
    linespacing: float = 1.12,
) -> None:
    """Draw side-column numbers inside alternating light colored rounded cards."""
    from matplotlib import patches as _patches

    faces = face_cycle or [C["band"], "#EEF3F6", "#F3F0EA", "#EEF3F6"]
    edge = "none" if edge is None else edge
    ax.set_xlim(0, 1)
    ax.axis("off")
    if header is not None:
        hy = header_y if header_y is not None else (float(max(y_positions)) + half_h + 0.22)
        ax.text(0.50, hy, header, va="bottom", ha="center", fontsize=header_fontsize, color=C["mute"], fontproperties=FP_SANS, clip_on=False)
    for i, (yi, text) in enumerate(zip(y_positions, texts)):
        face = faces[i % len(faces)]
        ax.add_patch(
            _patches.FancyBboxPatch(
                (x0, yi - half_h),
                width,
                2 * half_h,
                boxstyle="round,pad=0.012,rounding_size=0.06",
                facecolor=face,
                edgecolor=edge,
                lw=0.55,
                clip_on=False,
                zorder=1,
            )
        )
        ax.text(
            0.50,
            yi,
            text,
            va="center",
            ha="center",
            fontsize=fontsize,
            color=C["ink"],
            fontproperties=FP_SANS,
            linespacing=linespacing,
            zorder=2,
            clip_on=False,
        )


def note_card(ax, x, y, text, *, transform=None, fontsize=6.2, ha="left", va="top", face="#F3F6F7", edge=None) -> None:
    """Opaque padded note so stats do not collide with strokes/points."""
    transform = transform or ax.transAxes
    ax.text(
        x,
        y,
        text,
        transform=transform,
        va=va,
        ha=ha,
        fontsize=fontsize,
        color=C["ink"],
        fontproperties=FP_SANS,
        linespacing=1.15,
        clip_on=False,
        zorder=5,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": face,
            "edgecolor": edge or C["rule"],
            "linewidth": 0.45,
            "alpha": 0.97,
        },
    )


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def save_figure(
    fig,
    stem: str,
    *,
    axes,
    panel_ids,
    row_groups,
    column_groups,
    sources: list[Path],
    caption: str,
    boundary: str,
    notes: list[str] | None = None,
    exemptions: list[dict] | None = None,
) -> Path:
    ensure_dirs()
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig,
        axes=axes,
        panel_ids=panel_ids,
        row_groups=row_groups,
        column_groups=column_groups,
        exemptions=exemptions or [],
        json_out=QA / f"{stem}_alignment.json",
        overlay_svg=QA / f"{stem}_alignment_overlay.svg",
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

    legend = OUT / f"{stem}_legend.txt"
    legend.write_text(caption + "\n\nBoundary: " + boundary + "\n", encoding="utf-8")
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stem": stem,
        "caption": caption,
        "boundary": boundary,
        "notes": notes or [],
        "sources": [{"path": str(p), "sha256": sha256(p)} for p in sources if p.exists()],
        "outputs": [{"path": str(p), "sha256": sha256(p)} for p in paths.values()],
    }
    (OUT / f"{stem}_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return paths["pdf"]
