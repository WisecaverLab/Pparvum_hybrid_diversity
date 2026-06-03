#!/usr/bin/env python3
"""
plot_haplotypes.py

Visualize haplotype segments from a TSV file.
Columns expected: strain_id, start, stop, chromosome, haplotype

Usage:
    python plot_haplotypes.py input.tsv
    python plot_haplotypes.py input.tsv --output haplotypes.png
    python plot_haplotypes.py input.tsv --output haplotypes.pdf --dpi 300
"""

import sys
import argparse
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np


# ── Argument parsing ──────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tsv", help="Input TSV file")
    parser.add_argument("--output", "-o", default="haplotypes.png",
                        help="Output file (png/pdf/svg). Default: haplotypes.png")
    parser.add_argument("--dpi", type=int, default=150,
                        help="Resolution for raster output. Default: 150")
    parser.add_argument("--linewidth", type=float, default=2,
                        help="Width of each segment line. Default: 6")
    parser.add_argument("--col-strain", default=None,
                        help="Name of the strain-ID column (auto-detected if omitted)")
    parser.add_argument("--col-start", default=None,
                        help="Name of the start column (auto-detected if omitted)")
    parser.add_argument("--col-stop", default=None,
                        help="Name of the stop column (auto-detected if omitted)")
    parser.add_argument("--col-chrom", default=None,
                        help="Name of the chromosome column (auto-detected if omitted)")
    parser.add_argument("--col-haplotype", default=None,
                        help="Name of the haplotype column (auto-detected if omitted)")
    return parser.parse_args()


# ── Column auto-detection ─────────────────────────────────────────────────────

GUESSES = {
    "strain":     ["strain_id", "strain", "sample", "id", "name"],
    "start":      ["start", "begin", "pos_start", "chromstart"],
    "stop":       ["stop", "end", "pos_end", "chromend"],
    "chrom":      ["chromosome", "chrom", "chr", "contig"],
    "haplotype":  ["haplotype", "hap", "haplo", "allele"],
}

def autodetect_col(columns, role, override):
    if override:
        if override in columns:
            return override
        raise ValueError(f"Column '{override}' not found. Available: {list(columns)}")
    cols_lower = {c.lower(): c for c in columns}
    for guess in GUESSES[role]:
        if guess in cols_lower:
            return cols_lower[guess]
    raise ValueError(
        f"Cannot auto-detect the '{role}' column. "
        f"Use --col-{role} to specify it. Available columns: {list(columns)}"
    )


# ── Color palette ─────────────────────────────────────────────────────────────

# Colorblind-friendly qualitative palette (up to 20 haplotypes).
# Falls back to matplotlib tab20 for larger sets.
CB_PALETTE = [
    "#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE",
    "#AA3377", "#BBBBBB", "#000000", "#AA7733", "#44BB99",
    "#EEDD88", "#EE8866", "#DDDDDD", "#99DDFF", "#BB5566",
    "#004488", "#DDAA33", "#009988", "#CC3311", "#33BBEE",
]


# Base colors for each haplotype group prefix.
# Any haplotype whose string starts with a key gets that group's base color.
# Sub-types (A1-1, A1-2, …) are given lightened/darkened shades of the base.
HAP_GROUP_COLORS = {
    "A1": "#BB5566",
    "A2": "#004488",
}

def _adjust_lightness(hex_color, factor):
    """Return hex_color lightened (factor>1) or darkened (factor<1) in HLS space."""
    import colorsys
    hex_color = hex_color.lstrip("#")
    r, g, b = (int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    l = max(0.0, min(1.0, l * factor))
    r2, g2, b2 = colorsys.hls_to_rgb(h, l, s)
    return "#{:02X}{:02X}{:02X}".format(int(r2*255), int(g2*255), int(b2*255))

def build_palette(haplotypes):
    unique_haps = sorted(set(str(h) for h in haplotypes), key=natural_key)

    palette = {}
    # Group haplotypes by their matching prefix (longest match wins)
    groups = {}
    ungrouped = []
    for hap in unique_haps:
        matched = None
        for prefix in sorted(HAP_GROUP_COLORS, key=len, reverse=True):
            if hap.startswith(prefix):
                matched = prefix
                break
        if matched:
            groups.setdefault(matched, []).append(hap)
        else:
            ungrouped.append(hap)

    # Assign shaded variants within each group
    for prefix, members in groups.items():
        base = HAP_GROUP_COLORS[prefix]
        n = len(members)
        if n == 1:
            palette[members[0]] = base
        else:
            # Spread lightness: darker to lighter around the base (0.75 to 1.25)
            factors = [0.75 + 0.5 * (i / (n - 1)) for i in range(n)]
            for hap, factor in zip(members, factors):
                palette[hap] = _adjust_lightness(base, factor)

    # Fall back to CB_PALETTE for any unrecognised haplotypes
    cb_iter = iter(c for c in CB_PALETTE if c not in HAP_GROUP_COLORS.values())
    for hap in ungrouped:
        palette[hap] = next(cb_iter, "#AAAAAA")

    return palette


# ── Natural sort for strain IDs ───────────────────────────────────────────────

import re

def natural_key(s):
    """Sort '12A10' before '12A9' won't happen — numeric chunks are sorted as ints."""
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r"(\d+)", str(s))]


# ── Chromosome ordering ───────────────────────────────────────────────────────

def chrom_sort_key(c):
    s = str(c).lower().lstrip("chr")
    try:
        return (0, int(s))
    except ValueError:
        return (1, s)


# ── Sub-row offset layout ─────────────────────────────────────────────────────

def build_offsets(df, col_strain, col_hap):
    """
    For each strain, assign a vertical sub-row index to every distinct haplotype
    that appears for that strain.  Returns a dict keyed by (strain, hap) → float y.

    Haplotypes are sorted consistently across all strains so the same haplotype
    label always maps to the same relative sub-row position.

    Layout within a strain band of height 1.0:
        n sub-rows → offsets evenly spaced and centred around the strain's integer y.
        e.g. 2 haplotypes → -0.2 and +0.2
             3 haplotypes → -0.25, 0, +0.25
    The spacing shrinks for larger numbers of haplotypes so everything stays
    within the band.
    """
    # Global sorted haplotype list (consistent ordering across strains)
    all_haps = sorted(df[col_hap].unique(), key=natural_key)
    hap_rank = {h: i for i, h in enumerate(all_haps)}

    # Per-strain haplotype sets
    strain_haps = (
        df.groupby(col_strain)[col_hap]
        .apply(lambda s: sorted(s.unique(), key=lambda h: hap_rank[h]))
        .to_dict()
    )

    max_n = max(len(v) for v in strain_haps.values())
    # Max total spread: 0.7 of the row height (leaves 0.15 margin top & bottom)
    max_spread = 0.70

    offsets = {}
    strains_sorted = sorted(df[col_strain].unique(), key=natural_key)
    for strain_i, strain in enumerate(strains_sorted):
        haps = strain_haps[strain]
        n = len(haps)
        if n == 1:
            offsets[(strain, haps[0])] = strain_i
        else:
            spread = min(max_spread, 0.15 * (n - 1))   # tighten for large n
            step   = spread / (n - 1)
            start  = strain_i - spread / 2
            for j, hap in enumerate(haps):
                offsets[(strain, hap)] = start + j * step

    return offsets, strains_sorted


# ── Main plotting logic ───────────────────────────────────────────────────────

def plot_haplotypes(df, col_strain, col_start, col_stop, col_chrom, col_hap,
                    linewidth=6, output="haplotypes.png", dpi=150):

    # Cast numeric columns
    df[col_start] = pd.to_numeric(df[col_start])
    df[col_stop]  = pd.to_numeric(df[col_stop])
    df[col_hap]   = df[col_hap].astype(str)

    hap_colors = build_palette(df[col_hap])

    # Ordered axes
    chroms  = sorted(df[col_chrom].unique(), key=chrom_sort_key)
    offsets, strains = build_offsets(df, col_strain, col_hap)

    n_chroms  = len(chroms)
    n_strains = len(strains)
    strain_idx = {s: i for i, s in enumerate(strains)}

    # Figure sizing — subplot widths proportional to each chromosome's coordinate span
    chrom_spans = {
        c: df[df[col_chrom] == c][col_stop].max() - df[df[col_chrom] == c][col_start].min()
        for c in chroms
    }
    total_span = sum(chrom_spans.values())
    width_ratios = [chrom_spans[c] / total_span for c in chroms]

    # Fixed landscape letter size with margins (8.5x11 rotated, ~0.5in margins each side)
    PAGE_W, PAGE_H = 10.5, 7.5   # usable area after margins (landscape 11x8.5)
    fig_w, fig_h = PAGE_W, PAGE_H

    fig, axes = plt.subplots(
        1, n_chroms,
        figsize=(fig_w, fig_h),
        sharey=True,
        gridspec_kw={"wspace": 0.08, "width_ratios": width_ratios},
    )
    if n_chroms == 1:
        axes = [axes]

    fig.patch.set_facecolor("#F8F8F8")

    for ax, chrom in zip(axes, chroms):
        sub = df[df[col_chrom] == chrom]
        ax.set_facecolor("#FFFFFF")

        # Faint alternating band per strain for readability
        for i in range(n_strains):
            if i % 2 == 0:
                ax.axhspan(i - 0.5, i + 0.5, color="#F0F0F0", zorder=0)

        # Draw each segment at its sub-row y offset
        for _, row in sub.iterrows():
            key   = (row[col_strain], row[col_hap])
            y     = offsets.get(key, strain_idx[row[col_strain]])
            color = hap_colors[row[col_hap]]
            ax.plot(
                [row[col_start], row[col_stop]],
                [y, y],
                color=color,
                linewidth=linewidth,
                solid_capstyle="round",
                alpha=0.85,
            )

        # Chromosome axis formatting
        x_min = df[df[col_chrom] == chrom][col_start].min()
        x_max = df[df[col_chrom] == chrom][col_stop].max()
        # Snap limits to nearest Mbp outside the data range
        x_lo = np.floor(x_min / 1e6) * 1e6
        x_hi = np.ceil(x_max / 1e6) * 1e6
        ax.set_xlim(x_lo, x_hi)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1e6))
        ax.xaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, _: f"{int(x/1e6)}")
        )
        # Scale title/tick font to available row height
        row_px = fig_h / max(n_strains, 1)
        title_fs = max(5, min(9, row_px * 4))
        tick_fs  = max(4, min(7, row_px * 3))
        ax.set_title(str(chrom), fontsize=title_fs, fontweight="bold", pad=2, color="#333333")
        ax.set_xlabel("Mbp", fontsize=tick_fs, color="#555555", labelpad=1)
        ax.set_ylim(-0.8, n_strains - 0.2)
        ax.tick_params(axis="x", labelsize=tick_fs, rotation=45, pad=1)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.yaxis.set_visible(False)
        ax.grid(axis="x", linestyle=":", linewidth=0.4, color="#CCCCCC", zorder=1)

    # Y-axis strain labels centred at integer positions
    axes[0].yaxis.set_visible(True)
    axes[0].set_yticks(range(n_strains))
    axes[0].set_yticklabels(strains, fontsize=max(4, min(8, int(fig_h / n_strains * 5))))
    axes[0].tick_params(axis="y", length=0, pad=4)
    axes[0].spines["left"].set_visible(False)

    # Legend
    legend_handles = [
        mpatches.Patch(color=color, label=hap)
        for hap, color in sorted(hap_colors.items())
    ]
    fig.legend(
        handles=legend_handles,
        title="Haplotype",
        title_fontsize=7,
        fontsize=6,
        loc="lower center",
        ncol=min(len(hap_colors), 10),
        bbox_to_anchor=(0.5, -0.02),
        frameon=True,
        framealpha=0.9,
        edgecolor="#CCCCCC",
    )

    fig.suptitle("Single Cell Anomalies",
                 fontsize=10, fontweight="bold", y=1.005, color="#222222")

    plt.tight_layout(rect=[0, 0.07, 1, 1])
    fig.savefig(output, dpi=dpi, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    print(f"Saved: {output}  ({n_strains} strains × {n_chroms} chromosomes, "
          f"{len(hap_colors)} haplotypes)")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    print(f"Loaded {len(df):,} rows. Columns: {list(df.columns)}")

    col_strain = autodetect_col(df.columns, "strain",   args.col_strain)
    col_start  = autodetect_col(df.columns, "start",    args.col_start)
    col_stop   = autodetect_col(df.columns, "stop",     args.col_stop)
    col_chrom  = autodetect_col(df.columns, "chrom",    args.col_chrom)
    col_hap    = autodetect_col(df.columns, "haplotype",args.col_haplotype)

    print(f"Using columns → strain: '{col_strain}', start: '{col_start}', "
          f"stop: '{col_stop}', chrom: '{col_chrom}', haplotype: '{col_hap}'")

    plot_haplotypes(
        df, col_strain, col_start, col_stop, col_chrom, col_hap,
        linewidth=args.linewidth,
        output=args.output,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()