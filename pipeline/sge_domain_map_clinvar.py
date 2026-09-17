"""
SGE Domain Map Visualiser for SMC1A
====================================
For each input TSV file (one per targeton/exon), this script:
  1. Fetches SMC1A (P38531) domain annotations from UniProt REST API
  2. Fetches exon/CDS coordinates from Ensembl REST API to convert
     genomic positions -> protein amino acid positions
  3. Produces two plots per targeton:
       a) Lollipop plot  – variants as stems on a linear domain map
       b) Heatmap        – per-residue LFC coloured by domain
  4. Saves all plots to an output folder

Usage
-----
  python sge_domain_map.py --input *.tsv --outdir plots/
  ex: python code/sge_domain_map_clinvar_gaussian.py --input SMC1A_maveqc/all_screens/thesis/GMM_shrinkage/D4_ref/results/*.tsv --outdir SMC1A_maveqc/all_screens/thesis/GMM_shrinkage/D4_ref/plots/domain --clinvar SMC1A_maveqc/all_screens/clin_var/clinvar.vcf.gz


Requirements
------------
  pip install pandas matplotlib requests seaborn

Column assumptions (from your DESeq2/GMM results TSV):
  position                    : genomic position (chrX coordinates)
  consequence                 : variant consequence string
  pos_adj_log2FoldChange_raw  : LFC value to plot (pos_adj raw)
  anchor_tier                 : per-variant call from the GMM/shrinkage
                                 anchor pipelines ("enriched", "no impact",
                                 "weakly depleting", "strongly depleting").
                                 Used to flag variants as significant/
                                 impactful on the plots, with three levels
                                 of visual weight (strongly depleting /
                                 enriched = most prominent, weakly
                                 depleting = intermediate, no impact =
                                 faint background). Falls back to the
                                 older "GMM_status" column (depleted/no
                                 impact/enriched, pre-shrinkage pipelines),
                                 then to FDR-thresholding
                                 (pos_adj_fdr_raw < 0.05), if neither is
                                 present, with a warning, for backward
                                 compatibility with older result files.
"""

import argparse
import os
import sys
import requests
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import warnings
warnings.filterwarnings("ignore")

# ── Constants ────────────────────────────────────────────────────────────────

UNIPROT_ID   = "Q14683"          # SMC1A human
TRANSCRIPT   = "ENST00000322213" # canonical SMC1A transcript
LFC_COL      = "pos_adj_log2FoldChange_raw"
FDR_COL      = "pos_adj_fdr_raw"
FDR_THRESH   = 0.05
STATUS_COL   = "anchor_tier"      # preferred significance source (GMM/shrinkage anchor outputs)
LEGACY_STATUS_COL = "GMM_status"  # older pre-shrinkage pipelines, used if STATUS_COL absent
# Anything in this set is treated as significant/impactful; "no impact" is not.
# Includes both the new 4-level anchor_tier vocabulary ("weakly depleting",
# "strongly depleting", "enriched") and the older 3-level GMM_status
# vocabulary ("depleted", "enriched") for backward compatibility.
IMPACT_STATUSES = {"depleted", "enriched", "weakly depleting", "strongly depleting"}

# Per-tier marker styling for the lollipop plot: strongly depleting/enriched
# are the most visually prominent, weakly depleting is intermediate, and
# no impact (or anything unrecognised) fades into the background. Falls
# back to DEFAULT_TIER_STYLE for any status string not listed here.
TIER_STYLE = {
    "strongly depleting": {"s": 20, "alpha": 1.0,  "edge": "white", "lw": 0.5},
    "enriched":           {"s": 20, "alpha": 1.0,  "edge": "white", "lw": 0.5},
    "weakly depleting":   {"s": 13, "alpha": 0.75, "edge": "white", "lw": 0.3},
    "no impact":          {"s": 8,  "alpha": 0.5,  "edge": "none",  "lw": 0.0},
    "depleted":           {"s": 20, "alpha": 1.0,  "edge": "white", "lw": 0.5},  # legacy GMM_status
}
DEFAULT_TIER_STYLE = {"s": 8, "alpha": 0.5, "edge": "none", "lw": 0.0}
PROTEIN_LEN  = 1233              # SMC1A aa length (fallback if fetch fails)

# Consequence colour palette — keys are substrings matched case-insensitively
CONSEQUENCE_COLOURS = {
    "Missense_Variant":                    "#e74c3c",
    "Synonymous_Variant":                  "#3498db",
    "Nonsense_Variant":                    "#8e44ad",
    "LOF":                                 "#8e44ad",
    "Intronic_Variant":                    "#95a5a6",
    "Splice_Polypyrimidine_Tract_Variant": "#e67e22",
    "Splice_Variant":                      "#d35400",
    "Inframe_Deletion":                    "#1abc9c",
    "Inframe_Insertion":                   "#16a085",
    "Others":                              "#bdc3c7",
}
DEFAULT_COLOUR = "#bdc3c7"

# Domain colour palette (cycles if more domains than colours)
DOMAIN_PALETTE = [
    "#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f",
    "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac",
]

# ── UniProt fetch ─────────────────────────────────────────────────────────────

def fetch_uniprot_features(uniprot_id: str) -> tuple[list[dict], int]:
    """
    Returns (features, protein_length).
    features is a list of dicts with keys: type, description, start, end
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    print(f"  Fetching UniProt features for {uniprot_id}...")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    length = data.get("sequence", {}).get("length", PROTEIN_LEN)

    keep_types = {"Domain", "Region", "Coiled coil", "Motif", "Zinc finger",
                  "Repeat", "Compositional bias"}
    features = []
    for f in data.get("features", []):
        ftype = f.get("type", "")
        if ftype in keep_types:
            try:
                start = int(f["location"]["start"]["value"])
                end   = int(f["location"]["end"]["value"])
                desc  = f.get("description") or f.get("type", "Unknown")
                features.append({"type": ftype, "description": desc,
                                  "start": start, "end": end})
            except (KeyError, TypeError):
                continue

    print(f"  Found {len(features)} annotated features, protein length = {length} aa")
    return features, length


# ── AA -> Genomic coordinate mapping (hardcoded SMC1A / ENST00000322213) ──────
#
# SMC1A is on chrX, NEGATIVE strand.
# CDS exons in genomic order (low -> high), but transcribed right-to-left.
# Source: Ensembl ENST00000322213.9 / UCSC hg38
# Each tuple: (genomic_start, genomic_end) — both inclusive, 1-based.
#
SMC1A_CDS_EXONS_HG38 = [
    # CDS-only exon coordinates, trimmed to chrX:53380103-53422600
    # Source: Ensembl ENST00000322213.9 (MANE Select), hg38
    # Gene is on minus strand; exons listed high->low for readability
    (53422492, 53422600),  # exon 1,  109 nt
    (53414981, 53415169),  # exon 2,  189 nt
    (53414758, 53414870),  # exon 3,  113 nt
    (53413232, 53413435),  # exon 4,  204 nt
    (53412900, 53413138),  # exon 5,  239 nt
    (53411995, 53412253),  # exon 6,  259 nt
    (53411761, 53411902),  # exon 7,  142 nt
    (53409421, 53409503),  # exon 8,   83 nt
    (53409062, 53409269),  # exon 9,  208 nt
    (53405771, 53405956),  # exon 10, 186 nt
    (53405493, 53405672),  # exon 11, 180 nt
    (53405245, 53405391),  # exon 12, 147 nt
    (53405012, 53405149),  # exon 13, 138 nt
    (53403777, 53403893),  # exon 14, 117 nt
    (53403566, 53403672),  # exon 15, 107 nt
    (53399589, 53399730),  # exon 16, 142 nt
    (53396472, 53396617),  # exon 17, 146 nt
    (53396227, 53396380),  # exon 18, 154 nt
    (53394778, 53394888),  # exon 19, 111 nt
    (53383097, 53383253),  # exon 20, 157 nt
    (53382506, 53382660),  # exon 21, 155 nt
    (53382232, 53382383),  # exon 22, 152 nt
    (53381018, 53381087),  # exon 23,  70 nt
    (53380620, 53380730),  # exon 24, 111 nt
    (53380103, 53380186),  # exon 25,  84 nt
]


def build_aa_genomic_map() -> pd.DataFrame:
    """
    Build a DataFrame mapping aa_pos (1-based) -> genomic_pos for SMC1A.
    SMC1A is on the minus strand so we walk exons from high genomic coord
    to low genomic coord, assigning CDS nucleotide positions in order.
    """
    # Minus strand: walk exons from largest genomic coord downward
    exons_sorted = sorted(SMC1A_CDS_EXONS_HG38, key=lambda e: -e[0])
    rows = []
    cds_pos = 1
    for (e_start, e_end) in exons_sorted:
        for g in range(e_end, e_start - 1, -1):  # count down (minus strand)
            aa = (cds_pos - 1) // 3 + 1
            rows.append({"genomic_pos": g, "aa_pos": aa})
            cds_pos += 1
    df = pd.DataFrame(rows)
    print(f"  Built AA->genomic map: {len(df)} nt positions, "
          f"{df['aa_pos'].max()} amino acids")
    return df


def aa_to_genomic_blocks(aa_start: int, aa_end: int,
                          aa_map: pd.DataFrame) -> list[tuple[int, int]]:
    """
    Convert an AA range to a list of contiguous genomic blocks.
    Returns [(g_start, g_end), ...] — one per exon the domain spans.
    """
    sub = aa_map[(aa_map["aa_pos"] >= aa_start) &
                 (aa_map["aa_pos"] <= aa_end)]["genomic_pos"].sort_values()
    if sub.empty:
        return []
    # Split into contiguous blocks
    blocks = []
    block_start = prev = sub.iloc[0]
    for g in sub.iloc[1:]:
        if g != prev - 1 and g != prev + 1:  # gap detected
            blocks.append((min(block_start, prev), max(block_start, prev)))
            block_start = g
        prev = g
    blocks.append((min(block_start, prev), max(block_start, prev)))
    return blocks


def fetch_domain_genomic_coords(features: list[dict],
                                 transcript_id: str) -> list[dict]:
    """
    Map each UniProt domain's AA coordinates to genomic blocks using the
    hardcoded SMC1A exon structure. No network calls needed.
    """
    print(f"  Computing domain genomic coordinates from hardcoded "
          f"SMC1A exon structure (no Ensembl needed)...")
    aa_map = build_aa_genomic_map()
    updated = []
    for f in features:
        blocks = aa_to_genomic_blocks(f["start"], f["end"], aa_map)
        f = dict(f, genomic_blocks=blocks)
        print(f"    {f['description']} (AA {f['start']}-{f['end']}): "
              f"{len(blocks)} genomic block(s) -> "
              f"{[f'chrX:{b[0]}-{b[1]}' for b in blocks]}")
        updated.append(f)
    return updated


# ── Data loading ──────────────────────────────────────────────────────────────

def load_targeton(filepath: str, status_col: str = STATUS_COL) -> pd.DataFrame:
    """Load one targeton TSV.

    Significance/impact flagging prefers `status_col` (default "anchor_tier"):
    rows are treated as significant if their value is in IMPACT_STATUSES
    ("weakly depleting", "strongly depleting", "enriched" — or the legacy
    "depleted"/"enriched" from older GMM_status files); anything else (e.g.
    "no impact") is not.

    If `status_col` isn't present, falls back to LEGACY_STATUS_COL
    ("GMM_status", for older pre-shrinkage pipelines), then to FDR
    thresholding (pos_adj_fdr_raw < FDR_THRESH) for backward compatibility
    with pre-GMM result files, printing a warning at each fallback step so
    this doesn't happen silently.
    """
    df = pd.read_csv(filepath, sep="\t")

    required = ["position", "consequence", LFC_COL]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing column '{col}' in {filepath}. "
                             f"Available: {list(df.columns)}")

    optional_cols = []
    rename_map = {"position": "genomic_pos", LFC_COL: "lfc"}

    if status_col in df.columns:
        optional_cols.append(status_col)
        rename_map[status_col] = "status"
    elif LEGACY_STATUS_COL in df.columns:
        print(f"    NOTE: status column '{status_col}' not found in "
              f"{os.path.basename(filepath)} — falling back to legacy "
              f"'{LEGACY_STATUS_COL}' column for significance flagging.")
        optional_cols.append(LEGACY_STATUS_COL)
        rename_map[LEGACY_STATUS_COL] = "status"
    else:
        print(f"    WARNING: neither '{status_col}' nor '{LEGACY_STATUS_COL}' "
              f"found in {os.path.basename(filepath)} — falling back to FDR "
              f"thresholding for significance flagging.")
        if FDR_COL in df.columns:
            optional_cols.append(FDR_COL)
            rename_map[FDR_COL] = "fdr"

    df = df[required + optional_cols].copy()
    df.rename(columns=rename_map, inplace=True)

    df["aa_pos"] = np.nan

    df["lfc"] = pd.to_numeric(df["lfc"], errors="coerce")
    df.dropna(subset=["lfc"], inplace=True)

    if "status" in df.columns:
        status_norm = df["status"].astype(str).str.strip().str.lower()
        counts = status_norm.value_counts()
        counts_str = ", ".join(f"{k}: {v}" for k, v in counts.items())
        print(f"    Status counts — {counts_str}")

    return df


# ── Plotting helpers ──────────────────────────────────────────────────────────

def assign_domain_colours(features: list[dict]) -> dict:
    """Return {description: colour} for each unique domain."""
    seen = {}
    colour_idx = 0
    for f in features:
        desc = f["description"]
        if desc not in seen:
            seen[desc] = DOMAIN_PALETTE[colour_idx % len(DOMAIN_PALETTE)]
            colour_idx += 1
    return seen


def get_variant_colour(consequence: str) -> str:
    csq = str(consequence).strip()
    # Exact match first
    if csq in CONSEQUENCE_COLOURS:
        return CONSEQUENCE_COLOURS[csq]
    # Substring match (case-insensitive) as fallback
    csq_lower = csq.lower()
    for key, col in CONSEQUENCE_COLOURS.items():
        if key.lower() in csq_lower:
            return col
    return DEFAULT_COLOUR


def draw_domain_bar(ax, features, domain_colours, protein_len,
                    x_start=1, y=0, height=0.35, label=True,
                    use_genomic=False):
    """Draw the protein backbone and domain rectangles on ax.
    If use_genomic=True, uses f['genomic_blocks'] for coordinates."""
    # Backbone
    ax.add_patch(mpatches.FancyBboxPatch(
        (x_start, y - height / 2), protein_len, height,
        boxstyle="round,pad=0", linewidth=0.5,
        edgecolor="#555", facecolor="#e8e8e8", zorder=1))

    for f in features:
        col = domain_colours.get(f["description"], "#cccccc")
        if use_genomic:
            blocks = f.get("genomic_blocks", [])
        else:
            blocks = [(f["start"], f["end"])]

        for (b_start, b_end) in blocks:
            w = max(b_end - b_start, 1)
            ax.add_patch(mpatches.FancyBboxPatch(
                (b_start, y - height / 2), w, height,
                boxstyle="round,pad=0", linewidth=0.3,
                edgecolor="white", facecolor=col, alpha=0.85, zorder=2))
            if label and w > protein_len * 0.03:
                ax.text(b_start + w / 2, y, f["description"],
                        ha="center", va="center", fontsize=5.5,
                        fontweight="bold", color="white", zorder=3,
                        clip_on=True)


# ── ClinVar fetch ─────────────────────────────────────────────────────────────

CLINVAR_COLOURS = {
    "Pathogenic":           "#c0392b",
    "Likely pathogenic":    "#e67e22",
}

def fetch_clinvar_variants(vcf_path: str) -> pd.DataFrame:
    """
    Read Pathogenic and Likely pathogenic ClinVar variants for SMC1A
    from a local ClinVar VCF file (GRCh38, chrX).
    Returns DataFrame with columns: genomic_pos, significance.
    """
    import gzip
    REGION_START = 53374149
    REGION_END   = 53422654

    print(f"  Reading ClinVar VCF from {vcf_path}...")
    rows = []
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            # Filter to chrX / X
            chrom = parts[0]
            if chrom not in ("X", "chrX"):
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            if not (REGION_START <= pos <= REGION_END):
                continue
            info = parts[7]
            clnsig = ""
            for field in info.split(";"):
                if field.startswith("CLNSIG="):
                    clnsig = field[7:].replace("_", " ").replace("%2C", ",")
                    break
            sig = None
            clnsig_lower = clnsig.lower()
            if "likely pathogenic" in clnsig_lower:
                sig = "Likely pathogenic"
            elif "pathogenic" in clnsig_lower:
                sig = "Pathogenic"
            if sig is None:
                continue
            rows.append({"genomic_pos": pos, "significance": sig})
    df = pd.DataFrame(rows)
    print(f"  Retained {len(df)} P/LP variants in SMC1A region")
    return df


# ── Main plots ────────────────────────────────────────────────────────────────

def plot_lollipop(df: pd.DataFrame, features: list[dict],
                  domain_colours: dict, protein_len: int,
                  title: str, outpath: str,
                  clinvar_df: pd.DataFrame = None):
    """Lollipop plot with domain bar and optional ClinVar triangle track."""
    use_aa = df["aa_pos"].notna().any()
    x_col  = "aa_pos" if use_aa else "genomic_pos"
    x_label = "Amino acid position" if use_aa else "Genomic position (chrX)"

    fig, ax = plt.subplots(figsize=(14, 5))
    fig.patch.set_facecolor("white")

    # Zero line
    ax.axhline(0, color="#aaaaaa", lw=0.8, zorder=0)

    # Domain bar position
    lfc_range = df["lfc"].abs().max()
    bar_y = -(lfc_range + 0.6)
    bar_width = protein_len if use_aa else (df["genomic_pos"].max() - df["genomic_pos"].min())
    bar_start = 1 if use_aa else df["genomic_pos"].min()
    draw_domain_bar(ax, features, domain_colours,
                    bar_width, x_start=bar_start,
                    y=bar_y, height=0.35,
                    use_genomic=not use_aa)

    # ClinVar triangle track just below domain bar
    clinvar_y = bar_y - 0.45
    if clinvar_df is not None and not clinvar_df.empty:
        # Filter to variants within this targeton's x range
        if use_aa:
            x_lo_data, x_hi_data = 1, protein_len
        else:
            pad = (df["genomic_pos"].max() - df["genomic_pos"].min()) * 0.05 + 5
            x_lo_data = df["genomic_pos"].min() - pad
            x_hi_data = df["genomic_pos"].max() + pad
        cv_visible = clinvar_df[
            (clinvar_df["genomic_pos"] >= x_lo_data) &
            (clinvar_df["genomic_pos"] <= x_hi_data)
        ]
        for _, cv in cv_visible.iterrows():
            col = CLINVAR_COLOURS.get(cv["significance"], "#888888")
            ax.scatter(cv["genomic_pos"], clinvar_y,
                       marker="v", s=30, color=col,
                       zorder=5, linewidths=0, clip_on=True)

    # Lollipops
    for _, row in df.iterrows():
        x   = row[x_col]
        lfc = row["lfc"]
        if pd.isna(x):
            continue
        if "status" in df.columns:
            status_val = str(row.get("status", "")).strip().lower()
            style = TIER_STYLE.get(status_val, DEFAULT_TIER_STYLE)
        elif "fdr" in df.columns:
            sig = row.get("fdr", 1.0) < FDR_THRESH
            style = TIER_STYLE["strongly depleting"] if sig else DEFAULT_TIER_STYLE
        else:
            style = TIER_STYLE["strongly depleting"]
        col = get_variant_colour(row["consequence"])
        ax.vlines(x, 0, lfc, colors=col, linewidth=0.6, alpha=0.7, zorder=3)
        ax.scatter(x, lfc, color=col,
                   s=style["s"],
                   zorder=4,
                   edgecolors=style["edge"],
                   linewidths=style["lw"],
                   alpha=style["alpha"])

    # Domain legend
    domain_handles = [
        mpatches.Patch(facecolor=col, label=desc, alpha=0.85)
        for desc, col in domain_colours.items()
    ]
    # Consequence legend
    present_csq = df["consequence"].unique()
    csq_handles = [
        mpatches.Patch(facecolor=get_variant_colour(c), label=c)
        for c in present_csq
    ]
    # ClinVar legend
    cv_handles = [
        plt.scatter([], [], marker="v", color=col, s=30, label=sig)
        for sig, col in CLINVAR_COLOURS.items()
    ]

    leg1 = ax.legend(handles=domain_handles, title="Domains",
                     loc="upper left", fontsize=6, title_fontsize=7,
                     framealpha=0.8, ncol=2)
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=csq_handles, title="Consequence",
                     loc="upper right", fontsize=6, title_fontsize=7,
                     framealpha=0.8)
    ax.add_artist(leg2)
    leg3 = None
    if clinvar_df is not None and not clinvar_df.empty:
        leg3 = ax.legend(handles=cv_handles, title="ClinVar",
                         loc="lower right", fontsize=6, title_fontsize=7,
                         framealpha=0.8)

    sig_note = ("Largest, white-edged points = strongly depleting/enriched; "
                 "medium, white-edged points = weakly depleting; "
                 "small, faded points = no impact"
                 if "status" in df.columns else
                 "Larger, white-edged points = FDR-significant")
    sig_note_artist = ax.text(0.5, -0.16, sig_note, transform=ax.transAxes,
                              ha="center", va="top", fontsize=6.5, color="#555555")

    if use_aa:
        x_lo, x_hi = 1, protein_len
    else:
        pad = (df["genomic_pos"].max() - df["genomic_pos"].min()) * 0.05 + 5
        x_lo = df["genomic_pos"].min() - pad
        x_hi = df["genomic_pos"].max() + pad
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(clinvar_y - 0.3, lfc_range + 0.5)
    ax.set_xlabel(x_label, fontsize=9)
    ax.set_ylabel("LFC (pos_adj raw)", fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)

    plt.tight_layout()
    fig.subplots_adjust(bottom=0.22)  # reserve room for sig_note text below axes
    extra_artists = [a for a in (leg1, leg2, leg3, sig_note_artist) if a is not None]
    plt.savefig(outpath, dpi=180, bbox_inches="tight",
                bbox_extra_artists=extra_artists)
    plt.close()
    print(f"    Saved: {outpath}")


def plot_heatmap(df: pd.DataFrame, features: list[dict],
                 domain_colours: dict, protein_len: int,
                 title: str, outpath: str):
    """1-D heatmap of LFC across positions with domain boundaries."""
    use_aa = df["aa_pos"].notna().any()
    x_col  = "aa_pos" if use_aa else "genomic_pos"
    x_label = "Amino acid position" if use_aa else "Genomic position (chrX)"

    # Build per-position mean LFC array
    x_max = int(protein_len) if use_aa else int(df["genomic_pos"].max())
    x_min = 1 if use_aa else int(df["genomic_pos"].min())
    lfc_arr = np.full(x_max - x_min + 1, np.nan)

    grp = df.dropna(subset=[x_col]).groupby(x_col)["lfc"].mean()
    for pos, val in grp.items():
        idx = int(pos) - x_min
        if 0 <= idx < len(lfc_arr):
            lfc_arr[idx] = val

    # Figure height needs to flex with how many legend rows the domains need —
    # a fixed height here is what was causing the legend (and the row below it)
    # to sometimes run off the bottom of the saved image when a targeton has
    # many overlapping domains.
    n_domains  = max(len(domain_colours), 1)
    leg_ncol   = min(4, n_domains)
    leg_rows   = int(np.ceil(n_domains / leg_ncol))
    leg_ratio  = 0.35 * leg_rows
    fig_height = 3 + 0.5 + leg_ratio

    fig, axes = plt.subplots(
        3, 1, figsize=(16, fig_height),
        gridspec_kw={"height_ratios": [3, 0.5, leg_ratio]},
    )
    fig.patch.set_facecolor("white")

    # ── Heatmap row ──
    ax_heat = axes[0]
    vmax = np.nanmax(np.abs(lfc_arr)) if not np.all(np.isnan(lfc_arr)) else 1
    cmap = plt.get_cmap("RdBu_r")
    norm = Normalize(vmin=-vmax, vmax=vmax)

    ax_heat.imshow(
        lfc_arr[np.newaxis, :],
        aspect="auto", cmap=cmap, norm=norm,
        extent=[x_min - 0.5, x_max + 0.5, -0.5, 0.5],
        interpolation="nearest",
    )
    # Domain boundary lines
    for f in features:
        for boundary in (f["start"], f["end"]):
            if x_min <= boundary <= x_max:
                ax_heat.axvline(boundary, color="white", lw=0.8, alpha=0.7)

    ax_heat.set_yticks([])
    ax_heat.set_xlim(x_min, x_max)
    ax_heat.set_ylabel("LFC", fontsize=8)
    ax_heat.set_title(title, fontsize=10, fontweight="bold")
    ax_heat.tick_params(axis="x", labelsize=7)

    # Colorbar
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax_heat, orientation="vertical",
                        fraction=0.015, pad=0.01)
    cbar.set_label("LFC", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    # ── Domain bar row ──
    ax_dom = axes[1]
    ax_dom.set_xlim(x_min, x_max)
    ax_dom.set_ylim(-0.5, 0.5)
    ax_dom.set_yticks([])
    draw_domain_bar(ax_dom, features, domain_colours,
                    x_max - x_min, x_start=x_min,
                    y=0, height=0.8, label=True,
                    use_genomic=not use_aa)
    ax_dom.spines[:].set_visible(False)
    ax_dom.tick_params(axis="x", labelsize=7)
    ax_dom.set_xlabel(x_label, fontsize=8)

    # ── Domain colour legend row ──
    ax_leg = axes[2]
    ax_leg.axis("off")
    handles = [mpatches.Patch(facecolor=col, label=desc, alpha=0.85)
               for desc, col in domain_colours.items()]
    dom_legend = ax_leg.legend(handles=handles, loc="center", ncol=leg_ncol,
                               fontsize=6, frameon=False, title="Domains",
                               title_fontsize=7)

    plt.tight_layout(h_pad=0.3)
    plt.savefig(outpath, dpi=180, bbox_inches="tight",
                bbox_extra_artists=[dom_legend])
    plt.close()
    print(f"    Saved: {outpath}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="SGE domain map visualiser for SMC1A targetons")
    parser.add_argument("--input", nargs="+", required=True,
                        help="One or more DESeq2 result TSV files (one per targeton)")
    parser.add_argument("--outdir", default="sge_domain_plots",
                        help="Output directory for plots (default: sge_domain_plots/)")
    parser.add_argument("--uniprot", default=UNIPROT_ID,
                        help=f"UniProt ID (default: {UNIPROT_ID})")
    parser.add_argument("--transcript", default=TRANSCRIPT,
                        help=f"Ensembl transcript ID (default: {TRANSCRIPT})")
    parser.add_argument("--clinvar", default=None,
                        help="Path to local ClinVar VCF file (clinvar_chrX.vcf.gz)")
    parser.add_argument("--no-aa-map", action="store_true",
                        help="Skip genomic->AA mapping and plot in genomic coords")
    parser.add_argument("--status_col", default=STATUS_COL,
                        help=f"Column used to flag significant/impactful variants "
                             f"(default: {STATUS_COL}). Values in {sorted(IMPACT_STATUSES)} "
                             f"are treated as significant, with visual weight scaled by "
                             f"tier (strongly depleting/enriched > weakly depleting > "
                             f"no impact); falls back to '{LEGACY_STATUS_COL}' and then to "
                             f"FDR thresholding if this column is absent.")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ── Fetch annotations once ──
    print("\n[1/3] Fetching domain annotations from UniProt...")
    try:
        features, protein_len = fetch_uniprot_features(args.uniprot)
    except Exception as e:
        print(f"  ERROR fetching UniProt: {e}")
        features, protein_len = [], PROTEIN_LEN

    domain_colours = assign_domain_colours(features)

    print("\n[2/3] Mapping domain coordinates to genomic positions...")
    if not args.no_aa_map and features:
        try:
            features = fetch_domain_genomic_coords(features, args.transcript)
            print("  Domain genomic mapping complete.")
        except Exception as e:
            print(f"  WARNING: Domain mapping failed ({e}). "
                  "Domain bar will be unavailable.")

    print("\n[3/4] Reading ClinVar P/LP variants...")
    clinvar_df = pd.DataFrame()
    if args.clinvar:
        try:
            clinvar_df = fetch_clinvar_variants(args.clinvar)
        except Exception as e:
            print(f"  WARNING: ClinVar fetch failed ({e}). "
                  "ClinVar track will be omitted from plots.")
    else:
        print("  No --clinvar file provided, skipping ClinVar track.")

    # ── Process each targeton file ──
    print(f"\n[4/4] Processing {len(args.input)} targeton file(s)...\n")
    for filepath in args.input:
        name = os.path.splitext(os.path.basename(filepath))[0]
        print(f"  → {name}")

        try:
            df = load_targeton(filepath, status_col=args.status_col)
        except Exception as e:
            print(f"    ERROR loading file: {e}")
            continue

        if df.empty:
            print("    WARNING: No usable rows after filtering — skipping.")
            continue

        title_base = name.replace("_", " ")

        # Lollipop
        lollipop_path = os.path.join(args.outdir, f"{name}_lollipop.png")
        plot_lollipop(df, features, domain_colours, protein_len,
                      title=f"{title_base} — Lollipop",
                      outpath=lollipop_path,
                      clinvar_df=clinvar_df)

        # Heatmap
        heatmap_path = os.path.join(args.outdir, f"{name}_heatmap.png")
        plot_heatmap(df, features, domain_colours, protein_len,
                     title=f"{title_base} — Heatmap",
                     outpath=heatmap_path)

    print(f"\nDone! Plots saved to: {os.path.abspath(args.outdir)}/\n")

if __name__ == "__main__":
    main()
