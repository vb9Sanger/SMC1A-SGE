#!/usr/bin/env python3
"""
sge_clinvar_intersect_plot.py  (V2)

For each targeton, intersects SGE DESeq2 results (with GMM-derived depletion
calls) with ClinVar variants and produces:

  1. A "full" lollipop plot (all ClinVar variants for the targeton), showing:
       - Shape     : depletion status (anchor_tier: v strongly/weakly depleting |
                     o no impact | ^ enriched). Weakly depleting points are drawn
                     as a smaller downward triangle than strongly depleting/
                     enriched ones (same full opacity), so tier strength reads
                     at a glance without adding a 4th shape.
       - Colour    : disease condition (CdLS only | DEE only | Both | Other SMC1A |
                     Unclassified | No Condition Provided)
       - Border    : ClinVar classification (P/LP thick black border | other thin grey)

  2. A "condition-split" plot restricted to variants classified as DEE-only or
     CdLS-only (the two conditions of primary interest). Since condition is
     fixed per plot/panel here, points are coloured by variant consequence
     type instead (Missense_Variant, LOF, Synonymous_Variant, etc. - read
     from 'Summary_Plot', falling back to 'Summary_Consequence', falling
     back to 'Consequence'):
       - If a targeton has variants in only ONE of these categories, this is
         the same style of lollipop plot as (1), limited to that category.
       - If a targeton has variants in BOTH categories, this is a two-panel
         figure with a shared x-axis (genomic position): DEE-only variants in
         the top panel, CdLS-only variants in the bottom panel, so the two
         conditions can be visually compared side by side.

Depletion status is read from the 'anchor_tier' column if present (4 levels:
"strongly depleting", "weakly depleting", "no impact", "enriched" — the
output of the GMM/shrinkage anchor pipelines). Falls back to the older
'GMM_status' column (3-level: depleted/no impact/enriched) if 'anchor_tier'
is absent, then to 'stat_pos_raw' for older input files without either
anchor column.

Colour encoding (condition — used in the full plot and its own legend):
    DEE only        -> orange
    CdLS only       -> purple
    Both            -> green
    Other SMC1A     -> blue
    Unclassified    -> grey
    No Condition Provided  -> light grey (not plotted, background category only)

Colour encoding (variant consequence — used in the condition-split plots):
    Missense_Variant                    -> green
    LOF                                 -> red
    Synonymous_Variant                  -> dark green
    Inframe_Deletion                    -> grey
    Intronic_Variant                    -> blue
    Splice_Variant                      -> pink/magenta
    Splice_Polypyrimidine_Tract_Variant -> dark slate blue
    (anything else)                     -> light grey ("Other / unclassified")

Usage:
    python sge_clinvar_intersect_plot.py \\
        --deseq2_dir        all_deseq2_results.tsv \\
        --meta_dir      /path/to/meta_consequence/files \\
        --vcf_dir       /path/to/per_targeton_vcfs/ \\
        --regions       targeton_regions.tsv \\
        --outdir        clinvar_intersect_plots/

Arguments:
    --deseq2_dir    DESeq2/GMM-anchor results TSVs (all targetons, has 'position',
                    'oligo_name', 'pos_adj_log2FoldChange_raw', and one of
                    'anchor_tier' (preferred), 'GMM_status', or 'stat_pos_raw')
    --meta_dir      Directory containing *_meta_consequences.tsv files
                    (used to join ref/alt alleles onto DESeq2 results)
    --vcf_dir       Directory containing per-targeton ClinVar VCFs
                    (<Targeton_ID>_clinvar.vcf.gz)
    --regions       TSV produced by extract_targeton_regions.py
                    (columns: Targeton_ID, chrom, start, end)
    --outdir        Output directory for plots
    --de_novo_only  Restrict PLOTS to ClinVar variants with the de-novo bit
                    set in ORIGIN (does not affect the TSVs — see below)

Outputs per targeton (in --outdir):
    <tid>_clinvar_annotated.tsv     Original DESeq2 columns + meta annotation
                                     + ClinVar match info, incl. clinvar_id
                                     (ClinVar Allele ID from VCF ID column -
                                     paste into ClinVar's search box to look
                                     up a specific variant), origin_raw (raw
                                     ORIGIN bitmask string), origin_flags
                                     (decoded terms, e.g. "germline;de novo"),
                                     and is_de_novo (bool). These origin
                                     columns are always populated regardless
                                     of --de_novo_only, so you can filter on
                                     your own terms after the fact.
    <tid>_clinvar_intersect.png     Full lollipop plot (all conditions).
    <tid>_DEE_only.png              (if only DEE-only variants present)
    <tid>_CdLS_only.png             (if only CdLS-only variants present)
    <tid>_DEE_CdLS_stacked.png      (if both DEE-only and CdLS-only present)
    clinvar_variants_summary.tsv    ONE combined TSV (not per-targeton) listing
                                     every ClinVar-matched variant across all
                                     targetons, with Targeton_ID and exon
                                     (from TARGETON_EXON) prepended, so e.g.
                                     the proportion of CdLS-only variants that
                                     are depleted vs not can be checked in one
                                     place by filtering condition/anchor_tier.

Requirements:
    pip install pandas matplotlib numpy
    bcftools on PATH (for VCF reading)
"""

import argparse
import gzip
import os
import re
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# ── Disease condition classification ─────────────────────────────────────────

# Substrings to match in CLNDN field (case-insensitive)
CDLS_TERMS = [
    "congenital_muscular_hypertrophy-cerebral_syndrome",
    "cornelia_de_lange",
    "de_lange_syndrome",
]
DEE_TERMS = [
    "developmental_and_epileptic_encephalopathy",
    "epileptic_encephalopathy",
    "seizure",
]
# Terms that suggest SMC1A-specific but condition-unspecified
OTHER_SMC1A_TERMS = [
    "smc1a-related",
    "x-linked_complex_neurodevelopmental",
    "cohesinopathy",
]
# Terms to treat as uninformative
UNINFORMATIVE_TERMS = [
    "not_provided", "not_specified", "inborn_genetic_diseases",
    "see_cases", "history_of_neurodevelopmental_disorder",
]

# Visual encoding ─────────────────────────────────────────────────────────────

CONDITION_COLOURS = {
    "DEE only":          "#f5a623",   # orange
    "CdLS only":         "#8b4fcf",   # purple
    "Both":              "#4caf50",   # green
    "Other SMC1A":       "#4a90d9",   # blue
    "Unclassified":      "#9b9b9b",   # mid grey
    "No Condition Provided":    "#d4d4d4",   # light grey
}

# The two conditions of primary interest for the condition-split plot
SPLIT_CONDITIONS = ["DEE only", "CdLS only"]

# In the condition-split plots (DEE-only / CdLS-only), condition is already
# fixed per panel, so points there are coloured by variant consequence type
# instead. Read from 'Summary_Plot' (falling back to 'Summary_Consequence',
# then 'Consequence'). Note 'Summary_Plot' is preferred over
# 'Summary_Consequence' because the latter labels stop-gained variants as
# 'Nonsense_Variant' rather than 'LOF'.
CONSEQUENCE_COLOURS = {
    "Missense_Variant":                    "#5aa95a",   # green
    "LOF":                                 "#e05c4e",   # red
    "Synonymous_Variant":                  "#3f6b3f",   # dark green
    "Inframe_Deletion":                    "#b0b0b0",   # grey
    "Intronic_Variant":                    "#5b7fe0",   # blue
    "Splice_Variant":                      "#c66bb0",   # pink/magenta
    "Splice_Polypyrimidine_Tract_Variant": "#5b6b8c",   # dark slate blue
    "Others":                              "#dddddd",   # light grey (catch-all bucket in Summary_Plot)
}
CONSEQUENCE_OTHER_COLOUR = "#dddddd"  # any consequence label not in the map above (unexpected labels)

# Exon each targeton corresponds to (used for the cross-targeton ClinVar
# variant summary TSV).
TARGETON_EXON = {
    "SBPS": 1,  "APDY": 2,  "NBXL": 3,  "RSDC": 4,  "AXFG": 5,  "EXTP": 5,
    "SWAZ": 6,  "CQEJ": 6,  "IHMH": 7,  "KJQT": 8,  "YTUB": 9,  "RYFS": 10,
    "PSYW": 11, "NQBU": 12, "NUIB": 13, "SUHO": 14, "HDOI": 15, "MGTQ": 16,
    "RVMB": 17, "ELWX": 18, "OJTT": 19, "ZANL": 20, "RBMZ": 21, "CNJV": 22,
    "NLVE": 23, "FKKY": 24, "FBSZ": 25,
}

# Column set (and order) for the cross-targeton ClinVar variant summary TSV,
# in addition to the leading Targeton_ID and exon columns.
SUMMARY_COLUMNS = [
    "sequence", "oligo_name", "position", "pos_fit", "pos_fit_se",
    "pos_total_se_raw", "pos_total_se_shrunk", "pos_adj_log2FoldChange_raw",
    "pos_adj_score_raw", "pos_adj_pval_raw", "pos_adj_fdr_raw", "stat_pos_raw",
    # GMM/shrinkage anchor pipeline columns (replaces the older GMM_fit_set /
    # GMM_cluster / GMM_label / GMM_prob_cluster* / GMM_status* columns from
    # the pre-anchor classifier — kept out since they no longer exist in the
    # anchor-tier output and would just generate "missing column" warnings).
    "anchor_mu_lof", "anchor_sd_lof", "anchor_mu_lof_local", "anchor_mu_lof_global",
    "anchor_weight_lof_local", "anchor_mu_noimpact", "anchor_sd_noimpact",
    "anchor_mu_noimpact_local", "anchor_mu_noimpact_global",
    "anchor_weight_noimpact_local", "anchor_direction", "anchor_lof_threshold",
    "anchor_post_lof", "anchor_call", "anchor_tier",
    "HGVSc", "HGVSp", "Protein_position", "Consequence", "Summary_Consequence",
    "Summary_Plot", "vcf_ref", "vcf_alt", "in_clinvar", "condition",
    "clnsig_norm", "clndn_raw", "clnsig_raw", "clinvar_id",
    "origin_raw", "origin_flags", "is_de_novo",
    "any_submission_de_novo", "n_submissions", "n_submissions_de_novo",
    "de_novo_submitters",
]

def is_plp(clnsig_norm: str) -> bool:
    """Return True if the ClinVar classification is Pathogenic or Likely pathogenic."""
    s = clnsig_norm.lower()
    return "pathogenic" in s and "benign" not in s

# ClinVar's ORIGIN INFO field is a bitmask (submitter-reported allele origin;
# multiple bits may be set at once, e.g. 33 = germline + de-novo).
# Reference: ClinVar VCF header / README_VCF.txt.
ORIGIN_BITS = {
    1:   "germline",
    2:   "somatic",
    4:   "inherited",
    8:   "paternal",
    16:  "maternal",
    32:  "de novo",
    64:  "biparental",
    128: "uniparental",
    256: "not tested",
    512: "tested inconclusive",
    1073741824: "other",
}
DE_NOVO_BIT = 32

def decode_origin(origin_raw: str) -> list:
    """Decode ClinVar's ORIGIN bitmask into a list of human-readable terms.
    Returns [] for missing/unparseable/unknown (0) values."""
    if not origin_raw or origin_raw in (".", ""):
        return []
    try:
        val = int(origin_raw)
    except ValueError:
        return []
    return [label for bit, label in ORIGIN_BITS.items() if val & bit]

def is_de_novo_origin(origin_raw: str) -> bool:
    """True if the de-novo bit (32) is set anywhere in ORIGIN."""
    if not origin_raw or origin_raw in (".", ""):
        return False
    try:
        val = int(origin_raw)
    except ValueError:
        return False
    return bool(val & DE_NOVO_BIT)

DEPLETION_MARKERS = {
    "strongly depleting": "v",  # downward triangle (anchor_tier)
    "weakly depleting":   "v",  # same shape as strongly depleting — distinguished
                                # by size/alpha via DEPLETION_STYLE below, so tier
                                # strength reads at a glance without a 4th shape
    "depleted":           "v",  # downward triangle (legacy GMM_status/stat_pos_raw)
    "no impact":          "o",   # circle
    "enriched":           "^",   # upward triangle
}

# Per-tier size modulation layered on top of DEPLETION_MARKERS, so that
# "weakly depleting" reads as a visibly smaller downward triangle than
# "strongly depleting" even though they share a marker shape — both drawn
# at full opacity so weak-tier points don't get lost. Anything not listed
# (e.g. an unexpected status string) gets size_mult=1.0, alpha=1.0.
DEPLETION_STYLE = {
    "strongly depleting": {"size_mult": 1.15, "alpha": 1.0},
    "weakly depleting":   {"size_mult": 0.55, "alpha": 1.0},
    "depleted":           {"size_mult": 1.15, "alpha": 1.0},  # legacy
    "enriched":           {"size_mult": 1.15, "alpha": 1.0},
    "no impact":          {"size_mult": 0.8,  "alpha": 0.55},
}
DEFAULT_DEPLETION_STYLE = {"size_mult": 1.0, "alpha": 1.0}

# ── Helper functions ──────────────────────────────────────────────────────────

def extract_variant_key(oligo_name: str) -> str | None:
    """Extract chrX:POS_VARIANT key from oligo_name for joining."""
    m = re.search(r"(chrX:\d+_[^_]+(?:>[A-Z]+)?)", str(oligo_name))
    return m.group(1) if m else None


def classify_condition(clndn: str) -> str:
    """Map a raw CLNDN string to one of the five condition categories."""
    if not clndn or clndn in (".", ""):
        return "No Condition Provided"
    s = clndn.lower()
    has_cdls = any(t in s for t in CDLS_TERMS)
    has_dee  = any(t in s for t in DEE_TERMS)
    if has_cdls and has_dee:
        return "Both"
    if has_cdls:
        return "CdLS only"
    if has_dee:
        return "DEE only"
    if any(t in s for t in OTHER_SMC1A_TERMS):
        return "Other SMC1A"
    return "Unclassified"


def normalise_clnsig(raw: str) -> str:
    """Normalise CLNSIG string (may contain slashes and underscores)."""
    if not raw:
        return ""
    # Replace underscores with spaces for display, keep slashes
    return raw.replace("_", " ")


def read_clinvar_vcf(vcf_path: str) -> pd.DataFrame:
    """
    Parse a per-targeton ClinVar VCF (bgzipped or plain).
    Returns DataFrame with: chrom, pos, ref, alt, clnsig, clndn, condition,
    clnsig_norm, clinvar_id (VCF ID column), origin_raw (raw ORIGIN bitmask
    string), origin_flags (decoded terms, semicolon-joined, e.g.
    "germline;de novo"), is_de_novo (True if the de-novo bit is set).
    """
    rows = []
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
            info = parts[7]
            clnsig = clndn = origin_raw = ""
            for field in info.split(";"):
                if field.startswith("CLNSIG="):
                    clnsig = field[7:]
                elif field.startswith("CLNDN="):
                    clndn = field[6:]
                elif field.startswith("ORIGIN="):
                    origin_raw = field[7:]
            rows.append({
                "chrom": chrom, "pos": pos,
                "vcf_ref": ref, "vcf_alt": alt,
                "clnsig_raw": clnsig, "clndn_raw": clndn,
                "clinvar_id": vid if vid not in (".", "") else None,
                "origin_raw": origin_raw,
            })
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df["condition"]   = df["clndn_raw"].apply(classify_condition)
    df["clnsig_norm"] = df["clnsig_raw"].apply(normalise_clnsig)
    df["origin_flags"] = df["origin_raw"].apply(lambda v: ";".join(decode_origin(v)))
    df["is_de_novo"]   = df["origin_raw"].apply(is_de_novo_origin)
    # Build a join key matching the SGE var_key format (chrX:POS_REF>ALT)
    # ClinVar VCF uses plain contig names (X not chrX)
    def make_var_key(row):
        chrom_str = "chrX" if row["chrom"] in ("X", "chrX") else f"chr{row['chrom']}"
        ref, alt = row["vcf_ref"], row["vcf_alt"]
        if len(ref) == len(alt) == 1:
            return f"{chrom_str}:{row['pos']}_{ref}>{alt}"
        elif len(alt) < len(ref):  # deletion
            n_del = len(ref) - len(alt)
            return f"{chrom_str}:{row['pos'] + 1}_{n_del}del"
        else:  # insertion
            n_ins = len(alt) - len(ref)
            return f"{chrom_str}:{row['pos'] + 1}_{n_ins}ins"
    df["var_key"] = df.apply(make_var_key, axis=1)
    return df


# ── ClinVar submission-level data (submission_summary.txt.gz) ───────────────
#
# One row per submitted record (SCV), as opposed to the VCF's one row per
# variant. Lets us ask "did AT LEAST ONE submission for this variant report
# de novo?", which can catch cases the VCF's aggregated ORIGIN bitmask
# collapses away. Downloaded separately from:
#   https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz
# Genome-wide (not gene-specific), so we filter to a caller-supplied set of
# VariationIDs (the ones we already found in this gene's per-targeton VCFs)
# while streaming, rather than loading the whole file into memory.

def parse_origin_counts(s: str) -> dict:
    """Parse an OriginCounts string like 'germline:1' or 'germline:5;de novo:1'
    into {origin_lowercase: count_or_None}. Returns {} for 'na'/'-'/empty."""
    if not s or s in ("na", "-", "."):
        return {}
    out = {}
    for part in s.split(";"):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            origin, count = part.rsplit(":", 1)
            origin = origin.strip().lower()
            try:
                count = int(count)
            except ValueError:
                count = None
        else:
            origin, count = part.lower(), None
        out[origin] = count
    return out


def origin_counts_has_de_novo(s: str) -> bool:
    """True if any origin key in an OriginCounts string is 'de novo'
    (tolerant of 'de-novo' hyphenation)."""
    for origin in parse_origin_counts(s):
        if origin.replace("-", " ").strip() == "de novo":
            return True
    return False


def collect_clinvar_variation_ids(vcf_dir: str, regions: pd.DataFrame) -> set:
    """Pre-scan every <Targeton_ID>_clinvar.vcf.gz in vcf_dir (for targetons
    listed in regions) to build the full set of VariationIDs relevant to this
    gene, so submission_summary.txt.gz (genome-wide) can be filtered to just
    these while streaming instead of being loaded in full."""
    ids = set()
    for _, reg in regions.iterrows():
        tid = reg["Targeton_ID"]
        vcf_path = os.path.join(vcf_dir, f"{tid}_clinvar.vcf.gz")
        if not os.path.exists(vcf_path):
            continue
        clinvar = read_clinvar_vcf(vcf_path)
        if clinvar.empty:
            continue
        ids.update(
            str(int(v)) if isinstance(v, str) and v.isdigit() else str(v)
            for v in clinvar["clinvar_id"].dropna()
        )
    return ids


def read_submission_summary(path: str, variation_ids: set = None) -> pd.DataFrame:
    """Stream-parse submission_summary.txt.gz (one row per SCV). If
    variation_ids is given, only rows whose VariationID is in that set are
    kept — essential given this file covers all of ClinVar, not just one
    gene. The file has '##' description lines, then a '#VariationID...'
    header line, then data rows (data rows do NOT start with '#')."""
    opener = gzip.open if path.endswith(".gz") else open
    header = None
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#VariationID"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if line.startswith("#"):
                continue
            if header is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                continue
            row = dict(zip(header, parts))
            if variation_ids is not None and row["VariationID"] not in variation_ids:
                continue
            rows.append(row)
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def summarize_submissions_by_variation(sub_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse submission_summary rows (one per SCV) down to one row per
    VariationID: any_submission_de_novo (bool), n_submissions,
    n_submissions_de_novo, de_novo_submitters (';'-joined, for provenance)."""
    if sub_df.empty:
        return pd.DataFrame(columns=[
            "VariationID", "n_submissions", "n_submissions_de_novo",
            "any_submission_de_novo", "de_novo_submitters"])

    sub_df = sub_df.copy()
    sub_df["submission_is_de_novo"] = sub_df["OriginCounts"].apply(origin_counts_has_de_novo)

    out_rows = []
    for vid, g in sub_df.groupby("VariationID"):
        dn = g[g["submission_is_de_novo"]]
        out_rows.append({
            "VariationID": vid,
            "n_submissions": len(g),
            "n_submissions_de_novo": len(dn),
            "any_submission_de_novo": len(dn) > 0,
            "de_novo_submitters": ";".join(sorted(dn["Submitter"].dropna().unique())),
        })
    return pd.DataFrame(out_rows)


def load_meta_consequences(meta_dir: str, targeton_id: str) -> pd.DataFrame:
    """Find and load the meta_consequences TSV for a given targeton.

    Join key: strip trailing _START_END coords from unique_oligo_name to produce
    oligo_name_base, which matches DESeq2 oligo_name exactly.
    """
    for fname in os.listdir(meta_dir):
        if fname.endswith("_meta_consequences.tsv"):
            path = os.path.join(meta_dir, fname)
            try:
                peek = pd.read_csv(path, sep="\t", usecols=["Targeton_ID"], nrows=1)
                if peek["Targeton_ID"].iloc[0] == targeton_id:
                    base_cols = ["unique_oligo_name", "Targeton_ID",
                                 "vcf_pos", "vcf_ref", "vcf_alt",
                                 "HGVSc", "HGVSp", "Protein_position",
                                 "Consequence", "Summary_Consequence"]
                    try:
                        # Summary_Plot has the LOF/Missense_Variant/etc. category
                        # set used for consequence colouring; older files may lack it.
                        df = pd.read_csv(path, sep="\t", usecols=base_cols + ["Summary_Plot"])
                    except ValueError:
                        df = pd.read_csv(path, sep="\t", usecols=base_cols)
                    # Strip trailing _start_end from unique_oligo_name to get
                    # the base name that matches DESeq2 oligo_name
                    df["oligo_name"] = df["unique_oligo_name"].str.replace(
                        r"_\d+_\d+$", "", regex=True)
                    # Deduplicate: keep one row per oligo_name
                    # (multiple equivalent_codon rows share the same oligo_name)
                    df = df.drop_duplicates(subset=["oligo_name"])
                    return df
            except Exception:
                continue
    return pd.DataFrame()


def get_depletion_column(df: pd.DataFrame) -> str:
    """Prefer anchor_tier (GMM/shrinkage anchor pipeline, 4 levels including
    'weakly depleting'); fall back to the older GMM_status (3-level), then
    to stat_pos_raw for older input files with neither anchor column."""
    if "anchor_tier" in df.columns:
        return "anchor_tier"
    if "GMM_status" in df.columns:
        return "GMM_status"
    return "stat_pos_raw"


# ── Plotting ──────────────────────────────────────────────────────────────────

def get_consequence_label(row) -> str:
    """Resolve the consequence-type label for a row: prefer 'Summary_Plot'
    (matches the LOF/Missense_Variant/etc. category set used for colouring),
    falling back to 'Summary_Consequence' (note: this labels stop-gained as
    'Nonsense_Variant' rather than 'LOF'), then 'Consequence', then 'Unknown'."""
    for col in ("Summary_Plot", "Summary_Consequence", "Consequence"):
        val = row.get(col)
        if val is not None and not (isinstance(val, float) and pd.isna(val)) and str(val).strip() != "":
            return str(val)
    return "Unknown"


def draw_variants(ax, df: pd.DataFrame, colour_by: str = "condition"):
    """Scatter one point per row of df onto ax. df must have genomic_pos, lfc,
    depletion_status, clnsig_norm columns, plus 'condition' (colour_by=
    'condition') or 'Summary_Consequence'/'Consequence' (colour_by=
    'consequence')."""
    for _, row in df.iterrows():
        x   = row["genomic_pos"]
        lfc = row["lfc"]
        status = str(row["depletion_status"]).lower().strip()
        clnsig = row["clnsig_norm"]

        if colour_by == "consequence":
            label  = get_consequence_label(row)
            colour = CONSEQUENCE_COLOURS.get(label, CONSEQUENCE_OTHER_COLOUR)
        else:
            condition = row["condition"]
            colour = CONDITION_COLOURS.get(condition, CONDITION_COLOURS["Unclassified"])

        marker = DEPLETION_MARKERS.get(status, "o")
        tier_style = DEPLETION_STYLE.get(status, DEFAULT_DEPLETION_STYLE)

        plp = is_plp(clnsig)
        size = 60 * tier_style["size_mult"]
        lw   = 1.5 if plp else 0.3
        ec   = "black" if plp else "#bbbbbb"

        ax.scatter(x, lfc,
                   marker=marker, s=size,
                   color=colour, edgecolors=ec,
                   linewidths=lw, alpha=tier_style["alpha"],
                   zorder=3)


def condition_legend_handles(present_conditions, exclude=("No Condition Provided",)):
    return [
        mpatches.Patch(facecolor=col, label=label, alpha=0.85)
        for label, col in CONDITION_COLOURS.items()
        if label in present_conditions and label not in exclude
    ]


def consequence_legend_handles(df: pd.DataFrame):
    """Build legend patches for whichever consequence labels are present in df."""
    labels_present = set(df.apply(get_consequence_label, axis=1)) if not df.empty else set()
    handles = [
        mpatches.Patch(facecolor=col, label=label.replace("_", " "), alpha=0.85)
        for label, col in CONSEQUENCE_COLOURS.items()
        if label in labels_present
    ]
    if any(lbl not in CONSEQUENCE_COLOURS for lbl in labels_present):
        handles.append(mpatches.Patch(facecolor=CONSEQUENCE_OTHER_COLOUR,
                                      label="Other / unclassified", alpha=0.85))
    return handles


def depletion_legend_handles(present_statuses):
    handles = []
    for label, mk in DEPLETION_MARKERS.items():
        if label not in present_statuses:
            continue
        style = DEPLETION_STYLE.get(label, DEFAULT_DEPLETION_STYLE)
        handles.append(Line2D(
            [0], [0], marker=mk, color="none", markerfacecolor="#666",
            markersize=7 * (style["size_mult"] ** 0.5),
            markeredgewidth=0.5, markeredgecolor="#444",
            alpha=style["alpha"], label=label))
    return handles


PLP_LEGEND_HANDLES = [
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="#888888", markersize=8,
           markeredgewidth=1.5, markeredgecolor="black",
           label="Pathogenic / Likely pathogenic"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="#888888", markersize=8,
           markeredgewidth=0.3, markeredgecolor="#bbbbbb",
           label="Other ClinVar"),
]


def plot_lollipop(merged: pd.DataFrame, targeton_id: str, region: str, outpath: str,
                   restrict_conditions=None, title_suffix=None, depletion_col_label="anchor_tier",
                   colour_by="condition"):
    """
    Single-panel lollipop plot for one targeton.
    merged columns required: genomic_pos, lfc, depletion_status, condition,
                              clnsig_norm, in_clinvar
    restrict_conditions: optional list of condition labels to limit plotting to
                         (e.g. ["DEE only"]).  x/y ranges are still computed
                         from the full `merged` targeton data so plots for the
                         same targeton stay visually comparable.
    colour_by: "condition" (default, used for the full multi-condition plot)
               or "consequence" (used when restrict_conditions narrows the
               plot to a single condition, where a condition-colour legend
               would be redundant - one colour for every point).
    """
    fig, ax = plt.subplots(figsize=(14, 5.5))
    fig.patch.set_facecolor("#fafafa")
    ax.set_facecolor("#fafafa")

    ax.axhline(0, color="#cccccc", lw=0.8, zorder=0)
    lfc_max = merged["lfc"].abs().max()
    ax.axhspan(-0.5, 0.5, color="#f0f0f0", alpha=0.5, zorder=0)

    plot_df = merged[merged["in_clinvar"]].copy()
    if restrict_conditions is not None:
        plot_df = plot_df[plot_df["condition"].isin(restrict_conditions)]

    # x-range from the whole targeton (not just the plotted subset) so that
    # the full plot and any condition-split plot for the same targeton share
    # the same genomic-position axis.
    x_vals = merged["genomic_pos"]
    x_lo = x_vals.min() - (x_vals.max() - x_vals.min()) * 0.03
    x_hi = x_vals.max() + (x_vals.max() - x_vals.min()) * 0.03

    draw_variants(ax, plot_df, colour_by=colour_by)

    # ── Legends ──────────────────────────────────────────────────────────────
    dep_handles = depletion_legend_handles(
        plot_df["depletion_status"].astype(str).str.lower().values)

    if colour_by == "consequence":
        first_handles = consequence_legend_handles(plot_df)
        first_title = "Variant consequence"
    else:
        present_conditions = plot_df["condition"].unique()
        first_handles = condition_legend_handles(present_conditions)
        first_title = "Disease condition"

    if first_handles:
        leg1 = ax.legend(handles=first_handles, title=first_title,
                         loc="upper left", fontsize=6.5, title_fontsize=7.5,
                         framealpha=0.9, ncol=1,
                         bbox_to_anchor=(1.01, 1), borderaxespad=0)
        ax.add_artist(leg1)

    if dep_handles:
        leg2 = ax.legend(handles=dep_handles, title=f"Depletion status ({depletion_col_label})",
                         loc="center left", fontsize=6.5, title_fontsize=7.5,
                         framealpha=0.9,
                         bbox_to_anchor=(1.01, 0.45), borderaxespad=0)
        ax.add_artist(leg2)

    ax.legend(handles=PLP_LEGEND_HANDLES, title="ClinVar significance",
              loc="lower left", fontsize=6.5, title_fontsize=7.5,
              framealpha=0.9,
              bbox_to_anchor=(1.01, 0.0), borderaxespad=0)

    # ── Axes ─────────────────────────────────────────────────────────────────
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(-(lfc_max + 0.8), lfc_max + 0.8)
    ax.set_xlabel(f"Genomic position (chrX) — {region}", fontsize=9)
    ax.set_ylabel("LFC (pos_adj_log2FoldChange_raw)", fontsize=9)
    title = f"Targeton {targeton_id}  ·  {title_suffix}" if title_suffix else \
            f"Targeton {targeton_id}  ·  SGE vs ClinVar intersection"
    ax.set_title(title, fontsize=11, fontweight="bold", pad=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)

    plt.tight_layout(rect=[0, 0, 0.82, 1])
    plt.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"    Saved: {outpath}")


def plot_stacked_conditions(merged: pd.DataFrame, targeton_id: str, region: str,
                             outpath: str, depletion_col_label="anchor_tier"):
    """
    Two-panel stacked figure (shared x-axis): DEE-only variants on top,
    CdLS-only variants on the bottom. x/y ranges shared across panels and
    matched to the full targeton range for comparability with the main plot.
    """
    plot_df = merged[merged["in_clinvar"]].copy()
    dee_df  = plot_df[plot_df["condition"] == "DEE only"]
    cdls_df = plot_df[plot_df["condition"] == "CdLS only"]

    fig, (ax_top, ax_bottom) = plt.subplots(
        2, 1, figsize=(14, 9), sharex=True, sharey=True)
    fig.patch.set_facecolor("#fafafa")

    lfc_max = merged["lfc"].abs().max()
    x_vals = merged["genomic_pos"]
    x_lo = x_vals.min() - (x_vals.max() - x_vals.min()) * 0.03
    x_hi = x_vals.max() + (x_vals.max() - x_vals.min()) * 0.03

    panels = [
        (ax_top,    dee_df,  "DEE-only variants",  CONDITION_COLOURS["DEE only"]),
        (ax_bottom, cdls_df, "CdLS-only variants", CONDITION_COLOURS["CdLS only"]),
    ]

    for ax, df, label, colour in panels:
        ax.set_facecolor("#fafafa")
        ax.axhline(0, color="#cccccc", lw=0.8, zorder=0)
        ax.axhspan(-0.5, 0.5, color="#f0f0f0", alpha=0.5, zorder=0)
        draw_variants(ax, df, colour_by="consequence")
        ax.set_ylabel("LFC (pos_adj_log2FoldChange_raw)", fontsize=9)
        ax.set_title(label, loc="left", fontsize=10, fontweight="bold", color=colour)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=7)

    ax_top.set_xlim(x_lo, x_hi)
    ax_top.set_ylim(-(lfc_max + 0.8), lfc_max + 0.8)
    ax_bottom.set_xlabel(f"Genomic position (chrX) — {region}", fontsize=9)

    # Shared legends: condition is already indicated by each panel's
    # title/colour, so points here are coloured by variant consequence type
    # instead, plus the usual depletion-shape and P/LP-border legends.
    cons_handles = consequence_legend_handles(plot_df)
    dep_handles = depletion_legend_handles(
        plot_df["depletion_status"].astype(str).str.lower().values)

    if cons_handles:
        leg0 = fig.legend(handles=cons_handles, title="Variant consequence",
                          loc="upper left", fontsize=6.5, title_fontsize=7.5,
                          framealpha=0.9, bbox_to_anchor=(0.83, 0.97), borderaxespad=0)
        fig.add_artist(leg0)
    if dep_handles:
        leg1 = fig.legend(handles=dep_handles, title=f"Depletion status ({depletion_col_label})",
                          loc="upper left", fontsize=6.5, title_fontsize=7.5,
                          framealpha=0.9, bbox_to_anchor=(0.83, 0.55), borderaxespad=0)
        fig.add_artist(leg1)
    fig.legend(handles=PLP_LEGEND_HANDLES, title="ClinVar significance",
              loc="upper left", fontsize=6.5, title_fontsize=7.5,
              framealpha=0.9, bbox_to_anchor=(0.83, 0.30), borderaxespad=0)

    fig.suptitle(f"Targeton {targeton_id}  ·  DEE vs CdLS ClinVar variants",
                 fontsize=12, fontweight="bold")

    plt.tight_layout(rect=[0, 0, 0.82, 0.96])
    plt.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"    Saved: {outpath}")


def generate_condition_split_plot(merged: pd.DataFrame, targeton_id: str, region: str,
                                   outdir: str, depletion_col_label="anchor_tier"):
    """
    Decide which condition-split plot to make for this targeton:
      - both DEE-only and CdLS-only present -> stacked two-panel plot
      - only one of the two present         -> single lollipop plot, limited to it
      - neither present                     -> skip
    """
    plot_df = merged[merged["in_clinvar"]]
    has_dee  = (plot_df["condition"] == "DEE only").any()
    has_cdls = (plot_df["condition"] == "CdLS only").any()

    if has_dee and has_cdls:
        outpath = os.path.join(outdir, f"{targeton_id}_DEE_CdLS_stacked.png")
        plot_stacked_conditions(merged, targeton_id, region, outpath, depletion_col_label)
    elif has_dee:
        outpath = os.path.join(outdir, f"{targeton_id}_DEE_only.png")
        plot_lollipop(merged, targeton_id, region, outpath,
                      restrict_conditions=["DEE only"],
                      title_suffix="DEE variants only",
                      depletion_col_label=depletion_col_label,
                      colour_by="consequence")
    elif has_cdls:
        outpath = os.path.join(outdir, f"{targeton_id}_CdLS_only.png")
        plot_lollipop(merged, targeton_id, region, outpath,
                      restrict_conditions=["CdLS only"],
                      title_suffix="CdLS variants only",
                      depletion_col_label=depletion_col_label,
                      colour_by="consequence")
    else:
        print("    No DEE-only or CdLS-only ClinVar variants — skipping condition-split plot")


# ── Main ──────────────────────────────────────────────────────────────────────

def find_deseq2_files(deseq2_dir: str) -> dict:
    """
    Scan a directory for per-targeton DESeq2 result TSVs.
    Targeton_ID is extracted as everything before the first underscore.
    Returns {targeton_id: filepath}.
    """
    mapping = {}
    for fname in sorted(os.listdir(deseq2_dir)):
        if not fname.endswith(".tsv"):
            continue
        tid = fname.split("_")[0]
        if tid:
            mapping[tid] = os.path.join(deseq2_dir, fname)
    return mapping


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--deseq2_dir", required=True,
                        help="Directory of per-targeton DESeq2 TSVs "
                             "(filenames must start with Targeton_ID)")
    parser.add_argument("--meta_dir", required=True,
                        help="Directory with *_meta_consequences.tsv files")
    parser.add_argument("--vcf_dir",  required=True,
                        help="Directory with <Targeton_ID>_clinvar.vcf.gz files")
    parser.add_argument("--regions",  required=True,
                        help="targeton_regions.tsv from extract_targeton_regions.py")
    parser.add_argument("--outdir",   default="clinvar_intersect_plots",
                        help="Output directory (default: clinvar_intersect_plots/)")
    parser.add_argument("--de_novo_only", action="store_true",
                        help="Restrict the PLOTS (lollipop/stacked) to ClinVar "
                             "variants with the de-novo bit set in ORIGIN. Does "
                             "NOT affect the annotated TSV or the cross-targeton "
                             "summary TSV, which always include all matched "
                             "variants (with origin_raw/origin_flags/is_de_novo "
                             "columns) so you can filter on your own terms "
                             "downstream. Note ORIGIN reflects what ClinVar "
                             "submitters chose to report, so this is a subset "
                             "of true de-novo occurrences, not a complete count.")
    parser.add_argument("--submission_summary", default=None,
                        help="Path to ClinVar's submission_summary.txt.gz "
                             "(https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/"
                             "submission_summary.txt.gz). Genome-wide, so it is "
                             "pre-filtered to this gene's VariationIDs (collected "
                             "from the per-targeton VCFs in --vcf_dir) before "
                             "being fully loaded. Adds any_submission_de_novo "
                             "(True if AT LEAST ONE submission for that variant "
                             "reported de novo origin, which the VCF's aggregated "
                             "ORIGIN field can sometimes miss), n_submissions, "
                             "n_submissions_de_novo, and de_novo_submitters to "
                             "the annotated/summary TSVs. If provided, "
                             "--de_novo_only also honours this column (in "
                             "addition to the VCF-level is_de_novo) when "
                             "restricting plots. Optional — omit to skip "
                             "entirely.")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Discover per-targeton DESeq2 files
    print("Scanning DESeq2 directory...")
    deseq2_files = find_deseq2_files(args.deseq2_dir)
    print(f"  Found {len(deseq2_files)} targeton result file(s)")

    # Load targeton regions
    regions = pd.read_csv(args.regions, sep="\t")
    print(f"  {len(regions)} targetons in regions file")

    # Optional: submission-level ClinVar data (one row per SCV), for
    # any_submission_de_novo. Pre-scan this gene's per-targeton VCFs for the
    # VariationIDs we actually need before streaming the (genome-wide) file.
    submission_agg = pd.DataFrame()
    if args.submission_summary:
        print("\nScanning per-targeton ClinVar VCFs for VariationIDs "
              "(to filter submission_summary.txt.gz)...")
        variation_ids = collect_clinvar_variation_ids(args.vcf_dir, regions)
        print(f"  {len(variation_ids)} distinct VariationID(s) found across all targetons")
        print(f"Reading {args.submission_summary} (filtered to those IDs)...")
        sub_df = read_submission_summary(args.submission_summary, variation_ids=variation_ids)
        print(f"  {len(sub_df)} matching submission row(s) (one per SCV)")
        submission_agg = summarize_submissions_by_variation(sub_df)
        n_dn_variants = int(submission_agg["any_submission_de_novo"].sum()) if not submission_agg.empty else 0
        print(f"  {n_dn_variants} / {len(submission_agg)} VariationID(s) have "
              f"at least one submission reporting de novo origin")

    # Collects ClinVar-matched rows across all targetons for the combined
    # cross-targeton summary TSV written at the end.
    summary_rows = []

    # Process each targeton
    for _, reg in regions.iterrows():
        tid   = reg["Targeton_ID"]
        chrom = reg["chrom"]
        start = int(reg["start"])
        end   = int(reg["end"])
        region_str = f"{chrom}:{start}-{end}"

        print(f"\n→ {tid}  {region_str}")

        # Load this targeton's DESeq2 file
        if tid not in deseq2_files:
            print(f"  No DESeq2 result file found for {tid} — skipping")
            continue
        deseq_t = pd.read_csv(deseq2_files[tid], sep="\t")
        deseq_t["var_key"] = deseq_t["oligo_name"].apply(extract_variant_key)
        deseq_t = deseq_t.dropna(subset=["var_key"])  # drop PAM/ref controls
        if deseq_t.empty:
            print("  No variant rows after dropping controls — skipping")
            continue
        print(f"  {len(deseq_t)} DESeq2 variants")

        depletion_col = get_depletion_column(deseq_t)
        print(f"  Using '{depletion_col}' as depletion status")

        # Load meta_consequences for ref/alt and annotation
        meta = load_meta_consequences(args.meta_dir, tid)
        if meta.empty:
            print(f"  WARNING: no meta_consequences file found for {tid} — skipping")
            continue
        print(f"  {len(meta)} meta_consequence rows")

        # Join DESeq2 + meta on oligo_name to get HGVSc/HGVSp/Protein_position
        # meta["oligo_name"] = unique_oligo_name with trailing _start_end stripped
        meta_join_cols = ["oligo_name", "vcf_pos", "vcf_ref", "vcf_alt",
                          "HGVSc", "HGVSp", "Protein_position",
                          "Consequence", "Summary_Consequence"]
        if "Summary_Plot" in meta.columns:
            meta_join_cols.append("Summary_Plot")
        deseq_t = deseq_t.merge(meta[meta_join_cols], on="oligo_name", how="left")

        # Load per-targeton ClinVar VCF
        vcf_path = os.path.join(args.vcf_dir, f"{tid}_clinvar.vcf.gz")
        if not os.path.exists(vcf_path):
            print(f"  WARNING: VCF not found at {vcf_path} — ClinVar track will be empty")
            clinvar = pd.DataFrame()
        else:
            clinvar = read_clinvar_vcf(vcf_path)
            print(f"  {len(clinvar)} ClinVar variants in VCF")

        # Intersect DESeq2 variants with ClinVar on vcf_pos + vcf_ref + vcf_alt
        # (from meta join above) for exact allele matching
        if not clinvar.empty:
            clinvar_merge = clinvar[["pos", "vcf_ref", "vcf_alt",
                                     "condition", "clnsig_norm", "clndn_raw",
                                     "clnsig_raw", "clinvar_id",
                                     "origin_raw", "origin_flags", "is_de_novo"]].rename(
                                     columns={"pos": "vcf_pos"})
            deseq_t = deseq_t.merge(
                clinvar_merge,
                on=["vcf_pos", "vcf_ref", "vcf_alt"], how="left"
            )
            deseq_t["in_clinvar"] = deseq_t["condition"].notna()
            deseq_t["condition"]  = deseq_t["condition"].fillna("No Condition Provided")
            deseq_t["clnsig_norm"] = deseq_t["clnsig_norm"].fillna("")
            deseq_t["origin_raw"]   = deseq_t["origin_raw"].fillna("")
            deseq_t["origin_flags"] = deseq_t["origin_flags"].fillna("")
            # is_de_novo is False for both "not matched to ClinVar" and
            # "matched but no de-novo bit set" — same semantics as in_clinvar,
            # but note downstream that is_de_novo=False alone doesn't
            # distinguish "No Condition Provided" from "in ClinVar, origin unspecified".
            deseq_t["is_de_novo"] = deseq_t["is_de_novo"].fillna(False).astype(bool)
            # clinvar_id becomes float64 after the merge (NaN forces this even
            # though IDs are integer-like) - restore clean integer-looking strings
            deseq_t["clinvar_id"] = deseq_t["clinvar_id"].apply(
                lambda v: "" if pd.isna(v) else
                (str(int(v)) if float(v).is_integer() else str(v)))
        else:
            deseq_t["in_clinvar"]   = False
            deseq_t["condition"]    = "No Condition Provided"
            deseq_t["clnsig_norm"]  = ""
            deseq_t["clinvar_id"]   = ""
            deseq_t["origin_raw"]   = ""
            deseq_t["origin_flags"] = ""
            deseq_t["is_de_novo"]   = False

        n_cv = deseq_t["in_clinvar"].sum()
        print(f"  {n_cv} / {len(deseq_t)} DESeq2 variants matched to ClinVar")
        n_dn = deseq_t["is_de_novo"].sum()
        print(f"  {n_dn} / {n_cv} ClinVar-matched variants flagged de novo (ORIGIN bit 32)")

        # ── Merge submission-level de-novo data (optional) ──────────────────
        if not submission_agg.empty:
            deseq_t = deseq_t.merge(
                submission_agg.rename(columns={"VariationID": "clinvar_id"}),
                on="clinvar_id", how="left"
            )
            deseq_t["any_submission_de_novo"] = deseq_t["any_submission_de_novo"].fillna(False).astype(bool)
            deseq_t["n_submissions"] = deseq_t["n_submissions"].fillna(0).astype(int)
            deseq_t["n_submissions_de_novo"] = deseq_t["n_submissions_de_novo"].fillna(0).astype(int)
            deseq_t["de_novo_submitters"] = deseq_t["de_novo_submitters"].fillna("")
            n_sub_dn = deseq_t["any_submission_de_novo"].sum()
            print(f"  {n_sub_dn} / {n_cv} ClinVar-matched variants have >=1 "
                  f"submission reporting de novo (submission_summary.txt.gz)")
        elif args.submission_summary:
            # Flag was given but nothing matched/loaded — still add the
            # columns (all False/0/"") so the TSV schema stays consistent
            # across targetons regardless of whether this one had any hits.
            deseq_t["any_submission_de_novo"] = False
            deseq_t["n_submissions"] = 0
            deseq_t["n_submissions_de_novo"] = 0
            deseq_t["de_novo_submitters"] = ""

        # ── Write annotated TSV ──────────────────────────────────────────────
        # Columns added vs original DESeq2 file:
        #   from meta: HGVSc, HGVSp, Protein_position, Consequence,
        #              Summary_Consequence, Summary_Plot, vcf_ref, vcf_alt
        #   from ClinVar VCF: in_clinvar, condition, clnsig_norm, clndn_raw,
        #                     clnsig_raw, clinvar_id, origin_raw, origin_flags,
        #                     is_de_novo
        #   from submission_summary.txt.gz (if --submission_summary given):
        #     any_submission_de_novo, n_submissions, n_submissions_de_novo,
        #     de_novo_submitters
        tsv_cols = (
            # Original DESeq2 columns (all present)
            [c for c in deseq_t.columns if c in pd.read_csv(
                deseq2_files[tid], sep="\t", nrows=0).columns]
            # New annotation columns
            + [c for c in ["HGVSc", "HGVSp", "Protein_position",
                            "Consequence", "Summary_Consequence", "Summary_Plot",
                            "vcf_ref", "vcf_alt",
                            "in_clinvar", "condition",
                            "clnsig_norm", "clndn_raw", "clnsig_raw",
                            "clinvar_id", "origin_raw", "origin_flags", "is_de_novo",
                            "any_submission_de_novo", "n_submissions",
                            "n_submissions_de_novo", "de_novo_submitters"]
               if c in deseq_t.columns]
        )
        tsv_path = os.path.join(args.outdir, f"{tid}_clinvar_annotated.tsv")
        deseq_t[tsv_cols].to_csv(tsv_path, sep="\t", index=False)
        print(f"  Annotated TSV: {tsv_path}")

        # ── Collect ClinVar-matched rows for the cross-targeton summary ────────
        clinvar_only = deseq_t[deseq_t["in_clinvar"]].copy()
        if not clinvar_only.empty:
            exon = TARGETON_EXON.get(tid)
            if exon is None:
                print(f"  WARNING: no exon mapping found for targeton '{tid}' — "
                      f"'exon' will be blank in the summary TSV for these rows")
            missing_cols = [c for c in SUMMARY_COLUMNS if c No Condition Provided_only.columns]
            if missing_cols:
                print(f"  WARNING: summary TSV missing expected column(s): {missing_cols}")
            present_cols = [c for c in SUMMARY_COLUMNS if c in clinvar_only.columns]
            row_block = clinvar_only[present_cols].copy()
            row_block.insert(0, "Targeton_ID", tid)
            row_block.insert(1, "exon", exon)
            summary_rows.append(row_block)

        # Rename columns for plotting
        deseq_t = deseq_t.rename(columns={
            "position":                    "genomic_pos",
            "pos_adj_log2FoldChange_raw":  "lfc",
            depletion_col:                 "depletion_status",
        })
        deseq_t = deseq_t.dropna(subset=["lfc", "genomic_pos"])

        if deseq_t.empty:
            print("  No plottable rows — skipping")
            continue

        # --de_novo_only restricts the PLOTS only (not the TSVs written above,
        # which always retain every ClinVar match with origin_raw/origin_flags/
        # is_de_novo [and any_submission_de_novo, if --submission_summary was
        # given] so you can filter downstream on your own terms). Variants
        # that matched ClinVar but aren't de-novo-flagged by EITHER source
        # fall back to the "No Condition Provided" background category for plotting.
        deseq_plot = deseq_t.copy()
        if args.de_novo_only:
            de_novo_flag = deseq_plot["is_de_novo"]
            if "any_submission_de_novo" in deseq_plot.columns:
                de_novo_flag = de_novo_flag | deseq_plot["any_submission_de_novo"]
            keep_cv = deseq_plot["in_clinvar"] & de_novo_flag
            n_dropped = int((deseq_plot["in_clinvar"] & ~keep_cv).sum())
            print(f"  --de_novo_only: {n_dropped} ClinVar-matched but non-de-novo "
                  f"variant(s) moved to background for plotting")
            deseq_plot.loc[~keep_cv, "condition"] = "No Condition Provided"
            deseq_plot["in_clinvar"] = keep_cv

        outpath = os.path.join(args.outdir, f"{tid}_clinvar_intersect.png")
        plot_lollipop(deseq_plot, tid, region_str, outpath,
                      depletion_col_label=depletion_col)

        generate_condition_split_plot(deseq_plot, tid, region_str, args.outdir,
                                      depletion_col_label=depletion_col)

    if summary_rows:
        summary_df = pd.concat(summary_rows, ignore_index=True)
        summary_path = os.path.join(args.outdir, "clinvar_variants_summary.tsv")
        summary_df.to_csv(summary_path, sep="\t", index=False)
        print(f"\nClinVar variant summary across all targetons: {summary_path} "
              f"({len(summary_df)} variants)")
    else:
        print("\nNo ClinVar-matched variants found across any targeton — "
              "summary TSV not created")

    print(f"\nDone. Plots and annotated TSVs saved to: {os.path.abspath(args.outdir)}/")


if __name__ == "__main__":
    main()
