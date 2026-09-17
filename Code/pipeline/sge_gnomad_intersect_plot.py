#!/usr/bin/env python3
"""
sge_gnomad_intersect_plot.py

Intersects SGE DESeq2/GMM-anchor results with gnomAD (exomes + genomes) for
every assay variant across all targetons, and tests whether the GMM's
depleted/enriched calls are under-represented in the population relative to
"no impact" — the purifying-selection sanity check requested alongside the
in-silico concordance analysis.

Mirrors sge_clinvar_intersect_plot_gaussian.py's structure/conventions
(same find_deseq2_files / load_meta_consequences / get_depletion_column
helpers, same merge key: vcf_pos + vcf_ref + vcf_alt from the meta join),
with two differences suited to this analysis:

  1. TWO gnomAD sources per targeton (exomes + genomes), pooled per variant,
     rather than one ClinVar VCF — see combine_gnomad_sources().
  2. Unlike the ClinVar summary (matched variants only), the cross-targeton
     summary TSV here includes EVERY assayed variant (matched or not),
     because testing for absence requires the full denominator, not just
     the hits.

IMPORTANT — verify the INFO field names before trusting the output:
  This script assumes gnomAD v4.1 VCF INFO field names (AC, AN, AF,
  grpmax, AF_grpmax). These are documented gnomAD v4.1 fields but have NOT
  been verified against your actual downloaded VCFs. Before running in
  earnest:
      bcftools view -h <TID>_gnomad_exomes.vcf.gz | grep '^##INFO'
  and adjust --ac_field / --an_field / --af_field / --grpmax_af_field /
  --grpmax_pop_field if the names differ. The script prints a warning
  (once per file) if a configured field is never found.

Usage:
    python sge_gnomad_intersect_plot.py \\
        --deseq2_dir     all_deseq2_results/ \\
        --meta_dir       /path/to/meta_consequence/files \\
        --gnomad_vcf_dir gnomad_vcfs/ \\
        --regions        targeton_regions.tsv \\
        --outdir         gnomad_intersect_results/

Outputs (in --outdir):
    <TID>_gnomad_annotated.tsv   Every DESeq2/anchor row for that targeton,
                                  + gnomAD exomes/genomes/pooled columns.
    sge_gnomad_summary.tsv       ALL assay variants across all targetons,
                                  one row each, gnomAD + anchor_tier columns
                                  (this is the full denominator table the
                                  stats below are computed from).
    gnomad_absence_stats.tsv     Per-tier gnomAD-match rate + the
                                  depleted/enriched-vs-no-impact contingency
                                  test (Fisher's exact) results.
    gnomad_match_rate_by_tier.png     Bar chart: % gnomAD-matched per tier.
    gnomad_af_by_tier.png             Boxplot: log10(pooled AF) per tier,
                                       gnomAD-matched variants only.

Requirements:
    pip install pandas numpy scipy matplotlib
"""

import argparse
import gzip
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, chi2_contingency, mannwhitneyu

# ── Exon lookup (for the cross-targeton summary), same mapping as the
#    ClinVar script — keep in sync if the exon map changes.
TARGETON_EXON = {
    "SBPS": 1,  "APDY": 2,  "NBXL": 3,  "RSDC": 4,  "AXFG": 5,  "EXTP": 5,
    "SWAZ": 6,  "CQEJ": 6,  "IHMH": 7,  "KJQT": 8,  "YTUB": 9,  "RYFS": 10,
    "PSYW": 11, "NQBU": 12, "NUIB": 13, "SUHO": 14, "HDOI": 15, "MGTQ": 16,
    "RVMB": 17, "ELWX": 18, "OJTT": 19, "ZANL": 20, "RBMZ": 21, "CNJV": 22,
    "NLVE": 23, "FKKY": 24, "FBSZ": 25,
}

TIER_ORDER = ["no impact", "weakly depleting", "strongly depleting", "enriched"]
TIER_COLOURS = {
    "no impact":          "#9b9b9b",
    "weakly depleting":   "#f5a623",
    "strongly depleting": "#e05c4e",
    "enriched":           "#4caf50",
}


# ── gnomAD VCF reading ────────────────────────────────────────────────────────

def _parse_info(info_str: str) -> dict:
    """Split a VCF INFO string into a dict. Flag fields (no '=') map to True."""
    out = {}
    for field in info_str.split(";"):
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
            out[k] = v
        else:
            out[field] = True
    return out


def _first(value):
    """gnomAD INFO fields are normally already biallelic-split, but take the
    first comma-separated value defensively in case any slipped through."""
    if value is None:
        return None
    if isinstance(value, str) and "," in value:
        return value.split(",")[0]
    return value


def _to_float(value):
    v = _first(value)
    if v is None:
        return np.nan
    try:
        return float(v)
    except (ValueError, TypeError):
        return np.nan


def _to_int(value):
    v = _first(value)
    if v is None:
        return np.nan
    try:
        return int(float(v))
    except (ValueError, TypeError):
        return np.nan


def read_gnomad_vcf(vcf_path: str, source_label: str,
                    ac_field="AC", an_field="AN", af_field="AF",
                    grpmax_af_field="AF_grpmax", grpmax_pop_field="grpmax",
                    pass_only=True) -> pd.DataFrame:
    """
    Parse a per-targeton gnomAD sites VCF (bgzipped or plain).
    Returns DataFrame with: pos, vcf_ref, vcf_alt, filter, AC, AN, AF,
    AF_grpmax, grpmax_pop, source (source_label).
    Rows failing FILTER (i.e. FILTER != PASS) are dropped unless
    pass_only=False.
    """
    if not os.path.exists(vcf_path):
        return pd.DataFrame()

    rows = []
    field_hits = {ac_field: 0, an_field: 0, af_field: 0,
                  grpmax_af_field: 0, grpmax_pop_field: 0}
    n_total = 0
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            n_total += 1
            chrom, pos, vid, ref, alt, qual, flt, info = parts[:8]
            if pass_only and flt not in ("PASS", "."):
                continue
            info_d = _parse_info(info)
            for f in field_hits:
                if f in info_d:
                    field_hits[f] += 1
            rows.append({
                "pos": int(pos),
                "vcf_ref": ref,
                "vcf_alt": alt,
                "filter": flt,
                "AC": _to_int(info_d.get(ac_field)),
                "AN": _to_int(info_d.get(an_field)),
                "AF": _to_float(info_d.get(af_field)),
                "AF_grpmax": _to_float(info_d.get(grpmax_af_field)),
                "grpmax_pop": _first(info_d.get(grpmax_pop_field)),
                "source": source_label,
            })

    # Warn (once per file) if a configured field was never seen at all —
    # likely means the field name doesn't match this VCF's actual schema.
    if n_total > 0:
        for f, hits in field_hits.items():
            if hits == 0:
                print(f"    WARNING: INFO field '{f}' not found in any of "
                      f"{n_total} record(s) in {os.path.basename(vcf_path)} — "
                      f"check field name with 'bcftools view -h {vcf_path} | grep INFO'",
                      file=sys.stderr)

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def combine_gnomad_sources(exomes: pd.DataFrame, genomes: pd.DataFrame) -> pd.DataFrame:
    """
    Outer-merge exomes + genomes on (pos, vcf_ref, vcf_alt) and pool:
      - pooled_AC / pooled_AN = summed across whichever source(s) have the
        variant; pooled_AF = pooled_AC / pooled_AN.
      - combined_grpmax_af = max(exomes AF_grpmax, genomes AF_grpmax) — a
        conservative "most common in any single population group, in
        either callset" figure, NOT a single well-defined population
        (exomes and genomes grpmax can come from different population
        groups). Use pooled_AF for a single coherent frequency; use
        combined_grpmax_af only as a "how common has this ever been seen
        to be, anywhere" upper-bound check.
      - gnomad_source: "exomes", "genomes", or "exomes+genomes".
    Returns one row per distinct (pos, vcf_ref, vcf_alt); empty DataFrame
    if both inputs are empty.
    """
    if exomes.empty and genomes.empty:
        return pd.DataFrame(columns=[
            "pos", "vcf_ref", "vcf_alt", "in_gnomad", "gnomad_source",
            "gnomad_exomes_af", "gnomad_genomes_af", "pooled_ac", "pooled_an",
            "pooled_af", "combined_grpmax_af", "grpmax_pop_exomes",
            "grpmax_pop_genomes"])

    ex = exomes.rename(columns={c: f"{c}_exomes" for c in exomes.columns
                                if c not in ("pos", "vcf_ref", "vcf_alt")}) if not exomes.empty \
        else pd.DataFrame(columns=["pos", "vcf_ref", "vcf_alt"])
    gn = genomes.rename(columns={c: f"{c}_genomes" for c in genomes.columns
                                 if c not in ("pos", "vcf_ref", "vcf_alt")}) if not genomes.empty \
        else pd.DataFrame(columns=["pos", "vcf_ref", "vcf_alt"])

    merged = pd.merge(ex, gn, on=["pos", "vcf_ref", "vcf_alt"], how="outer")

    for col in ["AC_exomes", "AN_exomes", "AF_exomes", "AF_grpmax_exomes",
                "AC_genomes", "AN_genomes", "AF_genomes", "AF_grpmax_genomes"]:
        if col not in merged.columns:
            merged[col] = np.nan
    if "grpmax_pop_exomes" not in merged.columns:
        merged["grpmax_pop_exomes"] = None
    if "grpmax_pop_genomes" not in merged.columns:
        merged["grpmax_pop_genomes"] = None

    merged["pooled_ac"] = merged[["AC_exomes", "AC_genomes"]].sum(axis=1, skipna=True)
    merged["pooled_an"] = merged[["AN_exomes", "AN_genomes"]].sum(axis=1, skipna=True)
    merged["pooled_af"] = np.where(merged["pooled_an"] > 0,
                                   merged["pooled_ac"] / merged["pooled_an"], np.nan)
    merged["combined_grpmax_af"] = merged[["AF_grpmax_exomes", "AF_grpmax_genomes"]].max(axis=1, skipna=True)

    merged["in_gnomad_exomes"]  = merged["AC_exomes"].notna()
    merged["in_gnomad_genomes"] = merged["AC_genomes"].notna()
    merged["in_gnomad"] = merged["in_gnomad_exomes"] | merged["in_gnomad_genomes"]
    merged["gnomad_source"] = np.select(
        [merged["in_gnomad_exomes"] & merged["in_gnomad_genomes"],
         merged["in_gnomad_exomes"],
         merged["in_gnomad_genomes"]],
        ["exomes+genomes", "exomes", "genomes"],
        default="")

    out = merged.rename(columns={
        "AF_exomes": "gnomad_exomes_af", "AF_genomes": "gnomad_genomes_af",
        "AC_exomes": "gnomad_exomes_ac", "AC_genomes": "gnomad_genomes_ac",
        "AN_exomes": "gnomad_exomes_an", "AN_genomes": "gnomad_genomes_an",
    })
    keep = ["pos", "vcf_ref", "vcf_alt", "in_gnomad", "in_gnomad_exomes",
            "in_gnomad_genomes", "gnomad_source",
            "gnomad_exomes_ac", "gnomad_exomes_an", "gnomad_exomes_af",
            "gnomad_genomes_ac", "gnomad_genomes_an", "gnomad_genomes_af",
            "pooled_ac", "pooled_an", "pooled_af", "combined_grpmax_af",
            "grpmax_pop_exomes", "grpmax_pop_genomes"]
    return out[keep]


# ── Shared helpers (same as sge_clinvar_intersect_plot_gaussian.py) ─────────

def find_deseq2_files(deseq2_dir: str) -> dict:
    """Targeton_ID = everything before the first underscore in the filename."""
    mapping = {}
    for fname in sorted(os.listdir(deseq2_dir)):
        if not fname.endswith(".tsv"):
            continue
        tid = fname.split("_")[0]
        if tid:
            mapping[tid] = os.path.join(deseq2_dir, fname)
    return mapping


def load_meta_consequences(meta_dir: str, targeton_id: str) -> pd.DataFrame:
    """Identical to the ClinVar script's loader: gives vcf_pos/vcf_ref/vcf_alt
    (needed for the merge key) plus HGVSc/HGVSp/consequence annotation."""
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
                        df = pd.read_csv(path, sep="\t", usecols=base_cols + ["Summary_Plot"])
                    except ValueError:
                        df = pd.read_csv(path, sep="\t", usecols=base_cols)
                    df["oligo_name"] = df["unique_oligo_name"].str.replace(
                        r"_\d+_\d+$", "", regex=True)
                    df = df.drop_duplicates(subset=["oligo_name"])
                    return df
            except Exception:
                continue
    return pd.DataFrame()


def get_depletion_column(df: pd.DataFrame) -> str:
    if "anchor_tier" in df.columns:
        return "anchor_tier"
    if "GMM_status" in df.columns:
        return "GMM_status"
    return "stat_pos_raw"


def extract_variant_key(oligo_name: str):
    """Same filter as the ClinVar script: drop PAM/ref-seq control oligos."""
    import re
    m = re.search(r"(chrX:\d+_[^_]+(?:>[A-Z]+)?)", str(oligo_name))
    return m.group(1) if m else None


SUMMARY_COLUMNS = [
    "sequence", "oligo_name", "position", "pos_adj_log2FoldChange_raw",
    "pos_adj_score_raw", "pos_adj_pval_raw", "pos_adj_fdr_raw", "stat_pos_raw",
    "anchor_tier", "anchor_call",
    "HGVSc", "HGVSp", "Protein_position", "Consequence", "Summary_Consequence",
    "Summary_Plot", "vcf_ref", "vcf_alt",
    "in_gnomad", "in_gnomad_exomes", "in_gnomad_genomes", "gnomad_source",
    "gnomad_exomes_ac", "gnomad_exomes_an", "gnomad_exomes_af",
    "gnomad_genomes_ac", "gnomad_genomes_an", "gnomad_genomes_af",
    "pooled_ac", "pooled_an", "pooled_af", "combined_grpmax_af",
    "grpmax_pop_exomes", "grpmax_pop_genomes",
]


# ── Stats ─────────────────────────────────────────────────────────────────────

def run_absence_stats(summary_df: pd.DataFrame, depletion_col: str, outdir: str):
    """
    Per-tier gnomAD match rate + Fisher's exact test comparing
    (depleted or enriched) vs "no impact" for gnomAD presence.
    Also Mann-Whitney U on pooled_af between tiers, for gnomAD-matched
    variants only. Writes gnomad_absence_stats.tsv, returns nothing.
    """
    df = summary_df.dropna(subset=[depletion_col]).copy()
    df[depletion_col] = df[depletion_col].astype(str)

    lines = []
    lines.append("=== gnomAD match rate by tier ===")
    rate_rows = []
    for tier in TIER_ORDER:
        sub = df[df[depletion_col] == tier]
        n = len(sub)
        n_matched = int(sub["in_gnomad"].sum())
        pct = 100 * n_matched / n if n else float("nan")
        rate_rows.append({"tier": tier, "n_variants": n,
                          "n_in_gnomad": n_matched, "pct_in_gnomad": round(pct, 2)})
        lines.append(f"  {tier:<20s}  n={n:<5d}  in_gnomad={n_matched:<5d}  ({pct:.2f}%)")
    rate_df = pd.DataFrame(rate_rows)

    lines.append("")
    lines.append("=== Fisher's exact test: (weakly+strongly depleting, enriched) vs "
                 "no impact, gnomAD presence ===")
    fisher_rows = []
    no_impact = df[df[depletion_col] == "no impact"]
    for group_name, group_tiers in [
        ("weakly depleting",   ["weakly depleting"]),
        ("strongly depleting", ["strongly depleting"]),
        ("any depleting",      ["weakly depleting", "strongly depleting"]),
        ("enriched",           ["enriched"]),
    ]:
        grp = df[df[depletion_col].isin(group_tiers)]
        if grp.empty or no_impact.empty:
            lines.append(f"  {group_name}: skipped (no rows in one group)")
            continue
        table = [
            [int(grp["in_gnomad"].sum()), int((~grp["in_gnomad"]).sum())],
            [int(no_impact["in_gnomad"].sum()), int((~no_impact["in_gnomad"]).sum())],
        ]
        odds_ratio, p = fisher_exact(table)
        lines.append(f"  {group_name:<20s} vs no impact: OR={odds_ratio:.3f}  p={p:.3e}  "
                     f"(n={len(grp)} vs n={len(no_impact)})")
        fisher_rows.append({"comparison": f"{group_name} vs no impact",
                            "odds_ratio": odds_ratio, "p_value": p,
                            "n_group": len(grp), "n_no_impact": len(no_impact)})
    fisher_df = pd.DataFrame(fisher_rows)

    lines.append("")
    lines.append("=== Mann-Whitney U: pooled_af (gnomAD-matched variants only), "
                 "each tier vs no impact ===")
    mw_rows = []
    ni_af = no_impact.loc[no_impact["in_gnomad"], "pooled_af"].dropna()
    for tier in ["weakly depleting", "strongly depleting", "enriched"]:
        tier_af = df.loc[(df[depletion_col] == tier) & df["in_gnomad"], "pooled_af"].dropna()
        if len(tier_af) < 2 or len(ni_af) < 2:
            lines.append(f"  {tier}: skipped (too few gnomAD-matched variants with AF)")
            continue
        stat, p = mannwhitneyu(tier_af, ni_af, alternative="two-sided")
        lines.append(f"  {tier:<20s} vs no impact: U={stat:.1f}  p={p:.3e}  "
                     f"(n={len(tier_af)} vs n={len(ni_af)})")
        mw_rows.append({"tier": tier, "U_stat": stat, "p_value": p,
                        "n_tier": len(tier_af), "n_no_impact": len(ni_af)})
    mw_df = pd.DataFrame(mw_rows)

    print("\n".join(lines))
    stats_path = os.path.join(outdir, "gnomad_absence_stats.tsv")
    with open(stats_path, "w") as fh:
        fh.write("\n".join(lines) + "\n\n")
        fh.write("# match_rate_by_tier\n")
        rate_df.to_csv(fh, sep="\t", index=False)
        fh.write("\n# fisher_exact_gnomad_presence\n")
        fisher_df.to_csv(fh, sep="\t", index=False)
        fh.write("\n# mannwhitney_pooled_af\n")
        mw_df.to_csv(fh, sep="\t", index=False)
    print(f"\nWrote stats: {stats_path}")


# ── Plots ─────────────────────────────────────────────────────────────────────

def plot_match_rate_by_tier(summary_df: pd.DataFrame, depletion_col: str, outpath: str):
    df = summary_df.dropna(subset=[depletion_col])
    tiers = [t for t in TIER_ORDER if (df[depletion_col] == t).any()]
    rates, ns = [], []
    for t in tiers:
        sub = df[df[depletion_col] == t]
        rates.append(100 * sub["in_gnomad"].mean())
        ns.append(len(sub))

    fig, ax = plt.subplots(figsize=(7, 5))
    fig.patch.set_facecolor("#fafafa")
    ax.set_facecolor("#fafafa")
    bars = ax.bar(tiers, rates, color=[TIER_COLOURS.get(t, "#999") for t in tiers])
    for bar, n in zip(bars, ns):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
               f"n={n}", ha="center", fontsize=8)
    ax.set_ylabel("% of variants observed in gnomAD (exomes and/or genomes)")
    ax.set_title("gnomAD match rate by SGE anchor tier", fontsize=11, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"Saved: {outpath}")


def plot_af_by_tier(summary_df: pd.DataFrame, depletion_col: str, outpath: str):
    df = summary_df.dropna(subset=[depletion_col])
    df = df[df["in_gnomad"] & df["pooled_af"].notna() & (df["pooled_af"] > 0)]
    if df.empty:
        print("  No gnomAD-matched variants with a positive pooled AF — skipping AF boxplot")
        return
    tiers = [t for t in TIER_ORDER if (df[depletion_col] == t).any()]
    data = [np.log10(df.loc[df[depletion_col] == t, "pooled_af"]) for t in tiers]

    fig, ax = plt.subplots(figsize=(7, 5))
    fig.patch.set_facecolor("#fafafa")
    ax.set_facecolor("#fafafa")
    bp = ax.boxplot(data, labels=tiers, patch_artist=True, showfliers=True)
    for patch, t in zip(bp["boxes"], tiers):
        patch.set_facecolor(TIER_COLOURS.get(t, "#999"))
        patch.set_alpha(0.6)
    for i, t in enumerate(tiers):
        n = (df[depletion_col] == t).sum()
        ax.text(i + 1, ax.get_ylim()[1], f"n={n}", ha="center", fontsize=8, va="bottom")
    ax.set_ylabel("log10(pooled gnomAD allele frequency)")
    ax.set_title("gnomAD allele frequency by SGE anchor tier\n(gnomAD-matched variants only)",
                 fontsize=11, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"Saved: {outpath}")


# ── Main ──────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--deseq2_dir", required=True,
                        help="Directory of per-targeton DESeq2/anchor TSVs "
                             "(filenames must start with Targeton_ID)")
    parser.add_argument("--meta_dir", required=True,
                        help="Directory with *_meta_consequences.tsv files")
    parser.add_argument("--gnomad_vcf_dir", required=True,
                        help="Directory with <Targeton_ID>_gnomad_exomes.vcf.gz "
                             "and <Targeton_ID>_gnomad_genomes.vcf.gz "
                             "(from extract_gnomad_vcfs.sh)")
    parser.add_argument("--regions", required=True,
                        help="targeton_regions.tsv from extract_targeton_regions.py")
    parser.add_argument("--outdir", default="gnomad_intersect_results",
                        help="Output directory (default: gnomad_intersect_results/)")
    parser.add_argument("--ac_field", default="AC")
    parser.add_argument("--an_field", default="AN")
    parser.add_argument("--af_field", default="AF")
    parser.add_argument("--grpmax_af_field", default="AF_grpmax",
                        help="gnomAD v4.1 popmax-replacement field name. "
                             "VERIFY against your VCF header before trusting output.")
    parser.add_argument("--grpmax_pop_field", default="grpmax")
    parser.add_argument("--include_non_pass", action="store_true",
                        help="By default only FILTER=PASS gnomAD records are kept. "
                             "Set this to include non-PASS records too.")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    pass_only = not args.include_non_pass

    print("Scanning DESeq2 directory...")
    deseq2_files = find_deseq2_files(args.deseq2_dir)
    print(f"  Found {len(deseq2_files)} targeton result file(s)")

    regions = pd.read_csv(args.regions, sep="\t")
    print(f"  {len(regions)} targetons in regions file")

    summary_rows = []
    depletion_col_seen = None

    for _, reg in regions.iterrows():
        tid = reg["Targeton_ID"]
        print(f"\n→ {tid}")

        if tid not in deseq2_files:
            print(f"  No DESeq2 result file found for {tid} — skipping")
            continue
        deseq_t = pd.read_csv(deseq2_files[tid], sep="\t")
        deseq_t["var_key"] = deseq_t["oligo_name"].apply(extract_variant_key)
        deseq_t = deseq_t.dropna(subset=["var_key"])
        if deseq_t.empty:
            print("  No variant rows after dropping controls — skipping")
            continue
        print(f"  {len(deseq_t)} DESeq2 variants")

        depletion_col = get_depletion_column(deseq_t)
        depletion_col_seen = depletion_col

        meta = load_meta_consequences(args.meta_dir, tid)
        if meta.empty:
            print(f"  WARNING: no meta_consequences file found for {tid} — skipping")
            continue
        meta_join_cols = ["oligo_name", "vcf_pos", "vcf_ref", "vcf_alt",
                          "HGVSc", "HGVSp", "Protein_position",
                          "Consequence", "Summary_Consequence"]
        if "Summary_Plot" in meta.columns:
            meta_join_cols.append("Summary_Plot")
        deseq_t = deseq_t.merge(meta[meta_join_cols], on="oligo_name", how="left")

        exomes_path  = os.path.join(args.gnomad_vcf_dir, f"{tid}_gnomad_exomes.vcf.gz")
        genomes_path = os.path.join(args.gnomad_vcf_dir, f"{tid}_gnomad_genomes.vcf.gz")

        exomes  = read_gnomad_vcf(exomes_path,  "exomes",
                                  args.ac_field, args.an_field, args.af_field,
                                  args.grpmax_af_field, args.grpmax_pop_field, pass_only)
        genomes = read_gnomad_vcf(genomes_path, "genomes",
                                  args.ac_field, args.an_field, args.af_field,
                                  args.grpmax_af_field, args.grpmax_pop_field, pass_only)
        print(f"  gnomAD exomes: {len(exomes)} PASS record(s)  |  "
              f"genomes: {len(genomes)} PASS record(s)"
              f"{'  (including non-PASS)' if not pass_only else ''}")

        gnomad = combine_gnomad_sources(exomes, genomes)

        if not gnomad.empty:
            gnomad_merge = gnomad.rename(columns={"pos": "vcf_pos"})
            deseq_t = deseq_t.merge(gnomad_merge, on=["vcf_pos", "vcf_ref", "vcf_alt"], how="left")
        for col in ["in_gnomad", "in_gnomad_exomes", "in_gnomad_genomes"]:
            if col not in deseq_t.columns:
                deseq_t[col] = False
            deseq_t[col] = deseq_t[col].fillna(False).astype(bool)
        if "gnomad_source" not in deseq_t.columns:
            deseq_t["gnomad_source"] = ""
        deseq_t["gnomad_source"] = deseq_t["gnomad_source"].fillna("")

        n_gnomad = int(deseq_t["in_gnomad"].sum())
        print(f"  {n_gnomad} / {len(deseq_t)} DESeq2 variants matched to gnomAD")

        # ── Per-targeton annotated TSV ───────────────────────────────────────
        out_cols = [c for c in deseq_t.columns if c != "var_key"]
        tsv_path = os.path.join(args.outdir, f"{tid}_gnomad_annotated.tsv")
        deseq_t[out_cols].to_csv(tsv_path, sep="\t", index=False)
        print(f"  Annotated TSV: {tsv_path}")

        # ── Collect for cross-targeton summary (ALL variants, not just matched) ──
        exon = TARGETON_EXON.get(tid)
        if exon is None:
            print(f"  WARNING: no exon mapping found for targeton '{tid}' — "
                  f"'exon' will be blank in the summary TSV for these rows")
        present_cols = [c for c in SUMMARY_COLUMNS if c in deseq_t.columns]
        row_block = deseq_t[present_cols].copy()
        row_block.insert(0, "Targeton_ID", tid)
        row_block.insert(1, "exon", exon)
        summary_rows.append(row_block)

    if not summary_rows:
        sys.exit("\nNo targetons produced any rows — nothing to summarise. Exiting.")

    summary_df = pd.concat(summary_rows, ignore_index=True)
    summary_path = os.path.join(args.outdir, "sge_gnomad_summary.tsv")
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"\nWrote cross-targeton summary ({len(summary_df)} variants): {summary_path}")

    depletion_col = depletion_col_seen or "anchor_tier"
    if depletion_col not in summary_df.columns:
        sys.exit(f"\nDepletion column '{depletion_col}' not present in summary — "
                 f"cannot run stats/plots. Check your DESeq2/anchor input files.")

    print("\nRunning absence/enrichment stats...")
    run_absence_stats(summary_df, depletion_col, args.outdir)

    print("\nGenerating plots...")
    plot_match_rate_by_tier(summary_df, depletion_col,
                            os.path.join(args.outdir, "gnomad_match_rate_by_tier.png"))
    plot_af_by_tier(summary_df, depletion_col,
                    os.path.join(args.outdir, "gnomad_af_by_tier.png"))

    print(f"\nDone. Results in: {os.path.abspath(args.outdir)}/")


if __name__ == "__main__":
    main()
