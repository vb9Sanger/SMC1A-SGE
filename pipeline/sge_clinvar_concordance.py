#!/usr/bin/env python3
"""
sge_clinvar_concordance.py

Quantitatively compares any number of SGE classification approaches (e.g.
gene-wide anchor / per-targeton anchor / shrinkage anchor / the original
un-anchored per-exon GMM) against ClinVar, to answer: which approach's
"depleted" and "no impact" calls best agree with ClinVar's "Pathogenic/
Likely pathogenic" and "Benign/Likely benign" classifications?

Variant matching reuses the same logic as sge_clinvar_intersect_plot*.py:
chrX:POS_REF>ALT keys extracted directly from oligo_name (no separate
meta_consequences join needed for this).

For each approach, computes:
  - sensitivity, STRATIFIED BY DISEASE CONDITION, not just pooled. This
    matters specifically for SMC1A: CdLS-causing variants act via a
    dominant-negative/altered-function mechanism (the protein is made and
    incorporated into cohesin, just disrupts its function) rather than
    loss-of-function, while DEE-causing variants are predominantly true
    loss-of-function. A depletion-style functional assay is only expected
    to reliably catch the LOF mechanism -- so a CdLS pathogenic variant
    showing "no impact" is not evidence the classification is wrong, it's
    the expected result of an assay that isn't built to detect a
    dominant-negative effect. Pooling all P/LP variants together would
    dilute (or contradict) the real signal. Condition is parsed from
    ClinVar's CLNDN field using the same term lists as
    sge_clinvar_intersect_plot*.py: CdLS only / DEE only / Both / Other
    SMC1A / Unclassified.
      sensitivity_plp_all         -- pooled (kept for reference, but
                                      treat as the least meaningful number
                                      here given the above)
      sensitivity_plp_dee_only    -- the PRIMARY validation metric: DEE
                                      variants are expected to deplete
      sensitivity_plp_cdls_only   -- expected to be lower / uninformative
                                      by design, not a failure signal
      sensitivity_plp_both / _other_smc1a / _unclassified -- for completeness
  - specificity   = fraction of ClinVar B/LB variants called "no impact"
                    (pooled -- benign variants aren't subject to the same
                    mechanism confound, so not stratified)
  - balanced accuracy = mean(sensitivity_plp_all, specificity)
  - AUC, similarly stratified:
      auc_all             -- all P/LP vs all B/LB (as before)
      auc_dee_restricted  -- DEE-only P/LP vs all B/LB (the more
                              mechanistically appropriate validation AUC)
      auc_cdls_restricted -- CdLS-only P/LP vs all B/LB (expected near
                              0.5 by design, included for contrast)
    Uses each approach's own continuous "probability of depletion" score
    (anchor_post_lof for the anchor pipelines; reconstructed from
    GMM_prob_cluster1/2 + GMM_cluster + GMM_label for the original
    un-anchored per-exon GMM), computed via the rank-sum/Mann-Whitney
    formula (no sklearn dependency). This is what makes a POOLED,
    cross-exon AUC meaningful even when different targetons/approaches
    define "how depleted is depleted" on different local scales.
  - VUS/other-classified ClinVar variants are counted but excluded from
    the above (standard practice -- only unambiguous ground truth used).

Outputs (in --outdir):
  concordance_summary.tsv         one row per approach: the core metrics,
                                   sensitivity/AUC stratified by condition
  tier_breakdown_by_approach.tsv  long-format: approach x clinvar_call x
                                   condition x tier value x count (finer
                                   detail than the summary's binary
                                   sensitivity, e.g. do DEE P/LP variants
                                   skew toward "strongly" vs "weakly"
                                   depleting specifically)
  <approach>_matched_variants.tsv one row per ClinVar-matched SGE variant
                                   for that approach, incl. its condition,
                                   for manual inspection
  concordance_comparison.png      grouped bar chart: specificity /
                                   balanced accuracy / sensitivity
                                   (DEE-only vs CdLS-only vs all) / AUC
                                   (DEE-restricted vs CdLS-restricted vs
                                   all), side by side across approaches

Usage:
  python sge_clinvar_concordance.py \\
    --clinvar_vcf SMC1A_clinvar_aggregate.vcf.gz \\
    --approach "GeneWide=GMM_geneWide_results/*.tsv:anchor_tier" \\
    --approach "PerTargeton=GMM_perTargeton_results/*.tsv:anchor_tier" \\
    --approach "Shrinkage=GMM_shrinkage_results/*.tsv:anchor_tier" \\
    --approach "OriginalGMM=GMM_corrected/*_with_GMM.tsv:GMM_status" \\
    --outdir clinvar_concordance

  Or with per-targeton ClinVar VCFs instead of one aggregate file:
    --vcf_dir per_targeton_clinvar_vcfs/

Each --approach is repeatable: NAME=GLOB_PATTERN:TIER_COLUMN
  (glob pattern quoted so the shell doesn't expand it; split on the LAST
  ':' so Unix paths -- which don't normally contain literal colons --
  are safe)

Requirements: pandas, numpy, scipy, matplotlib
"""

import argparse
import glob as globmod
import gzip
import os
import re
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy.stats import rankdata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── Disease condition classification (ported from sge_clinvar_intersect_plot_V2-3.py) --
# Used to stratify sensitivity/AUC, since CdLS and DEE variants act through
# different mechanisms in SMC1A (see module docstring).

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
OTHER_SMC1A_TERMS = [
    "smc1a-related",
    "x-linked_complex_neurodevelopmental",
    "cohesinopathy",
]


def classify_condition(clndn):
    """Map a raw CLNDN string to one of five condition categories."""
    if not clndn or clndn in (".", ""):
        return "Not in ClinVar"
    s = str(clndn).lower()
    has_cdls = any(t in s for t in CDLS_TERMS)
    has_dee = any(t in s for t in DEE_TERMS)
    if has_cdls and has_dee:
        return "Both"
    if has_cdls:
        return "CdLS only"
    if has_dee:
        return "DEE only"
    if any(t in s for t in OTHER_SMC1A_TERMS):
        return "Other SMC1A"
    return "Unclassified"


# Conditions to report sensitivity for, in a stable order (used for both the
# summary table's column suffixes and the plot's bar ordering).
CONDITION_STRATA = ["DEE only", "CdLS only", "Both", "Other SMC1A", "Unclassified"]
CONDITION_SUFFIX = {
    "DEE only": "dee_only", "CdLS only": "cdls_only", "Both": "both",
    "Other SMC1A": "other_smc1a", "Unclassified": "unclassified",
}


# ── ClinVar VCF parsing (adapted from sge_clinvar_intersect_plot_V2-3.py) ────

def extract_variant_key(oligo_name):
    """Extract chrX:POS_VARIANT key from oligo_name for joining."""
    m = re.search(r"(chrX:\d+_[^_]+(?:>[A-Z]+)?)", str(oligo_name))
    return m.group(1) if m else None


def normalise_clnsig(raw):
    if not raw:
        return ""
    return raw.replace("_", " ")


def is_plp(clnsig_norm):
    s = clnsig_norm.lower()
    return "pathogenic" in s and "benign" not in s


def is_blb(clnsig_norm):
    s = clnsig_norm.lower()
    return "benign" in s and "pathogenic" not in s


def read_clinvar_vcf(vcf_path):
    """
    Parse a ClinVar VCF (bgzipped or plain, aggregate or per-targeton).
    Returns DataFrame with: chrom, pos, clinvar_id, vcf_ref, vcf_alt,
    clnsig_raw, clndn_raw, clnsig_norm, var_key, clinvar_call.
    """
    rows = []
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
            info = parts[7]
            clnsig = clndn = ""
            for field in info.split(";"):
                if field.startswith("CLNSIG="):
                    clnsig = field[7:]
                elif field.startswith("CLNDN="):
                    clndn = field[6:]
            rows.append({
                "chrom": chrom, "pos": pos, "clinvar_id": vid,
                "vcf_ref": ref, "vcf_alt": alt,
                "clnsig_raw": clnsig, "clndn_raw": clndn,
            })
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df["clnsig_norm"] = df["clnsig_raw"].apply(normalise_clnsig)

    def make_var_key(row):
        chrom_str = "chrX" if row["chrom"] in ("X", "chrX") else f"chr{row['chrom']}"
        ref, alt = row["vcf_ref"], row["vcf_alt"]
        if len(ref) == len(alt) == 1:
            return f"{chrom_str}:{row['pos']}_{ref}>{alt}"
        elif len(alt) < len(ref):
            n_del = len(ref) - len(alt)
            return f"{chrom_str}:{row['pos'] + 1}_{n_del}del"
        else:
            n_ins = len(alt) - len(ref)
            return f"{chrom_str}:{row['pos'] + 1}_{n_ins}ins"

    df["var_key"] = df.apply(make_var_key, axis=1)
    df["clinvar_call"] = np.where(
        df["clnsig_norm"].apply(is_plp), "P/LP",
        np.where(df["clnsig_norm"].apply(is_blb), "B/LB", "VUS/other")
    )
    df["condition"] = df["clndn_raw"].apply(classify_condition)
    return df


def load_clinvar(clinvar_vcf, vcf_dir):
    if clinvar_vcf:
        cv = read_clinvar_vcf(clinvar_vcf)
        print(f"Loaded {len(cv)} ClinVar variant(s) from aggregate VCF: {clinvar_vcf}")
    elif vcf_dir:
        dfs = []
        for fname in sorted(os.listdir(vcf_dir)):
            if fname.endswith(".vcf") or fname.endswith(".vcf.gz"):
                d = read_clinvar_vcf(os.path.join(vcf_dir, fname))
                if not d.empty:
                    dfs.append(d)
        cv = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
        print(f"Loaded {len(cv)} ClinVar variant(s) pooled from {len(dfs)} file(s) in {vcf_dir}")
    else:
        raise SystemExit("Must supply either --clinvar_vcf or --vcf_dir")
    if cv.empty:
        raise SystemExit("No ClinVar variants loaded -- check VCF path/contents.")
    n_before = len(cv)
    cv = cv.drop_duplicates(subset=["var_key"])
    if len(cv) < n_before:
        print(f"  ({n_before - len(cv)} duplicate var_key(s) dropped, e.g. overlapping targeton VCF windows)")
    n_plp = (cv["clinvar_call"] == "P/LP").sum()
    n_blb = (cv["clinvar_call"] == "B/LB").sum()
    n_vus = (cv["clinvar_call"] == "VUS/other").sum()
    print(f"  ClinVar breakdown: P/LP={n_plp}  B/LB={n_blb}  VUS/other={n_vus}")
    plp_by_cond = cv[cv["clinvar_call"] == "P/LP"]["condition"].value_counts()
    print(f"  P/LP by condition: {plp_by_cond.to_dict()}\n")
    return cv


# ── Depletion-call normalisation across different pipelines' vocabularies ────

def normalise_tier(value):
    """Map a tier/status string to one of: depleted / no impact / enriched / None."""
    if pd.isna(value):
        return None
    s = str(value).strip().lower()
    if "deplet" in s:
        return "depleted"
    if s == "no impact":
        return "no impact"
    if s == "enriched":
        return "enriched"
    return None


# ── Reconstruct a comparable depletion-probability score for the original ---
# ── un-anchored per-exon GMM output (GMM_prob_cluster1/2), so it can be   ---
# ── compared on equal footing with anchor_post_lof.                      ---

def add_gmm_prob_depleted(df):
    """
    Adds 'prob_depleted' = posterior probability of whichever GMM_cluster
    was labelled 'depleted' -- inferred PER SOURCE FILE (since the
    un-anchored GMM was fit independently per targeton, cluster 1 vs 2
    isn't consistently "depleted" across files).
    """
    needed = {"GMM_cluster", "GMM_label", "GMM_prob_cluster1", "GMM_prob_cluster2", "__source_file"}
    if not needed.issubset(df.columns):
        df["prob_depleted"] = np.nan
        return df

    def per_file(sub):
        depleted_rows = sub[sub["GMM_label"] == "depleted"]
        if depleted_rows.empty:
            sub = sub.copy()
            sub["prob_depleted"] = np.nan
            return sub
        depleted_cluster = depleted_rows["GMM_cluster"].mode()
        if depleted_cluster.empty:
            sub = sub.copy()
            sub["prob_depleted"] = np.nan
            return sub
        depleted_cluster = int(depleted_cluster.iloc[0])
        prob_col = f"GMM_prob_cluster{depleted_cluster}"
        sub = sub.copy()
        sub["prob_depleted"] = sub[prob_col] if prob_col in sub.columns else np.nan
        return sub

    return df.groupby("__source_file", group_keys=False).apply(per_file)


# ── AUC via rank-sum / Mann-Whitney U (no sklearn dependency) ────────────────

def compute_auc(scores, labels_binary):
    """
    scores: continuous scores, higher = more like the positive class (P/LP)
    labels_binary: 0/1 array, 1 = P/LP
    """
    scores = np.asarray(scores, dtype=float)
    labels_binary = np.asarray(labels_binary)
    mask = ~np.isnan(scores)
    scores = scores[mask]
    labels_binary = labels_binary[mask]
    n_pos = labels_binary.sum()
    n_neg = len(labels_binary) - n_pos
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = rankdata(scores)
    sum_ranks_pos = ranks[labels_binary == 1].sum()
    return (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


# ── Per-approach analysis ─────────────────────────────────────────────────────

def analyse_approach(name, glob_pattern, tier_col, clinvar_df, oligo_col):
    files = sorted(globmod.glob(glob_pattern))
    if not files:
        print(f"[{name}] WARNING: no files matched glob '{glob_pattern}' -- skipping")
        return None, None, None

    dfs = []
    for f in files:
        d = pd.read_csv(f, sep="\t", low_memory=False)
        d["__source_file"] = os.path.basename(f)
        dfs.append(d)
    all_df = pd.concat(dfs, ignore_index=True)

    if tier_col not in all_df.columns:
        print(f"[{name}] WARNING: tier column '{tier_col}' not found in these {len(files)} file(s) -- skipping")
        return None, None, None
    if oligo_col not in all_df.columns:
        print(f"[{name}] WARNING: oligo column '{oligo_col}' not found -- cannot match to ClinVar -- skipping")
        return None, None, None

    all_df["var_key"] = all_df[oligo_col].apply(extract_variant_key)
    n_before = len(all_df)
    all_df = all_df.dropna(subset=["var_key"])
    print(f"[{name}] {len(files)} file(s), {n_before} row(s), {len(all_df)} with an extractable variant key")

    if tier_col == "GMM_status":
        all_df = add_gmm_prob_depleted(all_df)
        score_col = "prob_depleted"
    elif "anchor_post_lof" in all_df.columns:
        score_col = "anchor_post_lof"
    else:
        score_col = None
        print(f"[{name}] NOTE: no recognised probability-of-depletion column found -- AUC will be NA for this approach")

    keep_cols = ["var_key", "clinvar_call", "condition", "clnsig_norm", "clnsig_raw", "clndn_raw", "clinvar_id"]
    keep_cols = [c for c in keep_cols if c in clinvar_df.columns]
    merged = all_df.merge(clinvar_df[keep_cols], on="var_key", how="inner")

    if merged.empty:
        print(f"[{name}] WARNING: 0 variants matched to ClinVar -- check var_key construction/--oligo_col")
        return None, None, None

    merged["tier_norm"] = merged[tier_col].apply(normalise_tier)
    if "condition" not in merged.columns:
        merged["condition"] = "Unclassified"

    plp = merged[merged["clinvar_call"] == "P/LP"]
    blb = merged[merged["clinvar_call"] == "B/LB"]
    vus = merged[merged["clinvar_call"] == "VUS/other"]

    n_plp, n_blb, n_vus = len(plp), len(blb), len(vus)
    sensitivity_all = (plp["tier_norm"] == "depleted").mean() if n_plp > 0 else np.nan
    specificity = (blb["tier_norm"] == "no impact").mean() if n_blb > 0 else np.nan
    balanced_acc = np.nanmean([sensitivity_all, specificity]) if (n_plp > 0 or n_blb > 0) else np.nan

    def auc_for(plp_subset):
        auc_df = pd.concat([plp_subset, blb])
        if score_col is None:
            return np.nan
        auc_df = auc_df.dropna(subset=[score_col])
        if auc_df.empty:
            return np.nan
        labels = (auc_df["clinvar_call"] == "P/LP").astype(int).values
        return compute_auc(auc_df[score_col].values, labels)

    auc_all = auc_for(plp)

    result = {
        "approach": name,
        "n_files": len(files),
        "n_total_variants": len(all_df),
        "n_matched_clinvar": len(merged),
        "n_plp": n_plp, "n_blb": n_blb, "n_vus_other": n_vus,
        "sensitivity_plp_all": sensitivity_all,
        "specificity_blb_no_impact": specificity,
        "balanced_accuracy": balanced_acc,
        "auc_all": auc_all,
        "score_col_used": score_col if score_col is not None else "",
    }

    # Condition-stratified sensitivity and AUC -- see module docstring for
    # why this matters more than the pooled numbers above for SMC1A.
    for cond in CONDITION_STRATA:
        suffix = CONDITION_SUFFIX[cond]
        plp_cond = plp[plp["condition"] == cond]
        n_cond = len(plp_cond)
        sens_cond = (plp_cond["tier_norm"] == "depleted").mean() if n_cond > 0 else np.nan
        result[f"n_plp_{suffix}"] = n_cond
        result[f"sensitivity_plp_{suffix}"] = sens_cond
        result[f"auc_{suffix}_restricted"] = auc_for(plp_cond) if n_cond > 0 else np.nan

    # Long-format tier breakdown: approach x clinvar_call x condition x tier value x count
    breakdown_rows = []
    for label, sub in [("P/LP", plp), ("B/LB", blb), ("VUS/other", vus)]:
        counts = sub.groupby(["condition", tier_col], dropna=False).size()
        for (cond, tier_val), n in counts.items():
            if n == 0:
                continue
            breakdown_rows.append({"approach": name, "clinvar_call": label,
                                    "condition": cond, "tier_value": tier_val, "n": n})
    breakdown_df = pd.DataFrame(breakdown_rows)

    return result, merged, breakdown_df


# ── Comparison plot ───────────────────────────────────────────────────────────

def make_comparison_plot(summary_df, outpath):
    # Deliberately leads with the condition-stratified numbers, since those
    # are the ones that are actually interpretable for SMC1A (see module
    # docstring) -- the pooled "_all" versions are included at the end for
    # reference, not as the headline metrics.
    metrics = [
        "sensitivity_plp_dee_only", "sensitivity_plp_cdls_only",
        "specificity_blb_no_impact", "balanced_accuracy",
        "auc_dee_only_restricted", "auc_cdls_only_restricted",
        "sensitivity_plp_all", "auc_all",
    ]
    labels = [
        "Sensitivity\n(DEE-only P/LP)", "Sensitivity\n(CdLS-only P/LP)",
        "Specificity\n(B/LB -> no impact)", "Balanced accuracy\n(pooled)",
        "AUC\n(DEE-restricted)", "AUC\n(CdLS-restricted)",
        "Sensitivity\n(all P/LP, pooled)", "AUC\n(all, pooled)",
    ]
    metrics = [m for m in metrics if m in summary_df.columns]
    labels = labels[:len(metrics)]

    n_approach = len(summary_df)
    x = np.arange(len(metrics))
    width = 0.8 / max(n_approach, 1)

    fig, ax = plt.subplots(figsize=(13, 6))
    for i, (_, row) in enumerate(summary_df.iterrows()):
        vals = [row[m] for m in metrics]
        ax.bar(x + i * width, vals, width, label=row["approach"])

    ax.axhline(0.5, color="grey", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axvline(1.5, color="grey", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.axvline(3.5, color="grey", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.axvline(5.5, color="grey", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.set_xticks(x + width * (n_approach - 1) / 2)
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Value")
    ax.set_title("SGE classification vs. ClinVar: concordance by approach\n"
                 "(DEE-only sensitivity/AUC is the primary validation metric for SMC1A -- CdLS acts via a\n"
                 "dominant-negative mechanism a depletion assay isn't expected to catch; pooled numbers shown for reference only)",
                 fontsize=9.5)
    ax.legend(fontsize=8, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


# ── Main ──────────────────────────────────────────────────────────────────────

def parse_approach_arg(raw):
    if "=" not in raw:
        raise argparse.ArgumentTypeError(f"--approach must be NAME=GLOB:TIER_COL, got: {raw}")
    name, rest = raw.split("=", 1)
    if ":" not in rest:
        raise argparse.ArgumentTypeError(f"--approach must be NAME=GLOB:TIER_COL, got: {raw}")
    glob_pattern, tier_col = rest.rsplit(":", 1)
    return name, glob_pattern, tier_col


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clinvar_vcf", default=None, help="Single aggregate ClinVar VCF (bgzipped or plain)")
    p.add_argument("--vcf_dir", default=None, help="Directory of per-targeton ClinVar VCFs (alternative to --clinvar_vcf)")
    p.add_argument("--approach", action="append", required=True,
                   help="Repeatable: NAME=GLOB_PATTERN:TIER_COLUMN")
    p.add_argument("--oligo_col", default="oligo_name", help="Column to extract chrX:POS_REF>ALT from [default: %(default)s]")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    if not args.clinvar_vcf and not args.vcf_dir:
        sys.exit("ERROR: must supply --clinvar_vcf or --vcf_dir")
    os.makedirs(args.outdir, exist_ok=True)

    clinvar_df = load_clinvar(args.clinvar_vcf, args.vcf_dir)

    summary_rows = []
    breakdown_dfs = []

    for raw in args.approach:
        try:
            name, glob_pattern, tier_col = parse_approach_arg(raw)
        except argparse.ArgumentTypeError as e:
            sys.exit(f"ERROR: {e}")
        print(f"=== {name} ===")
        result, merged, breakdown_df = analyse_approach(name, glob_pattern, tier_col, clinvar_df, args.oligo_col)
        if result is None:
            continue
        summary_rows.append(result)
        breakdown_dfs.append(breakdown_df)

        detail_path = os.path.join(args.outdir, f"{name}_matched_variants.tsv")
        merged.to_csv(detail_path, sep="\t", index=False)

        def fmt(v):
            return "NA" if pd.isna(v) else f"{v:.3f}"

        print(f"  specificity={fmt(result['specificity_blb_no_impact'])}  "
              f"sensitivity(DEE-only)={fmt(result['sensitivity_plp_dee_only'])} (n={result['n_plp_dee_only']})  "
              f"sensitivity(CdLS-only)={fmt(result['sensitivity_plp_cdls_only'])} (n={result['n_plp_cdls_only']})  "
              f"sensitivity(all)={fmt(result['sensitivity_plp_all'])}")
        print(f"  AUC(DEE-restricted)={fmt(result['auc_dee_only_restricted'])}  "
              f"AUC(CdLS-restricted)={fmt(result['auc_cdls_only_restricted'])}  "
              f"AUC(all)={fmt(result['auc_all'])}")
        print(f"  Wrote: {detail_path}\n")

    if not summary_rows:
        sys.exit("No approach produced results -- check glob patterns and column names.")

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(args.outdir, "concordance_summary.tsv")
    summary_df.to_csv(summary_path, sep="\t", index=False)

    breakdown_all = pd.concat(breakdown_dfs, ignore_index=True) if breakdown_dfs else pd.DataFrame()
    breakdown_path = os.path.join(args.outdir, "tier_breakdown_by_approach.tsv")
    breakdown_all.to_csv(breakdown_path, sep="\t", index=False)

    plot_path = os.path.join(args.outdir, "concordance_comparison.png")
    make_comparison_plot(summary_df, plot_path)

    print("=== Summary ===")
    print(summary_df.to_string(index=False))
    print(f"\nWrote: {summary_path}")
    print(f"Wrote: {breakdown_path}")
    print(f"Wrote: {plot_path}")


if __name__ == "__main__":
    main()
