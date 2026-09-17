#!/usr/bin/env python3
"""
nondepleting_plp_vs_synonymous.py

Question: for ClinVar Pathogenic/Likely-pathogenic (P/LP) SMC1A variants that
the SGE assay calls anchor_tier == "no impact", is their functional-score
signal genuinely indistinguishable from synonymous (no-impact control)
variants -- or is there a subtle residual shift toward depletion that just
didn't cross the calling threshold (i.e. a power/sensitivity issue rather
than a true absence of effect)?

Approach:
  - Non-depleting P/LP set: pulled from the ClinVar intersection TSV
    (needs the clnsig/condition annotations, which only live there).
  - Synonymous reference set: pulled from the FULL per-targeton
    anchor-tier output files (all consequence classes, not just
    ClinVar-matched ones) -- i.e. every synonymous variant across all 22
    targetons, not just the ~210 that happen to also be ClinVar hits.
    This gives a much larger, more representative no-impact reference
    (thousands of variants instead of ~210).
  - Because raw log2FoldChange (LFC) scale/noise differs targeton to
    targeton, both groups are converted to a per-targeton normalised score:
        z = (LFC - anchor_mu_noimpact) / anchor_sd_noimpact
    i.e. "how many no-impact-control SDs away from that targeton's own
    no-impact anchor". z=0 means dead-center of the no-impact distribution;
    more negative = further toward the depleting direction.
  - Compares the two z-score distributions with a two-sided Mann-Whitney U
    test, for all non-depleting P/LP variants and for the missense-only
    subset.
  - Also compares anchor_post_lof (posterior probability of being
    LOF-like/depleted) between groups -- a direct read on how close each
    variant sat to the depleting/no-impact decision boundary (50% cutoff).

Inputs:
  --input       ClinVar intersection TSV (e.g. clinvar_variants_summary.tsv):
                 needs Summary_Plot, clnsig_norm, anchor_tier,
                 pos_adj_log2FoldChange_raw, anchor_post_lof, Targeton_ID,
                 HGVSp, condition.
  --anchor_dir  Directory of per-targeton anchor-tier output TSVs (the full
                 gaussian_shrinkage_classifier.R output, one file per
                 targeton, e.g. APDY_exon2_all_deseq2_results_condition_
                 Day15_vs_Day4.tsv), needs consequence,
                 pos_adj_log2FoldChange_raw, anchor_mu_noimpact,
                 anchor_sd_noimpact, anchor_post_lof. Targeton_ID is
                 inferred from each filename (stem, minus the
                 "_all_deseq2_results_condition_*" suffix).

Usage:
    python nondepleting_plp_vs_synonymous.py \\
        --input clinvar_variants_summary.tsv \\
        --anchor_dir per_targeton_results/ \\
        --outdir results/
"""

import argparse
import glob
import os
import sys

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PLP_LABELS = ["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="ClinVar intersection TSV (e.g. clinvar_variants_summary.tsv)")
    p.add_argument(
        "--anchor_dir", required=True,
        help="Directory of per-targeton anchor-tier output TSVs, used as the full synonymous reference set",
    )
    p.add_argument("--outdir", default=".", help="Output directory")
    return p.parse_args()


def infer_targeton_name(path):
    bn = os.path.basename(path).replace(".tsv", "")
    return bn.split("_all_deseq2_results_condition_")[0]


def load_full_synonymous_set(anchor_dir):
    """Every synonymous variant across all per-targeton anchor-tier output files."""
    files = sorted(glob.glob(os.path.join(anchor_dir, "*.tsv")))
    files = [f for f in files if not os.path.basename(f).startswith("gmm_shrinkage_anchor_summary")]
    if not files:
        sys.exit(f"ERROR: no *.tsv files found in --anchor_dir {anchor_dir}")

    required = ["consequence", "pos_adj_log2FoldChange_raw", "anchor_mu_noimpact", "anchor_sd_noimpact", "anchor_post_lof"]
    frames = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        missing = [c for c in required if c not in df.columns]
        if missing:
            sys.exit(f"ERROR: {f} is missing required column(s): {missing}")
        syn = df[df["consequence"] == "Synonymous_Variant"].copy()
        syn["Targeton_ID"] = infer_targeton_name(f)
        frames.append(syn)

    out = pd.concat(frames, ignore_index=True)
    print(f"Loaded {len(out)} synonymous variants across {len(files)} targeton files (full reference set).")
    return out


def add_z(df):
    df = df.copy()
    df["z_noimpact"] = (df["pos_adj_log2FoldChange_raw"] - df["anchor_mu_noimpact"]) / df["anchor_sd_noimpact"]
    return df


def mwu(a, b):
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    return u, p


def summarize(name, s):
    return {
        "group": name,
        "n": len(s),
        "median": round(s.median(), 3) if len(s) else np.nan,
        "q25": round(s.quantile(0.25), 3) if len(s) else np.nan,
        "q75": round(s.quantile(0.75), 3) if len(s) else np.nan,
    }


def make_plot(syn, nondep, nondep_missense, outpath):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    # Panel A: z-score (LFC normalised to each targeton's no-impact anchor)
    ax = axes[0]
    groups = [("Synonymous\n(reference, n=%d)" % len(syn), syn["z_noimpact"]),
              ("Non-depleting P/LP\n(all, n=%d)" % len(nondep), nondep["z_noimpact"]),
              ("Non-depleting P/LP\n(missense only, n=%d)" % len(nondep_missense), nondep_missense["z_noimpact"])]
    data = [g[1].dropna().values for g in groups]
    labels = [g[0] for g in groups]
    parts = ax.violinplot(data, showmedians=True, widths=0.7)
    for pc in parts["bodies"]:
        pc.set_facecolor("#4285F4")
        pc.set_alpha(0.35)
    for i, d in enumerate(data, start=1):
        jitter = np.random.default_rng(0).normal(0, 0.04, size=len(d))
        ax.scatter(np.full(len(d), i) + jitter, d, s=10, color="#EA4335", alpha=0.6, zorder=3)
    ax.axhline(0, color="grey", linestyle="--", linewidth=1)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("z (LFC relative to no-impact anchor)", fontsize=10)
    ax.set_title("LFC distribution: non-depleting P/LP vs synonymous", fontsize=12, pad=12)

    # Clip extreme outliers so a handful of points don't compress the rest of the
    # distribution. Percentile-based clipping breaks down if outliers exceed 1% of
    # a group, so use Tukey fences (Q1/Q3 +/- 3*IQR on the pooled data) instead,
    # and annotate how many points are clipped per group.
    pooled = np.concatenate(data)
    q1, q3 = np.percentile(pooled, [25, 75])
    iqr = q3 - q1
    ylim = (q1 - 3 * iqr, q3 + 3 * iqr)
    n_clipped = [int(((d < ylim[0]) | (d > ylim[1])).sum()) for d in data]
    if any(n_clipped):
        yspan = ylim[1] - ylim[0]
        ax.set_ylim(ylim[0] - 0.08 * yspan, ylim[1])  # small extra headroom at the bottom for the annotation
        for i, n in enumerate(n_clipped, start=1):
            if n:
                ax.annotate(f"{n} off-scale (z<{ylim[0]:.1f})", xy=(i, ylim[0] - 0.04 * yspan),
                            ha="center", va="center", fontsize=7, color="grey")

    # Panel B: posterior probability of LOF/depleted call (log scale)
    ax2 = axes[1]
    groups_p = [("Synonymous\n(n=%d)" % len(syn), syn["anchor_post_lof"]),
                ("Non-depleting P/LP\n(all, n=%d)" % len(nondep), nondep["anchor_post_lof"]),
                ("Non-depleting P/LP\n(missense only, n=%d)" % len(nondep_missense), nondep_missense["anchor_post_lof"])]
    data_p = [g[1].dropna().values for g in groups_p]
    labels_p = [g[0] for g in groups_p]
    for i, d in enumerate(data_p, start=1):
        jitter = np.random.default_rng(1).normal(0, 0.06, size=len(d))
        ax2.scatter(np.full(len(d), i) + jitter, d, s=10, color="#4285F4", alpha=0.6)
        ax2.scatter([i], [np.median(d)], marker="_", s=800, color="black", linewidths=2, zorder=4)
    ax2.axhline(0.5, color="#EA4335", linestyle="--", linewidth=1, label="depleting call threshold (0.5)")
    ax2.set_yscale("symlog", linthresh=0.01)
    ax2.set_xticks(range(1, len(labels_p) + 1))
    ax2.set_xticklabels(labels_p, fontsize=9)
    ax2.set_ylabel("posterior P(LOF/depleted-like)")
    ax2.set_title("How close to the calling threshold?")
    ax2.legend(fontsize=8, loc="upper left", frameon=False)

    for a in axes:
        for spine in ["top", "right"]:
            a.spines[spine].set_visible(False)

    plt.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return ylim


def main():
    args = parse_args()
    df = pd.read_csv(args.input, sep="\t")

    required = [
        "Summary_Plot", "clnsig_norm", "anchor_tier",
        "pos_adj_log2FoldChange_raw", "anchor_mu_noimpact", "anchor_sd_noimpact",
        "anchor_post_lof",
    ]
    for col in required:
        if col not in df.columns:
            sys.exit(f"ERROR: expected column '{col}' not found in {args.input}")

    os.makedirs(args.outdir, exist_ok=True)

    plp_mask = df["clnsig_norm"].isin(PLP_LABELS)
    nondep_plp = add_z(df[(df["anchor_tier"] == "no impact") & plp_mask])
    nondep_missense = nondep_plp[nondep_plp["Summary_Plot"] == "Missense_Variant"]

    syn_raw = load_full_synonymous_set(args.anchor_dir)
    syn = add_z(syn_raw)

    # --- Statistics ---
    u_z_all, p_z_all = mwu(nondep_plp["z_noimpact"], syn["z_noimpact"])
    u_z_mis, p_z_mis = mwu(nondep_missense["z_noimpact"], syn["z_noimpact"])
    u_post_all, p_post_all = mwu(nondep_plp["anchor_post_lof"], syn["anchor_post_lof"])
    u_post_mis, p_post_mis = mwu(nondep_missense["anchor_post_lof"], syn["anchor_post_lof"])

    stats_rows = [
        {"comparison": "z_noimpact: non-depleting P/LP (all) vs synonymous", "n1": len(nondep_plp), "n2": len(syn), "U": u_z_all, "p_value": p_z_all},
        {"comparison": "z_noimpact: non-depleting P/LP (missense only) vs synonymous", "n1": len(nondep_missense), "n2": len(syn), "U": u_z_mis, "p_value": p_z_mis},
        {"comparison": "anchor_post_lof: non-depleting P/LP (all) vs synonymous", "n1": len(nondep_plp), "n2": len(syn), "U": u_post_all, "p_value": p_post_all},
        {"comparison": "anchor_post_lof: non-depleting P/LP (missense only) vs synonymous", "n1": len(nondep_missense), "n2": len(syn), "U": u_post_mis, "p_value": p_post_mis},
    ]
    stats_df = pd.DataFrame(stats_rows)

    summary_rows = [
        summarize("synonymous (z)", syn["z_noimpact"]),
        summarize("non-depleting P/LP - all (z)", nondep_plp["z_noimpact"]),
        summarize("non-depleting P/LP - missense only (z)", nondep_missense["z_noimpact"]),
        summarize("synonymous (post_lof)", syn["anchor_post_lof"]),
        summarize("non-depleting P/LP - all (post_lof)", nondep_plp["anchor_post_lof"]),
        summarize("non-depleting P/LP - missense only (post_lof)", nondep_missense["anchor_post_lof"]),
    ]
    summary_df = pd.DataFrame(summary_rows)

    detail_cols = [
        "Targeton_ID", "exon", "HGVSc", "HGVSp", "Summary_Plot", "condition", "clnsig_norm",
        "pos_adj_log2FoldChange_raw", "z_noimpact", "anchor_post_lof", "anchor_tier", "clinvar_id",
    ]
    detail_cols = [c for c in detail_cols if c in nondep_plp.columns]
    detail_df = nondep_plp[detail_cols].sort_values("z_noimpact")

    fig_path = os.path.join(args.outdir, "nondepleting_plp_vs_synonymous.png")
    ylim = make_plot(syn, nondep_plp, nondep_missense, fig_path)

    # Synonymous variants falling outside the plot's own off-scale range (same Tukey
    # fence used for the "N point(s) off-scale" annotation) -- worth a manual look,
    # since synonymous variants that score as depleted can be real splice-disrupting
    # artifacts rather than noise.
    syn_outlier_cols = [c for c in ["Targeton_ID", "position", "oligo_name", "pos_adj_log2FoldChange_raw",
                                     "z_noimpact", "anchor_post_lof", "anchor_tier"] if c in syn.columns]
    syn_outliers = syn[(syn["z_noimpact"] < ylim[0]) | (syn["z_noimpact"] > ylim[1])][syn_outlier_cols]
    syn_outliers = syn_outliers.sort_values("z_noimpact")

    # --- Write combined TSV (one file, multiple sections) ---
    combined_path = os.path.join(args.outdir, "nondepleting_plp_vs_synonymous_tables.tsv")
    with open(combined_path, "w") as f:
        f.write("# Table 1: group summary stats (median/IQR) for z-score and posterior P(LOF)\n")
        summary_df.to_csv(f, sep="\t", index=False)
        f.write("\n# Table 2: Mann-Whitney U tests, non-depleting P/LP vs synonymous reference\n")
        stats_df.to_csv(f, sep="\t", index=False)
        f.write("\n# Table 3: per-variant detail, non-depleting P/LP variants, sorted by z (most depletion-leaning first)\n")
        detail_df.to_csv(f, sep="\t", index=False)
        f.write(f"\n# Table 4: synonymous outliers outside the plot's Tukey-fence range (z < {ylim[0]:.2f} or z > {ylim[1]:.2f})\n")
        syn_outliers.to_csv(f, sep="\t", index=False)

    print(summary_df.to_string(index=False))
    print()
    print(stats_df.to_string(index=False))
    print()
    if len(syn_outliers):
        print(f"NOTE: {len(syn_outliers)} synonymous variant(s) fall outside the plot's shown range (see Table 4) -- "
              f"worth a manual look, e.g. via position, since synonymous variants near splice junctions can "
              f"genuinely disrupt splicing despite not changing the amino acid.")
        print()
    print("Wrote:")
    for p in [fig_path, combined_path]:
        print(f"  {p}")


if __name__ == "__main__":
    main()
