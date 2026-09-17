#!/usr/bin/env python3
"""
cassette_exon_consequence_analysis.py

Gene-wide, variant-level test of whether cassette-exon-annotated exons
behave differently from non-cassette exons, run in parallel for two
variant classes:

  1. Synonymous variants -- continuous z-score comparison (as before).
  2. Splice-region variants ("Splice_Variant" + "Splice_Polypyrimidine_
     Tract_Variant" by default) -- BOTH a continuous z-score comparison
     AND a categorical comparison of % depleted (anchor_tier in
     {weakly depleting, strongly depleting}) -- this directly tests the
     original observation ("more splice region variants depleting in
     certain exons") rather than using synonymous variants as a proxy.

z_noimpact = (LFC - anchor_mu_noimpact) / anchor_sd_noimpact
  i.e. how many no-impact-control SDs a variant's LFC sits from its own
  targeton's no-impact anchor. More negative = more depletion-like.

Inputs:
  --anchor_dir        Directory of per-targeton anchor-tier output TSVs (the
                       full gaussian_shrinkage_classifier.R output, one file
                       per targeton, needs consequence, anchor_tier,
                       pos_adj_log2FoldChange_raw, anchor_mu_noimpact,
                       anchor_sd_noimpact). Targeton_ID is inferred from
                       each filename.
  --exon_events_tsv   The exon-level annotation table (e.g.
                       SMC1A_exon_splice_events_full.tsv) with columns
                       Targeton_ID and events_within_exon -- a targeton
                       counts as "cassette" if "cassetteExon" appears in its
                       events_within_exon string.
  --splice_labels     Comma-separated consequence labels counted as
                       "splice-region" (default:
                       Splice_Variant,Splice_Polypyrimidine_Tract_Variant)

Outputs (single combined TSV + one figure):
  - cassette_exon_consequence_tables.tsv:
      Table 1: per-exon summary, synonymous variants
      Table 2: pooled group summary + Mann-Whitney, synonymous
      Table 3: per-exon summary, splice-region variants (n, %depleted, median z)
      Table 4: pooled group summary + Mann-Whitney, splice-region z-scores
      Table 5: 2x2 depletion counts + Fisher's exact test, splice-region
  - cassette_exon_consequence_analysis.png: 2x2 grid
      Row 1: synonymous (per-exon bar, pooled violin)
      Row 2: splice-region (per-exon %-depleted bar, pooled %-depleted bar)

Usage:
    python cassette_exon_consequence_analysis.py \
        --anchor_dir per_targeton_results/ \
        --exon_events_tsv SMC1A_exon_splice_events_full.tsv \
        --outdir results/
"""

import argparse
import glob
import os
import re
import sys

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

DEPLETED_TIERS = {"weakly depleting", "strongly depleting"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--anchor_dir", required=True, help="Directory of per-targeton anchor-tier output TSVs")
    p.add_argument("--exon_events_tsv", required=True, help="Exon-level annotation table")
    p.add_argument("--outdir", default=".", help="Output directory")
    p.add_argument("--consequence_col", default="consequence", help="Consequence column name (default: consequence)")
    p.add_argument("--tier_col", default="anchor_tier", help="Anchor tier column name (default: anchor_tier)")
    p.add_argument(
        "--splice_labels", default="Splice_Variant,Splice_Polypyrimidine_Tract_Variant",
        help="Comma-separated consequence labels counted as splice-region",
    )
    return p.parse_args()


def infer_targeton_name(path):
    bn = os.path.basename(path).replace(".tsv", "")
    return bn.split("_all_deseq2_results_condition_")[0]


def normalize_targeton_id(name):
    return re.sub(r"_exon\d+$", "", name)


def load_variants(anchor_dir, consequence_col, tier_col):
    files = sorted(glob.glob(os.path.join(anchor_dir, "*.tsv")))
    files = [f for f in files if not os.path.basename(f).startswith("gmm_shrinkage_anchor_summary")]
    if not files:
        sys.exit(f"ERROR: no *.tsv files found in --anchor_dir {anchor_dir}")

    required = [consequence_col, tier_col, "pos_adj_log2FoldChange_raw", "anchor_mu_noimpact", "anchor_sd_noimpact"]
    frames = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        missing = [c for c in required if c not in df.columns]
        if missing:
            sys.exit(f"ERROR: {f} is missing required column(s): {missing}")
        df["Targeton_ID"] = normalize_targeton_id(infer_targeton_name(f))
        frames.append(df)

    out = pd.concat(frames, ignore_index=True)
    print(f"Loaded {len(out)} variants (all consequence classes) across {len(files)} targeton files.")
    return out


def add_cassette_status(df, events_path):
    events = pd.read_csv(events_path, sep="\t")
    if "Targeton_ID" not in events.columns or "events_within_exon" not in events.columns:
        sys.exit("ERROR: --exon_events_tsv must have 'Targeton_ID' and 'events_within_exon' columns")
    events["is_cassette"] = events["events_within_exon"].fillna("").str.contains("cassetteExon")
    cassette_map = dict(zip(events["Targeton_ID"], events["is_cassette"]))
    exon_map = dict(zip(events["Targeton_ID"], events["EXON"])) if "EXON" in events.columns else {}

    df["is_cassette"] = df["Targeton_ID"].map(cassette_map)
    unmapped = df[df["is_cassette"].isna()]["Targeton_ID"].unique()
    if len(unmapped):
        print(f"NOTE: {len(unmapped)} targeton(s) have no entry in --exon_events_tsv, excluded: {sorted(unmapped)}")
    df = df.dropna(subset=["is_cassette"]).copy()
    df["EXON"] = df["Targeton_ID"].map(exon_map)
    return df, cassette_map, exon_map


def per_exon_summary_continuous(df, value_col, label_col="Targeton_ID"):
    rows = []
    for targ, sub in df.groupby(label_col):
        rows.append({
            "Targeton_ID": targ, "EXON": sub["EXON"].iloc[0], "is_cassette": sub["is_cassette"].iloc[0],
            "n": len(sub), "median_z": round(sub[value_col].median(), 3),
            "q25_z": round(sub[value_col].quantile(0.25), 3), "q75_z": round(sub[value_col].quantile(0.75), 3),
        })
    out = pd.DataFrame(rows)
    if out["EXON"].notna().any():
        out = out.sort_values("EXON")
    return out


def per_exon_summary_depletion(df, tier_col, label_col="Targeton_ID"):
    rows = []
    for targ, sub in df.groupby(label_col):
        n = len(sub)
        n_dep = sub[tier_col].isin(DEPLETED_TIERS).sum()
        rows.append({
            "Targeton_ID": targ, "EXON": sub["EXON"].iloc[0], "is_cassette": sub["is_cassette"].iloc[0],
            "n": n, "n_depleted": int(n_dep), "pct_depleted": round(100 * n_dep / n, 1) if n else 0.0,
        })
    out = pd.DataFrame(rows)
    if out["EXON"].notna().any():
        out = out.sort_values("EXON")
    return out


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    splice_labels = [s.strip() for s in args.splice_labels.split(",") if s.strip()]

    all_variants = load_variants(args.anchor_dir, args.consequence_col, args.tier_col)
    all_variants["z_noimpact"] = (
        (all_variants["pos_adj_log2FoldChange_raw"] - all_variants["anchor_mu_noimpact"])
        / all_variants["anchor_sd_noimpact"]
    )
    all_variants, cassette_map, exon_map = add_cassette_status(all_variants, args.exon_events_tsv)

    syn = all_variants[all_variants[args.consequence_col] == "Synonymous_Variant"].copy()
    splice = all_variants[all_variants[args.consequence_col].isin(splice_labels)].copy()

    print(f"Synonymous variants: n={len(syn)} | Splice-region variants ({'+'.join(splice_labels)}): n={len(splice)}")

    for name, d in [("synonymous", syn), ("splice-region", splice)]:
        if d["is_cassette"].sum() < 2 or (~d["is_cassette"]).sum() < 2:
            sys.exit(f"ERROR: not enough {name} variants in both cassette and non-cassette exons to compare.")

    # --- Synonymous: per-exon + pooled z-score comparison ---
    syn_per_exon = per_exon_summary_continuous(syn, "z_noimpact")
    syn_cassette_z = syn[syn["is_cassette"]]["z_noimpact"]
    syn_noncassette_z = syn[~syn["is_cassette"]]["z_noimpact"]
    u_syn, p_syn = mannwhitneyu(syn_cassette_z, syn_noncassette_z, alternative="two-sided")
    syn_group_df = pd.DataFrame([
        {"group": "cassette exons", "n": len(syn_cassette_z), "median_z": round(syn_cassette_z.median(), 3)},
        {"group": "non-cassette exons", "n": len(syn_noncassette_z), "median_z": round(syn_noncassette_z.median(), 3)},
    ])
    syn_stats_df = pd.DataFrame([{"comparison": "synonymous z_noimpact: cassette vs non-cassette",
                                   "n1": len(syn_cassette_z), "n2": len(syn_noncassette_z), "U": u_syn, "p_value": p_syn}])

    # --- Splice-region: per-exon + pooled z-score AND %-depleted comparison ---
    splice_per_exon_z = per_exon_summary_continuous(splice, "z_noimpact")
    splice_per_exon_dep = per_exon_summary_depletion(splice, args.tier_col)

    splice_cassette_z = splice[splice["is_cassette"]]["z_noimpact"]
    splice_noncassette_z = splice[~splice["is_cassette"]]["z_noimpact"]
    u_spl, p_spl = mannwhitneyu(splice_cassette_z, splice_noncassette_z, alternative="two-sided")
    splice_group_df = pd.DataFrame([
        {"group": "cassette exons", "n": len(splice_cassette_z), "median_z": round(splice_cassette_z.median(), 3)},
        {"group": "non-cassette exons", "n": len(splice_noncassette_z), "median_z": round(splice_noncassette_z.median(), 3)},
    ])
    splice_zstats_df = pd.DataFrame([{"comparison": "splice-region z_noimpact: cassette vs non-cassette",
                                       "n1": len(splice_cassette_z), "n2": len(splice_noncassette_z), "U": u_spl, "p_value": p_spl}])

    n_cas = len(splice[splice["is_cassette"]])
    n_cas_dep = splice[splice["is_cassette"]][args.tier_col].isin(DEPLETED_TIERS).sum()
    n_noncas = len(splice[~splice["is_cassette"]])
    n_noncas_dep = splice[~splice["is_cassette"]][args.tier_col].isin(DEPLETED_TIERS).sum()

    table_2x2 = pd.DataFrame(
        [[n_cas_dep, n_cas - n_cas_dep], [n_noncas_dep, n_noncas - n_noncas_dep]],
        index=["cassette exons", "non-cassette exons"], columns=["depleted", "not_depleted"],
    )
    odds, p_fisher = fisher_exact(table_2x2.values)
    fisher_df = table_2x2.reset_index().rename(columns={"index": "group"})
    fisher_df["pct_depleted"] = round(100 * fisher_df["depleted"] / (fisher_df["depleted"] + fisher_df["not_depleted"]), 1)
    fisher_stats_row = pd.DataFrame([{"comparison": "splice-region %% depleted: cassette vs non-cassette",
                                       "odds_ratio": round(odds, 3), "p_value": p_fisher}])

    # --- Combined TSV ---
    combined_path = os.path.join(args.outdir, "cassette_exon_consequence_tables.tsv")
    with open(combined_path, "w") as f:
        f.write("# Table 1: per-exon summary, synonymous variants\n")
        syn_per_exon.to_csv(f, sep="\t", index=False)
        f.write("\n# Table 2: pooled group summary + Mann-Whitney, synonymous z-scores\n")
        syn_group_df.to_csv(f, sep="\t", index=False)
        syn_stats_df.to_csv(f, sep="\t", index=False)
        f.write(f"\n# Table 3: per-exon summary, splice-region variants ({'+'.join(splice_labels)})\n")
        splice_per_exon_dep.merge(
            splice_per_exon_z[["Targeton_ID", "median_z"]], on="Targeton_ID", how="left"
        ).to_csv(f, sep="\t", index=False)
        f.write("\n# Table 4: pooled group summary + Mann-Whitney, splice-region z-scores\n")
        splice_group_df.to_csv(f, sep="\t", index=False)
        splice_zstats_df.to_csv(f, sep="\t", index=False)
        f.write("\n# Table 5: 2x2 depletion counts + Fisher's exact test, splice-region variants\n")
        fisher_df.to_csv(f, sep="\t", index=False)
        fisher_stats_row.to_csv(f, sep="\t", index=False)

    # --- Figure: 2x2 grid ---
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # Row 1: synonymous (unchanged from before)
    ax = axes[0, 0]
    colors = syn_per_exon["is_cassette"].map({True: "#EA4335", False: "#4285F4"})
    x = syn_per_exon["EXON"] if syn_per_exon["EXON"].notna().any() else range(len(syn_per_exon))
    ax.bar(x, syn_per_exon["median_z"], color=colors, width=0.7)
    ax.axhline(0, color="grey", linestyle="--", linewidth=1)
    ax.set_xlabel("exon number")
    ax.set_ylabel("median z (synonymous)")
    ax.set_title("Synonymous: per-exon median z")
    ax.legend(handles=[Patch(color="#EA4335", label="cassette"), Patch(color="#4285F4", label="non-cassette")],
              loc="best", fontsize=8, frameon=False)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)

    ax = axes[0, 1]
    data = [syn_cassette_z.values, syn_noncassette_z.values]
    parts = ax.violinplot(data, showmedians=True, widths=0.7)
    for pc in parts["bodies"]:
        pc.set_facecolor("#4285F4"); pc.set_alpha(0.35)
    for i, d in enumerate(data, start=1):
        jitter = np.random.default_rng(0).normal(0, 0.05, size=len(d))
        ax.scatter(np.full(len(d), i) + jitter, d, s=6, color="#EA4335", alpha=0.4, zorder=3)
    ax.axhline(0, color="grey", linestyle="--", linewidth=1)
    ax.set_xticks([1, 2]); ax.set_xticklabels([f"Cassette\n(n={len(syn_cassette_z)})", f"Non-cassette\n(n={len(syn_noncassette_z)})"])
    ax.set_ylabel("z (synonymous)")
    ax.set_title(f"Synonymous: pooled (Mann-Whitney p={p_syn:.3g})")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)

    # Row 2: splice-region, % depleted
    ax = axes[1, 0]
    colors = splice_per_exon_dep["is_cassette"].map({True: "#EA4335", False: "#4285F4"})
    x = splice_per_exon_dep["EXON"] if splice_per_exon_dep["EXON"].notna().any() else range(len(splice_per_exon_dep))
    ax.bar(x, splice_per_exon_dep["pct_depleted"], color=colors, width=0.7)
    ax.set_xlabel("exon number")
    ax.set_ylabel("% splice-region variants depleted")
    ax.set_title("Splice-region: per-exon % depleted")
    ax.legend(handles=[Patch(color="#EA4335", label="cassette"), Patch(color="#4285F4", label="non-cassette")],
              loc="best", fontsize=8, frameon=False)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)

    ax = axes[1, 1]
    bar_x = [0, 1]
    pct = fisher_df.set_index("group").loc[["cassette exons", "non-cassette exons"], "pct_depleted"]
    ax.bar(bar_x, pct.values, color=["#EA4335", "#4285F4"], width=0.5)
    ax.set_xticks(bar_x)
    ax.set_xticklabels([f"Cassette\n(n={n_cas})", f"Non-cassette\n(n={n_noncas})"])
    ax.set_ylabel("% splice-region variants depleted")
    ax.set_title(f"Splice-region: pooled (Fisher's p={p_fisher:.3g}, OR={odds:.2f})")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)

    plt.tight_layout()
    fig_path = os.path.join(args.outdir, "cassette_exon_consequence_analysis.png")
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print()
    print("=== Synonymous ===")
    print(syn_group_df.to_string(index=False))
    print(syn_stats_df.to_string(index=False))
    print()
    print("=== Splice-region ===")
    print(splice_group_df.to_string(index=False))
    print(splice_zstats_df.to_string(index=False))
    print()
    print(fisher_df.to_string(index=False))
    print(fisher_stats_row.to_string(index=False))
    print()
    print("Wrote:")
    for p_ in [fig_path, combined_path]:
        print(f"  {p_}")


if __name__ == "__main__":
    main()
