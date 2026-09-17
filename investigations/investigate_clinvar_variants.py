#!/usr/bin/env python3
"""
investigate_clinvar_variants.py

Pull submission-level detail (reviewer status, date last evaluated,
submitter, collection method, SCV, whether each SCV contributes to the
aggregate classification) from ClinVar's submission_summary.txt.gz for a
chosen subset of variants in clinvar_variants_summary.tsv (the
cross-targeton summary produced by sge_clinvar_intersect_plot.py).

Built for exactly this kind of question: "this cluster of variants looks
functionally/clinically inconsistent -- are the underlying ClinVar
submissions old / single-submitter / low-review-status (suggesting stale
labeling), or well-supported by multiple recent, high-confidence
submissions?"

Select the variant subset either with filters (--condition, --consequence,
--anchor_tier, --clnsig) or by passing explicit --clinvar_ids. Filters
combine with AND; leave any filter unset to not apply it.

Outputs two TSVs:
  <prefix>_variant_level.tsv    One row per variant: the usual annotation
                                 columns (Targeton_ID, HGVSc, HGVSp,
                                 condition, clnsig_norm, Summary_Plot,
                                 anchor_tier, clndn_raw, clinvar_id) plus
                                 aggregated submission stats (n_submissions
                                 found, earliest/latest DateLastEvaluated,
                                 unique review statuses, unique submitters,
                                 unique collection methods) -- meant as the
                                 quick "does this look stale/thin?" table.
  <prefix>_submission_level.tsv One row per individual SCV, full detail,
                                 for manual read-through of specific
                                 submissions.

Usage examples:
    # The 42 CdLS-only LOF variants from the earlier analysis
    python investigate_clinvar_variants.py \\
        --summary_tsv clinvar_variants_summary.tsv \\
        --submission_summary submission_summary.txt.gz \\
        --condition "CdLS only" --consequence LOF \\
        --prefix cdls_lof_investigation

    # The 3 discordant "no impact" Pathogenic LOF variants
    python investigate_clinvar_variants.py \\
        --summary_tsv clinvar_variants_summary.tsv \\
        --submission_summary submission_summary.txt.gz \\
        --condition "CdLS only" --consequence LOF --anchor_tier "no impact" \\
        --clnsig "Pathogenic,Pathogenic/Likely pathogenic" \\
        --prefix cdls_lof_noimpact

    # An explicit list of VariationIDs
    python investigate_clinvar_variants.py \\
        --summary_tsv clinvar_variants_summary.tsv \\
        --submission_summary submission_summary.txt.gz \\
        --clinvar_ids 489224,1069969,1806292 \\
        --prefix manual_check
"""
import argparse
import gzip
import os
import sys
import pandas as pd


def read_submission_summary(path: str, variation_ids: set) -> pd.DataFrame:
    """Stream-parse submission_summary.txt.gz (one row per SCV), keeping
    only rows whose VariationID is in variation_ids. Format: '##' comment
    lines, then a '#VariationID...' header line, then data rows (data rows
    do NOT start with '#')."""
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
            if row["VariationID"] in variation_ids:
                rows.append(row)
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--summary_tsv", required=True,
                         help="clinvar_variants_summary.tsv from sge_clinvar_intersect_plot.py")
    parser.add_argument("--submission_summary", required=True,
                         help="Path to ClinVar's submission_summary.txt.gz")
    parser.add_argument("--condition", default=None,
                         help="Filter to this 'condition' value exactly, e.g. 'CdLS only'")
    parser.add_argument("--consequence", default=None,
                         help="Filter to this 'Summary_Plot' value exactly, e.g. 'LOF'")
    parser.add_argument("--anchor_tier", default=None,
                         help="Comma-separated anchor_tier values to include, "
                              "e.g. 'strongly depleting,weakly depleting'")
    parser.add_argument("--clnsig", default=None,
                         help="Comma-separated clnsig_norm values to include, "
                              "e.g. 'Pathogenic,Pathogenic/Likely pathogenic'")
    parser.add_argument("--clinvar_ids", default=None,
                         help="Comma-separated explicit clinvar_id/VariationID list. "
                              "If given, all other filters above are ignored.")
    parser.add_argument("--outdir", default=".",
                         help="Output directory (default: current directory)")
    parser.add_argument("--prefix", default="clinvar_investigation",
                         help="Prefix for output filenames (default: clinvar_investigation)")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.summary_tsv, sep="\t")
    for col in ["clinvar_id", "condition", "clnsig_norm", "Summary_Plot", "anchor_tier"]:
        if col not in df.columns:
            sys.exit(f"ERROR: '{col}' column not found in {args.summary_tsv}")

    if args.clinvar_ids:
        ids_wanted = {v.strip() for v in args.clinvar_ids.split(",") if v.strip()}
        subset = df[df["clinvar_id"].astype(str).isin(ids_wanted)].copy()
        print(f"Selected {len(subset)} row(s) by explicit --clinvar_ids "
              f"({len(ids_wanted)} ID(s) requested).")
    else:
        subset = df.copy()
        if args.condition:
            subset = subset[subset["condition"] == args.condition]
        if args.consequence:
            subset = subset[subset["Summary_Plot"] == args.consequence]
        if args.anchor_tier:
            wanted = {v.strip() for v in args.anchor_tier.split(",")}
            subset = subset[subset["anchor_tier"].isin(wanted)]
        if args.clnsig:
            wanted = {v.strip() for v in args.clnsig.split(",")}
            subset = subset[subset["clnsig_norm"].isin(wanted)]
        print(f"Selected {len(subset)} row(s) matching filters "
              f"(condition={args.condition!r}, consequence={args.consequence!r}, "
              f"anchor_tier={args.anchor_tier!r}, clnsig={args.clnsig!r}).")

    if subset.empty:
        sys.exit("No matching variants -- nothing to investigate. Check your filters.")

    # clinvar_id can be blank ("") for non-ClinVar rows -- drop those, and
    # de-duplicate in case the same variant appears in overlapping targetons.
    subset = subset[subset["clinvar_id"].astype(str).str.strip() != ""]
    subset = subset.drop_duplicates(subset=["clinvar_id"])
    variation_ids = set(subset["clinvar_id"].astype(str))
    print(f"{len(variation_ids)} distinct VariationID(s) to look up in {args.submission_summary}")

    sub_df = read_submission_summary(args.submission_summary, variation_ids)
    if sub_df.empty:
        print("WARNING: no matching rows found in submission_summary.txt.gz for these IDs. "
              "Still writing the variant-level table (submission stats will be blank).")

    # ---- Per-variant summary: the "does this look stale/thin?" table ----
    variant_cols = ["Targeton_ID", "HGVSc", "HGVSp", "Protein_position",
                     "Summary_Plot", "condition", "clnsig_norm", "anchor_tier",
                     "pos_adj_log2FoldChange_raw", "clindn_raw" if "clindn_raw" in subset.columns else "clndn_raw",
                     "clinvar_id"]
    variant_cols = [c for c in variant_cols if c in subset.columns]
    variant_level = subset[variant_cols].copy()
    variant_level["clinvar_id"] = variant_level["clinvar_id"].astype(str)

    if not sub_df.empty:
        agg_rows = []
        for vid, g in sub_df.groupby("VariationID"):
            dates = pd.to_datetime(g["DateLastEvaluated"], errors="coerce")
            agg_rows.append({
                "clinvar_id": vid,
                "n_submissions_found": len(g),
                "earliest_date_last_evaluated": dates.min().date().isoformat() if dates.notna().any() else "",
                "latest_date_last_evaluated": dates.max().date().isoformat() if dates.notna().any() else "",
                "review_statuses": ";".join(sorted(g["ReviewStatus"].dropna().unique())),
                "collection_methods": ";".join(sorted(g["CollectionMethod"].dropna().unique())),
                "submitters": ";".join(sorted(g["Submitter"].dropna().unique())),
                "n_contributing": int((g["ContributesToAggregateClassification"].str.lower() == "yes").sum())
                                   if "ContributesToAggregateClassification" in g.columns else None,
            })
        agg_df = pd.DataFrame(agg_rows)
        variant_level = variant_level.merge(agg_df, on="clinvar_id", how="left")
        variant_level["n_submissions_found"] = variant_level["n_submissions_found"].fillna(0).astype(int)
    else:
        for c in ["n_submissions_found", "earliest_date_last_evaluated",
                  "latest_date_last_evaluated", "review_statuses",
                  "collection_methods", "submitters", "n_contributing"]:
            variant_level[c] = "" if c not in ("n_submissions_found", "n_contributing") else 0

    variant_out = os.path.join(args.outdir, f"{args.prefix}_variant_level.tsv")
    variant_level.to_csv(variant_out, sep="\t", index=False)
    print(f"Wrote per-variant summary: {variant_out}")

    # ---- Full submission-level detail for manual read-through ----
    if not sub_df.empty:
        # Attach the variant's consequence/condition/anchor_tier alongside
        # each submission row, so the read-through doesn't require a
        # separate lookup.
        annot_cols = [c for c in ["clinvar_id", "Targeton_ID", "HGVSc", "HGVSp",
                                   "Summary_Plot", "condition", "clnsig_norm",
                                   "anchor_tier"] if c in subset.columns]
        annot = subset[annot_cols].copy()
        annot["clinvar_id"] = annot["clinvar_id"].astype(str)
        sub_annotated = sub_df.merge(
            annot, left_on="VariationID", right_on="clinvar_id", how="left"
        )
        sub_out = os.path.join(args.outdir, f"{args.prefix}_submission_level.tsv")
        sub_annotated.to_csv(sub_out, sep="\t", index=False)
        print(f"Wrote per-submission detail: {sub_out} ({len(sub_annotated)} row(s))")
    else:
        print("No submission-level file written (no matches found).")

    # ---- Quick eyeball summary to the terminal ----
    if not sub_df.empty:
        print("\n--- Quick summary ---")
        print(f"Variants with 0 matching submissions found: "
              f"{(variant_level['n_submissions_found'] == 0).sum()} / {len(variant_level)}")
        print(f"Variants with only 1 matching submission: "
              f"{(variant_level['n_submissions_found'] == 1).sum()} / {len(variant_level)}")
        if variant_level["latest_date_last_evaluated"].astype(str).str.strip().ne("").any():
            years = pd.to_datetime(variant_level["latest_date_last_evaluated"], errors="coerce").dt.year
            print("Latest DateLastEvaluated year distribution:")
            print(years.value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
