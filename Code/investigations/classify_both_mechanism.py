#!/usr/bin/env python3
"""
classify_both_mechanism.py

For every variant labeled condition == "Both" (i.e. ClinVar's aggregated
CLNDN field contains both a CdLS-associated and a DEE85-associated term),
determines WHY it ended up "Both":

  - "single submission lists both conditions": at least one individual
    submission (SCV), on its own, cites both a CdLS term and a DEE term
    together (e.g. SubmittedPhenotypeInfo = "OMIM:300590;OMIM:301044") --
    i.e. one lab attached both conditions to one classification, rather
    than the "Both" label reflecting different labs disagreeing.

  - "different submitters disagree (one CdLS, one DEE)": no single
    submission mentions both terms, but across the variant's submissions,
    at least one cites CdLS only and a separate one cites DEE only.

  - "other/unclear": neither of the above (shouldn't normally occur for a
    variant already labeled "Both", but included for completeness/safety).

Input: a *_submission_level_flagged.tsv produced by
check_omim_coding_patterns.py (which itself expects a
*_submission_level.tsv from investigate_clinvar_variants.py -- ideally an
UNFILTERED run covering all ClinVar-matched variants, or at minimum a run
that includes the "Both" condition variants you want to check).

Usage:
    python classify_both_mechanism.py \\
        --input omim_coding_check_submission_level_flagged.tsv \\
        --output both_variants_mechanism.tsv
"""
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True,
                         help="*_submission_level_flagged.tsv from check_omim_coding_patterns.py "
                              "(must have VariationID, condition, mentions_dee, mentions_cdls columns)")
    parser.add_argument("--output", default="both_variants_mechanism.tsv",
                         help="Output TSV path (default: both_variants_mechanism.tsv)")
    parser.add_argument("--condition_value", default="Both",
                         help="Which value of the 'condition' column to investigate "
                              "(default: 'Both')")
    return parser.parse_args()


def main():
    args = parse_args()
    sl = pd.read_csv(args.input, sep="\t")

    required = ["VariationID", "condition", "mentions_dee", "mentions_cdls"]
    for col in required:
        if col not in sl.columns:
            raise SystemExit(f"ERROR: missing required column '{col}' in {args.input}. "
                              f"This script expects the *_submission_level_flagged.tsv "
                              f"output of check_omim_coding_patterns.py, not the raw "
                              f"submission_level.tsv.")

    target_ids = sl[sl["condition"] == args.condition_value]["VariationID"].unique()
    if len(target_ids) == 0:
        raise SystemExit(f"No rows found with condition == '{args.condition_value}' in {args.input}")

    rows = []
    for vid, g in sl[sl["VariationID"].isin(target_ids)].groupby("VariationID"):
        # A submission "individually" cites both if that single row's own
        # mentions_dee AND mentions_cdls flags are both True (i.e. that
        # one SCV's SubmittedPhenotypeInfo/ReportedPhenotypeInfo text
        # matched both term sets on its own).
        n_dual = int((g["mentions_dee"] & g["mentions_cdls"]).sum())
        n_cdls_only = int((g["mentions_cdls"] & ~g["mentions_dee"]).sum())
        n_dee_only = int((g["mentions_dee"] & ~g["mentions_cdls"]).sum())
        n_neither = int((~g["mentions_dee"] & ~g["mentions_cdls"]).sum())

        if n_dual > 0:
            mechanism = "single submission lists both conditions"
        elif n_cdls_only > 0 and n_dee_only > 0:
            mechanism = "different submitters disagree (one CdLS, one DEE)"
        else:
            mechanism = "other/unclear"

        rows.append({
            "VariationID": vid,
            "n_submissions_total": len(g),
            "n_submissions_citing_both_in_one_row": n_dual,
            "n_submissions_citing_cdls_only": n_cdls_only,
            "n_submissions_citing_dee_only": n_dee_only,
            "n_submissions_citing_neither": n_neither,
            "mechanism": mechanism,
        })

    out = pd.DataFrame(rows).sort_values(["mechanism", "n_submissions_total"])
    out.to_csv(args.output, sep="\t", index=False)

    print(f"{len(target_ids)} distinct variant(s) with condition == '{args.condition_value}' checked.")
    print(f"Wrote per-variant mechanism table: {args.output}\n")
    print("Mechanism breakdown:")
    print(out["mechanism"].value_counts().to_string())
    print("\nn_submissions_total distribution, split by mechanism:")
    print(pd.crosstab(out["n_submissions_total"], out["mechanism"]).to_string())


if __name__ == "__main__":
    main()
