#!/usr/bin/env python3
"""
check_omim_coding_patterns.py

Follow-up to investigate_clinvar_variants.py: checks, GENE-WIDE across all
ClinVar-matched variants (not just one filtered subset), how often
submissions actually cite a DEE85-specific code/term (OMIM:301044, or
'developmental and epileptic encephalopathy' in the resolved
ReportedPhenotypeInfo text) versus a CdLS-specific code/term (OMIM:300590,
MedGen:C1802395, 'Cornelia de Lange', 'Congenital muscular
hypertrophy-cerebral syndrome') -- and whether that tracks consequence
(Summary_Plot) and anchor_tier better than ClinVar's own aggregated
'condition' column does.

Input: a *_submission_level.tsv produced by investigate_clinvar_variants.py.
For a true gene-wide check, first (re)run investigate_clinvar_variants.py
with NO filters (--condition/--consequence/--anchor_tier/--clnsig all
omitted) so it processes every ClinVar-matched variant in
clinvar_variants_summary.tsv, e.g.:

    python investigate_clinvar_variants.py \\
        --summary_tsv clinvar_variants_summary.tsv \\
        --submission_summary submission_summary.txt.gz \\
        --prefix all_clinvar_variants

...then feed the resulting all_clinvar_variants_submission_level.tsv into
this script:

    python check_omim_coding_patterns.py \\
        --input all_clinvar_variants_submission_level.tsv \\
        --output_prefix omim_coding_check
"""
import argparse
import re
import pandas as pd

DEE_PATTERNS = [
    r"301044",
    r"developmental[_ ]and[_ ]epileptic",
    r"dee85",
    r"epileptic[_ ]encephalopathy",
]
CDLS_PATTERNS = [
    r"300590",
    r"c1802395",
    r"cornelia[_ ]de[_ ]lange",
    r"congenital[_ ]muscular[_ ]hypertrophy",
    r"de[_ ]lange[_ ]syndrome",
]


def matches_any(text: str, patterns: list) -> bool:
    if not isinstance(text, str) or not text:
        return False
    t = text.lower()
    return any(re.search(p, t) for p in patterns)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True,
                         help="*_submission_level.tsv from investigate_clinvar_variants.py "
                              "(ideally from an unfiltered/all-variants run)")
    parser.add_argument("--output_prefix", default="omim_coding_check",
                         help="Prefix for output TSVs (default: omim_coding_check)")
    return parser.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.input, sep="\t")

    required = ["VariationID", "SubmittedPhenotypeInfo", "ReportedPhenotypeInfo"]
    for col in required:
        if col not in df.columns:
            raise SystemExit(f"ERROR: missing required column '{col}' in {args.input}")

    # Check both the raw submitted code/text AND ClinVar's resolved
    # CUI:name text, since a bare code like 'OMIM:300590' in
    # SubmittedPhenotypeInfo won't textually match 'cornelia' unless we
    # also check ReportedPhenotypeInfo (or the numeric code itself).
    combined_text = (df["SubmittedPhenotypeInfo"].fillna("") + " | " +
                      df["ReportedPhenotypeInfo"].fillna(""))
    df["mentions_dee"] = combined_text.apply(lambda t: matches_any(t, DEE_PATTERNS))
    df["mentions_cdls"] = combined_text.apply(lambda t: matches_any(t, CDLS_PATTERNS))

    submission_out = f"{args.output_prefix}_submission_level_flagged.tsv"
    df.to_csv(submission_out, sep="\t", index=False)
    print(f"Wrote per-submission flags: {submission_out}")

    # ---- Per-variant aggregation ----
    annot_cols = [c for c in ["Summary_Plot", "condition", "clnsig_norm", "anchor_tier",
                               "Targeton_ID", "HGVSc", "HGVSp"] if c in df.columns]

    agg_rows = []
    for vid, g in df.groupby("VariationID"):
        row = {"VariationID": vid, "n_submissions": len(g),
               "any_mentions_dee": bool(g["mentions_dee"].any()),
               "any_mentions_cdls": bool(g["mentions_cdls"].any()),
               "n_mentions_dee": int(g["mentions_dee"].sum()),
               "n_mentions_cdls": int(g["mentions_cdls"].sum())}
        for c in annot_cols:
            vals = g[c].dropna().unique()
            row[c] = vals[0] if len(vals) else ""
        agg_rows.append(row)
    variant_level = pd.DataFrame(agg_rows)

    def coding_pattern(row):
        if row["any_mentions_dee"] and row["any_mentions_cdls"]:
            return "both cited"
        if row["any_mentions_dee"]:
            return "DEE cited only"
        if row["any_mentions_cdls"]:
            return "CdLS cited only"
        return "neither (uncoded/other)"

    variant_level["coding_pattern"] = variant_level.apply(coding_pattern, axis=1)

    variant_out = f"{args.output_prefix}_variant_level.tsv"
    variant_level.to_csv(variant_out, sep="\t", index=False)
    print(f"Wrote per-variant coding pattern summary: {variant_out}")

    # ---- Summary crosstabs to the terminal ----
    print(f"\n{len(variant_level)} distinct variants checked.\n")

    print("Overall coding_pattern counts:")
    print(variant_level["coding_pattern"].value_counts().to_string())

    if "Summary_Plot" in variant_level.columns:
        print("\ncoding_pattern x consequence (Summary_Plot):")
        print(pd.crosstab(variant_level["Summary_Plot"], variant_level["coding_pattern"]).to_string())

    if "condition" in variant_level.columns:
        print("\ncoding_pattern x ClinVar's own 'condition' column:")
        print(pd.crosstab(variant_level["condition"], variant_level["coding_pattern"]).to_string())

    if "anchor_tier" in variant_level.columns:
        print("\ncoding_pattern x anchor_tier:")
        print(pd.crosstab(variant_level["anchor_tier"], variant_level["coding_pattern"]).to_string())

    if "Summary_Plot" in variant_level.columns:
        print("\n--- Key check: among LOF-consequence variants specifically, "
              "how often is DEE85 ever cited vs. only CdLS? ---")
        lof = variant_level[variant_level["Summary_Plot"] == "LOF"]
        print(lof["coding_pattern"].value_counts().to_string())


if __name__ == "__main__":
    main()
