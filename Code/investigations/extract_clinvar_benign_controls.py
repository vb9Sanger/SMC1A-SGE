#!/usr/bin/env python3
"""
extract_clinvar_benign_controls.py

Prepares an additional benign-control set from ClinVar, to widen the curated
truth set's small Exclude_Benign group (n=25) for the sensitivity/specificity/
OddsPath calibration (09_calibrate_sensitivity_oddspath.py). Disease-condition
attribution -- the specific thing ClinVar was found unreliable for in this
gene (see repo README STEP TEN) -- does not apply to a Benign/Likely benign
call, so ClinVar's own classification can be used here without inheriting
that specific caveat. Review-status tier is still worth checking, since most
ClinVar submissions for this gene are single-submitter (see --min_review_tier).

Two independent inputs are needed because clinvar_variants_summary.tsv (from
sge_clinvar_intersect.py) does not itself carry ClinVar's review-status field
(CLNREVSTAT):
  1. --summary_tsv   clinvar_variants_summary.tsv -- supplies anchor_tier/
                      anchor_call and clnsig_norm, already joined to the assay.
  2. --clinvar_vcf    the same whole-gene(-region) ClinVar VCF used upstream --
                      supplies CLNREVSTAT via a direct bcftools query, since
                      the intersect script never extracted it.

Output is a plain (variant_key, anchor_call) TSV in the same shape
09_calibrate_sensitivity_oddspath.py's --extra_group_tsv expects, i.e. it can
be fed straight back into that script as a new group.

Requirements: pandas, bcftools on PATH.

Usage:
    python extract_clinvar_benign_controls.py \\
        --summary_tsv clinvar_variants_summary.tsv \\
        --clinvar_vcf clinvar.vcf.gz \\
        --region X:53374149-53422654 \\
        --min_review_tier any \\
        --output clinvar_benign_controls.tsv

--min_review_tier:
    any             -- every Benign/Likely benign/Benign_Likely_benign row
                       (default; matches the "All ClinVar B/LB" figure)
    multi_submitter -- criteria_provided,_multiple_submitters,_no_conflicts
                       or better only (matches the "high-confidence tier"
                       figure; expect far fewer rows)
"""

import argparse
import subprocess
import sys

import pandas as pd

BENIGN_VALUES = {"Benign", "Likely benign", "Benign/Likely benign"}

# Normalises clinvar_variants_summary.tsv's Summary_Plot vocabulary to
# assay_join_all.tsv's consequence_class vocabulary, so a --consequence_filter
# on 09_calibrate_sensitivity_oddspath.py can be applied consistently to both
# the curated-set groups and this ClinVar-derived one. Only the classes that
# actually matter for that filter are mapped explicitly; anything else passes
# through lowercased as a reasonable fallback (not expected to be filtered on).
SUMMARY_PLOT_TO_CONSEQUENCE_CLASS = {
    "Missense_Variant": "missense",
    "Inframe_Deletion": "inframe_deletion",
    "Inframe_Insertion": "inframe_insertion",
    "Synonymous_Variant": "synonymous",
    "Intronic_Variant": "intronic",
}

REVIEW_TIER_RANK = {
    "no_assertion_criteria_provided": 0,
    "criteria_provided,_conflicting_classifications": 0,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "reviewed_by_expert_panel": 3,
    "practice_guideline": 3,
}


def fetch_revstat(clinvar_vcf: str, region: str) -> dict:
    """Returns {(pos, ref, alt): CLNREVSTAT} for every record in region."""
    cmd = [
        "bcftools", "view", "-r", region, clinvar_vcf,
    ]
    view = subprocess.run(cmd, capture_output=True, text=True, check=True)
    query = subprocess.run(
        ["bcftools", "query", "-f", "%POS\t%REF\t%ALT\t%INFO/CLNREVSTAT\n"],
        input=view.stdout, capture_output=True, text=True, check=True,
    )
    lookup = {}
    for line in query.stdout.strip().splitlines():
        pos, ref, alt, revstat = line.split("\t")
        lookup[(pos, ref, alt)] = revstat
    return lookup


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary_tsv", required=True)
    ap.add_argument("--clinvar_vcf", required=True)
    ap.add_argument("--region", required=True,
                     help="e.g. X:53374149-53422654 (bcftools region syntax)")
    ap.add_argument("--min_review_tier", choices=["any", "multi_submitter"],
                     default="any")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.summary_tsv, sep="\t", dtype=str)
    df = df[df["clnsig_norm"].isin(BENIGN_VALUES)].copy()
    print(f"{len(df)} Benign/Likely benign row(s) in {args.summary_tsv}", file=sys.stderr)

    revstat_lookup = fetch_revstat(args.clinvar_vcf, args.region)
    print(f"{len(revstat_lookup)} ClinVar record(s) fetched for region {args.region}",
          file=sys.stderr)

    df["position_int"] = df["position"].astype(float).astype(int).astype(str)
    df["review_status"] = df.apply(
        lambda r: revstat_lookup.get((r["position_int"], r["vcf_ref"], r["vcf_alt"]), "."),
        axis=1,
    )
    n_unmatched = (df["review_status"] == ".").sum()
    if n_unmatched:
        print(f"NOTE: {n_unmatched} row(s) had no matching CLNREVSTAT record "
              f"(kept, tier treated as unknown/lowest)", file=sys.stderr)

    if args.min_review_tier == "multi_submitter":
        df["tier_rank"] = df["review_status"].map(REVIEW_TIER_RANK).fillna(0)
        before = len(df)
        df = df[df["tier_rank"] >= 2]
        print(f"Filtered to multi-submitter-or-better: {before} -> {len(df)} row(s)",
              file=sys.stderr)

    df["variant_key"] = df["position_int"] + ":" + df["vcf_ref"] + ":" + df["vcf_alt"]
    df["anchor_call"] = df["anchor_tier"].map(
        lambda t: "depleted" if t in ("strongly depleting", "weakly depleting") else t
    )
    df["consequence_class"] = df["Summary_Plot"].map(
        lambda s: SUMMARY_PLOT_TO_CONSEQUENCE_CLASS.get(s, s.lower())
    )

    out = df[["variant_key", "anchor_call", "anchor_tier", "clnsig_norm",
              "review_status", "consequence_class"]].drop_duplicates("variant_key")
    out.to_csv(args.output, sep="\t", index=False)
    print(f"\nWrote {args.output} ({len(out)} unique variant(s))", file=sys.stderr)
    print(out["anchor_tier"].value_counts().to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
