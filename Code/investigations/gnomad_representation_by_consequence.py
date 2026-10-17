#!/usr/bin/env python3
"""
gnomad_representation_by_consequence.py

Tests whether variants the assay classifies as depleting are under-represented
in gnomAD relative to those it classifies as no impact -- first pooled across
the gene, then **within each consequence class**.

The splice class is split into canonical acceptor and donor sites and splice
regions before stratifying. Summary_Consequence pools them, and the two carry
baselines an order of magnitude apart, so the pooled row averages over a
confound of exactly the kind this script exists to remove.

Why the stratified test is the one to quote
-------------------------------------------
The pooled comparison is confounded, and badly. Baseline gnomAD representation
varies by more than twentyfold between consequence classes: synonymous and
intronic variants are well represented, missense variants barely at all. The
assay's depleting calls are not spread evenly across those classes -- they are
concentrated in exactly the classes where baseline representation is lowest.
So a pooled depleting-vs-no-impact odds ratio partly measures the consequence
composition of the depleting set rather than the effect of depletion itself,
and it overstates the contrast.

Stratifying removes that confound: within a single consequence class, both
groups share the same baseline, and the comparison isolates the assay call.
The pooled figure is still computed and printed, so the size of the confound
is visible rather than hidden, but the per-class table is the result.

Two classes of caveat, both reported rather than silently handled
-----------------------------------------------------------------
**Classes with no comparison pool.** Nonsense, frameshift and start-lost
variants are classified as depleting almost without exception, so there is no
no-impact group within the class to compare against. These rows are computed
and marked `no_comparison_pool` rather than dropped, since their absence from
a published table is itself worth being able to justify.

**Design-incorporated variants inflate the intronic class.** Variants already
reported in gnomAD v2.1/v3.1 were deliberately incorporated into the library
at the design stage, so a fraction of the gnomAD matches arrived by that route
rather than by the assay independently sampling population variation. Those
variants are overwhelmingly intronic and no impact, which inflates the
intronic baseline specifically. `--systematic_only` repeats the whole analysis
over variants reached by the systematic mutators alone, which is the clean
comparison; the intronic row is the only one that moves materially.

Counts cover **every** assayed variant, matched or not -- testing for absence
requires the full denominator.

Usage:
    python gnomad_representation_by_consequence.py \\
        --gnomad_summary sge_gnomad_summary.tsv \\
        --outdir gnomad_intersect_results/
"""
import argparse
import sys

import pandas as pd
from scipy.stats import fisher_exact

# Classified as depleting at either strength. `enriched` is deliberately not
# pooled in: an enriched call is a different observation, not a weaker
# depleting one, and it is tested separately below.
DEPLETING = {"weakly depleting", "strongly depleting"}

# Classes where essentially every variant is classified as depleting, leaving
# no within-class no-impact pool. Reported, flagged, not dropped.
NO_POOL_THRESHOLD = 5   # minimum variants needed in *both* groups to test

# Acceptor and donor changes only. Summary_Consequence ranks every splice
# annotation above the protein-level one, so its Splice_Variant class pools
# these with splice-region changes that VEP scores below missense. Pooled, the
# class carries two baselines at once and the row means nothing; split, the
# splice-region row is the one with a usable comparison pool. Exact term
# matching, not substring: splice_donor_variant has to be told from
# splice_donor_region_variant and splice_donor_5th_base_variant.
CANONICAL_SPLICE_TERMS = {"splice_acceptor_variant", "splice_donor_variant"}


def consequence_class(df):
    """Summary_Consequence, with the splice class split on the raw VEP terms."""
    for col in ("Summary_Consequence", "Consequence"):
        if col not in df.columns:
            sys.exit(f"ERROR: --gnomad_summary has no '{col}' column")
    cls = df["Summary_Consequence"].copy()
    is_splice = cls == "Splice_Variant"
    canonical = df["Consequence"].fillna("").str.split(",").apply(
        lambda t: bool(CANONICAL_SPLICE_TERMS.intersection(t)))
    cls[is_splice & canonical] = "Canonical_Splice"
    cls[is_splice & ~canonical] = "Splice_Region"
    return cls

# Oligo-name suffixes marking a variant incorporated at the design stage
# because it was already reported, rather than reached by a systematic
# mutator. `--systematic_only` drops these and exits if it matches nothing,
# since a silently empty filter would look like a null sensitivity result.
DESIGN_TAGS = {"v2.1", "v3.1", "clinvar", "supplementary"}


def contingency(df):
    """(k_depleting, n_depleting, k_noimpact, n_noimpact) for one subset."""
    dep = df[df.anchor_tier.isin(DEPLETING)]
    nim = df[df.anchor_tier == "no impact"]
    return (int(dep.in_gnomad.sum()), len(dep),
            int(nim.in_gnomad.sum()), len(nim))


def test_row(label, df):
    k_d, n_d, k_n, n_n = contingency(df)
    row = {
        "consequence_class": label,
        "n_no_impact": n_n,
        "pct_gnomad_no_impact": round(100 * k_n / n_n, 2) if n_n else None,
        "n_depleting": n_d,
        "pct_gnomad_depleting": round(100 * k_d / n_d, 2) if n_d else None,
    }
    if n_d < NO_POOL_THRESHOLD or n_n < NO_POOL_THRESHOLD:
        row.update(odds_ratio=None, p_value=None, status="no_comparison_pool")
        return row
    # Two-sided. The direction is strongly expected, but these comparisons are
    # a check on whether the assay's calls track population constraint at all,
    # not a test of a pre-registered directional hypothesis, and a two-sided
    # test is the conservative choice for that. The calibration scripts, where
    # the direction genuinely is specified in advance, test one-sided instead.
    odds, p = fisher_exact([[k_d, n_d - k_d], [k_n, n_n - k_n]])
    row.update(odds_ratio=round(odds, 4), p_value=p, status="tested")
    return row


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gnomad_summary", required=True,
                    help="sge_gnomad_summary.tsv from sge_gnomad_intersect_plot.py")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--systematic_only", action="store_true",
                    help="restrict to variants reached by the systematic "
                         "mutators, excluding those incorporated at design "
                         "because they were already reported in gnomAD or "
                         "ClinVar (oligo-name suffixes in DESIGN_TAGS)")
    a = ap.parse_args()

    d = pd.read_csv(a.gnomad_summary, sep="\t", low_memory=False)
    d["in_gnomad"] = d["in_gnomad"].astype(bool)
    d = d[d.anchor_tier.notna()]

    if a.systematic_only:
        before = len(d)
        tag = d.oligo_name.str.rsplit("_", n=1).str[-1]
        d = d[~tag.isin(DESIGN_TAGS)]
        if before == len(d):
            sys.exit("--systematic_only matched nothing: oligo-name suffixes "
                     "are " + ", ".join(sorted(tag.unique())) + ". Update "
                     "DESIGN_TAGS rather than assuming the filter ran.")
        print(f"systematic-only: {before - len(d)} design-incorporated "
              f"variants excluded, {len(d)} remain\n")

    # --- pooled, reported for contrast rather than as the result -----------
    k_d, n_d, k_n, n_n = contingency(d)
    odds, p = fisher_exact([[k_d, n_d - k_d], [k_n, n_n - k_n]])
    print(f"POOLED across all consequence classes")
    print(f"  no impact  {k_n:>6}/{n_n:<6} {100*k_n/n_n:6.2f}%")
    print(f"  depleting  {k_d:>6}/{n_d:<6} {100*k_d/n_d:6.2f}%")
    print(f"  OR {odds:.4f}, p = {p:.2e}")
    print("  (confounded by consequence composition -- see module docstring)\n")

    enr = d[d.anchor_tier == "enriched"]
    if len(enr):
        k_e = int(enr.in_gnomad.sum())
        odds_e, p_e = fisher_exact(
            [[k_e, len(enr) - k_e], [k_n, n_n - k_n]])
        print(f"ENRICHED vs no impact: {k_e}/{len(enr)} "
              f"({100*k_e/len(enr):.2f}%), OR {odds_e:.2f}, p = {p_e:.3g}\n")

    # --- stratified: the reportable comparison ------------------------------
    d = d.assign(consequence_class=consequence_class(d))
    rows = [test_row(c, g) for c, g in d.groupby("consequence_class")]
    out = pd.DataFrame(rows).sort_values(
        "pct_gnomad_no_impact", ascending=False, na_position="last")

    print("BY CONSEQUENCE CLASS")
    print(out.to_string(index=False))

    tested = out[out.status == "tested"]
    print(f"\n{len(tested)} of {len(out)} classes had a comparison pool "
          f"(>= {NO_POOL_THRESHOLD} variants in each group).")

    import os
    os.makedirs(a.outdir, exist_ok=True)
    path = os.path.join(a.outdir, "gnomad_representation_by_consequence.tsv")
    out.to_csv(path, sep="\t", index=False)
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
