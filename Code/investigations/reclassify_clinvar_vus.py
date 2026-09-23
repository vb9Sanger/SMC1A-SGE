#!/usr/bin/env python3
"""
reclassify_clinvar_vus.py

Applies the calibrated sensitivity/OddsPath evidence strengths
(09_calibrate_sensitivity_oddspath.py) to ClinVar Variants of Uncertain
Significance (VUS), to see which ones the SGE assay's own depletion data
would support moving toward a PS3/BS3 call, and which it can't speak to.

The mapping is deliberately conservative, not automatic-for-everything:

  - PTV/splice-class VUS (LOF, Splice_Variant, Splice_Polypyrimidine_Tract_
    Variant) that deplete -> the DEE85_pathogenic-calibrated OddsPath tier,
    read from the FULL (consequence-unrestricted) calibration run. This is
    the only option available for this class: a matched PTV/splice-only
    benign reference is not constructible for this gene -- ClinVar has zero
    variants classified Benign/Likely benign in the actual PTV consequence
    classes (nonsense, frameshift, splice_acceptor, splice_donor), for the
    biologically expected reason that a true truncating variant essentially
    never gets classified benign in a haploinsufficiency-sensitive gene.
    Checked directly (2026-09-23): restricting the benign reference to
    these classes collapses it to 0 ClinVar rows and 7 curated rows, all
    the softer `splice_region` class rather than a true PTV, far too few to
    calibrate. So the unrestricted comparison is used here of necessity,
    not as an unexamined shortcut -- and it also happens to be a good match
    already, since DEE85_pathogenic's own calibration group is itself ~83%
    frameshift/nonsense/splice (Brnich et al. 2019's recommendation that
    "controls should also be relevant to... the type of variant under
    consideration" is satisfied on the pathogenic side even though the
    benign side can't be restricted the same way).
  - Missense VUS that deplete -> the CdLS_pathogenic-calibrated OddsPath
    tier, read from a SEPARATE, missense/in-frame-restricted calibration
    run (--missense_or_tsv), not the full-group run. Brnich et al. 2019 (the
    OddsPath framework this whole calibration follows) explicitly
    recommends matching reference-control consequence class to the variant
    class under interpretation ("missense controls for evaluating missense
    variants of uncertain significance") -- unlike the PTV/splice case, a
    matched missense/in-frame benign reference *is* constructible here
    (n=14: 3 curated + 11 ClinVar), so the matched figure is used rather
    than the full-group one. That comparison is markedly weaker than the
    full-group figure (OR=2.6, 95% CI 0.1-54.2, not significant, vs the
    full-group OR=26.1) -- reported plainly rather than silently upgraded,
    since misrepresenting an unmatched, more optimistic figure as if it
    were the properly-matched one would be exactly the confound Brnich et
    al.'s recommendation is meant to avoid.
  - Inframe deletion/insertion VUS -> NOT mapped to either tier. They were
    lumped into "not missense" in earlier quick cuts of the truth set, but
    mechanistically they don't clearly belong with true PTVs (frame is
    preserved) or with missense either. Flagged for manual review, not
    silently assigned to whichever bucket happens to be more convenient.
  - Any other consequence class (synonymous, intronic, "Others") that
    depletes -> also flagged for manual review, not mapped. A depleting
    synonymous/intronic VUS is exactly the kind of case the exon-level
    case-studies document investigates individually (e.g. splice-adjacent
    positional effects) -- it doesn't fit either calibrated mechanism class
    and deserves the same individual scrutiny, not a mechanical tier.
  - "no impact" VUS, any class -> NOT given a BS3 recommendation. Even
    DEE85's sensitivity (82.6%) is imperfect, and CdLS's (16.7%) means a
    non-depleting missense variant is common even among true CdLS
    pathogenic variants -- a single non-depleting call is not by itself
    safe BS3 evidence, per the write-up's own ACMG guidance (Section 3).
  - "enriched" VUS -> flagged for manual review; not evaluated by this
    calibration framework at all (enriched wasn't part of the P/LP-vs-
    benign comparisons the OddsPath tiers were computed from).

Inputs
------
--summary_tsv    clinvar_variants_summary.tsv (sge_clinvar_intersect.py).
                 Needs: clnsig_norm, Summary_Plot, anchor_tier,
                 pos_adj_log2FoldChange_raw, HGVSp, variant identifying
                 columns.
--or_tsv         odds_ratios.tsv from 09_calibrate_sensitivity_oddspath.py's
                 FULL (consequence-unrestricted) run. Used for the
                 PTV/splice mapping.
--missense_or_tsv
                 odds_ratios.tsv from a separate run of
                 09_calibrate_sensitivity_oddspath.py with
                 `--consequence_filter missense inframe_deletion
                 inframe_insertion`. Used for the missense mapping. Defaults
                 to --or_tsv if not given (not recommended -- see docstring
                 above on why the two should differ).
--ptv_splice_comparison
                 which row of --or_tsv to use for PTV/splice-class VUS,
                 given as "GROUP_A|GROUP_B" matching that row's group_a/
                 group_b columns exactly (order matters: the pathogenic
                 group should be group_a). Example:
                 "DEE85_pathogenic|Combined_Benign"
--missense_comparison
                 same, but a row of --missense_or_tsv, for missense-class
                 VUS. Example: "CdLS_pathogenic|Combined_Benign"
--output

Usage
-----
    python reclassify_clinvar_vus.py \\
        --summary_tsv clinvar_variants_summary.tsv \\
        --or_tsv sensitivity_results/odds_ratios.tsv \\
        --missense_or_tsv sensitivity_missense_only/odds_ratios.tsv \\
        --ptv_splice_comparison "DEE85_pathogenic|Combined_Benign" \\
        --missense_comparison "CdLS_pathogenic|Combined_Benign" \\
        --output vus_reclassification.tsv
"""

import argparse
import sys

import pandas as pd

PTV_SPLICE_CLASSES = {"LOF", "Splice_Variant", "Splice_Polypyrimidine_Tract_Variant"}
MISSENSE_CLASSES = {"Missense_Variant"}
NOT_MAPPED_CLASSES = {"Inframe_Deletion", "Inframe_Insertion", "Synonymous_Variant",
                       "Intronic_Variant", "Others"}

DEPLETED_TIERS = {"strongly depleting", "weakly depleting"}


def mechanism_class(consequence: str) -> str:
    if consequence in PTV_SPLICE_CLASSES:
        return "ptv_or_splice"
    if consequence in MISSENSE_CLASSES:
        return "missense"
    return "not_mapped"


def load_comparison(or_tsv: str, spec: str) -> dict:
    group_a, group_b = spec.split("|", 1)
    df = pd.read_csv(or_tsv, sep="\t")
    row = df[(df["group_a"] == group_a) & (df["group_b"] == group_b)]
    if row.empty:
        sys.exit(f"No row in {or_tsv} with group_a='{group_a}' group_b='{group_b}'. "
                  f"Available pairs:\n" +
                  df[["group_a", "group_b"]].to_string(index=False))
    r = row.iloc[0]
    return {
        "comparison": f"{group_a} vs {group_b}",
        "OR": r["OR"], "ci_lo": r["ci_lo"], "ci_hi": r["ci_hi"],
        "p": r["fisher_p"], "tier": r["oddspath_tier_point_estimate"],
        "ci_crosses_one": r["ci_lo"] <= 1.0 <= r["ci_hi"],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary_tsv", required=True)
    ap.add_argument("--or_tsv", required=True,
                     help="Full (consequence-unrestricted) odds_ratios.tsv, for PTV/splice VUS")
    ap.add_argument("--missense_or_tsv", default=None,
                     help="Missense/in-frame-restricted odds_ratios.tsv, for missense VUS "
                          "(defaults to --or_tsv if omitted -- not recommended, see docstring)")
    ap.add_argument("--ptv_splice_comparison", required=True, metavar="GROUP_A|GROUP_B")
    ap.add_argument("--missense_comparison", required=True, metavar="GROUP_A|GROUP_B")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    missense_or_tsv = args.missense_or_tsv or args.or_tsv
    if args.missense_or_tsv is None:
        print("WARNING: --missense_or_tsv not given, falling back to --or_tsv for the "
              "missense mapping too. This uses an unmatched (consequence-unrestricted) "
              "benign reference for missense VUS, which Brnich et al. 2019 recommends "
              "against -- pass --missense_or_tsv explicitly.", file=sys.stderr)

    ptv_cal = load_comparison(args.or_tsv, args.ptv_splice_comparison)
    missense_cal = load_comparison(missense_or_tsv, args.missense_comparison)

    for label, cal in [("PTV/splice", ptv_cal), ("missense", missense_cal)]:
        flag = " ** 95% CI CROSSES 1 -- do not treat as reliable evidence **" if cal["ci_crosses_one"] else ""
        print(f"{label} calibration: {cal['comparison']}  OR={cal['OR']:.1f} "
              f"(95% CI {cal['ci_lo']:.1f}-{cal['ci_hi']:.1f})  tier={cal['tier']}{flag}",
              file=sys.stderr)

    df = pd.read_csv(args.summary_tsv, sep="\t", dtype=str)
    vus = df[df["clnsig_norm"] == "Uncertain significance"].copy()
    print(f"\n{len(vus)} VUS row(s) in {args.summary_tsv}", file=sys.stderr)

    def classify(row):
        mclass = mechanism_class(row["Summary_Plot"])
        tier = row["anchor_tier"]
        if tier == "enriched":
            return mclass, "manual_review", "enriched call, not evaluated by this framework"
        if tier == "no impact":
            return mclass, "no_recommendation", "single non-depleting call is not safe BS3 evidence"
        if tier not in DEPLETED_TIERS:
            return mclass, "manual_review", f"unrecognised anchor_tier value: {tier!r}"
        # depleted (weakly or strongly)
        if mclass == "ptv_or_splice":
            cal = ptv_cal
        elif mclass == "missense":
            cal = missense_cal
        else:
            return mclass, "manual_review", "consequence class not mapped to either calibration (see docstring)"
        reliability = " (CI crosses 1 -- treat cautiously)" if cal["ci_crosses_one"] else ""
        return mclass, f"PS3 {cal['tier']}{reliability}", f"borrowed from {cal['comparison']} (OR={cal['OR']:.1f})"

    results = vus.apply(classify, axis=1, result_type="expand")
    results.columns = ["mechanism_class", "recommendation", "note"]
    out = pd.concat([
        vus[["Targeton_ID", "position", "HGVSc", "HGVSp", "Summary_Plot",
             "anchor_tier", "pos_adj_log2FoldChange_raw", "clnsig_norm", "clinvar_id"]],
        results,
    ], axis=1)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"\nWrote {args.output} ({len(out)} rows)", file=sys.stderr)

    print("\n=== mechanism_class x anchor_tier ===", file=sys.stderr)
    print(pd.crosstab(out["mechanism_class"], out["anchor_tier"]).to_string(), file=sys.stderr)
    print("\n=== recommendation counts ===", file=sys.stderr)
    print(out["recommendation"].value_counts().to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
