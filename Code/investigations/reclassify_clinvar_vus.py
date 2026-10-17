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
                 Needs: clnsig_norm, Summary_Consequence, anchor_tier,
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

# WHY DEPLETING MISSENSE VUS ARE NOT ROUTED TO A DISEASE ARM
# -----------------------------------------------------------
# Earlier versions mapped every depleting missense VUS through the CdLS
# calibration, on the grounds that missense is the predominant class of CdLS
# variants. That asserts an arm assignment the assay cannot make: DEE85 also
# produces missense variants (9 in the curated set), and the assay measures
# depletion, not mechanism.
#
# The assumption is not harmless. Routing the same 13 depleting missense VUS to
# DEE85 instead of CdLS moves 9 of them from "not carried forward" to PS3
# Strong -- the two arms' missense calibrations differ by a full tier. A
# conclusion that flips on an unsupported assumption should not be drawn.
#
# So by default they are reported as an observation and given no
# disease-specific tier. Depletion remains a specific signal in its own right:
# 0/290 benign in the matched reference and 0/11 ClinVar benign missense
# deplete. `--missense_policy route` restores the old behaviour for comparison.
#
# A NOTE ON THE ONE JUDGEMENT IN THIS SCRIPT
# -------------------------------------------
# `--min_k_path` is not a published rule. Brnich et al. 2019 requires a minimum
# of 11 total validation controls only "in the absence of rigorous statistical
# analysis", and lifts that explicitly once a formal OddsPath can be computed --
# which it can be here, against 290-1402 controls. So the framework offers no
# principled way to decline a tier resting on one or two pathogenic variants,
# and the CdLS depleting cells would otherwise be carried forward on 1-2
# variants with intervals spanning 2.1-2.8 orders of magnitude.
#
# The gate is therefore a stated judgement, kept visible as a parameter rather
# than buried in the logic, and every excluded variant records the counts and
# interval that caused the exclusion in its `note` column. --min_k_path 1
# reproduces the unfiltered Brnich reading.

import argparse
import sys

import pandas as pd

# Keyed on `Summary_Consequence`, the project's single consequence vocabulary
# (DECISIONS_LOG.md D156). The equivalent `Summary_Plot` column collapses
# Nonsense_Variant + Frameshift_Variant into a single `LOF` value, which is a
# misnomer -- a missense variant can also be loss-of-function, so what that
# label really means is "PTV". Both terms are listed explicitly here instead.
# The partition is identical either way; on the current ClinVar release the
# VUS set happens to contain no nonsense or frameshift variant at all, so the
# switch is results-neutral, but it will stay correct if a refresh adds one.
#
# Polypyrimidine-tract variants are deliberately NOT in the PTV/splice group.
# PVS1 does not reach them -- it covers nonsense, frameshift, canonical +/-1,2
# splice sites, the initiation codon and whole-exon deletions -- and, more to
# the point, the DEE85 truthset contains none: its 46 variants are 27 PTV,
# 10 canonical splice and 9 missense. Routing them to the DEE85 tier would
# extrapolate to a consequence class with no representation in the numerator,
# which is the same objection that keeps in-frame indels unmapped. Gene-wide
# they also behave unlike canonical splice variants (6.0% strongly depleting
# against 32.5%), and much more like missense (6.8%).
# `Summary_Consequence` pools canonical acceptor and donor sites with the
# surrounding splice region under one `Splice_Variant` label, and only 18% of
# that class is canonical. Membership of the PTV group is therefore taken from
# the VEP sequence ontology terms rather than from the summary class, matching
# `is_predicted_ptv` in `smc1a_lib.py`. PVS1 reaches canonical +/-1,2 sites and
# the DEE85 truthset's splice contribution is canonical donors, so a
# splice-region variant that is not a canonical site has no more claim on the
# DEE85 tier than a missense variant does, and is routed the same way.
PTV_SO_TERMS = {"stop_gained", "frameshift_variant", "splice_acceptor_variant",
                "splice_donor_variant", "start_lost"}
UNASSIGNABLE_CLASSES = {"Missense_Variant", "Splice_Variant"}
NOT_MAPPED_CLASSES = {"Inframe_Deletion", "Inframe_Insertion", "Synonymous_Variant",
                       "Intronic_Variant", "Others",
                       "Splice_Polypyrimidine_Tract_Variant"}

DEPLETED_TIERS = {"strongly depleting", "weakly depleting"}


def mechanism_class(summary: str, consequence: str) -> str:
    """Route a variant by the evidence its consequence is able to carry.

    `ptv_or_canonical_splice` can take the DEE85 tier. `arm_unassignable`
    covers consequence classes occurring in both disorders, for which no
    disease-specific tier applies however the variant behaves in the assay.
    """
    so = set((consequence or "").split(","))
    if PTV_SO_TERMS & so:
        return "ptv_or_canonical_splice"
    if summary in UNASSIGNABLE_CLASSES:
        return "arm_unassignable"
    return "not_mapped"


# Brnich et al. 2019 Table 3. The OddsPath point estimate is tiered against
# these, as the paper does.
PS3_THRESHOLDS = [("PS3 Very Strong", 350.0), ("PS3 Strong", 18.7),
                  ("PS3 Moderate", 4.3), ("PS3 Supporting", 2.1)]


def tier_of(oddspath) -> str:
    if oddspath is None or pd.isna(oddspath):
        return "indeterminate"
    for label, threshold in PS3_THRESHOLDS:
        if oddspath > threshold:
            return label
    return "indeterminate"


def load_tier_calibration(tsv: str, arm: str, reference: str) -> dict:
    """Per-tier calibration for one arm, keyed by assay call.

    Reads `10_calibrate_by_tier.py`'s output rather than the binary
    depleted-vs-not odds ratios. The binary figure pools strongly- and
    weakly-depleting calls against a near-empty benign comparator, which D154
    showed inflates the CdLS tier -- so mapping a VUS through it contradicts
    the by-tier reporting the calibration itself adopted. A VUS carries its own
    `anchor_tier`, so it can and should be mapped through the calibration for
    *that* tier.

    Each tier returns the significance gate (Fisher exact, correction-free) and
    the Brnich +1 OddsPath, which is what gets tiered.
    """
    df = pd.read_csv(tsv, sep="\t")
    sub = df[(df["arm"] == arm) & (df["reference"] == reference)
             & (df["mode"] == "by_tier")]
    if sub.empty:
        sys.exit(f"No by_tier rows in {tsv} for arm='{arm}' reference='{reference}'.\n"
                 "Available:\n" +
                 df[["arm", "reference", "mode"]].drop_duplicates().to_string(index=False))
    out = {}
    for _, r in sub.iterrows():
        out[r["result"]] = {
            "n_path": int(r["n_path"]), "k_path": int(r["k_path"]),
            "n_benign": int(r["n_benign"]), "k_benign": int(r["k_benign"]),
            "oddspath": r["LR_brnich"], "ci_lo": r["ci_lo"], "ci_hi": r["ci_hi"],
            "fisher_p": r["fisher_p"], "arm": arm,
        }
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary_tsv", required=True)
    ap.add_argument("--ptv_splice_tsv", required=True,
                    help="calibration_by_tier_*.tsv for the PTV-inclusive cut, "
                         "used for PTV/splice VUS")
    ap.add_argument("--missense_tsv", required=True,
                    help="calibration_by_tier_*.tsv for the matched "
                         "(missense/in-frame) cut, used for missense VUS")
    ap.add_argument("--ptv_splice_arm", default="DEE85_pathogenic")
    ap.add_argument("--missense_policy", choices=["observe", "route"],
                    default="observe",
                    help="what to do with a depleting missense VUS. 'observe' "
                         "(default) reports it as an observation and applies no "
                         "disease-specific tier, because the assay cannot say "
                         "which arm a missense variant belongs to. 'route' maps "
                         "it through --missense_arm, which asserts a mechanism "
                         "the data does not support -- see the note below.")
    ap.add_argument("--missense_arm", default="CdLS_pathogenic",
                    help="only used with --missense_policy route")
    ap.add_argument("--reference", default="curated+ClinVar+gnomAD")
    ap.add_argument("--alpha", type=float, default=0.05,
                    help="significance gate on the Fisher exact p of the "
                         "calibration cell a VUS is mapped through")
    ap.add_argument("--min_k_path", type=int, default=3,
                    help="minimum pathogenic variants in the calibration cell "
                         "before its tier is carried forward. Encodes a "
                         "JUDGEMENT, not a published rule -- see the note at "
                         "the top of this file and section 3.2 of the "
                         "calibration README. Set to 1 to disable.")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    ptv_cal = load_tier_calibration(args.ptv_splice_tsv, args.ptv_splice_arm,
                                    args.reference)
    mis_cal = load_tier_calibration(args.missense_tsv, args.missense_arm,
                                    args.reference)

    for label, cal in [("PTV/splice -> " + args.ptv_splice_arm, ptv_cal),
                       ("missense   -> " + args.missense_arm, mis_cal)]:
        print(f"{label}, by tier:", file=sys.stderr)
        for t in ("strongly depleting", "weakly depleting"):
            c = cal.get(t)
            if c is None:
                print(f"    {t:<20} absent from the calibration table", file=sys.stderr)
                continue
            if pd.isna(c["fisher_p"]) or c["fisher_p"] >= args.alpha:
                sig = "   ** not significant -- not carried forward **"
            elif c["k_path"] < args.min_k_path:
                sig = (f"   ** significant but only {c['k_path']} pathogenic "
                       f"variant(s) -- not carried forward **")
            else:
                sig = ""
            op = "--" if pd.isna(c["oddspath"]) else f"{c['oddspath']:.2f}"
            print(f"    {t:<20} {c['k_path']}/{c['n_path']} vs "
                  f"{c['k_benign']}/{c['n_benign']}  OddsPath {op:>8}  "
                  f"p={c['fisher_p']:.3g}  -> {tier_of(c['oddspath'])}{sig}",
                  file=sys.stderr)

    df = pd.read_csv(args.summary_tsv, sep="\t", dtype=str)
    vus = df[df["clnsig_norm"] == "Uncertain significance"].copy()
    print(f"\n{len(vus)} VUS row(s) in {args.summary_tsv}", file=sys.stderr)

    def classify(row):
        mclass = mechanism_class(row["Summary_Consequence"],
                                 row.get("Consequence", ""))
        tier = row["anchor_tier"]
        if tier == "enriched":
            return mclass, "manual_review", "enriched call, not evaluated by this framework"
        if tier == "no impact":
            return (mclass, "no_recommendation",
                    "a single non-depleting call is not safe BS3 evidence")
        if tier not in DEPLETED_TIERS:
            return mclass, "manual_review", f"unrecognised anchor_tier value: {tier!r}"
        if mclass == "arm_unassignable" and args.missense_policy == "observe":
            # No disease-specific tier. Both SMC1A disorders produce missense
            # variants, the assay measures depletion rather than mechanism, and
            # the two arms' calibrations differ by a full tier -- so routing
            # would assert an arm assignment the data cannot support. The
            # depletion itself is still a real and specific observation.
            return (mclass, "observation_functional_abnormality",
                    f"{tier}; benign variants of this class essentially never "
                    f"deplete (0/290 in the matched benign reference), so this "
                    f"is evidence of functional abnormality. No "
                    f"disease-specific PS3 tier is applied, the arm being "
                    f"unassignable from consequence class alone.")
        if mclass == "ptv_or_canonical_splice":
            cal = ptv_cal
        elif mclass == "arm_unassignable":
            cal = mis_cal
        else:
            return (mclass, "manual_review",
                    "consequence class not mapped to either calibration")
        c = cal.get(tier)
        if c is None or pd.isna(c["oddspath"]):
            return (mclass, "not_carried_forward",
                    f"no pathogenic variant occupies the {tier} tier in "
                    f"{c['arm'] if c else 'the'} calibration")
        # Significance gate first, then the tier. A cell that does not clear it
        # is reported, not applied -- see the calibration README section 3.2.
        if pd.isna(c["fisher_p"]) or c["fisher_p"] >= args.alpha:
            return (mclass, "not_carried_forward",
                    f"{c['arm']} {tier}: {c['k_path']}/{c['n_path']} vs "
                    f"{c['k_benign']}/{c['n_benign']}, Fisher p={c['fisher_p']:.3g} "
                    f"-- not distinguishable from no association")
        if c["k_path"] < args.min_k_path:
            return (mclass, "not_carried_forward",
                    f"{c['arm']} {tier}: significant (p={c['fisher_p']:.3g}, "
                    f"OddsPath {c['oddspath']:.2f}) but rests on only "
                    f"{c['k_path']} pathogenic variant(s); interval "
                    f"{c['ci_lo']:.2f}-{c['ci_hi']:.1f}. Reported, not applied "
                    f"-- see calibration README 3.2")
        note = (f"{c['arm']} {tier}: OddsPath {c['oddspath']:.2f} "
                f"(95% CI {c['ci_lo']:.2f}-{c['ci_hi']:.1f}, p={c['fisher_p']:.3g})")
        return mclass, tier_of(c["oddspath"]), note

    results = vus.apply(classify, axis=1, result_type="expand")
    results.columns = ["mechanism_class", "recommendation", "note"]
    out = pd.concat([
        vus[["Targeton_ID", "position", "HGVSc", "HGVSp", "Summary_Consequence",
             "anchor_tier", "pos_adj_log2FoldChange_raw", "clnsig_norm", "clinvar_id"]],
        results,
    ], axis=1)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"\nWrote {args.output} ({len(out)} rows)", file=sys.stderr)

    print("\n=== mechanism_class x anchor_tier ===", file=sys.stderr)
    print(pd.crosstab(out["mechanism_class"], out["anchor_tier"]).to_string(),
          file=sys.stderr)
    print("\n=== recommendation counts ===", file=sys.stderr)
    print(out["recommendation"].value_counts().to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
