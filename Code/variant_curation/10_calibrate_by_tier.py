#!/usr/bin/env python3
"""
10_calibrate_by_tier.py

Two extensions to 09_calibrate_sensitivity_oddspath.py, run together because
they interact.

**1. Benign reference composition.** The curated benign group is n=25 with
zero false positives, which makes LR+ depend almost entirely on a
Haldane-Anscombe correction rather than on data. This script recomputes the
calibration against four reference sets so the effect of composition is
visible rather than assumed:

    curated                      the original Exclude_Benign group
    curated + ClinVar            ClinVar B/LB, any review tier
    curated + gnomAD             population frequency only, no clinical claim
    curated + ClinVar + gnomAD   everything

The gnomAD-only variant matters because it carries no clinical assertion at
all: it is evidence against pathogenicity from allele frequency in a severe,
de novo, fitness-reducing disorder, and it is independent of whether anyone
classified the variant correctly. Running it with and without ClinVar answers
directly how much the result depends on trusting ClinVar.

**2. Calibration by tier, not by a binary call.** `09_...py` collapses
`anchor_tier` to depleted / not depleted, which puts "weakly depleting" and
"strongly depleting" in one bucket. They are not equivalent evidence: a weak
call sits near the anchor and carries real uncertainty about whether anything
happened. Pooling them lets the confident calls carry the uncertain ones.

This computes a separate likelihood ratio for each observable result:

    LR(tier) = P(tier | pathogenic) / P(tier | benign)

which is what Brnich et al. 2019's OddsPath reduces to for a given result --
the binary LR+ is the special case where the result is "depleted". Each tier
gets its own evidence strength, on the PS3 side where LR > 1 and the BS3 side
where LR < 1, so "no impact" is tiered as evidence AGAINST pathogenicity
rather than merely being the absence of evidence for it.
"""
import argparse
import math
import re

import pandas as pd

ARMS = ("CdLS_pathogenic", "DEE85_pathogenic")
TIERS = ("strongly depleting", "weakly depleting", "no impact")

# Brnich et al. 2019 Table 3 / ClinGen SVI. The BS3 thresholds are the exact
# reciprocals of the PS3 ones: the same Bayes factor read in the other
# direction.
PS3 = ((350, "Very Strong"), (18.7, "Strong"), (4.33, "Moderate"), (2.08, "Supporting"))
BS3 = ((1 / 350, "Very Strong"), (0.053, "Strong"), (0.23, "Moderate"), (0.48, "Supporting"))


def tier_of(lr):
    """Map a likelihood ratio to an ACMG evidence strength, either direction."""
    if lr is None or not math.isfinite(lr):
        return "n/a"
    for t, name in PS3:
        if lr >= t:
            return f"PS3 {name}"
    for t, name in BS3:
        if lr <= t:
            return f"BS3 {name}"
    return "indeterminate"


def lr_brnich(a, n1, b, n2):
    """OddsPath with the zero-cell handling Brnich et al. 2019 actually use.

    The paper addresses this directly in its Methods: "To circumvent posterior
    probabilities of zero or infinity, and to account for the possibility that
    the next variant tested in the assay might have a discordant result, we
    added exactly one misclassified variant to each set."

    So: one pathogenic control that reads normal, and one benign control that
    reads abnormal. That is a +1 correction applied to the discordant cell of
    each set, and it is materially more conservative than the +0.5
    Haldane-Anscombe correction applied to every cell -- for the matched
    DEE85 comparison it is the difference between Strong and Moderate.

    Applied unconditionally, as the paper does, not only when a cell is zero.
    """
    if n1 + 1 == 0 or n2 + 1 == 0:
        return None
    # As in lr_ci: with no pathogenic variant in this tier there is nothing to
    # estimate, and the formula returns a flat 0.00 that would be labelled
    # "BS3 Very Strong" for a result neither group exhibits.
    if a == 0:
        return None
    return (a / (n1 + 1)) / ((b + 1) / (n2 + 1))


def lr_ci(a, n1, b, n2, z=1.96):
    """LR = (a/n1)/(b/n2) with a log-scale 95% CI (Simel et al. 1991).

    Haldane-Anscombe: 0.5 is added to every cell only when a zero would make
    the ratio or its variance undefined, and the correction is reported so a
    reader can see which figures rest on it.
    """
    # No pathogenic variant landed in this tier, so there is nothing to
    # estimate. A Haldane correction here produces a large LR purely from
    # n2 >> n1 -- with 0/11 pathogenic and 0/291 benign it reports 24.08 and
    # a "PS3 Strong" label for a tier no pathogenic variant occupies. Refuse
    # rather than print it.
    if a == 0:
        return None, None, None, False
    corrected = (a == 0 or b == 0)
    if corrected:
        a, n1, b, n2 = a + 0.5, n1 + 1, b + 0.5, n2 + 1
    if b == 0 or n1 == 0 or n2 == 0:
        return None, None, None, corrected
    lr = (a / n1) / (b / n2)
    try:
        se = math.sqrt(1 / a - 1 / n1 + 1 / b - 1 / n2)
    except (ValueError, ZeroDivisionError):
        return lr, None, None, corrected
    return lr, lr * math.exp(-z * se), lr * math.exp(z * se), corrected


def fmt(v):
    """A likelihood ratio, or a dash where none is estimable."""
    return "-" if v is None or pd.isna(v) else f"{v:.2f}"


def norm_key(k):
    """Normalise a variant_key to `chrX:pos:ref:alt`.

    The three benign sources do not agree on format:
    `assay_join_all.tsv` and the gnomAD extractor write `chrX:53380084:C:A`,
    while `extract_clinvar_benign_controls.py` writes `53380118:G:A` with no
    contig. Without this, no key from ClinVar can ever match one from the
    other two, so the "combined" reference silently counts shared variants
    twice -- 13 between curated and ClinVar, 73 between ClinVar and gnomAD.
    """
    k = str(k).strip()
    if k.startswith("chrX:"):
        return k
    if re.match(r"^\d+:", k):
        return "chrX:" + k
    return k


def load_extra(path):
    d = pd.read_csv(path, sep="\t", dtype=str)
    if "anchor_tier" not in d.columns:
        # older files carry only the binary call; treat a depleted call as
        # unresolved between the two depleting tiers rather than guessing
        d["anchor_tier"] = d["anchor_call"].map(
            {"no impact": "no impact"}).fillna("depleted (tier unknown)")
    d["variant_key"] = d["variant_key"].map(norm_key)
    return d[["variant_key", "anchor_tier"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--join_tsv", required=True)
    ap.add_argument("--clinvar_benign", required=True)
    ap.add_argument("--gnomad_benign", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--consequence_filter", nargs="*", default=None,
                    help="restrict BOTH the pathogenic arms and the benign "
                         "reference to these consequence_class values. Brnich "
                         "et al. 2019 ask that controls match the variant "
                         "class under interpretation, so the restriction must "
                         "apply to both sides or the comparison is not matched.")
    ap.add_argument("--consequence_map", nargs="*", default=[],
                    help="extra TSVs supplying a consequence for benign "
                         "controls that are not in the curated join, e.g. "
                         "clinvar_variants_summary.tsv and "
                         "sge_gnomad_summary.tsv. Both carry oligo_name and "
                         "Summary_Consequence.")
    ap.add_argument("--label", default="",
                    help="tag written into the output, e.g. 'excl_PTV'")
    a = ap.parse_args()

    j = pd.read_csv(a.join_tsv, sep="\t", dtype=str, low_memory=False)
    j = j.drop_duplicates("variant_key")
    j["variant_key"] = j["variant_key"].map(norm_key)

    # Consequence class for every variant, so the benign sources -- which
    # carry only a key and a tier -- can be restricted on the same basis as
    # the arms. A control whose class is unknown is dropped rather than
    # assumed to match, so a matched run never silently includes an unmatched
    # control.
    cls_of = dict(zip(j.variant_key, j.consequence_class))

    # Most benign controls are NOT in the curated join, so they have no
    # consequence_class there and a matched run would silently drop them --
    # the first attempt collapsed the matched reference from 14 to 3 this way.
    # ClinVar and gnomAD both carry Summary_Consequence keyed by oligo_name,
    # in a different vocabulary, so it is mapped onto the curated one.
    SUMMARY_TO_CLASS = {
        "Missense_Variant": "missense",
        "Synonymous_Variant": "synonymous",
        "Nonsense_Variant": "nonsense",
        "Frameshift_Variant": "frameshift",
        "Inframe_Deletion": "inframe_deletion",
        "Inframe_Insertion": "inframe_insertion",
        "Intronic_Variant": "intronic",
        "Splice_Variant": "splice_region",
        "Splice_Polypyrimidine_Tract_Variant": "splice_region",
    }
    OLIGO_KEY = re.compile(r"(chrX):(\d+)_([ACGT]+)>([ACGT]+)")
    for path in a.consequence_map:
        m = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        if "oligo_name" not in m.columns or "Summary_Consequence" not in m.columns:
            continue
        e = m.oligo_name.str.extract(OLIGO_KEY)
        keys = (e[0] + ":" + e[1] + ":" + e[2] + ":" + e[3])
        for k, sc in zip(keys, m.Summary_Consequence):
            if pd.isna(k) or k in cls_of:
                continue                      # the curated class wins
            c = SUMMARY_TO_CLASS.get(str(sc))
            if c:
                cls_of[k] = c
    if a.consequence_filter:
        keep_cls = set(a.consequence_filter)
        print(f"CONSEQUENCE FILTER: {sorted(keep_cls)} "
              f"(applied to both arms and controls)\n")
        j = j[j.consequence_class.isin(keep_cls)]
    else:
        keep_cls = None
    cur_ben = j[j.curation_group == "Exclude_Benign"][["variant_key", "anchor_tier"]]
    cv_ben = load_extra(a.clinvar_benign)
    gn_ben = load_extra(a.gnomad_benign)

    def union(*frames):
        u = pd.concat(frames).drop_duplicates("variant_key")
        u = u[u.anchor_tier.isin(TIERS)]
        if keep_cls is not None:
            u = u[u.variant_key.map(lambda k: cls_of.get(k)).isin(keep_cls)]
        return u

    REFS = {
        "curated": union(cur_ben),
        "curated+ClinVar": union(cur_ben, cv_ben),
        "curated+gnomAD": union(cur_ben, gn_ben),
        "curated+ClinVar+gnomAD": union(cur_ben, cv_ben, gn_ben),
    }

    rows = []
    print("=" * 96)
    print("BENIGN REFERENCE COMPOSITION")
    print("=" * 96)
    for name, d in REFS.items():
        vc = d.anchor_tier.value_counts()
        dep = int(vc.get("strongly depleting", 0) + vc.get("weakly depleting", 0))
        print(f"  {name:<26} n={len(d):<5} strong={vc.get('strongly depleting',0):<4} "
              f"weak={vc.get('weakly depleting',0):<4} none={vc.get('no impact',0):<5} "
              f"specificity={1-dep/len(d):.4f}")

    for arm in ARMS:
        p = j[(j.curation_group == arm) & (j.anchor_tier.isin(TIERS))]
        n1 = len(p)
        pc = p.anchor_tier.value_counts()
        for ref_name, ben in REFS.items():
            n2 = len(ben)
            bc = ben.anchor_tier.value_counts()

            # --- binary, for comparability with 09_...py
            a_dep = int(pc.get("strongly depleting", 0) + pc.get("weakly depleting", 0))
            b_dep = int(bc.get("strongly depleting", 0) + bc.get("weakly depleting", 0))
            lr, lo, hi, corr = lr_ci(a_dep, n1, b_dep, n2)
            rows.append(dict(label=a.label, arm=arm, reference=ref_name, mode="binary",
                             result="depleted", n_path=n1, k_path=a_dep,
                             n_benign=n2, k_benign=b_dep, LR=lr, ci_lo=lo,
                             ci_hi=hi, haldane=corr, evidence=tier_of(lr),
                             LR_brnich=lr_brnich(a_dep, n1, b_dep, n2),
                             evidence_brnich=tier_of(lr_brnich(a_dep, n1, b_dep, n2))))

            # --- by tier
            for t in TIERS:
                ka, kb = int(pc.get(t, 0)), int(bc.get(t, 0))
                lr, lo, hi, corr = lr_ci(ka, n1, kb, n2)
                rows.append(dict(label=a.label, arm=arm, reference=ref_name, mode="by_tier",
                                 result=t, n_path=n1, k_path=ka, n_benign=n2,
                                 k_benign=kb, LR=lr, ci_lo=lo, ci_hi=hi,
                                 haldane=corr, evidence=tier_of(lr),
                                 LR_brnich=lr_brnich(ka, n1, kb, n2),
                                 evidence_brnich=tier_of(lr_brnich(ka, n1, kb, n2))))

    out = pd.DataFrame(rows)
    suffix = f"_{a.label}" if a.label else ""
    out.to_csv(f"{a.outdir}/calibration_by_tier{suffix}.tsv", sep="\t", index=False)

    for mode, label in (("binary", "BINARY (depleted vs not) -- comparable to 09_...py"),
                        ("by_tier", "BY TIER -- separate evidence per observable result")):
        print()
        print("=" * 96)
        print(label)
        print("=" * 96)
        s = out[out["mode"] == mode]
        for arm in ARMS:
            print(f"\n### {arm}")
            print(f"  {'reference':<26}{'result':<20}{'path':<10}{'benign':<11}"
                  f"{'LR':>9}  {'95% CI':<22}{'evidence (Haldane)':<18}"
                  f"{'LR':>8}  {'evidence (Brnich +1)'}")
            for _, r in s[s.arm == arm].iterrows():
                ci = (f"{r.ci_lo:.2f}-{r.ci_hi:.1f}"
                      if pd.notna(r.ci_lo) and pd.notna(r.ci_hi) else "-")
                star = "*" if r.haldane else " "
                print(f"  {r.reference:<26}{r.result:<20}"
                      f"{str(r.k_path)+'/'+str(r.n_path):<10}"
                      f"{str(r.k_benign)+'/'+str(r.n_benign):<11}"
                      f"{fmt(r.LR):>8}{star} {ci:<22}{r.evidence:<18}"
                      f"{fmt(r.LR_brnich):>8}  {r.evidence_brnich}")
    print("\n  * LR rests on a Haldane-Anscombe correction for a zero cell.")
    print("  The right-hand pair applies instead the correction Brnich et al. 2019")
    print("  themselves use -- one misclassified variant added to each set -- which")
    print("  is more conservative and is applied unconditionally, as in the paper.")
    print(f"\nWrote {a.outdir}/calibration_by_tier{suffix}.tsv")


if __name__ == "__main__":
    main()
