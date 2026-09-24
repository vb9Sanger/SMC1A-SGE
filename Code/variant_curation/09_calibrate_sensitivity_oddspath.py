#!/usr/bin/env python3
"""
09_calibrate_sensitivity_oddspath.py

Formal sensitivity/specificity/OddsPath calibration of the SGE assay against
the curated truth set, following the Brnich et al. 2019 / ClinGen SVI
OddsPath framework -- the item STEP SEVEN's README section previously listed
as "planned, not yet done".

Input is the output of `06_join_assay.py` (assay_join_all.tsv): one row per
matched (curated variant, oligo) pair, with `curation_group`, `anchor_call`,
`variant_key`, and `match_type` columns.

Question this answers: given the truth set's own CdLS/DEE85 mechanism split
(missense/dominant-negative vs protein-truncating/loss-of-function, see
STEP SEVEN background), how well does the assay's depletion call track
disease-causing status -- and does that differ by mechanism group? Computes,
per curation group (CdLS_pathogenic, DEE85_pathogenic, Exclude_Benign by
default -- override with --groups):
  - sensitivity/specificity with Clopper-Pearson 95% CIs
  - for every pair of groups supplied, with group_a treated as the
    pathogenic/case reference and group_b as the benign/control reference:
      * OddsPath itself (Brnich et al. 2019's own quantity, not an
        approximation of it -- see note below), separately for the two
        results their Table 3 actually calibrates against:
          - LR+ = sensitivity_a / (1 - specificity_b), for an ABNORMAL
            (depleting) result -> compared against the PS3-side thresholds
            (>=2.08 / >=4.33 / >=18.7 / >=350)
          - LR- = (1 - sensitivity_a) / specificity_b, for a NORMAL
            (non-depleting) result -> compared against the BS3-side
            thresholds (<=0.48 / <=0.23 / <=0.053 / <=1/350), the exact
            reciprocals of the PS3 ones
        each with a log-scale 95% CI (Simel et al. 1991) and Fisher's exact
        p on the underlying 2x2 table
      * the case-control/diagnostic odds ratio (DOR = LR+ / LR-), reported
        alongside as a supplementary omnibus statistic, NOT tiered against
        Table 3 -- see note below on why it isn't the same quantity

  Why LR+/LR- rather than the pooled odds ratio: Brnich's OddsPath is
  literally a Bayes factor, OddsPath = posterior odds / prior odds. Writing
  the posterior via Bayes' rule from sensitivity/specificity and a prior
  P(pathogenic), the prior cancels out algebraically and OddsPath for an
  abnormal result reduces exactly to LR+ = sensitivity/(1-specificity); for
  a normal result it reduces to LR- = (1-sensitivity)/specificity. The
  case-control odds ratio computed from the full 2x2 table is a different,
  related quantity -- the diagnostic odds ratio, DOR = LR+/LR- -- which
  conflates the abnormal-result and normal-result evidence into one pooled
  number instead of treating them separately, the way Table 3's own
  (non-reciprocal-looking, because they ARE reciprocals of each other:
  0.48=1/2.08, 0.23=1/4.33, 0.053=1/18.7) PS3/BS3 thresholds require. An
  earlier version of this script tiered the pooled DOR directly against
  Table 3, which is only a good approximation of the true OddsPath when
  specificity is high (LR- close to 1); this version computes the paper's
  own quantities directly instead of relying on that approximation holding.

De-duplicates by variant_key within each group before computing anything
(a variant matched by more than one oligo -- see the data-quality note
below -- would otherwise be double-counted or silently resolved by whichever
row happens to sort first).

Note on `position_only_verify_manually` rows: assay_join_all.tsv's fallback
match (06_join_assay.py, load_and_match(), the `pos_only` merge) joins on
`span_start` alone, not on span length/end, by design -- allele agreement
can't be established automatically for these, so the label exists to route
them to manual review rather than trusting a position-only match silently.
Any variant carrying this label should be manually verified against its
matched oligo(s) before being used downstream. One such case was checked
this way for the truth set's CdLS/DEE85/Benign groups (chrX:53382339:GCCAT:G,
c.3326_3329del): it had matched two oligos on position, one a genuine allele
match, the other an unrelated different-length edit -- the spurious
`..._53382339_53382341_inframe` row was removed from all three
assay_join_all.tsv copies once verified. A full check of the other 25 indels
and all other multi-oligo variant_keys in these groups confirmed this was
the only row carrying the position_only_verify_manually label in that set.
This script's de-duplication (`--on_duplicate`) makes the resolution choice
explicit and auditable for any future such case, rather than silent.

Usage:
    python 09_calibrate_sensitivity_oddspath.py \\
        --join_tsv assay_join_all.tsv \\
        --outdir sensitivity_results/ \\
        --pathogenic_groups CdLS_pathogenic DEE85_pathogenic \\
        --benign_group Exclude_Benign

Additional benign controls (e.g. from ClinVar -- see
extract_clinvar_benign_controls.py, which produces exactly the (variant_key,
anchor_call) TSV shape --extra_group_tsv expects) can be added as new named
groups without touching assay_join_all.tsv itself:
    --extra_group_tsv "ClinVar_Benign=clinvar_benign_controls.tsv"

Pooled comparisons (e.g. "how strong is the evidence if CdLS and DEE85 are
combined into one P/LP group", or "curated + ClinVar benign combined") are
computed from the already-counted groups, not by re-reading files, so both
the pooled and unpooled numbers are always reported side by side rather than
the pooled figure silently replacing the unpooled ones:
    --pool "Combined_PLP=CdLS_pathogenic,DEE85_pathogenic" \\
    --pool "Combined_Benign=Exclude_Benign,ClinVar_Benign"
"""

import argparse
import math
import os
import sys
from collections import defaultdict
from itertools import combinations

import pandas as pd
from scipy.stats import beta, fisher_exact

PS3_TIERS = [
    (350.0, "Very Strong"),
    (18.7, "Strong"),
    (4.33, "Moderate"),
    (2.08, "Supporting"),
]
BS3_TIERS = [(1.0 / t, label) for t, label in PS3_TIERS]  # exact reciprocals


def clopper_pearson(k, n, alpha=0.05):
    if n == 0:
        return (float("nan"), float("nan"))
    lo = 0.0 if k == 0 else beta.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - alpha / 2, k + 1, n - k)
    return (lo, hi)


Z = 1.959963984540054  # 97.5th percentile, standard normal


def odds_ratio_with_ci(a, b, c, d, alpha=0.05):
    """2x2 table [[a, c], [b, d]] -- a,c = 'depleted' counts, b,d =
    'not depleted' counts, for the numerator and denominator group
    respectively. Haldane-Anscombe (+0.5 to every cell) applied only if any
    cell is zero, matching Brnich et al. 2019's own convention for OddsPath
    at the extremes.

    This is the case-control/diagnostic odds ratio (DOR), a supplementary
    omnibus statistic -- NOT itself OddsPath. See `likelihood_ratios_with_ci`
    for the two quantities Table 3 actually calibrates against."""
    table = [[a, c], [b, d]]
    _, p = fisher_exact(table)
    if 0 in (a, b, c, d):
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_est = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lo = math.exp(math.log(or_est) - Z * se)
    hi = math.exp(math.log(or_est) + Z * se)
    return or_est, lo, hi, p


def likelihood_ratios_with_ci(a, n1, c, n2):
    """a/n1 = group_a (pathogenic-role) depletion rate (sensitivity);
    c/n2 = group_b (benign-role) depletion rate (false-positive rate,
    1 - specificity). b = n1 - a, d = n2 - c.

    Returns (lr_plus, lr_plus_lo, lr_plus_hi, lr_minus, lr_minus_lo,
    lr_minus_hi). LR+ = sensitivity / (1 - specificity) is OddsPath for an
    ABNORMAL (depleting) result -- the PS3-side quantity. LR- =
    (1 - sensitivity) / specificity is OddsPath for a NORMAL (non-depleting)
    result -- the BS3-side quantity. CIs use the standard log-scale
    likelihood-ratio formula (Simel, Samsa & Matchar 1991); Haldane-Anscombe
    (+0.5 to all four cells) applied only if any cell is zero, same
    convention as `odds_ratio_with_ci`."""
    b = n1 - a
    d = n2 - c
    if 0 in (a, b, c, d):
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        n1, n2 = a + b, c + d

    lr_plus = (a / n1) / (c / n2)
    se_plus = math.sqrt(1 / a - 1 / n1 + 1 / c - 1 / n2)
    lr_plus_lo = math.exp(math.log(lr_plus) - Z * se_plus)
    lr_plus_hi = math.exp(math.log(lr_plus) + Z * se_plus)

    lr_minus = (b / n1) / (d / n2)
    se_minus = math.sqrt(1 / b - 1 / n1 + 1 / d - 1 / n2)
    lr_minus_lo = math.exp(math.log(lr_minus) - Z * se_minus)
    lr_minus_hi = math.exp(math.log(lr_minus) + Z * se_minus)

    return lr_plus, lr_plus_lo, lr_plus_hi, lr_minus, lr_minus_lo, lr_minus_hi


def ps3_tier(lr_plus):
    for threshold, label in PS3_TIERS:
        if lr_plus >= threshold:
            return label
    return "none (below Supporting)"


def bs3_tier(lr_minus):
    for threshold, label in BS3_TIERS:
        if lr_minus <= threshold:
            return label
    return "none (below Supporting)"


def dedupe_group(df, on_duplicate):
    """One row per variant_key. `on_duplicate` resolves any variant_key with
    more than one distinct anchor_call among its matched oligos:
      - 'conservative_depleted': treat as depleted if ANY oligo says so
      - 'first': keep whichever row sorts first (fragile -- avoid for real use)
      - 'drop': exclude the variant entirely (safest once you've audited why
                 it's discordant, since a real discordant case and a bad
                 join look identical from this table alone)
    Prints every discordant variant_key it encounters, regardless of mode,
    since silently resolving one is exactly the failure mode this script
    exists to avoid repeating.
    """
    by_key = defaultdict(list)
    for _, row in df.iterrows():
        by_key[row["variant_key"]].append(row["anchor_call"])

    resolved = {}
    for k, calls in by_key.items():
        uniq = set(calls)
        if len(uniq) > 1:
            print(f"  ! discordant variant_key {k}: calls={sorted(uniq)} "
                  f"-- resolving via --on_duplicate={on_duplicate} "
                  f"(verify this is intentional, not a join artifact)",
                  file=sys.stderr)
        if on_duplicate == "drop" and len(uniq) > 1:
            continue
        if on_duplicate == "conservative_depleted" and "depleted" in uniq:
            resolved[k] = "depleted"
        else:
            resolved[k] = calls[0]
    return resolved


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--join_tsv", required=True,
                     help="assay_join_all.tsv from 06_join_assay.py")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--pathogenic_groups", nargs="+",
                     default=["CdLS_pathogenic", "DEE85_pathogenic"],
                     help="curation_group values treated as disease-causing "
                          "(default: %(default)s)")
    ap.add_argument("--benign_group", default="Exclude_Benign",
                     help="curation_group value treated as the benign/control "
                          "reference (default: %(default)s)")
    ap.add_argument("--on_duplicate", choices=["conservative_depleted", "first", "drop"],
                     default="conservative_depleted",
                     help="how to resolve a variant_key matched by more than "
                          "one oligo with discordant anchor_call (default: "
                          "%(default)s) -- see this script's docstring")
    ap.add_argument("--extra_group_tsv", action="append", default=[],
                     metavar="NAME=PATH",
                     help="repeatable. Adds a new group read from a plain "
                          "(variant_key, anchor_call) TSV instead of "
                          "--join_tsv's curation_group column -- e.g. "
                          "additional benign controls from ClinVar. Subject "
                          "to the same --on_duplicate de-duplication.")
    ap.add_argument("--pool", action="append", default=[],
                     metavar="NAME=GROUP1,GROUP2,...",
                     help="repeatable. Adds a virtual group whose "
                          "depletion count/n is the sum of the named "
                          "groups' own counts (which must already exist, "
                          "from --pathogenic_groups/--benign_group/"
                          "--extra_group_tsv). Reported alongside, not "
                          "instead of, the groups it pools.")
    ap.add_argument("--consequence_filter", nargs="+", default=None,
                     metavar="CLASS",
                     help="restrict every group (including --benign_group "
                          "and any --extra_group_tsv) to these "
                          "consequence_class values before counting -- e.g. "
                          "--consequence_filter missense inframe_deletion "
                          "inframe_insertion, to answer 'how good is this "
                          "assay specifically where annotation alone can't "
                          "already call pathogenicity', excluding PTVs "
                          "whose classification often doesn't need "
                          "functional evidence in the first place (PVS1). "
                          "Applied to the benign reference too, not just "
                          "the pathogenic groups, so the comparison stays "
                          "consequence-matched. Every --extra_group_tsv "
                          "must carry a consequence_class column when this "
                          "is set (extract_clinvar_benign_controls.py's "
                          "output does).")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.join_tsv, sep="\t", dtype=str)
    if args.consequence_filter:
        before = len(df)
        df = df[df["consequence_class"].isin(args.consequence_filter)]
        print(f"--consequence_filter {args.consequence_filter}: "
              f"{args.join_tsv} {before} -> {len(df)} row(s)", file=sys.stderr)

    groups = list(args.pathogenic_groups) + [args.benign_group]
    # Tracks which groups play the "benign/control" role for reporting
    # purposes (false-positive rate) vs "pathogenic" (sensitivity). Seeded
    # with --benign_group; an --extra_group_tsv is inferred as benign-role
    # if "benign" appears in its name (case-insensitive) -- name it
    # accordingly (e.g. "ClinVar_Benign") if that's what it is. A --pool is
    # benign-role only if every group it pools is.
    benign_role_groups = {args.benign_group}
    counts = {}
    for g in groups:
        sub = df[df["curation_group"] == g]
        if sub.empty:
            print(f"WARNING: curation_group '{g}' has no rows in {args.join_tsv}",
                  file=sys.stderr)
        resolved = dedupe_group(sub, args.on_duplicate)
        n = len(resolved)
        dep = sum(1 for c in resolved.values() if c == "depleted")
        counts[g] = (dep, n)

    for spec in args.extra_group_tsv:
        name, path = spec.split("=", 1)
        extra = pd.read_csv(path, sep="\t", dtype=str)
        missing = {"variant_key", "anchor_call"} - set(extra.columns)
        if missing:
            sys.exit(f"--extra_group_tsv {path} missing column(s): {missing}")
        if args.consequence_filter:
            if "consequence_class" not in extra.columns:
                sys.exit(f"--consequence_filter is set but {path} has no "
                          f"consequence_class column -- cannot apply the "
                          f"same filter to this group, which would make it "
                          f"an inconsistent (unfiltered) comparison")
            before = len(extra)
            extra = extra[extra["consequence_class"].isin(args.consequence_filter)]
            print(f"--consequence_filter {args.consequence_filter}: "
                  f"{path} {before} -> {len(extra)} row(s)", file=sys.stderr)
        resolved = dedupe_group(extra, args.on_duplicate)
        n = len(resolved)
        dep = sum(1 for c in resolved.values() if c == "depleted")
        counts[name] = (dep, n)
        groups.append(name)
        if "benign" in name.lower():
            benign_role_groups.add(name)
        print(f"Loaded extra group '{name}' from {path}: n={n}, depleted={dep}",
              file=sys.stderr)

    for spec in args.pool:
        name, members_str = spec.split("=", 1)
        members = members_str.split(",")
        unknown = [m for m in members if m not in counts]
        if unknown:
            sys.exit(f"--pool {name}: unknown group(s) {unknown}; must be "
                      f"one of {list(counts)} (define with --pathogenic_groups/"
                      f"--benign_group/--extra_group_tsv first)")
        dep = sum(counts[m][0] for m in members)
        n = sum(counts[m][1] for m in members)
        counts[name] = (dep, n)
        groups.append(name)
        if all(m in benign_role_groups for m in members):
            benign_role_groups.add(name)
        print(f"Pooled group '{name}' = {'+'.join(members)}: n={n}, depleted={dep}",
              file=sys.stderr)

    rows = []
    print("\n=== Per-group depletion rate ===")
    for g in groups:
        dep, n = counts[g]
        lo, hi = clopper_pearson(dep, n)
        role = "benign" if g in benign_role_groups else "pathogenic"
        label = "false-positive rate" if role == "benign" else "sensitivity"
        print(f"{g:<20} n={n:<5} depleted={dep:<5} {label}={dep/n:.1%} "
              f"(95% CI {lo:.1%}-{hi:.1%})")
        rows.append({"group": g, "role": role, "n": n, "n_depleted": dep,
                     "rate": dep / n if n else float("nan"),
                     "ci_lo": lo, "ci_hi": hi})
    pd.DataFrame(rows).to_csv(os.path.join(args.outdir, "depletion_rates.tsv"),
                               sep="\t", index=False)

    print("\n=== Pairwise OddsPath (group_a = pathogenic/case reference, "
          "group_b = benign/control reference) ===")
    or_rows = []
    for g1, g2 in combinations(groups, 2):
        a, n1 = counts[g1]
        b = n1 - a
        c, n2 = counts[g2]
        d = n2 - c
        or_est, or_lo, or_hi, p = odds_ratio_with_ci(a, b, c, d)
        (lr_plus, lr_plus_lo, lr_plus_hi,
         lr_minus, lr_minus_lo, lr_minus_hi) = likelihood_ratios_with_ci(a, n1, c, n2)
        ps3 = ps3_tier(lr_plus)
        bs3 = bs3_tier(lr_minus)
        print(f"{g1} vs {g2}: "
              f"LR+={lr_plus:.1f} (95% CI {lr_plus_lo:.1f}-{lr_plus_hi:.1f}) -> PS3 tier: {ps3}  |  "
              f"LR-={lr_minus:.2f} (95% CI {lr_minus_lo:.2f}-{lr_minus_hi:.2f}) -> BS3 tier: {bs3}  |  "
              f"DOR={or_est:.1f} (95% CI {or_lo:.1f}-{or_hi:.1f}), p={p:.2e}")
        or_rows.append({
            "group_a": g1, "group_b": g2,
            "LR_plus": lr_plus, "LR_plus_ci_lo": lr_plus_lo, "LR_plus_ci_hi": lr_plus_hi,
            "oddspath_tier_ps3": ps3,
            "LR_minus": lr_minus, "LR_minus_ci_lo": lr_minus_lo, "LR_minus_ci_hi": lr_minus_hi,
            "oddspath_tier_bs3": bs3,
            "fisher_p": p,
            "DOR": or_est, "DOR_ci_lo": or_lo, "DOR_ci_hi": or_hi,
        })
    pd.DataFrame(or_rows).to_csv(os.path.join(args.outdir, "odds_ratios.tsv"),
                                  sep="\t", index=False)

    print(f"\nWritten: {args.outdir}/depletion_rates.tsv, "
          f"{args.outdir}/odds_ratios.tsv")


if __name__ == "__main__":
    main()
