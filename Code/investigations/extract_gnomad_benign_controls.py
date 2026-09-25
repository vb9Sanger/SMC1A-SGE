#!/usr/bin/env python3
"""
extract_gnomad_benign_controls.py

Builds a population-frequency benign control set for the OddsPath calibration,
from the assay's own gnomAD intersection (`sge_gnomad_summary.tsv`).

Rationale: SMC1A's two disorders are severe, early-onset and essentially
always de novo, so a pathogenic allele has close to zero reproductive fitness
and should not persist in a population database at appreciable frequency. A
variant seen repeatedly in gnomAD is therefore evidence against pathogenicity
that is independent of any clinical assertion -- which is the point, since the
curated benign group is small (n=25) and the ClinVar benign group inherits
whatever ClinVar's submitters believed.

**Which frequency rule.** For SMC1A/CdLS the Whiffin et al. 2017 maximum
credible population AF is 5.0e-7, and that sits BELOW the resolution of the
database: at the observed AN of ~1.1M a single allele is already 9.2e-7. So
every variant observed in gnomAD at all exceeds the maximum credible
pathogenic frequency, and the interesting question is not the AF cut but
whether to additionally require the CONFIDENCE BOUND to clear it. Two modes:

  --proxy_clinical         every observed variant is a control. CanVIG-UK
                           (Allen/Rowlands et al. 2026) Recommendation 4,
                           which extends the maximum-tolerated-allele-frequency
                           principle this far for very rare, highly penetrant,
                           early-onset disease. The default for this work.

  --max_credible_af 5e-7   Whiffin's rigorous form: filter on the FILTERING
                           allele frequency (the CI lower bound), not the point
                           estimate. The required allele count is derived, not
                           asserted -- here AC >= 3. Robustness comparator.

  --min_af / --min_ac      the original hand-set rule (grpmax AF > 1e-4 and
                           AC >= 2). Retired: the AF cut sat 200x above maxAF
                           for no stated reason and AC >= 2 is in fact more
                           permissive than Whiffin's own faf95 requirement.

Two further filters apply in every mode, and each was added after looking at
what the data does without it:

1. **Deduplicate by genomic variant.** `sge_gnomad_summary.tsv` is per oligo,
   and SMC1A's targeton windows overlap, so 1,313 gnomAD-observed rows are
   only 1,272 distinct variants. Left as rows, the overlapping windows would
   count the same variant twice in the control group.

2. **Exclude anything asserted pathogenic elsewhere.** A variant in the
   curated CdLS/DEE85 arms, or called P/LP by ClinVar, is removed regardless
   of frequency: a control group must not contain the cases. Worth noting that
   under the faf95 rule this removes nothing -- no variant asserted pathogenic
   anywhere reaches AC >= 3 in gnomAD, so that criterion excludes every known
   pathogenic variant without being told to.

Output is a (variant_key, anchor_call, anchor_tier) TSV in the shape
`10_calibrate_by_tier.py --gnomad_benign` and
`09_calibrate_sensitivity_oddspath.py --extra_group_tsv` expect.

Usage:
    python extract_gnomad_benign_controls.py \\
        --gnomad_summary sge_gnomad_summary.tsv \\
        --curated_tsv smc1a_variants_all.tsv \\
        --clinvar_summary clinvar_variants_summary.tsv \\
        --proxy_clinical \\
        --output gnomad_benign_controls.tsv
"""
import argparse
import re
import sys

import pandas as pd

ARM_GROUPS = ("CdLS_pathogenic", "DEE85_pathogenic")
PLP = ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic")
KEY_RE = re.compile(r"(chrX):(\d+)_([ACGT]+)>([ACGT]+)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gnomad_summary", required=True)
    ap.add_argument("--curated_tsv", required=True,
                    help="smc1a_variants_all.tsv, to exclude the pathogenic arms")
    ap.add_argument("--clinvar_summary", default=None,
                    help="clinvar_variants_summary.tsv, to exclude ClinVar P/LP")
    ap.add_argument("--min_af", type=float, default=1e-4,
                    help="minimum gnomAD grpmax allele frequency [%(default)s]")
    ap.add_argument("--min_ac", type=int, default=2,
                    help="minimum total allele count; guards against a single "
                         "allele in a small ancestry group producing a high "
                         "grpmax AF [%(default)s]")
    ap.add_argument("--af_field", default="combined_grpmax_af",
                    choices=["combined_grpmax_af", "pooled_af"])
    ap.add_argument("--proxy_clinical", action="store_true",
                    help="CanVIG-UK Recommendation 4 mode: take EVERY variant "
                         "observed in gnomAD as a proxy-clinical benign "
                         "control, ignoring --min_af and --min_ac. Allen, "
                         "Rowlands et al. 2026 extend the maximum tolerated "
                         "allele frequency principle so that, for a very rare "
                         "highly-penetrant early-onset disorder, mere presence "
                         "in a population database is the evidence -- the "
                         "allele count is not the point, because an allele of "
                         "near-zero reproductive fitness should not be there "
                         "at all. Much larger reference, weaker per-variant "
                         "claim; report alongside the stricter set, not "
                         "instead of it.")
    ap.add_argument("--max_credible_af", type=float, default=None,
                    help="Whiffin et al. 2017 maximum credible population "
                         "allele frequency for this gene and disease. When "
                         "given, --min_ac is DERIVED rather than asserted: the "
                         "script finds the smallest allele count whose 95%% "
                         "lower confidence bound on the allele frequency (the "
                         "filtering allele frequency, gnomAD faf95) exceeds "
                         "this value, and requires at least that many alleles. "
                         "This is what the paper actually recommends filtering "
                         "on, and it is the reason a singleton carries no "
                         "benign evidence however high its point-estimate "
                         "grpmax AF looks. Overrides --min_ac and --min_af.")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    g = pd.read_csv(a.gnomad_summary, sep="\t", dtype=str, low_memory=False)
    g["_in"] = g["in_gnomad"].astype(str).str.lower().isin(("true", "1", "yes"))
    g["_af"] = pd.to_numeric(g[a.af_field], errors="coerce")
    for c in ("gnomad_exomes_ac", "gnomad_genomes_ac"):
        g["_" + c] = pd.to_numeric(g.get(c), errors="coerce").fillna(0)
    g["_ac"] = g["_gnomad_exomes_ac"] + g["_gnomad_genomes_ac"]
    for c in ("gnomad_exomes_an", "gnomad_genomes_an"):
        g["_" + c] = pd.to_numeric(g.get(c), errors="coerce").fillna(0)
    g["_an"] = g["_gnomad_exomes_an"] + g["_gnomad_genomes_an"]

    m = g.oligo_name.str.extract(KEY_RE)
    g["variant_key"] = m[0] + ":" + m[1] + ":" + m[2] + ":" + m[3]

    obs = g[g["_in"] & g.variant_key.notna()].copy()
    print(f"gnomAD-observed assayed rows: {len(obs)}  "
          f"({obs.variant_key.nunique()} distinct variants)")

    # 1. one row per genomic variant. Overlapping targeton windows otherwise
    #    contribute the same variant twice. Where a variant appears in two
    #    windows its anchor_tier is taken from the first; disagreement is
    #    reported rather than silently resolved.
    disagree = (obs.groupby("variant_key").anchor_tier.nunique() > 1)
    if disagree.any():
        print(f"  NOTE: {int(disagree.sum())} variant(s) carry different "
              f"anchor_tier in two overlapping windows; keeping the first and "
              f"listing them:")
        for k in disagree[disagree].index[:10]:
            print(f"    {k}: {sorted(set(obs[obs.variant_key == k].anchor_tier))}")
    obs = obs.drop_duplicates("variant_key")

    # 2. frequency and allele-count filters
    if a.max_credible_af:
        from scipy.stats import chi2
        an = float(obs["_an"].median())
        # One-sided 95% lower bound on a Poisson rate given k observed alleles.
        # gnomAD computes faf95 per ancestry group on a smaller AN, so the
        # count derived here from the total AN is a floor, not the exact
        # gnomAD-equivalent threshold.
        k = 1
        while k < 1000 and (0.5 * chi2.ppf(0.05, 2 * k)) / an <= a.max_credible_af:
            k += 1
        print(f"WHIFFIN FAF MODE: maxAF={a.max_credible_af:g}, median AN="
              f"{an:,.0f}\n  smallest AC whose faf95 exceeds maxAF: {k} "
              f"(faf95={0.5 * chi2.ppf(0.05, 2 * k) / an:.2e}); "
              f"AC={k - 1} gives {0.5 * chi2.ppf(0.05, 2 * (k - 1)) / an:.2e}, "
              "below maxAF")
        keep = obs[obs["_ac"] >= k].copy()
        print(f"  kept {len(keep)} variants with AC >= {k}")
    elif a.proxy_clinical:
        keep = obs.copy()
        print("PROXY-CLINICAL MODE (CanVIG-UK Rec 4): every gnomAD-observed "
              f"variant kept, AF and AC filters skipped: {len(keep)}")
    else:
        keep = obs[(obs["_af"] > a.min_af) & (obs["_ac"] >= a.min_ac)].copy()
        print(f"after grpmax AF > {a.min_af:g} and AC >= {a.min_ac}: {len(keep)}")

    # 3. remove anything asserted pathogenic elsewhere
    cur = pd.read_csv(a.curated_tsv, sep="\t", dtype=str, low_memory=False)
    path_keys = set(cur[cur.curation_group.isin(ARM_GROUPS)].variant_key)
    n0 = len(keep)
    keep = keep[~keep.variant_key.isin(path_keys)]
    print(f"  removed {n0 - len(keep)} in a curated pathogenic arm")

    if a.clinvar_summary:
        cv = pd.read_csv(a.clinvar_summary, sep="\t", dtype=str, low_memory=False)
        if "clnsig_norm" in cv.columns:
            cvp = cv[cv.clnsig_norm.isin(PLP)]
            cvm = cvp.oligo_name.str.extract(KEY_RE) if "oligo_name" in cvp.columns else None
            if cvm is not None:
                cvk = set((cvm[0] + ":" + cvm[1] + ":" + cvm[2] + ":" + cvm[3]).dropna())
                n0 = len(keep)
                keep = keep[~keep.variant_key.isin(cvk)]
                print(f"  removed {n0 - len(keep)} called P/LP by ClinVar")

    keep = keep[keep.anchor_tier.notna()]
    out = keep[["variant_key", "anchor_tier"]].copy()
    # the calibration script consumes a binary anchor_call; anchor_tier is kept
    # alongside so the same file also serves the by-tier ("bucket") analysis
    out["anchor_call"] = out.anchor_tier.map(
        lambda t: "no impact" if t == "no impact" else
                  ("enriched" if t == "enriched" else "depleted"))
    out = out[["variant_key", "anchor_call", "anchor_tier"]]
    out.to_csv(a.output, sep="\t", index=False)

    print(f"\nWrote {a.output} ({len(out)} variants)")
    print(out.anchor_tier.value_counts().to_string())


if __name__ == "__main__":
    main()
