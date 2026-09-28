#!/usr/bin/env python3
"""
build_excalibr_input.py

Builds ExCALIBR `BasicScoreset` input from the SMC1A SGE results, one file per
disease arm.

ExCALIBR (Zeiberg, Stewart et al. 2026, bioRxiv 2025.04.29.651326;
https://github.com/rosstewart/exCALIBR) calibrates an assay by jointly
modelling the SCORE DISTRIBUTIONS of four samples as skew-normal mixtures, and
returns a per-variant posterior probability of pathogenicity expressed as ACMG
evidence points from -8 to +8. That is a different object from the Brnich
OddsPath calibration in `09_/10_calibrate_*.py`, which discards the score and
tiers a categorical call: there, every "strongly depleting" variant receives
the same evidence regardless of how far past the threshold it sits.

The BasicScoreset format is two columns:

    score,sample_assignments

where `sample_assignments` indexes 0=Pathogenic, 1=Benign, 2=gnomAD,
3=Synonymous, comma-separated where a variant belongs to more than one, and
blank for variants that inform nothing but still receive an evidence
assignment (ExCALIBR treats those as VUS). The gnomAD sample is mandatory in
every mode -- it is what the prior probability of pathogenicity is estimated
from, via label-shift correction, because clinical controls are ascertainment
biased and cannot supply a population prior.

**Why one file per arm.** SMC1A causes two disorders by two mechanisms -- CdLS2
(MIM 300590, missense/dominant-negative) and DEE85 (MIM 301044,
PTV/loss-of-function) -- and the assay detects one and not the other. Pooling
them into a single pathogenic sample would model a bimodal pathogenic
distribution as though it were one population, and the skew-normal mixture
would fit the average of two mechanisms that behave differently. Running the
arms separately is also the extension the paper explicitly flags as
unevaluated: "the calibration approach has not been evaluated on genes with
both gain- and loss-of-function disease mechanisms". The benign, gnomAD and
synonymous samples are shared between the two runs; only sample 0 changes.

**Three judgement calls, each recorded because none is forced by the format.**

1. **Splice-region synonymous variants are excluded from the synonymous
   sample** (--keep_splice_synonymous to retain them). ExCALIBR uses the
   synonymous sample as a functionally-normal reference, "under the assumption
   that synonymous variants are predominantly functionally normal". In this
   assay that assumption fails for exactly this subset: synonymous variants
   that also touch a splice region deplete at 10/88 (11.4%) against 23/2545
   (0.9%) for synonymous variants elsewhere, a 12-fold enrichment. SGE edits
   the endogenous locus and therefore reports splicing, so those depletions are
   real rather than noise -- which is precisely why they must not sit in the
   normal reference. The paper's own pipeline drops splice-consequence rows
   when an assay cannot measure splicing; the same reasoning inverts for an
   assay that can.

2. **No SpliceAI filtering.** The paper removes variants with SpliceAI > 0.2
   for assays that do not detect splice effects, since a pathogenic effect
   there cannot be attributed to protein function. 64 of their 80 datasets were
   cDNA-based and needed this; SGE is one of the 19 that does not.

3. **Curated-pathogenic variants that also appear in gnomAD keep both labels.**
   A pathogenic allele can be present in a population database at very low
   frequency, and the prior estimate depends on the gnomAD sample being a
   population sample rather than a filtered one. ExCALIBR resolves the overlap
   itself (assigning to one sample at random, synonymous taking precedence), so
   the honest thing is to declare both memberships and let it apply its rule.

Usage:
    python build_excalibr_input.py \\
        --gnomad_summary sge_gnomad_summary.tsv \\
        --join_tsv assay_join_all.tsv \\
        --clinvar_summary clinvar_variants_summary.tsv \\
        --outdir excalibr_calibration/
"""
import argparse
import os

import pandas as pd

ARMS = ("CdLS_pathogenic", "DEE85_pathogenic")
BLB = ("Benign", "Likely benign", "Benign/Likely benign")
# ExCALIBR BasicScoreset sample indices
PATHOGENIC, BENIGN, GNOMAD, SYNONYMOUS = 0, 1, 2, 3


def variant_key(df):
    """`chrX:pos:ref:alt` from the VCF columns.

    Deliberately NOT parsed out of `oligo_name`, the way the other scripts in
    this repo do it. That pattern (`chrX:pos_REF>ALT`) only matches SNVs --
    deletions are named `chrX:53380079_1del` and `chrX:53380081_53380082_2del1`
    -- so keying on it silently drops 25,279 of 39,135 assayed rows and takes
    the CdLS arm from 41 variants to 36. The VCF triple covers every variant
    type.
    """
    pos = df["position"].astype(str).str.replace(r"\.0$", "", regex=True)
    return "chrX:" + pos + ":" + df["vcf_ref"] + ":" + df["vcf_alt"]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gnomad_summary", required=True,
                    help="sge_gnomad_summary.tsv -- carries every assayed "
                         "variant, its score, its consequence and whether it "
                         "is in gnomAD")
    ap.add_argument("--join_tsv", required=True,
                    help="assay_join_all.tsv from the all-variants scope: "
                         "supplies curation_group and, via oligo_name, the "
                         "link to the assay table")
    ap.add_argument("--clinvar_summary", required=True)
    ap.add_argument("--score_field", default="pos_adj_log2FoldChange_raw",
                    help="the continuous assay score [%(default)s]")
    ap.add_argument("--keep_splice_synonymous", action="store_true",
                    help="retain splice-region synonymous variants in the "
                         "synonymous sample; see note 1 in the docstring")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    g = pd.read_csv(a.gnomad_summary, sep="\t", dtype=str, low_memory=False)
    g["variant_key"] = variant_key(g)
    g["score"] = pd.to_numeric(g[a.score_field], errors="coerce")
    g = g[g.variant_key.notna() & g.score.notna()]
    oligo_to_key = dict(zip(g.oligo_name, g.variant_key))
    # One row per genomic variant: the table is per oligo and the targeton
    # windows overlap, so a variant in two windows would enter a sample twice.
    n_rows = len(g)
    g = g.drop_duplicates("variant_key").set_index("variant_key")
    print(f"assayed variants: {n_rows} rows -> {len(g)} distinct")

    # Group membership comes from the assay join rather than from matching
    # curated keys against the assay table directly. 06_join_assay.py already
    # resolved HGVS 3'-shifting against VCF left-alignment for this
    # minus-strand gene; re-deriving it here by string equality loses 2 CdLS
    # and 3 DEE85 variants. A curated variant matched by several oligos is
    # collapsed to one, preferring the oligo whose own key equals the curated
    # key.
    j = pd.read_csv(a.join_tsv, sep="\t", dtype=str, low_memory=False)
    j["assay_key"] = j.oligo_name.map(oligo_to_key)
    j = j[j.assay_key.notna()]
    j["_exact"] = (j.assay_key == j.variant_key).astype(int)
    j = (j.sort_values("_exact", ascending=False)
           .drop_duplicates("variant_key", keep="first"))

    def group_keys(name):
        return set(j[j.curation_group == name].assay_key) & set(g.index)

    arms = {arm: group_keys(arm) for arm in ARMS}
    benign = group_keys("Exclude_Benign")

    cv = pd.read_csv(a.clinvar_summary, sep="\t", dtype=str, low_memory=False)
    cv["variant_key"] = variant_key(cv)
    benign |= set(cv[cv.clnsig_norm.isin(BLB)].variant_key.dropna())
    benign &= set(g.index)

    in_gnomad = set(g.index[g.in_gnomad.astype(str).str.lower()
                            .isin(("true", "1", "yes"))])

    syn_mask = g.Consequence.str.contains("synonymous", na=False)
    if not a.keep_splice_synonymous:
        splice = g.Consequence.str.contains("splice", na=False)
        dropped = int((syn_mask & splice).sum())
        syn_mask &= ~splice
        print(f"synonymous sample: excluded {dropped} splice-region "
              f"synonymous variants (see docstring note 1)")
    synonymous = set(g.index[syn_mask])

    # A variant asserted pathogenic must never sit in the benign or synonymous
    # reference, whatever its consequence says.
    all_path = arms[ARMS[0]] | arms[ARMS[1]]
    for name, s in (("benign", benign), ("synonymous", synonymous)):
        clash = s & all_path
        if clash:
            print(f"  removed {len(clash)} curated-pathogenic variant(s) from "
                  f"the {name} sample")
            s -= clash
    benign -= all_path
    synonymous -= all_path

    print(f"\nshared samples: benign={len(benign)}  gnomAD={len(in_gnomad)}  "
          f"synonymous={len(synonymous)}")
    print(f"overlaps: benign&gnomAD={len(benign & in_gnomad)}  "
          f"synonymous&gnomAD={len(synonymous & in_gnomad)}  "
          f"benign&synonymous={len(benign & synonymous)}")

    for arm in ARMS:
        path = arms[arm]
        rows = []
        for k, score in g.score.items():
            m = []
            if k in path:
                m.append(PATHOGENIC)
            if k in benign:
                m.append(BENIGN)
            if k in in_gnomad:
                m.append(GNOMAD)
            if k in synonymous:
                m.append(SYNONYMOUS)
            rows.append((score, ",".join(str(i) for i in m)))
        out = pd.DataFrame(rows, columns=["score", "sample_assignments"])
        dest = os.path.join(a.outdir, f"smc1a_{arm}.csv")
        out.to_csv(dest, index=False)
        labelled = (out.sample_assignments != "").sum()
        print(f"\n{arm}: pathogenic={len(path)}  "
              f"{len(out)} rows ({labelled} labelled, "
              f"{len(out) - labelled} unlabelled/VUS) -> {dest}")


if __name__ == "__main__":
    main()
