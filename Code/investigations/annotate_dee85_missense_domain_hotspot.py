#!/usr/bin/env python3
"""
annotate_dee85_missense_domain_hotspot.py

Adds three columns to smc1a_variants_curated.tsv recording whether a
DEE85_pathogenic missense variant sits in the N-/C-terminal ATPase-head
domain that the literature (Baranano et al. 2022; Bozarth et al. 2023; Di
Nardo et al. 2026) reports as where non-loss-of-function (missense/in-frame)
DEE85 variants specifically cluster -- see
classifier_run/cdls_dee85_domain_literature_review/ for the full literature
review this is drawn from.

Only the curated project's 11 DEE85_pathogenic missense variants are
annotated (the specific claim this literature makes); every other row gets
literature_domain_hotspot=False with the domain/source columns blank, since
the claim was not evaluated for -- and does not obviously apply to -- CdLS
or non-DEE85 variants.

A variant is marked True only where a specific paper explicitly places it
in the N-terminal or C-terminal ATPase head domain (the two clustering
regions reported); variants individually discussed in the literature but
placed in the coiled-coil arm (outside the reported cluster) are marked
False, with the discussion still recorded in *_source for reference.

Usage:
    python annotate_dee85_missense_domain_hotspot.py \\
        --curated_tsv smc1a_variants_curated.tsv \\
        --output smc1a_variants_curated.tsv
"""

import argparse

import pandas as pd

# hgvs_p_mane -> (literature_domain_hotspot, domain, source)
ANNOTATIONS = {
    "p.Gly32Glu": (True, "N-terminal ATPase head",
                   "Di Nardo et al. 2026 (MD-modeled: distorts Lys38-ATP interaction)"),
    "p.Phe47Cys": (True, "N-terminal ATPase head (aa 4-148)",
                   "Baranano et al. 2022 (1 of 3 head-domain missense variants in their cohort)"),
    "p.Arg96Pro": (True, "N-terminal ATPase head (aa 4-148)",
                   "Baranano et al. 2022 (1 of 3 head-domain missense variants in their cohort)"),
    "p.Arg96Cys": (True, "N-terminal ATPase head",
                   "Di Nardo et al. 2026 (MD-modeled; same residue as Baranano's Arg96Pro)"),
    "p.Ile1081Phe": (True, "C-terminal ATPase head (adjacent)",
                     "Bozarth et al. 2023 (\"close to P2 head\")"),
    "p.Val1186Ala": (False, "", ""),
    "p.Glu313Lys": (False, "coiled-coil arm (outside reported cluster)", ""),
    "p.Val772Leu": (False, "coiled-coil arm (outside reported cluster)", ""),
    "p.Arg807His": (False, "coiled-coil arm (outside reported cluster)", ""),
    "p.Arg895Gly": (False, "second coiled-coil domain (outside reported cluster)",
                    "Kruszka et al. 2019 (individually reported, hypothesized dominant-negative"
                    " or LOF; not in the head/hinge cluster)"),
    "p.Tyr983Cys": (False, "coiled-coil (outside reported cluster)",
                    "Di Nardo et al. 2026 (MD-modeled; no substantial structural deviation found)"),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--curated_tsv", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.curated_tsv, sep="\t", dtype=str)

    df["literature_domain_hotspot"] = False
    df["literature_domain_hotspot_domain"] = ""
    df["literature_domain_hotspot_source"] = ""

    n_annotated = 0
    for hgvs_p, (is_hotspot, domain, source) in ANNOTATIONS.items():
        mask = df["hgvs_p_mane"] == hgvs_p
        n = mask.sum()
        if n == 0:
            print(f"WARNING: {hgvs_p} not found in {args.curated_tsv}")
            continue
        df.loc[mask, "literature_domain_hotspot"] = is_hotspot
        df.loc[mask, "literature_domain_hotspot_domain"] = domain
        df.loc[mask, "literature_domain_hotspot_source"] = source
        n_annotated += n

    df.to_csv(args.output, sep="\t", index=False)
    n_true = (df["literature_domain_hotspot"] == True).sum()  # noqa: E712
    print(f"Annotated {n_annotated} row(s) across {len(ANNOTATIONS)} variant(s); "
          f"{n_true} marked literature_domain_hotspot=True")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
