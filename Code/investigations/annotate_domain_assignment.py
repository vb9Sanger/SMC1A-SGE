#!/usr/bin/env python3
"""
annotate_domain_assignment.py

Assigns every curated variant to a protein domain from the structural
annotation, and flags whether it falls in the region the DEE85 literature
reports as the cluster for non-truncating variants: the N-/C-terminal
ATPase head and the central hinge (Bozarth et al. 2023; Baranano et al.
2022; Di Nardo et al. 2026 -- see
classifier_run/cdls_dee85_domain_literature_review/).

**This replaces an earlier publication-based approach**, which marked a
variant as being in the cluster only where a named publication explicitly
placed it there. That approach had a coverage bias it could not escape: a
variant could only be marked True if some paper had individually discussed
it, so a DEE85 missense variant sitting squarely inside the ATPase head but
never written about was marked False. Testing whether the literature's
claim holds in this dataset, using a method that can only recover variants
the literature already named, is close to circular.

Domain boundaries come from the same annotation used for the per-residue
constraint analysis (`smc1a_domain_annotation.json`, built from UniProt and
InterPro), so the two analyses speak the same structural language:

    ATPase head   InterPro IPR028468 "Smc1, ATP-binding cassette domain",
                  residues 4-148 and 1117-1220
    SMC hinge     UniProt `Domain` "SMC hinge", residues 515-629
                  (Pfam PF06470 and SMART SM00968 agree to within 2 aa)

The literature annotation is retained, but demoted to **provenance**: which
variants were individually reported, by whom, and where that publication
placed them. Disagreements between the two are reported rather than
reconciled, as they are informative in both directions -- a publication may
place a variant "close to" a domain it does not formally fall within, and a
variant no publication discussed may still fall inside one.

Usage:
    python annotate_domain_assignment.py \\
        --curated_tsv smc1a_variants_curated.tsv \\
        --domain_json smc1a_domain_annotation.json \\
        --output smc1a_variants_curated.tsv
"""

import argparse
import json

import pandas as pd

# The consequence classes the literature claim concerns: it contrasts
# non-truncating (missense, in-frame) DEE85 variants, which it reports as
# clustered, against truncating ones, which it reports as gene-wide.
NON_TRUNCATING = {"missense", "inframe_deletion", "inframe_insertion"}

HEAD_ACCESSION = "IPR028468"
HINGE_DESCRIPTION = "SMC hinge"


def load_domains(path):
    """ATPase head and hinge intervals, from the shared domain annotation."""
    d = json.load(open(path))
    head = [(f["start"], f["end"]) for f in d["pfam"]
            if f["accession"] == HEAD_ACCESSION]
    hinge = [(f["start"], f["end"]) for f in d["features"]
             if f["type"] == "Domain" and f["description"] == HINGE_DESCRIPTION]
    if not head or not hinge:
        raise SystemExit(f"could not find head/hinge features in {path}")
    return head, hinge


def in_any(pos, intervals):
    return any(a <= pos <= b for a, b in intervals)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--curated_tsv", required=True)
    ap.add_argument("--domain_json", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    head, hinge = load_domains(a.domain_json)
    print(f"ATPase head : {head}")
    print(f"SMC hinge   : {hinge}\n")

    df = pd.read_csv(a.curated_tsv, sep="\t", dtype=str)
    pos = pd.to_numeric(df.protein_position, errors="coerce")

    df["domain_assignment"] = ""
    df["in_atpase_head"] = ""
    df["in_hinge"] = ""
    df["in_reported_cluster"] = ""

    has_pos = pos.notna()
    for i in df.index[has_pos]:
        p = int(pos[i])
        h, g = in_any(p, head), in_any(p, hinge)
        df.at[i, "in_atpase_head"] = str(h)
        df.at[i, "in_hinge"] = str(g)
        df.at[i, "in_reported_cluster"] = str(h or g)
        df.at[i, "domain_assignment"] = (
            "ATPase head" if h else "SMC hinge" if g else "outside head and hinge")

    df.to_csv(a.output, sep="\t", index=False)
    print(f"Wrote {a.output} ({int(has_pos.sum())} of {len(df)} rows carry a "
          f"protein position)\n")

    # ---- summary: the claim concerns non-truncating variants, per arm
    nt = df[has_pos & df.consequence_class.isin(NON_TRUNCATING)]
    print("Non-truncating curated variants, by arm and domain")
    print(pd.crosstab(nt.curation_group, nt.domain_assignment).to_string(), "\n")

    for arm in sorted(nt.curation_group.dropna().unique()):
        s = nt[nt.curation_group == arm]
        k = (s.in_reported_cluster == "True").sum()
        print(f"  {arm:<20} {k}/{len(s)} in the reported cluster "
              f"({100 * k / len(s):.1f}%)")

    # ---- where the structural and literature assignments disagree
    #
    # Compared only within the scope the literature annotation attempted:
    # DEE85 missense variants. It never addressed CdLS or truncating
    # variants, so a "disagreement" there would only restate that fact.
    if {"literature_domain_hotspot",
            "literature_domain_hotspot_source"} <= set(df.columns):
        scope = (has_pos & df.curation_group.eq("DEE85_pathogenic")
                 & df.consequence_class.eq("missense"))
        lit = df.literature_domain_hotspot.astype(str).str.lower().eq("true")
        struct = df.in_reported_cluster.eq("True")
        clash = df[scope & (lit != struct)]
        print(f"\nStructural vs literature assignment, within the "
              f"{int(scope.sum())} DEE85 missense variants the literature "
              f"annotation covered: {len(clash)} disagreement(s)")
        for _, r in clash.iterrows():
            print(f"  {r.hgvs_p_mane:<16} residue {int(float(r.protein_position)):<5} "
                  f"literature={str(r.literature_domain_hotspot):<6} "
                  f"structural={r.domain_assignment}")
            placed = r.literature_domain_hotspot_domain
            if isinstance(placed, str) and placed.strip():
                print(f"{'':>18}literature placed it: {placed}")


if __name__ == "__main__":
    main()
