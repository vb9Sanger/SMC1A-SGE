#!/usr/bin/env python3
"""
smc1a_mii.py

Per-residue Mutation Intolerance Index (MII) for SMC1A from the SGE data,
stratified by protein domain.

MII is the fraction of amino acid substitutions at a protein position that the
assay calls functionally abnormal: 0.00 = every substitution tolerated,
1.00 = every substitution abnormal. It converts a per-variant functional map
into a per-residue measure of mutational constraint, which can then be
stratified by domain to ask whether the assay recapitulates known
structure-function relationships.

The index follows the definition in Jaramillo-Calle, D., PhD thesis
("Mutation intolerance indices", §5.3.4.6), applied there to CTCF. Two
deliberate departures, both forced by the SMC1A library being more complete:

**1. Substitutions, not variants.** The original counts missense *variants*
per position. This library includes multi-nucleotide codon replacements, so
30% of amino acid substitutions are reachable by more than one genomic variant
(median 1.32 encodings each, up to 25 variants per position encoding 19
substitutions). Counting variants would weight a substitution by how many
codon paths reach it, which is a property of the genetic code and the library
design, not of the residue. This collapses to distinct `p.Xxx###Yyy` first, so
the denominator is the number of amino acid changes actually tested -- 19 of
the 19 possible at 94% of positions, which makes MII directly comparable
between positions.

**2. No minimum-observations filter by default.** The original requires >=5
distinct missense variants per position and drops the rest as unstable. Here
1,166 of 1,170 covered positions clear that threshold anyway, so it excludes
nothing; `--min_subs` is kept for symmetry and set to 5.

**Redundant encodings that disagree.** Where several genomic variants encode
one substitution and they do not agree on depleting/not, `--on_discordant`
decides: `conservative` (default; any depleting call makes the substitution
depleting), `majority`, or `drop`. The concordance rate across the 6,565
multiply-encoded substitutions is reported either way and is worth quoting on
its own -- it measures assay reproducibility across independent nucleotide
changes with an identical protein consequence, which is the internal control
the original also uses.

**Domain annotation** is fetched from UniProt and InterPro rather than
hard-coded, cached to `--domain_cache`, and written into the output so every
position carries both its coarse label and the source that label came from.
The ATPase head is taken from InterPro IPR028468, which annotates it as two
sequence fragments because the domain is discontinuous -- see
`build_domain_map` for why UniProt has no equivalent feature.

Usage:
    python smc1a_mii.py \\
        --gnomad_summary sge_gnomad_summary.tsv \\
        --outdir mutation_intolerance/
"""
import argparse
import json
import os
import sys
import urllib.request

import numpy as np
import pandas as pd

UNIPROT = "Q14683"
# Full record rather than a field subset: the first version of this script
# requested only ft_domain/ft_coiled/ft_region/ft_binding and so never saw the
# secondary-structure records, which is part of why the ATPase head looked
# unannotated.
UNIPROT_URL = f"https://rest.uniprot.org/uniprotkb/{UNIPROT}.json"
PFAM_URL = ("https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/"
            f"{UNIPROT}/?page_size=100")
# InterPro's integrated entry for the SMC1 ATPase head, sourced from CDD
# cd03275. It is annotated as TWO fragments because the domain is
# discontinuous in sequence -- see build_domain_map.
ATPASE_ENTRY = "IPR028468"
DEPLETING = ("strongly depleting", "weakly depleting")


def fetch_annotation(cache):
    """UniProt features plus Pfam domains, cached as JSON."""
    if cache and os.path.exists(cache):
        with open(cache) as fh:
            return json.load(fh)
    out = {"uniprot": UNIPROT, "features": [], "pfam": []}
    with urllib.request.urlopen(UNIPROT_URL, timeout=60) as r:
        u = json.load(r)
    out["length"] = u["sequence"]["length"]
    for f in u.get("features", []):
        loc = f["location"]
        out["features"].append({"type": f["type"],
                                "start": loc["start"]["value"],
                                "end": loc["end"]["value"],
                                "description": f.get("description", "")})
    try:
        with urllib.request.urlopen(PFAM_URL, timeout=60) as r:
            p = json.load(r)
        for res in p.get("results", []):
            md = res["metadata"]
            for prot in res.get("proteins", []):
                for l in prot.get("entry_protein_locations", []):
                    for frag in l["fragments"]:
                        out["pfam"].append({"accession": md["accession"],
                                            "name": md.get("name") or "",
                                            "db": md.get("source_database", ""),
                                            "start": frag["start"],
                                            "end": frag["end"]})
    except Exception as e:                      # Pfam is supplementary
        print(f"WARNING: Pfam lookup failed ({e}); continuing on UniProt only",
              file=sys.stderr)
    if cache:
        os.makedirs(os.path.dirname(cache) or ".", exist_ok=True)
        with open(cache, "w") as fh:
            json.dump(out, fh, indent=2)
    return out


def build_domain_map(ann):
    """One coarse domain label per residue, plus the fine annotations.

    **The ATPase head is annotated, not derived.** UniProt has no `Domain`
    feature for it, which is a limitation of the data model rather than missing
    knowledge: a UniProt `Domain` carries one contiguous start-end, and the SMC
    ATPase head is *discontinuous in sequence*. The nucleotide-binding domain is
    formed when the two long coiled-coil arms fold back and bring the
    N-terminal segment (carrying the Walker A / P-loop, UniProt `Binding site`
    32-39, ligand ATP) into apposition with the C-terminal segment (carrying
    Walker B and the ABC signature motif). Expressing that needs two fragments.
    Pfam cannot either -- its SMC_N model (PF02463) matches 3-514 and 630-1208,
    i.e. arm plus head together, so it does not isolate the head.

    The databases that *do* support discontinuous matches all annotate it:

        InterPro IPR028468  "Smc1, ATP-binding cassette domain"   4-148; 1117-1220
        CDD      cd03275    ATP-binding cassette domain of SMC1   4-148; 1117-1220
        SUPERFAM SSF52540   P-loop NTP hydrolase                  7-171; 1116-1203
        CATH     3.40.50.300 P-loop NTP hydrolase                 3-221; 966-1225

    IPR028468 is used here: it is the InterPro integrated entry, it is specific
    to SMC1 rather than to P-loop ATPases generally, and its two fragments are
    the tightest of the four. The others are recorded in the annotation string
    so a reader can see the spread -- the head boundary is somewhere in
    148-221 on the N side and 966-1117 on the C side depending on whether you
    take the ABC cassette or the whole P-loop superfamily fold.
    """
    n = ann["length"]
    hinge = [(f["start"], f["end"]) for f in ann["features"]
             if f["type"] == "Domain" and "hinge" in f["description"].lower()]
    if not hinge:
        sys.exit("no SMC hinge domain in the UniProt annotation -- aborting "
                 "rather than guessing its position")
    h_start, h_end = hinge[0]

    head = sorted((p["start"], p["end"]) for p in ann["pfam"]
                  if p["accession"] == ATPASE_ENTRY)
    if len(head) != 2:
        sys.exit(f"expected 2 fragments for {ATPASE_ENTRY} (the ATPase head is "
                 f"discontinuous); got {len(head)}. Refusing to fall back on a "
                 f"derived boundary silently -- check the InterPro response.")
    (n_lo, n_hi), (c_lo, c_hi) = head
    provenance = {}

    coarse = {}
    for pos in range(1, n + 1):
        if n_lo <= pos <= n_hi:
            coarse[pos] = "ATPase head, N lobe"
            provenance[pos] = f"InterPro {ATPASE_ENTRY} (CDD cd03275)"
        elif c_lo <= pos <= c_hi:
            coarse[pos] = "ATPase head, C lobe"
            provenance[pos] = f"InterPro {ATPASE_ENTRY} (CDD cd03275)"
        elif h_start <= pos <= h_end:
            coarse[pos] = "SMC hinge"
            provenance[pos] = "UniProt Domain / Pfam PF06470"
        elif pos < n_lo or pos > c_hi:
            coarse[pos] = "terminus, unassigned"
            provenance[pos] = "outside every annotated domain"
        elif pos < h_start:
            coarse[pos] = "coiled-coil arm, N-side"
            provenance[pos] = "between the N lobe and the hinge"
        else:
            coarse[pos] = "coiled-coil arm, C-side"
            provenance[pos] = "between the hinge and the C lobe"

    fine = {pos: [] for pos in range(1, n + 1)}
    for f in ann["features"]:
        if f["type"] in ("Chain", "Natural variant", "Sequence conflict"):
            continue
        for pos in range(f["start"], f["end"] + 1):
            label = f["type"] + (f": {f['description']}" if f["description"] else "")
            fine[pos].append(label)
    for p in ann["pfam"]:
        for pos in range(p["start"], min(p["end"], n) + 1):
            fine[pos].append(f"{p['db']}:{p['accession']} {p['name']}".strip())
    return coarse, {k: "; ".join(sorted(set(v))) for k, v in fine.items()}, provenance


def collapse_substitutions(mis, policy):
    """One row per (position, amino acid substitution).

    Returns the collapsed frame and the concordance statistics across
    substitutions reachable by more than one genomic variant.
    """
    grp = mis.groupby(["aa_pos", "aa_sub"])
    n_enc = grp.size().rename("n_encodings")
    n_dep = grp.depleting.sum().rename("n_encodings_depleting")
    # Mean across redundant encodings of one substitution, so a substitution
    # reachable by three codon paths still contributes once -- the same
    # collapse MII uses, applied to the continuous score.
    score = grp.score.mean().rename("LFC")
    tab = pd.concat([n_enc, n_dep, score], axis=1).reset_index()

    multi = tab[tab.n_encodings > 1]
    disc = multi[(multi.n_encodings_depleting > 0)
                 & (multi.n_encodings_depleting < multi.n_encodings)]
    stats = {"n_substitutions": len(tab),
             "n_multi_encoded": len(multi),
             "n_discordant": len(disc),
             "concordance": (1 - len(disc) / len(multi)) if len(multi) else float("nan")}

    if policy == "conservative":
        tab["depleting"] = tab.n_encodings_depleting > 0
    elif policy == "majority":
        tab["depleting"] = tab.n_encodings_depleting * 2 > tab.n_encodings
    elif policy == "drop":
        keep = ~((tab.n_encodings_depleting > 0)
                 & (tab.n_encodings_depleting < tab.n_encodings))
        tab = tab[keep].copy()
        tab["depleting"] = tab.n_encodings_depleting > 0
    else:
        sys.exit(f"unknown --on_discordant policy: {policy}")
    return tab, stats


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gnomad_summary", required=True,
                    help="sge_gnomad_summary.tsv -- every assayed variant with "
                         "its Protein_position, HGVSp and anchor_tier")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--domain_cache", default=None,
                    help="JSON cache for the UniProt/Pfam annotation "
                         "[<outdir>/smc1a_domain_annotation.json]")
    ap.add_argument("--score_field", default="pos_adj_log2FoldChange_raw",
                    help="continuous assay score, for the mean-LFC index "
                         "[%(default)s]")
    ap.add_argument("--min_subs", type=int, default=5,
                    help="minimum distinct substitutions for a position to be "
                         "scored [%(default)s]")
    ap.add_argument("--on_discordant", default="conservative",
                    choices=["conservative", "majority", "drop"])
    ap.add_argument("--abnormal", default="depleting",
                    choices=["depleting", "strongly_only"],
                    help="what counts as functionally abnormal: either "
                         "depletion tier, or strongly depleting alone "
                         "[%(default)s]")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    cache = a.domain_cache or os.path.join(a.outdir, "smc1a_domain_annotation.json")

    ann = fetch_annotation(cache)
    coarse, fine, provenance = build_domain_map(ann)
    print(f"annotation: UniProt {ann['uniprot']}, {ann['length']} aa, "
          f"{len(ann['features'])} features, {len(ann['pfam'])} Pfam domains")

    g = pd.read_csv(a.gnomad_summary, sep="\t", dtype=str, low_memory=False)
    pos = g.position.astype(str).str.replace(r"\.0$", "", regex=True)
    g["variant_key"] = "chrX:" + pos + ":" + g.vcf_ref + ":" + g.vcf_alt
    g = g.drop_duplicates("variant_key")

    # Summary_Consequence, not the raw VEP Consequence and not Summary_Plot.
    # `Consequence.startswith("missense")` also catches the 1,038
    # `missense_variant,splice_region_variant` calls that both summary columns
    # class as Splice_Variant, and 28.5% of those deplete -- consistent with a
    # real splice effect rather than a missense one. Summary_Consequence rather
    # than Summary_Plot because the two differ only in that Summary_Plot
    # collapses nonsense and frameshift into "LOF", which is a mechanism label
    # rather than a consequence class; Summary_Consequence maps 1:1 onto
    # smc1a_lib.consequence_class.
    if "Summary_Consequence" not in g.columns:
        sys.exit("--gnomad_summary has no Summary_Consequence column; that is the "
                 "field this project classes variants by, so refusing to fall "
                 "back on the raw VEP Consequence")
    mis = g[g.Summary_Consequence == "Missense_Variant"].copy()
    mis["aa_pos"] = pd.to_numeric(mis.Protein_position, errors="coerce")
    mis["aa_sub"] = mis.HGVSp.str.extract(r"p\.([A-Za-z]{3}\d+[A-Za-z]{3})$")[0]
    mis = mis[mis.aa_pos.notna() & mis.aa_sub.notna()].copy()
    mis["aa_pos"] = mis.aa_pos.astype(int)
    abnormal = DEPLETING if a.abnormal == "depleting" else ("strongly depleting",)
    mis["depleting"] = mis.anchor_tier.isin(abnormal)
    mis["score"] = pd.to_numeric(mis[a.score_field], errors="coerce")
    n_before = len(mis)
    mis = mis[mis.score.notna()]
    if len(mis) < n_before:
        print(f"  dropped {n_before - len(mis)} variant(s) with no "
              f"{a.score_field}", file=sys.stderr)
    print(f"{len(mis)} distinct missense variants at {mis.aa_pos.nunique()} positions")

    subs, cstats = collapse_substitutions(mis, a.on_discordant)
    print(f"collapsed to {cstats['n_substitutions']} distinct amino acid "
          f"substitutions\n  {cstats['n_multi_encoded']} reachable by >1 "
          f"genomic variant; {cstats['n_discordant']} of those disagree "
          f"=> concordance {cstats['concordance']:.1%}")
    print(f"  discordant handling: --on_discordant {a.on_discordant}")

    per_pos = subs.groupby("aa_pos").agg(
        n_subs=("aa_sub", "nunique"),
        n_abnormal=("depleting", "sum"),
        mean_LFC=("LFC", "mean"),
        median_LFC=("LFC", "median")).reset_index()
    per_pos["MII"] = per_pos.n_abnormal / per_pos.n_subs
    scored = per_pos[per_pos.n_subs >= a.min_subs].copy()
    scored["domain"] = scored.aa_pos.map(coarse)
    scored["domain_source"] = scored.aa_pos.map(provenance)
    scored["annotation"] = scored.aa_pos.map(fine)
    dropped = len(per_pos) - len(scored)
    print(f"\n{len(scored)} positions scored ({dropped} dropped for "
          f"< {a.min_subs} substitutions)")
    print(f"  MII      (fraction abnormal): median "
          f"{scored.MII.median():.3f}, mean {scored.MII.mean():.3f}")
    print(f"  mean LFC (mean substitution effect): median "
          f"{scored.mean_LFC.median():.3f}, mean {scored.mean_LFC.mean():.3f} "
          f"(more negative = less tolerant)")
    try:
        from scipy import stats as _st
        rho, prho = _st.spearmanr(scored.MII, scored.mean_LFC)
        print(f"  the two indices correlate at Spearman rho = {rho:.3f} "
              f"(p = {prho:.2g})")
    except ImportError:
        pass
    print(f"  fully dispensable (MII = 0): {(scored.MII == 0).sum()} "
          f"({(scored.MII == 0).mean():.1%})")
    print(f"  fully intolerant  (MII = 1): {(scored.MII == 1).sum()} "
          f"({(scored.MII == 1).mean():.1%})")

    pos_path = os.path.join(a.outdir, "smc1a_mii_per_position.tsv")
    scored.to_csv(pos_path, sep="\t", index=False)
    sub_path = os.path.join(a.outdir, "smc1a_substitution_calls.tsv")
    subs.to_csv(sub_path, sep="\t", index=False)

    print("\n=== MII by domain ===")
    rows = []
    for dom, sub in scored.groupby("domain"):
        rows.append(dict(domain=dom, n_positions=len(sub),
                         median_mean_LFC=round(sub.mean_LFC.median(), 3),
                         mean_mean_LFC=round(sub.mean_LFC.mean(), 3),
                         median_MII=round(sub.MII.median(), 3),
                         mean_MII=round(sub.MII.mean(), 3),
                         frac_dispensable=round((sub.MII == 0).mean(), 3),
                         frac_MII_over_0_5=round((sub.MII > 0.5).mean(), 3)))
    dom_df = pd.DataFrame(rows).sort_values("median_mean_LFC")
    print(dom_df.to_string(index=False))
    dom_path = os.path.join(a.outdir, "smc1a_mii_by_domain.tsv")
    dom_df.to_csv(dom_path, sep="\t", index=False)

    try:
        from scipy import stats as st
        names = [d for d, _ in scored.groupby("domain")]
        for index_name, col in (("mean LFC", "mean_LFC"), ("MII", "MII")):
            gs = [s[col].values for _, s in scored.groupby("domain")]
            if len(gs) > 2:
                h, p = st.kruskal(*gs)
                print(f"\nKruskal-Wallis across {len(gs)} domains, "
                      f"{index_name}: H={h:.1f}, p={p:.3g}")
        groups = [s.MII.values for _, s in scored.groupby("domain")]
        print("\npairwise Mann-Whitney on MII (two-sided), BH-adjusted:")
        pairs, raw = [], []
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                u, p = st.mannwhitneyu(groups[i], groups[j], alternative="two-sided")
                pairs.append((names[i], names[j])); raw.append(p)
        order = np.argsort(raw)
        adj = np.empty(len(raw)); m = len(raw); prev = 1.0
        for rank, idx in enumerate(order[::-1]):
            prev = min(prev, raw[idx] * m / (m - rank))
            adj[idx] = prev
        for (x, y), p, q in sorted(zip(pairs, raw, adj), key=lambda t: t[2]):
            print(f"  {x:<34} vs {y:<34} p={p:.2g}  q={q:.2g}")
    except ImportError:
        print("\n(scipy not available -- skipping significance tests)")

    with open(os.path.join(a.outdir, "smc1a_mii_summary.txt"), "w") as fh:
        fh.write(f"UniProt {ann['uniprot']} ({ann['length']} aa)\n")
        fh.write(f"abnormal = {a.abnormal}; on_discordant = {a.on_discordant}; "
                 f"min_subs = {a.min_subs}\n")
        fh.write(f"redundant-encoding concordance = {cstats['concordance']:.4f} "
                 f"({cstats['n_discordant']}/{cstats['n_multi_encoded']} discordant)\n\n")
        fh.write(dom_df.to_string(index=False) + "\n")

    print(f"\nWrote {pos_path}\n      {sub_path}\n      {dom_path}")


if __name__ == "__main__":
    main()
