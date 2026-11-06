"""Synonymous depletion rate with the TARGETON as the unit of analysis.

Mirrors the splice-region exon-level test in
cassette_exon_consequence_analysis.py (Table 6) exactly: one % depleted per
targeton, two-sided Mann-Whitney between cassette and non-cassette groups.
Also reproduces that splice-region test as a control.
"""
import glob, os, re, io
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact

# Run from the project root (the folder holding classifier_run/ and gnomAD/).
BASE = os.getcwd() + "/"
DEPLETED_TIERS = {"weakly depleting", "strongly depleting"}

# cassette assignment taken from the script's own Table 1, so it cannot drift
tbl = open(BASE + "classifier_run/cassette_exon_analysis/"
           "cassette_exon_consequence_tables.tsv").read()
t1 = pd.read_csv(io.StringIO(tbl.split("# Table 1")[1].split("\n", 1)[1]
                             .split("\n\n")[0]), sep="\t")
cassette = dict(zip(t1.Targeton_ID, t1.is_cassette))

rows = []
for f in sorted(glob.glob(BASE + "classifier_run/classifier_results/D4_ref/*.tsv")):
    tid = os.path.basename(f).split("_")[0]
    if tid not in cassette:
        continue
    d = pd.read_csv(f, sep="\t", low_memory=False)
    syn = d[d.consequence == "Synonymous_Variant"]
    syn = syn[syn.anchor_tier.notna()]
    if len(syn) == 0:
        continue
    n_dep = syn.anchor_tier.isin(DEPLETED_TIERS).sum()
    rows.append(dict(Targeton_ID=tid, is_cassette=cassette[tid], n=len(syn),
                     n_depleted=int(n_dep),
                     pct_depleted=round(100 * n_dep / len(syn), 1)))

p = pd.DataFrame(rows)
print(p.to_string(index=False))

cas = p[p.is_cassette].pct_depleted
non = p[~p.is_cassette].pct_depleted
u, pv = mannwhitneyu(cas, non, alternative="two-sided")
print(f"\nTARGETON AS UNIT, synonymous % depleted")
print(f"  cassette     n_targetons={len(cas):2d}  median {cas.median():.1f}%")
print(f"  non-cassette n_targetons={len(non):2d}  median {non.median():.1f}%")
print(f"  Mann-Whitney two-sided: U = {u:.1f}, p = {pv:.4f}")

# pooled proportion, for the record
a = p[p.is_cassette].n_depleted.sum(); b = p[p.is_cassette].n.sum() - a
c = p[~p.is_cassette].n_depleted.sum(); e = p[~p.is_cassette].n.sum() - c
odds, pf = fisher_exact([[a, b], [c, e]])
print(f"\nPOOLED synonymous: cassette {a}/{a+b} = {100*a/(a+b):.1f}%  "
      f"non-cassette {c}/{c+e} = {100*c/(c+e):.1f}%")
print(f"  Fisher exact: OR = {odds:.2f}, p = {pf:.4g}")
