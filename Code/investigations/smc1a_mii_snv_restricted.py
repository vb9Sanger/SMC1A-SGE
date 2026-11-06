"""Reproduce smc1a_mii.py baseline, then repeat restricted to
single-nucleotide-reachable substitutions, under two candidate definitions."""
import os

import pandas as pd, numpy as np

# Run from the project root (the folder holding classifier_run/ and gnomAD/).
BASE = os.getcwd() + "/"
DEPLETING = ("strongly depleting", "weakly depleting")

g = pd.read_csv(BASE + "gnomAD/gnomad_intersect_results/sge_gnomad_summary.tsv",
                sep="\t", low_memory=False)
mis = g[g.Summary_Consequence == "Missense_Variant"].copy()
mis["aa_pos"] = pd.to_numeric(mis.Protein_position, errors="coerce")
mis["aa_sub"] = mis.HGVSp.str.extract(r"p\.([A-Za-z]{3}\d+[A-Za-z]{3})$")[0]
mis = mis[mis.aa_pos.notna() & mis.aa_sub.notna()].copy()
mis["aa_pos"] = mis.aa_pos.astype(int)
mis["depleting"] = mis.anchor_tier.isin(DEPLETING)
mis["score"] = pd.to_numeric(mis["pos_adj_log2FoldChange_raw"], errors="coerce")
mis = mis[mis.score.notna()]

# domain map from the published per-position output
pp = pd.read_csv(BASE + "classifier_run/mutation_intolerance/smc1a_mii_per_position.tsv",
                 sep="\t")
dom = dict(zip(pp.aa_pos, pp.domain))

def run(m, label, min_subs=5):
    grp = m.groupby(["aa_pos", "aa_sub"])
    tab = pd.concat([grp.size().rename("n_encodings"),
                     grp.depleting.sum().rename("n_dep"),
                     grp.score.mean().rename("LFC")], axis=1).reset_index()
    tab["depleting"] = tab.n_dep > 0                      # conservative policy
    per = tab.groupby("aa_pos").agg(n_subs=("aa_sub", "nunique"),
                                    n_abn=("depleting", "sum"),
                                    mean_LFC=("LFC", "mean")).reset_index()
    per["MII"] = per.n_abn / per.n_subs
    sc = per[per.n_subs >= min_subs].copy()
    sc["domain"] = sc.aa_pos.map(dom)
    print(f"\n### {label}")
    print(f"  substitutions {len(tab)}   positions scored {len(sc)}")
    print(f"  mean MII overall {sc.MII.mean():.4f}")
    for d in ["ATPase head, N lobe", "ATPase head, C lobe", "SMC hinge",
              "coiled-coil arm, N-side", "coiled-coil arm, C-side",
              "terminus, unassigned"]:
        s = sc[sc.domain == d]
        if len(s):
            print(f"    {d:<26} n={len(s):4d}  mean MII {s.MII.mean():.4f}")
    return sc

run(mis, "BASELINE (all encodings)")

snv_route = mis[mis.oligo_name.str.endswith("_snv", na=False)]
run(snv_route, "RESTRICTED (a) oligo route == _snv")

is_snv = (mis.vcf_ref.str.len() == 1) & (mis.vcf_alt.str.len() == 1)
run(mis[is_snv], "RESTRICTED (b) vcf_ref/vcf_alt both single base")
