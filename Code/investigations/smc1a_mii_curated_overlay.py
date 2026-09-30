#!/usr/bin/env python3
"""
smc1a_mii_curated_overlay.py

Overlays the curated CdLS2 / DEE85 / benign variants onto the per-residue
Mutation Intolerance Index (`smc1a_mii.py`), to ask whether clinically
established pathogenic variants sit at residues the assay independently finds
constrained -- and whether the two disease arms differ in that respect.

**The circularity, and how it is handled.** MII is computed from the same
depletion calls the curated join uses, so a curated variant contributes to the
MII at its own position. At a median of 19 substitutions per position that is
~5.3% of the index, small but not nothing, and it biases in exactly the
direction the analysis is testing. Every comparison here therefore uses a
**leave-one-out MII**: for each curated variant, the MII at its position
recomputed with that variant's own amino acid substitution removed from both
numerator and denominator. The naive value is reported alongside so the size
of the effect is visible, but the statistics use leave-one-out.

Note this removes the variant's *own* contribution, not the contribution of
other curated variants at the same residue. Where two curated variants share a
position the residual dependence remains; `--report_shared` lists them.

**Scope.** The primary comparison is missense-only, because MII is a missense
index: CdLS 35, DEE85 11 variants. The curated benign group contributes only 3
missense variants -- the same scarcity that made the matched benign reference
unusable in the OddsPath calibration -- so a population proxy-benign set is
built from the assay's own gnomAD intersection (missense, observed in gnomAD,
not in a curated pathogenic arm), following the rule adopted in
`DECISIONS_LOG.md` D153. PTV-class curated variants are reported separately as
domain context only; a nonsense variant's MII describes how tolerant its
residue is to *missense* change, which is a different question.

Usage:
    python smc1a_mii_curated_overlay.py \\
        --mii smc1a_mii_per_position.tsv \\
        --substitution_calls smc1a_substitution_calls.tsv \\
        --join_tsv assay_join_all.tsv \\
        --gnomad_summary sge_gnomad_summary.tsv \\
        --outdir mutation_intolerance/
"""
import argparse
import os
import re

import numpy as np
import pandas as pd

ARMS = ("CdLS_pathogenic", "DEE85_pathogenic", "Exclude_Benign")
AA3 = ("Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|"
       "Trp|Tyr|Val")
SUB_RE = re.compile(rf"p\.(({AA3})\d+({AA3}))$")


def leave_one_out(subs, pos, sub, col="depleting"):
    """The index at `pos` with substitution `sub` removed.

    `col="depleting"` gives leave-one-out MII (fraction abnormal);
    `col="LFC"` gives leave-one-out mean LFC (mean substitution effect).
    Returns NaN when removal leaves nothing to average.
    """
    at = subs[subs.aa_pos == pos]
    rest = at[at.aa_sub != sub]
    if len(rest) == 0:
        return np.nan
    return rest[col].mean()


def describe(name, vals, frac_over=None):
    """One summary line. `frac_over` adds a tail fraction where that is
    meaningful (MII, bounded 0-1); omitted for mean LFC, where a 0.5 cut has
    no interpretation."""
    v = pd.Series(vals).dropna()
    if len(v) == 0:
        return f"  {name:<34} n=0"
    out = (f"  {name:<34} n={len(v):<5} median={v.median():>7.3f}  "
           f"mean={v.mean():>7.3f}")
    if frac_over is not None:
        out += f"  frac>{frac_over}={(v > frac_over).mean():.3f}"
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mii", required=True)
    ap.add_argument("--substitution_calls", required=True)
    ap.add_argument("--join_tsv", required=True)
    ap.add_argument("--gnomad_summary", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--report_shared", action="store_true",
                    help="list curated variants sharing a residue with another")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    mii = pd.read_csv(a.mii, sep="\t")
    subs = pd.read_csv(a.substitution_calls, sep="\t")
    subs["depleting"] = subs.depleting.astype(bool)
    mii_by_pos = dict(zip(mii.aa_pos, mii.MII))
    lfc_by_pos = dict(zip(mii.aa_pos, mii.mean_LFC))
    dom_by_pos = dict(zip(mii.aa_pos, mii.domain))

    j = pd.read_csv(a.join_tsv, sep="\t", dtype=str, low_memory=False)
    j = j.drop_duplicates("variant_key")
    j = j[j.curation_group.isin(ARMS)].copy()
    j["aa_pos"] = pd.to_numeric(j.protein_position, errors="coerce")
    j["aa_sub"] = j.hgvs_p_mane.str.extract(SUB_RE)[0]
    j = j[j.aa_pos.notna()].copy()
    j["aa_pos"] = j.aa_pos.astype(int)
    j["MII_naive"] = j.aa_pos.map(mii_by_pos)
    j["domain"] = j.aa_pos.map(dom_by_pos)
    j["MII_loo"] = [
        leave_one_out(subs, p, s) if isinstance(s, str) else mii_by_pos.get(p, np.nan)
        for p, s in zip(j.aa_pos, j.aa_sub)]
    j["meanLFC_naive"] = j.aa_pos.map(lfc_by_pos)
    j["meanLFC_loo"] = [
        leave_one_out(subs, p, s, "LFC") if isinstance(s, str)
        else lfc_by_pos.get(p, np.nan)
        for p, s in zip(j.aa_pos, j.aa_sub)]

    mis = j[j.consequence_class == "missense"].copy()

    # population proxy-benign missense, per D153
    g = pd.read_csv(a.gnomad_summary, sep="\t", dtype=str, low_memory=False)
    pos = g.position.astype(str).str.replace(r"\.0$", "", regex=True)
    g["variant_key"] = "chrX:" + pos + ":" + g.vcf_ref + ":" + g.vcf_alt
    g = g.drop_duplicates("variant_key")
    gm = g[g.Consequence.str.startswith("missense", na=False)
           & g.in_gnomad.astype(str).str.lower().isin(("true", "1", "yes"))].copy()
    gm["aa_pos"] = pd.to_numeric(gm.Protein_position, errors="coerce")
    gm["aa_sub"] = gm.HGVSp.str.extract(SUB_RE)[0]
    path_keys = set(j[j.curation_group != "Exclude_Benign"].variant_key)
    gm = gm[~gm.variant_key.isin(path_keys) & gm.aa_pos.notna()].copy()
    gm["aa_pos"] = gm.aa_pos.astype(int)
    gm["MII_loo"] = [
        leave_one_out(subs, p, s) if isinstance(s, str) else mii_by_pos.get(p, np.nan)
        for p, s in zip(gm.aa_pos, gm.aa_sub)]
    gm["meanLFC_loo"] = [
        leave_one_out(subs, p, s, "LFC") if isinstance(s, str)
        else lfc_by_pos.get(p, np.nan)
        for p, s in zip(gm.aa_pos, gm.aa_sub)]
    gm["domain"] = gm.aa_pos.map(dom_by_pos)

    index_cols = [("mean LFC (primary)", "meanLFC_loo", "mean_LFC", None),
                  ("MII (fraction abnormal)", "MII_loo", "MII", 0.5)]
    all_groups = {}
    for label, loo_col, bg_col, fo in index_cols:
        print(f"\n=== leave-one-out {label} at curated variant positions "
              f"(missense only)")
        groups = {}
        for arm in ARMS:
            v = mis[mis.curation_group == arm]
            groups[arm] = v[loo_col].dropna().values
            print(describe(arm, v[loo_col], fo))
        groups["gnomAD_proxy_benign"] = gm[loo_col].dropna().values
        print(describe("gnomAD proxy-benign (D153 rule)", gm[loo_col], fo))
        groups["_background"] = mii[bg_col].dropna().values
        print(describe("all assayed positions (background)",
                       groups["_background"], fo))
        all_groups[label] = groups
    groups = all_groups["MII (fraction abnormal)"]
    bg = groups["_background"]

    print("\n  naive vs leave-one-out, to show the circularity's size:")
    for arm in ARMS:
        v = mis[mis.curation_group == arm]
        if len(v):
            print(f"    {arm:<20} MII {v.MII_naive.mean():.3f} -> "
                  f"{v.MII_loo.mean():.3f} ({v.MII_loo.mean() - v.MII_naive.mean():+.3f})"
                  f"   meanLFC {v.meanLFC_naive.mean():.3f} -> "
                  f"{v.meanLFC_loo.mean():.3f} "
                  f"({v.meanLFC_loo.mean() - v.meanLFC_naive.mean():+.3f})")

    try:
        from scipy import stats as st
        for label, gset in all_groups.items():
            print(f"\n=== pairwise Mann-Whitney, leave-one-out {label}")
            keys = [k for k in ["CdLS_pathogenic", "DEE85_pathogenic",
                                "gnomAD_proxy_benign"] if len(gset.get(k, []))]
            for i in range(len(keys)):
                for k2 in keys[i + 1:]:
                    u, p = st.mannwhitneyu(gset[keys[i]], gset[k2],
                                           alternative="two-sided")
                    n1, n2 = len(gset[keys[i]]), len(gset[k2])
                    r = 1 - 2 * u / (n1 * n2)
                    print(f"  {keys[i]:<22} vs {k2:<22} "
                          f"p={p:.3g}  rank-biserial r={r:+.3f}")
            for k in keys:
                u, p = st.mannwhitneyu(gset[k], gset["_background"],
                                       alternative="two-sided")
                print(f"  {k:<22} vs background          p={p:.3g}")
    except ImportError:
        print("\n(scipy unavailable -- no tests)")

    print("\n=== domain distribution of curated variants (all consequence classes)")
    print(pd.crosstab(j.domain, j.curation_group).to_string())

    if a.report_shared:
        shared = (mis.groupby("aa_pos").size()
                  .loc[lambda s: s > 1].index.tolist())
        if shared:
            print(f"\ncurated missense variants sharing a residue: "
                  f"{len(shared)} residue(s)")
            print(mis[mis.aa_pos.isin(shared)]
                  [["aa_pos", "hgvs_p_mane", "curation_group", "MII_loo"]]
                  .sort_values("aa_pos").to_string(index=False))

    cols = ["variant_key", "hgvs_c_mane", "hgvs_p_mane", "curation_group",
            "consequence_class", "aa_pos", "aa_sub", "anchor_tier",
            "domain", "meanLFC_naive", "meanLFC_loo", "MII_naive", "MII_loo"]
    out = os.path.join(a.outdir, "smc1a_mii_curated_overlay.tsv")
    j[[c for c in cols if c in j.columns]].sort_values(
        ["curation_group", "aa_pos"]).to_csv(out, sep="\t", index=False)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
