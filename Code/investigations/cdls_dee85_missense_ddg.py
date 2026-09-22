#!/usr/bin/env python3
"""
cdls_dee85_missense_ddg.py

Follow-up to cdls_vs_dee85_missense_predictors.R: that script found no
in-silico predictor (SIFT/PolyPhen/CADD/REVEL) separates CdLS_pathogenic
from DEE85_pathogenic missense variants, even though the assay's own
depletion calls do. This script asks a more specific mechanistic question --
raised as "if I trust the assay's mechanism separation, then get a delta
delta G" -- using a predicted protein-stability change (ddG) rather than a
generic pathogenicity score: do assay-depleting DEE85 missense variants look
structurally destabilizing (consistent with a loss-of-function mechanism),
while the assay-depleting CdLS missense variant does not (consistent with a
dominant-negative/altered-function mechanism that need not unfold the
protein)?

ddG is not computed by this script. It expects a pre-computed
site-saturation-mutagenesis (SSM) table from ThermoMPNN
(https://github.com/Kuhlman-Lab/ThermoMPNN, MIT license), run once against
an AlphaFold structure of human SMC1A (UniProt Q14683):

    curl -L -o thermompnn.zip \\
        https://github.com/Kuhlman-Lab/ThermoMPNN/archive/refs/heads/main.zip
    unzip thermompnn.zip
    cd ThermoMPNN-main
    uv sync --extra cpu   # installs its own Python 3.10 + CPU PyTorch, no admin needed

    curl -o structures/SMC1A_Q14683.pdb \\
        https://alphafold.ebi.ac.uk/files/AF-Q14683-F1-model_v6.pdb

    .venv/bin/python3 analysis/custom_inference.py \\
        --pdb structures/SMC1A_Q14683.pdb --chain A \\
        --out_dir <ddg_outdir>

That produces ThermoMPNN_inference_SMC1A_Q14683.csv: one row per
(position, substitution) pair for every residue in the protein (all 20 amino
acids scanned at every position), with a `position` column that is
zero-indexed against the AlphaFold/UniProt sequence. This script converts
each curated variant's HGVS protein notation (e.g. p.Glu502Lys) to that same
zero-indexed (position, wildtype, mutation) key to look up its predicted ddG,
verifying the AlphaFold sequence's residue at that position matches the
variant's own wild-type call before accepting the match -- this also
transparently handles any off-by-one/indexing differences between numbering
conventions, rather than assuming they agree.

Usage:
    python cdls_dee85_missense_ddg.py \\
        --assay_join_tsv assay_join_all.tsv \\
        --ddg_csv ThermoMPNN_inference_SMC1A_Q14683.csv \\
        --outdir ddg_analysis/
"""

import argparse
import re
import sys

import pandas as pd
from scipy import stats

AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}

HGVS_MISSENSE_RE = re.compile(r"^p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})$")


def parse_missense_hgvs(hgvs_p: str):
    """Return (position, wt_1letter, mut_1letter) for a p.XxxNNNYyy string,
    or None if it isn't a simple substitution (e.g. a deletion/indel)."""
    m = HGVS_MISSENSE_RE.match(hgvs_p.strip())
    if not m:
        return None
    wt3, pos, mut3 = m.group(1), int(m.group(2)), m.group(3)
    if wt3 not in AA3TO1 or mut3 not in AA3TO1:
        return None
    return pos, AA3TO1[wt3], AA3TO1[mut3]


def load_curated_missense(assay_join_tsv: str) -> pd.DataFrame:
    df = pd.read_csv(assay_join_tsv, sep="\t", dtype=str)
    df = df[
        df["curation_group"].isin(["CdLS_pathogenic", "DEE85_pathogenic"])
        & (df["consequence_class"] == "missense")
    ].copy()

    parsed = df["hgvs_p_mane"].map(parse_missense_hgvs)
    df = df[parsed.notna()].copy()
    df[["position", "wt", "mut"]] = pd.DataFrame(parsed[parsed.notna()].tolist(), index=df.index)
    df["position"] = df["position"].astype(int)

    df["anchor_call"] = df["anchor_tier"].map(
        lambda t: "depleted" if t in ("strongly depleting", "weakly depleting") else "no impact"
    )

    keep = ["variant_key", "hgvs_p_mane", "curation_group", "anchor_call",
            "anchor_tier", "position", "wt", "mut"]
    return df[keep].drop_duplicates("variant_key")


def merge_ddg(curated: pd.DataFrame, ddg_csv: str) -> pd.DataFrame:
    ddg = pd.read_csv(ddg_csv)
    curated = curated.copy()
    curated["position0"] = curated["position"] - 1

    merged = curated.merge(
        ddg[["position", "wildtype", "mutation", "ddG_pred"]],
        left_on=["position0", "mut"], right_on=["position", "mutation"],
        suffixes=("", "_struct"),
    )

    mismatched = merged[merged["wt"] != merged["wildtype"]]
    if len(mismatched):
        print(f"WARNING: dropping {len(mismatched)} variant(s) whose curated "
              f"wild-type residue does not match the structure at that "
              f"position (see --keep_mismatches to inspect):",
              file=sys.stderr)
        print(mismatched[["hgvs_p_mane", "wt", "wildtype"]].to_string(index=False),
              file=sys.stderr)
        merged = merged[merged["wt"] == merged["wildtype"]]

    out = merged[["variant_key", "hgvs_p_mane", "curation_group", "anchor_call",
                  "anchor_tier", "position", "wt", "mut", "ddG_pred"]]
    return out.sort_values(["curation_group", "position"]).reset_index(drop=True)


def summarise(df: pd.DataFrame):
    lines = []

    def group_stats(sub, label):
        if len(sub) == 0:
            return f"{label}: n=0"
        return (f"{label}: n={len(sub)}, median ddG={sub.median():.3f}, "
                f"mean ddG={sub.mean():.3f} kcal/mol")

    cdls = df[df.curation_group == "CdLS_pathogenic"]["ddG_pred"]
    dee85 = df[df.curation_group == "DEE85_pathogenic"]["ddG_pred"]
    lines.append("=== CdLS vs DEE85, all missense ===")
    lines.append(group_stats(cdls, "CdLS_pathogenic"))
    lines.append(group_stats(dee85, "DEE85_pathogenic"))
    if len(cdls) and len(dee85):
        u, p = stats.mannwhitneyu(cdls, dee85, alternative="two-sided")
        lines.append(f"Mann-Whitney U={u:.1f}, p={p:.4f}")

    lines.append("")
    lines.append("=== depleting vs non-depleting, within each disease group ===")
    for grp in ["CdLS_pathogenic", "DEE85_pathogenic"]:
        sub = df[df.curation_group == grp]
        dep = sub[sub.anchor_call == "depleted"]["ddG_pred"]
        nodep = sub[sub.anchor_call == "no impact"]["ddG_pred"]
        lines.append(group_stats(dep, f"{grp}, depleting"))
        lines.append(group_stats(nodep, f"{grp}, no impact"))
        if len(dep) >= 2 and len(nodep) >= 2:
            u, p = stats.mannwhitneyu(dep, nodep, alternative="two-sided")
            lines.append(f"  Mann-Whitney U={u:.1f}, p={p:.4f}")
        lines.append("")

    lines.append("=== pooled missense (both diseases): depleting vs no impact ===")
    dep = df[df.anchor_call == "depleted"]["ddG_pred"]
    nodep = df[df.anchor_call == "no impact"]["ddG_pred"]
    lines.append(group_stats(dep, "depleting"))
    lines.append(group_stats(nodep, "no impact"))
    if len(dep) and len(nodep):
        u, p = stats.mannwhitneyu(dep, nodep, alternative="two-sided")
        auc = u / (len(dep) * len(nodep))
        lines.append(f"Mann-Whitney U={u:.1f}, p={p:.4f}, AUC={auc:.3f}")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--assay_join_tsv", required=True,
                         help="assay_join_all.tsv (curated truth set joined to assay results)")
    parser.add_argument("--ddg_csv", required=True,
                         help="ThermoMPNN custom_inference.py SSM output CSV")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    import os
    os.makedirs(args.outdir, exist_ok=True)

    curated = load_curated_missense(args.assay_join_tsv)
    print(f"{len(curated)} curated CdLS/DEE85 missense variant(s) loaded", file=sys.stderr)

    merged = merge_ddg(curated, args.ddg_csv)
    out_tsv = os.path.join(args.outdir, "cdls_dee85_missense_ddg.tsv")
    merged.to_csv(out_tsv, sep="\t", index=False)
    print(f"Wrote {out_tsv} ({len(merged)} variant(s) with matched ddG)", file=sys.stderr)

    summary = summarise(merged)
    summary_path = os.path.join(args.outdir, "cdls_dee85_missense_ddg_summary.txt")
    with open(summary_path, "w") as fh:
        fh.write(summary + "\n")
    print(f"Wrote {summary_path}", file=sys.stderr)
    print("\n" + summary)


if __name__ == "__main__":
    main()
