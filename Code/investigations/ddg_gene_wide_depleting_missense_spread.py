#!/usr/bin/env python3
"""
ddg_gene_wide_depleting_missense_spread.py

Follow-up to ddg_missense_mechanism_discriminator.py, widened from the small
curated CdLS/DEE85/benign/VUS comparison groups (n in the tens) to every
missense variant the assay actually screened gene-wide, regardless of
clinical curation status -- to test whether the ddG-based mechanism split
proposed there (destabilizing -> DEE85-consistent LOF-like; non-
destabilizing -> no confident disease-mechanism-specific call) holds up as
a general pattern, and at what scale the "non-discriminating" minority
actually occurs. Framed as a candidate methodological layer for ACMG
evidence attribution generally, not only for VUS reclassification.

Input: `merged_variant_data.tsv` (../insilico_concordance/), the complete
per-oligo annotation table (every consequence class, every construct type,
all 25 targetons -- the same row count as the raw classifier output, unlike
the ClinVar- or gnomAD-restricted intersect tables). Every depleting
(weakly or strongly depleting) missense-consequence oligo is looked up
against the ThermoMPNN site-saturation scan, with the same wild-type-residue
verification used throughout this analysis chain.

Usage:
    python ddg_gene_wide_depleting_missense_spread.py \\
        --merged_variant_data_tsv merged_variant_data.tsv \\
        --ddg_csv ThermoMPNN_inference_SMC1A_Q14683.csv \\
        --outdir ddg_analysis/
"""

import argparse
import os
import re
import sys

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}
HGVS_RE = re.compile(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})")
DEPLETING_TIERS = {"weakly depleting", "strongly depleting"}
DESTABILIZING_THRESHOLD = 1.0
BENIGN_CEILING = 0.5  # matches the highest benign-missense ddG observed (0.97) rounded down


def parse_missense_hgvs(hgvs_p: str):
    m = HGVS_RE.search(hgvs_p.strip())
    if not m:
        return None
    wt3, pos, mut3 = m.group(1), int(m.group(2)), m.group(3)
    if wt3 not in AA3TO1 or mut3 not in AA3TO1:
        return None
    return pos, AA3TO1[wt3], AA3TO1[mut3]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--merged_variant_data_tsv", required=True)
    parser.add_argument("--ddg_csv", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--destabilizing_threshold", type=float, default=DESTABILIZING_THRESHOLD)
    parser.add_argument("--benign_ceiling", type=float, default=BENIGN_CEILING)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.merged_variant_data_tsv, sep="\t", dtype=str, low_memory=False)
    dep = df[(df["consequence"] == "Missense_Variant") & (df["anchor_tier"].isin(DEPLETING_TIERS))]
    dep = dep.drop_duplicates("oligo_name")
    print(f"{len(dep)} distinct depleting missense oligo(s) gene-wide", file=sys.stderr)

    ddg = pd.read_csv(args.ddg_csv)
    ddg_lookup = {(r.position, r.mutation): (r.wildtype, r.ddG_pred) for r in ddg.itertuples()}

    rows, n_no_hgvs, n_no_match, n_wt_mismatch = [], 0, 0, 0
    for r in dep.itertuples():
        hp = getattr(r, "HGVSp")
        if pd.isna(hp):
            n_no_hgvs += 1
            continue
        parsed = parse_missense_hgvs(hp.split(":")[-1])
        if parsed is None:
            n_no_match += 1
            continue
        pos, wt1, mut1 = parsed
        key = (pos - 1, mut1)
        if key not in ddg_lookup:
            n_no_match += 1
            continue
        struct_wt, ddg_pred = ddg_lookup[key]
        if struct_wt != wt1:
            n_wt_mismatch += 1
            continue
        rows.append({"oligo_name": r.oligo_name, "hgvs_p": hp.split(":")[-1],
                     "anchor_tier": getattr(r, "anchor_tier"),
                     "position": pos, "wt": wt1, "mut": mut1, "ddG_pred": ddg_pred})

    out = pd.DataFrame(rows)
    print(f"no HGVSp: {n_no_hgvs}, no ddG lookup: {n_no_match}, "
          f"wild-type mismatch: {n_wt_mismatch}, final n={len(out)}", file=sys.stderr)

    out_tsv = os.path.join(args.outdir, "gene_wide_depleting_missense_ddg.tsv")
    out.to_csv(out_tsv, sep="\t", index=False)
    print(f"Wrote {out_tsv}", file=sys.stderr)

    below_benign = (out.ddG_pred <= args.benign_ceiling).sum()
    above_destab = (out.ddG_pred > args.destabilizing_threshold).sum()
    summary_lines = [
        f"n = {len(out)}",
        f"median ddG = {out.ddG_pred.median():.3f}, mean = {out.ddG_pred.mean():.3f}",
        f"25th/75th percentile = {out.ddG_pred.quantile(0.25):.3f} / {out.ddG_pred.quantile(0.75):.3f}",
        f"fraction > {args.destabilizing_threshold} (destabilizing): "
        f"{above_destab}/{len(out)} = {above_destab/len(out):.1%}",
        f"fraction <= {args.benign_ceiling} (non-discriminating / benign-like band): "
        f"{below_benign}/{len(out)} = {below_benign/len(out):.1%}",
    ]
    summary = "\n".join(summary_lines)
    with open(os.path.join(args.outdir, "gene_wide_depleting_missense_ddg_summary.txt"), "w") as fh:
        fh.write(summary + "\n")
    print("\n" + summary)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(out.ddG_pred, bins=40, color="#4C72B0", edgecolor="white")
    ax.axvline(args.benign_ceiling, color="#999999", linestyle="--", linewidth=1,
               label=f"benign ceiling ({args.benign_ceiling})")
    ax.axvline(args.destabilizing_threshold, color="#C44E52", linestyle="--", linewidth=1,
               label=f"destabilizing threshold ({args.destabilizing_threshold})")
    ax.set_xlabel("Predicted ΔΔG (kcal/mol)")
    ax.set_ylabel("Number of depleting missense variants")
    ax.set_title(f"Gene-wide depleting missense variants (n={len(out)}):\n"
                 f"predicted structural destabilization")
    ax.legend(frameon=False)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    plt.tight_layout()
    fig_path = os.path.join(args.outdir, "gene_wide_depleting_missense_ddg_histogram.png")
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)
    print(f"Wrote {fig_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
