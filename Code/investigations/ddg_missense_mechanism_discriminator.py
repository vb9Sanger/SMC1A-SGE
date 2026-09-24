#!/usr/bin/env python3
"""
ddg_missense_mechanism_discriminator.py

Follow-up to cdls_dee85_missense_ddg.py, prompted by a specific objection to
the VUS reclassification (09_calibrate_sensitivity_oddspath.py /
reclassify_clinvar_vus.py, STEP SEVEN): mapping a depleting missense VUS to
the CdLS-calibrated OddsPath tier is questionable on mechanistic grounds,
since depletion is fundamentally a loss-of-function-consistent readout --
closer to DEE85's mechanism than to CdLS's dominant-negative one -- and only
one confirmed depleting CdLS missense variant exists to calibrate that
mapping against.

This script asks whether predicted protein-stability change (ddG, from the
same ThermoMPNN site-saturation scan used by cdls_dee85_missense_ddg.py) can
discriminate, for a depleting missense variant, whether it is behaving in a
DEE85-consistent (destabilizing) or CdLS-consistent/non-discriminating
(non-destabilizing) way -- and specifically whether "non-destabilizing" is
actually CdLS-specific, or just what an ordinary missense variant (benign or
pathogenic) looks like structurally.

Adds a benign-missense reference group (ClinVar B/LB missense variants, from
clinvar_variants_summary.tsv) to the CdLS/DEE85 comparison already in
cdls_dee85_missense_ddg.tsv, and scores the depleting ClinVar VUS missense
set (from vus_reclassification.tsv) on the same scale.

Usage:
    python ddg_missense_mechanism_discriminator.py \\
        --cdls_dee85_ddg_tsv cdls_dee85_missense_ddg.tsv \\
        --clinvar_summary_tsv clinvar_variants_summary.tsv \\
        --vus_reclassification_tsv vus_reclassification.tsv \\
        --ddg_csv ThermoMPNN_inference_SMC1A_Q14683.csv \\
        --outdir ddg_analysis/
"""

import argparse
import os
import re
import sys

import pandas as pd

AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}
HGVS_MISSENSE_RE = re.compile(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})")
BENIGN_LABELS = {"Benign", "Likely benign", "Benign/Likely benign"}
DESTABILIZING_THRESHOLD = 1.0  # kcal/mol, matches cdls_dee85_missense_ddg.py's usage


def parse_missense_hgvs(hgvs_p: str):
    m = HGVS_MISSENSE_RE.search(hgvs_p.strip())
    if not m:
        return None
    wt3, pos, mut3 = m.group(1), int(m.group(2)), m.group(3)
    if wt3 not in AA3TO1 or mut3 not in AA3TO1:
        return None
    return pos, AA3TO1[wt3], AA3TO1[mut3]


def lookup_ddg(ddg: pd.DataFrame, pos: int, wt1: str, mut1: str):
    match = ddg[(ddg.position == pos - 1) & (ddg.mutation == mut1)]
    if not len(match):
        return None, None
    row = match.iloc[0]
    return row.wildtype, row.ddG_pred


def load_benign_missense(clinvar_summary_tsv: str, ddg: pd.DataFrame) -> pd.DataFrame:
    df = pd.read_csv(clinvar_summary_tsv, sep="\t", dtype=str)
    df = df[(df["Summary_Plot"] == "Missense_Variant") & (df["clnsig_norm"].isin(BENIGN_LABELS))]
    rows = []
    for _, r in df.drop_duplicates("HGVSp").iterrows():
        parsed = parse_missense_hgvs(r["HGVSp"])
        if parsed is None:
            continue
        pos, wt1, mut1 = parsed
        struct_wt, ddg_pred = lookup_ddg(ddg, pos, wt1, mut1)
        if struct_wt is None or struct_wt != wt1:
            continue
        rows.append({"hgvs_p": r["HGVSp"].split(":")[-1], "group": "Benign_missense",
                     "anchor_call": r["anchor_call"], "position": pos,
                     "wt": wt1, "mut": mut1, "ddG_pred": ddg_pred})
    return pd.DataFrame(rows)


def load_depleting_vus_missense(vus_tsv: str, ddg: pd.DataFrame) -> pd.DataFrame:
    df = pd.read_csv(vus_tsv, sep="\t", dtype=str)
    df = df[df["recommendation"] == "PS3 Strong"]
    rows = []
    for _, r in df.iterrows():
        parsed = parse_missense_hgvs(r["HGVSp"])
        if parsed is None:
            continue
        pos, wt1, mut1 = parsed
        struct_wt, ddg_pred = lookup_ddg(ddg, pos, wt1, mut1)
        if struct_wt is None or struct_wt != wt1:
            continue
        rows.append({"hgvs_p": r["HGVSp"].split(":")[-1], "group": "Depleting_missense_VUS",
                     "anchor_tier": r["anchor_tier"], "position": pos,
                     "wt": wt1, "mut": mut1, "ddG_pred": ddg_pred,
                     "mechanism_call": ("DEE85-consistent" if ddg_pred > DESTABILIZING_THRESHOLD
                                        else "no disease-specific tier")})
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cdls_dee85_ddg_tsv", required=True,
                         help="Output of cdls_dee85_missense_ddg.py")
    parser.add_argument("--clinvar_summary_tsv", required=True,
                         help="clinvar_variants_summary.tsv (from sge_clinvar_intersect.py)")
    parser.add_argument("--vus_reclassification_tsv", required=True,
                         help="Output of reclassify_clinvar_vus.py")
    parser.add_argument("--ddg_csv", required=True,
                         help="ThermoMPNN custom_inference.py SSM output CSV")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--destabilizing_threshold", type=float, default=DESTABILIZING_THRESHOLD)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ddg = pd.read_csv(args.ddg_csv)
    cdls_dee85 = pd.read_csv(args.cdls_dee85_ddg_tsv, sep="\t")

    benign = load_benign_missense(args.clinvar_summary_tsv, ddg)
    print(f"{len(benign)} ClinVar benign/likely-benign missense variant(s) matched", file=sys.stderr)

    vus = load_depleting_vus_missense(args.vus_reclassification_tsv, ddg)
    print(f"{len(vus)} depleting missense VUS (PS3 Strong candidates) matched", file=sys.stderr)

    benign.to_csv(os.path.join(args.outdir, "benign_missense_ddg.tsv"), sep="\t", index=False)
    vus.to_csv(os.path.join(args.outdir, "vus_missense_ddg_mechanism_call.tsv"), sep="\t", index=False)

    def summarize(label, series):
        if len(series) == 0:
            return f"{label}: n=0"
        return (f"{label}: n={len(series)}, median={series.median():.3f}, "
                f"mean={series.mean():.3f}, frac>{args.destabilizing_threshold}="
                f"{(series > args.destabilizing_threshold).mean():.2f}")

    lines = ["=== ddG by group (benign, CdLS, DEE85, depleting VUS) ==="]
    lines.append(summarize("Benign_missense", benign.ddG_pred))
    for grp in ["CdLS_pathogenic", "DEE85_pathogenic"]:
        for call in ["no impact", "depleted"]:
            sub = cdls_dee85[(cdls_dee85.curation_group == grp) & (cdls_dee85.anchor_call == call)]
            lines.append(summarize(f"{grp}, {call}", sub.ddG_pred))
    lines.append(summarize("Depleting_missense_VUS", vus.ddG_pred))
    lines.append("")
    lines.append(f"VUS mechanism-call breakdown (threshold={args.destabilizing_threshold} kcal/mol):")
    lines.append(vus.mechanism_call.value_counts().to_string())

    summary = "\n".join(lines)
    with open(os.path.join(args.outdir, "ddg_mechanism_discriminator_summary.txt"), "w") as fh:
        fh.write(summary + "\n")
    print("\n" + summary)


if __name__ == "__main__":
    main()
