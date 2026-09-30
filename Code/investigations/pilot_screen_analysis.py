#!/usr/bin/env python3
"""
pilot_screen_analysis.py

Tabulates the SMC1A pilot SGE screen (exons 3 and 4, two independent sgRNA HDR
libraries each), which was run ahead of the full 25-targeton screen to confirm
that the assay discriminates variant classes in SMC1A/HAP1.

Reports, per HDR library: depletion rate and median LFC by consequence class,
the median LFC of protein-truncating against synonymous variants at each
timepoint, and a flag for any library returning no significance calls.

Consequence classes use `Summary_Consequence` (`DECISIONS_LOG.md` D156). "PTV"
means `Nonsense_Variant` + `Frameshift_Variant`.

Note: pilot LFCs are not adjusted for positional bias, which the full screen's
`pos_adj_log2FoldChange_raw` does account for, so the two are not on the same
footing and are never combined. This script takes no full-screen input, and
nothing it produces feeds another analysis.

Usage:
    python pilot_screen_analysis.py \\
        --pilot_xlsx "SMC1A_Pilot/SMC1A Pilot LFC.xlsx" \\
        --outdir classifier_run/pilot_screen/
"""
import argparse
import os

import pandas as pd

# sheet in the pilot workbook -> library label
PILOT_LIBRARIES = {
    "SMC1A_exon3_sg1": "exon3_sg1",
    "SMC1A_exon3_sg2": "exon3_sg2",
    "SMC1A_exon4_sg1": "exon4_sg1",
    "SMC1A_exon4_sg2": "exon4_sg2",
}

PTV_CLASSES = ("Nonsense_Variant", "Frameshift_Variant")
TIMEPOINTS = ("D7", "D11", "D15", "D21")

# most to least disruptive, so the tables read as a gradient
CSQ_ORDER = ["Nonsense_Variant", "Frameshift_Variant", "Inframe_Deletion",
             "Splice_Variant", "Missense_Variant",
             "Splice_Polypyrimidine_Tract_Variant",
             "Synonymous_Variant", "Intronic_Variant"]


def load_pilot(xlsx):
    """One row per (library, oligo)."""
    xl = pd.ExcelFile(xlsx)
    missing = [s for s in PILOT_LIBRARIES if s not in xl.sheet_names]
    if missing:
        raise SystemExit(f"pilot workbook is missing sheet(s): {missing}")
    frames = []
    for sheet, label in PILOT_LIBRARIES.items():
        d = xl.parse(sheet)
        d["library"] = label
        pos = d.vcf_pos.astype(str).str.replace(r"\.0$", "", regex=True)
        d["variant_key"] = (d.chrom.astype(str) + ":" + pos + ":"
                            + d.vcf_ref.astype(str) + ":" + d.vcf_alt.astype(str))
        frames.append(d)
    p = pd.concat(frames, ignore_index=True)
    p["is_depleted"] = p.Day15_significance.eq("Depleted")
    return p


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pilot_xlsx", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    p = load_pilot(a.pilot_xlsx)

    say("=== libraries")
    for lab, g in p.groupby("library"):
        say(f"  {lab:<11} n_oligos={len(g):<6} sgRNA={g.sgRNA.iloc[0]} "
            f"targeton={g.Targeton.iloc[0]:<7} depleted={int(g.is_depleted.sum())}")
    say(f"  distinct variants: {p.variant_key.nunique()} "
        f"(from {len(p)} oligo-level rows)")
    say()

    say("=== depletion by consequence class, within each library (Day 15)")
    rows = [{"library": lab, "Summary_Consequence": csq, "n": len(d),
             "n_depleted": int(d.is_depleted.sum()),
             "pct_depleted": round(100 * d.is_depleted.mean(), 1),
             "median_LFC": round(d.adj_lfc_D15.median(), 3)}
            for lab, g in p.groupby("library")
            for csq, d in g.groupby("Summary_Consequence")]
    dep = pd.DataFrame(rows)
    dep["_o"] = dep.Summary_Consequence.apply(
        lambda c: CSQ_ORDER.index(c) if c in CSQ_ORDER else len(CSQ_ORDER))
    dep = dep.sort_values(["library", "_o"]).drop(columns="_o")
    for lab, g in dep.groupby("library"):
        say(f"  -- {lab}")
        say(g.drop(columns="library").to_string(index=False))
        say()
    dep.to_csv(os.path.join(a.outdir, "pilot_depletion_by_consequence.tsv"),
               sep="\t", index=False)

    say("=== median LFC by timepoint: PTV vs synonymous")
    traj = []
    for lab, g in p.groupby("library"):
        ptv = g[g.Summary_Consequence.isin(PTV_CLASSES)]
        syn = g[g.Summary_Consequence.eq("Synonymous_Variant")]
        r = {"library": lab}
        for tp in TIMEPOINTS:
            r[f"PTV_{tp}"] = round(ptv[f"adj_lfc_{tp}"].median(), 2)
            r[f"syn_{tp}"] = round(syn[f"adj_lfc_{tp}"].median(), 2)
        traj.append(r)
    traj = pd.DataFrame(traj)
    say(traj.to_string(index=False))
    traj.to_csv(os.path.join(a.outdir, "pilot_timepoint_trajectory.tsv"),
                sep="\t", index=False)
    say()

    say("=== per-library QC")
    for lab, g in p.groupby("library"):
        ptv = g[g.Summary_Consequence.isin(PTV_CLASSES)]
        n_dep = int(g.is_depleted.sum())
        flag = ("  <-- no significant calls anywhere, despite a PTV median LFC "
                f"of {ptv.adj_lfc_D15.median():+.2f}" if n_dep == 0 else "")
        say(f"  {lab:<11} depleted calls={n_dep:<5} "
            f"PTV median LFC={ptv.adj_lfc_D15.median():+.3f}  "
            f"D15->D21 PTV shift="
            f"{ptv.adj_lfc_D21.median() - ptv.adj_lfc_D15.median():+.2f}{flag}")
    say()

    cols = ["variant_key", "library", "sgRNA", "HGVSc", "HGVSp",
            "Consequence", "Summary_Consequence", "Summary_Plot",
            "adj_lfc_D7", "adj_lfc_D11", "adj_lfc_D15", "adj_lfc_D21",
            "BH_FDR_D15", "Day15_significance", "is_depleted"]
    path = os.path.join(a.outdir, "pilot_variant_table.tsv")
    p[cols].to_csv(path, sep="\t", index=False)
    say(f"Wrote {path} ({len(p)} rows)")

    with open(os.path.join(a.outdir, "pilot_screen_summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
