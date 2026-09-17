#!/usr/bin/env python3
"""
Quick diagnostic: compares is_de_novo (VCF ORIGIN bit) against
any_submission_de_novo (submission_summary.txt.gz) across all annotated
TSVs from sge_clinvar_intersect_plot.py, to check whether they ever
disagree (they should, at least sometimes, if the two sources are adding
independent information) and to surface example rows where they do.

Usage:
    python check_de_novo_agreement.py /path/to/clinvar_intersect_plots/
"""
import sys
import glob
import os
import pandas as pd

outdir = sys.argv[1] if len(sys.argv) > 1 else "."
files = sorted(glob.glob(os.path.join(outdir, "*_clinvar_annotated.tsv")))

if not files:
    print(f"No *_clinvar_annotated.tsv files found in {outdir}")
    sys.exit(1)

total_in_clinvar = 0
total_origin_dn = 0
total_sub_dn = 0
total_disagree = 0
disagree_examples = []

for f in files:
    df = pd.read_csv(f, sep="\t")
    if "is_de_novo" not in df.columns or "any_submission_de_novo" not in df.columns:
        print(f"  MISSING COLUMNS in {f} — skipping")
        continue
    cv = df[df["in_clinvar"]]
    total_in_clinvar += len(cv)
    total_origin_dn += int(cv["is_de_novo"].sum())
    total_sub_dn += int(cv["any_submission_de_novo"].sum())

    disagree = cv[cv["is_de_novo"] != cv["any_submission_de_novo"]]
    if not disagree.empty:
        total_disagree += len(disagree)
        cols = [c for c in ["clinvar_id", "oligo_name", "is_de_novo",
                             "any_submission_de_novo", "origin_flags",
                             "n_submissions", "n_submissions_de_novo",
                             "de_novo_submitters"] if c in disagree.columns]
        disagree_examples.append((f, disagree[cols]))

print(f"Files checked: {len(files)}")
print(f"Total ClinVar-matched rows: {total_in_clinvar}")
print(f"Total is_de_novo (ORIGIN):        {total_origin_dn}")
print(f"Total any_submission_de_novo:      {total_sub_dn}")
print(f"Total disagreeing rows:            {total_disagree}")
print()

if disagree_examples:
    print("Example disagreements:")
    for f, ex in disagree_examples[:8]:
        print(f"--- {f} ---")
        print(ex.to_string(index=False))
        print()
else:
    print("No disagreements found — the two columns are identical everywhere "
          "they overlap in these files.")

# Cross-check against the gene-wide pre-scan number the main script printed
# ("N / M VariationID(s) have at least one submission reporting de novo").
# If total_sub_dn (summed across targetons, with overlap double-counting) is
# noticeably LESS than that gene-wide N, some flagged variants are being
# lost between the pre-scan and the per-targeton merge.
print()
print("If total_sub_dn is smaller than the gene-wide pre-scan count printed "
      "at the start of the main script's run (the '43 / 991' line), that's "
      "evidence some flagged VariationIDs never make it into any targeton's "
      "merged output. Check whether their clinvar_id shows up at all in any "
      "*_clinvar_annotated.tsv (in_clinvar column) vs. only in the VCF.")
