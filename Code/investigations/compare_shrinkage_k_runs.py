#!/usr/bin/env python3
"""
compare_shrinkage_k_runs.py

Sensitivity check for gaussian_shrinkage_classifier.R's shrinkage strength
(--shrinkage_k_lof / --shrinkage_k_noimpact). The auto-selected k is a
heuristic (median per-targeton n), not fit from the data's actual
between-targeton variability -- this script quantifies how much that
choice actually matters by diffing anchor_tier calls across 2+ full runs
of the classifier made with different k values.

Recommended runs to generate on the farm (k_lof_auto=68.5, k_ctl_auto=175.5
for the current SMC1A dataset -- confirm your own auto values from your
run's console output or gmm_shrinkage_anchor_summary.tsv):

    # half k
    Rscript Code/gaussian_shrinkage_classifier.R --input "..." \
        --out_dir results_k_half --shrinkage_k_lof 34.25 --shrinkage_k_noimpact 87.75

    # current / auto (baseline -- already have this)
    Rscript Code/gaussian_shrinkage_classifier.R --input "..." \
        --out_dir results_k_auto

    # double k
    Rscript Code/gaussian_shrinkage_classifier.R --input "..." \
        --out_dir results_k_double --shrinkage_k_lof 137 --shrinkage_k_noimpact 351

Then:
    python compare_shrinkage_k_runs.py \
        --run results_k_half:half \
        --run results_k_auto:auto \
        --run results_k_double:double \
        --outdir sensitivity_results/

Each --run points at one output directory from a classifier run (containing
the per-targeton *.tsv files with an anchor_tier column) plus a short label.
The FIRST --run given is treated as the baseline all others are compared to.

Matches rows across runs by (Targeton_ID-equivalent inferred from filename)
+ position + oligo_name if present, falling back to row order within a
matched targeton file if no stable ID column is found -- row order is safe
here because every run classifies the exact same input rows, just with
different anchor parameters.

Outputs:
  - anchor_tier_flip_summary.tsv: per targeton, how many variants changed
    tier between baseline and each comparison run
  - anchor_tier_flip_detail.tsv: every variant whose tier changed in any
    comparison run, with old/new tier, consequence, and identifying columns
    (position, oligo_name by default -- pass --id_cols to change) so flips
    can be joined against downstream files like clinvar_variants_summary.tsv
  - anchor_tier_flip_by_consequence.tsv: flip counts broken down by
    consequence class (are flips concentrated in missense? LOF? etc.)
"""

import argparse
import glob
import os
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--run", action="append", required=True,
        help="outdir:label, e.g. results_k_auto:auto. First one given is the baseline.",
    )
    p.add_argument("--outdir", default=".", help="Where to write comparison outputs")
    p.add_argument("--consequence_col", default="consequence", help="Consequence column name (default: consequence)")
    p.add_argument("--tier_col", default="anchor_tier", help="Tier column name (default: anchor_tier)")
    p.add_argument(
        "--id_cols", default="position,oligo_name",
        help="Comma-separated columns to carry into the flip-detail output for downstream joins "
             "(e.g. against a ClinVar summary file). Default: position,oligo_name. "
             "Silently skips any that aren't present in the baseline file.",
    )
    return p.parse_args()


def load_run(outdir):
    """Load every per-targeton TSV in a classifier output dir, keyed by targeton (filename stem)."""
    files = sorted(glob.glob(os.path.join(outdir, "*.tsv")))
    files = [f for f in files if not os.path.basename(f).startswith("gmm_shrinkage_anchor_summary")]
    if not files:
        sys.exit(f"ERROR: no per-targeton *.tsv files found in {outdir}")
    frames = {}
    for f in files:
        targeton = os.path.basename(f).replace(".tsv", "")
        frames[targeton] = pd.read_csv(f, sep="\t")
    return frames


def main():
    args = parse_args()
    runs = []
    for r in args.run:
        if ":" not in r:
            sys.exit(f"--run must be outdir:label, got '{r}'")
        outdir, label = r.split(":", 1)
        runs.append((label, load_run(outdir)))

    baseline_label, baseline_frames = runs[0]
    other_runs = runs[1:]
    if not other_runs:
        sys.exit("Need at least 2 --run entries (a baseline plus at least one comparison).")

    os.makedirs(args.outdir, exist_ok=True)

    flip_summary_rows = []
    flip_detail_rows = []

    common_targetons = set(baseline_frames.keys())
    for _, frames in other_runs:
        common_targetons &= set(frames.keys())
    if not common_targetons:
        sys.exit("No targeton files in common across the given runs -- check filenames match across output dirs.")

    for targeton in sorted(common_targetons):
        base_df = baseline_frames[targeton].reset_index(drop=True)
        n_rows = len(base_df)
        conseq = base_df[args.consequence_col] if args.consequence_col in base_df.columns else pd.Series(["NA"] * n_rows)
        base_tier = base_df[args.tier_col]

        id_cols_present = [c for c in [c.strip() for c in args.id_cols.split(",") if c.strip()] if c in base_df.columns]

        for label, frames in other_runs:
            cmp_df = frames[targeton].reset_index(drop=True)
            if len(cmp_df) != n_rows:
                print(f"WARNING: {targeton} has {n_rows} rows in '{baseline_label}' but {len(cmp_df)} in '{label}' -- skipping (row order comparison unsafe).")
                continue
            cmp_tier = cmp_df[args.tier_col]

            changed = base_tier != cmp_tier
            n_changed = int(changed.sum())
            flip_summary_rows.append({
                "targeton": targeton,
                "comparison": f"{baseline_label} -> {label}",
                "n_variants": n_rows,
                "n_tier_changed": n_changed,
                "pct_changed": round(100 * n_changed / n_rows, 2) if n_rows else 0.0,
            })

            if n_changed:
                idx = changed[changed].index
                for i in idx:
                    row = {
                        "targeton": targeton,
                        "row_index": i,
                        "consequence": conseq.iloc[i] if i < len(conseq) else "NA",
                        "comparison": f"{baseline_label} -> {label}",
                        f"tier_{baseline_label}": base_tier.iloc[i],
                        f"tier_{label}": cmp_tier.iloc[i],
                    }
                    for c in id_cols_present:
                        row[c] = base_df[c].iloc[i]
                    flip_detail_rows.append(row)

    flip_summary = pd.DataFrame(flip_summary_rows)
    flip_detail = pd.DataFrame(flip_detail_rows)

    summary_path = os.path.join(args.outdir, "anchor_tier_flip_summary.tsv")
    flip_summary.to_csv(summary_path, sep="\t", index=False)

    detail_path = os.path.join(args.outdir, "anchor_tier_flip_detail.tsv")
    flip_detail.to_csv(detail_path, sep="\t", index=False)

    conseq_path = os.path.join(args.outdir, "anchor_tier_flip_by_consequence.tsv")
    if not flip_detail.empty:
        by_conseq = flip_detail.groupby(["comparison", "consequence"]).size().reset_index(name="n_flipped")
        by_conseq.to_csv(conseq_path, sep="\t", index=False)
    else:
        pd.DataFrame(columns=["comparison", "consequence", "n_flipped"]).to_csv(conseq_path, sep="\t", index=False)

    print(f"Baseline: {baseline_label}\n")
    print(flip_summary.to_string(index=False))
    print()
    total_variants = flip_summary.drop_duplicates("targeton")["n_variants"].sum()
    for label, _ in other_runs:
        sub = flip_summary[flip_summary["comparison"].str.endswith(f"-> {label}")]
        total_changed = sub["n_tier_changed"].sum()
        print(f"TOTAL, {baseline_label} -> {label}: {total_changed}/{total_variants} variants changed tier ({100*total_changed/total_variants:.2f}%)")

    print()
    print("Wrote:")
    for p in [summary_path, detail_path, conseq_path]:
        print(f"  {p}")


if __name__ == "__main__":
    main()
