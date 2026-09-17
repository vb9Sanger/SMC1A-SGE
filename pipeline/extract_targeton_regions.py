#!/usr/bin/env python3
"""
extract_targeton_regions.py

Scans a directory of SGE meta files and builds a lookup table mapping
each Targeton_ID to its chromosomal region.

For each targeton the script expects two files:
  <prefix>_meta_consequences.tsv   — contains Targeton_ID column
  <prefix>_meta.csv                — contains ref_chr, ref_start, ref_end columns

Output: a TSV file with columns:
  Targeton_ID | chrom | start | end | source_prefix

Usage:
  python extract_targeton_regions.py --input_dir /path/to/meta/files \
                                     --output targeton_regions.tsv

Optional:
  --meta_suffix         suffix for meta files          (default: _meta.csv)
  --consequences_suffix suffix for consequence files   (default: _meta_consequences.tsv)
"""

import argparse
import os
import sys
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input_dir", required=True,
                        help="Directory containing meta and meta_consequences files")
    parser.add_argument("--output", default="targeton_regions.tsv",
                        help="Output TSV path (default: targeton_regions.tsv)")
    parser.add_argument("--meta_suffix", default="_meta.csv",
                        help="Suffix identifying meta files (default: _meta.csv)")
    parser.add_argument("--consequences_suffix", default="_meta_consequences.tsv",
                        help="Suffix identifying consequence files (default: _meta_consequences.tsv)")
    return parser.parse_args()


def find_file_pairs(input_dir, meta_suffix, cons_suffix):
    """
    Return a list of (prefix, meta_path, consequences_path) tuples.
    Pairs are identified by finding all *_meta_consequences.tsv files
    and looking for the matching *_meta.csv alongside them.
    """
    pairs = []
    missing_meta = []

    for fname in sorted(os.listdir(input_dir)):
        if not fname.endswith(cons_suffix):
            continue

        prefix = fname[: -len(cons_suffix)]
        cons_path = os.path.join(input_dir, fname)
        meta_path = os.path.join(input_dir, prefix + meta_suffix)

        if not os.path.isfile(meta_path):
            missing_meta.append(fname)
            continue

        pairs.append((prefix, meta_path, cons_path))

    if missing_meta:
        print(f"WARNING: {len(missing_meta)} consequence file(s) had no matching meta file:",
              file=sys.stderr)
        for f in missing_meta:
            print(f"  {f}", file=sys.stderr)

    return pairs


def extract_targeton_region(prefix, meta_path, cons_path):
    """
    Extract unique Targeton_ID(s) and the region (chrom, start, end)
    for one targeton pair.

    Returns a list of dicts (one per unique Targeton_ID found).
    """
    # --- Targeton_ID from consequences file ---
    try:
        cons_df = pd.read_csv(cons_path, sep="\t", usecols=["Targeton_ID"])
    except Exception as e:
        print(f"ERROR reading {cons_path}: {e}", file=sys.stderr)
        return []

    targeton_ids = cons_df["Targeton_ID"].dropna().unique().tolist()
    if not targeton_ids:
        print(f"WARNING: no Targeton_ID values found in {cons_path}", file=sys.stderr)
        return []

    # --- Region from meta file ---
    try:
        meta_df = pd.read_csv(meta_path, usecols=["ref_chr", "ref_start", "ref_end"])
    except Exception as e:
        print(f"ERROR reading {meta_path}: {e}", file=sys.stderr)
        return []

    meta_df = meta_df.dropna(subset=["ref_chr", "ref_start", "ref_end"])
    if meta_df.empty:
        print(f"WARNING: no region data found in {meta_path}", file=sys.stderr)
        return []

    # Use the first row's region (all rows in a targeton share the same region)
    chrom = meta_df["ref_chr"].iloc[0]
    start = int(meta_df["ref_start"].iloc[0])
    end   = int(meta_df["ref_end"].iloc[0])

    # Sanity check: warn if multiple distinct regions appear (shouldn't happen)
    if meta_df[["ref_chr", "ref_start", "ref_end"]].nunique().max() > 1:
        print(f"WARNING: multiple distinct regions found in {meta_path} — using first row",
              file=sys.stderr)

    records = []
    for tid in targeton_ids:
        records.append({
            "Targeton_ID":   tid,
            "chrom":         chrom,
            "start":         start,
            "end":           end,
            "source_prefix": prefix,
        })
    return records


def main():
    args = parse_args()

    if not os.path.isdir(args.input_dir):
        sys.exit(f"ERROR: --input_dir does not exist: {args.input_dir}")

    pairs = find_file_pairs(args.input_dir, args.meta_suffix, args.consequences_suffix)
    if not pairs:
        sys.exit("ERROR: no file pairs found. Check --input_dir and suffix arguments.")

    print(f"Found {len(pairs)} targeton file pair(s) in {args.input_dir}")

    all_records = []
    for prefix, meta_path, cons_path in pairs:
        records = extract_targeton_region(prefix, meta_path, cons_path)
        all_records.extend(records)

    if not all_records:
        sys.exit("ERROR: no records extracted. Check file formats.")

    out_df = pd.DataFrame(all_records, columns=["Targeton_ID", "chrom", "start", "end",
                                                 "source_prefix"])

    # Deduplicate in case the same Targeton_ID appears across multiple files
    n_before = len(out_df)
    out_df = out_df.drop_duplicates(subset=["Targeton_ID", "chrom", "start", "end"])
    n_dupes = n_before - len(out_df)
    if n_dupes:
        print(f"Removed {n_dupes} duplicate row(s)")

    out_df = out_df.sort_values(["chrom", "start", "Targeton_ID"]).reset_index(drop=True)
    out_df.to_csv(args.output, sep="\t", index=False)

    print(f"Written {len(out_df)} targeton(s) → {args.output}")
    print(out_df.to_string(index=False))


if __name__ == "__main__":
    main()
