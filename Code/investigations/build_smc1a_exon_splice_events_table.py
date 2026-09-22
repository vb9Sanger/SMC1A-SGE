#!/usr/bin/env python3
"""
build_smc1a_exon_splice_events_table.py

Builds the exon-level alternative-splicing annotation table that
cassette_exon_consequence_analysis.py needs as its --exon_events_tsv input
(Targeton_ID, EXON, events_within_exon, events_adjacent_boundary,
design_status) -- cross-referencing UCSC's "Alt Events" track (table
knownAlt: altPromoter, altFivePrime, altThreePrime, altFinish, cassetteExon,
retainedIntron, bleedingExon -- derived from comparing GENCODE transcript
structures) against every SMC1A exon's genomic coordinates.

For each exon, every UCSC event overlapping it is classified as:
  - "within exon"       if the event lies entirely inside the exon's
                         [exon_g_start, exon_g_end] genomic span.
  - "adjacent boundary"  if it overlaps the exon but extends past one of its
                         boundaries (e.g. a retainedIntron or bleedingExon
                         spanning the exon-intron junction).

A targeton counts as "cassette" (for cassette_exon_consequence_analysis.py)
only via events_within_exon containing "cassetteExon" -- an
events_adjacent_boundary annotation alone does not count, since that reflects
a neighbouring intron/exon event rather than the assayed exon itself being
skippable.

design_status is "assayed, results available" if a per-targeton anchor-tier
result file exists for that targeton in --results_dir, else
"designed, results pending".

Inputs:
  --targetons_tsv   curation_attempt2/data/reference/smc1a_targetons.tsv
                     (needs targeton_id, exon, exon_g_start, exon_g_end)
  --results_dir     Directory of per-targeton anchor-tier result TSVs (e.g.
                     classifier_run/classifier_results/D4_ref/) -- used only
                     to check which targetons have a result file.
  --genome          UCSC genome assembly (default: hg38)
  --region_chrom / --region_start / --region_end
                     Genomic region to query the knownAlt track over
                     (default: the SMC1A locus with a 20kb margin).

Usage:
    python build_smc1a_exon_splice_events_table.py \\
        --targetons_tsv smc1a_targetons.tsv \\
        --results_dir classifier_results/D4_ref/ \\
        --output SMC1A_exon_splice_events_full.tsv
"""

import argparse
import json
import sys
import urllib.request

import pandas as pd


def fetch_known_alt(genome, chrom, start, end):
    url = (f"https://api.genome.ucsc.edu/getData/track?genome={genome};"
           f"track=knownAlt;chrom={chrom};start={start};end={end}")
    with urllib.request.urlopen(url) as resp:
        data = json.load(resp)
    events = data.get("knownAlt", [])
    if not events:
        sys.exit(f"ERROR: no knownAlt events returned for {chrom}:{start}-{end} on {genome}")
    return events


def classify_events(exon_start, exon_end, events):
    """exon_start/exon_end are 1-based inclusive (as in smc1a_targetons.tsv).
    UCSC's knownAlt coordinates are BED-style 0-based half-open, so a
    feature's 1-based inclusive span is (chromStart + 1) to chromEnd."""
    within, adjacent = [], []
    for ev in events:
        ev_start, ev_end, name = ev["chromStart"] + 1, ev["chromEnd"], ev["name"]
        if ev_end < exon_start or ev_start > exon_end:
            continue  # no overlap
        if ev_start >= exon_start and ev_end <= exon_end:
            within.append(name)
        else:
            adjacent.append(name)
    # de-duplicate, keep first-seen order
    within = list(dict.fromkeys(within))
    adjacent = list(dict.fromkeys(adjacent))
    return ", ".join(within), ", ".join(adjacent)


def has_results(targeton_id, exon, results_dir):
    import glob
    import os
    pattern = os.path.join(results_dir, f"{targeton_id}_exon{exon}_*.tsv")
    return len(glob.glob(pattern)) > 0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targetons_tsv", required=True)
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--genome", default="hg38")
    parser.add_argument("--region_chrom", default="chrX")
    parser.add_argument("--region_start", type=int, default=53370000)
    parser.add_argument("--region_end", type=int, default=53440000)
    args = parser.parse_args()

    targetons = pd.read_csv(args.targetons_tsv, sep="\t")
    events = fetch_known_alt(args.genome, args.region_chrom, args.region_start, args.region_end)
    print(f"Fetched {len(events)} knownAlt event(s) from UCSC ({args.genome}, "
          f"{args.region_chrom}:{args.region_start}-{args.region_end})", file=sys.stderr)

    rows = []
    for _, r in targetons.iterrows():
        within, adjacent = classify_events(int(r["exon_g_start"]), int(r["exon_g_end"]), events)
        status = ("assayed, results available"
                   if has_results(r["targeton_id"], int(r["exon"]), args.results_dir)
                   else "designed, results pending")
        rows.append({
            "EXON": int(r["exon"]),
            "Targeton_ID": r["targeton_id"],
            "events_within_exon": within,
            "events_adjacent_boundary": adjacent,
            "design_status": status,
        })

    out = pd.DataFrame(rows).sort_values(["EXON", "Targeton_ID"])
    out.to_csv(args.output, sep="\t", index=False)

    n_cassette = out["events_within_exon"].str.contains("cassetteExon").sum()
    n_none = ((out["events_within_exon"] == "") & (out["events_adjacent_boundary"] == "")).sum()
    print(f"Wrote {args.output} ({len(out)} targeton row(s))", file=sys.stderr)
    print(f"{n_cassette}/{len(out)} exons carry a cassetteExon (within-exon) annotation", file=sys.stderr)
    print(f"{n_none}/{len(out)} exons have no alt-splicing annotation at all", file=sys.stderr)


if __name__ == "__main__":
    main()
