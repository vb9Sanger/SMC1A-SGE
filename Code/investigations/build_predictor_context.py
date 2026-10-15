#!/usr/bin/env python3
"""
build_predictor_context.py

Assembles per-variant orthogonal annotations for comparison against the
assay's own functional calls: pathogenicity prediction (CADD), evolutionary
constraint (phyloP), structural burial (relative solvent accessibility) and
predicted folding stability (ddG).

Four sources, four different kinds of evidence:

  CADD    already in the merged in-silico table -- a trained pathogenicity
          predictor, i.e. partly circular with clinical labels
  phyloP  per-base 100-way vertebrate conservation, fetched from the UCSC API
          and cached. Evolutionary, independent of any clinical annotation
  RSA     relative solvent accessibility, retrieved per residue from AlphaSync
          (Lang & Babu), which computes it with DSSP from the AlphaFold model
          and normalises by the residue's theoretical maximum (Tien et al.
          2013). Structural, and likewise independent of clinical labels
  ddG     ThermoMPNN's predicted change in folding free energy, from the
          site-saturation scan of the same AlphaFold model. Structural, and
          the only annotation here that covers the whole missense library

The point of pairing them: CADD agreeing with the assay is weak evidence
because both are tuned on overlapping notions of damage, whereas phyloP, RSA
and ddG know nothing about pathogenicity at all.

Coordinate notes. phyloP is per genomic base, so it applies to every variant.
RSA is per residue, so it applies only to variants with a protein position.
pLDDT is carried through per residue, because RSA from a low-confidence region
of a predicted model is not trustworthy -- 12.4% of SMC1A residues sit below
pLDDT 70 and can be excluded or flagged downstream.

A note on which structure each annotation comes from. AlphaSync's stored
residue features for Q14683 were computed from an earlier AlphaFold DB release
than the model_v6 one its own download endpoint now points at, and the two
differ materially across residues 240-268 (mean |dRSA| 0.21 there against 0.04
elsewhere). That block holds 713 of 27,937 missense variants and does not move
the result -- the depleting-vs-no-impact AUC is 0.704 on AlphaSync values
against 0.710 on the current model -- but RSA and ddG are therefore not
guaranteed to describe identical coordinates, and the RSA column should not be
read as a measurement on the v6 model.

Usage:
    python build_predictor_context.py \\
        --insilico classifier_run/insilico_concordance/merged_variant_data.tsv \\
        --assay_summary gnomAD/gnomad_intersect_results/sge_gnomad_summary.tsv \\
        --ddg classifier_run/ddg_analysis/ThermoMPNN_inference_SMC1A_Q14683.csv \\
        --uniprot Q14683 \\
        --region chrX:53374149-53422654 \\
        --outdir classifier_run/plots/predictor_context
"""
import argparse
import json
import os
import re
import ssl
import subprocess
import time
import urllib.error
import urllib.request

import pandas as pd

UCSC = ("https://api.genome.ucsc.edu/getData/track"
        "?genome=hg38;track={track};chrom={chrom};start={start};end={end}")

ALPHASYNC = "https://alphasync.stjude.org/api/v1/protein/{acc}"

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def fetch_phylop(chrom, start, end, track, cache, chunk=8000, pause=0.4):
    """Per-base phyloP over [start, end), cached to disk.

    Chunked because the API caps response size; the cache means the figure can
    be regenerated offline and that a rerun does not hammer UCSC.
    """
    if os.path.exists(cache):
        df = pd.read_csv(cache, sep="\t")
        print(f"phyloP: {len(df):,} bases from cache {cache}")
        return df
    rows = []
    for s in range(start, end, chunk):
        e = min(s + chunk, end)
        url = UCSC.format(track=track, chrom=chrom, start=s, end=e)
        with urllib.request.urlopen(url, timeout=120) as r:
            payload = json.load(r)
        vals = payload.get(track, [])
        if isinstance(vals, dict):            # some tracks nest by chrom
            vals = vals.get(chrom, [])
        for v in vals:
            rows.append((int(v["start"]), float(v["value"])))
        print(f"  {chrom}:{s:,}-{e:,}  {len(vals):,} bases")
        time.sleep(pause)
    df = pd.DataFrame(rows, columns=["pos0", "phyloP"]).drop_duplicates("pos0")
    # the API is 0-based half-open; variant tables are 1-based
    df["position"] = df.pos0 + 1
    df = df[["position", "phyloP"]]
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    df.to_csv(cache, sep="\t", index=False)
    print(f"phyloP: {len(df):,} bases fetched, cached to {cache}")
    return df


def get_json(url, timeout=180):
    """GET a JSON document, falling back to curl.

    AlphaSync negotiates a TLS version that Python builds older than 3.10
    refuse, which is what the project's curation virtualenv happens to be;
    curl on the same machine completes the same request. Trying the standard
    library first keeps the script dependency-free where it can be.
    """
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.load(r)
    except (urllib.error.URLError, ssl.SSLError) as exc:
        print(f"  urllib failed ({exc}); retrying with curl")
        out = subprocess.run(
            ["curl", "-sS", "-L", "-m", str(timeout),
             "-H", "Accept: application/json", url],
            capture_output=True, text=True, check=True)
        return json.loads(out.stdout)


def fetch_alphasync_rsa(acc, cache):
    """Per-residue RSA and pLDDT for one UniProt accession, from AlphaSync.

    Taken from the published resource rather than recomputed here so the
    values are citable and the normalisation is someone else's documented
    choice: AlphaSync runs DSSP 4.2.2.1 on the AlphaFold model and divides by
    the Tien et al. 2013 theoretical maximum per residue type.

    Values are passed through exactly as served, including the occasional
    residue whose RSA exceeds 1 -- a published number that has been altered is
    no longer the published number.
    """
    if os.path.exists(cache):
        df = pd.read_csv(cache, sep="\t")
        print(f"RSA: {len(df)} residues from cache {cache}")
    else:
        payload = get_json(ALPHASYNC.format(acc=acc))
        df = pd.DataFrame(payload)[["site", "aa", "asa", "relasa", "plddt"]]
        df = df.rename(columns={"site": "aa_pos", "aa": "resname",
                                "asa": "sasa", "relasa": "rsa"})
        df = df.sort_values("aa_pos").reset_index(drop=True)
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        df.to_csv(cache, sep="\t", index=False)
        print(f"RSA: {len(df)} residues fetched from AlphaSync, cached to {cache}")
    n_over = int((df.rsa > 1).sum())
    print(f"  median {df.rsa.median():.3f}, {n_over} residue(s) above 1.0, "
          f"{100 * (df.plddt < 70).mean():.1f}% below pLDDT 70")
    return df


def load_ddg(ddg_csv):
    """ThermoMPNN site-saturation scan, keyed by (residue, wild type, mutant).

    The scan's `position` is zero-indexed against the model sequence; every
    other table here is 1-based, so it is shifted once on load rather than at
    each use.
    """
    d = pd.read_csv(ddg_csv)
    d["aa_pos"] = d.position.astype(int) + 1
    d = d[["aa_pos", "wildtype", "mutation", "ddG_pred"]]
    print(f"ddG: {len(d):,} (position, substitution) pairs from {ddg_csv}")
    return d


def attach_ddg(v, ddg):
    """Join predicted ddG onto variants by parsing their protein HGVS.

    Only simple single-residue substitutions are matched. The wild-type
    residue is carried into the key, so a variant whose annotated wild type
    disagrees with the model sequence fails to join and arrives as missing
    rather than silently taking another substitution's score.
    """
    p = v.HGVSp.fillna("").str.replace(r"^.*p\.", "", regex=True)
    m = p.str.extract(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")
    v = v.copy()
    v["_wt"] = m[0].map(AA3_TO_1)
    v["_pos"] = pd.to_numeric(m[1], errors="coerce").astype("Int64")
    v["_mut"] = m[2].map(AA3_TO_1)
    v = v.merge(ddg.rename(columns={"aa_pos": "_pos", "wildtype": "_wt",
                                    "mutation": "_mut", "ddG_pred": "ddg"}),
                on=["_pos", "_wt", "_mut"], how="left")
    mis = v[v.Summary_Consequence == "Missense_Variant"]
    print(f"  matched {mis.ddg.notna().sum():,} of {len(mis):,} missense variants "
          f"({100 * mis.ddg.notna().mean():.1f}%)")
    return v.drop(columns=["_wt", "_pos", "_mut"])


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--insilico", required=True)
    ap.add_argument("--assay_summary", required=True)
    ap.add_argument("--ddg", required=True)
    ap.add_argument("--uniprot", default="Q14683")
    ap.add_argument("--region", required=True, help="chrX:start-end, 1-based")
    ap.add_argument("--track", default="phyloP100way")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    chrom, span = a.region.split(":")
    start, end = (int(x) for x in span.split("-"))

    phylop = fetch_phylop(chrom, start - 1, end, a.track,
                          os.path.join(a.outdir, f"{a.track}_{chrom}_{start}_{end}.tsv"))
    rsa = fetch_alphasync_rsa(a.uniprot,
                              os.path.join(a.outdir, "rsa_per_residue.tsv"))
    ddg = load_ddg(a.ddg)

    # ---- variants: calls and consequence from the assay summary, CADD from
    # the merged in-silico table. Joined on oligo_name, which is what the
    # in-silico merge itself keys on.
    g = pd.read_csv(a.assay_summary, sep="\t", dtype=str, low_memory=False)
    g["position"] = pd.to_numeric(g.position, errors="coerce").astype("Int64")
    g["aa_pos"] = pd.to_numeric(g.Protein_position, errors="coerce").astype("Int64")

    ins = pd.read_csv(a.insilico, sep="\t", dtype=str, low_memory=False)
    ins["CADD"] = pd.to_numeric(ins.CADD_num, errors="coerce")
    ins["REVEL"] = pd.to_numeric(ins.REVEL_num, errors="coerce")
    ins = ins[["oligo_name", "CADD", "REVEL"]].drop_duplicates("oligo_name")

    v = g[["oligo_name", "Targeton_ID", "position", "aa_pos", "anchor_tier",
           "Summary_Consequence", "HGVSp", "pos_adj_log2FoldChange_raw"]].copy()
    v["lfc"] = pd.to_numeric(v.pos_adj_log2FoldChange_raw, errors="coerce")
    v.drop(columns="pos_adj_log2FoldChange_raw", inplace=True)
    v = v.merge(ins, on="oligo_name", how="left")
    v = v.merge(phylop, on="position", how="left")
    v = v.merge(rsa[["aa_pos", "rsa", "plddt", "resname"]], on="aa_pos", how="left")
    v = attach_ddg(v, ddg)

    # SNV flag. CADD scores single-nucleotide substitutions; this library is
    # mostly multi-nucleotide codon replacements, which CADD does not cover, so
    # the CADD panel necessarily describes a subset. Recorded per variant so
    # the figure can state its own denominator instead of implying full
    # coverage.
    g2 = pd.read_csv(a.assay_summary, sep="\t", dtype=str, low_memory=False)
    g2["is_snv"] = (g2.vcf_ref.str.len() == 1) & (g2.vcf_alt.str.len() == 1)
    v = v.merge(g2[["oligo_name", "is_snv"]].drop_duplicates("oligo_name"),
                on="oligo_name", how="left")

    # Per-residue conservation, as well as per-base. phyloP at the variant's
    # own base is dominated by which codon position it hits -- a third-base
    # change sits at a far less conserved site than a second-base one -- which
    # is noise when the comparison is against a residue-level functional call.
    # Averaging over the codon's assayed bases gives the residue's constraint.
    by_res = (v.dropna(subset=["aa_pos", "phyloP"])
              .drop_duplicates(["aa_pos", "position"])
              .groupby("aa_pos")
              .agg(phylop_residue=("phyloP", "mean"),
                   n_codon_bases=("position", "nunique")).reset_index())
    v = v.merge(by_res, on="aa_pos", how="left")
    print(f"per-residue phyloP: {len(by_res)} residues, "
          f"median {by_res.n_codon_bases.median():.0f} assayed bases each")

    out = os.path.join(a.outdir, "variant_context.tsv")
    v.to_csv(out, sep="\t", index=False)
    print(f"\nWrote {out} ({len(v):,} rows)")

    print("\n--- annotation coverage, by consequence class")
    cov = v.groupby("Summary_Consequence").agg(
        n=("oligo_name", "size"),
        CADD=("CADD", lambda s: f"{100 * s.notna().mean():.0f}%"),
        phyloP=("phyloP", lambda s: f"{100 * s.notna().mean():.0f}%"),
        RSA=("rsa", lambda s: f"{100 * s.notna().mean():.0f}%"),
        ddG=("ddg", lambda s: f"{100 * s.notna().mean():.0f}%"))
    print(cov.to_string())
    print("\n--- missense only, by assay call")
    mis = v[v.Summary_Consequence == "Missense_Variant"]
    print(mis.groupby("anchor_tier").agg(
        n=("oligo_name", "size"),
        n_CADD=("CADD", "count"), CADD=("CADD", "median"),
        phyloP_base=("phyloP", "median"),
        phyloP_residue=("phylop_residue", "median"),
        RSA=("rsa", "median"),
        n_ddG=("ddg", "count"), ddG=("ddg", "median")).round(3).to_string())
    print("\n--- does restricting to SNVs (CADD's own scope) change RSA/phyloP?")
    for lab, sub in [("all missense", mis), ("missense SNVs", mis[mis.is_snv == True])]:
        t = sub.groupby("anchor_tier").agg(
            n=("oligo_name", "size"), RSA=("rsa", "median"),
            phyloP_residue=("phylop_residue", "median")).round(3)
        print(f"  {lab}:")
        print("    " + t.to_string().replace("\n", "\n    "))


if __name__ == "__main__":
    main()
