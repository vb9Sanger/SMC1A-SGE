#!/usr/bin/env python3
"""
merge_genedx_clinvar.py

Combines GeneDx's de novo SMC1A cohort data (patient-level HPO phenotype
terms + diagnostic status) with the ClinVar-derived per-variant summary
(clinvar_variants_summary.tsv, from sge_clinvar_intersect.py), giving a
patient-phenotype-informed view of each variant alongside its ClinVar
condition label and SGE anchor_tier functional call.

Input 1 --clinvar_summary: clinvar_variants_summary.tsv
    Needs: position, vcf_ref, vcf_alt (join key), plus whatever annotation
    columns you want carried through (condition, clnsig_norm, anchor_tier,
    Summary_Plot, HGVSc, HGVSp, Targeton_ID, clinvar_id, ...).

Input 2 --genedx: GeneDx TSV, one row per patient
    Needs: chrom, pos, ref, alt, SAMPLE_ID_anon, hpos (';'-separated HPO
    IDs), diagnostic_status (e.g. positive / possible / negative).

IMPORTANT — diagnostic_status handling:
    'negative' means this variant was found in that patient but was NOT
    deemed diagnostic/causal for their phenotype. By default their HPO
    terms are EXCLUDED from the pooled phenotype set for a variant (pass
    --include_negative_in_pool to override), since including them would
    misattribute an unrelated patient's presentation to the variant. They
    are still counted (genedx_n_negative) and their terms kept in a
    separate genedx_negative_only_hpos column for reference.

Output: one row per distinct variant (by position/ref/alt), covering the
UNION of variants in both inputs (outer join), with in_clinvar/in_genedx
flags -- so GeneDx variants not yet in ClinVar are visible too, not
silently dropped.

Optional --resolve_hpo_names queries the live HPO REST API
(https://ontology.jax.org/api/hp/terms/{id}) to resolve each HPO code to
its human-readable name (requires internet access; results are cached
locally so repeat runs don't re-query known codes), and applies a SMALL
STARTER keyword classifier (CDLS_KEYWORDS / DEE_KEYWORDS below) to flag
each variant's pooled phenotype as CdLS-suggestive / DEE-suggestive /
mixed / no match. This keyword list only contains terms independently
verified against HPO/OMIM sources at the time of writing -- review and
expand it yourself via https://hpo.jax.org/app before relying on it for
anything beyond a rough first pass.

Usage:
    python merge_genedx_clinvar.py \\
        --clinvar_summary clinvar_variants_summary.tsv \\
        --genedx gdx_dnm_180k_SMC1A_HPOs_diagnostic_status.tsv \\
        --output clinvar_genedx_merged.tsv \\
        --resolve_hpo_names
"""
import argparse
import json
import os
import sys
import time
import urllib.request
import urllib.error
import pandas as pd

HPO_API_URL = "https://ontology.jax.org/api/hp/terms/{}"

# ---- STARTER phenotype keyword lists ---------------------------------------
# Matched case-insensitively against RESOLVED HPO term NAMES (only populated
# if --resolve_hpo_names is used). This is a small illustrative starting
# point, NOT a validated phenotype classifier -- only terms independently
# verified against HPO/OMIM sources are included below. Everything else in
# your 600+ distinct HPO terms is unclassified until you expand this list
# using the HPO browser (https://hpo.jax.org/app) and the clinical
# literature on SMC1A-CdLS vs SMC1A-DEE85.
CDLS_KEYWORDS = [
    "synophrys",             # HP:0000664 -- verified (classic CdLS hallmark, OMIM #300590)
]
DEE_KEYWORDS = [
    "seizure",                # HP:0001250 -- verified
    "epileptic",
    "epilepsy",
    "developmental regression",  # HP:0002376 -- verified, present in this dataset
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clinvar_summary", required=True,
                    help="clinvar_variants_summary.tsv (from sge_clinvar_intersect.py)")
    p.add_argument("--genedx", required=True,
                    help="GeneDx patient-level TSV (chrom, pos, ref, alt, SAMPLE_ID_anon, "
                         "hpos, diagnostic_status)")
    p.add_argument("--output", default="clinvar_genedx_merged.tsv",
                    help="Output TSV path (default: clinvar_genedx_merged.tsv)")
    p.add_argument("--include_negative_in_pool", action="store_true",
                    help="Also fold diagnostic_status=='negative' patients' HPO terms into "
                         "the pooled phenotype set (default: excluded -- see docstring).")
    p.add_argument("--resolve_hpo_names", action="store_true",
                    help="Query the live HPO REST API to resolve HPO codes to names, and "
                         "apply the starter CdLS/DEE keyword classifier. Requires internet "
                         "access from wherever this script runs.")
    p.add_argument("--hpo_cache", default="hpo_term_cache.json",
                    help="Local cache file for resolved HPO id->name mappings "
                         "(default: hpo_term_cache.json, reused across runs)")
    return p.parse_args()


def aggregate_genedx(gdx: pd.DataFrame, include_negative: bool) -> pd.DataFrame:
    """One row per distinct (chrom,pos,ref,alt), aggregated across patients."""
    rows = []
    pool_statuses = {"positive", "possible", "negative"} if include_negative else {"positive", "possible"}

    for (chrom, pos, ref, alt), g in gdx.groupby(["chrom", "pos", "ref", "alt"]):
        status_lower = g["diagnostic_status"].str.lower()
        status_counts = status_lower.value_counts()

        pooled_hpos = set()
        for h in g.loc[status_lower.isin(pool_statuses), "hpos"].dropna():
            pooled_hpos.update(x for x in h.split(";") if x)

        negative_only_hpos = set()
        for h in g.loc[status_lower == "negative", "hpos"].dropna():
            negative_only_hpos.update(x for x in h.split(";") if x)

        rows.append({
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "genedx_n_patients":  len(g),
            "genedx_n_positive":  int(status_counts.get("positive", 0)),
            "genedx_n_possible":  int(status_counts.get("possible", 0)),
            "genedx_n_negative":  int(status_counts.get("negative", 0)),
            "genedx_sample_ids":  ";".join(g["SAMPLE_ID_anon"].astype(str)),
            "genedx_pooled_hpos": ";".join(sorted(pooled_hpos)),
            "genedx_negative_only_hpos": ";".join(sorted(negative_only_hpos)),
        })
    return pd.DataFrame(rows)


def resolve_hpo_names(hpo_ids: set, cache_path: str) -> dict:
    """Resolve a set of HPO IDs to human-readable names via the live HPO
    REST API, caching results locally. Returns {hpo_id: name}; unresolved
    IDs (network failure, 404, etc.) map to "" rather than raising."""
    cache = {}
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            cache = json.load(f)

    to_fetch = sorted(h for h in hpo_ids if h not in cache)
    if not to_fetch:
        print(f"All {len(hpo_ids)} HPO term(s) already in cache ({cache_path}).")
        return cache

    print(f"Resolving {len(to_fetch)} new HPO term name(s) via ontology.jax.org "
          f"({len(cache)} already cached)...")
    n_failed = 0
    for i, hid in enumerate(to_fetch):
        try:
            req = urllib.request.Request(
                HPO_API_URL.format(hid), headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                cache[hid] = data.get("name", "")
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError,
                json.JSONDecodeError) as e:
            cache[hid] = ""
            n_failed += 1
        if (i + 1) % 50 == 0:
            print(f"  ...{i + 1}/{len(to_fetch)}")
        time.sleep(0.05)  # be polite to a free public API

    with open(cache_path, "w") as f:
        json.dump(cache, f, indent=2, sort_keys=True)

    if n_failed:
        print(f"  WARNING: {n_failed}/{len(to_fetch)} lookups failed (network issue or "
              f"unrecognised ID) -- those will show as blank names.")
    return cache


def score_phenotype(hpo_ids: list, name_lookup: dict):
    """Return (cdls_keyword_hits, dee_keyword_hits) lists for one variant's
    pooled HPO ids, using the STARTER keyword lists (see module docstring)."""
    names_blob = " | ".join(name_lookup.get(h, "").lower() for h in hpo_ids if h)
    cdls_hits = [kw for kw in CDLS_KEYWORDS if kw in names_blob]
    dee_hits  = [kw for kw in DEE_KEYWORDS if kw in names_blob]
    return cdls_hits, dee_hits


def main():
    args = parse_args()

    cvs = pd.read_csv(args.clinvar_summary, sep="\t")
    for c in ["position", "vcf_ref", "vcf_alt"]:
        if c not in cvs.columns:
            sys.exit(f"ERROR: missing column '{c}' in {args.clinvar_summary}")

    gdx = pd.read_csv(args.genedx, sep="\t")
    for c in ["chrom", "pos", "ref", "alt", "SAMPLE_ID_anon", "hpos", "diagnostic_status"]:
        if c not in gdx.columns:
            sys.exit(f"ERROR: missing column '{c}' in {args.genedx}")

    print(f"ClinVar summary: {len(cvs)} row(s)")
    print(f"GeneDx: {len(gdx)} patient-row(s), {gdx['SAMPLE_ID_anon'].nunique()} distinct patients")

    gdx_agg = aggregate_genedx(gdx, args.include_negative_in_pool)
    print(f"GeneDx aggregated to {len(gdx_agg)} distinct variant(s)")

    # De-dup the ClinVar summary by variant (position/vcf_ref/vcf_alt): the
    # same variant can appear multiple times if overlapping targetons both
    # cover it.
    cvs_dedup = cvs.drop_duplicates(subset=["position", "vcf_ref", "vcf_alt"]).copy()
    if len(cvs_dedup) < len(cvs):
        print(f"  ({len(cvs) - len(cvs_dedup)} duplicate row(s) dropped from overlapping targetons)")

    merged = cvs_dedup.merge(
        gdx_agg,
        left_on=["position", "vcf_ref", "vcf_alt"], right_on=["pos", "ref", "alt"],
        how="outer", indicator=True,
    )
    merged["in_clinvar"] = merged["_merge"].isin(["left_only", "both"])
    merged["in_genedx"]  = merged["_merge"].isin(["right_only", "both"])
    merged.drop(columns=["_merge"], inplace=True)

    for c in ["genedx_n_patients", "genedx_n_positive", "genedx_n_possible", "genedx_n_negative"]:
        merged[c] = merged[c].fillna(0).astype(int)
    for c in ["genedx_sample_ids", "genedx_pooled_hpos", "genedx_negative_only_hpos"]:
        merged[c] = merged[c].fillna("")

    n_both         = int((merged["in_clinvar"] & merged["in_genedx"]).sum())
    n_clinvar_only = int((merged["in_clinvar"] & ~merged["in_genedx"]).sum())
    n_genedx_only  = int((~merged["in_clinvar"] & merged["in_genedx"]).sum())
    print(f"\nMatched (in both ClinVar and GeneDx): {n_both}")
    print(f"ClinVar-only: {n_clinvar_only}")
    print(f"GeneDx-only (NOT currently in your ClinVar summary): {n_genedx_only}")

    if args.resolve_hpo_names:
        all_ids = set()
        for h in merged["genedx_pooled_hpos"]:
            if h:
                all_ids.update(h.split(";"))
        name_lookup = resolve_hpo_names(all_ids, args.hpo_cache)

        def names_for(h):
            return "" if not h else ";".join(name_lookup.get(x, x) for x in h.split(";"))
        merged["genedx_pooled_hpo_names"] = merged["genedx_pooled_hpos"].apply(names_for)

        cdls_col, dee_col, signal_col = [], [], []
        for h in merged["genedx_pooled_hpos"]:
            if not h:
                cdls_col.append(""); dee_col.append(""); signal_col.append("no GeneDx phenotype data")
                continue
            cdls_hits, dee_hits = score_phenotype(h.split(";"), name_lookup)
            cdls_col.append(";".join(cdls_hits))
            dee_col.append(";".join(dee_hits))
            if cdls_hits and dee_hits:
                signal_col.append("mixed signal")
            elif cdls_hits:
                signal_col.append("CdLS-suggestive")
            elif dee_hits:
                signal_col.append("DEE-suggestive")
            else:
                signal_col.append("no keyword match (expand CDLS_KEYWORDS/DEE_KEYWORDS)")
        merged["cdls_keyword_hits"] = cdls_col
        merged["dee_keyword_hits"]  = dee_col
        merged["genedx_phenotype_signal"] = signal_col

        print("\ngenedx_phenotype_signal counts (variants with >=1 GeneDx patient only):")
        has_gdx = merged[merged["in_genedx"]]
        print(has_gdx["genedx_phenotype_signal"].value_counts().to_string())

        if "condition" in merged.columns:
            both = merged[merged["in_clinvar"] & merged["in_genedx"]]
            if not both.empty:
                print("\ngenedx_phenotype_signal x ClinVar 'condition' (matched variants only):")
                print(pd.crosstab(both["condition"], both["genedx_phenotype_signal"]).to_string())

    merged.to_csv(args.output, sep="\t", index=False)
    print(f"\nWrote: {args.output}")


if __name__ == "__main__":
    main()
