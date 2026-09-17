#!/usr/bin/env python3
"""Ingest the DDD/DECIPHER de novo SMC1A cohort into the observation ledger.

Input
-----
`genedx/DDD_SMC1A_variants_with_phenotype.tsv` -- 15 de novo SMC1A observations
formed by joining the DECIPHER variant extract to the DECIPHER phenotype
extract on `patient_id` (one-to-one, identical identifier sets, no duplicates).

Why this cohort is straightforward on coordinates
-------------------------------------------------
Unlike every literature source, DECIPHER supplies `assembly`, `chr`, `start`,
`ref_allele` and `alt_allele` directly, all 15 rows on **GRCh38**, with
VCF-style anchored alleles for indels (`GCCAT>G` for a 4 bp deletion). That is
already the convention the registry stores, so `variant_key` is built from the
supplied coordinate rather than derived from a cDNA description, and the
HGVS-versus-VCF normalisation hazard of METHODS 4A.12 does not arise.

`user_hgvs` is empty for all 15 rows, so c. and p. descriptions are obtained
from VEP on the MANE transcript rather than taken from the source. One row is
reported against `ENST00000375340` rather than MANE `ENST00000322213`; because
the join key is genomic this affects only the source's own annotation, not the
variant's identity, and the MANE description is derived here regardless.

Disease attribution
-------------------
DECIPHER supplies a pathogenicity call but no syndrome-level diagnosis. The
`associated_condition` column carries the attribution assigned from each
patient's DECIPHER clinical notes under the criteria recorded in
docs/ASSOCIATED_CONDITION_LABELLING.md section 4. An empty cell is a real value
meaning the notes do not name either syndrome; such rows keep a classification
and stay in the all-variants list, but carry no disease and cannot enter either
analysis arm.

Verification checks
-------------------
  V1  every row is SMC1A and `sequence_variant`
  V2  assembly is GRCh38 for every row, so no liftover is required
  V3  the supplied reference allele matches the GRCh38 reference base(s)
  V4  VEP returns a consequence on the MANE transcript
  V5  inheritance is de novo for every row (the cohort is DNM-ascertained)
  V6  every `associated_condition` value is one this script knows how to map
  V7  any DEE85 attribution is female (the premise in the labelling document)
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import urllib.parse
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smc1a_lib as L  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, "data", "reference")
INTERIM = os.path.join(ROOT, "data", "interim")
DOCS = os.path.join(ROOT, "docs", "verification")
SRC = os.path.join(ROOT, "genedx", "DDD_SMC1A_variants_with_phenotype.tsv")

http = L.CachedHTTP(os.path.join(ROOT, "data", "cache", "ensembl"))
MODEL = L.TranscriptModel.from_json(os.path.join(REF, "smc1a_transcript_model.json"))
REFSEQ = L.ReferenceSequence(os.path.join(REF, "smc1a_chrX_region.fa"))
TARGETONS = L.TargetonMap.from_tsv(os.path.join(REF, "smc1a_targetons.tsv"))

report: list[str] = []


def say(m: str = "") -> None:
    print(m)
    report.append(m)


# --------------------------------------------------------------------------
# DECIPHER pathogenicity -> registry classification
#
# DECIPHER's values are curator-assigned clinical calls. They map directly onto
# the registry vocabulary; `uncertain` maps to VUS and is therefore ineligible
# for the curated set regardless of any disease attribution.
# --------------------------------------------------------------------------
PATH_MAP = {
    "pathogenic": ("P", "DECIPHER curator classification: pathogenic"),
    "likely_pathogenic": ("LP", "DECIPHER curator classification: likely pathogenic"),
    "uncertain": ("VUS", "DECIPHER curator classification: uncertain significance"),
}

# --------------------------------------------------------------------------
# Disease attribution from `associated_condition` (see the labelling document).
#
#   CdLS   <- DECIPHER clinical notes name Cornelia de Lange syndrome
#   DEE85  <- the notes instead name holoprosencephaly or Rett syndrome,
#             and the patient is female
#
# Nothing is inferred here; this maps the column rather than re-deriving it.
# --------------------------------------------------------------------------
ASSOC_MAP = {
    "CdLS": ("Cornelia de Lange syndrome", "CdLS2",
             "assigned from the patient's DECIPHER clinical notes, which name "
             "Cornelia de Lange syndrome "
             "(docs/ASSOCIATED_CONDITION_LABELLING.md section 4)"),
    "DEE85": ("SMC1A-related developmental and epileptic encephalopathy", "DEE85",
              "assigned from the patient's DECIPHER clinical notes, which name "
              "holoprosencephaly or a Rett-like presentation, with female sex "
              "confirmed (docs/ASSOCIATED_CONDITION_LABELLING.md section 4)"),
}

SEX_MAP = {"female": "female", "male": "male"}


def vep_region(pos: int, ref: str, alt: str):
    """Consequence, c. and p. on the MANE transcript from a plus-strand allele.

    Submitted as an HGVS genomic description for consistency with
    03_resolve_variants.py and 04_ingest_genedx.py.
    """
    ref, alt = ref.upper(), alt.upper()
    end = pos + len(ref) - 1
    if len(ref) == 1 and len(alt) == 1:
        desc = f"{L.GENOMIC_ACC}:g.{pos}{ref}>{alt}"
    elif len(alt) == 0 or alt == "-":
        desc = (f"{L.GENOMIC_ACC}:g.{pos}del" if len(ref) == 1
                else f"{L.GENOMIC_ACC}:g.{pos}_{end}del")
    else:
        desc = (f"{L.GENOMIC_ACC}:g.{pos}delins{alt}" if len(ref) == 1
                else f"{L.GENOMIC_ACC}:g.{pos}_{end}delins{alt}")
    url = (f"{L.ENSEMBL_REST}/vep/human/hgvs/{urllib.parse.quote(desc)}"
           "?content-type=application/json;numbers=1;hgvs=1;mane=1")
    d, rec = http.get_json(url)
    if not d or not isinstance(d, list) or not d:
        return {}, f"vep failed HTTP {rec['status']} for {desc}"
    tc = [t for t in d[0].get("transcript_consequences", [])
          if t.get("transcript_id") == L.TRANSCRIPT_NOVER]
    if not tc:
        return {}, f"no consequence on the MANE transcript for {desc}"
    t = tc[0]
    terms = t.get("consequence_terms", [])
    hc, hp = L.clean(t.get("hgvsc")), L.clean(t.get("hgvsp"))
    return {
        "consequence_terms": ";".join(terms),
        "consequence_class": L.consequence_class(terms),
        "hgvs_c_mane": hc.split(":", 1)[1] if ":" in hc else "",
        "hgvs_p_mane": hp.split(":", 1)[1] if ":" in hp else "",
        "vep_exon": L.clean(t.get("exon")).split("/")[0],
        "vep_intron": L.clean(t.get("intron")).split("/")[0],
        "protein_position": L.clean(t.get("protein_start")),
    }, ""


def main():
    argparse.ArgumentParser().parse_args()
    rows = list(csv.DictReader(open(SRC, encoding="utf-8-sig"), delimiter="\t"))

    say("# DDD / DECIPHER de novo SMC1A cohort -- ingest report")
    say()
    say(f"Source: `{os.path.relpath(SRC, ROOT)}`")
    say(f"Rows read: **{len(rows)}**")
    say()

    v1, v2, v3, v4, v5, v6, v7 = [], [], [], [], [], [], []
    obs = []

    for r in rows:
        pid = L.clean(r["patient_id"])
        pos = int(L.clean(r["start"]))
        ref = L.clean(r["ref_allele"]).upper()
        alt = L.clean(r["alt_allele"]).upper()

        # ---- V1 gene and variant class ------------------------------
        if L.clean(r["gene_name"]) != "SMC1A":
            v1.append((pid, L.clean(r["gene_name"])))
        if L.clean(r["variant_class"]) != "sequence_variant":
            v1.append((pid, L.clean(r["variant_class"])))

        # ---- V2 assembly --------------------------------------------
        if L.clean(r["assembly"]) != "GRCh38":
            v2.append((pid, L.clean(r["assembly"])))

        # ---- V3 reference allele matches GRCh38 ---------------------
        obs_ref = REFSEQ.get(pos, pos + len(ref) - 1).upper()
        if obs_ref != ref:
            v3.append((pid, pos, ref, obs_ref))

        # ---- V4 VEP on MANE -----------------------------------------
        v, verr = vep_region(pos, ref, alt)
        if verr:
            v4.append((pid, pos, verr))

        # ---- V5 inheritance -----------------------------------------
        inh = L.clean(r["inheritance"])
        if "de novo" not in inh.lower():
            v5.append((pid, inh))

        # ---- V6/V7 disease attribution ------------------------------
        assoc = L.clean(r["associated_condition"])
        if assoc and assoc not in ASSOC_MAP:
            v6.append((pid, assoc))
        sex = SEX_MAP.get(L.clean(r["gender"]).lower(), "")
        if assoc == "DEE85" and sex != "female":
            v7.append((pid, sex))

        d_reported, d_harmonised, d_basis = ASSOC_MAP.get(
            assoc, ("", "",
                    "none -- the DECIPHER clinical notes do not name either "
                    "syndrome, so no disease is attributed"))

        cls, cls_note = PATH_MAP.get(L.clean(r["pathogenicity"]).lower(),
                                     ("not_classified", ""))
        eligible = cls in ("P", "LP")

        g_end = pos + len(ref) - 1
        design, core, amp = TARGETONS.assign(pos, g_end)

        hpo_ids = [x.strip() for x in L.clean(r["child_hpo"]).split(";") if x.strip()]
        hpo_names = [x.strip() for x in L.clean(r["child_terms"]).split(";") if x.strip()]

        obs.append({
            "source": "DDD_DECIPHER",
            "citation": "DDD / DECIPHER (unpublished cohort extract)",
            "pmid": "",
            "individual_key": f"DDD:{pid}",
            "patient_id_source": pid,
            "decipher_variant_id": L.clean(r["variant_id"]),

            "chrom": L.CHROM,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "variant_key": L.variant_key(L.CHROM, pos, ref, alt),
            "coordinate_basis": ("GRCh38 plus-strand VCF-style allele supplied "
                                 "directly by DECIPHER; not derived from a cDNA "
                                 "description"),
            # DECIPHER supplies no HGVS description (`user_hgvs` is empty for
            # every row), so the MANE c. derived here from the supplied
            # coordinate is offered as the resolution input. This routes the
            # cohort through the same guards (G1-G6) and the same downstream
            # join as every other source, rather than bypassing them.
            "hgvs_c_input": v.get("hgvs_c_mane", ""),
            "hgvs_p_input": "",
            "source_transcript": L.clean(r["user_transcript"]),
            "source_transcript_is_mane": (
                L.clean(r["user_transcript"]) == L.TRANSCRIPT_NOVER),
            **{k: v.get(k, "") for k in (
                "consequence_terms", "consequence_class", "hgvs_c_mane",
                "hgvs_p_mane", "vep_exon", "vep_intron", "protein_position")},
            "source_consequence": L.clean(r["vep_consequence"]),
            "in_design_window": bool(design),
            "exon": v.get("vep_exon", ""),
            "intron": v.get("vep_intron", ""),
            "phenotype_detail": "; ".join(dict.fromkeys(hpo_names)),
            "is_independent_proband": "assumed",
            "family_key": f"DDD:{pid}",
            "source_record_id": L.clean(r["variant_id"]),
            "targetons_design": "|".join(design),
            "targetons_exon_core": "|".join(core),
            "targetons_amplicon": "|".join(amp),

            "pathogenicity_reported_raw": L.clean(r["pathogenicity"]),
            "pathogenicity_reported": cls,
            "classification_method_raw": "DECIPHER curator classification",
            "classification_note": cls_note,

            "disease_reported_raw": d_reported,
            "disease_harmonised": d_harmonised,
            "attribution_basis": d_basis,
            "associated_condition_raw": assoc,
            "associated_condition_notes": L.clean(
                r.get("associated_condition_notes")),
            "eligible_for_curated_set": eligible,

            "sex": sex,
            "zygosity_reported": L.clean(r["genotype"]),
            "inheritance": "de_novo",
            "inheritance_detail": inh,
            "hpo_ids": ";".join(dict.fromkeys(hpo_ids)),
            "hpo_names": ";".join(dict.fromkeys(hpo_names)),
            "n_hpo_terms": len(dict.fromkeys(hpo_ids)),
            "has_clinical_detail": bool(hpo_ids),
            "proband_independence_basis": (
                "each row is a distinct DECIPHER patient_id and the identifiers "
                "are unique, so the 15 rows are 15 individuals; relatedness "
                "between DECIPHER records is not exposed, so independence is "
                "assumed rather than established"),
        })

    # ------------------------------------------------------------------
    say("## Verification")
    say()
    for name, fails, desc in (
        ("V1", v1, "every row is SMC1A and a sequence variant"),
        ("V2", v2, "assembly is GRCh38 (no liftover required)"),
        ("V3", v3, "supplied reference allele matches GRCh38"),
        ("V4", v4, "VEP returns a consequence on the MANE transcript"),
        ("V5", v5, "inheritance is de novo"),
        ("V6", v6, "every associated_condition value is mappable"),
        ("V7", v7, "any DEE85 attribution is female"),
    ):
        say(f"- **{name}** {desc}: "
            + ("PASS" if not fails else f"**{len(fails)} FAILURE(S)** {fails[:5]}"))
    say()

    nonmane = [o for o in obs if not o["source_transcript_is_mane"]]
    say(f"- rows reported against a non-MANE transcript: **{len(nonmane)}**"
        + (f" ({', '.join(o['patient_id_source'] + ' -> ' + o['source_transcript'] for o in nonmane)})"
           if nonmane else ""))
    say("  These affect the source's own annotation only. The join key is genomic "
        "and the MANE c./p. description is derived here, so the variant's "
        "identity is unaffected.")
    say()

    say("## Disease attribution")
    say()
    say("| associated_condition | variants | eligible for curated set |")
    say("|---|---|---|")
    for k, n in Counter(o["associated_condition_raw"] or "(none)"
                        for o in obs).most_common():
        el = len([o for o in obs if (o["associated_condition_raw"] or "(none)") == k
                  and o["eligible_for_curated_set"]])
        say(f"| `{k}` | {n} | {el} |")
    say()
    say("An empty `associated_condition` is a real value: the DECIPHER notes do "
        "not name either syndrome. Those rows keep their classification and "
        "remain in the all-variants list, but carry no disease and cannot enter "
        "either analysis arm.")
    say()

    say("## Classification")
    say()
    for k, n in Counter(o["pathogenicity_reported"] for o in obs).most_common():
        say(f"- `{k}`: {n}")
    say()

    say("## Sex and zygosity")
    say()
    for k, n in Counter(o["sex"] or "(not supplied)" for o in obs).most_common():
        say(f"- {k}: {n}")
    for k, n in Counter(o["zygosity_reported"] for o in obs).most_common():
        say(f"- zygosity `{k}`: {n}")
    say()

    say("## Targeton coverage")
    say()
    for k, n in Counter(o["targetons_design"] or "(none)" for o in obs).most_common():
        say(f"- `{k}`: {n}")
    say()

    cols = list(dict.fromkeys(k for o in obs for k in o))
    out = os.path.join(INTERIM, "ddd_observations.tsv")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for o in obs:
            # escape tabs and newlines -- an unescaped newline splits a logical
            # row across physical lines (the defect fixed project-wide earlier)
            fh.write("\t".join(str(o.get(c, "")).replace("\t", " ")
                                .replace("\n", " ").replace("\r", " ")
                                for c in cols) + "\n")
    say(f"Wrote `data/interim/ddd_observations.tsv` "
        f"({len(obs)} rows, {len(cols)} columns).")
    print(f"HTTP cache: {http.n_hits} hits, {http.n_misses} misses.")

    os.makedirs(DOCS, exist_ok=True)
    with open(os.path.join(DOCS, "ddd_ingest.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(report) + "\n")
    print("Wrote docs/verification/ddd_ingest.md")

    if any((v1, v2, v3, v4, v5, v6, v7)):
        sys.exit(1)


if __name__ == "__main__":
    main()
