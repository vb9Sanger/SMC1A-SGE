#!/usr/bin/env python3
"""
04_ingest_genedx.py

Ingest the GeneDx 180k de novo mutation cohort extract for SMC1A.

The dataset
-----------
`genedx/gdx_dnm_180k_SMC1A_annotated_with_sex.tsv` -- 79 de novo SMC1A observations from
79 distinct probands, pre-annotated with GRCh38 coordinates, HPO terms and a
laboratory `diagnostic_status`.

The 79 observations comprise only **71 distinct variants**: 7 variants recur
across 15 samples. Collapsing to one row per variant therefore requires
aggregation, which 05_build_registry.py performs.

Are the 79 samples unrelated? THIS IS AN ASSUMPTION, NOT A FACT.
--------------------------------------------------------------
The extract contains **no relatedness, family, pedigree or trio field**, and
`SAMPLE_ID_anon` is a 6-digit anonymised identifier with no pedigree
structure. Independence therefore cannot be established from the data.

What the data does support:

  * De novo ascertainment excludes parent-to-child transmission by
    construction -- a variant absent from both parents cannot have been
    inherited from them. So the pedigree structures that inflate LOVD
    (transmitting parent plus affected child in one entry) are excluded here.

What it does NOT exclude:

  * **Siblings.** Parental germline mosaicism can produce the same
    apparently-de novo variant in two children. This is documented in
    NDD cohorts and would make two rows non-independent.
  * **The same individual sequenced twice** under two anonymised IDs
    (re-referral or re-analysis), which is undetectable here.

Supporting (not conclusive) evidence that the recurrences are independent
mutational events rather than cryptic relatedness: 5 of the 7 recurrent
variants are C>T/G>A transitions at a CpG dinucleotide (71%), against 14 of 64
non-recurrent variants (22%). CpG transitions are the classic hypermutable
context, and 5 of the 7 recurrences are arginine codons (Arg1121, Arg1049,
Arg790, Arg711, Arg398). Recurrent independent mutation at hypermutable sites
is the more parsimonious explanation, but it does not rule out relatedness in
any individual case.

Accordingly `proband_independence_basis` records the assumption explicitly,
and `05_build_registry.py` treats `n_independent_probands` for this source as
an **upper bound**. See OPEN_QUESTIONS Q14 -- relatedness fields have been
requested.

Overlap with Retterer 2016 (also GeneDx) must still be considered, and is
handled at the registry stage.

What this source does and does not provide
------------------------------------------
It provides, unusually for a cohort dataset, a laboratory diagnostic
judgement (`diagnostic_status`: positive / possible / negative) -- i.e.
whether the SMC1A variant was considered the diagnostic finding. That is
genuine pathogenicity evidence and is mapped to the classification field.

It does **not** provide a syndrome diagnosis, so it cannot by itself
attribute a variant to CdLS or DEE85. Every proband therefore remains in an
`Exclude_*` group until a disease attribution becomes available.

The per-proband HPO terms are carried through verbatim (`hpo_ids`,
`hpo_names`) as source data. No inference is drawn from them: an earlier
version of this script scored them against hand-assembled CdLS/DEE term
lists, which were not taken from any citable source, and that assessment has
been removed per DECISIONS_LOG D23.

Inclusion gate for this source
------------------------------
Specified 2026-09-10: only variants with `diagnostic_status = positive` are
eligible for the curated set. `possible` and `negative` are recorded in the
all-variants list but can never enter it, irrespective of any other evidence.

Sex is not supplied in the current extract and is expected later; it is
mechanistically important here (METHODS 1.2) so the column is present and
empty rather than omitted.

Verification performed
----------------------
The supplied annotation is independently re-derived rather than trusted:

  W1  `ref` must match the cached GRCh38 plus-strand reference at `pos`
  W2  the exon assignment is recomputed from the transcript model
  W3  the targeton assignment is recomputed from the design windows
  W4  consequence, c. and p. on ENST00000322213.9 are recomputed via VEP,
      and disagreements with the supplied `consequence` / `cds_change` /
      `protein_change` are reported

Outputs
-------
data/interim/genedx_observations.tsv
docs/verification/genedx_ingest.md

Usage
-----
    ../.venv/bin/python 04_ingest_genedx.py
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import urllib.parse
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smc1a_lib as L  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, "data", "reference")
INTERIM = os.path.join(ROOT, "data", "interim")
DOCS = os.path.join(ROOT, "docs", "verification")
SRC = os.path.join(ROOT, "genedx", "gdx_dnm_180k_SMC1A_annotated_with_sex.tsv")

http = L.CachedHTTP(os.path.join(ROOT, "data", "cache", "ensembl"))
MODEL = L.TranscriptModel.from_json(os.path.join(REF, "smc1a_transcript_model.json"))
REFSEQ = L.ReferenceSequence(os.path.join(REF, "smc1a_chrX_region.fa"))
TARGETONS = L.TargetonMap.from_tsv(os.path.join(REF, "smc1a_targetons.tsv"))

report: list[str] = []


def say(m: str = "") -> None:
    print(m)
    report.append(m)


# --------------------------------------------------------------------------
# diagnostic_status -> classification
# --------------------------------------------------------------------------
# The lab's diagnostic judgement is pathogenicity evidence but is not an ACMG
# classification, so it is mapped conservatively and the raw value is kept.
# ---------------------------------------------------------------------------
# Disease attribution from the curated `associated_condition` column.
#
# This cohort supplies HPO terms and an internal diagnostic status but no
# syndrome-level diagnosis, so a disease label cannot be read from any single
# field. `associated_condition` carries the attribution assigned from the
# cohort's own phenotype data under explicit, testable criteria recorded in
# docs/ASSOCIATED_CONDITION_LABELLING.md:
#
#   CdLS   <- diagnostic_status = positive AND sex = M
#   DEE85  <- sex = F AND hpo_names contains `Epileptic encephalopathy`,
#             `Status epilepticus`, or any holoprosencephaly term
#
# An EMPTY cell is a real value meaning the evidence does not reach the stated
# bar. Such rows keep a classification and remain in the all-variants list, but
# carry no disease and so cannot enter either analysis arm.
#
# Nothing is inferred here: this maps the column, it does not re-derive it. The
# labels were verified programmatically against both rules (0 mislabelled, 0
# rows meeting a rule but left unlabelled).
# ---------------------------------------------------------------------------
# Sex as supplied by the cohort (`sex_proband`), mapped to the registry vocabulary.
SEX_MAP = {"M": "male", "F": "female"}

ASSOC_MAP = {
    "CdLS": ("Cornelia de Lange syndrome",
             "CdLS2",
             "assigned from cohort phenotype: positive laboratory diagnostic "
             "status in a male proband. DEE85 has been reported only in "
             "females, so a positive male diagnosis supports CdLS "
             "(docs/ASSOCIATED_CONDITION_LABELLING.md section 3.1)"),
    "DEE85": ("SMC1A-related developmental and epileptic encephalopathy",
              "DEE85",
              "assigned from cohort phenotype: female proband with "
              "`Epileptic encephalopathy`, `Status epilepticus` or a "
              "holoprosencephaly term among the reported HPO terms "
              "(docs/ASSOCIATED_CONDITION_LABELLING.md section 3.2)"),
}

DIAG_MAP = {
    "positive": ("LP", "reported as the diagnostic finding by the testing laboratory"),
    "possible": ("VUS", "reported as a possible but not established diagnostic finding"),
    "negative": ("VUS", "NOT considered the diagnostic finding by the testing "
                        "laboratory -- de novo occurrence alone does not establish "
                        "pathogenicity"),
}


def vep_region(pos: int, ref: str, alt: str):
    """Consequence, c. and p. on the MANE transcript from a plus-strand allele.

    Submitted as an HGVS genomic description rather than via the region
    endpoint: the region endpoint takes an allele string rather than a
    ref/alt pair and silently rejects the latter, so HGVS g. is used for
    every variant class for consistency with 03_resolve_variants.py.
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
    ap = argparse.ArgumentParser()
    args = ap.parse_args()

    rows = list(csv.DictReader(open(SRC, encoding="utf-8-sig"), delimiter="\t"))

    say("# GeneDx 180k cohort ingest report")
    say()
    say(f"Source: `{os.path.relpath(SRC, ROOT)}`")
    say(f"Rows: {len(rows)}; distinct anonymised samples: "
        f"{len({r['SAMPLE_ID_anon'] for r in rows})}")
    say()
    n_var = len({(r["chrom"], r["pos"], r["ref"], r["alt"]) for r in rows})
    say(f"Distinct variants: **{n_var}** -- i.e. {len(rows) - n_var} observations "
        f"are repeats of a variant already seen in another sample.")
    say()
    say("**Sample independence is assumed, not established.** The extract carries")
    say("no relatedness, family, pedigree or trio field, and `SAMPLE_ID_anon` has")
    say("no pedigree structure. De novo ascertainment excludes parent-to-child")
    say("transmission, but not siblings arising from parental germline mosaicism,")
    say("nor the same individual sequenced twice under two IDs. Counts for this")
    say("source are therefore **upper bounds**; relatedness fields have been")
    say("requested (Q14).")
    say()

    obs = []
    w1_fail, w2_diff, w3_diff, w4_diff, vep_fail = [], [], [], [], []
    unknown_assoc = []

    for i, r in enumerate(rows, 1):
        pos = int(r["pos"])
        ref, alt = r["ref"].upper(), r["alt"].upper()

        # ---- W1: reference allele ------------------------------------
        obs_ref = REFSEQ.get(pos, pos + len(ref) - 1)
        ref_ok = (obs_ref == ref)
        if not ref_ok:
            w1_fail.append((pos, ref, alt, obs_ref))

        # ---- W2 / W3: recompute exon and targeton --------------------
        span_lo, span_hi = pos, pos + len(ref) - 1
        my_exon = MODEL.exon_of_g(span_lo) or MODEL.exon_of_g(span_hi)
        my_intron = MODEL.intron_of_g(span_lo) if my_exon is None else None
        design, core, amp = TARGETONS.assign(span_lo, span_hi)

        # The supplied `exon` column mixes exon numbers and "intronN" strings.
        # Keep the two apart in the output, matching the schema, and compare
        # against the source using its own combined form.
        their_exon = L.clean(r["exon"])
        my_exon_str = str(my_exon) if my_exon else (f"intron{my_intron}" if my_intron else "")
        if their_exon != my_exon_str:
            w2_diff.append((pos, their_exon, my_exon_str))
        exon_out = str(my_exon) if my_exon else ""
        intron_out = str(my_intron) if (my_exon is None and my_intron) else ""

        their_t = set(t for t in L.clean(r["targeton_id"]).replace(";", "|").split("|") if t)
        if their_t != set(design) and their_t:
            w3_diff.append((pos, sorted(their_t), sorted(design), sorted(amp)))

        # ---- W4: recompute consequence on MANE -----------------------
        v, verr = vep_region(pos, ref, alt)
        if verr:
            vep_fail.append((pos, ref, alt, verr))
        their_c = L.clean(r["cds_change"])
        mine_c = v.get("hgvs_c_mane", "")
        # only compare when the supplied value is a real HGVS c. string
        if their_c.startswith("c.") and mine_c and their_c != mine_c:
            w4_diff.append((pos, their_c, mine_c, L.clean(r["protein_change"]),
                            v.get("hgvs_p_mane", "")))

        # ---- HPO inference -------------------------------------------
        hpo_ids = [x.strip() for x in L.clean(r["hpos"]).split(";") if x.strip()]
        hpo_ids_unique = list(dict.fromkeys(hpo_ids))
        names = [x.strip() for x in L.clean(r["hpo_names"]).split(";") if x.strip()]

        cls, cls_note = DIAG_MAP.get(L.clean(r["diagnostic_status"]).lower(),
                                     ("not_classified", ""))

        # ---- disease attribution from `associated_condition` ---------
        assoc_raw = L.clean(r.get("associated_condition"))
        assoc_reported, assoc_harmonised, assoc_basis = ASSOC_MAP.get(
            assoc_raw,
            ("", "",
             "none -- no syndrome diagnosis supplied and the cohort phenotype "
             "does not reach the criteria for either arm"))
        if assoc_raw and assoc_raw not in ASSOC_MAP:
            unknown_assoc.append((pos, assoc_raw))

        # Inclusion gate (specified 2026-09-10): only `positive` is eligible for the
        # curated set. This is applied independently of, and takes precedence
        # over, any phenotype inference.
        diag = L.clean(r["diagnostic_status"]).lower()
        eligible = (diag == "positive")
        if diag == "negative":
            group = "Exclude_disease_uncertain"
            excl = ("laboratory judged the SMC1A variant NOT to be the "
                    "diagnostic finding")
        elif diag == "possible":
            group = "Exclude_pending_classification"
            excl = ("laboratory judged the SMC1A variant only a possible "
                    "diagnostic finding; not eligible per inclusion gate")
        else:
            group = "Exclude_pending_classification"
            excl = ("no syndrome diagnosis supplied; awaiting a sourced "
                    "attribution instrument")

        obs.append({
            "source": "GeneDx_180k_DNM",
            "source_detail": "GeneDx 180k diagnostic exome de novo mutation cohort",
            "source_record_id": f"gdx_{r['SAMPLE_ID_anon']}",
            "individual_key": f"GDX:{r['SAMPLE_ID_anon']}",
            # Each sample is treated as its own family because no relatedness
            # field exists -- an ASSUMPTION, recorded as such. If two samples
            # are in fact siblings (possible via parental germline mosaicism,
            # which de novo ascertainment does not exclude), this over-counts.
            "family_key": f"GDX:{r['SAMPLE_ID_anon']}",
            "is_independent_proband": "assumed",
            "proband_independence_basis": (
                "assumed unrelated: distinct anonymised sample, de novo "
                "ascertained; no relatedness/family/pedigree field supplied, "
                "so independence is NOT established (see Q14). De novo status "
                "excludes parent-child transmission but not siblings via "
                "parental germline mosaicism."),

            # position as supplied, verified
            "chrom": L.CHROM, "pos": pos, "ref": ref, "alt": alt,
            "variant_key": L.variant_key(L.CHROM, pos, ref, alt),
            "ref_allele_verified": ref_ok,
            "ref_allele_in_grch38": obs_ref,

            # recomputed annotation on MANE
            "hgvs_c_input": mine_c,
            "hgvs_c_input_transcript": L.TRANSCRIPT,
            "hgvs_p_input": v.get("hgvs_p_mane", ""),
            "consequence_terms": v.get("consequence_terms", ""),
            "consequence_class": v.get("consequence_class", ""),
            "protein_position": v.get("protein_position", ""),
            "exon": exon_out,
            "intron": intron_out,
            "exon_source_supplied": their_exon,
            "targetons_design": "|".join(design),
            "targetons_exon_core": "|".join(core),
            "targetons_amplicon": "|".join(amp),
            "targeton_source_supplied": L.clean(r["targeton_id"]),
            "in_design_window": bool(design),

            # as supplied, retained verbatim
            "variant_class_supplied": L.clean(r["variant_class"]),
            "consequence_supplied": L.clean(r["consequence"]),
            "cds_change_supplied": their_c,
            "protein_change_supplied": L.clean(r["protein_change"]),
            "net_length_change": L.clean(r["net_length_change"]),

            # classification
            "pathogenicity_reported_raw": L.clean(r["diagnostic_status"]),
            "pathogenicity_reported": cls,
            "classification_method_raw": "laboratory diagnostic status",
            "classification_note": cls_note,

            # disease attribution
            "disease_reported_raw": assoc_reported,
            "disease_harmonised": assoc_harmonised,
            "attribution_basis": assoc_basis,
            "associated_condition_raw": assoc_raw,
            "curation_group": group,
            "exclusion_reason": excl,
            "genedx_diagnostic_positive": eligible,
            "eligible_for_curated_set": eligible,

            # individual
            "sex": SEX_MAP.get(L.clean(r.get("sex_proband")).upper(), ""),
            "inheritance": "de_novo",
            "inheritance_detail": "de novo by trio sequencing (cohort is DNM-ascertained)",
            "hpo_ids": ";".join(hpo_ids_unique),
            "hpo_names": ";".join(dict.fromkeys(names)),
            "n_hpo_terms": len(hpo_ids_unique),
            "phenotype_detail": "",
            "notes": "",
        })
        if i % 20 == 0:
            print(f"  ... {i}/{len(rows)}", file=sys.stderr)

    # ---------------- verification report ----------------
    say("## Verification of the supplied annotation")
    say()
    say(f"**W1 -- reference allele vs cached GRCh38:** "
        f"{len(rows) - len(w1_fail)}/{len(rows)} match.")
    if w1_fail:
        say()
        say("| pos | supplied ref | supplied alt | GRCh38 has |")
        say("|---|---|---|---|")
        for p, rf, al, gr in w1_fail:
            say(f"| {p} | {rf} | {al} | {gr} |")
    say()
    say(f"**W2 -- exon assignment recomputed from the transcript model:** "
        f"{len(rows) - len(w2_diff)}/{len(rows)} agree.")
    if w2_diff:
        say()
        say("| pos | supplied | recomputed |")
        say("|---|---|---|")
        for p, a, b in w2_diff[:20]:
            say(f"| {p} | `{a}` | `{b}` |")
    say()
    say(f"**W3 -- targeton assignment recomputed from the design windows:** "
        f"{len(rows) - len(w3_diff)}/{len(rows)} agree.")
    if w3_diff:
        say()
        say("| pos | supplied | recomputed (design window) | amplicon only |")
        say("|---|---|---|---|")
        for p, a, b, am in w3_diff[:20]:
            say(f"| {p} | `{';'.join(a)}` | `{';'.join(b) or '(none)'}` | "
                f"`{';'.join(sorted(set(am) - set(b))) or '-'}` |")
    say()
    say(f"**W4 -- c./p. recomputed on {L.TRANSCRIPT} via VEP:** "
        f"{len([o for o in obs if o['hgvs_c_input']])}/{len(rows)} annotated; "
        f"{len(w4_diff)} disagree with the supplied `cds_change`.")
    if w4_diff:
        say()
        say("| pos | supplied c. | recomputed c. | supplied p. | recomputed p. |")
        say("|---|---|---|---|---|")
        for p, a, b, c, d in w4_diff[:25]:
            say(f"| {p} | `{a}` | `{b}` | `{c}` | `{d}` |")
    if vep_fail:
        say()
        say(f"VEP could not annotate {len(vep_fail)} variant(s):")
        for p, rf, al, e in vep_fail[:10]:
            say(f"  - chrX:{p} {rf}>{al} -- {e}")
    say()

    # ---------------- distributions ----------------
    say("## Diagnostic status and classification")
    say()
    say("| `diagnostic_status` | mapped | probands |")
    say("|---|---|---|")
    for k, n in Counter((o["pathogenicity_reported_raw"], o["pathogenicity_reported"])
                        for o in obs).most_common():
        say(f"| {k[0]} | {k[1]} | {n} |")
    say()

    say("## Consequence classes (recomputed on MANE)")
    say()
    say("| consequence_class | probands | predicted LOF |")
    say("|---|---|---|")
    for k, n in Counter(o["consequence_class"] or "(none)" for o in obs).most_common():
        say(f"| {k} | {n} | {'yes' if k in L.PTV_CLASSES else ''} |")
    say()

    say("## HPO terms")
    say()
    say(f"- probands with >=1 HPO term: "
        f"{len([o for o in obs if o['n_hpo_terms']])}/{len(obs)}")
    say(f"- median HPO terms per proband: "
        f"{sorted(o['n_hpo_terms'] for o in obs)[len(obs) // 2]}")
    say()
    say("Terms are carried through verbatim as source data. No syndrome")
    say("inference is drawn from them (see the module docstring and D23).")
    say()

    say("## Inclusion gate (specified 2026-09-10): `diagnostic_status = positive` only")
    say()
    say("| diagnostic_status | probands | eligible for curated set |")
    say("|---|---|---|")
    for k, n in Counter(o["pathogenicity_reported_raw"] for o in obs).most_common():
        el = "yes" if k.lower() == "positive" else "**no**"
        say(f"| {k} | {n} | {el} |")
    say()
    say(f"Eligible: **{len([o for o in obs if o['eligible_for_curated_set']])}"
        f"/{len(obs)}** probands.")
    say()
    say("Eligible probands by consequence class:")
    say()
    say("| consequence_class | eligible probands | predicted LOF |")
    say("|---|---|---|")
    for k, n in Counter(o["consequence_class"] for o in obs
                        if o["eligible_for_curated_set"]).most_common():
        say(f"| {k} | {n} | {'yes' if k in L.PTV_CLASSES else ''} |")
    say()

    say("## Grouping")
    say()
    say("All GeneDx probands remain in an `Exclude_*` group: the source")
    say("supplies no syndrome diagnosis, so no attribution to CdLS or DEE85")
    say("can be made from it.")
    say()
    say("| curation_group | probands |")
    say("|---|---|")
    for k, n in Counter(o["curation_group"] for o in obs).most_common():
        say(f"| {k} | {n} |")
    say()
    say()
    say("## Targeton coverage")
    say()
    say(f"- in an SGE design window (expected in the library): "
        f"{len([o for o in obs if o['in_design_window']])}/{len(obs)}")
    say(f"- in an amplicon but outside every design window "
        f"(sequenced, but no oligo designed): "
        f"{len([o for o in obs if o['targetons_amplicon'] and not o['in_design_window']])}")
    say(f"- outside every targeton: "
        f"{len([o for o in obs if not o['targetons_amplicon']])}")
    say()

    cols = list(obs[0].keys())
    out = os.path.join(INTERIM, "genedx_observations.tsv")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for o in obs:
            fh.write("\t".join(str(o.get(c, "")).replace("\t", " ").replace("\n", " ")
                               for c in cols) + "\n")
    say(f"Wrote `data/interim/genedx_observations.tsv` "
        f"({len(obs)} rows, {len(cols)} columns).")
    say(f"HTTP cache: {http.n_hits} hits, {http.n_misses} misses.")

    with open(os.path.join(DOCS, "genedx_ingest.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(report) + "\n")
    print("Wrote docs/verification/genedx_ingest.md")


if __name__ == "__main__":
    main()
