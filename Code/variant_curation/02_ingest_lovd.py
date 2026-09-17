#!/usr/bin/env python3
"""
02_ingest_lovd.py

Ingest the LOVD Global Variome full-gene export for SMC1A and flatten its
relational structure into one row per *observation*, where an observation is
one (variant x individual/report) pair.

Why observations rather than variants
-------------------------------------
LOVD holds 179 variant entries linked to 92 individuals via screenings. The
same variant can appear as several entries (different submitters, different
individuals, or a submitter-level "CLASSIFICATION record" with no individual
attached). Collapsing straight to one row per variant would destroy exactly the
information needed to (a) count independent probands for hotspot analysis and
(b) avoid over-weighting a variant that appears three times because three
submitters deposited the same family. The observation table is therefore the
primary evidence ledger; 05_build_registry.py collapses it to one row per
variant while retaining the provenance.

Relational path walked
----------------------
    Variants_On_Genome  (genomic description, classification, dbSNP, ClinVar)
      |-- Variants_On_Transcripts (c./p./exon on NM_006306.2)
      |-- Screenings_To_Variants -> Screenings -> Individuals
                                                   |-- Individuals_To_Diseases -> Diseases
                                                   |-- Phenotypes

Input
-----
data/raw/lovd/LOVD_SMC1A_full_download_<date>.txt
    Obtained from the public endpoint
    https://databases.lovd.nl/shared/download/all/gene/SMC1A
    (no login required; this is the same content as the "Full data view" web
    page but complete and machine-readable, so the HTML/PDF page captures in
    variant_curation/sources/ are retained only as provenance snapshots).

Output
------
data/interim/lovd_observations.tsv
data/interim/lovd_variant_entries.tsv     one row per LOVD variant entry
docs/verification/lovd_ingest.md

Usage
-----
    ../.venv/bin/python 02_ingest_lovd.py [--export <path>]
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smc1a_lib as L  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INTERIM = os.path.join(ROOT, "data", "interim")
DOCS = os.path.join(ROOT, "docs", "verification")
os.makedirs(INTERIM, exist_ok=True)
os.makedirs(DOCS, exist_ok=True)

report: list[str] = []


def say(m: str = "") -> None:
    print(m)
    report.append(m)


# --------------------------------------------------------------------------
# Disease harmonisation
# --------------------------------------------------------------------------
# LOVD disease codes are free-ish text maintained per-submitter. They are
# mapped onto the four analysis categories here, deliberately conservatively:
# anything that does not unambiguously indicate CdLS or a DEE/EIEE-85
# presentation becomes 'unspecified_NDD' or 'other', which routes the variant
# to an Exclude_* group rather than into the high-confidence set.
DISEASE_MAP = {
    "CDLS": ("CdLS", "Cornelia de Lange syndrome (unspecified type)"),
    "CDLS2": ("CdLS", "Cornelia de Lange syndrome type 2 (SMC1A, MIM 300590)"),
    "EIEE85": ("DEE85", "Developmental and epileptic encephalopathy 85 (MIM 301044)"),
    "EIEE": ("DEE_unspecified", "Early infantile epileptic encephalopathy, type not stated"),
    "DEE": ("DEE_unspecified", "Developmental and epileptic encephalopathy, type not stated"),
    "epilepsy": ("epilepsy_only", "Epilepsy, no encephalopathy stated"),
    "seizures": ("epilepsy_only", "Seizures, no encephalopathy stated"),
    "ID": ("unspecified_NDD", "Intellectual disability, syndrome not stated"),
    "NDD": ("unspecified_NDD", "Neurodevelopmental disorder, syndrome not stated"),
    "?": ("unknown", "Disease not recorded"),
    "CF": ("other", "Cystic fibrosis -- SMC1A variant is incidental to the ascertained condition"),
    "MYOP": ("other", "Myopathy -- SMC1A variant is incidental to the ascertained condition"),
    "MD": ("other", "Muscular dystrophy -- SMC1A variant is incidental to the ascertained condition"),
    "CHTE": ("other", "Central hypothyroidism/testicular enlargement (IGSF1-region) -- incidental"),
}

# LOVD 'effectid' encodes the submitter's and curator's view of the effect on
# function as a two-digit code (reported/concluded). Decoded for transparency.
EFFECT_CODE = {
    "0": "not classified", "1": "no effect", "3": "probably no effect",
    "5": "effect unknown", "7": "probably affects function",
    "9": "affects function", "": "not classified",
}


def decode_effect(eid: str) -> str:
    eid = L.clean(eid)
    if len(eid) == 2:
        return (f"reported={EFFECT_CODE.get(eid[0], eid[0])}; "
                f"concluded={EFFECT_CODE.get(eid[1], eid[1])}")
    return EFFECT_CODE.get(eid, eid)


# LOVD clinical classification strings -> the 5-tier ACMG vocabulary.
CLASS_MAP = {
    "pathogenic": "P",
    "likely pathogenic": "LP",
    "pathogenic (dominant)": "P",
    "likely pathogenic (dominant)": "LP",
    "pathogenic (recessive)": "P",
    "vus": "VUS",
    "variant of unknown significance": "VUS",
    "likely benign": "LB",
    "benign": "B",
    "": "not_classified",
    "-": "not_classified",
}


def map_class(s: str) -> str:
    s = L.clean(s).lower()
    if s in CLASS_MAP:
        return CLASS_MAP[s]
    if "likely patho" in s:
        return "LP"
    if "patho" in s:
        return "P"
    if "likely benign" in s:
        return "LB"
    if "benign" in s:
        return "B"
    if "unknown" in s or "uncertain" in s or s == "vus":
        return "VUS"
    return "not_classified"


# --------------------------------------------------------------------------
# Genetic origin / inheritance harmonisation
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# LOVD `allele` codebook
# --------------------------------------------------------------------------
# The web "Full data view" renders this field as a label ('Maternal
# (confirmed)'), but the download stores the underlying numeric code. The
# mapping below was confirmed by matching the authenticated HTML capture
# against the export record-by-record (see 02b_verify_lovd_capture.py):
# code 21 -> 'Maternal (confirmed)', 20 -> 'Maternal (inferred)',
# 1 -> 'Parent #1', 0 -> 'Unknown'.
#
# This matters: without decoding, every parental-allele assignment in LOVD is
# invisible and all 14 maternally-inherited SMC1A entries would be harmonised
# to 'unknown'.
ALLELE_CODE = {
    "0":  "Unknown",
    "1":  "Parent #1",
    "2":  "Parent #2",
    "3":  "Both alleles",
    "10": "Paternal (inferred)",
    "11": "Paternal (confirmed)",
    "20": "Maternal (inferred)",
    "21": "Maternal (confirmed)",
    "":   "",
}


def decode_allele(code: str) -> str:
    code = L.clean(code)
    if code in ALLELE_CODE:
        return ALLELE_CODE[code]
    return f"unrecognised_allele_code_{code}"


# --------------------------------------------------------------------------
# Family grouping
# --------------------------------------------------------------------------
# Recurrence is evidence of a mutational hotspot only when the observations are
# independent, so siblings must not be counted as separate probands. LOVD
# records kinship in three places of differing reliability, and all three are
# used differently:
#
#   1. `panelid` -- a hard, structured link between two individual records.
#      Authoritative. Only 3 individuals carry one.
#   2. `Individual/Individual_ID` -- a submitter label. Several submitters use
#      a "<PMID>.<family><P|S>" convention (e.g. 16604071.7P / 16604071.7S =
#      proband and sibling of family 7 from PMID 16604071), which is
#      deterministic and parseable.
#   3. `Individual/Remarks` -- free text. NOT parsed for family assignment.
#      Most remarks describe a de novo singleton ("2-generation family, 1
#      affected, unaffected non-carrier parents") and imply independence, not
#      relatedness; only some imply an affected relative. Deciding which is a
#      per-record judgement, so these are FLAGGED for manual adjudication
#      rather than resolved by regex.
#
# Rules 1 and 2 are combined with a union-find so that transitive links are
# captured.

LABEL_FAMILY = re.compile(r"^(\d{6,8})\.(\d+)([PS])$")

# Remark wording that indicates an affected relative or an entry standing for
# more than one person -- i.e. where proband independence is questionable.
KINSHIP_SUSPECT = re.compile(
    r"sister|brother|sibling|twin|cousin|affected mother|affected father|"
    r"mother/child|affected parent|affected parant|inherited from|"
    r"\d+\s*affected(?!\s*,?\s*unaffected non-carrier)|"
    r"[23]-generation family,\s*[2-9]", re.I)


# --------------------------------------------------------------------------
# Reference / publication parsing
# --------------------------------------------------------------------------
# LOVD stores citations as structured markup in `VariantOnGenome/Reference` and
# `Individual/Reference`, e.g.
#     {PMID:Musio 2006:16604071}, {OMIM300040:0001}
#     {PMID:DDDS 2015:25533962}, {DOI:DDDS 2015:10.1038/nature14135}
#     PMID: 30158690
#     Redeker (unpublished)
# The PMID is extracted to its own field so every LOVD observation can be
# traced to its primary publication. This also drives the literature step: a
# paper already cited by LOVD needs verifying, not rediscovering.
#
# The OMIM allelic-variant ids (OMIM300040:000N) are kept too -- OMIM 300040 is
# the SMC1A gene entry, and these identify the specific allelic variants OMIM
# has curated, which is an independent attribution source.

REF_PMID = re.compile(r"\{PMID:([^:}]*):(\d{6,9})\}")
REF_PMID_BARE = re.compile(r"\bPMID:?\s*(\d{6,9})\b")
REF_DOI = re.compile(r"\{DOI:[^:}]*:([^}]+)\}")
REF_OMIM = re.compile(r"\{(OMIM\d+:[\dA-Za-z.]+)\}")
REF_UNPUB = re.compile(r"unpublished|personal communication", re.I)


def parse_references(*fields):
    """Extract PMIDs, citation labels, DOIs and OMIM allelic-variant ids."""
    text = " , ".join(L.clean(f) for f in fields if L.clean(f))
    pmids, cites = [], []
    for label, pmid in REF_PMID.findall(text):
        if pmid not in pmids:
            pmids.append(pmid)
            cites.append(L.clean(label))
    for pmid in REF_PMID_BARE.findall(text):
        if pmid not in pmids:
            pmids.append(pmid)
            cites.append("")
    return {
        "pmids": ";".join(pmids),
        "citations": ";".join(c for c in cites if c),
        "dois": ";".join(dict.fromkeys(REF_DOI.findall(text))),
        "omim_allelic_variants": ";".join(dict.fromkeys(REF_OMIM.findall(text))),
        "reference_is_unpublished": str(bool(REF_UNPUB.search(text))),
        "reference_raw": text,
    }


class _UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def build_families(individuals: dict):
    """Assign each LOVD individual to a family.

    Returns (family_of, basis_of): individual id -> family key, and the rule
    that assigned it.
    """
    uf = _UnionFind()
    basis = {}

    for iid in individuals:
        uf.find(iid)

    # rule 1: explicit panel links
    for iid, ind in individuals.items():
        pid = L.clean(ind.get("panelid"))
        if pid and pid in individuals:
            uf.union(pid, iid)
            basis[iid] = f"lovd_panelid link to {pid}"
            basis.setdefault(pid, f"lovd_panelid link with {iid}")

    # rule 2: "<PMID>.<family><P|S>" submitter labels
    by_label_family = defaultdict(list)
    for iid, ind in individuals.items():
        m = LABEL_FAMILY.match(L.clean(ind.get("Individual/Individual_ID")))
        if m:
            by_label_family[(m.group(1), m.group(2))].append(iid)
    for (study, fam), members in by_label_family.items():
        if len(members) > 1:
            for m_ in members[1:]:
                uf.union(members[0], m_)
            for m_ in members:
                basis[m_] = (f"submitter label convention "
                             f"<PMID {study}>.<family {fam}><P|S>")

    family_of = {}
    for iid in individuals:
        root = uf.find(iid)
        lab = L.clean(individuals[root].get("Individual/Individual_ID"))
        m = LABEL_FAMILY.match(lab)
        if m:
            family_of[iid] = f"LOVD_fam:PMID{m.group(1)}.fam{m.group(2)}"
        else:
            family_of[iid] = f"LOVD_fam:{root}"
    return family_of, basis


def harmonise_inheritance(genetic_origin: str, allele: str, segregation: str):
    """Return (inheritance, de_novo_confidence).

    LOVD records inheritance in two places that frequently disagree:
    `VariantOnGenome/Genetic_origin` (e.g. 'De novo', 'Germline') and the
    variant `allele` field (e.g. 'Maternal (confirmed)'). Both are kept in the
    output; this function derives a single harmonised value and flags the
    disagreements, because de novo status is the main driver of confidence in
    pathogenicity for a dominant NDD gene.
    """
    go = L.clean(genetic_origin).lower()
    al = decode_allele(allele).lower()

    if "parent #" in al:
        # LOVD's 'Parent #1'/'Parent #2' deliberately does not say which
        # parent, so it establishes inheritance without establishing the
        # transmitting parent. For an X-linked gene that distinction matters,
        # so it is kept separate from the maternal/paternal categories.
        return "inherited_unspecified_parent", f"LOVD allele={decode_allele(allele)}"
    if "both alleles" in al:
        return "unknown", (f"LOVD allele={decode_allele(allele)} "
                           f"(homozygous/both-allele assignment; implausible for a "
                           f"dominant X-linked variant -- flagged for review)")
    if al.startswith("unrecognised"):
        return "unknown", al

    if "de novo" in go:
        # 'De novo' in Genetic_origin but a parental allele assignment is a
        # genuine internal contradiction in the record.
        if "maternal" in al or "paternal" in al:
            return "conflicting", "conflict: origin=de novo but allele=" + L.clean(allele)
        return "de_novo", "stated de novo (parental testing not itemised in LOVD)"
    if "maternal" in al:
        return ("inherited_maternal",
                "confirmed" if "confirmed" in al else "inferred")
    if "paternal" in al:
        return ("inherited_paternal",
                "confirmed" if "confirmed" in al else "inferred")
    if "inherited" in go:
        return "inherited_unspecified", ""
    if "germline" in go:
        # 'Germline' only states it is not somatic; it says nothing about
        # whether it was inherited or arose de novo.
        return "unknown", "origin recorded only as germline"
    if "somatic" in go:
        return "somatic", ""
    if "classification record" in go:
        return "not_applicable", "submitter classification record, no individual"
    return "unknown", ""


# --------------------------------------------------------------------------
# Parse the LOVD export
# --------------------------------------------------------------------------

def parse_export(path: str):
    sections: dict[str, list[dict]] = {}
    meta: list[str] = []
    cur = None
    hdr: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            m = re.match(r"^## (\w+) ## Do not remove", line)
            if m:
                cur, hdr = m.group(1), []
                sections[cur] = []
                continue
            if line.startswith("###") or line.startswith("## ") or line.startswith("# "):
                meta.append(line)
                continue
            if not line.strip() or cur is None:
                continue
            cells = [c.strip('"').replace('""', '"') for c in line.split("\t")]
            if not hdr:
                hdr = [c.strip("{}") for c in cells]
                continue
            sections[cur].append(dict(zip(hdr, cells)))
    return sections, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--export", default=None)
    args = ap.parse_args()

    path = args.export
    if path is None:
        cands = sorted(glob.glob(os.path.join(ROOT, "data", "raw", "lovd",
                                              "LOVD_SMC1A_full_download_*.txt")))
        if not cands:
            sys.exit("FATAL: no LOVD export found in data/raw/lovd/")
        path = cands[-1]

    S, meta = parse_export(path)

    say("# LOVD ingest report")
    say()
    say(f"Source export: `{os.path.relpath(path, ROOT)}`")
    for m in meta:
        if "LOVD-version" in m or "Filter" in m:
            say(f"    {m}")
    say()
    say("Section row counts as parsed:")
    say()
    say("| section | rows |")
    say("|---|---|")
    for k, v in S.items():
        say(f"| {k} | {len(v)} |")
    say()

    # ---- index the relational tables ----
    vog = {r["id"]: r for r in S["Variants_On_Genome"]}
    vot = {r["id"]: r for r in S["Variants_On_Transcripts"]}
    diseases = {r["id"]: r for r in S["Diseases"]}
    individuals = {r["id"]: r for r in S["Individuals"]}
    screenings = {r["id"]: r for r in S["Screenings"]}

    scr_by_variant = defaultdict(list)
    for r in S["Screenings_To_Variants"]:
        scr_by_variant[r["variantid"]].append(r["screeningid"])

    dis_by_individual = defaultdict(list)
    for r in S["Individuals_To_Diseases"]:
        dis_by_individual[r["individualid"]].append(r["diseaseid"])

    phen_by_individual = defaultdict(list)
    for r in S["Phenotypes"]:
        phen_by_individual[r["individualid"]].append(r)

    # ---- family grouping ----------------------------------------------
    # LOVD 'panelid' links members of a submitted panel/family; 'panel_size'>1
    # marks an entry that represents several individuals at once. Both are used
    # to decide whether an observation counts as an independent proband.
    say("## Family and panel structure")
    say()
    panels = defaultdict(list)
    for iid, ind in individuals.items():
        pid = L.clean(ind.get("panelid"))
        if pid:
            panels[pid].append(iid)
    multi = {k: v for k, v in panels.items() if len(v) > 1}
    say(f"- individuals with a parent panel: {sum(len(v) for v in panels.values())}")
    say(f"- panels containing >1 individual: {len(multi)}")
    big = [i for i, ind in individuals.items()
           if L.clean(ind.get("panel_size")) not in ("", "1")]
    say(f"- entries with panel_size > 1 (an entry standing for several people): {len(big)}")
    for i in big:
        say(f"    individual {i}: panel_size={individuals[i]['panel_size']}, "
            f"remarks={L.clean(individuals[i].get('Individual/Remarks')) or '-'}")
    say()

    # ---- family assignment --------------------------------------------
    family_of, family_basis = build_families(individuals)
    fam_sizes = Counter(family_of.values())
    multi_fams = {f: n for f, n in fam_sizes.items() if n > 1}
    say(f"- families containing >1 LOVD individual (after union-find over "
        f"panelid + submitter labels): {len(multi_fams)}")
    for f, n in sorted(multi_fams.items()):
        members = sorted(i for i, fk in family_of.items() if fk == f)
        say(f"    {f}: {n} individuals -- "
            + "; ".join(f"{i} ({L.clean(individuals[i].get('Individual/Individual_ID'))})"
                        for i in members))
        say(f"        basis: {family_basis.get(members[0], 'n/a')}")
    say()

    # ---- build observations -------------------------------------------
    obs = []
    n_no_individual = 0
    for vid, g in vog.items():
        t = vot.get(vid, {})
        screen_ids = scr_by_variant.get(vid, [])
        ind_ids = []
        for sid in screen_ids:
            s = screenings.get(sid)
            if s and L.clean(s.get("individualid")):
                ind_ids.append(s["individualid"])
        if not ind_ids:
            ind_ids = [None]
            n_no_individual += 1

        for iid in ind_ids:
            ind = individuals.get(iid, {}) if iid else {}
            dis_ids = dis_by_individual.get(iid, []) if iid else []
            dis_syms = [L.clean(diseases[d]["symbol"]) for d in dis_ids if d in diseases]
            dis_names = [L.clean(diseases[d]["name"]) for d in dis_ids if d in diseases]
            harm = [DISEASE_MAP.get(s, ("other", s))[0] for s in dis_syms]

            phens = phen_by_individual.get(iid, []) if iid else []
            phen_bits = []
            for p in phens:
                for k, v in p.items():
                    if k.startswith("Phenotype/") and not L.blank(v):
                        phen_bits.append(f"{k.split('/', 1)[1]}={L.clean(v)}")

            inh, dnconf = harmonise_inheritance(
                g.get("VariantOnGenome/Genetic_origin"),
                g.get("allele"),
                g.get("VariantOnGenome/Segregation"))

            sids = [sid for sid in screen_ids]
            techs = L.join_unique(
                [L.clean(screenings[s].get("Screening/Technique")) for s in sids
                 if s in screenings])
            tmpl = L.join_unique(
                [L.clean(screenings[s].get("Screening/Template")) for s in sids
                 if s in screenings])

            obs.append({
                # --- provenance -------------------------------------
                "source": "LOVD",
                "source_detail": "Global Variome shared LOVD (databases.lovd.nl)",
                "source_record_id": vid,
                "lovd_dbid": L.clean(g.get("VariantOnGenome/DBID")),
                "lovd_individual_id": L.clean(iid) if iid else "",
                "lovd_panel_id": L.clean(ind.get("panelid")),
                "lovd_panel_size": L.clean(ind.get("panel_size")),
                "individual_key": f"LOVD:{iid}" if iid else "",
                "family_key": family_of.get(iid, ""),
                "family_key_basis": (family_basis.get(iid)
                                     or ("singleton: no panel link or family "
                                         "label" if iid else "")),
                "kinship_note": (L.clean(ind.get("Individual/Remarks"))
                                 if iid and KINSHIP_SUSPECT.search(
                                     L.clean(ind.get("Individual/Remarks")) or "")
                                 else ""),
                "kinship_needs_adjudication": str(bool(
                    iid and (KINSHIP_SUSPECT.search(
                        L.clean(ind.get("Individual/Remarks")) or "")
                        or L.clean(ind.get("panel_size")) not in ("", "1")))),
                "is_independent_proband": "yes" if iid else "no_individual",
                "lovd_owner": L.clean(g.get("owned_by")),
                **parse_references(g.get("VariantOnGenome/Reference"),
                                   ind.get("Individual/Reference")),
                "individual_label_raw": L.clean(ind.get("Individual/Individual_ID")),
                "data_availability": L.clean(ind.get("Individual/Data_av")),

                # --- variant description as submitted ---------------
                "hgvs_c_input": L.clean(t.get("VariantOnTranscript/DNA")),
                "hgvs_c_input_transcript": "NM_006306.2",
                "hgvs_p_input": L.clean(t.get("VariantOnTranscript/Protein")),
                "hgvs_r_input": L.clean(t.get("VariantOnTranscript/RNA")),
                "exon_input": L.clean(t.get("VariantOnTranscript/Exon")),
                "hgvs_g_hg19_input": L.clean(g.get("VariantOnGenome/DNA")),
                "hgvs_g_hg38_input": L.clean(g.get("VariantOnGenome/DNA/hg38")),
                "lovd_pos_g_start_hg19": L.clean(g.get("position_g_start")),
                "lovd_pos_g_end_hg19": L.clean(g.get("position_g_end")),
                "lovd_variant_type": L.clean(g.get("type")),
                "published_as": L.clean(g.get("VariantOnGenome/Published_as")),
                "iscn": L.clean(g.get("VariantOnGenome/ISCN")),

                # --- classification ----------------------------------
                "pathogenicity_reported_raw": L.clean(
                    g.get("VariantOnGenome/ClinicalClassification")),
                "pathogenicity_reported": map_class(
                    g.get("VariantOnGenome/ClinicalClassification")),
                "classification_method_raw": L.clean(
                    g.get("VariantOnGenome/ClinicalClassification/Method")),
                "lovd_effect_code": L.clean(g.get("effectid")),
                "lovd_effect_decoded": decode_effect(g.get("effectid")),

                # --- disease -----------------------------------------
                "disease_reported_raw": L.join_unique(dis_syms),
                "disease_reported_name": L.join_unique(dis_names),
                "disease_harmonised": L.join_unique(harm),

                # --- individual --------------------------------------
                "sex": {"M": "M", "F": "F"}.get(
                    L.clean(ind.get("Individual/Gender")), ""),
                "sex_raw": L.clean(ind.get("Individual/Gender")),
                "country": L.clean(ind.get("Individual/Origin/Geographic")),
                "population": L.clean(ind.get("Individual/Origin/Population")),
                "consanguinity": L.clean(ind.get("Individual/Consanguinity")),
                "age_of_death": L.clean(ind.get("Individual/Age_of_death")),
                "individual_remarks": L.clean(ind.get("Individual/Remarks")),

                # --- inheritance -------------------------------------
                "inheritance": inh,
                "inheritance_detail": dnconf,
                "genetic_origin_raw": L.clean(g.get("VariantOnGenome/Genetic_origin")),
                "allele_raw": decode_allele(g.get("allele")),
                "allele_code_raw": L.clean(g.get("allele")),
                "segregation_raw": L.clean(g.get("VariantOnGenome/Segregation")),

                # --- phenotype / assays ------------------------------
                "phenotype_detail": "; ".join(phen_bits),
                "methylation_raw": L.clean(g.get("VariantOnGenome/Methylation")),
                "screening_technique": techs,
                "screening_template": tmpl,
                "variant_remarks": L.clean(g.get("VariantOnGenome/Remarks")),
                "frequency_raw": L.clean(g.get("VariantOnGenome/Frequency")),

                # --- cross-references --------------------------------
                "clinvar_id": L.clean(g.get("VariantOnGenome/ClinVar")),
                "dbsnp_id": L.clean(g.get("VariantOnGenome/dbSNP")),
            })

    say(f"## Observations built")
    say()
    say(f"- LOVD variant entries: {len(vog)}")
    say(f"- entries with no linked individual (submitter classification records): "
        f"{n_no_individual}")
    say(f"- observations (variant x individual/report): {len(obs)}")
    say(f"- distinct submitted c. descriptions: "
        f"{len({o['hgvs_c_input'] for o in obs if o['hgvs_c_input']})}")
    say()

    # ---- integrity checks ---------------------------------------------
    say("## Integrity checks")
    say()
    missing_c = [o for o in obs if not o["hgvs_c_input"]]
    say(f"- observations with no c. description: {len(missing_c)}")
    for o in missing_c:
        say(f"    {o['source_record_id']} dbid={o['lovd_dbid']} "
            f"g(hg19)={o['hgvs_g_hg19_input']} type={o['lovd_variant_type']} "
            f"iscn={o['iscn'] or '-'}")

    unparseable = []
    for o in obs:
        if o["hgvs_c_input"] and L.parse_hgvs_c(o["hgvs_c_input"]) is None:
            unparseable.append(o["hgvs_c_input"])
    say(f"- c. descriptions this pipeline cannot parse: {len(set(unparseable))}")
    for u in sorted(set(unparseable)):
        say(f"    {u}")

    # Variants that are not confined to SMC1A. The LOVD DB-ID names the gene
    # database that *hosts* the entry, not the gene the variant affects: the
    # RIBC1_* entries are ordinary point variants in SMC1A intron 1 (RIBC1 is
    # the neighbouring gene). Span is therefore the correct discriminator, and
    # the exclusion criterion "variants only in the gene SMC1A" is applied on
    # span, with the DB-ID reported alongside for transparency.
    say(f"- entries hosted under another gene's LOVD DB-ID: "
        f"{len([o for o in obs if o['lovd_dbid'] and not o['lovd_dbid'].startswith('SMC1A')])}")
    LARGE = 1000  # bp; above this an event cannot be SMC1A-only at this locus
    multigene = []
    for o in obs:
        try:
            span = abs(int(o["lovd_pos_g_end_hg19"]) - int(o["lovd_pos_g_start_hg19"])) + 1
        except (ValueError, TypeError):
            span = None
        o["genomic_span_bp"] = span if span is not None else ""
        if span is not None and span > LARGE:
            multigene.append((o, span))
    say(f"- entries spanning >{LARGE} bp (cannot be SMC1A-only; flagged "
        f"`not_smc1a_only`): {len(multigene)}")
    for o, span in multigene:
        say(f"    {o['source_record_id']} dbid={o['lovd_dbid']} "
            f"c={o['hgvs_c_input']} span={span:,} bp iscn={o['iscn'] or '-'}")
    for o in obs:
        o["not_smc1a_only"] = any(o is x for x, _ in multigene)
    say(f"- entries hosted under a non-SMC1A DB-ID but confined to the SMC1A "
        f"locus (retained as SMC1A variants):")
    for o in obs:
        if (o["lovd_dbid"] and not o["lovd_dbid"].startswith("SMC1A")
                and not o["not_smc1a_only"]):
            say(f"    {o['source_record_id']} dbid={o['lovd_dbid']} "
                f"c={o['hgvs_c_input']} span={o['genomic_span_bp']} bp")

    # ---- audit LOVD's own transcript conversions ----------------------
    # `Published_as` records the description as it appeared in the original
    # source. Where that source used a different transcript, LOVD converted it
    # to NM_006306.2 numbering. That conversion is audited here rather than
    # trusted, because a systematic isoform offset would shift a variant by
    # whole codons.
    #
    # NM_001281463.1 is a different SMC1A ISOFORM (verified in step 01): CDS
    # 3636 nt vs 3702, sharing MANE's 3' end but with a distinct N-terminus.
    # Conversion rule for its c. >= 43:  c_MANE = c_alt + 66  (p + 22 codons).
    # Its c.1-42 lies in MANE intron 1 and has no MANE CDS equivalent.
    ISOFORM_OFFSET = {"NM_001281463.1": (66, 43)}

    say()
    say("## Audit of LOVD's transcript conversions (`Published_as`)")
    say()
    pub_tx = defaultdict(int)
    for o in obs:
        m = re.search(r"\((NM_[\d.]+)\)", o["published_as"] or "")
        o["published_as_transcript"] = m.group(1) if m else ""
        if o["published_as_transcript"]:
            pub_tx[o["published_as_transcript"]] += 1
    say("| transcript cited by the original source | observations | c. numbering vs MANE |")
    say("|---|---|---|")
    for tx, n in sorted(pub_tx.items(), key=lambda kv: -kv[1]):
        note = ("**different isoform -- offset +66 nt / +22 codons**"
                if tx in ISOFORM_OFFSET else "identical (verified in step 01)")
        say(f"| {tx} | {n} | {note} |")
    say()

    n_ok = n_bad = n_na = 0
    mismatches = []
    for o in obs:
        tx = o["published_as_transcript"]
        o["alt_isoform_numbering"] = tx in ISOFORM_OFFSET
        o["conversion_audit"] = ""
        if tx not in ISOFORM_OFFSET:
            continue
        off, min_c = ISOFORM_OFFSET[tx]
        pm = re.search(r"c\.([-*]?\d+)", o["published_as"])
        pl = L.parse_hgvs_c(o["hgvs_c_input"])
        if not pm or not pl:
            n_na += 1
            continue
        try:
            c_alt = int(pm.group(1))
        except ValueError:
            n_na += 1
            continue
        if c_alt < min_c:
            # isoform-specific first exon: expected to land in MANE intron 1
            intronic = "+" in pl["pos1"] or "-" in pl["pos1"][1:]
            o["conversion_audit"] = (
                f"published c.{c_alt} on {tx} is in the isoform-specific first "
                f"exon; expected to map into MANE intron 1 -- LOVD gives "
                f"{o['hgvs_c_input']} ({'intronic, consistent' if intronic else 'NOT intronic, INCONSISTENT'})")
            if intronic:
                n_ok += 1
            else:
                n_bad += 1
                mismatches.append(o)
            continue
        expected = c_alt + off
        try:
            got = int(pl["pos1"])
        except ValueError:
            n_na += 1
            continue
        if got == expected:
            n_ok += 1
            o["conversion_audit"] = f"consistent: c.{c_alt}+{off} = c.{got}"
        else:
            n_bad += 1
            o["conversion_audit"] = (f"INCONSISTENT: published c.{c_alt} on {tx} "
                                     f"implies MANE c.{expected}, LOVD stores c.{got} "
                                     f"(difference {got - expected})")
            mismatches.append(o)

    say(f"- conversions consistent with the derived offset: {n_ok}")
    say(f"- conversions inconsistent: {n_bad}")
    say(f"- not auditable (no comparable c. position): {n_na}")
    if mismatches:
        say()
        say("| LOVD c. (MANE numbering) | Published_as | audit |")
        say("|---|---|---|")
        for o in mismatches[:25]:
            say(f"| `{o['hgvs_c_input']}` | `{o['published_as'][:46]}` | "
                f"{o['conversion_audit'][:130]} |")
    say()
    say("Observations on an alternative isoform are flagged "
        "`alt_isoform_numbering` so that any variant whose only source cites "
        "`NM_001281463.1` is converted rather than taken at face value.")

    say()
    say("## Publications cited by LOVD")
    say()
    say("Every observation is traced to its primary publication where LOVD")
    say("records one. These papers need **verifying** during literature")
    say("extraction, not rediscovering, and the PMIDs give an exact join to the")
    say("PDFs in `variant_curation/sources/`.")
    say()
    pmid_cite, pmid_n = {}, defaultdict(int)
    for o in obs:
        cs = o["citations"].split(";")
        for i, pm in enumerate([x for x in o["pmids"].split(";") if x]):
            pmid_n[pm] += 1
            if pm not in pmid_cite and i < len(cs) and cs[i]:
                pmid_cite[pm] = cs[i]
    say("| PMID | citation | observations | PDF in sources/ (match basis) |")
    say("|---|---|---|---|")
    src_dir = os.path.join(ROOT, "variant_curation", "sources")

    def fold(t: str) -> str:
        """Accent- and case-insensitive fold, so 'Pie' matches 'Pie\u0301'."""
        import unicodedata
        return "".join(c for c in unicodedata.normalize("NFKD", t.lower())
                       if not unicodedata.combining(c))

    pdfs = [fold(f) for f in os.listdir(src_dir)] if os.path.isdir(src_dir) else []
    missing = []
    for pm, n in sorted(pmid_n.items(), key=lambda kv: -kv[1]):
        cite = pmid_cite.get(pm, "")
        parts = cite.split()
        author = fold(parts[0].rstrip(",")) if parts else ""
        year = next((w for w in parts if w.isdigit() and len(w) == 4), "")
        # Filenames in variant_curation/sources/ are hand-typed and contain
        # transliteration and letter-transposition errors -- Tzschach 2015 is
        # filed as "Tszchach_2015.pdf" -- so a prefix match alone both misses
        # real files and, for short surnames, risks false positives. Four
        # tests are applied in order of decreasing confidence, and for very
        # short surnames the year is required as well.
        import difflib

        def find_pdf(au: str, yr: str):
            if len(au) < 3:
                return None
            need_year = len(au) <= 3          # e.g. "Liu", "Sun"
            for f in pdfs:
                if need_year and yr and yr not in f:
                    continue
                if au in f:
                    return "substring"
            for f in pdfs:
                if need_year and yr and yr not in f:
                    continue
                for tok in re.split(r"[^a-z0-9]+", f):
                    if len(tok) < 3:
                        continue
                    if sorted(tok) == sorted(au):
                        return "anagram"     # catches letter transposition
                    if difflib.SequenceMatcher(None, au, tok).ratio() >= 0.85:
                        return "fuzzy"
            if yr and any(yr in f and au[:3] in f for f in pdfs):
                return "year+prefix"
            return None

        basis = find_pdf(author, year)
        have = basis is not None
        if not have:
            missing.append((pm, cite, n))
        say(f"| {pm} | {cite or '-'} | {n} | "
            f"{basis if have else '**NOT FOUND**'} |")
    say()
    say(f"- distinct PMIDs cited: {len(pmid_n)}")
    say(f"- observations with a PMID: "
        f"{len([o for o in obs if o['pmids']])}/{len(obs)}")
    say(f"- observations citing unpublished data: "
        f"{len([o for o in obs if o['reference_is_unpublished'] == 'True'])}")
    say(f"- observations with an OMIM allelic-variant id: "
        f"{len([o for o in obs if o['omim_allelic_variants']])} "
        f"({L.join_unique([o['omim_allelic_variants'] for o in obs])})")
    say()
    if missing:
        say(f"**{len(missing)} publication(s) cited by LOVD have no matching PDF "
            f"in `variant_curation/sources/`.** Filename matching is approximate, "
            f"so verify before sourcing:")
        say()
        for pm, cite, n in missing:
            say(f"  - PMID {pm} -- {cite or '(no label)'} ({n} observation(s)) "
                f"-- https://pubmed.ncbi.nlm.nih.gov/{pm}/")
        say()

    say("## Family grouping outcome")
    say()
    say("| basis for family assignment | observations |")
    say("|---|---|")
    for k, n in Counter(o["family_key_basis"] for o in obs if o["family_key_basis"]).most_common():
        say(f"| {k} | {n} |")
    say()
    adj = [o for o in obs if o["kinship_needs_adjudication"] == "True"]
    say(f"### Records needing manual kinship adjudication ({len(adj)})")
    say()
    say("Free-text remarks are **not** parsed for family assignment. Most")
    say("describe a de novo singleton (\"1 affected, unaffected non-carrier")
    say("parents\") and imply independence; only some imply an affected relative")
    say("or an entry standing for several people. Those flagged below are listed")
    say("for a human decision rather than resolved by regex.")
    say()
    say("| individual | family_key | panel_size | remark |")
    say("|---|---|---|---|")
    seen = set()
    for o in sorted(adj, key=lambda x: x["lovd_individual_id"]):
        if o["lovd_individual_id"] in seen:
            continue
        seen.add(o["lovd_individual_id"])
        say(f"| {o['lovd_individual_id']} | `{o['family_key']}` | "
            f"{o['lovd_panel_size'] or '-'} | {o['kinship_note'] or o['individual_remarks'][:70]} |")
    say()

    say("## Disease code usage")
    say()
    counts = defaultdict(int)
    for o in obs:
        counts[(o["disease_reported_raw"] or "(none)", o["disease_harmonised"] or "(none)")] += 1
    say("| LOVD code | harmonised | observations |")
    say("|---|---|---|")
    for (raw, harm), n in sorted(counts.items(), key=lambda kv: -kv[1]):
        say(f"| {raw} | {harm} | {n} |")

    say()
    say("## Classification distribution")
    say()
    c2 = defaultdict(int)
    for o in obs:
        c2[(o["pathogenicity_reported_raw"] or "(none)", o["pathogenicity_reported"])] += 1
    say("| LOVD classification | mapped | observations |")
    say("|---|---|---|")
    for (raw, m), n in sorted(c2.items(), key=lambda kv: -kv[1]):
        say(f"| {raw} | {m} | {n} |")

    say()
    say("## Inheritance distribution")
    say()
    c3 = defaultdict(int)
    for o in obs:
        c3[o["inheritance"]] += 1
    say("| harmonised inheritance | observations |")
    say("|---|---|")
    for k, n in sorted(c3.items(), key=lambda kv: -kv[1]):
        say(f"| {k} | {n} |")

    say()
    say("## Records carrying methylation / episignature data")
    say()
    meth = [o for o in obs if o["methylation_raw"]]
    say(f"- LOVD `VariantOnGenome/Methylation` populated: {len(meth)}")
    dav = [o for o in obs if o["data_availability"]]
    say(f"- LOVD `Individual/Data_av` populated: {len(dav)} "
        f"({L.join_unique([o['data_availability'] for o in dav])})")
    say()
    say("LOVD carries no episignature results for SMC1A. Episignature evidence is")
    say("therefore sourced from the dedicated publications (Aref-Eshghi 2020,")
    say("Levy 2022) in step 04.")

    # ---- write ---------------------------------------------------------
    cols = list(obs[0].keys())
    out = os.path.join(INTERIM, "lovd_observations.tsv")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for o in obs:
            fh.write("\t".join(str(o.get(c, "")).replace("\t", " ").replace("\n", " ")
                               for c in cols) + "\n")
    say()
    say(f"Wrote `data/interim/lovd_observations.tsv` ({len(obs)} rows, {len(cols)} columns).")

    with open(os.path.join(DOCS, "lovd_ingest.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(report) + "\n")
    print("Wrote docs/verification/lovd_ingest.md")


if __name__ == "__main__":
    main()
