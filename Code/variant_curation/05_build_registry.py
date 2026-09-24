#!/usr/bin/env python3
"""
05_build_registry.py

Collapse the per-source observation ledgers into the three deliverable tables,
assign each variant to an analysis group, and apply the inclusion criteria.

Outputs
-------
data/processed/smc1a_observations.tsv       the evidence ledger: one row per
                                            (variant x source x individual)
data/processed/smc1a_variants_all.tsv       one row per distinct variant,
                                            EVERY variant encountered
data/processed/smc1a_variants_curated.tsv   the high-confidence subset
docs/verification/registry_build.md

Design
------
The observation ledger is primary and lossless. The two variant-level tables
are derived from it, so nothing is ever discarded: a variant that fails every
inclusion criterion still appears in `..._all.tsv` with the reason recorded.

Counting rules (METHODS 4.4)
----------------------------
Three counts are kept deliberately separate, because recurrence is evidence of
a hotspot only when the observations are independent:

  n_observations          rows in the ledger. Inflated by re-reporting and by
                          multiple submitters depositing the same case.
  n_individuals           distinct individuals named by any source.
  n_independent_probands  ONE PER FAMILY -- the hotspot-safe count.
  n_families              distinct families (equal to the above by
                          construction; both are emitted for clarity).
  n_siblings_collapsed    n_individuals - n_families, i.e. how much the
                          family grouping reduced the naive count.

Cross-source individual linkage
-------------------------------
LOVD individuals and GeneDx probands cannot be linked: the GeneDx samples are
anonymised with no external identifier, and LOVD's individual identifiers are
submitter-local. No fuzzy matching on phenotype or variant identity is
attempted -- inferring that two records are the same person from a shared
variant would be circular, and a wrong merge is as damaging as a wrong split.
The consequence is stated explicitly in the report: for a variant seen in both
sources, `n_independent_probands` is an UPPER bound, and the affected variants
are flagged `possible_cross_source_duplicate` so the ambiguity is visible
rather than hidden.

Group assignment
----------------
Deterministic and applied in a fixed order (see `assign_group`). Classification
consensus is resolved first, then disease attribution. Any internal conflict
produces an explicit `conflicting` value rather than being silently resolved by
source precedence.

Usage
-----
    ../.venv/bin/python 05_build_registry.py
"""

from __future__ import annotations

import csv
import datetime
import re
import glob
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smc1a_lib as L  # noqa: E402
import smc1a_schema as S  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, "data", "reference")
INTERIM = os.path.join(ROOT, "data", "interim")
PROC = os.path.join(ROOT, "data", "processed")
DOCS = os.path.join(ROOT, "docs", "verification")
os.makedirs(PROC, exist_ok=True)

# Needed to normalise allele pairs onto one key form, which requires reading
# the reference to the left of an indel.
REFSEQ = L.ReferenceSequence(os.path.join(REF, "smc1a_chrX_region.fa"))


# The two analysis arms. Named once so that code asking "is this an arm?" does
# not re-spell the pair at each site.
ARM_GROUPS = ("CdLS_pathogenic", "DEE85_pathogenic")


# --------------------------------------------------------------------------
# Manual attribution overrides
# --------------------------------------------------------------------------
# Variants where the automated gate reaches an arm that patient-level review
# does not support. Each entry names the variant, the group it is moved to, and
# the reason, and every one is listed in the build report so an override cannot
# become invisible. A key that matches no variant ABORTS the build, so an
# override cannot silently go stale after a re-key or a source change.
#
# Overrides may only ever REMOVE a variant from an analysis arm, never admit
# one: an automated gate that is too permissive can be corrected here, but a
# variant that fails the stated criteria must not be let in by hand.
ATTRIBUTION_OVERRIDES = {
    # c.1756C>T p.Arg586Trp -- D137
    "chrX:53405648:G:A": {
        "group": "Exclude_disease_uncertain",
        "exclusion_reason": "conflicting_disease_attribution",
        "reason": ("neither arm is supported at patient level: the DEE85 "
                   "attribution rests on holoprosencephaly alone with no "
                   "seizure term of any kind, in a record whose other terms "
                   "include CdLS-compatible dysmorphology, while the only CdLS "
                   "claim is cohort ascertainment in Huisman 2017, whose "
                   "non-Dutch participants were CdLS-suspected before testing "
                   "(D137)"),
    },
}


def apply_attribution_overrides(variants, say_fn):
    """Apply the manual overrides, abort on a stale key, and report all of them."""
    by_key = {v["variant_key"]: v for v in variants}
    stale = [k for k in ATTRIBUTION_OVERRIDES if k not in by_key]
    if stale:
        raise SystemExit(
            f"ABORT: attribution override(s) match no variant: {stale}. An "
            f"override that no longer applies must be removed or re-keyed "
            f"deliberately, not left to do nothing.")

    applied = []
    for key, o in ATTRIBUTION_OVERRIDES.items():
        v = by_key[key]
        if o["group"] in ARM_GROUPS:
            raise SystemExit(
                f"ABORT: override for {key} targets the analysis arm "
                f"'{o['group']}'. Overrides may only remove a variant from an "
                f"arm, never admit one.")
        applied.append((v["hgvs_c_mane"], v["hgvs_p_mane"],
                        v["curation_group"], o["group"]))
        v["curation_group"] = o["group"]
        v["curation_group_reason"] = "manual override: " + o["reason"]
        v["exclusion_reason"] = o["exclusion_reason"]
        v["include_high_confidence"] = "False"
        v["notes"] = L.join_unique(
            [v.get("notes"), "curation_group set by manual override (D137)"])

    say_fn(f"Manual attribution overrides applied: {len(applied)}. These are "
           f"variants whose automated group patient-level review does not "
           f"support; an override can only remove a variant from an arm.")
    say_fn("")
    say_fn("| variant | protein | gate assigned | overridden to |")
    say_fn("|---|---|---|---|")
    for hg, p, was, now in applied:
        say_fn(f"| `{hg}` | {p} | `{was}` | `{now}` |")
    say_fn("")


def _coords_from_key(key, res, first):
    """Split a canonical variant key back into its coordinate columns.

    Guarantees `variant_key == chrom:pos_grch38:ref:alt` on every delivered
    row. An `UNRESOLVED::` key has no coordinates to split, so those rows fall
    back to whatever the resolver or the observation carried, which is how
    quarantined rows retain a best-effort position.
    """
    if not key.startswith("UNRESOLVED::"):
        parts = key.split(":")
        if len(parts) == 4:
            return {"chrom": parts[0], "pos_grch38": parts[1],
                    "ref": parts[2], "alt": parts[3]}
    return {"chrom": L.clean(res.get("chrom")) or L.clean(first.get("chrom")),
            "pos_grch38": L.clean(res.get("pos")) or L.clean(first.get("pos")),
            "ref": L.clean(res.get("ref")) or L.clean(first.get("ref")),
            "alt": L.clean(res.get("alt")) or L.clean(first.get("alt"))}


def check_key_coordinate_agreement(variants):
    """Abort the build if any row's key disagrees with its own coordinates, or
    is not canonical.

    A fifth build-time guard, added for the same reason as the other four: the
    defect it catches produced a plausible-looking table rather than an error.
    Two representations of one duplication coexisted for 113 curated rows
    without anything complaining, and the duplicate pair they created was only
    found when the assay join could not match a variant that was demonstrably
    in the library.
    """
    bad_cols, bad_canon = [], []
    for v in variants:
        key = v["variant_key"]
        if key.startswith("UNRESOLVED::"):
            continue
        rebuilt = L.variant_key(v["chrom"], int(v["pos_grch38"]),
                                v["ref"], v["alt"])
        if rebuilt != key:
            bad_cols.append(f"{key} versus columns {rebuilt}")
        canon = L.canonical_variant_key(v["chrom"], int(v["pos_grch38"]),
                                        v["ref"], v["alt"], REFSEQ)
        if canon != key:
            bad_canon.append(f"{key} normalises to {canon}")
    if bad_cols or bad_canon:
        for m in (bad_cols + bad_canon)[:20]:
            print(f"  {m}")
        raise SystemExit(
            f"ABORT: {len(bad_cols)} row(s) whose variant_key disagrees with "
            f"their own coordinate columns and {len(bad_canon)} whose key is "
            f"not canonical. Variant identity is not well defined, so nothing "
            f"downstream would be trustworthy.")
    say(f"Key/coordinate agreement: all {len(variants)} rows carry a "
        f"canonical `variant_key` equal to `chrom:pos_grch38:ref:alt`.")
    say()

# Provenance stamps written onto every delivered row. Both are DERIVED, not
# hardcoded: they were previously literals and went stale, so the delivered
# files claimed a build date five days earlier than the actual build and a
# pipeline version 34 releases out of date -- predating most of the curation
# they described. A provenance stamp that can drift silently is worse than none,
# because it is read as authoritative.
def _pipeline_version():
    """Latest version heading in the changelog."""
    path = os.path.join(ROOT, "docs", "CHANGELOG.md")
    best = None
    try:
        for line in open(path, encoding="utf-8"):
            m = re.match(r"^##\s+(v\d+\.\d+)\s*$", line.strip())
            if m:
                key = tuple(int(x) for x in m.group(1)[1:].split("."))
                if best is None or key > best[0]:
                    best = (key, m.group(1))
    except OSError:
        pass
    return best[1] if best else "unversioned"


PIPELINE_VERSION = _pipeline_version()
CURATION_DATE = datetime.date.today().isoformat()

report: list[str] = []


def say(m: str = "") -> None:
    print(m)
    report.append(m)


def truthy(v) -> bool:
    """Coerce a flag to bool regardless of whether it is a Python bool or the
    string form read back from a TSV.

    Exists because comparing a computed `bool` against the string "True"
    silently evaluates false, which previously caused every variant to fail
    inclusion criterion 4.
    """
    if isinstance(v, bool):
        return v
    return str(v).strip().lower() in ("true", "1", "yes")


def read_tsv(path: str) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --------------------------------------------------------------------------
# Consensus helpers
# --------------------------------------------------------------------------

CLASS_RANK = {"P": 5, "LP": 4, "VUS": 3, "LB": 2, "B": 1, "not_classified": 0}
PATHOGENIC = {"P", "LP"}
BENIGN = {"B", "LB"}


def classification_consensus(classes: list[str]) -> tuple[str, str]:
    """Resolve several reported classifications into one, plus a note.

    A variant reported both pathogenic and benign is `conflicting`; that is a
    finding, not something to average away.
    """
    s = {c for c in classes if c and c != "not_classified"}
    if not s:
        return "not_classified", ""
    if s & PATHOGENIC and s & BENIGN:
        return "conflicting", (f"reported as both pathogenic and benign: "
                               f"{'/'.join(sorted(s))}")
    best = max(s, key=lambda c: CLASS_RANK.get(c, 0))
    note = f"reported as {'/'.join(sorted(s))}" if len(s) > 1 else ""
    return best, note


# Which harmonised disease values constitute a usable attribution.
# Canonical disease tokens. Sources were registered with two spellings for the
# CdLS arm -- `CdLS` (LOVD) and `CdLS2` (several literature registrations, the
# GeneDx cohort and the DDD cohort) -- because the vocabulary was never enforced
# in one place. Both are accepted here and normalised to `CdLS` for reporting.
#
# The silent-failure mode this caused is the reason for `KNOWN_DISEASE_TOKENS`
# below: an unrecognised token matched neither arm and fell through to
# "unresolved", which is indistinguishable in the output from a genuine
# inability to resolve a recorded disease. Sources whose attributions were
# affected therefore looked correctly excluded rather than mis-mapped.
CDLS_SET = {"CdLS", "CdLS2"}
DEE_SET = {"DEE85"}
DISEASE_CANON = {"CdLS2": "CdLS"}
# Values that indicate a disease was recorded but cannot be resolved to
# either condition.
UNRESOLVED = {"DEE_unspecified", "epilepsy_only", "unspecified_NDD",
              "unknown", "other", ""}


# Definite diseases that are neither CdLS nor DEE85. Such a disease must not be
# absorbed into either arm, nor filed as "uncertain" -- the disease would be
# certain, merely outside the contrast -- so it is routed to its own exclusion
# reason, and to `conflicting` where another source asserts an in-scope disease.
#
# The set is EMPTY. Holoprosencephaly was assessed against it and belongs to
# DEE85, not outside it: MIM 301044 is "developmental and epileptic
# encephalopathy 85 WITH OR WITHOUT MIDLINE BRAIN DEFECTS", holoprosencephaly is
# the canonical midline brain defect, and OMIM cites the reporting series as a
# DEE85 primary. The routing is retained because the situation is real in
# principle; with an empty set it excludes nothing.
OUT_OF_SCOPE_DISEASES: set[str] = set()


def disease_consensus(values: list[str]) -> tuple[str, str]:
    """Resolve harmonised disease values across observations."""
    s = set()
    for v in values:
        for tok in str(v).split("|"):
            tok = tok.strip()
            if tok:
                s.add(DISEASE_CANON.get(tok, tok))
    # Diseases that are neither CdLS nor DEE85 but ARE definite. These must not
    # be silently absorbed into either arm, and must not be filed as "uncertain"
    # either -- the disease is certain, it is simply out of scope. The ruling
    # (D105): keep the variant in the all list with the disease recorded
    # verbatim, and exclude it from BOTH analysis arms.
    oos = sorted(t for t in s if t in OUT_OF_SCOPE_DISEASES)
    if oos:
        also = sorted((s & CDLS_SET) | (s & DEE_SET))
        if also:
            return "conflicting", (
                f"attributed to {'/'.join(oos)} by one source and "
                f"{'/'.join(also)} by another -- out-of-scope disease in conflict "
                f"with an in-scope attribution; requires manual adjudication")
        return "out_of_scope", f"recorded as {'/'.join(oos)}, a definite disease outside the CdLS/DEE85 contrast"

    cd, de = bool(s & CDLS_SET), bool(s & DEE_SET)
    if cd and de:
        return "conflicting", ("attributed to CdLS by one source and DEE85 by "
                               "another -- requires manual adjudication")
    if cd:
        extra = s - CDLS_SET - UNRESOLVED
        return "CdLS", (f"also carries {'/'.join(sorted(extra))}" if extra else "")
    if de:
        extra = s - DEE_SET - UNRESOLVED
        return "DEE85", (f"also carries {'/'.join(sorted(extra))}" if extra else "")
    if s - {""}:
        return "unresolved", f"recorded as {'/'.join(sorted(s - {''}))}"
    return "none", "no disease recorded"



# Every disease token a source may emit. An unknown token is a registration
# error, not a variant with an unresolvable disease, and must fail loudly.
KNOWN_DISEASE_TOKENS = (CDLS_SET | DEE_SET | UNRESOLVED | OUT_OF_SCOPE_DISEASES
                        | set(DISEASE_CANON))


def suppress_redundant_database_attributions(observations, say_fn):
    """A database deposit of a publication we hold contributes no disease claim.

    LOVD records frequently ARE a publication, deposited: 69 of its 92
    individuals carry a submitter label and PMID naming a specific paper, and 58
    of those name a paper ingested separately here. Where that is so, the
    deposit is not a second, independent observation of the disease -- it is the
    same evidence arriving by a second route. Counting it as an independent
    source lets a deposit outvote, or manufacture a conflict with, the very
    paper it was copied from.

    The concrete failure: Yuan 2019 records per case whether CdLS was even a
    differential diagnosis, and for 13 of 14 it was NOT, so D69 correctly takes
    no CdLS attribution from that paper. Its LOVD deposit nonetheless labels
    every patient CDLS -- `Diagnosis/Definite=CDLS2` appears verbatim on 17
    individuals in one submission block, alongside identical boilerplate
    remarks, so it is a submission-level field rather than a per-patient
    adjudication. It is applied to protein-truncating variants too, which is the
    CdLS/DEE conflation this curation exists to avoid. The result was 12
    variants carrying a CdLS token that the source publication explicitly
    declines to make, four of them excluded from the DEE85 arm on the conflict
    that created.

    Only the **disease** claim is suppressed. Everything else the deposit
    carries -- classification, sex, inheritance, phenotype fields -- is kept: a
    deposit may legitimately record details the paper omitted, and those are not
    duplicative claims about which disorder the patient has.

    Identity is a separate matter and is NOT changed here: these deposits still
    count as distinct individuals from the published patient, so proband counts
    remain upper bounds (Q37, Q38).
    """
    held = {str(o.get("source", "")) for o in observations}
    suppressed = []
    skipped_no_determination = []
    for o in observations:
        if o.get("source") != "LOVD":
            continue
        tok = L.clean(o.get("disease_harmonised"))
        if not tok:
            continue
        cite = L.clean(o.get("citations"))
        if not cite:
            continue
        # The literature sources are named with underscores for spaces.
        src = cite.replace(" ", "_")
        if src not in held:
            continue
        paper_tokens = {L.clean(x.get("disease_harmonised"))
                        for x in observations if x.get("source") == src}
        paper_tokens = {t for t in paper_tokens if t}
        # The rule is conditional on the publication having made its own
        # determination. Where a held paper records no disease token at all it
        # has said nothing to govern with, and a deposit's claim may be the only
        # evidence rather than a duplicate of it -- suppressing it would then
        # destroy information instead of removing redundancy. Every publication
        # linked in the present data does make a determination, so this branch
        # is not currently reached; it exists so that the code enforces the rule
        # as stated rather than a broader one that happens to coincide with it.
        if not paper_tokens:
            skipped_no_determination.append((o.get("individual_key"), cite, tok))
            continue
        suppressed.append((o.get("individual_key"), cite, tok,
                           "|".join(sorted(paper_tokens)) or "(none)"))
        o["disease_harmonised"] = ""
        o["disease_attribution_suppressed"] = (
            f"LOVD deposit of {cite}, which is ingested separately and records "
            f"its own disease determination ({'|'.join(sorted(paper_tokens))}); "
            f"the publication governs and the deposit adds no independent claim")

    say_fn(f"Redundant database disease claims suppressed: **{len(suppressed)}** "
           f"LOVD observation(s) that are deposits of a publication ingested "
           f"separately. The publication governs; the deposit adds no "
           f"independent disease claim. Classification, sex, inheritance and "
           f"phenotype from the deposit are retained.")
    say_fn("")
    if skipped_no_determination:
        say_fn(f"**{len(skipped_no_determination)} deposit(s) NOT suppressed** "
               f"because the publication they came from records no disease "
               f"determination of its own, so the deposit's claim is not a "
               f"duplicate of anything:")
        for ind, cite, tok in skipped_no_determination:
            say_fn(f"- `{ind}` (deposit of {cite}) keeps its `{tok}` claim")
        say_fn("")
    if suppressed:
        by_paper = Counter(c for _, c, _, _ in suppressed)
        say_fn("| publication | deposits suppressed | deposit said | publication says |")
        say_fn("|---|---|---|---|")
        for paper in sorted(by_paper):
            rows = [s for s in suppressed if s[1] == paper]
            dep = "|".join(sorted({r[2] for r in rows}))
            pub = "|".join(sorted({r[3] for r in rows}))
            say_fn(f"| {paper} | {len(rows)} | `{dep}` | `{pub}` |")
        say_fn("")
    return suppressed


def check_disease_vocabulary(observations):
    """Fail if any source emits a disease token this module cannot interpret."""
    seen = {}
    for o in observations:
        for tok in str(o.get("disease_harmonised", "")).split("|"):
            tok = tok.strip()
            if tok and tok not in KNOWN_DISEASE_TOKENS:
                seen.setdefault(tok, set()).add(o.get("source", "?"))
    return seen

# Threshold separating small indels from structural events. The observed span
# distribution is bimodal with a wide gap -- small indels reach 15 bp
# (c.173_187del, a genuine CdLS in-frame deletion) and the next largest event
# is 142 bp (a whole-exon deletion) -- so the exact value is not sensitive.
# ---------------------------------------------------------------------------
# Per-source criterion coverage
#
# An inclusion criterion written against one source family can silently reject
# another that stores the same information in a different field. Such a failure
# is invisible in aggregate output, because it surfaces as a per-variant
# exclusion carrying a legitimate-looking reason rather than as an error.
#
# This reports, per source, how many of its variants satisfy each criterion. A
# source at 0/N on a criterion its records plainly satisfy is a field-mapping
# defect, not a curation result. It is reported rather than enforced, because a
# source can legitimately score zero -- a classification-only submitter really
# does have no clinical detail -- so the numbers are surfaced for inspection.
# ---------------------------------------------------------------------------
def report_source_criterion_coverage(variants, say_fn):
    """Report, per source, how many of its variants satisfy each criterion."""
    crit = {
        "clinical detail": lambda v: truthy(v.get("has_clinical_detail")),
        "disease resolved": lambda v: v.get("disease_attribution") in ("CdLS", "DEE85"),
        "P/LP": lambda v: v.get("pathogenicity_consensus") in ("P", "LP"),
        "in design window": lambda v: truthy(v.get("in_design_window")),
    }
    srcs = sorted({x for v in variants
                   for x in str(v.get("sources", "")).split("|") if x})
    say_fn("| source | variants | " + " | ".join(crit) + " |")
    say_fn("|---" * (len(crit) + 2) + "|")
    zeroes = []
    for src in srcs:
        sub = [v for v in variants
               if src in str(v.get("sources", "")).split("|")]
        cells = []
        for name, fn in crit.items():
            n = len([v for v in sub if fn(v)])
            cells.append(f"{n}/{len(sub)}")
            if n == 0 and len(sub) > 2:
                zeroes.append((src, name, len(sub)))
        say_fn(f"| `{src}` | {len(sub)} | " + " | ".join(cells) + " |")
    say_fn("")
    if zeroes:
        say_fn("**Sources scoring zero on a criterion** -- check whether the "
               "criterion reads the field this source populates, rather than "
               "assuming the records genuinely fail it:")
        for src, name, n in zeroes:
            say_fn(f"- `{src}`: 0 of {n} variants pass **{name}**")
    else:
        say_fn("No source scores zero on any criterion across more than two "
               "variants.")
    say_fn("")
    return zeroes


# ---------------------------------------------------------------------------
# Contradicted CdLS attributions
#
# A CdLS label can arrive from two very different kinds of evidence: a
# publication that clinically diagnosed the patient, or a bare gene-level
# disease tag on a database record. Where a source additionally states that
# CdLS was NOT in the differential diagnosis for its patient, a label resting
# only on the second kind is not supportable -- it is disease inherited from
# the gene, which is the failure mode that disqualified aggregate database
# attribution for this project (section 1 of METHODS).
#
# The test therefore has three parts, and all must hold:
#   1. the consensus disease is CdLS;
#   2. at least one observation records that CdLS was NOT the differential, or
#      that the CdLS presentation was absent;
#   3. NO publication-derived observation asserts CdLS -- i.e. every CdLS
#      assertion comes from a database tag.
#
# Part 3 is what stops this from over-firing. A per-patient "CdLS was not the
# differential" is a fact about ONE patient, not about the variant, so a variant
# with independent clinical CdLS reports keeps its attribution and only the
# contradiction is recorded.
# ---------------------------------------------------------------------------
DB_SOURCES = {"LOVD"}


def cdls_attribution_contradicted(v, observations_for_variant):
    """Is a CdLS consensus contradicted, with no publication supporting it?"""
    if v.get("disease_attribution") != "CdLS":
        return False, ""
    contradicting = [o for o in observations_for_variant
                     if L.clean(o.get("cdls_differential_reported")).lower() == "no"
                     or truthy(o.get("cdls_presentation_absent_in_a_source"))]
    if not contradicting:
        return False, ""
    pub_cdls = [o for o in observations_for_variant
                if L.clean(o.get("source")) not in DB_SOURCES
                and "CdLS" in str(o.get("disease_harmonised", ""))]
    if pub_cdls:
        srcs = sorted({L.clean(o.get("source")) for o in pub_cdls})
        return False, (f"a source records that CdLS was not the differential, but the "
                       f"CdLS attribution is independently supported by "
                       f"{', '.join(srcs)}, so it stands")
    db = sorted({L.clean(o.get("source")) for o in observations_for_variant
                 if L.clean(o.get("source")) in DB_SOURCES
                 and "CdLS" in str(o.get("disease_harmonised", ""))})
    cs = sorted({L.clean(o.get("source")) for o in contradicting})
    return True, (f"CdLS is asserted only by database tag ({', '.join(db) or 'none'}) "
                  f"while {', '.join(cs)} records that CdLS was NOT in the differential "
                  f"diagnosis for its patient; a gene-level tag cannot carry a disease "
                  f"attribution against a primary source's own statement")


# ---------------------------------------------------------------------------
# Attribution provenance per arm
#
# An arm assignment is only as good as the evidence that supplies its disease.
# Three weaknesses are worth surfacing at every build rather than found by
# occasional audit:
#
#   * database-only  -- the disease is asserted only by a database record and by
#                       no publication. Not automatically wrong: a named
#                       submitter recording `Diagnosis/Definite` is case-level
#                       evidence. But a bare gene-level tag is not, and the two
#                       are indistinguishable once harmonised to a token.
#   * overlay-only   -- the disease is asserted only by a source registered
#                       overlay_only, which by definition may enrich but not
#                       establish. This should always be zero.
#   * single-source  -- only one source asserts the disease, so nothing
#                       corroborates it.
#
# Reported, not enforced. The CdLS contradiction gate is the enforcing rule;
# this is the instrument that would have surfaced the need for it.
# ---------------------------------------------------------------------------
def report_attribution_provenance(variants, obs_by_variant, say_fn):
    """Per arm, how many variants have weakly-provenanced disease attributions."""
    rows = []
    for arm, tok in (("CdLS_pathogenic", "CdLS"), ("DEE85_pathogenic", "DEE85")):
        sub = [v for v in variants if v.get("curation_group") == arm]
        db_only = ov_only = single = 0
        detail = []
        for v in sub:
            obs = obs_by_variant.get(v["variant_key"], [])
            asserts = [o for o in obs
                       if tok in str(o.get("disease_harmonised", ""))]
            if not asserts:
                continue
            srcs = {L.clean(o.get("source")) for o in asserts}
            pub = {x for x in srcs if x not in DB_SOURCES}
            nonov = {L.clean(o.get("source")) for o in asserts
                     if not truthy(o.get("is_overlay_source"))}
            if not pub:
                db_only += 1
                detail.append((v["hgvs_c_mane"], "database-only", sorted(srcs)))
            if not nonov:
                ov_only += 1
                detail.append((v["hgvs_c_mane"], "overlay-only", sorted(srcs)))
            if len(srcs) == 1:
                single += 1
        rows.append((arm, len(sub), db_only, ov_only, single, detail))
    say_fn("| arm | variants | database-only | overlay-only | single-source |")
    say_fn("|---|---|---|---|---|")
    for arm, n, d, o, sg, _ in rows:
        say_fn(f"| `{arm}` | {n} | {d} | {o} | {sg} |")
    say_fn("")
    for arm, n, d, o, sg, detail in rows:
        if detail:
            say_fn(f"**`{arm}` -- attributions with no publication or no "
                   f"non-overlay support:**")
            for hg, kind, srcs in detail:
                say_fn(f"- `{hg}` ({kind}): {', '.join(s for s in srcs if s)}")
            say_fn("")
    if any(o for _, _, _, o, _, _ in rows):
        say_fn("**An overlay-only attribution is a defect** -- an overlay source may "
               "enrich a variant but may not establish its disease.")
        say_fn("")
    return rows


# Sources write protein descriptions with varying rigour: some omit the `p.`
# prefix, and one uses a bare position (`p1095`) as shorthand for "protein
# position 1095", which is not a protein description at all. Where several
# sources describe the same variant, prefer a well-formed HGVS value over a
# malformed one instead of taking whichever happens to come first; and reject a
# value that is not a protein description rather than carrying it into the
# delivered file, per METHODS 4A.20.
P_WELLFORMED = re.compile(r"^p\.[A-Z][a-z]{2}\d+")
P_PLAUSIBLE = re.compile(r"^p?\.?\(?[A-Z][a-z]{2}\d+")


def pick_protein(vep_value, observations):
    """Choose the best protein description for a variant.

    Preference order, and each step exists because of an observed failure:

    1. the VEP-derived description from the resolution step. It is computed from
       the coordinate against MANE, so it cannot inherit a source's
       transcription error.
    2. a well-formed HGVS value from a NON-overlay source.
    3. an interpretable value from a non-overlay source, with the `p.` prefix
       supplied.
    4. the same from an overlay source, last.

    A source-stated value must never outrank a primary's or VEP's: one overlay
    review renders a primary's `Glu205fs*` as `Glu2055*`, and letting source
    order decide put that impossible description into the protein column.

    Any value implying a codon beyond the 1233-residue protein is rejected
    outright, which is what catches that class of error.
    """
    def usable(c):
        c = L.clean(c)
        if not c:
            return None
        m = re.match(r"^p?\.?\(?([A-Z][a-z]{2})(\d+)", c)
        if not m:
            return None
        if int(m.group(2)) > L.PROTEIN_LENGTH:
            return None          # impossible codon -- a transcription error
        return c

    v = usable(vep_value)
    if v:
        return (v if v.startswith("p.") else "p." + v.lstrip("p").lstrip(".")), ""

    tiers = [[], []]             # non-overlay first, overlay second
    for o in observations:
        bucket = 1 if truthy(o.get("is_overlay_source")) else 0
        for key in ("hgvs_p_mane", "hgvs_p_input"):
            c = usable(o.get(key))
            if c:
                tiers[bucket].append((c, L.clean(o.get("source"))))

    for bucket in tiers:
        for c, src in bucket:                       # well-formed first
            if re.match(r"^p\.[A-Z][a-z]{2}\d+", c):
                return c, ""
        for c, src in bucket:                       # prefix missing only
            return "p." + c.lstrip("p").lstrip("."), (
                f"protein description normalised from {c!r} as given by {src}: "
                f"the source omitted the `p.` prefix")

    rejected = sorted({L.clean(o.get(k)) for o in observations
                       for k in ("hgvs_p_mane", "hgvs_p_input")
                       if L.clean(o.get(k))})
    if rejected:
        return "", ("no usable protein description; source value(s) "
                    + ", ".join(repr(x) for x in rejected)
                    + " are not valid HGVS protein descriptions for a "
                    f"{L.PROTEIN_LENGTH}-residue protein, so the column is left "
                    "empty rather than carrying an unusable value")
    return "", ""


# variant_key -> its contributing observations, populated in main().
OBS_BY_VARIANT: dict[str, list] = {}

# Global individual -> family map, populated in main() from all observations.
FAMILY_OF: dict[str, str] = {}

LARGE_EVENT_BP = 50


def is_large_event(v: dict) -> tuple[bool, str]:
    """Is this a structural event rather than a small indel?

    Specified 2026-09-14: such variants are kept in the all-variants list but
    excluded from the curated set, because the SGE library contains designed
    SNVs and small indels only and cannot represent them. They are also the
    variants whose breakpoints are reported at inconsistent precision, which
    is what made one exon-16 deletion appear as four rows (D53, Q20).
    """
    ref, alt = L.clean(v.get("ref")), L.clean(v.get("alt"))
    if ref and alt:
        span = abs(len(ref) - len(alt))
        if span > LARGE_EVENT_BP:
            return True, f"deletion/insertion of {span} bp (> {LARGE_EVENT_BP} bp)"
    c = L.clean(v.get("hgvs_c_mane"))
    if "-?" in c or "+?" in c:
        return True, "breakpoints undetermined in the source description"
    if L.clean(v.get("consequence_class")).startswith("cnv"):
        return True, "annotated as a copy-number event"
    return False, ""


def assign_group(v: dict) -> tuple[str, str, bool, str]:
    """Return (curation_group, reason, include_high_confidence, exclusion_reason).

    Applied in a fixed order so the outcome is reproducible.
    """
    cls = v["pathogenicity_consensus"]
    dis = v["disease_attribution"]

    # 1. events that are not confined to SMC1A cannot attribute a phenotype to
    #    SMC1A regardless of how well classified they are
    if truthy(v["not_smc1a_only"]):
        return ("Exclude_disease_uncertain",
                "multi-gene event: phenotype cannot be attributed to SMC1A alone",
                False, "not_smc1a_only")

    # 2. coordinates that could not be resolved on the MANE transcript
    if v["resolution_status"] == "quarantined":
        return ("Exclude_pending_classification",
                f"coordinate resolution failed ({v['guard_failed']})",
                False, "unresolved_coordinates")

    # 3. classification
    if cls == "conflicting":
        return ("Exclude_VUS",
                "conflicting classifications across sources", False,
                "conflicting_classification")
    if cls in BENIGN:
        return ("Exclude_Benign", f"classified {cls}", False, "benign")
    if cls == "VUS":
        return ("Exclude_VUS", "classified VUS", False, "vus")
    if cls == "not_classified":
        return ("Exclude_pending_classification",
                "no classification available", False, "not_classified")

    # 4. classification is P/LP -- now the disease attribution decides
    if dis == "conflicting":
        return ("Exclude_disease_uncertain",
                "P/LP but sources disagree on the condition", False,
                "conflicting_disease_attribution")
    # The two exclusion groups are distinguished by whether the ambiguity is
    # irreducible or merely not yet resolved -- the distinction specified for this study:
    #   * a disease WAS recorded but cannot be mapped to CdLS or DEE85
    #     (e.g. '?', 'ID', 'NDD', 'epilepsy', 'seizures', 'EIEE') -> the
    #     evidence is ambiguous as it stands -> Exclude_disease_uncertain
    #   * NO disease was recorded at all -> nothing to be ambiguous about, the
    #     attribution is simply outstanding and is expected to become
    #     available (chiefly the GeneDx probands, who have HPO terms but no
    #     syndrome label) -> Exclude_pending_classification
    if dis == "unresolved":
        # Distinguish "no source has said yet" from "a source considered the
        # diagnosis and declined it". Both leave the variant unattributable,
        # but only the first is upgradeable, and `Exclude_pending_classification`
        # promises upgradeability. A paper that records CdLS was NOT the
        # differential diagnosis has already looked (D140).
        declined = truthy(v.get("disease_attribution_declined"))
        return ("Exclude_disease_uncertain",
                (f"P/LP, and a source explicitly considered a CdLS/DEE85 "
                 f"attribution and declined it, so the condition cannot be "
                 f"resolved to either arm and is NOT pending further detail "
                 f"({v['disease_attribution_note']})" if declined else
                 f"P/LP but the recorded condition cannot be resolved to CdLS "
                 f"or DEE85 ({v['disease_attribution_note']})"),
                False,
                "disease_attribution_declined_by_source" if declined
                else "disease_recorded_but_unresolvable")
    if dis == "out_of_scope":
        return ("Exclude_disease_uncertain",
                f"P/LP but the recorded disease is outside the CdLS/DEE85 contrast "
                f"({v['disease_attribution_note']}); retained in the all list, "
                f"excluded from both arms",
                False, "disease_out_of_scope")
    if dis == "none":
        return ("Exclude_pending_classification",
                "P/LP with no disease attribution recorded by any source; "
                "attribution outstanding rather than ambiguous",
                False, "no_disease_attribution_yet")

    # 5. sufficient clinical detail (inclusion criterion 4)
    if not truthy(v["has_clinical_detail"]):
        return ("Exclude_disease_uncertain",
                f"{dis} attribution rests on a disease label with no supporting "
                f"clinical detail", False, "insufficient_clinical_detail")

    group = "CdLS_pathogenic" if dis == "CdLS" else "DEE85_pathogenic"
    # A large structural event keeps its biological group -- it IS a pathogenic
    # CdLS/DEE85 variant -- but is gated out of the curated set because the SGE
    # library cannot contain it.
    big, why = is_large_event(v)
    if big:
        return (group,
                f"{cls} and attributed to {dis}, but excluded from the curated "
                f"set as a structural event: {why}",
                False, "large_structural_event_not_assayable")
    return group, f"{cls} and attributed to {dis} with clinical detail", True, ""



# ---------------------------------------------------------------------------
# Residue-level mutational hotspots
#
# Distinct from both variant-identity overlap (same HGVS_c reported twice) and
# patient-identity overlap (same individual reported twice). This is the third
# category: the SAME amino acid residue altered by DIFFERENT variants in
# different individuals. It is not a duplicate -- the HGVS_c differs -- and by
# definition the individuals differ. Its value is as weak, independent evidence
# that the residue matters to protein function, separate from any one variant's
# own pathogenicity evidence.
#
# Computed by a systematic scan over every variant carrying a parseable protein
# position, so it does not depend on a source volunteering the connection.
# `residue_hotspot_groups` records which analysis arms the co-hitting variants
# fall into, because a residue hit by a PTV on the DEE85 side and a missense on
# the CdLS side is directly informative for the mechanism question.
# ---------------------------------------------------------------------------
P_POS = re.compile(r"^p\.\(?([A-Z][a-z]{2})(\d+)")


def annotate_residue_hotspots(variants):
    """Flag residues altered by more than one distinct variant."""
    by_res = {}
    for v in variants:
        m = P_POS.match(v.get("hgvs_p_mane") or "")
        v["_res"] = (m.group(1), int(m.group(2))) if m else None
        if v["_res"]:
            by_res.setdefault(v["_res"], []).append(v)

    n_hot = 0
    for v in variants:
        res = v.pop("_res", None)
        v["residue_hotspot"] = "False"
        v["residue_hotspot_n_variants"] = ""
        v["residue_hotspot_residue"] = ""
        v["residue_hotspot_groups"] = ""
        v["residue_hotspot_partners"] = ""
        if not res:
            continue
        peers = by_res[res]
        keys = sorted({p["variant_key"] for p in peers})
        if len(keys) < 2:
            continue
        n_hot += 1
        others = sorted({p["hgvs_c_mane"] for p in peers
                         if p["variant_key"] != v["variant_key"]})
        v["residue_hotspot"] = "True"
        v["residue_hotspot_n_variants"] = str(len(keys))
        v["residue_hotspot_residue"] = f"{res[0]}{res[1]}"
        v["residue_hotspot_groups"] = L.join_unique(
            [p.get("curation_group", "") for p in peers])
        v["residue_hotspot_partners"] = "|".join(others)
    return n_hot

def main():
    # ------------------------------------------------------------------
    # Load observations from every source
    # ------------------------------------------------------------------
    obs_files = sorted(glob.glob(os.path.join(INTERIM, "*_observations.tsv")))
    if not obs_files:
        sys.exit("FATAL: no *_observations.tsv found -- run steps 02 and 04 first")

    observations = []
    for f in obs_files:
        src = os.path.basename(f).replace("_observations.tsv", "")
        for r in read_tsv(f):
            r["_source_file"] = src
            observations.append(r)

    # Fail loudly on an unrecognised disease token. Silently treating one as
    # "unresolved" hides a registration error as a curation outcome.
    # Global individual -> family resolution. An explicit family assignment
    # from ANY source wins over the default of "each individual is its own
    # family", so a relationship stated by one source is honoured even where
    # another source citing the same patient does not state it.
    for o in observations:
        ik = L.clean(o.get("individual_key"))
        fk = L.clean(o.get("family_key"))
        if ik and fk and fk != ik:
            prev = FAMILY_OF.get(ik)
            if prev and prev != fk:
                print(f"WARNING: individual {ik} assigned to two families "
                      f"({prev}, {fk}) -- keeping {prev}")
                continue
            FAMILY_OF[ik] = fk

    # Runs before anything reads a disease token, so no downstream step can
    # see a claim this rule removes.
    suppress_redundant_database_attributions(observations, say)

    unknown = check_disease_vocabulary(observations)
    if unknown:
        for tok, srcs in sorted(unknown.items()):
            print(f"FATAL: unknown disease_harmonised token {tok!r} "
                  f"from {sorted(srcs)}")
        sys.exit("FATAL: unrecognised disease token(s) -- fix the source "
                 "registration or add the token to the vocabulary")

    # ------------------------------------------------------------------
    # Attach resolved coordinates
    # ------------------------------------------------------------------
    resolved, quarantined = {}, {}
    p = os.path.join(INTERIM, "resolved_variants.tsv")
    if os.path.exists(p):
        for r in read_tsv(p):
            resolved[r["hgvs_c_input"]] = r
    p = os.path.join(INTERIM, "resolution_quarantine.tsv")
    if os.path.exists(p):
        for r in read_tsv(p):
            quarantined[r["hgvs_c_input"]] = r

    targeton_status = {}
    for t in read_tsv(os.path.join(REF, "smc1a_targetons.tsv")):
        targeton_status[t["targeton_id"]] = t["screening_status"]

    say("# Registry build report")
    say()
    say(f"Pipeline `{PIPELINE_VERSION}`, curated {CURATION_DATE}.")
    say()
    say("| source ledger | observations |")
    say("|---|---|")
    for f in obs_files:
        say(f"| `{os.path.basename(f)}` | {len(read_tsv(f))} |")
    say(f"| **total** | **{len(observations)}** |")
    say()

    # ------------------------------------------------------------------
    # Normalise each observation onto a variant key
    # ------------------------------------------------------------------
    unkeyed = []
    reconciled = []
    for o in observations:
        c = L.clean(o.get("hgvs_c_input"))
        res = resolved.get(c) or quarantined.get(c)

        if L.clean(o.get("variant_key")):
            # The cohort sources supply plus-strand coordinates directly
            # (verified in 04 and 04d).
            coords = (o.get("chrom"), o.get("pos"), o.get("ref"), o.get("alt"))
            provenance = "source_vcf"
            o["_res_status"] = "resolved" if c in resolved else (
                "quarantined" if c in quarantined else "direct_coordinates")
            supplied = o["variant_key"]
        elif res and L.clean(res.get("variant_key")):
            coords = (res.get("chrom"), res.get("pos"),
                      res.get("ref"), res.get("alt"))
            provenance = "variant_recoder"
            o["_res_status"] = res["status"]
            supplied = res["variant_key"]
        else:
            o["_key"] = f"UNRESOLVED::{c or o.get('source_record_id')}"
            o["_coord_provenance"] = "none"
            o["_res_status"] = "quarantined"
            unkeyed.append(o)
            o["_res"] = res or {}
            continue

        # Reconcile the two coordinate provenances onto one key. A source's own
        # VCF coordinates are left-aligned as a caller produces them, whereas
        # the Variant Recoder returns the 3'-shifted HGVS position, so the same
        # variant reached by the two routes carried two different keys and did
        # not collapse: c.2202dup was present twice, once from GeneDx as
        # chrX:53403887:A:AT and once from LOVD as chrX:53403888:T:TT. Keying
        # on the normalised allele pair removes the distinction. This is what
        # `normalise_vcf` was written for; it had never been wired in.
        chrom, pos, ref, alt = coords
        if L.clean(pos) and L.clean(ref) and L.clean(alt):
            o["_key"] = L.canonical_variant_key(
                L.clean(chrom) or L.CHROM, int(pos), ref, alt, REFSEQ)
        else:
            o["_key"] = supplied
        o["_coord_provenance"] = provenance
        o["_res"] = res or {}
        if o["_key"] != supplied:
            reconciled.append((supplied, o["_key"], c, provenance))

    say(f"Observations that could not be placed on a genomic coordinate: "
        f"{len(unkeyed)} "
        f"(retained under an `UNRESOLVED::` key so nothing is dropped).")
    say()

    prov = Counter(o.get("_coord_provenance") for o in observations)
    say("Coordinate provenance of the keys, before reconciliation:")
    say()
    say("| provenance | observations |")
    say("|---|---|")
    for k, v in sorted(prov.items(), key=lambda kv: -kv[1]):
        say(f"| {k} | {v} |")
    say()
    say(f"Keys changed by normalisation: {len(reconciled)} observation(s). "
        f"Both provenances are keyed on the normalised allele pair, so a "
        f"variant supplied as VCF coordinates by one source and derived from "
        f"HGVS by another collapses to one row.")
    if reconciled:
        say()
        say("| as supplied | normalised | submitted HGVS c. | provenance |")
        say("|---|---|---|---|")
        for supplied, key, c, p in sorted(set(reconciled)):
            say(f"| `{supplied}` | `{key}` | {c} | {p} |")
    say()

    # ------------------------------------------------------------------
    # Collapse to variants
    # ------------------------------------------------------------------
    groups = defaultdict(list)
    for o in observations:
        groups[o["_key"]].append(o)

    # ---- overlay gate ------------------------------------------------
    # A source marked `overlay_only` (currently the Astorino 2025 review) may
    # attach attribution and clinical detail to a variant established
    # elsewhere, but may never introduce one. A variant supported ONLY by
    # overlay observations is therefore dropped from both output tables and
    # reported, so a review's transcription errors cannot create variants.
    overlay_dropped = []
    for key in list(groups):
        obs_ = groups[key]
        if obs_ and all(truthy(o.get("is_overlay_source")) for o in obs_):
            overlay_dropped.append((key, obs_[0].get("hgvs_c_input", ""),
                                    obs_[0].get("source", "")))
            del groups[key]

    variants = []
    for key, os_ in sorted(groups.items()):
        res = next((o["_res"] for o in os_ if o["_res"]), {})
        first = os_[0]

        def vals(field):
            return [L.clean(o.get(field)) for o in os_ if L.clean(o.get(field))]

        # ---- independent probands and families -----------------------
        # An observation counts as an independent proband only if it names an
        # individual. LOVD submitter classification records name none, so they
        # contribute evidence of classification but not a case.
        ind_keys, fam_keys, no_individual = set(), set(), 0
        needs_adj = False
        for o in os_:
            ik = L.clean(o.get("individual_key"))
            if not ik:
                no_individual += 1
                continue
            ind_keys.add(ik)
            # `family_key` is assigned upstream: step 02 by union-find over the
            # LOVD panelid links and the "<PMID>.<family><P|S>" submitter
            # labels; step 04 gives each GeneDx sample its own family. An empty
            # value means a singleton family.
            # Resolve the family through a GLOBAL individual -> family map.
            # The same individual can be described by two sources, only one of
            # which states a family relationship; taking each observation's own
            # family_key would then keep both values and count the individual
            # twice, inflating the family count rather than collapsing it.
            fam_keys.add(FAMILY_OF.get(ik) or L.clean(o.get("family_key")) or ik)
            if truthy(o.get("kinship_needs_adjudication")):
                needs_adj = True

        srcs = sorted({o["_source_file"] for o in os_})
        cross_source_dup = len(srcs) > 1
        # Independence is "assumed" where the source supplies no relatedness
        # field (currently GeneDx). Such counts are upper bounds.
        assumed = sorted({L.clean(o.get("proband_independence_basis"))
                          for o in os_ if L.clean(o.get("is_independent_proband")) == "assumed"})

        cls_consensus, cls_note = classification_consensus(
            vals("pathogenicity_reported"))
        dis_attr, dis_note = disease_consensus(vals("disease_harmonised"))
        # True where a source explicitly considered a CdLS/DEE85 attribution
        # and recorded that it did not apply, as distinct from never addressing
        # it. Read from the source's own per-case field, not inferred.
        declined_by = sorted({L.clean(o.get("source")) for o in os_
                              if L.clean(o.get("cdls_differential_reported")) == "no"})

        # ---- clinical detail (inclusion criterion 4) ------------------
        # Satisfied when at least one observation carries phenotype detail, a
        # recorded definite/initial diagnosis, or an episignature result --
        # i.e. something beyond a bare disease label.
        # HPO terms count. The two unpublished cohorts carry their phenotype as
        # structured HPO rather than free text, so a test limited to the
        # free-text fields silently failed every cohort observation on a
        # criterion they in fact satisfy -- with deep phenotyping, at that.
        detail = any(L.clean(o.get("phenotype_detail"))
                     or L.clean(o.get("episignature_result"))
                     or L.clean(o.get("individual_remarks"))
                     or L.clean(o.get("hpo_ids"))
                     or L.clean(o.get("hpo_names"))
                     for o in os_)

        # Sources emit either single letters (LOVD) or words (the cohorts).
        # Normalise before counting, or one vocabulary is silently discarded.
        SEXN = {"M": "M", "F": "F", "male": "M", "female": "F"}
        sexes = Counter(x for x in (SEXN.get(v) for v in vals("sex")) if x)
        inh = Counter(vals("inheritance"))

        _prot, _prot_note = pick_protein(res.get("hgvs_p_mane"), os_)

        v = {
            # ---- identity ----
            "variant_key": key,
            # Taken from the key, not from whichever provenance happened to
            # supply this group's first observation. Previously the key came
            # from the grouping key and these came from the resolver, so a row
            # keyed on a source's VCF coordinates carried the recoder's
            # coordinates in its own columns and the two disagreed.
            **_coords_from_key(key, res, first),
            "hgvs_g": L.clean(res.get("hgvs_g")),
            "hgvs_c_mane": (L.clean(res.get("hgvs_c_mane"))
                            or L.clean(first.get("hgvs_c_input"))),
            "hgvs_p_mane": _prot,
            "hgvs_p_note": _prot_note,
            "transcript": L.TRANSCRIPT,
            "protein": L.PROTEIN,
            "assembly": L.ASSEMBLY,
            "hgvs_c_as_submitted": L.join_unique(vals("hgvs_c_input")),
            "protein_position": (L.clean(res.get("protein_position"))
                                 or L.clean(first.get("protein_position"))),
            "ref_aa": L.clean(res.get("ref_aa")),
            "edit_type": L.clean(res.get("edit_type")),
            "consequence_class": (L.clean(res.get("consequence_class"))
                                  or L.clean(first.get("consequence_class"))),
            "consequence_terms": (L.clean(res.get("consequence_terms"))
                                  or L.clean(first.get("consequence_terms"))),
            "is_predicted_ptv": str((L.clean(res.get("consequence_class"))
                                     or L.clean(first.get("consequence_class"))
                                     ) in L.PTV_CLASSES),
            "exon": L.clean(res.get("exon")) or L.clean(first.get("exon")),
            "intron": L.clean(res.get("intron")),

            # ---- assay ----
            "targetons_design": (L.clean(res.get("targetons_design"))
                                 or L.clean(first.get("targetons_design"))),
            "targetons_amplicon": (L.clean(res.get("targetons_amplicon"))
                                   or L.clean(first.get("targetons_amplicon"))),
            "in_design_window": (L.clean(res.get("in_design_window"))
                                 or L.clean(first.get("in_design_window"))),

            # ---- evidence counts ----
            "n_observations": len(os_),
            # ONE PER FAMILY: the hotspot-safe count, so siblings cannot
            # inflate recurrence (METHODS 4.4). n_individuals is kept beside
            # it so the difference is visible rather than hidden.
            "n_independent_probands": len(fam_keys),
            "n_individuals": len(ind_keys),
            "n_families": len(fam_keys),
            "n_siblings_collapsed": len(ind_keys) - len(fam_keys),
            "kinship_needs_adjudication": str(needs_adj),
            "n_observations_without_individual": no_individual,
            "is_recurrent": str(len(fam_keys) > 1),
            "proband_count_is_upper_bound": str(bool(assumed) or cross_source_dup),
            "proband_independence_caveat": (
                "; ".join(assumed)[:400] if assumed else ""),
            "n_sources": len(srcs),
            "sources": "|".join(srcs),
            "source_records": L.join_unique(vals("source_record_id")),
            "individual_keys": "|".join(sorted(ind_keys)),
            "lovd_dbids": L.join_unique(vals("lovd_dbid")),
            # Parsed PMIDs, not the raw citation markup. `reference_raw` is
            # kept separately so the source wording is still recoverable.
            "pmids": ";".join(dict.fromkeys(
                pm for o in os_ for pm in L.clean(o.get("pmids")).split(";") if pm)),
            "citations": L.join_unique(vals("citations")),
            "dois": L.join_unique(vals("dois")),
            "omim_allelic_variants": L.join_unique(vals("omim_allelic_variants")),
            "reference_raw": L.join_unique(vals("reference_raw"))[:1000],
            "clinvar_ids": L.join_unique(vals("clinvar_id")),
            "dbsnp_id": L.join_unique(vals("dbsnp_id")),
            "published_as": L.join_unique(vals("published_as")),

            # ---- classification ----
            "pathogenicity_consensus": cls_consensus,
            "pathogenicity_all_reported": L.join_unique(
                vals("pathogenicity_reported_raw")),
            "pathogenicity_note": cls_note,
            "classification_methods": L.join_unique(vals("classification_method_raw")),

            # ---- disease ----
            "disease_attribution": dis_attr,
            "disease_attribution_declined": str(bool(declined_by)),
            "disease_attribution_declined_by": L.join_unique(declined_by),
            "disease_attribution_note": dis_note,
            "disease_all_reported": L.join_unique(vals("disease_reported_raw")),
            "has_clinical_detail": str(detail),
            # Yuan 2019 records per case whether CdLS was a differential
            # diagnosis. Where a source says it was NOT, that is surfaced here
            # rather than allowed to disappear behind another source's CdLS
            # attribution (D69).
            "cdls_differential_reported": L.join_unique(vals("cdls_differential_reported")),
            "cdls_presentation_absent_in_a_source": str(
                "no" in vals("cdls_differential_reported")),
            "dual_molecular_diagnosis": L.join_unique(vals("dual_molecular_diagnosis"))[:400],
            "zygosity_reported": L.join_unique(vals("zygosity_reported")),

            # ---- individual-level ----
            "sex_M": sexes.get("M", 0),
            "sex_F": sexes.get("F", 0),
            "sex_unknown": len(ind_keys) - sexes.get("M", 0) - sexes.get("F", 0),
            "inheritance_summary": "|".join(f"{k}:{n}" for k, n in inh.most_common() if k),
            "n_de_novo": inh.get("de_novo", 0),
            "phenotype_detail": L.join_unique(vals("phenotype_detail"))[:2000],
            "hpo_ids": L.join_unique(vals("hpo_ids")),
            "individual_remarks": L.join_unique(vals("individual_remarks"))[:1000],

            # ---- episignature (populated at the literature step) ----
            "episignature_available": str(bool(vals("episignature_result"))),
            "episignature_result": L.join_unique(vals("episignature_result")),
            "episignature_source": L.join_unique(vals("episignature_source")),

            # ---- functional evidence ----
            "functional_evidence": L.join_unique(vals("functional_evidence")),
            # provenance: was any contributing cDNA back-derived from a protein
            # change under the unique-route rule (D95)? Kept visible so a
            # back-derived coordinate is never mistaken for an extracted one.
            "hgvs_c_derived_from_protein": L.join_unique(
                vals("hgvs_c_derived_from_protein")),

            # ---- flags ----
            "not_smc1a_only": str(any(L.clean(o.get("not_smc1a_only")) == "True"
                                      for o in os_)),
            "alt_isoform_numbering": str(any(
                L.clean(o.get("alt_isoform_numbering")) == "True" for o in os_)),
            "possible_cross_source_duplicate": str(cross_source_dup),
            "resolution_status": next((o["_res_status"] for o in os_
                                       if o["_res_status"] != "resolved"),
                                      "resolved"),
            "guard_failed": L.clean(res.get("guard_failed")),
            "resolution_notes": L.clean(res.get("notes")),
            "notes": L.join_unique(vals("variant_remarks"))[:1000],
        }

        # targeton screening status: distinguishes "awaiting results" from
        # "absent from the library"
        ts = {targeton_status.get(t, "unknown")
              for t in v["targetons_design"].split("|") if t}
        v["targeton_screening_status"] = "|".join(sorted(ts)) if ts else "not_in_design_window"
        v["assay_testable_now"] = str(bool(ts) and "screened" in ts)

        grp, reason, include, excl = assign_group(v)

        # Contradicted CdLS attribution overrides an arm assignment. Applied
        # after assign_group so it can see the resolved group, and only ever
        # removes a variant from an arm -- it never admits one.
        contradicted, contra_note = cdls_attribution_contradicted(v, os_)
        v["cdls_attribution_contradicted"] = str(contradicted)
        v["cdls_attribution_note"] = contra_note
        if contradicted:
            grp, include, excl = "Exclude_disease_uncertain", False, \
                "cdls_attribution_contradicted"
            reason = ("P/LP but the CdLS attribution is contradicted: " + contra_note)
        v["curation_group"] = grp
        v["curation_group_reason"] = reason
        v["include_high_confidence"] = str(include)
        v["exclusion_reason"] = excl
        v["curation_date"] = CURATION_DATE
        v["pipeline_version"] = PIPELINE_VERSION
        variants.append(v)
        OBS_BY_VARIANT[v["variant_key"]] = os_

    check_key_coordinate_agreement(variants)
    apply_attribution_overrides(variants, say)

    # Q24: flag excluded-but-arguable variants for sensitivity analysis.
    # Runs AFTER assign_group so it can read exclusion_reason, and it never
    # alters curation_group or include_high_confidence.
    n_hot = annotate_residue_hotspots(variants)
    hot_res = sorted({(v["residue_hotspot_residue"], v["residue_hotspot_n_variants"],
                       v["residue_hotspot_groups"])
                      for v in variants if v["residue_hotspot"] == "True"},
                     key=lambda t: (-int(t[1]), t[0]))
    print(f"\nResidue-level hotspots: {len(hot_res)} residues altered by >1 distinct "
          f"variant ({n_hot} variant rows flagged)")
    for res, n, grps in hot_res:
        print(f"  {res:<9} {n} distinct variants   "
              f"{grps.replace('_pathogenic', '')}")

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------
    say(f"## Collapse")
    say()
    say(f"- observations: {len(observations)}")
    say(f"- distinct variants: {len(variants)}")
    say(f"- variants seen in >1 independent proband: "
        f"{len([v for v in variants if truthy(v['is_recurrent'])])}")
    say(f"- variants seen in >1 source: "
        f"{len([v for v in variants if v['n_sources'] > 1])}")
    say()

    say("## Analysis groups")
    say()
    say("| curation_group | variants | independent probands |")
    say("|---|---|---|")
    gc = Counter(v["curation_group"] for v in variants)
    gp = defaultdict(int)
    for v in variants:
        gp[v["curation_group"]] += v["n_independent_probands"]
    for g, n in gc.most_common():
        say(f"| `{g}` | {n} | {gp[g]} |")
    say()

    say("### Why variants were excluded")
    say()
    say("| exclusion_reason | variants |")
    say("|---|---|")
    for r, n in Counter(v["exclusion_reason"] for v in variants
                        if v["exclusion_reason"]).most_common():
        say(f"| `{r}` | {n} |")
    say()

    curated = [v for v in variants if truthy(v["include_high_confidence"])]
    say("## Attribution provenance per arm")
    say()
    report_attribution_provenance(variants, OBS_BY_VARIANT, say)

    say("## Per-source criterion coverage")
    say()
    report_source_criterion_coverage(variants, say)

    say("## High-confidence curated set")
    say()
    say(f"**{len(curated)} variants** meeting all inclusion criteria "
        f"(METHODS 5.2).")
    say()
    if curated:
        # PTV is reported split into `truncating` (frameshift/nonsense) and
        # `splice`, because lumping them misrepresents the splice class here: of
        # the splice variants in the curated set, one has a predicted IN-FRAME
        # protein consequence (`c.2974-2A>G` -> `p.Asp992_Gln994del`) and one is
        # synonymous at the protein level (`c.1731G>A` -> `p.Glu577=`), where
        # "splice-affecting" is a prediction rather than an observation. Neither
        # is established loss of function, so counting them alongside nonsense
        # and frameshift overstates the PTV burden of an arm.
        SPLICE = {"splice_donor", "splice_acceptor", "splice_region"}
        TRUNC = {"frameshift", "nonsense"}
        say("| group | variants | probands | truncating | splice | other | "
            "in design window |")
        say("|---|---|---|---|---|---|---|")
        for g in ("CdLS_pathogenic", "DEE85_pathogenic"):
            sub = [v for v in curated if v["curation_group"] == g]
            if not sub:
                continue
            say(f"| `{g}` | {len(sub)} | "
                f"{sum(v['n_independent_probands'] for v in sub)} | "
                f"{len([v for v in sub if v['consequence_class'] in TRUNC])} | "
                f"{len([v for v in sub if v['consequence_class'] in SPLICE])} | "
                f"{len([v for v in sub if v['consequence_class'] not in TRUNC | SPLICE])} | "
                f"{len([v for v in sub if truthy(v['in_design_window'])])} |")
        say()
        say("Consequence breakdown:")
        say()
        say("| group | " + " | ".join(
            sorted({v["consequence_class"] for v in curated if v["consequence_class"]}))
            + " |")
        cs = sorted({v["consequence_class"] for v in curated if v["consequence_class"]})
        say("|---" * (len(cs) + 1) + "|")
        for g in ("CdLS_pathogenic", "DEE85_pathogenic"):
            sub = [v for v in curated if v["curation_group"] == g]
            if not sub:
                continue
            c = Counter(v["consequence_class"] for v in sub)
            say(f"| `{g}` | " + " | ".join(str(c.get(k, 0)) for k in cs) + " |")
    say()

    say("## Assay testability of the curated set")
    say()
    say("| targeton_screening_status | curated variants |")
    say("|---|---|")
    for s_, n in Counter(v["targeton_screening_status"] for v in curated).most_common():
        say(f"| {s_} | {n} |")
    say()
    say("`not_screened` means CQEJ (exon 6, second tiling window) or NLVE "
        "(exon 23): designed and present in the library, but not yet "
        "screened, so no assay result exists. A variant carrying it is not "
        "excluded from the curated set -- it is simply not assay-testable "
        "yet, and `assay_testable_now` is False.")
    say()

    say("## Cross-source duplication")
    say()
    dups = [v for v in variants if truthy(v["possible_cross_source_duplicate"])]
    say(f"{len(dups)} variants are reported by more than one source. LOVD "
        f"individuals and GeneDx probands cannot be linked (GeneDx samples are "
        f"anonymised; LOVD identifiers are submitter-local), and no fuzzy "
        f"matching is attempted. For these variants "
        f"`n_independent_probands` is therefore an **upper bound**.")
    if dups:
        say()
        say("| variant | sources | probands (upper bound) |")
        say("|---|---|---|")
        for v in sorted(dups, key=lambda x: -x["n_independent_probands"])[:25]:
            say(f"| `{v['hgvs_c_mane']}` | {v['sources']} | "
                f"{v['n_independent_probands']} |")
    say()

    say("## Effect of family grouping on recurrence counts")
    say()
    coll = [v for v in variants if int(v["n_siblings_collapsed"]) > 0]
    say(f"{len(coll)} variants had sibling observations collapsed, so their "
        f"recurrence count is lower than the number of individuals reported.")
    if coll:
        say()
        say("| variant | individuals | independent probands (one per family) | collapsed |")
        say("|---|---|---|---|")
        for v in sorted(coll, key=lambda x: -int(x["n_siblings_collapsed"])):
            say(f"| `{v['hgvs_c_mane']}` | {v['n_individuals']} | "
                f"{v['n_independent_probands']} | {v['n_siblings_collapsed']} |")
    say()
    adj = [v for v in variants if truthy(v["kinship_needs_adjudication"])]
    say(f"{len(adj)} variants involve at least one record whose kinship is stated "
        f"only in free text and so **needs manual adjudication** "
        f"(`kinship_needs_adjudication`); their proband counts may fall further.")
    say()

    say("## Recurrence (independent probands, one per family)")
    say()
    rec = sorted([v for v in variants if v["n_independent_probands"] > 1],
                 key=lambda v: -v["n_independent_probands"])
    say(f"{len(rec)} variants seen in more than one independent proband.")
    if rec:
        say()
        say("| variant | protein | probands | individuals | adjudication? | consequence | group |")
        say("|---|---|---|---|---|---|---|")
        for v in rec[:30]:
            say(f"| `{v['hgvs_c_mane']}` | {v['hgvs_p_mane'] or '-'} | "
                f"{v['n_independent_probands']} | {v['n_individuals']} | "
                f"{'**yes**' if truthy(v['kinship_needs_adjudication']) else 'no'} | "
                f"{v['consequence_class']} | `{v['curation_group']}` |")
    say()

    # ------------------------------------------------------------------
    # Write
    # ------------------------------------------------------------------
    # Column order comes from smc1a_schema, which also generates the data
    # dictionary -- so the two cannot drift apart. Columns absent from the
    # schema are appended rather than dropped, and reported.
    for o in observations:
        o["variant_key_resolved"] = o["_key"]
    present = {k for o in observations for k in o if not k.startswith("_")}
    ledger_cols, extra_obs = S.order_columns(present, S.OBSERVATIONS)
    with open(os.path.join(PROC, "smc1a_observations.tsv"), "w",
              encoding="utf-8") as fh:
        fh.write("\t".join(ledger_cols) + "\n")
        for o in observations:
            fh.write("\t".join(
                str(o.get(c, "")).replace("\t", " ").replace("\n", " ")
                for c in ledger_cols) + "\n")

    vcols, extra_var = S.order_columns(set(variants[0].keys()), S.VARIANTS)
    for name, data in (("smc1a_variants_all.tsv", variants),
                       ("smc1a_variants_curated.tsv", curated)):
        with open(os.path.join(PROC, name), "w", encoding="utf-8") as fh:
            fh.write("\t".join(vcols) + "\n")
            for v in data:
                fh.write("\t".join(
                    str(v.get(c, "")).replace("\t", " ").replace("\n", " ")
                    for c in vcols) + "\n")

    say("## Schema coverage")
    say()
    say("Column order and `docs/DATA_DICTIONARY.md` are both generated from "
        "`code/smc1a_schema.py`, so the documentation cannot drift from the "
        "files.")
    say()
    if extra_var or extra_obs:
        say("Columns present in the output but **not documented in the "
            "schema** (appended at the end of the file, never dropped):")
        say()
        for name, ex in (("variant tables", extra_var),
                         ("observation ledger", extra_obs)):
            if ex:
                say(f"- {name}: " + ", ".join(f"`{c}`" for c in sorted(ex)))
        say()
    else:
        say("All output columns are documented in the schema.")
        say()

    # ---- structural integrity assertion --------------------------------
    # A note field containing a newline once split one logical row across
    # several physical lines, so the quarantine file held 50 lines for 16
    # variants. Any reader would have mis-parsed it silently. This check makes
    # that class of corruption fail loudly instead of passing unnoticed.
    say("## Overlay-source gate")
    say()
    say("Sources marked `overlay_only` may attach attribution and clinical "
        "detail to variants established elsewhere, but may not introduce a "
        "variant or a coordinate. A variant supported only by overlay "
        "observations is dropped and listed here.")
    say()
    say(f"- variants dropped as overlay-only: **{len(overlay_dropped)}**")
    if overlay_dropped:
        say()
        say("| submitted c. | source |")
        say("|---|---|")
        for k, c, src in sorted(overlay_dropped, key=lambda x: x[1]):
            say(f"| `{c or k}` | {src} |")
    say()

    say("## Structural integrity of the written files")
    say()
    say("| file | rows | columns | ragged rows |")
    say("|---|---|---|---|")
    integrity_ok = True
    for name in ("smc1a_observations.tsv", "smc1a_variants_all.tsv",
                 "smc1a_variants_curated.tsv"):
        with open(os.path.join(PROC, name), encoding="utf-8") as fh:
            lines = fh.read().rstrip("\n").split("\n")
        ncol = lines[0].count("\t") + 1
        ragged = [i + 2 for i, l in enumerate(lines[1:])
                  if l.count("\t") + 1 != ncol]
        integrity_ok = integrity_ok and not ragged
        say(f"| `{name}` | {len(lines) - 1} | {ncol} | "
            f"{len(ragged) if ragged else '0'} |")
        if ragged:
            say(f"  FAIL: lines {ragged[:10]} have the wrong column count")
    say()
    if integrity_ok:
        say("All files are rectangular: every row has exactly the header's "
            "column count.")
    say()

    say("## Files written")
    say()
    say(f"| file | rows | columns |")
    say(f"|---|---|---|")
    say(f"| `data/processed/smc1a_observations.tsv` | {len(observations)} | "
        f"{len(ledger_cols)} |")
    say(f"| `data/processed/smc1a_variants_all.tsv` | {len(variants)} | {len(vcols)} |")
    say(f"| `data/processed/smc1a_variants_curated.tsv` | {len(curated)} | {len(vcols)} |")

    with open(os.path.join(DOCS, "registry_build.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(report) + "\n")
    print("\nWrote docs/verification/registry_build.md")
    if not integrity_ok:
        sys.exit("FATAL: a written file is not rectangular -- see "
                 "docs/verification/registry_build.md. Almost certainly an "
                 "unescaped newline or tab in a text field.")


if __name__ == "__main__":
    main()
