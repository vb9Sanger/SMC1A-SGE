#!/usr/bin/env python3
"""
06_join_assay.py -- join the curated variant set to the SGE screen results.

Amended from `variant_curation/overlap/truthset_sge_overlap.py`, which was
written against the historical spreadsheet truth set. The curated set produced
by step 05 differs from that spreadsheet in ways that change the matching
problem, so the matching logic here is not a port of the old one:

  * the old truth set stored variants as HGVS c. on the transcript strand, and
    the old script therefore complemented every allele before matching, SMC1A
    being on the minus strand. The curated set stores `pos_grch38`/`ref`/`alt`
    in plus-strand GRCh38 VCF form -- the same convention as the oligo names --
    so no complementation happens here at all, and the `--strand` switch is
    gone. Complementing now would silently produce wrong keys for every
    variant.
  * the old script could only match single-nucleotide substitutions exactly,
    because the spreadsheet carried no genomic alleles for indels; deletions,
    duplications and delins all fell back to "position_only_verify_manually".
    The curated set carries normalised alleles, so oligo names are decoded
    against the reference and left-aligned, so deletion and insertion oligos
    match exactly too and position-only matching is a last resort rather than
    the routine outcome for indels.
  * variants absent from an oligo file are separated by cause -- a targeton
    whose results were not supplied, versus a variant genuinely not built into
    the library -- instead of a single `no_match`.

Two screen-design properties are checked empirically rather than assumed,
following the pattern set by `smc1a_lib.TargetonMap`:

  1. **Oligo orientation and alignment.** Oligo sequences are aligned to the
     plus-strand reference to locate the amplicon, rather than trusting a
     filename or a coordinate column.
  2. **Fixed background edits.** Every library member in a targeton carries the
     same PAM/protospacer-disrupting edits, which are not the variant under
     test. In APDY these are chrX:53,415,051 G>A, chrX:53,415,054 G>A and
     chrX:53,415,057 A>G, present in 1870/1871 oligos that cover them. Two
     consequences are handled:
       - at such a position the oligo *name* is written against the edited
         background, so its stated reference base is not the genomic reference
         and a key derived from it would be wrong;
       - a curated variant falling on one of these positions needs care. A
         *substitution* there is still exactly interpretable, since what the
         cell carries is the oligo's own base, and is re-expressed against the
         genome; an *indel* there cannot be separated from the edit and is
         withheld from the analysis set.
     Both are detected per targeton from the oligo sequences themselves and
     reported; no position list is hard-coded.
  3. **Positional coverage by oligo class.** Substitutions are near-saturated
     in every design window (61-100% of positions), but single-base deletions
     reach 64% only in APDY and 5-24% elsewhere, and duplications and
     insertions are not designed at all. Since the DEE85 arm is
     protein-truncating and so indel-heavy while the CdLS arm is missense, this
     screens 80% of the CdLS arm against 59% of DEE85, for reasons of library
     design rather than biology. Reported per targeton and cross-tabulated by
     arm and edit type.

Matching
--------
`exact_allele`
    The curated `variant_key` equals the oligo's canonical key, after decoding
    the oligo name against the reference and left-normalising it. This is the
    only match type admitted to the analysis tables.

`position_only_verify_manually`
    The curated variant starts at the position of an oligo whose own allele
    could not be established -- an undecodable name, or a name written against
    a background edit -- so allele agreement can be neither confirmed nor
    ruled out. Written to the full table, flagged, and excluded from the
    analysis tables.

    This is deliberately much narrower than in the old script, where every
    indel fell here. Position-only matching against an oligo whose allele *is*
    known and demonstrably different is not a weak match but a wrong one: it
    attaches an unrelated oligo's depletion score to the curated variant. In
    APDY the 15 bp in-frame deletion c.173_187del shares a start coordinate
    with four well-formed oligos, one of them strongly depleting, and none of
    them that deletion; reading across such a row would attribute a depletion
    result to a variant the library never contained. Candidate oligos with a
    known, non-matching allele are therefore not reported as matches.

`no_oligo_in_supplied_file` / `no_oligo_targeton_not_supplied`
    Not matched. The first means the targeton's results were read and the
    variant is not in them: the expected outcome for duplications and
    insertions, which are not designed, and common for deletions, which are
    built at only a minority of positions (point 3 above). `site_built_in_library`
    separates an allele not designed at a covered site from a site the library
    does not reach. The second means the variant cannot yet be assessed because
    that targeton's results are absent.

Inputs
------
data/processed/smc1a_variants_curated.tsv   (step 05; --variants to override)
data/reference/smc1a_chrX_region.fa         (step 01)
data/reference/smc1a_targetons.tsv
variant_curation/sge_results/*_all_deseq2_results_condition_*.tsv
                                            (22 targetons; --sge_glob to
                                            override. Read-only.)

Outputs (--outdir, default data/assay_join)
-------------------------------------------
assay_join_all.tsv          one row per curated-variant x oligo match, any
                            confidence, with clinical and assay columns
assay_join_summary.tsv      one row per curated variant, matched or not, with
                            its match outcome -- the "was this variant
                            screened" table
assay_join_analysis.tsv     the exact_allele subset carrying a curation_group,
                            restricted by --class_filter, one row per variant
assay_join_summary_stats.tsv  per-group depletion counts and Fisher's exact
assay_join_depletion.png    stripplot and tier breakdown, CdLS vs DEE85
                            (written only if matplotlib is installed)
docs/verification/assay_join.md

Requirements: pandas, numpy. matplotlib is optional -- every table is written
without it and only the figure is skipped. scipy is optional and only the
Fisher's exact test depends on it.

Usage
-----
    ../.venv/bin/python 06_join_assay.py
"""

from __future__ import annotations

import argparse
import datetime
import glob as globmod
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smc1a_lib as L  # noqa: E402

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False

try:
    from scipy.stats import fisher_exact
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, "data", "reference")
DOCS = os.path.join(ROOT, "docs", "verification")

GROUPS = ["CdLS_pathogenic", "DEE85_pathogenic"]

# Columns carried through from the curated set. Deliberately a fixed list: a
# silent schema change upstream should surface as a missing-column error here
# rather than as a join that quietly drops a field.
VARIANT_COLS = [
    "variant_key", "chrom", "pos_grch38", "ref", "alt",
    "hgvs_c_mane", "hgvs_p_mane", "protein_position",
    "edit_type", "consequence_class", "is_predicted_ptv", "exon", "intron",
    "targetons_design", "targeton_screening_status", "assay_testable_now",
    "in_design_window",
    "curation_group", "curation_group_reason", "include_high_confidence",
    "pathogenicity_consensus", "disease_attribution",
    "n_individuals", "n_independent_probands", "n_families", "is_recurrent",
    "sex_M", "sex_F", "sex_unknown",
    "residue_hotspot", "residue_hotspot_residue",
    "literature_domain_hotspot", "literature_domain_hotspot_domain",
    "literature_domain_hotspot_source",
    "functional_evidence", "episignature_result",
    "n_sources", "sources", "pmids",
    "resolution_status", "notes",
]

# Assay columns kept from the DESeq2 results. `pos_adj_*` are the
# position-adjusted statistics used as the project's readout; the unadjusted
# and shrunk columns are kept for reference but are not the default score.
SGE_KEEP_COLS = [
    "oligo_name", "position", "consequence",
    "baseMean_raw", "log2FoldChange_raw", "pvalue_raw", "padj_raw",
    "log2FoldChange_shrunk", "pvalue_shrunk", "padj_shrunk",
    "adj_log2FoldChange_raw", "adj_pval_raw", "adj_fdr_raw", "stat_adj_raw",
    "pos_adj_log2FoldChange_raw", "pos_adj_pval_raw", "pos_adj_fdr_raw",
    "stat_pos_raw",
    "pos_adj_log2FoldChange_shrunk", "pos_adj_fdr_shrunk", "stat_pos_shrunk",
    "anchor_post_lof", "anchor_call", "anchor_tier", "anchor_direction",
    "anchor_lof_threshold",
]

SCORE_PREFERENCE = ["pos_adj_log2FoldChange_raw", "adj_log2FoldChange_raw",
                    "log2FoldChange_raw", "anchor_post_lof"]
TIER_PREFERENCE = ["anchor_tier", "anchor_call", "stat_adj_raw"]

report: list[str] = []


def say(m: str = "") -> None:
    print(m)
    report.append(m)


# --------------------------------------------------------------------------
# Allele normalisation
# --------------------------------------------------------------------------

def left_align(pos: int, ref: str, alt: str, refseq):
    """Left-align and trim a plus-strand allele pair.

    Delegates to `smc1a_lib.normalise_vcf`, which is also what step 05 keys on,
    so the join and the variant table cannot drift apart on allele
    representation. Kept as a named function here because both sides of the
    join go through it and the intent is worth stating at the call sites.

    Canonicalising both sides is what makes indel matching work at all. A
    duplication reaches the oligo list as an insertion at whichever position
    the library designer chose and reaches the variant table at whichever
    position its source used; only after left-alignment do the two agree. The
    same applies to single-base deletions in a homopolymer, where the oligo
    name gives the designed base and the variant table an anchored pair four
    bases away.
    """
    return L.normalise_vcf(pos, ref, alt, refseq)


def canonical_key(pos: int, ref: str, alt: str, refseq) -> str:
    return L.canonical_variant_key(L.CHROM, pos, ref, alt, refseq)


# --------------------------------------------------------------------------
# Oligo name decoding
# --------------------------------------------------------------------------

_DEL_TAG = re.compile(r"^(\d+)del\d*$")


def parse_oligo_name(name: str):
    """Decode the `chrX:...` segment of an oligo name into a structured dict.

    The grammar observed in the SGE results, with counts from APDY:

        POS_REF>ALT_snv                654   single-nucleotide substitution
        POS1_POS2_REF>ALT_aa           814   codon replacement
        POS1_POS2_REF>ALT_snvre        369   multi-nucleotide substitution
        POS1_POS2_REF>ALT_stop          51   stop-introducing codon replacement
        POS_1del                       155   single-base deletion
        POS1_POS2_2del0 / _2del1        13   two-base deletion
        POS1_POS2_inframe               61   three-base in-frame deletion
        POS[_POS2]_{clinvar,gnomAD_*}     9   seeded known variants

    Returns None when there is no decodable `chrX:` segment. `kind` is one of
    `sub`, `del` or `unknown`; `unknown` carries a span but no allele and falls
    through to position-only matching.
    """
    i = name.find("chrX:")
    if i == -1:
        return None
    parts = name[i + len("chrX:"):].split("_")
    if not parts or not parts[0].isdigit():
        return None

    pos1 = int(parts[0])
    rest = parts[1:]
    pos2 = None
    if rest and rest[0].isdigit():
        pos2 = int(rest[0])
        rest = rest[1:]
    tag = rest[-1] if rest else ""

    allele = next((t for t in rest if ">" in t), None)
    if allele is not None:
        ref, alt = allele.split(">", 1)
        return {"kind": "sub", "pos1": pos1, "pos2": pos2,
                "ref": ref, "alt": alt, "tag": tag}

    # Deletion families. The deleted span is pos1..pos2 inclusive, or pos1
    # alone where no second coordinate is given. `del_len` is checked against
    # the oligo sequence length by verify_deletion_lengths().
    if _DEL_TAG.match(tag) or tag == "inframe":
        end = pos2 if pos2 is not None else pos1
        return {"kind": "del", "pos1": pos1, "pos2": end,
                "ref": None, "alt": None, "tag": tag}

    return {"kind": "unknown", "pos1": pos1,
            "pos2": pos2 if pos2 is not None else pos1,
            "ref": None, "alt": None, "tag": tag}


def oligo_key_from_name(p: dict, refseq, background: set):
    """Canonical key implied by the oligo *name*, used only as a cross-check on
    the sequence-derived key. Returns (key_or_None, reason).

    Substitution names are verified against the reference before keying,
    because a stated reference base that does not match the genome means the
    name is written against something other than GRCh38 -- in practice a fixed
    background edit -- so a key built from it would be wrong.
    """
    if p is None:
        return None, "no_chrX_segment"

    if p["kind"] == "sub":
        end = p["pos2"] if p["pos2"] is not None else p["pos1"] + len(p["ref"]) - 1
        if end - p["pos1"] + 1 != len(p["ref"]):
            return None, "name_span_disagrees_with_allele_length"
        observed = refseq.get(p["pos1"], end)
        if not observed:
            return None, "outside_cached_reference"
        if observed != p["ref"]:
            span = set(range(p["pos1"], end + 1))
            return None, ("named_on_background_edit" if span & background
                          else "reference_base_mismatch")
        return canonical_key(p["pos1"], p["ref"], p["alt"], refseq), "ok"

    if p["kind"] == "del":
        deleted = refseq.get(p["pos1"], p["pos2"])
        if not deleted:
            return None, "outside_cached_reference"
        return canonical_key(p["pos1"], deleted, "", refseq), "ok"

    return None, "oligo_name_not_decodable"


def oligo_key_from_sequence(seq, bg_seq, start, background, refseq):
    """Canonical key for an oligo derived from its own sequence.

    This, not the name, is the primary route. Oligo names follow a systematic
    grammar only for the systematically designed members; variants seeded from
    ClinVar and gnomAD carry names that do not encode their allele at all. In
    APDY two curated variants -- the 15 bp in-frame deletion c.173_187del and
    the duplication c.157dup -- are present in the library as
    `53415092_53415106_clinvar` and `53415122_T_clinvar`, and neither name is
    decodable. Matching on names alone reported both as absent from the
    library, one of them with four unrelated oligos at the same coordinate.
    Sequences are unambiguous and are read the same way for every member.

    The oligo is compared against the *background* sequence -- the reference
    with that targeton's fixed edits applied -- so that the PAM/protospacer
    edits carried by every member are not mistaken for the variant under test.
    What remains after trimming the shared prefix and suffix is the variant,
    which is then put through `left_align()`, so that an indel placed anywhere
    within a repeat converges on the same key.
    """
    if seq == bg_seq:
        return None, "no_variant_relative_to_background"

    i = 0
    while i < min(len(seq), len(bg_seq)) and seq[i] == bg_seq[i]:
        i += 1
    j = 0
    while j < min(len(seq), len(bg_seq)) - i and seq[-1 - j] == bg_seq[-1 - j]:
        j += 1

    ref = bg_seq[i:len(bg_seq) - j]
    alt = seq[i:len(seq) - j]
    pos = start + i

    span = set(range(pos, pos + len(ref))) if ref else {pos - 1, pos}
    if span & background:
        # The designed variant sits on a background-edited base, so `ref` above
        # is the edited base rather than the genome's. For a **substitution**
        # this is still unambiguous: what the cell carries at that position is
        # simply the oligo's base, and the variant relative to the genome is
        # that base against the genomic reference. CNJV's c.3321C>A is exactly
        # this case -- the background carries A at chrX:53,382,348 where the
        # genome has G, and the oligo `53382348_A>T_snv` carries T, so the
        # oligo is chrX:53382348:G:T, which is the curated variant. Refusing to
        # key it lost a curated DEE85 nonsense variant that had in fact been
        # screened.
        #
        # Re-expressed against the genome, therefore, whenever the change is a
        # pure substitution. For an **indel** overlapping an edited base there
        # is no equivalent reading -- the edit falls inside or beside the
        # altered span and the two cannot be separated -- so those stay
        # flagged for manual inspection.
        if len(ref) != len(alt) or not ref:
            return None, "indel_overlaps_background_edit"
        genomic_ref = refseq.get(pos, pos + len(ref) - 1)
        if not genomic_ref:
            return None, "outside_cached_reference"
        if genomic_ref == alt:
            # The oligo restores the reference base: it tests the background,
            # not a variant.
            return None, "reverts_background_edit"
        p, r, a = left_align(pos, genomic_ref, alt, refseq)
        return (L.variant_key(L.CHROM, p, r, a),
                "ok_re_expressed_against_genome")

    if ref and refseq.get(pos, pos + len(ref) - 1) != ref:
        return None, "background_reference_disagrees"

    p, r, a = left_align(pos, ref, alt, refseq)
    # Left-alignment can shift the variant onto a background-edited base that
    # the pre-normalisation span did not touch, so the test is repeated after.
    if set(range(p, p + len(r))) & background:
        return None, "variant_overlaps_background_edit"
    return L.variant_key(L.CHROM, p, r, a), "ok"


# --------------------------------------------------------------------------
# Empirical checks on each oligo file
# --------------------------------------------------------------------------

def align_amplicon(sequences, refseq, search_lo, search_hi):
    """Locate the plus-strand start of the modal-length oligo by direct
    alignment, and return (start, modal_length, median_mismatches).

    Orientation is established here rather than assumed. Oligo sequences for
    this library are plus-strand: reverse-complemented alignment of an APDY
    oligo gives 181 mismatches over 290 bases against 4 for the forward
    orientation, so a wrong assumption is not subtle, but it is still checked
    so that a future library built the other way round fails loudly.
    """
    lens = Counter(len(s) for s in sequences)
    modal = lens.most_common(1)[0][0]
    pool = [s for s in sequences if len(s) == modal][:40]
    if not pool:
        return None, modal, None

    def best_offset(queries):
        votes = Counter()
        scores = []
        for q in queries:
            best = None
            for start in range(search_lo, search_hi + 1):
                ref = refseq.get(start, start + modal - 1)
                if len(ref) != modal:
                    continue
                nm = sum(1 for a, b in zip(ref, q) if a != b)
                if best is None or nm < best[1]:
                    best = (start, nm)
            if best is not None:
                votes[best[0]] += 1
                scores.append(best[1])
        if not votes:
            return None, None
        return votes.most_common(1)[0][0], int(np.median(scores))

    fwd_start, fwd_nm = best_offset(pool)
    rev_start, rev_nm = best_offset([L.revcomp(s) for s in pool])
    if fwd_nm is None:
        return None, modal, None
    if rev_nm is not None and rev_nm < fwd_nm:
        raise SystemExit(
            "ERROR: oligo sequences align better reverse-complemented "
            f"({rev_nm} vs {fwd_nm} median mismatches). This code assumes "
            "plus-strand oligo sequences; the decoding would be wrong.")
    return fwd_start, modal, fwd_nm


def detect_background_edits(names, sequences, refseq, start, modal,
                            min_fraction=0.9):
    """Positions where nearly every oligo differs from the reference.

    Each oligo's own named span is excluded from its votes, so the variant
    under test cannot be mistaken for a background edit. Positions between 5%
    and `min_fraction` are returned separately: a background edit should be
    present in essentially all oligos or none, and an intermediate frequency
    means something about the library is not as assumed.
    """
    ref = refseq.get(start, start + modal - 1)
    mismatch, coverage = Counter(), Counter()
    votes = defaultdict(Counter)
    for name, seq in zip(names, sequences):
        if len(seq) != modal:
            continue
        p = parse_oligo_name(name)
        span = range(p["pos1"], (p["pos2"] or p["pos1"]) + 1) if p else range(0)
        lo, hi = (min(span), max(span)) if len(span) else (None, None)
        for i in range(modal):
            pos = start + i
            if lo is not None and lo <= pos <= hi:
                continue
            coverage[pos] += 1
            if ref[i] != seq[i]:
                mismatch[pos] += 1
                votes[pos][seq[i]] += 1

    fixed, ambiguous = {}, {}
    for pos, n in mismatch.items():
        frac = n / coverage[pos]
        if frac >= min_fraction:
            fixed[pos] = (refseq.base(pos), votes[pos].most_common(1)[0][0],
                          n, coverage[pos])
        elif frac > 0.05:
            ambiguous[pos] = (refseq.base(pos), votes[pos].most_common(1)[0][0],
                              n, coverage[pos])
    return fixed, ambiguous


def report_library_coverage(metas, say_fn):
    """Per-targeton positional coverage by oligo class.

    Coverage is not uniform across classes, and the difference is large enough
    to bias an arm-level comparison. Substitutions are near-saturated in every
    design window; single-base deletions are built at 64% of positions in APDY
    but only 12-24% elsewhere, and duplications and insertions are not designed
    at all. Since the DEE85 arm is protein-truncating and therefore
    indel-heavy, while the CdLS arm is missense and therefore substitution-
    heavy, unmatched variants fall disproportionately on DEE85 for reasons of
    library design rather than biology. Reported so that a difference in how
    much of each arm was screened cannot be mistaken for a difference in how
    each arm behaves.
    """
    say_fn("### Positional coverage by oligo class")
    say_fn()
    say_fn("Fraction of positions in the substitution span carrying at least "
           "one oligo of each class.")
    say_fn()
    say_fn("| targeton | span (bp) | substitution | 1 bp deletion | "
           "2 bp deletion | 3 bp in-frame |")
    say_fn("|---|---|---|---|---|---|")
    for m in metas:
        c = m.get("coverage")
        if not c:
            continue
        say_fn(f"| {m['targeton']} | {c['span']} | {c['sub']:.0%} | "
               f"{c['del1']:.0%} | {c['del2']:.0%} | {c['inframe']:.0%} |")
    say_fn()


def measure_coverage(names):
    """Positions covered by each oligo class, as a fraction of the span over
    which substitution oligos were designed."""
    sub, del1, del2, inframe = set(), set(), set(), set()
    for n in names:
        p = parse_oligo_name(str(n))
        if p is None:
            continue
        if p["kind"] == "sub":
            end = p["pos2"] if p["pos2"] is not None else p["pos1"]
            sub.update(range(p["pos1"], end + 1))
        elif p["kind"] == "del":
            span = range(p["pos1"], p["pos2"] + 1)
            n_del = p["pos2"] - p["pos1"] + 1
            (del1 if n_del == 1 else del2 if n_del == 2 else inframe).update(span)
    if not sub:
        return None
    lo, hi = min(sub), max(sub)
    span = hi - lo + 1
    return {"span": span, "sub": len(sub) / span, "del1": len(del1) / span,
            "del2": len(del2) / span, "inframe": len(inframe) / span}


def report_unmatched_by_class(summary, say_fn):
    """Cross-tabulate the match outcome by disease arm and edit type, so a
    coverage asymmetry between the arms is visible as a number."""
    sub = summary[summary["curation_group"].isin(GROUPS)]
    if sub.empty:
        return
    say_fn("Match outcome by arm and edit type, which is where the coverage "
           "asymmetry above shows up:")
    say_fn()
    ct = pd.crosstab([sub["curation_group"], sub["edit_type"]],
                     sub["match_type"])
    say_fn("```")
    say_fn(ct.to_string())
    say_fn("```")
    say_fn()
    for g in GROUPS:
        arm = sub[sub["curation_group"] == g]
        screened = (arm["match_type"] == "exact_allele").sum()
        say_fn(f"- {g}: {screened}/{len(arm)} screened "
               f"({screened / len(arm):.0%})")
    say_fn()


def verify_deletion_lengths(df, modal):
    """Cross-check each deletion oligo's decoded span against its sequence
    length: a span of n bases must shorten the oligo by exactly n.

    This validates the span rule independently of the naming grammar. In APDY
    it holds for all 229 deletion oligos (155 x 1 bp, 13 x 2 bp, 61 x 3 bp).
    """
    ok = bad = 0
    examples = []
    for name, seq in zip(df["oligo_name"], df["sequence"]):
        p = parse_oligo_name(str(name))
        if p is None or p["kind"] != "del":
            continue
        expected = p["pos2"] - p["pos1"] + 1
        observed = modal - len(str(seq))
        if expected == observed:
            ok += 1
        else:
            bad += 1
            if len(examples) < 5:
                examples.append(f"{name}: span {expected} bp, "
                                f"sequence shorter by {observed} bp")
    return ok, bad, examples


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------

def load_variants(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    missing = [c for c in VARIANT_COLS if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: {os.path.basename(path)} is missing expected "
                 f"column(s): {missing}")
    df = df[VARIANT_COLS].copy()
    df["pos_grch38"] = pd.to_numeric(df["pos_grch38"], errors="coerce").astype("Int64")

    # The curated set is fully resolved, but the full variant table also holds
    # quarantined descriptions with no coordinate. Those cannot be joined to
    # anything, so they are dropped here and counted rather than being carried
    # through as silent non-matches.
    unresolved = df[df["pos_grch38"].isna() | df["ref"].isna() | df["alt"].isna()]
    if not unresolved.empty:
        df = df.drop(unresolved.index)
        say(f"Dropped {len(unresolved)} row(s) with no resolved coordinate, "
            f"which cannot be joined: "
            + ", ".join(f"`{k}`" for k in unresolved["variant_key"].head(12)))
        say()

    # Span of the variant on the plus strand, used for position-only matching
    # and for the background-edit overlap check. For a VCF-style allele pair
    # the span is the reference allele.
    df["span_start"] = df["pos_grch38"]
    df["span_end"] = df["pos_grch38"] + df["ref"].fillna("").str.len() - 1

    df["targeton_list"] = df["targetons_design"].fillna("").apply(
        lambda v: [t for t in (x.strip() for x in v.split("|")) if t])
    return df


def add_join_keys(df, refseq):
    """Add `join_key`, the canonically left-aligned key used for matching.

    `variant_key` as delivered is used for reporting and for tracing a row
    back to the curated table, but not for matching: see `left_align()` for
    why some delivered keys are not canonical. `join_key` is rebuilt from the
    `pos_grch38`/`ref`/`alt` columns, which are the authoritative coordinates.
    """
    df = df.copy()
    df["join_key"] = [
        canonical_key(int(p), str(r), str(a), refseq)
        for p, r, a in zip(df["pos_grch38"], df["ref"], df["alt"])]
    df["key_is_canonical"] = df["join_key"] == df["variant_key"]
    df["key_matches_columns"] = df["variant_key"] == [
        L.variant_key(L.CHROM, int(p), str(r), str(a))
        for p, r, a in zip(df["pos_grch38"], df["ref"], df["alt"])]
    return df


def check_key_canonicalisation(df):
    """Report delivered keys that are not canonical, and rows whose key
    disagrees with its own coordinate columns.

    Step 05 now keys on the normalised allele pair and asserts key/coordinate
    agreement itself, so both counts should be zero. This is the independent
    restatement of that at the point of use: if either is non-zero, the variant
    table was built by an older pipeline version than this join expects, and
    the matching would silently miss duplications."""
    non_canonical = df[~df["key_is_canonical"]]
    key_vs_cols = df[~df["key_matches_columns"]]
    if non_canonical.empty and key_vs_cols.empty:
        say("Key canonicalisation: every delivered `variant_key` is "
            "left-aligned and agrees with its own coordinate columns.")
        say()
        return

    if not non_canonical.empty:
        say(f"**{len(non_canonical)} delivered `variant_key`(s) are not "
            f"left-aligned** and are matched on a rebuilt canonical key "
            f"instead. All are duplications, reached through the Variant "
            f"Recoder rather than from a source's own VCF coordinates; see "
            f"`left_align()`. Matching on the delivered key would report "
            f"these as absent from the library.")
        say()
        say("| delivered key | canonical key | HGVS c. | group |")
        say("|---|---|---|---|")
        for _, r in non_canonical.iterrows():
            say(f"| `{r['variant_key']}` | `{r['join_key']}` | "
                f"{r['hgvs_c_mane']} | {r['curation_group']} |")
        say()
    if not key_vs_cols.empty:
        say(f"**{len(key_vs_cols)} row(s) whose `variant_key` disagrees with "
            f"their own `pos_grch38`/`ref`/`alt` columns.** Step 05 takes the "
            f"key from the grouping key, which for a cohort source is that "
            f"source's own VCF coordinate, but takes the coordinate columns "
            f"from the resolver. Both representations of the same variant "
            f"therefore coexist, and grouping on the key string cannot "
            f"collapse a variant reached by both routes. The coordinate "
            f"columns are treated as authoritative here.")
        say()
        for _, r in key_vs_cols.iterrows():
            say(f"- `{r['variant_key']}` versus columns "
                f"chrX:{int(r['pos_grch38'])}:{r['ref']}:{r['alt']} "
                f"({r['hgvs_c_mane']}, sources {r['sources']})")
        say()


def verify_variant_ref_alleles(df, refseq):
    """Confirm each curated reference allele against the plus-strand genome.

    The curated keys are produced by step 03 and audited by step 08, so this is
    a cheap independent restatement at the point of use: if it fails, the
    curated table and this script disagree about the coordinate convention and
    nothing downstream can be trusted.
    """
    bad = []
    for _, r in df.iterrows():
        ref = str(r["ref"] or "")
        observed = refseq.get(int(r["span_start"]), int(r["span_end"]))
        if observed != ref.upper():
            bad.append(f"{r['variant_key']}: table says {ref}, "
                       f"reference has {observed or '<outside cached span>'}")
    return bad


def load_sge_file(path, refseq, background_min_fraction):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    for required in ("oligo_name", "sequence"):
        if required not in d.columns:
            say(f"  WARNING: {os.path.basename(path)} has no '{required}' "
                f"column -- skipping file")
            return None, None
    names = d["oligo_name"].astype(str).tolist()
    seqs = d["sequence"].astype(str).tolist()

    positions = [p["pos1"] for p in (parse_oligo_name(n) for n in names) if p]
    if not positions:
        say(f"  WARNING: no decodable oligo names in "
            f"{os.path.basename(path)} -- skipping file")
        return None, None
    lo, hi = min(positions) - 400, max(positions) + 400

    start, modal, nm = align_amplicon(seqs, refseq, lo, hi)
    fixed, ambiguous = detect_background_edits(
        names, seqs, refseq, start, modal, background_min_fraction)
    del_ok, del_bad, del_examples = verify_deletion_lengths(d, modal)

    meta = {"file": os.path.basename(path), "n_oligos": len(d),
            "coverage": measure_coverage(names),
            "amplicon_start": start, "modal_length": modal,
            "median_mismatch": nm, "background": fixed,
            "background_ambiguous": ambiguous,
            "del_ok": del_ok, "del_bad": del_bad, "del_examples": del_examples}

    bg_positions = set(fixed)
    bg_seq = list(refseq.get(start, start + modal - 1))
    for pos, (_ref, obs, _n, _cov) in fixed.items():
        bg_seq[pos - start] = obs
    bg_seq = "".join(bg_seq)

    keys, reasons, spans, name_keys = [], [], [], []
    disagree = []
    for n, s in zip(names, seqs):
        p = parse_oligo_name(n)
        key, why = oligo_key_from_sequence(s, bg_seq, start, bg_positions, refseq)
        name_key, _name_why = oligo_key_from_name(p, refseq, bg_positions)
        keys.append(key)
        reasons.append(why)
        name_keys.append(name_key)
        if key and name_key and key != name_key:
            if len(disagree) < 10:
                disagree.append(f"{n}: sequence says {key}, name says {name_key}")
        spans.append((p["pos1"], p["pos2"] if p and p["pos2"] is not None
                      else (p["pos1"] if p else None)) if p else (None, None))

    keep = [c for c in SGE_KEEP_COLS if c in d.columns]
    out = d[keep].copy()
    out["oligo_key"] = keys
    out["oligo_key_status"] = reasons
    out["oligo_key_from_name"] = name_keys
    out["oligo_span_start"] = [s[0] for s in spans]
    out["oligo_span_end"] = [s[1] for s in spans]
    out["source_file"] = os.path.basename(path)

    meta["key_status"] = Counter(reasons)
    meta["name_disagreements"] = disagree
    meta["n_name_disagreements"] = sum(
        1 for k, nk in zip(keys, name_keys) if k and nk and k != nk)
    meta["n_keyed_by_sequence_only"] = sum(
        1 for k, nk in zip(keys, name_keys) if k and not nk)
    return out, meta


TARGETON_PREFIX_RE = re.compile(r"^([A-Za-z0-9]+?)_")


def infer_targeton_from_filename(path):
    m = TARGETON_PREFIX_RE.match(os.path.basename(path))
    return m.group(1) if m else None


# --------------------------------------------------------------------------
# Matching
# --------------------------------------------------------------------------

def match(variants, sge_files, refseq, known_targetons, background_min_fraction):
    matches = []
    outcome = {}           # variant_key -> best match type seen
    matched_targeton = defaultdict(list)
    metas = []
    seen_targetons = set()
    oligo_positions: set[int] = set()

    for path in sorted(sge_files):
        targeton = infer_targeton_from_filename(path)
        if targeton is None:
            say(f"  WARNING: cannot infer a targeton from "
                f"'{os.path.basename(path)}' -- skipping")
            continue
        if known_targetons and targeton not in known_targetons:
            say(f"  WARNING: '{targeton}' from "
                f"'{os.path.basename(path)}' is not in the targeton map -- "
                f"skipping, since its variants could not be attributed")
            continue

        sge, meta = load_sge_file(path, refseq, background_min_fraction)
        if sge is None:
            continue
        seen_targetons.add(targeton)
        meta["targeton"] = targeton
        metas.append(meta)
        # Positions spanned by at least one oligo, so that a variant absent
        # from the library can be read as "this site was not built" or "this
        # site was built but not this allele".
        for a, b in zip(sge["oligo_span_start"], sge["oligo_span_end"]):
            if pd.notna(a) and pd.notna(b):
                oligo_positions.update(range(int(a), int(b) + 1))

        here = variants[variants["targeton_list"].apply(lambda ts: targeton in ts)]
        meta["n_curated_in_targeton"] = len(here)
        say(f"  [{targeton}] {meta['file']}: {meta['n_oligos']} oligos, "
            f"amplicon start chrX:{meta['amplicon_start']}, "
            f"{len(meta['background'])} fixed background edit(s), "
            f"{len(here)} curated variant(s) in this design window")
        if here.empty:
            continue

        keyed = sge[sge["oligo_key"].notna()]
        exact = here.merge(keyed, left_on="join_key", right_on="oligo_key",
                           how="inner", suffixes=("", "_sge"))
        if not exact.empty:
            exact = exact.copy()
            exact["match_type"] = "exact_allele"
            exact["matched_targeton"] = targeton
            matches.append(exact)
            for k in exact["variant_key"]:
                outcome[k] = "exact_allele"
                matched_targeton[k].append(targeton)

        # Every curated variant carries a normalised allele, so an oligo with a
        # derived key either matches it or does not; a shared start coordinate
        # adds nothing. The fallback is therefore restricted to oligos whose
        # own allele could not be derived, which are the only ones for which
        # allele agreement is genuinely undecidable.
        got = set(exact["variant_key"]) if not exact.empty else set()
        rest = here[~here["variant_key"].isin(got)]
        undecidable = sge[sge["oligo_key"].isna()]
        if not rest.empty and not undecidable.empty:
            pos_only = rest.merge(undecidable, left_on="span_start",
                                  right_on="oligo_span_start",
                                  how="inner", suffixes=("", "_sge"))
            if not pos_only.empty:
                pos_only = pos_only.copy()
                pos_only["match_type"] = "position_only_verify_manually"
                pos_only["matched_targeton"] = targeton
                matches.append(pos_only)
                for k in pos_only["variant_key"]:
                    if outcome.get(k) != "exact_allele":
                        outcome[k] = "position_only_verify_manually"
                        matched_targeton[k].append(targeton)

    all_matches = (pd.concat(matches, ignore_index=True) if matches
                   else pd.DataFrame())
    return (all_matches, outcome, matched_targeton, metas, seen_targetons,
            oligo_positions)


def build_summary(variants, outcome, matched_targeton, seen_targetons,
                  oligo_positions):
    """One row per curated variant, with its match outcome and, for misses, the
    reason -- distinguishing a library that does not contain the variant from a
    targeton whose results have not been supplied."""
    s = variants.drop(columns=["targeton_list"]).copy()

    def classify(row):
        k = row["variant_key"]
        if k in outcome:
            return outcome[k]
        ts = [t for t in row["targetons_design"].split("|") if t] \
            if isinstance(row["targetons_design"], str) else []
        if not ts:
            return "no_design_targeton"
        if any(t in seen_targetons for t in ts):
            return "no_oligo_in_supplied_file"
        return "no_oligo_targeton_not_supplied"

    s["match_type"] = variants.apply(classify, axis=1)
    s["matched_targeton"] = s["variant_key"].map(
        lambda k: "|".join(sorted(set(matched_targeton.get(k, [])))))
    s["targeton_results_supplied"] = variants["targeton_list"].apply(
        lambda ts: "|".join(t for t in ts if t in seen_targetons))
    s["site_built_in_library"] = [
        any(p in oligo_positions
            for p in range(int(a), int(b) + 1))
        for a, b in zip(variants["span_start"], variants["span_end"])]
    return s


def flag_background_confounding(summary, metas):
    """Mark curated variants whose span overlaps a fixed background edit in the
    targeton that measured them. Their readout reflects the variant on an
    edited background, so it is not comparable with the rest."""
    bg = {}
    for m in metas:
        for pos, info in m["background"].items():
            bg.setdefault(m["targeton"], {})[pos] = info

    flags, detail = [], []
    for _, r in summary.iterrows():
        hit = []
        for t in (r["matched_targeton"].split("|") if r["matched_targeton"] else []):
            for pos, (ref, obs, n, cov) in bg.get(t, {}).items():
                if int(r["span_start"]) <= pos <= int(r["span_end"]):
                    hit.append(f"{t}:chrX:{pos} {ref}>{obs}")
        flags.append(bool(hit))
        detail.append("; ".join(hit))
    summary = summary.copy()
    summary["on_background_edit"] = flags
    summary["background_edit_detail"] = detail
    return summary


def collapse_to_one_row_per_variant(matched, score_col, tier_col):
    """Reduce exact matches to one row per curated variant.

    Distinct designed oligos can normalise to the same allele -- in APDY two
    pairs of single-base deletions in a repeated dinucleotide do, for example
    chrX:53,415,049 CG>C. These are independent measurements of the same
    variant, so the score is averaged, `n_oligos` records how many contributed,
    and the tier is kept only where the oligos agree; disagreement is recorded
    rather than resolved by picking one.
    """
    rows = []
    for key, grp in matched.groupby("variant_key", sort=False):
        base = grp.iloc[0].to_dict()
        scores = pd.to_numeric(grp[score_col], errors="coerce").dropna()
        tiers = sorted(set(grp[tier_col].dropna().astype(str)))
        base["n_oligos"] = len(grp)
        base["oligo_names"] = "|".join(grp["oligo_name"].astype(str))
        base[score_col] = scores.mean() if len(scores) else np.nan
        base[score_col + "_range"] = (
            f"{scores.min():.3f}..{scores.max():.3f}" if len(scores) > 1 else "")
        base[tier_col] = tiers[0] if len(tiers) == 1 else ""
        base["tier_agreement"] = ("single_oligo" if len(grp) == 1
                                  else "agree" if len(tiers) == 1
                                  else "disagree:" + "|".join(tiers))
        rows.append(base)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Summary statistics
# --------------------------------------------------------------------------

def is_depleted_label(x):
    return isinstance(x, str) and "deplet" in x.lower()


def compute_summary_stats(df, score_col, tier_col):
    rows, table = [], []
    for g in GROUPS:
        sub = df[df["curation_group"] == g]
        n_tier = int(sub[tier_col].replace("", np.nan).notna().sum())
        n_depl = int(sub[tier_col].apply(is_depleted_label).sum())
        scores = pd.to_numeric(sub[score_col], errors="coerce").dropna()
        rows.append({
            "group": g, "n": len(sub), "n_with_tier_call": n_tier,
            "n_depleted": n_depl,
            "pct_depleted": (n_depl / n_tier * 100) if n_tier else np.nan,
            f"median_{score_col}": scores.median() if len(scores) else np.nan,
            f"mean_{score_col}": scores.mean() if len(scores) else np.nan,
        })
        table.append([n_depl, n_tier - n_depl])

    stats = pd.DataFrame(rows)
    fisher = None
    if HAVE_SCIPY and all(sum(r) > 0 for r in table):
        odds, p = fisher_exact(table)
        fisher = {"odds_ratio": odds, "p_value": p, "table": table}
    return stats, fisher


# --------------------------------------------------------------------------
# Plot
# --------------------------------------------------------------------------

TIER_COLOR_MAP = {
    "strongly depleting": "#7B241C",
    "depleted": "#C0392B",
    "weakly depleting": "#E67E22",
    "no impact": "#2ECC71",
    "enriched": "#2E86AB",
}
DEFAULT_TIER_COLOR = "#B0B0B0"
GROUP_COLORS = {"CdLS_pathogenic": "#7A5195", "DEE85_pathogenic": "#EF5675"}


def make_plot(df, outpath, score_col, tier_col, class_label,
              stats=None, fisher=None):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))

    ax = axes[0]
    rng = np.random.default_rng(0)
    for i, g in enumerate(GROUPS):
        vals = pd.to_numeric(df.loc[df["curation_group"] == g, score_col],
                             errors="coerce").dropna()
        if vals.empty:
            continue
        ax.scatter(np.full(len(vals), i) + rng.normal(0, 0.05, len(vals)),
                   vals, alpha=0.6, color=GROUP_COLORS[g], label=g)
        ax.scatter([i], [vals.median()], color="black", marker="_",
                   s=800, linewidths=2, zorder=5)
    ref_line = 0.5 if ("post_lof" in score_col.lower()
                       or "prob" in score_col.lower()) else 0
    ax.axhline(ref_line, color="grey", linestyle="--", linewidth=0.8)
    ax.set_xticks(range(len(GROUPS)))
    ax.set_xticklabels([f"CdLS2\n({class_label})", f"DEE85\n({class_label})"])
    ax.set_ylabel(score_col)
    ax.set_title("Depletion score by disease group")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=8)

    ax2 = axes[1]
    preferred = ["strongly depleting", "depleted", "weakly depleting",
                 "no impact", "enriched"]
    present = [t for t in df[tier_col].replace("", np.nan).dropna().unique()]
    order = ([t for t in preferred if t in present]
             + [t for t in present if t not in preferred])
    bottom = np.zeros(len(GROUPS))
    for tier in order:
        fracs = []
        for g in GROUPS:
            sub = df.loc[df["curation_group"] == g, tier_col].replace("", np.nan)
            n = sub.notna().sum()
            fracs.append((sub == tier).sum() / n if n else 0)
        ax2.bar(GROUPS, fracs, bottom=bottom, label=tier,
                color=TIER_COLOR_MAP.get(tier, DEFAULT_TIER_COLOR))
        bottom += np.array(fracs)
    ax2.set_ylabel("Fraction of matched variants")
    ax2.set_xticks(range(len(GROUPS)))
    ax2.set_xticklabels([f"CdLS2\n({class_label})", f"DEE85\n({class_label})"])
    ax2.set_title(f"SGE call breakdown ({tier_col})")
    ax2.legend(fontsize=8)
    ax2.spines[["top", "right"]].set_visible(False)

    fig.suptitle("Curated CdLS2 vs DEE85 variants: SGE depletion behaviour\n"
                 f"({class_label}; exact allele matches only)", fontsize=11)

    if stats is not None and not stats.empty:
        lines = []
        for _, r in stats.iterrows():
            label = "CdLS2" if "CdLS" in r["group"] else "DEE85"
            lines.append(f"{label}: n={int(r['n'])}, "
                         f"depleted={int(r['n_depleted'])}/"
                         f"{int(r['n_with_tier_call'])}"
                         + (f" ({r['pct_depleted']:.0f}%)"
                            if pd.notna(r["pct_depleted"]) else ""))
        if fisher is not None:
            p = fisher["p_value"]
            lines.append("Fisher's exact (depleted vs not): "
                         f"OR={fisher['odds_ratio']:.2g}, "
                         f"p={p:.2e}" if p < 0.001 else
                         "Fisher's exact (depleted vs not): "
                         f"OR={fisher['odds_ratio']:.2g}, p={p:.3f}")
        elif not HAVE_SCIPY:
            lines.append("(scipy not installed: no Fisher's exact test)")
        fig.text(0.5, 0.005, "   |   ".join(lines), ha="center", va="bottom",
                 fontsize=8.5, color="#333333")

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(outpath, dpi=180)
    plt.close()


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

CLASS_FILTERS = {
    "missense": lambda d: d["consequence_class"] == "missense",
    "missense_inframe": lambda d: d["consequence_class"].isin(
        ["missense", "inframe_deletion", "inframe_insertion"]),
    "all": lambda d: pd.Series(True, index=d.index),
}


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--variants",
                   default=os.path.join(ROOT, "data", "processed",
                                        "smc1a_variants_curated.tsv"),
                   help="curated variant table from step 05 "
                        "[default: data/processed/smc1a_variants_curated.tsv]")
    p.add_argument("--sge_glob",
                   default=os.path.join(ROOT, "variant_curation", "sge_results",
                                        "*_all_deseq2_results_condition_*.tsv"),
                   help="glob for per-targeton DESeq2 results; the targeton is "
                        "taken from the filename prefix")
    p.add_argument("--outdir", default=os.path.join(ROOT, "data", "assay_join"))
    p.add_argument("--score_col", default=None,
                   help="continuous depletion score. Default: the first of "
                        + ", ".join(SCORE_PREFERENCE) + " present in the data. "
                        "`anchor_post_lof` is a 0-1 posterior probability, not "
                        "a fold change, and saturates once the classifier is "
                        "confident, so it is the last resort.")
    p.add_argument("--tier_col", default=None,
                   help="categorical call. Default: the first of "
                        + ", ".join(TIER_PREFERENCE) + " present. The binary "
                        "depleted-vs-not statistics treat both 'strongly "
                        "depleting' and 'weakly depleting' as depleted.")
    p.add_argument("--class_filter", default="missense",
                   choices=sorted(CLASS_FILTERS),
                   help="consequence classes admitted to the analysis table. "
                        "'missense' compares the two groups at matched variant "
                        "class, which is the comparison that is not confounded "
                        "by the class difference between the groups; 'all' "
                        "includes every matched variant [default: %(default)s]")
    p.add_argument("--background_min_fraction", type=float, default=0.9,
                   help="fraction of covering oligos that must differ from the "
                        "reference at a position for it to be called a fixed "
                        "background edit [default: %(default)s]")
    p.add_argument("--allow_background_confounded", action="store_true",
                   help="keep variants that overlap a fixed background edit in "
                        "the analysis table. They are excluded by default "
                        "because they were measured on an edited background.")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(DOCS, exist_ok=True)

    refseq = L.ReferenceSequence(os.path.join(REF, "smc1a_chrX_region.fa"))
    tmap = L.TargetonMap.from_tsv(os.path.join(REF, "smc1a_targetons.tsv"))
    known = {t["targeton_id"] for t in tmap.targetons}

    say(f"# Assay join report")
    say()
    say(f"Generated {datetime.date.today().isoformat()}")
    say()
    say(f"Variant table: `{os.path.relpath(args.variants, ROOT)}`")

    variants = add_join_keys(load_variants(args.variants), refseq)
    say(f"Variants loaded: {len(variants)}")
    say()
    check_key_canonicalisation(variants)

    bad_ref = verify_variant_ref_alleles(variants, refseq)
    if bad_ref:
        for b in bad_ref[:10]:
            say(f"  {b}")
        sys.exit(f"ERROR: {len(bad_ref)} curated reference allele(s) disagree "
                 "with the plus-strand reference. The variant table and this "
                 "script do not share a coordinate convention; nothing "
                 "downstream would be meaningful.")
    say("Reference-allele check: all curated `ref` alleles match the "
        "plus-strand GRCh38 reference, so both sides of the join use the same "
        "convention and no complementation is applied.")
    say()

    files = sorted(globmod.glob(args.sge_glob))
    if not files:
        sys.exit(f"ERROR: no files matched --sge_glob '{args.sge_glob}'")
    say(f"## Oligo files ({len(files)})")
    say()
    (matched, outcome, matched_targeton, metas, seen,
     oligo_positions) = match(variants, files, refseq, known,
                              args.background_min_fraction)
    say()

    # ---- per-file empirical checks -------------------------------------
    say("## Empirical checks on each oligo file")
    say()
    for m in metas:
        say(f"### {m['targeton']} -- `{m['file']}`")
        say()
        say(f"- oligos: {m['n_oligos']}; modal length {m['modal_length']} bp")
        say(f"- plus-strand alignment: amplicon starts chrX:"
            f"{m['amplicon_start']}, median {m['median_mismatch']} mismatch(es) "
            f"per oligo against the reference")
        if m["background"]:
            say(f"- fixed background edits ({len(m['background'])}), present in "
                f"essentially every library member and not the variant under "
                f"test:")
            for pos in sorted(m["background"]):
                ref, obs, n, cov = m["background"][pos]
                say(f"    - chrX:{pos} {ref}>{obs} in {n}/{cov} oligos "
                    f"({n / cov:.1%})")
        else:
            say("- fixed background edits: none detected")
        if m["background_ambiguous"]:
            say(f"- **positions differing from the reference in an "
                f"intermediate fraction of oligos "
                f"({len(m['background_ambiguous'])})** -- neither a fixed "
                f"background edit nor absent, which the library design does "
                f"not predict:")
            for pos in sorted(m["background_ambiguous"]):
                ref, obs, n, cov = m["background_ambiguous"][pos]
                say(f"    - chrX:{pos} {ref}>{obs} in {n}/{cov} oligos "
                    f"({n / cov:.1%})")
        say(f"- deletion span cross-check: {m['del_ok']} deletion oligo(s) "
            f"shorten the sequence by exactly the decoded span, {m['del_bad']} "
            f"do not")
        for e in m["del_examples"]:
            say(f"    - {e}")
        undecoded = {k: v for k, v in m["key_status"].items()
                     if not k.startswith("ok")}
        n_keyed = sum(v for k, v in m["key_status"].items()
                      if k.startswith("ok"))
        n_reexp = m["key_status"].get("ok_re_expressed_against_genome", 0)
        say(f"- oligos keyed from their sequence: {n_keyed}/{m['n_oligos']}"
            + (f", of which {n_reexp} re-expressed against the genome because "
               f"the designed substitution sits on a background-edited base"
               if n_reexp else ""))
        if undecoded:
            say("- oligos deliberately not keyed:")
            for k, v in sorted(undecoded.items(), key=lambda kv: -kv[1]):
                say(f"    - {k}: {v}")
        say(f"- name cross-check: {m['n_name_disagreements']} oligo(s) where a "
            f"decodable name implies a different key from the sequence; "
            f"{m['n_keyed_by_sequence_only']} keyed from sequence alone "
            f"because the name does not encode an allele")
        for e in m["name_disagreements"]:
            say(f"    - {e}")
        say(f"- curated variants in this design window: "
            f"{m['n_curated_in_targeton']}")
        say()

    report_library_coverage(metas, say)

    summary = build_summary(variants, outcome, matched_targeton, seen,
                            oligo_positions)
    summary = flag_background_confounding(summary, metas)
    # A substitution landing on a background-edited base is still exactly
    # interpretable (see oligo_key_from_sequence), so it is not confounded even
    # though it overlaps an edit. Only overlaps that could not be resolved that
    # way are withheld from the analysis set.
    resolved_keys = set()
    if not matched.empty and "oligo_key_status" in matched.columns:
        resolved_keys = set(
            matched.loc[matched["oligo_key_status"].astype(str)
                        .str.startswith("ok_re_expressed"), "variant_key"])
    summary["background_edit_resolved"] = summary["variant_key"].isin(resolved_keys)
    summary["background_edit_unresolved"] = (
        summary["on_background_edit"] & ~summary["background_edit_resolved"])

    # ---- resolve score / tier columns ----------------------------------
    available = set(matched.columns) if not matched.empty else set()
    score_col = args.score_col or next(
        (c for c in SCORE_PREFERENCE if c in available), SCORE_PREFERENCE[0])
    tier_col = args.tier_col or next(
        (c for c in TIER_PREFERENCE if c in available), TIER_PREFERENCE[0])
    say(f"Readout columns: score `{score_col}`, tier `{tier_col}`")
    say()

    all_path = os.path.join(args.outdir, "assay_join_all.tsv")
    summary_path = os.path.join(args.outdir, "assay_join_summary.tsv")
    matched.to_csv(all_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)

    say("## Match outcome, one row per curated variant")
    say()
    say("| outcome | variants |")
    say("|---|---|")
    for k, v in summary["match_type"].value_counts().items():
        say(f"| {k} | {v} |")
    say()
    say("By disease group:")
    say()
    ct = pd.crosstab(summary["curation_group"], summary["match_type"])
    say("```")
    say(ct.to_string())
    say("```")
    say()

    screened = summary[summary["match_type"] == "exact_allele"]
    missing_in_file = summary[summary["match_type"] == "no_oligo_in_supplied_file"]
    if not missing_in_file.empty:
        say("Curated variants whose targeton results were read but which are "
            "not in the library. Substitutions are near-saturated across each "
            "design window, but deletions are built at only a minority of "
            "positions outside APDY, and no duplications or insertions are "
            "designed at all, so indels dominate this list:")
        say()
        say("`site built` distinguishes an allele that was not designed at a "
            "site the library does cover, from a site the library does not "
            "reach at all.")
        say()
        say("| variant | HGVS c. | edit | consequence | targeton | site built |")
        say("|---|---|---|---|---|---|")
        for _, r in missing_in_file.iterrows():
            say(f"| `{r['variant_key']}` | {r['hgvs_c_mane']} | "
                f"{r['edit_type']} | {r['consequence_class']} | "
                f"{r['targetons_design']} | "
                f"{'yes' if r['site_built_in_library'] else 'no'} |")
        say()
        say("Edit types among them: "
            + ", ".join(f"{k} {v}" for k, v in
                        missing_in_file["edit_type"].value_counts().items()))
        say()
        report_unmatched_by_class(summary, say)

    confounded = summary[summary["on_background_edit"]]
    if not confounded.empty:
        say(f"**{len(confounded)} matched variant(s) coincide with a fixed "
            f"background edit.** A substitution there is still exactly "
            f"interpretable, because what the cell carries at that position is "
            f"the oligo's own base; an indel there is not, and is withheld "
            f"from the analysis set.")
        say()
        say("| variant | HGVS c. | edit | background edit | resolved |")
        say("|---|---|---|---|---|")
        for _, r in confounded.iterrows():
            say(f"| `{r['variant_key']}` | {r['hgvs_c_mane']} | "
                f"{r['edit_type']} | {r['background_edit_detail']} | "
                f"{'yes, re-expressed against the genome' if r['background_edit_resolved'] else 'no, withheld'} |")
        say()
    else:
        say("No matched curated variant overlaps a fixed background edit.")
        say()

    say(f"Wrote `{os.path.relpath(all_path, ROOT)}` "
        f"({len(matched)} variant x oligo match rows)")
    say(f"Wrote `{os.path.relpath(summary_path, ROOT)}` "
        f"({len(summary)} rows, one per curated variant)")
    say()

    # ---- analysis table -------------------------------------------------
    say("## Analysis set")
    say()
    if matched.empty:
        say("No matches, so no analysis table.")
        write_report(args.outdir)
        return

    exact = matched[matched["match_type"] == "exact_allele"].copy()
    if score_col not in exact.columns or tier_col not in exact.columns:
        say(f"Readout column(s) absent from the matched data "
            f"(score `{score_col}`: {score_col in exact.columns}, "
            f"tier `{tier_col}`: {tier_col in exact.columns}); "
            f"no analysis table.")
        write_report(args.outdir)
        return

    one_per = collapse_to_one_row_per_variant(exact, score_col, tier_col)
    multi = one_per[one_per["n_oligos"] > 1]
    if not multi.empty:
        say(f"{len(multi)} curated variant(s) matched more than one oligo, "
            f"distinct designed edits normalising to the same allele. Scores "
            f"are averaged and the tier kept only where the oligos agree:")
        for _, r in multi.iterrows():
            say(f"- `{r['variant_key']}`: {int(r['n_oligos'])} oligos, "
                f"{r['tier_agreement']}")
        say()

    analysis = one_per[one_per["curation_group"].isin(GROUPS)].copy()
    n_before = len(analysis)
    if not args.allow_background_confounded:
        confounded_keys = set(
            summary.loc[summary["background_edit_unresolved"], "variant_key"])
        analysis = analysis[~analysis["variant_key"].isin(confounded_keys)]
        if len(analysis) != n_before:
            say(f"Excluded {n_before - len(analysis)} background-confounded "
                f"variant(s); --allow_background_confounded keeps them.")
    analysis = analysis[CLASS_FILTERS[args.class_filter](analysis)]
    say(f"Class filter `{args.class_filter}`: {len(analysis)} variant(s) "
        f"admitted.")
    say()

    analysis_path = os.path.join(args.outdir, "assay_join_analysis.tsv")
    analysis.to_csv(analysis_path, sep="\t", index=False)
    say(f"Wrote `{os.path.relpath(analysis_path, ROOT)}`")
    say()

    if analysis.empty:
        say("Analysis table is empty, so no statistics and no figure.")
        write_report(args.outdir)
        return

    say("```")
    say(pd.crosstab(analysis["curation_group"],
                    analysis[tier_col].replace("", "<no call>")).to_string())
    say("```")
    say()

    stats, fisher = compute_summary_stats(analysis, score_col, tier_col)
    stats_path = os.path.join(args.outdir, "assay_join_summary_stats.tsv")
    stats.to_csv(stats_path, sep="\t", index=False)
    say("```")
    say(stats.to_string(index=False))
    say("```")
    say()
    if fisher is not None:
        say(f"Fisher's exact test, depleted versus not, CdLS2 versus DEE85: "
            f"OR={fisher['odds_ratio']:.3g}, p={fisher['p_value']:.3g} "
            f"(2x2 {fisher['table']}).")
    elif not HAVE_SCIPY:
        say("scipy is not installed, so no Fisher's exact test. The counts "
            "above are sufficient to compute one separately.")
    else:
        say("A group had no tier calls, so no Fisher's exact test.")
    n_min = int(stats["n"].min()) if not stats.empty else 0
    if n_min < 10:
        say()
        say(f"The smaller group has n={n_min}. Any test on this is severely "
            f"underpowered and the counts should be read as descriptive until "
            f"more targetons are joined.")
    say()
    say(f"Wrote `{os.path.relpath(stats_path, ROOT)}`")

    if HAVE_MPL:
        plot_path = os.path.join(args.outdir, "assay_join_depletion.png")
        make_plot(analysis, plot_path, score_col, tier_col,
                  args.class_filter, stats, fisher)
        say(f"Wrote `{os.path.relpath(plot_path, ROOT)}`")
    else:
        say("matplotlib is not installed, so the figure was skipped. Every "
            "table above is written regardless; `pip install matplotlib` to "
            "produce it.")

    write_report(args.outdir)


DEFAULT_OUTDIR = os.path.join(ROOT, "data", "assay_join")


def write_report(outdir):
    """Write the run report beside its own outputs, and additionally to
    docs/verification/ for the default configuration only.

    Every configuration used to write the same docs/verification/assay_join.md,
    so running a non-default one -- a different class filter, or the full
    variant table -- silently replaced the canonical report for the delivered
    join. The report now always sits next to the outputs it describes, and the
    canonical copy is written only by the run that produced the delivered
    tables.
    """
    body = "\n".join(report).rstrip() + "\n"
    local = os.path.join(outdir, "assay_join_report.md")
    with open(local, "w", encoding="utf-8") as fh:
        fh.write(body)
    print(f"\nWrote {os.path.relpath(local, ROOT)}")
    if os.path.abspath(outdir) == os.path.abspath(DEFAULT_OUTDIR):
        with open(os.path.join(DOCS, "assay_join.md"), "w",
                  encoding="utf-8") as fh:
            fh.write(body)
        print("Wrote docs/verification/assay_join.md")


if __name__ == "__main__":
    main()
