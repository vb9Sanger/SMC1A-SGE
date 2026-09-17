"""
smc1a_lib.py -- shared utilities for SMC1A variant curation.

Central definitions used by every step of the pipeline:
  * the MANE Select transcript that all variants are recalibrated against
  * an offline c. <-> GRCh38 coordinate mapper built from the Ensembl exon
    structure (used to independently validate every API-derived coordinate)
  * targeton / mutagenised-region assignment
  * a disk-cached HTTP client so that every external API response is stored
    verbatim, making the pipeline reproducible without re-querying live APIs

SMC1A is on the MINUS strand of chrX. Two conventions are used throughout and
are never mixed:
  * "genomic" / "plus-strand": chrX GRCh38 forward-strand coordinates and
    alleles. This is what the SGE library oligo names use, what VCF uses, and
    what every `pos/ref/alt` column in the output tables uses.
  * "transcript" / "c.": HGVS coding coordinates on ENST00000322213.9, which
    run in the direction of decreasing genomic coordinate.

Author: Vanessa Burns
"""

from __future__ import annotations

import hashlib
import re
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request

# --------------------------------------------------------------------------
# Project constants
# --------------------------------------------------------------------------

GENE = "SMC1A"
CHROM = "chrX"
STRAND = -1

# MANE Select v1.4. The SGE library was designed against this transcript, so it
# is the single reference for all coordinate recalibration.
TRANSCRIPT = "ENST00000322213.9"
TRANSCRIPT_NOVER = "ENST00000322213"
PROTEIN = "ENSP00000323421.3"
REFSEQ_MANE = "NM_006306.3"
GENOMIC_ACC = "NC_000023.11"  # chrX, GRCh38
ASSEMBLY = "GRCh38"

# RefSeq transcript versions seen across the source literature. CDS numbering
# (c.1-c.3702) is IDENTICAL in all of these -- verified by
# 01_build_transcript_reference.py -- so c. positions are portable between
# them without an offset. NM_006306.1 carries a 2 bp reference difference at
# c.489_490 (p.163/164) but the same numbering.
EQUIVALENT_CDS_TRANSCRIPTS = [
    "NM_006306.1",
    "NM_006306.2",
    "NM_006306.3",
    "NM_006306.4",
    "ENST00000322213.9",
]
CDS_LENGTH = 3702
PROTEIN_LENGTH = 1233
UTR5_LENGTH = 54

ENSEMBL_REST = "https://rest.ensembl.org"

# Reverse complement / complement
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def complement(base: str) -> str:
    return _COMP[base.upper()]


def revcomp(seq: str) -> str:
    return "".join(_COMP[b] for b in reversed(seq.upper()))


# --------------------------------------------------------------------------
# Cached HTTP client
# --------------------------------------------------------------------------

class CachedHTTP:
    """Disk-cached GET client.

    Every response is written to `cache_dir` keyed by a hash of the URL, so a
    re-run of the pipeline is deterministic and does not depend on the remote
    service still being available or unchanged. Delete the cache directory to
    force a refresh.
    """

    # Ensembl asks that REST clients identify themselves with a contact address.
    # Taken from the SMC1A_CONTACT environment variable rather than hardcoded, so
    # a published copy of this code does not send requests identifying its
    # original author, and anyone running it supplies their own contact.
    DEFAULT_UA = "SMC1A-curation/1.0 (research)"

    def __init__(self, cache_dir: str, delay: float = 0.34,
                 user_agent: str | None = None):
        contact = os.environ.get("SMC1A_CONTACT", "").strip()
        user_agent = user_agent or (
            f"SMC1A-curation/1.0 (research; {contact})" if contact
            else self.DEFAULT_UA)
        self.cache_dir = cache_dir
        self.delay = delay
        self.user_agent = user_agent
        os.makedirs(cache_dir, exist_ok=True)
        self._last = 0.0
        self.n_hits = 0
        self.n_misses = 0

    def _path(self, url: str) -> str:
        return os.path.join(self.cache_dir, hashlib.sha256(url.encode()).hexdigest() + ".json")

    def get(self, url: str, tries: int = 4):
        """Return {'ok':bool, 'status':int, 'body':str, 'url':str}."""
        p = self._path(url)
        if os.path.exists(p):
            self.n_hits += 1
            with open(p, encoding="utf-8") as fh:
                return json.load(fh)

        rec = None
        for attempt in range(tries):
            gap = self.delay - (time.time() - self._last)
            if gap > 0:
                time.sleep(gap)
            req = urllib.request.Request(url, headers={
                "User-Agent": self.user_agent, "Accept": "application/json"})
            try:
                with urllib.request.urlopen(req, timeout=90) as resp:
                    rec = {"ok": True, "status": resp.status,
                           "body": resp.read().decode("utf-8"), "url": url}
                self._last = time.time()
                break
            except urllib.error.HTTPError as e:
                body = e.read().decode("utf-8", "replace")
                self._last = time.time()
                # 4xx are deterministic (bad HGVS) -> cache them; 5xx/429 retry
                if e.code in (429, 500, 502, 503, 504) and attempt < tries - 1:
                    time.sleep(2 ** attempt)
                    continue
                rec = {"ok": False, "status": e.code, "body": body, "url": url}
                break
            except Exception as e:  # network-level; retry then give up uncached
                self._last = time.time()
                if attempt < tries - 1:
                    time.sleep(2 ** attempt)
                    continue
                return {"ok": False, "status": -1, "body": f"{type(e).__name__}: {e}", "url": url}

        self.n_misses += 1
        # NEVER cache a transient server-side failure: a momentary outage would
        # otherwise become a permanent "unresolvable variant". 4xx responses are
        # deterministic (genuinely bad HGVS) and are cached.
        if rec["ok"] or rec["status"] not in (429, 500, 502, 503, 504, -1):
            with open(p, "w", encoding="utf-8") as fh:
                json.dump(rec, fh)
        else:
            self.n_transient_failures = getattr(self, "n_transient_failures", 0) + 1
        return rec

    def get_json(self, url: str):
        rec = self.get(url)
        if not rec["ok"]:
            return None, rec
        try:
            return json.loads(rec["body"]), rec
        except json.JSONDecodeError:
            return None, rec


# --------------------------------------------------------------------------
# Offline transcript coordinate model
# --------------------------------------------------------------------------

class TranscriptModel:
    """Offline c. <-> genomic mapper for a minus-strand transcript.

    Built from the exon table produced by 01_build_transcript_reference.py.
    Exists to provide a second, independent implementation that every
    API-derived coordinate is checked against, so that a silent API
    mis-mapping cannot propagate into the curated table.
    """

    def __init__(self, exons, cds_t_start: int, cds_t_end: int):
        # exons: list of dicts with exon_number, g_start, g_end, t_start, t_end
        # ordered 5'->3' along the transcript (decreasing genomic coordinate)
        self.exons = sorted(exons, key=lambda e: e["t_start"])
        self.cds_t_start = cds_t_start
        self.cds_t_end = cds_t_end
        self.t_len = self.exons[-1]["t_end"]

    # ---- construction -------------------------------------------------
    @classmethod
    def from_json(cls, path: str) -> "TranscriptModel":
        d = json.load(open(path, encoding="utf-8"))
        return cls(d["exons"], d["cds_t_start"], d["cds_t_end"])

    # ---- genomic <-> transcript --------------------------------------
    def g_to_t(self, g: int):
        """Genomic position -> transcript position (None if intronic/outside)."""
        for e in self.exons:
            if e["g_start"] <= g <= e["g_end"]:
                return e["t_start"] + (e["g_end"] - g)  # minus strand
        return None

    def t_to_g(self, t: int):
        for e in self.exons:
            if e["t_start"] <= t <= e["t_end"]:
                return e["g_end"] - (t - e["t_start"])
        return None

    # ---- transcript <-> c. -------------------------------------------
    def t_to_c(self, t: int) -> str:
        if t < self.cds_t_start:
            return f"-{self.cds_t_start - t}"
        if t > self.cds_t_end:
            return f"*{t - self.cds_t_end}"
        return str(t - self.cds_t_start + 1)

    def c_to_t(self, c: str):
        c = c.strip()
        if c.startswith("-"):
            return self.cds_t_start - int(c[1:])
        if c.startswith("*"):
            return self.cds_t_end + int(c[1:])
        return self.cds_t_start + int(c) - 1

    # ---- full c. position parsing (incl. intronic offsets) -----------
    def c_pos_to_g(self, cpos: str):
        """Parse one HGVS c. position, e.g. '1951', '616-2', '109+108', '-19',
        '*23', and return (genomic_position, region) where region is
        'exonic' or 'intronic'. Returns (None, reason) if unmappable.

        Minus-strand aware: a '+n' offset moves to LOWER genomic coordinate,
        a '-n' offset to HIGHER.
        """
        cpos = cpos.strip()
        # split anchor from intronic offset. Careful: leading '-' is a UTR
        # marker, not an offset, so scan from position 1.
        anchor, offset = cpos, 0
        for i in range(1, len(cpos)):
            if cpos[i] in "+-":
                anchor, offset = cpos[:i], int(cpos[i:])
                break
        try:
            t = self.c_to_t(anchor)
        except ValueError:
            return None, f"unparseable c. position {cpos!r}"
        if t is None or not (1 <= t <= self.t_len):
            return None, f"c. position {cpos!r} outside transcript"
        g_anchor = self.t_to_g(t)
        if g_anchor is None:
            return None, f"c. position {cpos!r} not on an exon"
        if offset == 0:
            return g_anchor, "exonic"

        # Validate that the anchor really is an exon boundary in the direction
        # implied by the offset -- this is the check that catches invalid HGVS
        # such as 'c.587-2' (c.587 is mid-exon), which some tools silently
        # coerce rather than reject.
        ex = self._exon_of_t(t)
        if offset > 0 and t != ex["t_end"]:
            return None, (f"c.{cpos}: +offset requires an exon 3' end; "
                          f"c.{anchor} is not the last base of exon {ex['exon_number']}")
        if offset < 0 and t != ex["t_start"]:
            return None, (f"c.{cpos}: -offset requires an exon 5' start; "
                          f"c.{anchor} is not the first base of exon {ex['exon_number']}")
        # minus strand: transcript 3' == lower genomic coordinate
        return g_anchor - offset, "intronic"

    def _exon_of_t(self, t: int):
        for e in self.exons:
            if e["t_start"] <= t <= e["t_end"]:
                return e
        return None

    def exon_of_g(self, g: int):
        """Return exon number if `g` is exonic, else None."""
        for e in self.exons:
            if e["g_start"] <= g <= e["g_end"]:
                return e["exon_number"]
        return None

    def intron_of_g(self, g: int):
        """Return intron number (n = intron between exon n and n+1) or None."""
        for a, b in zip(self.exons, self.exons[1:]):
            # minus strand: exon a is at HIGHER coordinates than exon b
            lo, hi = b["g_end"] + 1, a["g_start"] - 1
            if lo <= g <= hi:
                return a["exon_number"]
        return None

    def protein_position_of_c(self, c: int):
        """c. coding position -> 1-based codon number."""
        if not (1 <= c <= CDS_LENGTH):
            return None
        return (c - 1) // 3 + 1


# --------------------------------------------------------------------------
# Targeton assignment
# --------------------------------------------------------------------------

class TargetonMap:
    """Assign variants to SGE targetons.

    THREE nested regions per targeton matter and are kept separate, because
    conflating them gives the wrong answer to "could this variant be in my
    library?":

      * amplicon      (ref_start..ref_end)      the sequenced amplicon
      * design window ("Output targeton coordinates")
                      the window within which library variants were
                      systematically designed. This is WIDER than the exonic
                      core: it extends into the flanking intron so that splice
                      donor/acceptor and polypyrimidine-tract variants are
                      covered.
      * exon core     (r2_start..r2_end)        the exonic portion only

    Verified empirically against the APDY results: oligos span
    chrX:53,414,941-53,415,183, whereas r2 is only 53,414,981-53,415,169. The
    122 oligos outside r2 are exactly the intronic, splice and
    polypyrimidine-tract variants. Assigning on r2 would therefore wrongly
    declare every splice variant absent from the library.

    A handful of library members lie outside even the design window: known
    variants seeded from ClinVar and gnomAD (recognisable by the `_clinvar` /
    `_gnomAD_*` suffix in the oligo name). For that reason `in_sge_library` is
    ALWAYS determined empirically by joining to the real oligo list
    (06_join_assay.py). These region flags are the a-priori expectation and
    serve to explain misses, never to substitute for the join.
    """

    REGION_KEYS = {
        "design": ("design_start", "design_end"),
        "exon_core": ("mut_start", "mut_end"),
        "amplicon": ("amplicon_start", "amplicon_end"),
    }

    def __init__(self, targetons):
        self.targetons = targetons  # list of dicts

    @classmethod
    def from_tsv(cls, path: str) -> "TargetonMap":
        rows = []
        with open(path, encoding="utf-8-sig") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                if not line.strip():
                    continue
                d = dict(zip(hdr, line.rstrip("\n").split("\t")))
                for k in ("exon", "amplicon_start", "amplicon_end",
                          "mut_start", "mut_end", "design_start", "design_end"):
                    if d.get(k) not in (None, "", "NA"):
                        d[k] = int(d[k])
                rows.append(d)
        return cls(rows)

    def assign(self, g_start: int, g_end: int):
        """Return (design_hits, exon_core_hits, amplicon_hits): lists of
        targeton IDs whose respective region overlaps the variant span.

        `design_hits` is the operative one for "expected in the library".
        """
        out = {k: [] for k in self.REGION_KEYS}
        for t in self.targetons:
            for region, (ks, ke) in self.REGION_KEYS.items():
                if ks not in t or ke not in t or t[ks] == "":
                    continue
                if g_start <= t[ke] and g_end >= t[ks]:
                    out[region].append(t["targeton_id"])
        return out["design"], out["exon_core"], out["amplicon"]


# --------------------------------------------------------------------------
# Variant key normalisation
# --------------------------------------------------------------------------

class ReferenceSequence:
    """Plus-strand GRCh38 reference for the SMC1A locus, read from the FASTA
    cached by 01_build_transcript_reference.py.

    Exists so that every coordinate the pipeline produces can be checked
    against the actual reference base without a network call.
    """

    def __init__(self, path: str):
        with open(path, encoding="utf-8") as fh:
            hdr = fh.readline().strip()
            span = hdr.split()[0].lstrip(">").split(":")[1]
            self.start, self.end = (int(x) for x in span.split("-"))
            self.seq = "".join(l.strip() for l in fh).upper()
        if len(self.seq) != self.end - self.start + 1:
            raise ValueError(f"FASTA length {len(self.seq)} != span "
                             f"{self.end - self.start + 1}")

    def get(self, start: int, end: int) -> str:
        """Inclusive 1-based plus-strand slice; '' if outside the cached span."""
        if start < self.start or end > self.end or end < start:
            return ""
        return self.seq[start - self.start:end - self.start + 1]

    def base(self, pos: int) -> str:
        return self.get(pos, pos)


def normalise_vcf(pos: int, ref: str, alt: str, refseq: ReferenceSequence):
    """Left-align and parsimoniously trim a plus-strand VCF-style allele pair,
    so the same biological variant from different sources collapses to one key.

    Implements the standard vt/bcftools normalisation:
      1. while ref and alt both end in the same base, drop the last base of
         each -- extending one base leftwards first whenever an allele would
         otherwise be trimmed away entirely;
      2. while ref and alt both start with the same base and both are >1 long,
         drop the first base of each and advance pos;
      3. while either allele is empty, extend one base to the left.

    Step 1 previously required both alleles to be longer than one base, which
    made the function a no-op on an *anchored* indel: `T>TT` was returned
    unchanged instead of being shifted to the 5'-most equivalent position.
    Duplications arrive in exactly that form, so two representations of one
    duplication -- `chrX:53415122:T:TT` and `chrX:53415121:G:GT` for c.157dup
    -- both survived, defeating the function's whole purpose. The left
    extension inside step 1 is what allows an anchored indel to shift.
    """
    ref, alt = ref.upper(), alt.upper()

    while ref and alt and ref[-1] == alt[-1]:
        if len(ref) == 1 or len(alt) == 1:
            prev = refseq.base(pos - 1)
            if not prev:
                break
            ref, alt, pos = prev + ref, prev + alt, pos - 1
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0] and len(ref) > 1 and len(alt) > 1:
        ref, alt, pos = ref[1:], alt[1:], pos + 1

    # a pure indel represented without an anchor base: shift left to the most
    # 5' equivalent position, then add the anchor
    if not ref or not alt:
        while True:
            prev = refseq.base(pos - 1)
            if not prev:
                break
            # can we shift one base left and stay equivalent?
            if ref and ref[-1] == prev:
                ref, pos = prev + ref[:-1], pos - 1
            elif alt and alt[-1] == prev:
                alt, pos = prev + alt[:-1], pos - 1
            else:
                break
        anchor = refseq.base(pos - 1)
        if anchor:
            ref, alt, pos = anchor + ref, anchor + alt, pos - 1
    return pos, ref, alt


def variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Join key, in the same form as the SGE oligo naming.

    Takes the allele pair as given. Use `canonical_variant_key` wherever keys
    from different sources have to agree.
    """
    return f"{chrom}:{pos}:{ref.upper()}:{alt.upper()}"


def canonical_variant_key(chrom: str, pos: int, ref: str, alt: str,
                          refseq: ReferenceSequence) -> str:
    """Join key built from the normalised allele pair.

    The single key form that sources with different coordinate conventions can
    be compared on: a variant supplied as VCF coordinates by one source and
    derived from HGVS by another reduces to the same string here. Always
    prefer this over `variant_key` for identity, grouping and joining.
    """
    return variant_key(chrom, *normalise_vcf(pos, ref, alt, refseq))


# --------------------------------------------------------------------------
# HGVS c. parsing
# --------------------------------------------------------------------------

import re as _re

_HGVS_C = _re.compile(
    r"^(?:(?P<ref>[A-Za-z0-9_.]+):)?"          # optional reference accession
    r"c\."                                      # coding coordinate prefix
    r"(?P<p1>[-*]?\d+(?:[+-]\d+)?)"             # first position
    r"(?:_(?P<p2>[-*]?\d+(?:[+-]\d+)?))?"       # optional second position
    r"(?P<edit>.*)$"
)

_EDIT = [
    ("sub",    _re.compile(r"^(?P<r>[ACGT])>(?P<a>[ACGT])$", _re.I)),
    ("delins", _re.compile(r"^del(?P<r>[ACGT]*)ins(?P<a>[ACGT]+)$", _re.I)),
    ("del",    _re.compile(r"^del(?P<r>[ACGT]*)$", _re.I)),
    ("dup",    _re.compile(r"^dup(?P<r>[ACGT]*)$", _re.I)),
    ("ins",    _re.compile(r"^ins(?P<a>[ACGT]+)$", _re.I)),
    ("inv",    _re.compile(r"^inv(?P<r>[ACGT]*)$", _re.I)),
    ("ident",  _re.compile(r"^=$")),
]


def parse_hgvs_c(desc: str):
    """Parse an HGVS c. description into its components.

    Returns a dict with keys: accession, pos1, pos2, edit_type, ref_c, alt_c,
    normalised (the canonical 'c....' string with the accession stripped).
    Returns None if the description cannot be parsed, which is treated as a
    curation failure rather than being silently guessed at.
    """
    desc = str(desc).strip().replace(" ", "")
    # tolerate the common 'NM_006306.2(SMC1A):c.123A>G' and 'SMC1A:c.' forms
    desc = _re.sub(r"\([^)]*\)", "", desc)
    m = _HGVS_C.match(desc)
    if not m:
        return None
    edit = m.group("edit")
    etype = ref_c = alt_c = None
    for name, rx in _EDIT:
        em = rx.match(edit)
        if em:
            etype = name
            gd = em.groupdict()
            ref_c = (gd.get("r") or "").upper()
            alt_c = (gd.get("a") or "").upper()
            break
    if etype is None:
        return None
    p1, p2 = m.group("p1"), m.group("p2")
    return {
        "accession": m.group("ref") or "",
        "pos1": p1,
        "pos2": p2 or p1,
        "edit_type": etype,
        "ref_c": ref_c,
        "alt_c": alt_c,
        "normalised": f"c.{p1}" + (f"_{p2}" if p2 else "") + edit,
    }


# --------------------------------------------------------------------------
# Codon / amino-acid utilities (used to independently verify protein changes)
# --------------------------------------------------------------------------

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*", "Sec": "U", "Xaa": "X",
}
AA1 = {v: k for k, v in AA3.items()}


def aa3_to_1(s: str) -> str:
    """'p.(Val651Met)' / 'p.Val651Met' -> 'V651M'. Returns '' if unparseable."""
    s = str(s).strip()
    s = _re.sub(r"^[A-Za-z0-9_.]+:", "", s)
    s = s.replace("p.", "").strip("()")
    m = _re.match(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|\*)", s)
    if not m:
        return ""
    a, n, b = m.groups()
    if a not in AA3:
        return ""
    if b == "=":
        return f"{AA3[a]}{n}{AA3[a]}"
    return f"{AA3[a]}{n}{AA3.get(b, b)}"


# --------------------------------------------------------------------------
# Consequence harmonisation
# --------------------------------------------------------------------------

# Map VEP SO terms onto the coarse classes used for analysis, chosen to line up
# with the `consequence` column of the SGE DESeq2 output.
CONSEQUENCE_CLASS = {
    "missense_variant": "missense",
    "synonymous_variant": "synonymous",
    "stop_gained": "nonsense",
    "stop_lost": "stop_lost",
    "start_lost": "start_lost",
    "frameshift_variant": "frameshift",
    "inframe_deletion": "inframe_deletion",
    "inframe_insertion": "inframe_insertion",
    "protein_altering_variant": "protein_altering",
    "splice_acceptor_variant": "splice_acceptor",
    "splice_donor_variant": "splice_donor",
    "splice_donor_5th_base_variant": "splice_region",
    "splice_region_variant": "splice_region",
    "splice_polypyrimidine_tract_variant": "splice_region",
    "splice_donor_region_variant": "splice_region",
    "intron_variant": "intronic",
    "5_prime_UTR_variant": "utr5",
    "3_prime_UTR_variant": "utr3",
    "upstream_gene_variant": "upstream",
    "downstream_gene_variant": "downstream",
    "coding_sequence_variant": "coding_unknown",
    "stop_retained_variant": "synonymous",
    "transcript_ablation": "cnv_whole_gene",
    "feature_truncation": "cnv_partial",
    "feature_elongation": "cnv_partial",
    "transcript_amplification": "cnv_gain",
}

# Predicted loss-of-function classes (the classes SGE depletion is expected to
# report on directly).
PTV_CLASSES = {"nonsense", "frameshift", "splice_acceptor", "splice_donor",
               "start_lost", "cnv_whole_gene", "cnv_partial"}


def consequence_class(so_terms) -> str:
    """Most severe harmonised class from a list of VEP SO terms."""
    order = ["cnv_whole_gene", "cnv_partial", "frameshift", "nonsense",
             "splice_acceptor", "splice_donor", "start_lost", "stop_lost",
             "inframe_deletion", "inframe_insertion", "protein_altering",
             "missense", "splice_region", "synonymous", "coding_unknown",
             "utr5", "utr3", "intronic", "upstream", "downstream", "cnv_gain"]
        # unknown terms fall through to 'other'
    classes = {CONSEQUENCE_CLASS.get(t, "other") for t in so_terms}
    for o in order:
        if o in classes:
            return o
    return "other"


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------

def blank(v) -> bool:
    """LOVD and hand-curated sheets use several tokens for 'no data'.

    Note that '?' is deliberately NOT treated as blank. In LOVD it is a
    meaningful disease symbol ('unclassified / mixed'), and conflating it with
    a missing value would erase the distinction between "explicitly
    unclassified" and "no disease recorded" -- two different curation states.
    The null token in the export is the empty string; '-' is the null used by
    the web view.
    """
    return v is None or str(v).strip() in ("", "-", "NA", "na", "N/A", ".", "none", "None")


def clean(v) -> str:
    return "" if blank(v) else str(v).strip()


def join_unique(values, sep="|"):
    """Order-preserving unique join, used for the aggregated list columns."""
    out = []
    for v in values:
        v = clean(v)
        if v and v not in out:
            out.append(v)
    return sep.join(out)
AA1_TO_AA3 = {v: k for k, v in AA3.items()}

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


# ---------------------------------------------------------------------------
# Back-derivation of a cDNA description from a protein description
#
# Some sources report only a protein change (legacy `R496H` style, or a
# modelling target such as Di Nardo's Arg96Cys) with no cDNA anywhere in the
# paper. Joining such a source to a variant by protein-change match is unsafe in
# general, because the same amino acid substitution can arise from different
# nucleotide changes -- which is why the standing rule forbids it.
#
# There is however a decidable subset. If the reference codon admits EXACTLY ONE
# single-nucleotide change producing the stated amino acid, the cDNA description
# is uniquely determined and no inference is involved. That condition is what
# this function tests; it returns a result only when the route is unique.
#
# Validated across the registry against variants whose cDNA was independently
# extracted: 99/111 missense and 22/27 nonsense have a unique route, and in
# 121/121 of those the back-derived description reproduced the extracted cDNA
# exactly. The ~12% ambiguous cases are refused, not guessed.
#
# Frameshifts are refused outright: `p.(Gln361GlyfsTer42)` fixes the first
# affected codon but says nothing about which bases moved or where, so many
# distinct indels share one annotation.
# ---------------------------------------------------------------------------
_P_SIMPLE = re.compile(
    r"^p\.\(?([A-Z][a-z]{2})(\d+)(Ter|\*|[A-Z][a-z]{2})\)?$")


def back_derive_c_from_p(protein: str, tm, refseq):
    """Derive a cDNA description from a protein change, or refuse.

    Returns a dict with `status` one of:
      unique             -- exactly one single-nt route; `hgvs_c` is determined
      ambiguous          -- more than one route; refused
      no_single_nt_route -- needs >1 nucleotide change; refused
      not_substitution   -- frameshift/indel/unparseable; refused
      reference_mismatch -- stated reference residue is not what MANE has
    """
    out = {"status": "not_substitution", "hgvs_c": "", "routes": [],
           "codon": None, "ref_codon": "", "note": ""}
    m = _P_SIMPLE.match((protein or "").strip())
    if not m:
        out["note"] = "not a simple single-residue substitution"
        return out
    ref3, num, alt3 = m.group(1), int(m.group(2)), m.group(3)
    if not 1 <= num <= 1233:
        out["status"] = "not_substitution"
        out["note"] = f"codon {num} outside the 1233-residue protein"
        return out
    out["codon"] = num

    def base_at(cpos):
        g = tm.c_pos_to_g(str(cpos))
        if isinstance(g, tuple):
            g = g[0]
        return COMPLEMENT[refseq.base(int(g)).upper()]

    ref_codon = "".join(base_at((num - 1) * 3 + 1 + i) for i in range(3))
    out["ref_codon"] = ref_codon
    observed = CODON_TABLE.get(ref_codon, "?")
    if AA1_TO_AA3.get(observed) != ref3:
        out["status"] = "reference_mismatch"
        out["note"] = (f"source states {ref3}{num} but MANE codon {num} "
                       f"({ref_codon}) encodes {AA1_TO_AA3.get(observed, observed)}")
        return out

    target = "*" if alt3 in ("Ter", "*") else AA3.get(alt3)
    if not target:
        out["note"] = f"unrecognised target residue {alt3}"
        return out

    routes = []
    for i in range(3):
        for b in "ACGT":
            if b == ref_codon[i]:
                continue
            mut = ref_codon[:i] + b + ref_codon[i + 1:]
            if CODON_TABLE.get(mut) == target:
                routes.append(f"c.{(num - 1) * 3 + 1 + i}{ref_codon[i]}>{b}")
    out["routes"] = routes
    if len(routes) == 1:
        out["status"] = "unique"
        out["hgvs_c"] = routes[0]
        out["note"] = (f"codon {num} reference {ref_codon} admits exactly one "
                       f"single-nucleotide change producing {alt3}, so the cDNA "
                       f"description is uniquely determined, not inferred")
    elif routes:
        out["status"] = "ambiguous"
        out["note"] = (f"{len(routes)} single-nucleotide routes give {alt3} "
                       f"({', '.join(routes)}); refused")
    else:
        out["status"] = "no_single_nt_route"
        out["note"] = (f"no single-nucleotide change of codon {num} ({ref_codon}) "
                       f"gives {alt3}; refused")
    return out
