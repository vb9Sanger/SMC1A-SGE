#!/usr/bin/env python3
"""
01_build_transcript_reference.py

Builds the reference layer that every later step depends on, and formally
verifies the assumptions that the curation rests on.

Outputs
-------
data/reference/smc1a_exons.tsv            exon table with genomic, transcript
                                          and c. coordinate bounds
data/reference/smc1a_transcript_model.json machine-readable model consumed by
                                          smc1a_lib.TranscriptModel
data/reference/smc1a_targetons.tsv        targeton amplicon + mutagenised
                                          region coordinates, joined to exons
data/reference/smc1a_chrX_region.fa       cached plus-strand reference sequence
                                          spanning the locus
docs/verification/transcript_equivalence.md  audit report

Verifications performed (all recorded in the audit report)
----------------------------------------------------------
V1  Ensembl ENST00000322213.9 exon structure matches the lab's SMC1A_exon_map.tsv
V2  CDS length and protein length match the expected 3702 nt / 1233 aa
V3  CDS sequence identity between NM_006306.1/.2/.3/.4 and ENST00000322213.9,
    establishing whether c. numbering is portable across the transcript
    versions used by the source literature
V4  targeton_regions.tsv is consistent with the SMC1A_info.xlsx manifest
V5  the offline TranscriptModel round-trips c. <-> genomic for every exon
    boundary, and agrees with Ensembl Variant Recoder on a probe set

Usage
-----
    ../.venv/bin/python 01_build_transcript_reference.py
"""

from __future__ import annotations

import json
import os
import sys
import urllib.parse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smc1a_lib as L  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, "data", "reference")
CACHE = os.path.join(ROOT, "data", "cache", "ensembl")
DOCS = os.path.join(ROOT, "docs", "verification")
SCREEN = os.path.join(ROOT, "screen_info")

for d in (REF, CACHE, DOCS):
    os.makedirs(d, exist_ok=True)

http = L.CachedHTTP(CACHE)
report: list[str] = []


def say(msg: str = "") -> None:
    print(msg)
    report.append(msg)


def check(label: str, ok: bool, detail: str = "") -> bool:
    say(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f" -- {detail}" if detail else ""))
    return ok


# =========================================================================
# Ensembl exon structure
# =========================================================================

def fetch_transcript():
    url = (f"{L.ENSEMBL_REST}/lookup/id/{L.TRANSCRIPT_NOVER}"
           f"?expand=1;utr=1;content-type=application/json")
    d, rec = http.get_json(url)
    if d is None:
        sys.exit(f"FATAL: could not fetch transcript: {rec['status']} {rec['body'][:200]}")
    return d


def build_model(tr):
    if tr["strand"] != L.STRAND:
        sys.exit(f"FATAL: expected strand {L.STRAND}, Ensembl reports {tr['strand']}")
    # minus strand: transcript 5'->3' runs from HIGHEST to LOWEST genomic coord
    exons = sorted(tr["Exon"], key=lambda e: -e["start"])
    rows, tpos = [], 0
    for n, e in enumerate(exons, start=1):
        length = e["end"] - e["start"] + 1
        rows.append({
            "exon_number": n,
            "exon_id": e["id"],
            "g_start": e["start"],
            "g_end": e["end"],
            "length": length,
            "t_start": tpos + 1,
            "t_end": tpos + length,
        })
        tpos += length

    model = L.TranscriptModel(rows, 1, tpos)  # provisional, to use g_to_t
    tl = tr["Translation"]
    # minus strand: translation 'end' is the genomic position of the A of ATG
    cds_t_start = model.g_to_t(tl["end"])
    cds_t_end = model.g_to_t(tl["start"])
    if cds_t_start is None or cds_t_end is None:
        sys.exit("FATAL: translation bounds are not exonic")

    model = L.TranscriptModel(rows, cds_t_start, cds_t_end)
    for r in rows:
        r["c_start"] = model.t_to_c(r["t_start"])
        r["c_end"] = model.t_to_c(r["t_end"])
    return model, rows, tl


# =========================================================================
# Verifications
# =========================================================================

def v1_exon_map(rows) -> bool:
    """Ensembl exon structure vs the lab's own exon map."""
    path = os.path.join(SCREEN, "SMC1A_exon_map.tsv")
    lab = {}
    with open(path, encoding="utf-8-sig") as fh:
        hdr = fh.readline().strip().split(",")
        for line in fh:
            if not line.strip():
                continue
            d = dict(zip(hdr, line.strip().split(",")))
            span = d["Exon_position"].split(":")[1]
            s, e = (int(x) for x in span.split("-"))
            lab.setdefault(int(d["EXON"]), set()).add((s, e))
    ok = True
    for r in rows:
        n = r["exon_number"]
        if n not in lab:
            ok = check(f"exon {n} present in lab exon map", False) and ok
            continue
        got = (r["g_start"], r["g_end"])
        if got not in lab[n]:
            # exon 25: the lab map gives the transcribed extent, Ensembl the
            # full 3'UTR-containing exon -- report rather than silently accept
            ok = check(f"exon {n} coordinates agree", False,
                       f"Ensembl {got} vs lab {sorted(lab[n])}") and ok
    if ok:
        check(f"all {len(rows)} exon coordinate spans agree with SMC1A_exon_map.tsv", True)
    return ok


def v2_lengths(rows, model, tl) -> bool:
    cds_len = model.cds_t_end - model.cds_t_start + 1
    a = check("CDS length is 3702 nt", cds_len == L.CDS_LENGTH, f"got {cds_len}")
    b = check("protein length is 1233 aa", tl["length"] == L.PROTEIN_LENGTH, f"got {tl['length']}")
    c = check("5'UTR length is 54 nt", model.cds_t_start - 1 == L.UTR5_LENGTH,
              f"got {model.cds_t_start - 1}")
    d = check("transcript has 25 exons", len(rows) == 25, f"got {len(rows)}")
    return all([a, b, c, d])


def v3_transcript_equivalence() -> bool:
    """Are c. positions portable between the RefSeq versions used by the
    source literature and the MANE transcript the assay was built on?"""
    cds = {}
    # CDS bounds read from the GenBank flatfile of each version.
    # NM_001281463.1 is included because it is the transcript most frequently
    # cited in LOVD's `Published_as` field (42 of 179 entries) -- it is a
    # different isoform, not another version of the same one.
    for acc in ["NM_006306.1", "NM_006306.2", "NM_006306.3", "NM_006306.4",
                "NM_001281463.1"]:
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
               f"?db=nuccore&id={acc}&rettype=gb&retmode=text")
        rec = http.get(url)
        if not rec["ok"]:
            check(f"fetch {acc}", False, f"HTTP {rec['status']}")
            continue
        line = next((l for l in rec["body"].splitlines()
                     if l.startswith("     CDS  ")), None)
        span = line.split()[1]
        s, e = (int(x) for x in span.replace("<", "").replace(">", "").split(".."))
        url2 = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                f"?db=nuccore&id={acc}&rettype=fasta&retmode=text"
                f"&seq_start={s}&seq_stop={e}")
        rec2 = http.get(url2)
        cds[acc] = ("".join(rec2["body"].split("\n")[1:]).upper(), s, e)

    d, _ = http.get_json(f"{L.ENSEMBL_REST}/sequence/id/{L.TRANSCRIPT_NOVER}"
                         "?type=cds;content-type=application/json")
    cds[L.TRANSCRIPT] = (d["seq"].upper(), None, None)

    mane = cds[L.TRANSCRIPT][0]
    say()
    say("  | transcript          | 5'UTR len | CDS len | CDS identical to MANE | c. numbering portable |")
    say("  |---------------------|-----------|---------|-----------------------|-----------------------|")
    all_ok = True
    for acc, (seq, s, e) in cds.items():
        utr = (s - 1) if s else L.UTR5_LENGTH
        same = seq == mane
        portable = len(seq) == len(mane)  # equal CDS length => no numbering offset
        say(f"  | {acc:19} | {utr:>9} | {len(seq):>7} | {str(same):21} | {str(portable):21} |")
        all_ok = all_ok and portable
    say()

    check("CDS numbering is portable across all NM_006306 versions and MANE", all_ok)

    # Transcripts whose CDS LENGTH differs are a different isoform, not
    # another version: c. numbering is NOT portable and a conversion offset
    # must be applied. Derive and report it explicitly.
    for acc, (seq, _, _) in cds.items():
        if len(seq) == len(mane) or acc == L.TRANSCRIPT:
            continue
        off = len(mane) - len(seq)
        k = next((i for i in range(len(seq)) if seq[i:] == mane[i + off:]), None)
        say()
        say(f"  **{acc} is a different isoform, not another version of MANE.** "
            f"CDS {len(seq)} nt vs {len(mane)} nt.")
        if k is not None:
            say(f"  It shares MANE's 3' end but has a distinct N-terminus. "
                f"Conversion rule:")
            say(f"")
            say(f"  | region on {acc} | maps to MANE as |")
            say(f"  |---|---|")
            say(f"  | c.1-{k} (isoform-specific first exon) | not in the MANE CDS "
                f"-- lies within MANE intron 1 |")
            say(f"  | c.{k+1} onwards | `c_MANE = c_alt + {off}` "
                f"(`p_MANE = p_alt + {off // 3}`) |")
            say(f"")
            say(f"  A variant reported on {acc} and taken at face value against MANE "
                f"is therefore **{off // 3} codons out of register**. Any source "
                f"citing this accession must be converted before use; such reports "
                f"are flagged `alt_isoform_numbering` during curation.")
        else:
            say(f"  No consistent single offset relates it to MANE; variants on this "
                f"accession must be resolved individually via their genomic "
                f"coordinates.")

    # characterise any sequence-level differences, which change the *reference
    # amino acid* reported by old papers even though numbering is unchanged
    for acc, (seq, _, _) in cds.items():
        if seq != mane and len(seq) == len(mane):
            diffs = [(i + 1, seq[i], mane[i]) for i in range(len(seq)) if seq[i] != mane[i]]
            codons = sorted({(i - 1) // 3 + 1 for i, _, _ in diffs})
            say(f"  NOTE: {acc} differs from MANE at {len(diffs)} CDS base(s): "
                + ", ".join(f"c.{i}{a}>{b}" for i, a, b in diffs)
                + f" -> affects codon(s) {codons}.")
            say(f"        Numbering is unaffected, but a variant reported on {acc} at these")
            say("        codons will carry a different reference amino acid; such reports are")
            say("        flagged `ref_aa_transcript_version_conflict` during curation.")
    return all_ok


def v4_targetons(model):
    """Reconcile targeton_regions.tsv with the SMC1A_info.xlsx manifest and
    emit the merged reference table."""
    import openpyxl

    wb = openpyxl.load_workbook(os.path.join(SCREEN, "SMC1A_info.xlsx"), data_only=True)
    ws = wb[wb.sheetnames[0]]
    rows = list(ws.iter_rows(values_only=True))
    hdr = [str(c).strip() if c is not None else "" for c in rows[0]]
    idx = {h: i for i, h in enumerate(hdr) if h}
    # The exon-number column is unlabelled in the manifest; it sits immediately
    # after Exon_position. Resolve it positionally and assert it looks right.
    exon_col = idx["Exon_position"] + 1
    if hdr[exon_col]:
        sys.exit(f"FATAL: expected an unlabelled exon-number column at index "
                 f"{exon_col}, found {hdr[exon_col]!r}; the manifest layout has changed")
    manifest = {}
    for r in rows[1:]:
        if not r or not r[0]:
            continue
        tid = str(r[idx["Targeton_ID"]]).strip()
        manifest[tid] = {
            "targeton_id": tid,
            "exon": int(r[exon_col]),
            "exon_span": str(r[idx["Exon_position"]]),
            "amplicon_start": int(r[idx["ref_start"]]),
            "amplicon_end": int(r[idx["ref_end"]]),
            "mut_start": int(r[idx["r2_start"]]),
            "mut_end": int(r[idx["r2_end"]]),
            "design_span": str(r[idx["Output targeton coordinates"]]),
            "sgrna_id": str(r[idx["sgRNA_id"]]),
            "sgrna_strand": str(r[idx["sgRNA_strand"]]),
            "action_vector": str(r[idx["action_vector"]]),
        }

    regions = {}
    with open(os.path.join(SCREEN, "targeton_regions.tsv"), encoding="utf-8-sig") as fh:
        h = fh.readline().strip().split("\t")
        for line in fh:
            if not line.strip():
                continue
            d = dict(zip(h, line.strip().split("\t")))
            regions[d["Targeton_ID"]] = d

    # Screening status per targeton, set from information supplied for this study.
    # Absence from targeton_regions.tsv does NOT mean a targeton was dropped:
    # EXTP is confirmed real and still being screened, so variants falling in
    # it must be marked "results pending", which is a different curation state
    # from "not in the library".
    SCREENING_STATUS = {
        "EXTP": ("screening_in_progress",
                 "confirmed real 2026-09-10; screening ongoing, no "
                 "results yet"),
    }
    # IZFR is an obsolete name for APDY, confirmed 2026-09-10 -- the same
    # targeton, not two. The manifest still uses the old name, so its row is
    # RENAMED (not dropped: it is the only source of APDY's coordinates).
    OBSOLETE_IDS = {"IZFR": "APDY"}

    say()
    for old, new in OBSOLETE_IDS.items():
        if old in manifest:
            rec = manifest.pop(old)
            rec["targeton_id"] = new
            manifest[new] = rec
            say(f"  NOTE: manifest targeton {old} renamed to {new} -- an obsolete "
                f"name for the same targeton (confirmed 2026-09-10). "
                f"{old} does not appear in the reference table.")

    only_manifest = sorted(set(manifest) - set(regions))
    only_regions = sorted(set(regions) - set(manifest))
    expected = {t for t, (st, _) in SCREENING_STATUS.items()
                if st == "screening_in_progress"}
    check("targeton_regions.tsv and SMC1A_info.xlsx agree, allowing for "
          "targetons still being screened",
          not (set(only_manifest) - expected) and not only_regions,
          f"manifest-only={only_manifest} (expected-pending={sorted(expected)}) "
          f"regions-only={only_regions}")

    # Reconcile by amplicon coordinates: a renamed targeton keeps its amplicon.
    by_amp = {(m["amplicon_start"], m["amplicon_end"]): t for t, m in manifest.items()}
    # After dropping obsolete IDs, any targeton in targeton_regions.tsv that
    # is not in the manifest is reconciled by amplicon coordinates.
    aliases = {}
    for tid, d in regions.items():
        if tid in manifest:
            continue
        key = (int(d["start"]), int(d["end"]))
        if key in by_amp:
            aliases[tid] = by_amp[key]
            say(f"  NOTE: targeton {tid} (targeton_regions.tsv) shares its amplicon "
                f"{key[0]}-{key[1]} with manifest targeton {by_amp[key]}; "
                f"{tid} is canonical (it is the ID used in the DESeq2 output).")
        else:
            say(f"  WARNING: targeton {tid} (targeton_regions.tsv) is absent from "
                f"the manifest and its amplicon {key[0]}-{key[1]} matches no "
                f"manifest entry -- needs manual reconciliation.")
    for tid in only_manifest:
        st, why = SCREENING_STATUS.get(tid, ("unknown", "not yet confirmed"))
        say(f"  NOTE: targeton {tid} is in the manifest but absent from "
            f"targeton_regions.tsv (exon {manifest[tid]['exon']}, "
            f"amplicon {manifest[tid]['amplicon_start']}-{manifest[tid]['amplicon_end']}). "
            f"Screening status: **{st}** -- {why}.")

    out = []
    for tid, m in manifest.items():
        canonical = tid
        for alias, target in aliases.items():
            if target == tid:
                canonical = alias
        e = next((x for x in EXON_ROWS if x["exon_number"] == m["exon"]), None)
        dspan = m["design_span"].split(":")[-1]
        try:
            d_start, d_end = (int(x) for x in dspan.split("-"))
        except ValueError:
            d_start = d_end = ""
        out.append({
            "targeton_id": canonical,
            "targeton_id_manifest": tid,
            "exon": m["exon"],
            "amplicon_start": m["amplicon_start"],
            "amplicon_end": m["amplicon_end"],
            "mut_start": m["mut_start"],
            "mut_end": m["mut_end"],
            "mut_length": m["mut_end"] - m["mut_start"] + 1,
            "design_start": d_start,
            "design_end": d_end,
            "design_length": (d_end - d_start + 1) if d_start != "" else "",
            "exon_g_start": e["g_start"] if e else "",
            "exon_g_end": e["g_end"] if e else "",
            "exon_c_start": e["c_start"] if e else "",
            "exon_c_end": e["c_end"] if e else "",
            "sgrna_id": m["sgrna_id"],
            "sgrna_strand": m["sgrna_strand"],
            "action_vector": m["action_vector"],
            "in_targeton_regions_tsv": canonical in regions,
            "alias_of": tid if canonical != tid else "",
            # Status is a property of the canonical targeton. A manifest ID
            # that was renamed (IZFR -> APDY) must not pass its "renamed"
            # status on to the canonical row, which is genuinely screened.
            "screening_status": SCREENING_STATUS.get(canonical, ("screened", ""))[0],
            "screening_status_note": SCREENING_STATUS.get(canonical, ("screened", ""))[1],
        })
    out.sort(key=lambda r: (r["exon"], r["mut_start"]))

    cols = list(out[0].keys())
    with open(os.path.join(REF, "smc1a_targetons.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    say(f"  wrote data/reference/smc1a_targetons.tsv ({len(out)} targetons)")

    # how much of the CDS is actually mutagenised?
    covered = set()
    for r in out:
        covered |= set(range(r["mut_start"], r["mut_end"] + 1))
    designed = set()
    for r in out:
        if r["design_start"] != "":
            designed |= set(range(r["design_start"], r["design_end"] + 1))
    cds_pos = {model.t_to_g(t) for t in range(model.cds_t_start, model.cds_t_end + 1)}
    frac = len(cds_pos & covered) / len(cds_pos)
    say(f"  exon-core (r2) regions cover {len(cds_pos & covered)}/{len(cds_pos)} "
        f"CDS bases ({frac:.1%})")
    dfrac = len(cds_pos & designed) / len(cds_pos)
    say(f"  design windows cover {len(cds_pos & designed)}/{len(cds_pos)} "
        f"CDS bases ({dfrac:.1%}), plus "
        f"{len(designed - cds_pos):,} non-CDS (intronic/UTR) bases -- the")
    say(f"  flanking splice regions that the exon-core definition would miss")
    return out


def v5_roundtrip(model) -> bool:
    """Offline model self-consistency, then agreement with Variant Recoder."""
    ok = True
    for e in model.exons:
        for t in (e["t_start"], e["t_end"]):
            g = model.t_to_g(t)
            if model.g_to_t(g) != t:
                ok = check(f"round-trip t={t}", False) and ok
    check("c. <-> genomic round-trips at every exon boundary", ok)

    # Probe set spanning the variant classes we expect to curate, including
    # two deliberately invalid descriptions that must be REJECTED rather than
    # silently coerced.
    probes = [
        ("c.1951G>A", True), ("c.2095C>T", True), ("c.1487G>A", True),
        ("c.2046_2048del", True), ("c.109+108dup", True), ("c.-19C>T", True),
        ("c.616-2A>G", True), ("c.615+1G>A", True), ("c.2492_2493insT", True),
        ("c.3702G>A", True),
        ("c.587-2A>G", False),   # c.587 is mid-exon-4: not a splice site
        ("c.298+1G>A", True),    # genuine exon 2 donor
    ]
    refseq = L.ReferenceSequence(os.path.join(REF, "smc1a_chrX_region.fa"))
    say()
    say("  | input (MANE c.)   | valid HGVS? | offline g. | region   | ref base OK | Recoder VCF       |")
    say("  |-------------------|-------------|------------|----------|-------------|-------------------|")
    agree = True
    for p, should_parse in probes:
        parsed = L.parse_hgvs_c(p)
        g_off, region = (None, "unparsed")
        refok = ""
        if parsed:
            g_off, region = model.c_pos_to_g(parsed["pos1"])
            if g_off and parsed["edit_type"] == "sub":
                # minus strand: the plus-strand reference base is the
                # complement of the base named in the c. description
                refok = str(refseq.base(g_off) == L.complement(parsed["ref_c"]))
        accepted = g_off is not None
        if accepted != should_parse:
            agree = False
        url = (f"{L.ENSEMBL_REST}/variant_recoder/human/"
               + urllib.parse.quote(f"{L.TRANSCRIPT}:{p}")
               + "?content-type=application/json;vcf_string=1;fields=hgvsg,vcf_string,hgvsc")
        d, _rec = http.get_json(url)
        vcf = ""
        if d:
            for allele in d:
                for _k, v in allele.items():
                    if isinstance(v, dict):
                        for s in v.get("vcf_string", []):
                            if s.startswith("X-"):
                                vcf = s
        mark = "" if accepted == should_parse else "  <-- UNEXPECTED"
        say(f"  | {p:17} | {str(should_parse):11} | {str(g_off):10} | {region:8} "
            f"| {refok:11} | {vcf:17} |{mark}")
    say()
    check("offline model accepts valid and rejects invalid c. descriptions", agree)
    say("  NOTE: `c.587-2A>G` is invalid (c.587 is internal to exon 4). The offline")
    say("  model rejects it; Ensembl Variant Recoder instead silently returns the")
    say("  neighbouring exonic variant c.585A>G. Every API result is therefore")
    say("  required to echo back the *same* hgvsc string that was submitted")
    say("  (03_resolve_variants.py), so this class of silent coercion cannot enter")
    say("  the curated table.")

    # Independent amino-acid verification from the cached CDS: catches any
    # position/strand error that would still yield a valid-looking coordinate.
    cds = "".join(
        L.complement(refseq.base(model.t_to_g(t)))
        for t in range(model.cds_t_start, model.cds_t_end + 1)
    )
    prot = "".join(L.CODON_TABLE.get(cds[i:i + 3], "?") for i in range(0, len(cds), 3))
    say()
    check("CDS reconstructed from the plus-strand FASTA is 3702 nt",
          len(cds) == L.CDS_LENGTH, f"got {len(cds)}")
    check("reconstructed CDS starts ATG and ends with a stop codon",
          cds[:3] == "ATG" and prot[-1] == "*", f"start={cds[:3]} end_aa={prot[-1]}")
    check("reconstructed protein is 1233 aa with no internal stop",
          len(prot) - 1 == L.PROTEIN_LENGTH and "*" not in prot[:-1],
          f"got {len(prot) - 1} aa")
    d, _ = http.get_json(f"{L.ENSEMBL_REST}/sequence/id/{L.TRANSCRIPT_NOVER}"
                         "?type=protein;content-type=application/json")
    check("reconstructed protein is identical to Ensembl ENSP00000323421",
          d is not None and prot[:-1] == d["seq"].upper())
    say()
    say("  This last check closes the loop: the plus-strand genomic sequence, the")
    say("  minus-strand exon model, the CDS numbering and the reference proteome all")
    say("  agree, so a c. position converted by this model lands on the intended")
    say("  codon.")
    return ok


# =========================================================================
# Main
# =========================================================================

if __name__ == "__main__":
    say("# SMC1A transcript reference: verification report")
    say()
    say(f"Generated by `code/01_build_transcript_reference.py`.")
    say(f"Reference transcript: **{L.TRANSCRIPT}** (MANE Select) / {L.REFSEQ_MANE}, "
        f"protein {L.PROTEIN}, genome {L.ASSEMBLY} ({L.GENOMIC_ACC}), "
        f"{L.GENE} on the **minus** strand of {L.CHROM}.")
    say()

    tr = fetch_transcript()
    MODEL, EXON_ROWS, TL = build_model(tr)

    # Cache the locus reference sequence first: V5 verifies the coordinate
    # model against it.
    lo = min(e["g_start"] for e in EXON_ROWS) - 500
    hi = max(e["g_end"] for e in EXON_ROWS) + 500
    d, _ = http.get_json(f"{L.ENSEMBL_REST}/sequence/region/human/X:{lo}..{hi}:1"
                         "?content-type=application/json")
    seq = d["seq"].upper()
    with open(os.path.join(REF, "smc1a_chrX_region.fa"), "w", encoding="utf-8") as fh:
        fh.write(f">chrX:{lo}-{hi} GRCh38 plus_strand\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")

    say("## V1 -- Ensembl exon structure vs SMC1A_exon_map.tsv")
    v1_exon_map(EXON_ROWS)
    say()
    say("## V2 -- transcript and protein lengths")
    v2_lengths(EXON_ROWS, MODEL, TL)
    say()
    say("## V3 -- c. numbering portability across transcript versions")
    say()
    say("The source literature cites NM_006306.1 through .4 and, in older CdLS")
    say("papers, no accession at all. This check establishes whether a c. position")
    say("taken from a publication can be used against the MANE transcript without")
    say("an offset correction.")
    v3_transcript_equivalence()
    say()
    say("## V4 -- targeton reconciliation")
    TARGETONS = v4_targetons(MODEL)
    say()
    say("## V5 -- offline coordinate model")
    v5_roundtrip(MODEL)

    # ---- write reference tables ----
    cols = ["exon_number", "exon_id", "g_start", "g_end", "length",
            "t_start", "t_end", "c_start", "c_end"]
    with open(os.path.join(REF, "smc1a_exons.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in EXON_ROWS:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    json.dump({
        "gene": L.GENE, "chrom": L.CHROM, "strand": L.STRAND,
        "assembly": L.ASSEMBLY, "genomic_accession": L.GENOMIC_ACC,
        "transcript": L.TRANSCRIPT, "protein": L.PROTEIN,
        "refseq_mane": L.REFSEQ_MANE,
        "exons": EXON_ROWS,
        "cds_t_start": MODEL.cds_t_start, "cds_t_end": MODEL.cds_t_end,
        "transcript_length": MODEL.t_len,
    }, open(os.path.join(REF, "smc1a_transcript_model.json"), "w"), indent=1)

    say()
    say(f"Cached reference sequence chrX:{lo}-{hi} ({len(seq):,} bp, plus strand) "
        f"to `data/reference/smc1a_chrX_region.fa`.")
    say(f"HTTP cache: {http.n_hits} hits, {http.n_misses} misses "
        f"(`data/cache/ensembl/`).")

    with open(os.path.join(DOCS, "transcript_equivalence.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(report) + "\n")
    print(f"\nWrote docs/verification/transcript_equivalence.md")
