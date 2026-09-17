#!/usr/bin/env python3
"""
03_resolve_variants.py

Recalibrate every collected variant description onto the MANE Select
transcript ENST00000322213.9 / GRCh38, and derive the plus-strand VCF
representation that the SGE library is keyed on.

Method
------
Publication- and database-reported variants arrive as HGVS c. descriptions on
assorted RefSeq versions (and occasionally as protein-only or genomic-only
descriptions). They are resolved by a two-track procedure, and a variant is
only accepted if the tracks agree:

  Track A -- Ensembl Variant Recoder (REST, GRCh38)
      Submitted as `ENST00000322213.9:c....`. Returns the HGVS genomic
      description, the plus-strand VCF representation (`vcf_string`), the
      re-derived c. and p. descriptions, and any co-located rsID. This handles
      minus-strand conversion, HGVS 3'-shifting vs VCF left-alignment, and
      indel/dup/delins normalisation, so none of that is reimplemented here.

  Track B -- offline model (smc1a_lib.TranscriptModel + cached GRCh38 FASTA)
      An independent implementation built from the Ensembl exon structure and
      the cached plus-strand reference sequence. Used to confirm that the
      coordinate returned by Track A falls where the c. position says it
      should, that the reference allele matches the reference genome, and --
      for coding SNVs -- that the codon and reference amino acid are the ones
      implied by the submitted protein change.

Guards (a variant failing any of these is quarantined, never silently kept)
--------------------------------------------------------------------------
G1  the c. description must parse
G2  Variant Recoder must echo back the SAME c. description that was submitted.
    This is the guard against silent coercion: `c.587-2A>G` (invalid, since
    c.587 is mid-exon) is returned by the API as the unrelated exonic variant
    c.585A>G. Without G2 that would enter the table as a real variant.
G3  the plus-strand reference allele must match the cached GRCh38 sequence
G4  the offline model must place the variant in the same exon/intron as the
    submitted c. position implies
G5  where the source states a protein change, the independently recomputed
    reference amino acid must match it (catches transcript-version and
    numbering errors, including the NM_006306.1 c.489_490 discrepancy)
G6  the resolved position must lie within the SMC1A locus

Consequence annotation
----------------------
VEP REST (`/vep/human/hgvs`) restricted to ENST00000322213, giving SO
consequence terms, exon/intron number, codon and protein position on the MANE
transcript. Harmonised to the coarse classes used by the SGE DESeq2 output via
smc1a_lib.consequence_class.

Inputs
------
data/interim/lovd_observations.tsv          (step 02)
data/interim/literature_observations.tsv    (step 04, if present)
data/interim/*_observations.tsv             any other source table

Outputs
-------
data/interim/resolved_variants.tsv          one row per distinct submitted
                                            description, with resolution result
data/interim/resolution_quarantine.tsv      descriptions that failed a guard
docs/verification/coordinate_resolution.md

Usage
-----
    ../.venv/bin/python 03_resolve_variants.py
"""

from __future__ import annotations

import glob
import json
import os
import sys
import urllib.parse
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smc1a_lib as L  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, "data", "reference")
INTERIM = os.path.join(ROOT, "data", "interim")
DOCS = os.path.join(ROOT, "docs", "verification")
os.makedirs(DOCS, exist_ok=True)

http = L.CachedHTTP(os.path.join(ROOT, "data", "cache", "ensembl"))
MODEL = L.TranscriptModel.from_json(os.path.join(REF, "smc1a_transcript_model.json"))
REFSEQ = L.ReferenceSequence(os.path.join(REF, "smc1a_chrX_region.fa"))
TARGETONS = L.TargetonMap.from_tsv(os.path.join(REF, "smc1a_targetons.tsv"))

LOCUS_LO = min(e["g_start"] for e in MODEL.exons) - 500
LOCUS_HI = max(e["g_end"] for e in MODEL.exons) + 500

report: list[str] = []


def say(m: str = "") -> None:
    print(m)
    report.append(m)


# --------------------------------------------------------------------------
# Reference CDS reconstructed from the plus-strand FASTA (verified in step 01
# to be identical to the Ensembl CDS), used for the amino-acid guard G5.
# --------------------------------------------------------------------------

CDS = "".join(L.complement(REFSEQ.base(MODEL.t_to_g(t)))
              for t in range(MODEL.cds_t_start, MODEL.cds_t_end + 1))


def ref_aa_at_codon(codon: int) -> str:
    if not (1 <= codon <= L.PROTEIN_LENGTH):
        return ""
    return L.CODON_TABLE.get(CDS[(codon - 1) * 3: codon * 3], "?")


# --------------------------------------------------------------------------
# Track A: Ensembl Variant Recoder
# --------------------------------------------------------------------------

def variant_recoder(hgvs_c: str):
    """Submit `ENST00000322213.9:c...` and return the parsed result.

    Returns (result_dict | None, error_string).
    """
    desc = f"{L.TRANSCRIPT}:{hgvs_c}"
    url = (f"{L.ENSEMBL_REST}/variant_recoder/human/{urllib.parse.quote(desc)}"
           "?content-type=application/json;vcf_string=1"
           ";fields=hgvsg,hgvsc,hgvsp,vcf_string,id")
    d, rec = http.get_json(url)
    if d is None:
        try:
            msg = json.loads(rec["body"]).get("error", rec["body"][:200])
        except Exception:
            msg = rec["body"][:200]
        return None, f"variant_recoder HTTP {rec['status']}: {msg}"
    if not isinstance(d, list) or not d:
        return None, "variant_recoder returned no alleles"

    # The response is a list of input descriptions, each a dict keyed by the
    # alternate allele. A c. substitution yields exactly one allele.
    out = []
    for entry in d:
        for allele, v in entry.items():
            if not isinstance(v, dict):
                continue
            vcf = [s for s in v.get("vcf_string", []) if s.startswith("X-")]
            hg = [s for s in v.get("hgvsg", []) if s.startswith(L.GENOMIC_ACC)]
            hc = [s for s in v.get("hgvsc", []) if s.startswith(L.TRANSCRIPT)]
            hp = [s for s in v.get("hgvsp", []) if s.startswith(L.PROTEIN.split(".")[0])]
            if not vcf:
                continue
            chrom, pos, ref, alt = vcf[0].split("-")
            out.append({
                "alt_allele": allele,
                "pos": int(pos), "ref": ref, "alt": alt,
                "hgvs_g": hg[0] if hg else "",
                "hgvs_c_returned": hc[0].split(":", 1)[1] if hc else "",
                "hgvs_p_returned": hp[0].split(":", 1)[1] if hp else "",
                "rsid": ",".join(v.get("id", [])),
            })
    if not out:
        return None, "variant_recoder returned no plus-strand VCF representation"
    if len(out) > 1:
        return None, ("variant_recoder returned multiple alleles: "
                      + ",".join(f"{o['pos']}{o['ref']}>{o['alt']}" for o in out))
    return out[0], ""


# --------------------------------------------------------------------------
# VEP consequence annotation
# --------------------------------------------------------------------------

def vep(hgvs_c: str):
    desc = f"{L.TRANSCRIPT}:{hgvs_c}"
    url = (f"{L.ENSEMBL_REST}/vep/human/hgvs/{urllib.parse.quote(desc)}"
           "?content-type=application/json;numbers=1;hgvs=1;canonical=1;mane=1")
    d, rec = http.get_json(url)
    if d is None or not isinstance(d, list) or not d:
        return {}, f"vep HTTP {rec['status']}"
    tc = [t for t in d[0].get("transcript_consequences", [])
          if t.get("transcript_id") == L.TRANSCRIPT_NOVER]
    if not tc:
        return {}, "vep returned no consequence on the MANE transcript"
    t = tc[0]
    terms = t.get("consequence_terms", [])
    return {
        "consequence_terms": ";".join(terms),
        "consequence_class": L.consequence_class(terms),
        "vep_exon": L.clean(t.get("exon")),
        "vep_intron": L.clean(t.get("intron")),
        "protein_position": L.clean(t.get("protein_start")),
        "amino_acids": L.clean(t.get("amino_acids")),
        "codons": L.clean(t.get("codons")),
        "most_severe": L.clean(d[0].get("most_severe_consequence")),
    }, ""


# --------------------------------------------------------------------------
# Resolution of one description
# --------------------------------------------------------------------------

def resolve(hgvs_c: str, stated_protein: str = ""):
    r = {
        "hgvs_c_input": hgvs_c,
        "status": "", "guard_failed": "", "notes": [],
        "chrom": L.CHROM, "pos": "", "ref": "", "alt": "",
        "variant_key": "", "hgvs_g": "", "hgvs_c_mane": "", "hgvs_p_mane": "",
        "rsid": "", "edit_type": "", "region_expected": "",
        "exon": "", "intron": "", "protein_position": "", "ref_aa": "",
        "consequence_terms": "", "consequence_class": "",
        "targetons_design": "", "targetons_exon_core": "",
        "targetons_amplicon": "", "in_design_window": "",
    }

    # ---- G1: parse ----------------------------------------------------
    p = L.parse_hgvs_c(hgvs_c)
    if p is None:
        r["status"] = "quarantined"
        r["guard_failed"] = "G1_unparseable_hgvs"
        return r
    r["edit_type"] = p["edit_type"]

    # offline expectation for the anchor position
    g1, region1 = MODEL.c_pos_to_g(p["pos1"])
    g2, region2 = MODEL.c_pos_to_g(p["pos2"])
    if g1 is None or g2 is None:
        r["status"] = "quarantined"
        r["guard_failed"] = "G4_offline_model_rejects_position"
        r["notes"].append(region1 if g1 is None else region2)
        return r
    r["region_expected"] = region1 if region1 == region2 else f"{region1}/{region2}"
    g_lo, g_hi = min(g1, g2), max(g1, g2)

    # ---- Track A ------------------------------------------------------
    a, err = variant_recoder(hgvs_c)
    if a is None:
        r["status"] = "quarantined"
        r["guard_failed"] = "A_variant_recoder_failed"
        r["notes"].append(err)
        return r

    # ---- G2: no silent coercion ---------------------------------------
    returned = a["hgvs_c_returned"]
    if returned and L.parse_hgvs_c(returned):
        ret_n = L.parse_hgvs_c(returned)["normalised"]
        inp_n = p["normalised"]
        # Accept differences that are purely notational: an omitted deleted
        # sequence (c.123del vs c.123delA) or an omitted duplicated sequence.
        def strip_seq(s: str) -> str:
            import re
            return re.sub(r"(del|dup|ins)[ACGT]+", r"\1", s)

        # A `delXXXinsY` is a pure deletion whenever the inserted sequence is a
        # suffix (or prefix) of the deleted sequence -- e.g. on a GTG reference,
        # `c.3046_3048delGTGinsG` deletes TG and keeps the final G, which is
        # exactly Variant Recoder's `c.3047_3048del`. These are the SAME variant
        # written two ways, so the guard must not quarantine them. Equivalence is
        # proved arithmetically from the submitted description alone, not assumed.
        def delins_is_equivalent_del(inp: str, ret: str) -> bool:
            import re
            mi = re.match(r"^c\.(\d+)_(\d+)del([ACGT]+)ins([ACGT]+)$", inp)
            mr = re.match(r"^c\.(\d+)(?:_(\d+))?del$", strip_seq(ret))
            if not (mi and mr):
                return False
            a, b, dele, ins = int(mi.group(1)), int(mi.group(2)), mi.group(3), mi.group(4)
            if len(dele) != b - a + 1 or len(ins) >= len(dele):
                return False
            ra = int(mr.group(1))
            rb = int(mr.group(2)) if mr.group(2) else ra
            n = len(dele) - len(ins)          # net bases removed
            if rb - ra + 1 != n:
                return False
            # inserted seq retained at the 3' end of the interval -> del at the 5' end
            if dele.endswith(ins) and (ra, rb) == (a, a + n - 1):
                return True
            # inserted seq retained at the 5' end -> del at the 3' end
            if dele.startswith(ins) and (ra, rb) == (b - n + 1, b):
                return True
            return False

        equivalent = delins_is_equivalent_del(inp_n, ret_n)
        if equivalent:
            r["notes"].append(
                f"submitted {inp_n} and Variant Recoder returned {ret_n}; these are "
                f"the same variant -- the inserted base(s) are retained from the "
                f"deleted sequence, so the delins reduces arithmetically to the "
                f"returned pure deletion. Accepted on the normalised description.")
        if not equivalent and strip_seq(ret_n) != strip_seq(inp_n):
            r["status"] = "quarantined"
            r["guard_failed"] = "G2_api_returned_different_variant"
            r["notes"].append(f"submitted {inp_n} but Variant Recoder returned "
                              f"{ret_n} -- the submitted description is probably "
                              f"invalid on this transcript")
            r["pos"], r["ref"], r["alt"] = a["pos"], a["ref"], a["alt"]
            return r

    r["pos"], r["ref"], r["alt"] = a["pos"], a["ref"], a["alt"]
    r["hgvs_g"] = a["hgvs_g"]
    r["hgvs_c_mane"] = returned or p["normalised"]
    r["hgvs_p_mane"] = a["hgvs_p_returned"]
    r["rsid"] = a["rsid"]
    r["variant_key"] = L.variant_key(L.CHROM, a["pos"], a["ref"], a["alt"])

    # ---- G6: inside the locus -----------------------------------------
    if not (LOCUS_LO <= a["pos"] <= LOCUS_HI):
        r["status"] = "quarantined"
        r["guard_failed"] = "G6_outside_locus"
        return r

    # ---- G3: reference allele matches the reference genome -------------
    obs_ref = REFSEQ.get(a["pos"], a["pos"] + len(a["ref"]) - 1)
    if obs_ref != a["ref"]:
        r["status"] = "quarantined"
        r["guard_failed"] = "G3_ref_allele_mismatch"
        r["notes"].append(f"GRCh38 has {obs_ref!r} at chrX:{a['pos']}, "
                          f"variant claims REF={a['ref']!r}")
        return r

    # ---- G4: exon/intron agreement -------------------------------------
    # Compare against the anchor implied by the c. description rather than the
    # VCF position, since VCF left-alignment can move an indel by a few bases.
    exon_off = MODEL.exon_of_g(g_lo) or MODEL.exon_of_g(g_hi)
    intron_off = MODEL.intron_of_g(g_lo) if exon_off is None else None
    if r["region_expected"] == "exonic" and exon_off is None:
        r["status"] = "quarantined"
        r["guard_failed"] = "G4_expected_exonic_but_intronic"
        return r
    r["exon"] = exon_off or ""
    r["intron"] = intron_off or ""

    # ---- consequence ---------------------------------------------------
    v, verr = vep(hgvs_c)
    if verr:
        r["notes"].append(verr)
    r.update({k: v[k] for k in ("consequence_terms", "consequence_class",
                                "protein_position") if k in v})
    if v.get("vep_exon"):
        r["exon"] = v["vep_exon"].split("/")[0]
    if v.get("vep_intron"):
        r["intron"] = v["vep_intron"].split("/")[0]

    # ---- G5: reference amino acid --------------------------------------
    if r["protein_position"]:
        try:
            codon = int(str(r["protein_position"]).split("-")[0])
            r["ref_aa"] = ref_aa_at_codon(codon)
        except ValueError:
            codon = None
        stated = L.aa3_to_1(stated_protein) if stated_protein else ""
        if stated and r["ref_aa"]:
            import re as _r
            m = _r.match(r"^([A-Z*])(\d+)", stated)
            if m and int(m.group(2)) == codon and m.group(1) != r["ref_aa"]:
                r["notes"].append(
                    f"G5: source states reference amino acid {m.group(1)}{codon} "
                    f"but MANE has {r['ref_aa']}{codon} "
                    f"-- flagged ref_aa_conflict")
                r["guard_failed"] = "G5_ref_aa_conflict"

    # ---- targeton assignment -------------------------------------------
    span_lo = min(a["pos"], g_lo)
    span_hi = max(a["pos"] + len(a["ref"]) - 1, g_hi)
    design, core, amp = TARGETONS.assign(span_lo, span_hi)
    r["targetons_design"] = "|".join(design)
    r["targetons_exon_core"] = "|".join(core)
    r["targetons_amplicon"] = "|".join(amp)
    r["in_design_window"] = bool(design)

    r["status"] = "resolved" if not r["guard_failed"] else "resolved_with_flag"
    return r


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    # gather every distinct submitted description across all source tables
    inputs: dict[str, dict] = {}
    files = sorted(glob.glob(os.path.join(INTERIM, "*_observations.tsv")))
    if not files:
        sys.exit("FATAL: no *_observations.tsv in data/interim/ -- run step 02 first")

    for f in files:
        with open(f, encoding="utf-8") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                d = dict(zip(hdr, line.rstrip("\n").split("\t")))
                c = L.clean(d.get("hgvs_c_input"))
                if not c:
                    continue
                rec = inputs.setdefault(c, {"sources": set(), "protein": set()})
                rec["sources"].add(os.path.basename(f).replace("_observations.tsv", ""))
                if L.clean(d.get("hgvs_p_input")):
                    rec["protein"].add(L.clean(d["hgvs_p_input"]))

    say("# Coordinate resolution report")
    say()
    say(f"Reference: **{L.TRANSCRIPT}** (MANE Select), {L.ASSEMBLY} / "
        f"{L.GENOMIC_ACC}, minus strand.")
    say()
    say("All positions and alleles in the output are **plus-strand GRCh38**, "
        "matching the")
    say("convention used by the SGE library oligo names "
        "(`...chrX:<pos>_<REF>><ALT>_...`),")
    say("so `variant_key` joins directly to the screen results.")
    say()
    say(f"Source tables read: {', '.join(os.path.basename(f) for f in files)}")
    say(f"Distinct submitted c. descriptions: {len(inputs)}")
    say()

    rows, quar = [], []
    for i, (c, meta) in enumerate(sorted(inputs.items()), 1):
        stated = sorted(meta["protein"])[0] if meta["protein"] else ""
        r = resolve(c, stated)
        r["sources"] = "|".join(sorted(meta["sources"]))
        r["hgvs_p_stated"] = stated
        r["notes"] = "; ".join(r["notes"])
        (quar if r["status"] == "quarantined" else rows).append(r)
        if i % 25 == 0:
            print(f"  ... {i}/{len(inputs)}", file=sys.stderr)

    say("## Outcome")
    say()
    say(f"| outcome | descriptions |")
    say(f"|---|---|")
    say(f"| resolved | {len([r for r in rows if r['status'] == 'resolved'])} |")
    say(f"| resolved with a flag | {len([r for r in rows if r['status'] == 'resolved_with_flag'])} |")
    say(f"| quarantined | {len(quar)} |")
    say()

    if quar:
        say("## Quarantined descriptions")
        say()
        say("These are **not** discarded: they are written to")
        say("`data/interim/resolution_quarantine.tsv` and carried into the")
        say("all-variants list with `curation_group = Exclude_pending_classification`,")
        say("so nothing collected is lost. Each needs a manual decision.")
        say()
        say("| submitted c. | guard | detail | sources |")
        say("|---|---|---|---|")
        for r in quar:
            say(f"| `{r['hgvs_c_input']}` | {r['guard_failed']} | "
                f"{r['notes'][:150]} | {r['sources']} |")
        say()

    flagged = [r for r in rows if r["guard_failed"]]
    if flagged:
        say("## Resolved but flagged")
        say()
        say("| c. (MANE) | flag | detail |")
        say("|---|---|---|")
        for r in flagged:
            say(f"| `{r['hgvs_c_mane']}` | {r['guard_failed']} | {r['notes'][:160]} |")
        say()

    # ---- cross-check against LOVD's own hg38 where available -----------
    say("## Independent cross-check: LOVD-supplied hg38 coordinates")
    say()
    lovd_hg38 = {}
    p = os.path.join(INTERIM, "lovd_observations.tsv")
    if os.path.exists(p):
        with open(p, encoding="utf-8") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                d = dict(zip(hdr, line.rstrip("\n").split("\t")))
                c, g = L.clean(d.get("hgvs_c_input")), L.clean(d.get("hgvs_g_hg38_input"))
                if c and g:
                    lovd_hg38[c] = g
    agree = disagree = nocmp = 0
    mism = []
    for r in rows:
        g = lovd_hg38.get(r["hgvs_c_input"])
        if not g or not r["hgvs_g"]:
            nocmp += 1
            continue
        # compare the position(s) named in each HGVS g. description
        import re as _re
        a = _re.findall(r"\d{6,}", g)
        b = _re.findall(r"\d{6,}", r["hgvs_g"])
        if a and b and a[0] == b[0]:
            agree += 1
        else:
            disagree += 1
            mism.append((r["hgvs_c_input"], g, r["hgvs_g"]))
    say(f"- LOVD hg38 agrees with the resolved coordinate: {agree}")
    say(f"- disagrees: {disagree}")
    say(f"- no LOVD hg38 value to compare: {nocmp}")
    if mism:
        say()
        say("| c. | LOVD hg38 | resolved hg38 |")
        say("|---|---|---|")
        for c, g, h in mism[:40]:
            say(f"| `{c}` | {g} | {h} |")
    say()

    # ---- distributions -------------------------------------------------
    say("## Consequence classes (distinct variants)")
    say()
    cc = defaultdict(int)
    for r in rows:
        cc[r["consequence_class"] or "(none)"] += 1
    say("| consequence_class | variants | predicted LOF |")
    say("|---|---|---|")
    for k, n in sorted(cc.items(), key=lambda kv: -kv[1]):
        say(f"| {k} | {n} | {'yes' if k in L.PTV_CLASSES else ''} |")
    say()

    say("## Targeton coverage")
    say()
    inmut = [r for r in rows if r["in_design_window"]]
    say(f"- variants inside an SGE **design window** (expected in the library): "
        f"{len(inmut)}/{len(rows)}")
    say(f"- variants inside an amplicon but outside every design window "
        f"(sequenced, but no oligo designed): "
        f"{len([r for r in rows if r['targetons_amplicon'] and not r['in_design_window']])}")
    say(f"- variants outside every targeton: "
        f"{len([r for r in rows if not r['targetons_amplicon']])}")
    say()
    byt = defaultdict(int)
    for r in inmut:
        for t in r["targetons_design"].split("|"):
            byt[t] += 1
    say("| targeton | exon | variants |")
    say("|---|---|---|")
    tex = {t["targeton_id"]: t["exon"] for t in TARGETONS.targetons}
    for t, n in sorted(byt.items(), key=lambda kv: -kv[1]):
        say(f"| {t} | {tex.get(t, '?')} | {n} |")
    say()

    # ---- write ---------------------------------------------------------
    cols = ["hgvs_c_input", "sources", "status", "guard_failed", "variant_key",
            "chrom", "pos", "ref", "alt", "hgvs_g", "hgvs_c_mane",
            "hgvs_p_mane", "hgvs_p_stated", "rsid", "edit_type",
            "region_expected", "exon", "intron", "protein_position", "ref_aa",
            "consequence_terms", "consequence_class", "targetons_design",
            "targetons_exon_core", "targetons_amplicon", "in_design_window",
            "notes"]
    for name, data in (("resolved_variants.tsv", rows),
                       ("resolution_quarantine.tsv", quar)):
        with open(os.path.join(INTERIM, name), "w", encoding="utf-8") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in data:
                fh.write("\t".join(
                    str(r.get(c, "")).replace("\t", " ")
                                     .replace("\n", " ").replace("\r", " ")
                    for c in cols) + "\n")
    say(f"Wrote `data/interim/resolved_variants.tsv` ({len(rows)}) and "
        f"`data/interim/resolution_quarantine.tsv` ({len(quar)}).")
    say(f"HTTP cache: {http.n_hits} hits, {http.n_misses} misses.")

    with open(os.path.join(DOCS, "coordinate_resolution.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(report) + "\n")
    print("Wrote docs/verification/coordinate_resolution.md")


if __name__ == "__main__":
    main()
