#!/usr/bin/env python3
"""Verify every literature extraction behind the curated set against its PDF.

For each source, checks that the cDNA and protein descriptions the extraction
records actually appear in the paper. Comparison is notation-insensitive
(3-letter vs 1-letter residues, whitespace, HGVS punctuation) because sources
write the same change many ways; what is being tested is whether the VARIANT
is in the paper, not whether the string is.

A miss is not automatically an error -- a paper may state protein only and
leave the cDNA to be derived, or describe the change in prose (Musio 2006).
Misses are listed for inspection, not failed.
"""
import csv, glob, logging, os, re, sys
logging.disable(logging.CRITICAL)
import pypdf

R = "/Users/vb9/Documents/Claude/SMC1A/thesis_workflow/curation_attempt2"
SRC = f"{R}/variant_curation/sources"

# Sources whose PDF filename does not follow <Author>_<year>
# Supplementary material held alongside the main paper. Several variants are
# reported only there, so a check against the main PDF alone under-reports.
SUPPLEMENTARY = {
    "DiNardo_2026": ["DiNardo_supp.pdf"],
    "DiNardo_2026_CdLS": ["DiNardo_supp.pdf"],
    "Kruszka_2019": ["Kruszka_2019_supplementary/Supplementary_Data1.pdf",
                     "Kruszka_2019_supplementary/Supplementary_Data2.pdf",
                     "Kruszka_2019_supplementary/Supplementary_Data3.pdf"],
    "Chuan_2022": ["Chuan_supplementary/Table 1.XLSX",
                   "Chuan_supplementary/Table 2.XLSX",
                   "Chuan_supplementary/Table 3.XLSX",
                   "Chuan_supplementary/Table 4.XLSX",
                   "Chuan_supplementary/Table 5-2.XLSX"],
    "Yuan_2019": ["Yuan_supplementary/NIHMS973104-supplement-Table_S1.docx",
                  "Yuan_supplementary/NIHMS973104-supplement-Table_S2.xlsx"],
    "Huisman_2017": ["Huisman_supplementary.doc"],
    "Chinen_2019": ["Chinen_supplementary.xlsx"],
}

PDF_OVERRIDE = {
    "Kruszka_2019": "Kruszka_2018_DEE",
    "Jansen_2016": "Jansen_2015_DEE",
    "Gorman_2017": "Gorman",
    "Amllal_2025": "Amllal_2024",
    "Yang_2025_table1": "Yang_2025",
    "DiNardo_2026_CdLS": "DiNardo_2026_DEE_therapeutic",
    "ArefEshghi_2020": "Aref-Eshghi_2020",
    "Tzschach_2015": "Tszchach_2015",   # filename misspells the author
}
# Not papers: databases and unpublished cohorts, verified by other means
NO_PAPER = {"LOVD", "GeneDx_180k_DNM", "DDD_DECIPHER"}

A3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
      "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
      "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
      "Tyr": "Y", "Val": "V", "Ter": "X"}


# Artefacts of extracting text out of typeset PDFs: papers print an arrow for
# the substitution operator, and the plus of an intronic offset and the star of
# a stop routinely come back as mojibake.
PDF_ARTEFACT = {"\u2192": ">", "\u2190": ">", "\u00fe": "+", "\u00de": "+",
                "/C3": "*", "\u2217": "*", "\u00b7": "", "\ufb01": "fi"}


def norm(s: str) -> str:
    s = str(s)
    for k, v in PDF_ARTEFACT.items():
        s = s.replace(k, v)
    for k, v in A3.items():
        s = s.replace(k, v)
    s = s.replace("*", "X")
    # Papers routinely omit the `p.` prefix in tables. Stripping it matters:
    # without this, `p.Gly32Glu` normalises to `pg32e` and does not match the
    # paper's `Gly32Glu` -> `g32e`, reporting a miss that is purely notation.
    s = re.sub(r"^\s*p\s*\.\s*", "", s)
    # A ">" between two nucleotide letters sometimes extracts as "4"
    # (Tzschach 2015 prints c.1937T>C, which comes back as "c.1937T4C").
    s = re.sub(r"(?<=[ACGT])4(?=[ACGT])", ">", s)
    # Papers space out the "c." prefix ("c. C1495T"), so collapse whitespace
    # before the ANNOVAR rewrite can see the pattern.
    s = re.sub(r"\bc\.\s+", "c.", s)
    # ANNOVAR-style descriptions put the reference base BEFORE the position
    # (`c.C121T` for `c.121C>T`), which no amount of separator-stripping will
    # reconcile. Yuan Table S2 and Chinen 2019 both use it. Rewrite to the
    # HGVS order before comparing.
    s = re.sub(r"\bc\.([ACGT])(\d+)([ACGT])\b", r"c.\2\1>\3", s)
    # keep only alphanumerics: every separator, operator and bracket a paper
    # might use is noise for the question being asked
    return re.sub(r"[^a-z0-9]", "", s.lower())


def read_any(path):
    """Text out of a pdf, xlsx or docx; empty string if unreadable."""
    try:
        if path.lower().endswith(".pdf"):
            rd = pypdf.PdfReader(path)
            return "\n".join(p.extract_text() or "" for p in rd.pages)
        if path.lower().endswith((".xlsx", ".xlsm")):
            import openpyxl
            wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
            return "\n".join(str(c) for ws in wb for row in ws.iter_rows()
                              for c in (x.value for x in row) if c is not None)
        if path.lower().endswith(".docx"):
            import docx
            d = docx.Document(path)
            out = [p.text for p in d.paragraphs]
            for tb in d.tables:
                for r in tb.rows:
                    out += [c.text for c in r.cells]
            return "\n".join(out)
    except Exception as e:
        print(f"    (could not read {os.path.basename(path)}: {e})", file=sys.stderr)
    return ""


pdfs = {os.path.basename(p)[:-4] for p in glob.glob(f"{SRC}/*.pdf")}


def find_pdf(src):
    if src in PDF_OVERRIDE:
        return PDF_OVERRIDE[src]
    m = re.match(r"([A-Za-z]+)_?(\d{4})", src)
    if not m:
        return None
    a, y = m.group(1).lower(), m.group(2)
    for p in pdfs:
        pl = p.lower().replace("_", "").replace("-", "")
        if a[:6] in pl and y in pl:
            return p
    return None


cur = {r["variant_key"] for r in
       csv.DictReader(open(f"{R}/data/processed/smc1a_variants_all.tsv"),
                      delimiter="\t")}
obs = [o for o in csv.DictReader(
    open(f"{R}/data/processed/smc1a_observations.tsv"), delimiter="\t")
    if o["variant_key_resolved"] in cur]

by_src = {}
for o in obs:
    by_src.setdefault(o["source"], []).append(o)

cache = {}
print(f"{'source':<22}{'vars':>5} {'cDNA':>9} {'protein':>9}  notes")
print("-" * 78)
problems = []
for src in sorted(by_src, key=lambda s: -len({o['variant_key_resolved'] for o in by_src[s]})):
    rows = by_src[src]
    nvar = len({o["variant_key_resolved"] for o in rows})
    if src in NO_PAPER:
        print(f"{src:<22}{nvar:>5} {'-':>9} {'-':>9}  database/cohort, no paper")
        continue
    pdf = find_pdf(src)
    if not pdf:
        print(f"{src:<22}{nvar:>5} {'?':>9} {'?':>9}  *** NO PDF FOUND ***")
        problems.append((src, "no pdf", []))
        continue
    key = (pdf, src)
    if key not in cache:
        try:
            rd = pypdf.PdfReader(f"{SRC}/{pdf}.pdf")
            txt = "\n".join(p.extract_text() or "" for p in rd.pages)
            for extra in SUPPLEMENTARY.get(src, []):
                txt += "\n" + read_any(f"{SRC}/{extra}")
            cache[key] = norm(txt)
        except Exception as e:
            cache[key] = ""
            print(f"{src:<22}{nvar:>5}   PDF READ FAILED: {e}")
            continue
    t = cache[key]
    cs = sorted({o["hgvs_c_input"] for o in rows if o.get("hgvs_c_input")})
    ps = sorted({o["hgvs_p_input"] for o in rows if o.get("hgvs_p_input")})
    mc = [v for v in cs if norm(v) not in t]
    mp = [v for v in ps if norm(v) not in t]
    note = ""
    if mc or mp:
        note = f"MISS c={len(mc)} p={len(mp)}"
        problems.append((src, pdf, {"cdna": mc, "protein": mp}))
    print(f"{src:<22}{nvar:>5} {len(cs)-len(mc):>4}/{len(cs):<4} "
          f"{len(ps)-len(mp):>4}/{len(ps):<4}  {note}")

print()
print("=" * 78)
print("SOURCES WITH MISSES -- for inspection")
print("=" * 78)
for src, pdf, m in problems:
    if not isinstance(m, dict):
        print(f"\n### {src}: {pdf}")
        continue
    print(f"\n### {src}  (pdf: {pdf})")
    if m["cdna"]:
        print(f"   cDNA not found   : {m['cdna']}")
    if m["protein"]:
        print(f"   protein not found: {m['protein']}")
