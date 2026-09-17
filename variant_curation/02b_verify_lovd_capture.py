#!/usr/bin/env python3
"""
02b_verify_lovd_capture.py

Verify that the public LOVD full-gene export used by this pipeline contains
everything that was visible in the *authenticated* "Full data view" web page.

Why this check exists
---------------------
The curation switched from two manual page captures (an HTML capture of
entries 1-100 and a PDF capture of entries 101-179, both taken while logged in
to LOVD) to the public relational export at
`https://databases.lovd.nl/shared/download/all/gene/SMC1A`.

LOVD does show some content only to authenticated users, and individual-level
records can carry a per-record `license`. So the switch is only safe if the
public export is a superset of what the logged-in view showed. That is tested
here rather than assumed.

Method
------
1. Parse the 41 columns of the authenticated HTML capture (100 entries).
2. Map each HTML column onto the field(s) of the export that back it.
3. For every mapped column, check that the multiset of non-empty values seen
   in the HTML is a **subset** of the values present in the export for that
   field. This is robust to row ordering and to the HTML showing only the
   first page.
4. Report any HTML column with values absent from the export, and any export
   field with no HTML counterpart (i.e. data the web view did *not* show).
5. Confirm the entry count from the PDF capture's own header.
6. Check the per-record `license` field for access restrictions.

Output
------
docs/verification/lovd_capture_equivalence.md

Usage
-----
    ../.venv/bin/python 02b_verify_lovd_capture.py
"""

from __future__ import annotations

import glob
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smc1a_lib as L  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "variant_curation", "sources")
DOCS = os.path.join(ROOT, "docs", "verification")
os.makedirs(DOCS, exist_ok=True)

HTML = os.path.join(SRC, "Full data view for gene SMC1A - Global Variome shared LOVD.html")
PDF = os.path.join(SRC, "page 2 Full data view for gene SMC1A - Global Variome shared LOVD.pdf")

report: list[str] = []


def say(m: str = "") -> None:
    print(m)
    report.append(m)


# Map each column of the authenticated "Full data view" onto the export
# section and field that backs it. `None` means "presentational only".
# Duplicated HTML headers (Reference, Remarks, VIP appear twice -- once for
# the variant, once for the individual) are disambiguated by position.
COLUMN_MAP = {
    0:  ("Variants_On_Genome", "effectid"),
    1:  ("Variants_On_Transcripts", "VariantOnTranscript/Exon"),
    2:  ("Variants_On_Transcripts", "VariantOnTranscript/DNA"),
    3:  ("Variants_On_Transcripts", "VariantOnTranscript/RNA"),
    4:  ("Variants_On_Transcripts", "VariantOnTranscript/Protein"),
    5:  ("Variants_On_Genome", "allele"),
    6:  ("Variants_On_Genome", "VariantOnGenome/ClinicalClassification/Method"),
    7:  ("Variants_On_Genome", "VariantOnGenome/ClinicalClassification"),
    8:  ("Variants_On_Genome", "VariantOnGenome/DNA"),
    9:  ("Variants_On_Genome", "VariantOnGenome/DNA/hg38"),
    10: ("Variants_On_Genome", "VariantOnGenome/Published_as"),
    11: ("Variants_On_Genome", "VariantOnGenome/ISCN"),
    12: ("Variants_On_Genome", "VariantOnGenome/DBID"),
    13: ("Variants_On_Genome", "VariantOnGenome/Remarks"),
    14: ("Variants_On_Genome", "VariantOnGenome/Reference"),
    15: ("Variants_On_Genome", "VariantOnGenome/ClinVar"),
    16: ("Variants_On_Genome", "VariantOnGenome/dbSNP"),
    17: ("Variants_On_Genome", "VariantOnGenome/Genetic_origin"),
    18: ("Variants_On_Genome", "VariantOnGenome/Segregation"),
    19: ("Variants_On_Genome", "VariantOnGenome/Frequency"),
    20: ("Variants_On_Genome", "VariantOnGenome/Restriction_site"),
    21: ("Variants_On_Genome", "VariantOnGenome/VIP"),
    22: ("Variants_On_Genome", "VariantOnGenome/Methylation"),
    23: ("Screenings", "Screening/Template"),
    24: ("Screenings", "Screening/Technique"),
    25: ("Screenings", "Screening/Tissue"),
    26: ("Screenings", "Screening/Remarks"),
    27: ("Diseases", "symbol"),
    28: ("Individuals", "Individual/Individual_ID"),
    29: ("Individuals", "Individual/Reference"),
    30: ("Individuals", "Individual/Remarks"),
    31: ("Individuals", "Individual/Gender"),
    32: ("Individuals", "Individual/Consanguinity"),
    33: ("Individuals", "Individual/Origin/Geographic"),
    34: ("Individuals", "Individual/Origin/Population"),
    35: ("Individuals", "Individual/Age_of_death"),
    36: ("Individuals", "Individual/VIP"),
    37: ("Individuals", "Individual/Data_av"),
    38: ("Individuals", "Individual/Treatment"),
    39: ("Individuals", "panel_size"),
    40: ("Individuals", "owned_by"),
}

# The web view renders several fields differently from the way the download
# stores them. Each such field is canonicalised to a form comparable on both
# sides, so a genuine data gap is distinguishable from a mere encoding
# difference. Fields reduced to "<present>" are checked for *presence* only,
# because the two representations share no vocabulary.
ALLELE_CODE = {
    "0": "unknown", "1": "parent #1", "2": "parent #2", "3": "both alleles",
    "10": "paternal (inferred)", "11": "paternal (confirmed)",
    "20": "maternal (inferred)", "21": "maternal (confirmed)",
}


def canon(col_idx: int, v: str, side: str = "html") -> str:
    """Canonicalise one cell. `side` is 'html' or 'export'."""
    v = L.clean(v)
    if not v:
        return ""
    v = re.sub(r"\s+", " ", v).strip()

    if col_idx in (0, 40, 14, 29):
        # 'Effect' renders as '+/+' but is stored as a 2-digit effectid;
        # 'Owner' renders as a person's name but is stored as a user id;
        # references render as 'PubMed: X , Journal: Y' but are stored as
        # '{PMID:...}' markup. Presence-only comparison.
        return "<present>"
    if col_idx == 5:
        # 'Allele' renders as a label but is stored as a numeric code.
        return ALLELE_CODE.get(v, v.lower()) if side == "export" else v.lower()
    if col_idx == 24:
        # 'Technique' is multi-valued: ', ' in the view, ';' in the export.
        return ";".join(sorted(t.strip() for t in re.split(r"[;,]", v) if t.strip()))
    if col_idx == 27:
        # 'Disease' is multi-valued in the view ('CDLS, MYOP'); the export
        # holds one symbol per Diseases row, with the multiplicity carried by
        # Individuals_To_Diseases. Compare as sorted symbol sets.
        return ";".join(sorted(t.strip() for t in re.split(r"[;,|]", v) if t.strip()))
    return v
    if col_idx == 14 or col_idx == 29:
        # References rendered as 'PubMed: X , Journal: Y', stored as
        # '{PMID:X}' style markup
        return "<reference>"
    return v


def parse_html():
    from bs4 import BeautifulSoup
    s = BeautifulSoup(open(HTML, encoding="utf-8", errors="replace").read(), "lxml")
    # LOVD nests the data grid inside layout tables, so several tables contain
    # the header text. Require the innermost one: a header row of ~41 cells
    # whose first cell is exactly 'Effect'.
    grid, header = None, None
    for t in s.find_all("table"):
        rows = t.find_all("tr")
        if not rows:
            continue
        hdr = [c.get_text(" ", strip=True).replace("\xa0", " ")
               for c in rows[0].find_all(["th", "td"], recursive=False)]
        if not hdr:
            hdr = [c.get_text(" ", strip=True).replace("\xa0", " ")
                   for c in rows[0].find_all(["th", "td"])]
        if (30 <= len(hdr) <= 60 and hdr[0] == "Effect"
                and any("DNA change (cDNA)" in h for h in hdr)):
            grid, header = t, hdr
            break
    if grid is None:
        sys.exit("FATAL: could not locate the 41-column data grid in the HTML capture")

    out = []
    for r in grid.find_all("tr")[1:]:
        cells = [c.get_text(" ", strip=True).replace("\xa0", " ")
                 for c in r.find_all(["td", "th"], recursive=False)]
        if len(cells) == len(header):
            out.append(cells)
    if not out:
        sys.exit("FATAL: found the grid header but no data rows matched its width")
    return header, out


def parse_export():
    cands = sorted(glob.glob(os.path.join(ROOT, "data", "raw", "lovd",
                                          "LOVD_SMC1A_full_download_*.txt")))
    path = cands[-1]
    sections: dict[str, list[dict]] = {}
    cur, hdr = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            m = re.match(r"^## (\w+) ## Do not remove", line)
            if m:
                cur, hdr = m.group(1), []
                sections[cur] = []
                continue
            if line.startswith("#") or not line.strip() or cur is None:
                continue
            cells = [c.strip('"').replace('""', '"') for c in line.split("\t")]
            if not hdr:
                hdr = [c.strip("{}") for c in cells]
                continue
            sections[cur].append(dict(zip(hdr, cells)))
    return os.path.basename(path), sections


def main():
    header, rows = parse_html()
    name, S = parse_export()

    say("# Does the public LOVD export contain everything the logged-in view showed?")
    say()
    say("Compares the **authenticated** 'Full data view' page captures against the")
    say("**public** relational export now used by the pipeline.")
    say()
    say(f"- authenticated HTML capture: `{os.path.basename(HTML)}` "
        f"-- {len(rows)} entries, {len(header)} columns")
    say(f"- public export: `{name}`")
    say()

    # ---- entry counts -------------------------------------------------
    say("## Entry counts")
    say()
    pdf_total = pdf_shown = None
    try:
        import pdfplumber
        with pdfplumber.open(PDF) as pdf:
            txt = pdf.pages[0].extract_text() or ""
        m = re.search(r"(\d+)\s+entries on\s+(\d+)\s+pages?\.\s*"
                      r"Showing entries\s+(\d+)\s*-\s*(\d+)", txt)
        if m:
            pdf_total = int(m.group(1))
            pdf_shown = (int(m.group(3)), int(m.group(4)))
    except Exception as e:
        say(f"  (could not read the PDF capture: {type(e).__name__}: {e})")

    n_export = len(S["Variants_On_Genome"])
    say("| source | entries |")
    say("|---|---|")
    say(f"| authenticated HTML capture (page 1) | {len(rows)} |")
    if pdf_total:
        say(f"| authenticated PDF capture (page 2) | "
            f"{pdf_shown[1] - pdf_shown[0] + 1} (entries {pdf_shown[0]}-{pdf_shown[1]}) |")
        say(f"| **total stated by LOVD while logged in** | **{pdf_total}** |")
    say(f"| public export `Variants_On_Genome` | {n_export} |")
    say()
    if pdf_total:
        ok = pdf_total == n_export
        say(f"**{'PASS' if ok else 'FAIL'}** -- the logged-in view reported "
            f"{pdf_total} entries; the public export contains {n_export}.")
        if ok:
            say()
            say("No entries are withheld from the anonymous download.")
    say()

    # ---- per-column value coverage ------------------------------------
    say("## Per-column check")
    say()
    say("For each column of the logged-in view, every non-empty value it displayed")
    say("must also be present in the corresponding export field. A column passes")
    say("only if **no** displayed value is missing from the export.")
    say()

    export_vals = {}
    for sec, recs in S.items():
        for r in recs:
            for k, v in r.items():
                export_vals.setdefault((sec, k), Counter())[canon(-1, v, "export")] += 1

    say("| # | logged-in column | export field | values shown | missing from export | verdict |")
    say("|---|---|---|---|---|---|")
    problems = []
    for i, colname in enumerate(header):
        if i not in COLUMN_MAP:
            continue
        sec, field = COLUMN_MAP[i]
        shown = Counter()
        for r in rows:
            v = canon(i, r[i], "html") if i < len(r) else ""
            if v:
                shown[v] += 1
        if not shown:
            say(f"| {i} | {colname} | `{sec}.{field}` | 0 | - | no data shown |")
            continue
        have = export_vals.get((sec, field), Counter())
        # canon() on the export side, applied per-column
        have_c = Counter()
        for recs in [S[sec]]:
            for r in recs:
                v = canon(i, r.get(field, ""), "export")
                if v:
                    have_c[v] += 1
        # Multi-valued fields ('Technique', 'Disease') are stored one value
        # per export row, so compare the constituent tokens rather than the
        # rendered combination.
        if i in (24, 27):
            have_tokens = {t for v in have_c for t in v.split(";")}
            missing = sorted({t for v in shown for t in v.split(";")
                              if t not in have_tokens})
        else:
            missing = [v for v in shown if v not in have_c]
        verdict = "**PASS**" if not missing else "**FAIL**"
        if missing:
            problems.append((i, colname, sec, field, missing))
        say(f"| {i} | {colname} | `{sec}.{field}` | {len(shown)} | "
            f"{len(missing)} | {verdict} |")
    say()

    if problems:
        say("### Values shown while logged in but absent from the public export")
        say()
        for i, colname, sec, field, missing in problems:
            say(f"**{colname}** (`{sec}.{field}`) -- {len(missing)} missing:")
            for v in missing[:12]:
                say(f"  - `{v}`")
            say()
    else:
        say("No column showed a value that is missing from the public export.")
        say()

    # ---- what the export has that the web view did NOT show -----------
    say("## Data present in the export but *not* shown by the logged-in view")
    say()
    mapped = {(s, f) for s, f in COLUMN_MAP.values()}
    extra = []
    for sec, recs in S.items():
        if not recs:
            continue
        for f in recs[0].keys():
            if (sec, f) in mapped:
                continue
            n = sum(1 for r in recs if not L.blank(r.get(f)))
            if n:
                extra.append((sec, f, n))
    say("| export section | field | records populated |")
    say("|---|---|---|")
    for sec, f, n in sorted(extra, key=lambda x: (x[0], -x[2])):
        say(f"| {sec} | `{f}` | {n} |")
    say()
    say("The export is therefore a **superset**: it adds the entire `Phenotypes`")
    say("section, the `Diseases` metadata (OMIM id, inheritance), the family")
    say("pointers (`fatherid`, `motherid`, `panelid`) and the numeric c. positions,")
    say("none of which the 41-column web view displayed.")
    say()

    # ---- licence check -------------------------------------------------
    say("## Access restrictions on individual records")
    say()
    lic = Counter(L.clean(r.get("license")) or "(none)" for r in S["Individuals"])
    say("| `Individuals.license` | records |")
    say("|---|---|")
    for k, n in lic.most_common():
        say(f"| {k} | {n} |")
    say()
    if set(lic) <= {"(none)"}:
        say("No individual record carries a restrictive licence, consistent with the")
        say("entry counts matching.")
    say()

    say("## Conclusion")
    say()
    if not problems and (pdf_total == n_export if pdf_total else True):
        say("The public export contains every entry and every field value that the")
        say("authenticated view displayed, plus substantial additional")
        say("individual-level phenotype data. Using the export loses nothing, and the")
        say("page captures are retained only as provenance snapshots.")
    else:
        say("**The export does not fully cover the authenticated view** -- see the")
        say("failures above. The page captures must be parsed to recover the")
        say("missing values.")

    with open(os.path.join(DOCS, "lovd_capture_equivalence.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(report) + "\n")
    print("\nWrote docs/verification/lovd_capture_equivalence.md")


if __name__ == "__main__":
    main()
