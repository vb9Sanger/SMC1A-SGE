#!/usr/bin/env python3
"""
07_make_data_dictionary.py

Generate `docs/DATA_DICTIONARY.md` from `code/smc1a_schema.py` and the actual
delivered files.

The dictionary is generated rather than hand-written so that it cannot drift
out of step with the tables: the same schema drives the column ordering in
05_build_registry.py. Every documented column is checked against the real
file header, and any mismatch in either direction is reported as a defect
rather than quietly ignored.

For each column the dictionary also reports the observed fill rate and, for
low-cardinality fields, the observed value distribution -- so the reader can
see what a column actually contains, not only what it is meant to contain.

Usage
-----
    ../.venv/bin/python 07_make_data_dictionary.py
"""

from __future__ import annotations

import csv
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smc1a_lib as L  # noqa: E402
import smc1a_schema as S  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed")
OUT = os.path.join(ROOT, "docs", "DATA_DICTIONARY.md")

MAX_VALUES = 12       # show the value distribution up to this cardinality
MAX_EXAMPLE = 60      # truncate example values to this many characters


def load(path):
    with open(path, encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        return r.fieldnames or [], list(r)


def summarise(rows, col):
    """(fill_rate_string, values_or_example_string)"""
    vals = [L.clean(r.get(col)) for r in rows]
    filled = [v for v in vals if v]
    if not rows:
        return "-", ""
    pct = f"{len(filled)}/{len(rows)} ({len(filled) / len(rows):.0%})"
    if not filled:
        return pct, "*(always empty)*"
    uniq = Counter(filled)
    if len(uniq) <= MAX_VALUES:
        return pct, ", ".join(
            f"`{v}` ({n})" for v, n in uniq.most_common())
    ex = sorted(uniq, key=lambda v: -uniq[v])[0]
    if len(ex) > MAX_EXAMPLE:
        ex = ex[:MAX_EXAMPLE] + "..."
    return pct, f"{len(uniq)} distinct; e.g. `{ex}`"


def main():
    out = []
    problems = []

    out.append("# Data dictionary")
    out.append("")
    out.append("**Generated** by `code/07_make_data_dictionary.py` from "
               "`code/smc1a_schema.py` and the delivered files. Do not edit by "
               "hand -- edit the schema and regenerate, so that the column "
               "order of the tables and this document stay in step.")
    out.append("")
    out.append("## Conventions that apply throughout")
    out.append("")
    out.append("| | |")
    out.append("|---|---|")
    out.append(f"| Reference transcript | `{L.TRANSCRIPT}` (MANE Select) = "
               f"`{L.REFSEQ_MANE}` |")
    out.append(f"| Reference protein | `{L.PROTEIN}` (1233 aa) |")
    out.append(f"| Genome build | {L.ASSEMBLY} (`{L.GENOMIC_ACC}`) |")
    out.append(f"| Gene strand | **minus** |")
    out.append("| Coordinates and alleles | **plus-strand** GRCh38 throughout. "
               "Matches VCF and the SGE oligo names, so `variant_key` joins "
               "directly to the screen results. Because *SMC1A* is on the "
               "minus strand, `ref`/`alt` are the **complement** of the bases "
               "named in the `c.` description (e.g. `c.1951G>A` is "
               "`chrX:53405352:C:T`). |")
    out.append("| `c.` numbering | Portable across NM_006306.1-.4 and MANE "
               "(verified). **Not** portable from NM_001281463.1, a different "
               "isoform offset by +66 nt / +22 codons. |")
    out.append("| List fields | `\\|`-separated, except `pmids`, "
               "`consequence_terms` and `hpo_ids` which are `;`-separated. |")
    out.append("| Booleans | The strings `True` / `False`. |")
    out.append("| Empty | Empty string means not recorded. Note `?` is a "
               "*value* in LOVD disease codes (\"unclassified / mixed\"), not a "
               "null. |")
    out.append("")
    out.append("### The three tables")
    out.append("")
    out.append("| File | Grain | Use it for |")
    out.append("|---|---|---|")
    out.append("| `smc1a_observations.tsv` | one row per **(variant x source x "
               "individual/report)** | auditing any variant-level value back to "
               "its source record |")
    out.append("| `smc1a_variants_all.tsv` | one row per **variant**, every "
               "variant encountered | seeing what was found and why anything "
               "was excluded |")
    out.append("| `smc1a_variants_curated.tsv` | one row per **variant**, "
               "high-confidence subset | the set to intersect with the screen |")
    out.append("")
    out.append("The observation ledger is primary and lossless; the two variant "
               "tables are derived from it. `variant_key_resolved` in the "
               "ledger is the foreign key to `variant_key` in the variant "
               "tables.")
    out.append("")
    out.append("### Counting rule")
    out.append("")
    out.append("Use **`n_independent_probands`** (one per family) for any "
               "recurrence or hotspot inference -- never `n_observations`, "
               "which is inflated by re-reporting. Always check "
               "`proband_count_is_upper_bound` first: it is `True` where "
               "independence is not established.")
    out.append("")

    for fname, schema in S.TABLES.items():
        path = os.path.join(PROC, fname)
        if not os.path.exists(path):
            problems.append(f"{fname}: file not found")
            continue
        header, rows = load(path)
        out.append("---")
        out.append("")
        out.append(f"## `{fname}`")
        out.append("")
        out.append(f"{len(rows)} rows, {len(header)} columns.")
        out.append("")

        documented = {c for _t, _n, cols in schema for c, _ty, _d in cols}
        missing = [c for c in header if c not in documented]
        absent = [c for c in documented if c not in header]
        if missing:
            problems.append(f"{fname}: columns present but undocumented: "
                            + ", ".join(missing))
        if absent:
            problems.append(f"{fname}: columns documented but absent: "
                            + ", ".join(absent))

        for title, note, cols in schema:
            cols_here = [(c, t, d) for c, t, d in cols if c in header]
            if not cols_here:
                continue
            out.append(f"### {title}")
            out.append("")
            if note:
                out.append(note)
                out.append("")
            out.append("| column | type | filled | description | observed values |")
            out.append("|---|---|---|---|---|")
            for c, t, d in cols_here:
                pct, vals = summarise(rows, c)
                d = d.replace("|", "\\|")
                vals = vals.replace("|", "\\|")
                out.append(f"| `{c}` | {t} | {pct} | {d} | {vals} |")
            out.append("")

    out.append("---")
    out.append("")
    out.append("## Schema integrity")
    out.append("")
    if problems:
        out.append("**Defects found:**")
        out.append("")
        for p in problems:
            out.append(f"- {p}")
    else:
        out.append("Every column in every delivered file is documented here, "
                   "and every documented column exists in its file. No "
                   "discrepancies.")
    out.append("")

    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")

    print(f"Wrote docs/DATA_DICTIONARY.md "
          f"({len(out)} lines, {len(S.TABLES)} tables)")
    if problems:
        print("\nSchema problems:")
        for p in problems:
            print("  -", p)
        sys.exit(1)
    print("Schema integrity: OK")


if __name__ == "__main__":
    main()
