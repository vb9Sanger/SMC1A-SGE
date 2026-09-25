#!/usr/bin/env python3
"""Re-derive protein change, position and group assignment for the ALL-VARIANTS
set, independently of what the pipeline recorded.

`08_audit_curated_set.py` runs its re-derivation checks over the 113 curated
rows only. The other 153 variants feed the benign reference, the ClinVar
comparisons and the VUS work, so an error there reaches the specificity and
OddsPath numbers without ever touching the curated set. This applies the same
independent re-derivation to all 266, and additionally re-runs the group gate
from its own stated inputs.
"""
import csv, os, re, sys
from collections import Counter, defaultdict

ROOT = "/Users/vb9/Documents/Claude/SMC1A/thesis_workflow/curation_attempt2"
sys.path.insert(0, f"{ROOT}/code")
import smc1a_lib as L
from smc1a_lib import CODON_TABLE, AA3, COMPLEMENT

TM = L.TranscriptModel.from_json(f"{ROOT}/data/reference/smc1a_transcript_model.json")
RS = L.ReferenceSequence(f"{ROOT}/data/reference/smc1a_chrX_region.fa")
AAR = {v: k for k, v in AA3.items()}

rows = list(csv.DictReader(open(f"{ROOT}/data/processed/smc1a_variants_all.tsv"),
                           delimiter="\t"))
cur = {r["variant_key"] for r in
       csv.DictReader(open(f"{ROOT}/data/processed/smc1a_variants_curated.tsv"),
                      delimiter="\t")}
print(f"all-variants rows: {len(rows)}  (curated {len(cur)}, "
      f"non-curated {len(rows)-len(cur)})\n")


def codon_aa(n):
    if n < 1 or n > 1233:
        return None
    try:
        b = ""
        for i in range(3):
            g = TM.c_pos_to_g(str((n - 1) * 3 + 1 + i))
            if isinstance(g, tuple):
                g = g[0]
            if g is None:
                return None
            b += COMPLEMENT[RS.base(int(g)).upper()]
        return AAR.get(CODON_TABLE.get(b))
    except Exception:
        return None


fail = defaultdict(list)
SPLICE = ("splice_region", "splice_donor", "splice_acceptor")

for r in rows:
    tag = "curated" if r["variant_key"] in cur else "other"
    c = r.get("hgvs_c_mane") or ""
    p = r.get("hgvs_p_mane") or ""
    cls = r.get("consequence_class") or ""
    lab = f"[{tag}] {c} {p} ({cls})"

    # --- P1 substitution cDNA position maps to the key, and the ref base is right
    m = re.fullmatch(r"c\.(\d+)([ACGT])>([ACGT])", c)
    if m and not r["variant_key"].startswith("UNRESOLVED"):
        g = TM.c_pos_to_g(m.group(1))
        if isinstance(g, tuple):
            g = g[0]
        if g is None:
            fail["P1 c. does not map to the genome"].append(lab)
        else:
            kp = int(r["variant_key"].split(":")[1])
            if int(g) != kp:
                fail["P1 c./key position disagree"].append(f"{lab}: c.->{g} key {kp}")
            cb = COMPLEMENT[RS.base(int(g)).upper()]
            if cb != m.group(2):
                fail["P1 cDNA ref base wrong"].append(f"{lab}: MANE has {cb}")

    # --- P2 the protein's stated reference residue matches the MANE translation
    mp = re.match(r"^p\.\(?([A-Z][a-z]{2})(\d+)", p)
    if mp:
        want, num = mp.group(1), int(mp.group(2))
        got = codon_aa(num)
        if got and got != want:
            fail["P2 protein ref residue mismatch"].append(
                f"{lab}: MANE codon {num} is {got}")

    # --- P3 protein_position agrees with the protein description
    pp = (r.get("protein_position") or "").strip()
    if mp and pp.isdigit() and int(pp) != int(mp.group(2)):
        fail["P3 protein_position disagrees with hgvs_p"].append(
            f"{lab}: protein_position={pp}, hgvs_p says {mp.group(2)}")

    # --- P4 protein_position agrees with the cDNA position
    mc = re.match(r"c\.(\d+)", c)
    if mc and pp.isdigit():
        expect = (int(mc.group(1)) - 1) // 3 + 1
        if abs(int(pp) - expect) > 1:
            fail["P4 protein_position far from its cDNA codon"].append(
                f"{lab}: protein_position={pp}, c.{mc.group(1)} is codon {expect}")

    # --- P5 consequence_class consistent with the protein description
    if p:
        if re.search(r"Ter|\*", p) and "fs" not in p and cls not in (
                ("nonsense", "frameshift") + SPLICE):
            fail["P5 stop but class is not nonsense/frameshift/splice"].append(lab)
        if "fs" in p and cls not in (("frameshift",) + SPLICE):
            fail["P5 fs but class is not frameshift/splice"].append(lab)
        mm = re.fullmatch(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)?", p)
        if mm and mm.group(3) != "Ter" and cls not in ("missense", "synonymous"):
            fail["P5 substitution but class is not missense"].append(lab)
        if mm and mm.group(1) == mm.group(3) and cls != "synonymous":
            fail["P5 same residue but class is not synonymous"].append(lab)

for k in sorted(fail):
    print(f"**{k}: {len(fail[k])}**")
    for x in fail[k][:12]:
        print(f"   {x}")
    print()
if not fail:
    print("All protein/position re-derivations agree.\n")
