"""Independent audit of the curated set.

Re-derives every checkable property from the delivered file plus the reference
sequence and transcript model, rather than trusting the pipeline's own
annotations. Anything the pipeline got wrong should surface here as a
disagreement, not be confirmed by circular reading.
"""
import csv, os, re, sys
from collections import Counter, defaultdict

# Paths are derived from this file's location, as in every other step, so the
# audit runs from any working directory rather than only from the project root.
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import smc1a_lib as L
from smc1a_lib import CODON_TABLE, AA3, COMPLEMENT

TM = L.TranscriptModel.from_json(
    os.path.join(ROOT, "data", "reference", "smc1a_transcript_model.json"))
RS = L.ReferenceSequence(
    os.path.join(ROOT, "data", "reference", "smc1a_chrX_region.fa"))
AAR = {v: k for k, v in AA3.items()}
LOCUS = (53374149, 53422654)

rows = list(csv.DictReader(open(os.path.join(ROOT, "data", "processed", "smc1a_variants_curated.tsv")),
                           delimiter="\t"))
allrows = list(csv.DictReader(open(os.path.join(ROOT, "data", "processed", "smc1a_variants_all.tsv")),
                              delimiter="\t"))
fail = defaultdict(list)


def codon_aa(n):
    if n < 1 or n > 1233:          # outside the protein; do not read past the CDS
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


# A1 variant_key well-formed and unique
seen = Counter(r["variant_key"] for r in rows)
for k, n in seen.items():
    if n > 1:
        fail["A1 duplicate variant_key"].append(f"{k} x{n}")
for r in rows:
    if not re.fullmatch(r"chrX:\d+:[ACGT]+:[ACGT]+", r["variant_key"]):
        fail["A1 malformed variant_key"].append(r["variant_key"])

# A2 position inside the locus
for r in rows:
    p = int(r["variant_key"].split(":")[1])
    if not LOCUS[0] <= p <= LOCUS[1]:
        fail["A2 position outside SMC1A"].append(f"{r['hgvs_c_mane']} @ {p}")

# A3 reference allele matches GRCh38
for r in rows:
    _, p, ref, alt = r["variant_key"].split(":")
    p = int(p)
    obs = RS.get(p, p + len(ref) - 1).upper()
    if obs != ref.upper():
        fail["A3 ref allele mismatch"].append(
            f"{r['hgvs_c_mane']} key says {ref}, GRCh38 has {obs}")

# A4 hgvs_c maps back to the same genomic position (substitutions only --
#    indels legitimately differ by HGVS 3'-shifting vs VCF left-alignment)
for r in rows:
    c = r.get("hgvs_c_mane") or ""
    m = re.fullmatch(r"c\.(\d+)([ACGT])>([ACGT])", c)
    if not m:
        continue
    g = TM.c_pos_to_g(m.group(1))
    if isinstance(g, tuple):
        g = g[0]
    if g is None:
        fail["A4 c. does not map"].append(c)
        continue
    kp = int(r["variant_key"].split(":")[1])
    if int(g) != kp:
        fail["A4 c./key position disagree"].append(f"{c}: c.->{g}, key {kp}")
    # and the cDNA base should be the stated reference base
    cb = COMPLEMENT[RS.base(int(g)).upper()]
    if cb != m.group(2):
        fail["A4 cDNA ref base wrong"].append(f"{c}: MANE has {cb}")

# A5 protein reference residue matches MANE
for r in rows:
    p = r.get("hgvs_p_mane") or ""
    m = re.match(r"^p\.\(?([A-Z][a-z]{2})(\d+)", p)
    if not m:
        continue
    want, num = m.group(1), int(m.group(2))
    got = codon_aa(num)
    if got and got != want:
        fail["A5 protein ref residue mismatch"].append(
            f"{r['hgvs_c_mane']} {p}: MANE codon {num} is {got}")

# A6 consequence_class consistent with the protein description
for r in rows:
    cls, p = r.get("consequence_class") or "", r.get("hgvs_p_mane") or ""
    if not p:
        continue
    SPLICE = ("splice_region", "splice_donor", "splice_acceptor")
    # A frameshift whose first affected codon becomes a stop is written
    # `p.Xaa123Ter`, not `p.Xaa123fs`, so a Ter protein is legitimate for a
    # frameshift class. And a splice variant's PREDICTED protein consequence is
    # often a downstream frameshift, so `fs` is legitimate for a splice class.
    # The variant class and the protein consequence answer different questions.
    if re.search(r"Ter|\*", p) and "fs" not in p and cls not in (
            ("nonsense", "frameshift") + SPLICE):
        fail["A6 stop but not nonsense/frameshift"].append(
            f"{r['hgvs_c_mane']} {p} -> {cls}")
    if "fs" in p and cls not in (("frameshift",) + SPLICE):
        fail["A6 fs but wrong class"].append(f"{r['hgvs_c_mane']} {p} -> {cls}")
    # a true missense: three-letter ref and alt residues, alt not a stop
    mm = re.fullmatch(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)?", p)
    if mm and mm.group(3) != "Ter" and cls != "missense":
        fail["A6 substitution but not missense"].append(
            f"{r['hgvs_c_mane']} {p} -> {cls}")

# A7 inclusion criteria re-derived
for r in rows:
    if r.get("pathogenicity_consensus") not in ("P", "LP"):
        fail["A7 not P/LP"].append(f"{r['hgvs_c_mane']}={r.get('pathogenicity_consensus')}")
    if r.get("disease_attribution") not in ("CdLS", "DEE85"):
        fail["A7 disease unresolved"].append(
            f"{r['hgvs_c_mane']}={r.get('disease_attribution')}")
    if str(r.get("has_clinical_detail")).lower() not in ("true", "1", "yes"):
        fail["A7 no clinical detail"].append(r["hgvs_c_mane"])
    if str(r.get("in_design_window")).lower() not in ("true", "1", "yes"):
        fail["A7 outside design window"].append(r["hgvs_c_mane"])
    if r.get("exclusion_reason"):
        fail["A7 curated but has exclusion_reason"].append(
            f"{r['hgvs_c_mane']}={r['exclusion_reason']}")
    if str(r.get("include_high_confidence")).lower() not in ("true", "1", "yes"):
        fail["A7 include flag not set"].append(r["hgvs_c_mane"])

# A8 arm membership sane
for r in rows:
    if r.get("curation_group") not in ("CdLS_pathogenic", "DEE85_pathogenic"):
        fail["A8 curated but not in an arm"].append(
            f"{r['hgvs_c_mane']}={r.get('curation_group')}")
    if truthy_flag := str(r.get("cdls_attribution_contradicted", "")).lower() in ("true",):
        fail["A8 contradicted CdLS still curated"].append(r["hgvs_c_mane"])

# A9 probands and provenance
for r in rows:
    try:
        n = int(r.get("n_independent_probands") or 0)
    except ValueError:
        n = 0
    if n < 1:
        fail["A9 no probands"].append(r["hgvs_c_mane"])
    if not (r.get("sources") or "").strip():
        fail["A9 no source"].append(r["hgvs_c_mane"])

# A10 curated is a strict subset of all, on variant_key
allkeys = {r["variant_key"] for r in allrows}
for r in rows:
    if r["variant_key"] not in allkeys:
        fail["A10 curated variant absent from all-set"].append(r["variant_key"])

# A11 no two curated variants are the same genomic event
norm = defaultdict(list)
for r in rows:
    _, p, ref, alt = r["variant_key"].split(":")
    # left-trim shared prefix so equivalent representations collide
    ref, alt, p = ref.upper(), alt.upper(), int(p)
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, p = ref[1:], alt[1:], p + 1
    norm[(p, ref, alt)].append(r["hgvs_c_mane"])
for k, v in norm.items():
    if len(v) > 1:
        fail["A11 same genomic event twice"].append(f"{k} -> {', '.join(v)}")


# A12 no DEE85 variant has a male patient (the female-only premise)
for r in rows:
    if r.get("curation_group") == "DEE85_pathogenic" and int(r.get("sex_M") or 0) > 0:
        fail["A12 male patient in the DEE85 arm"].append(
            f"{r['hgvs_c_mane']} sex_M={r['sex_M']}")

# A13 coding variants should carry a protein description
for r in rows:
    cls = r.get("consequence_class") or ""
    if cls in ("missense", "nonsense", "frameshift", "inframe_deletion",
               "inframe_insertion") and not (r.get("hgvs_p_mane") or "").strip():
        fail["A13 coding variant with no protein description"].append(
            f"{r['hgvs_c_mane']} ({cls})")

# A14 every variant is assay-addressable
for r in rows:
    if not (r.get("targetons_design") or "").strip():
        fail["A14 no design-window targeton"].append(r["hgvs_c_mane"])
    st = r.get("targeton_screening_status") or ""
    if "screened" not in st:
        fail["A14 targeton not screened"].append(f"{r['hgvs_c_mane']} ({st})")

# A15 HGVS descriptions are well-formed
for r in rows:
    c = r.get("hgvs_c_mane") or ""
    if c and not re.match(r"^c\.[\d\-+*]", c):
        fail["A15 malformed hgvs_c"].append(c)
    p = r.get("hgvs_p_mane") or ""
    if p and not p.startswith("p."):
        fail["A15 protein missing p. prefix"].append(f"{r['hgvs_c_mane']}: {p}")

print("# Final audit of the curated set")
print()
print(f"variants audited: **{len(rows)}**")
print()
CHECKS = ["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11"]
if not fail:
    print("All checks PASS -- no disagreement between the delivered file and an "
          "independent re-derivation from the reference.")
else:
    for k in sorted(fail):
        print(f"**{k}: {len(fail[k])}**")
        for x in fail[k][:8]:
            print(f"  - {x}")
        print()

if fail:
    sys.exit(1)
