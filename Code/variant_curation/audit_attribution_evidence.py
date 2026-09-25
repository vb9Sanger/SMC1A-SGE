#!/usr/bin/env python3
"""Evidence audit v2 of the curated SMC1A disease attributions.

v1 flagged cross-arm text from free prose and over-fired: the main driver was
Huisman 2017's deliberately undecided label "SMC1A-related disorder (CdLS
and/or epilepsy)", which the curation already declines to take an attribution
from. v2 uses the harmonised disease code per observation instead, so a
"contradiction" means another source positively claimed the OTHER arm.
"""
import re
import pandas as pd

ROOT = "/Users/vb9/Documents/Claude/SMC1A/thesis_workflow/curation_attempt2"
SC = ("/private/tmp/claude-503/-Users-vb9-Documents-Claude-SMC1A-thesis-workflow/"
      "3868e5c0-820a-4d63-9a8f-9da30bba6c9e/scratchpad")

cur = pd.read_csv(f"{ROOT}/data/processed/smc1a_variants_curated.tsv",
                  sep="\t", dtype=str, low_memory=False)
obs = pd.read_csv(f"{ROOT}/data/processed/smc1a_observations.tsv",
                  sep="\t", dtype=str, low_memory=False)

CDLS = {"CdLS", "CdLS2"}
DEE = {"DEE85"}
obs["_unsuppressed"] = obs.disease_attribution_suppressed.isna()
obs["_cdls"] = obs.disease_harmonised.isin(CDLS) & obs._unsuppressed
obs["_dee"] = obs.disease_harmonised.isin(DEE) & obs._unsuppressed

INFERRED = re.compile(r"assigned from cohort phenotype|assigned from the patient's DECIPHER", re.I)
GENOTYPE_ASC = re.compile(r"genetically-ascertained|genetically.defined", re.I)
CAVEATED = re.compile(r"\bNOTE\b|weakens the attribution|also carries a NIPBL|carries BOTH", re.I)

rows = []
for _, v in cur.iterrows():
    arm = v.curation_group
    vo = obs[obs.variant_key_resolved == v.variant_key]
    support = vo[vo._cdls] if arm == "CdLS_pathogenic" else vo[vo._dee]
    oppose = vo[vo._dee] if arm == "CdLS_pathogenic" else vo[vo._cdls]

    basis = " || ".join(support.attribution_basis.dropna().astype(str))
    detail = support.phenotype_detail.dropna().astype(str)
    hpo = vo.n_hpo_terms.dropna()
    rows.append(dict(
        hgvs_c=v.hgvs_c_mane, hgvs_p=v.hgvs_p_mane, arm=arm,
        consequence=v.consequence_class, variant_key=v.variant_key,
        n_support_src=len(set(support.source)),
        support_src="|".join(sorted(set(support.source))),
        n_oppose_src=len(set(oppose.source)),
        oppose_src="|".join(sorted(set(oppose.source))),
        inferred=bool(INFERRED.search(basis)),
        genotype_asc=bool(GENOTYPE_ASC.search(basis)),
        caveated=bool(CAVEATED.search(basis)),
        percase_detail=bool((len(detail) and detail.str.len().max() > 40) or (len(hpo) and hpo.astype(float).max() > 0)),
        n_hpo=int(hpo.astype(float).max()) if len(hpo) else 0,
        episig=str(v.get("episignature_available")) == "True",
        n_probands=int(v.n_independent_probands),
        dup=str(v.get("possible_cross_source_duplicate")) == "True",
    ))

df = pd.DataFrame(rows)
df["contradicted"] = df.n_oppose_src > 0
# insufficient: a single supporting source whose claim is rule-inferred, from a
# genotype-ascertained cohort, self-caveated, or carries no per-case detail
df["insufficient"] = (df.n_support_src <= 1) & (
    df.inferred | df.genotype_asc | df.caveated | ~df.percase_detail)
df["flagged"] = df.contradicted | df.insufficient | (df.n_support_src == 0)
df.to_csv(f"{SC}/attribution_evidence_audit.tsv", sep="\t", index=False)

for arm, g in df.groupby("arm"):
    print(f"=== {arm} (n={len(g)}) ===")
    print(f"  0 supporting sources            : {(g.n_support_src==0).sum()}")
    print(f"  1 supporting source             : {(g.n_support_src==1).sum()}")
    print(f"  >=2 supporting sources          : {(g.n_support_src>=2).sum()}")
    print(f"  CONTRADICTED (other arm claimed): {g.contradicted.sum()}")
    print(f"  INSUFFICIENT                    : {g.insufficient.sum()}")
    print(f"    of which rule-inferred        : {(g.insufficient & g.inferred).sum()}")
    print(f"    of which genotype-ascertained : {(g.insufficient & g.genotype_asc).sum()}")
    print(f"    of which source self-caveated : {(g.insufficient & g.caveated).sum()}")
    print(f"    of which no per-case detail   : {(g.insufficient & ~g.percase_detail).sum()}")
    print(f"  ** TOTAL FLAGGED **             : {g.flagged.sum()}  "
          f"(clean: {(~g.flagged).sum()})")
    print()
