#!/usr/bin/env python3
"""
04b_ingest_literature.py

Turn the per-paper extraction files in `data/raw/literature/` into observation
rows in the common ledger schema.

Extraction model
----------------
Extraction from a PDF is a judgement task, not an automated one, so it is kept
deliberately separate from ingestion:

  * `data/raw/literature/<paper>_curated.tsv` is the **extraction**: what the
    paper says, transcribed with any hand corrections documented in
    DECISIONS_LOG. It is reviewable line by line against the PDF.
  * this script is the **ingestion**: it validates those files, normalises them
    onto the ledger schema, and derives the individual/family keys.

Nothing here interprets a paper; it only reshapes an extraction.

Individual identity across papers
---------------------------------
Every observation is keyed on the **primary** report of that patient, not on
the paper it was read from:

    individual_key = PUB:<CitationCompact>:<patient>      e.g. PUB:Baranano2022:P6

Consequently, when Bozarth 2023 re-reports Barañano 2022 patient P6 and
Barañano 2022 is later read directly, both produce the same key and collapse to
**one** individual in the registry. De-duplication is therefore structural
rather than a post-hoc cleanup step, and `is_rereport` records only which paper
a row was read from -- it is not what prevents double counting.

LOVD individuals whose submitter label encodes a PMID and patient id (D38) are
mapped into the same namespace by `05_build_registry.py`, so a publication
patient already deposited in LOVD also collapses.

Classification
--------------
DEE cohort papers typically assert pathogenicity without stating ACMG
criteria. Bozarth 2023, for instance, uses "pathogenic" 15 times but never
mentions ACMG, "likely pathogenic" or "VUS". Such an assertion is mapped
conservatively to **LP**, never P, with the basis recorded in
`classification_method_raw`, so an author's assertion is never silently
promoted to a formal classification. See OPEN_QUESTIONS Q16.

Outputs
-------
data/interim/literature_observations.tsv
docs/verification/literature_ingest.md

Usage
-----
    ../.venv/bin/python 04b_ingest_literature.py
"""

from __future__ import annotations

import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smc1a_lib as L  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw", "literature")
INTERIM = os.path.join(ROOT, "data", "interim")
DOCS = os.path.join(ROOT, "docs", "verification")

report: list[str] = []


def say(m: str = "") -> None:
    print(m)
    report.append(m)


# --------------------------------------------------------------------------
# Paper registry
# --------------------------------------------------------------------------
# One entry per extraction file. `pmid` is recorded only where it is known
# from an authoritative source (LOVD's citation markup, or the paper itself);
# it is never guessed, because a wrong PMID would create a false cross-source
# link. Citation strings are the join key instead.
PAPERS = {
    "bozarth2023_table3_curated.tsv": {
        "source": "Bozarth_2023",
        "citation": "Bozarth 2023",
        "pmid": "",                     # to be filled once confirmed
        "design": "primary case series with compiled literature review",
        "disease_reported": "SMC1A-DEE",
        "disease_harmonised": "DEE85",
        "disease_basis": ("all 41 variants presented as SMC1A-related "
                          "developmental and epileptic encephalopathy, with "
                          "seizure-onset age, speech, ID, walking and MRI per "
                          "patient"),
        "classification_raw": "pathogenic (author assertion; no ACMG criteria stated)",
        "classification": "LP",
        "classification_method": ("author assertion in a peer-reviewed DEE "
                                  "cohort; no ACMG criteria stated, so mapped "
                                  "conservatively to LP rather than P"),
        "table": "Table 3",
    },
    "elwan2021_case_curated.tsv": {
        "source": "Elwan_2021",
        "citation": "Elwan 2021",
        "pmid": "",
        "design": "primary single-case report",
        "disease_reported": "SMC1A truncating variant with late-onset clustering seizures",
        "disease_harmonised": "DEE85",
        "disease_basis": ("primary case report of a de novo heterozygous "
                          "truncating SMC1A variant presenting with clustering "
                          "seizures; unusual in that onset was age 12 and "
                          "development was normal into adulthood"),
        "classification_raw": "pathogenic (author assertion; absent from population databases)",
        "classification": "LP",
        "classification_method": ("author assertion; de novo, truncating, absent "
                                  "from population databases. No ACMG criteria "
                                  "stated, so mapped conservatively to LP"),
        "table": "case description",
    },
    "oguni2019_case_curated.tsv": {
        "source": "Oguni_2019",
        "citation": "Oguni 2019",
        "pmid": "",
        "design": "primary case report",
        "disease_reported": "SMC1A-related cluster seizures",
        "disease_harmonised": "DEE85",
        "disease_basis": ("primary report of periodic cluster seizures with "
                          "long-term video-EEG. NOTE: the variant is "
                          "maternally inherited, not de novo, which weakens "
                          "the attribution relative to the de novo cases"),
        "classification_raw": "reported as the causative variant (author assertion)",
        "classification": "VUS",
        "classification_method": ("maternally inherited missense in an X-linked "
                                  "dominant condition, with the mother's "
                                  "phenotype not established in the extracted "
                                  "text. Mapped to VUS rather than LP pending "
                                  "that check -- see Q19"),
        "table": "case description",
    },
    "jansen2015_table1_curated.tsv": {
        "source": "Jansen_2016",
        "citation": "Jansen 2016",
        "pmid": "",
        "design": "primary WGS study, 2 patients (multi-gene de novo table)",
        "disease_reported": "de novo loss-of-function SMC1A, early-onset epilepsy",
        "disease_harmonised": "DEE85",
        "disease_basis": ("primary WGS report of de novo loss-of-function SMC1A "
                          "variants in X-linked early-onset epilepsy; both rated "
                          "'Likely' clinically relevant by the authors"),
        "classification_raw": "de novo loss-of-function, rated clinically relevant",
        "classification": "LP",
        "classification_method": ("de novo LOF in an established gene, author-rated "
                                  "clinically relevant; no ACMG criteria stated, so "
                                  "mapped conservatively to LP"),
        "table": "Table 1",
    },
    "goldstein2015_table1_curated.tsv": {
        "source": "Goldstein_2015",
        "citation": "Goldstein 2015",
        "pmid": "",
        "design": "primary case series, 2 patients (multi-gene table)",
        "disease_reported": "SMC1A-related epileptic encephalopathy",
        "disease_harmonised": "DEE85",
        "disease_basis": ("primary report of two patients with novel heterozygous "
                          "SMC1A frameshift variants and epileptic "
                          "encephalopathy; both parents negative in each case"),
        "classification_raw": "novel heterozygous frameshift, both parents negative",
        "classification": "LP",
        "classification_method": ("de novo frameshift (parents negative); no ACMG "
                                  "criteria stated, so mapped conservatively to LP"),
        "table": "Table 1",
    },
    "lebrun2015_case_curated.tsv": {
        "pmid": '',
        "classification": 'LP',
        "classification_raw": 'pathogenic / novel de novo (author assertion; no ACMG criteria stated)',
        "classification_method": 'author assertion in a primary report; no ACMG criteria stated, so mapped conservatively to LP rather than P',
        "source": 'Lebrun_2015',
        "citation": 'Lebrun 2015',
        "design": 'primary single-case report with functional studies',
        "disease_reported": 'early-onset encephalopathy with epilepsy AND clinical CdLS',
        "disease_harmonised": 'DEE85',
        "disease_basis": "primary report titled 'Early-Onset Encephalopathy with Epilepsy Associated with a Novel Splice Site Mutation in SMC1A'. The patient carries BOTH a clinical CdLS diagnosis and severe early-onset epileptic encephalopathy, with RT-PCR evidence of aberrant splicing and reduced transcript",
        "table": 'case description',
    },
    "wenger2016_case_curated.tsv": {
        "pmid": '',
        "classification": 'LP',
        "classification_raw": 'pathogenic / novel de novo (author assertion; no ACMG criteria stated)',
        "classification_method": 'author assertion in a primary report; no ACMG criteria stated, so mapped conservatively to LP rather than P',
        "source": 'Wenger_2016',
        "citation": 'Wenger 2016',
        "design": 'primary single-case report',
        "disease_reported": 'CdLS with left ventricular non-compaction cardiomyopathy',
        "disease_harmonised": 'CdLS',
        "disease_basis": 'primary report of novel CdLS features -- left ventricular non-compaction cardiomyopathy, microform cleft lip, poor growth -- with a de novo in-frame SMC1A deletion',
        "table": 'case description',
    },
    "fang2020_case_curated.tsv": {
        "pmid": '',
        "classification": 'LP',
        "classification_raw": 'pathogenic / novel de novo (author assertion; no ACMG criteria stated)',
        "classification_method": 'author assertion in a primary report; no ACMG criteria stated, so mapped conservatively to LP rather than P',
        "source": 'Fang_2020',
        "citation": 'Fang 2021',
        "design": 'primary case report (letter)',
        "disease_reported": 'SMC1A-related epilepsy with CdLS features',
        "disease_harmonised": 'DEE85',
        "disease_basis": 'primary report of a de novo nonsense variant with frequent daily seizures, EEG abnormalities and skeletal/cardiac findings',
        "table": 'case description',
    },
    "gorman_case_curated.tsv": {
        "pmid": '',
        "classification": 'LP',
        "classification_raw": 'pathogenic / novel de novo (author assertion; no ACMG criteria stated)',
        "classification_method": 'author assertion in a primary report; no ACMG criteria stated, so mapped conservatively to LP rather than P',
        "source": 'Gorman_2017',
        "citation": 'Gorman 2017',
        "design": 'primary case report',
        "disease_reported": 'SMC1A-related early-onset epilepsy',
        "disease_harmonised": 'DEE85',
        "disease_basis": 'primary report of a novel de novo frameshift with focal seizures from 17 weeks and multifocal EEG abnormalities',
        "table": 'case description',
    },
    "mannini2010_table1_curated.tsv": {
        "source": "Mannini_2010",
        "citation": "Mannini 2010",
        "pmid": "",
        "design": "SMC1A mutation update / review with functional work",
        "disease_reported": "Cornelia de Lange syndrome",
        "disease_harmonised": "CdLS",
        "disease_basis": ("Table 1 is explicitly titled 'Mutational spectrum of "
                          "the SMC1A gene in CdLS'. Tables 2 (colorectal cancer, "
                          "somatic) and 3 (polymorphisms) are excluded"),
        "classification_raw": "reported in the CdLS mutational spectrum",
        "classification": "LP",
        "classification_method": ("listed in a curated CdLS mutational spectrum "
                                  "with a first-description citation per variant; "
                                  "no ACMG criteria stated, so mapped "
                                  "conservatively to LP"),
        "table": "Table 1",
    },
    "jang2015_case_curated.tsv": {
        "source": "Jang_2015",
        "citation": "Jang 2015",
        "pmid": "",
        "design": "primary family report (3 generations, 4 affected)",
        "disease_reported": "Cornelia de Lange syndrome",
        "disease_harmonised": "CdLS",
        "disease_basis": ("primary report titled 'Novel Pathogenic Variant "
                          "(c.3178G>A) in the SMC1A Gene in a Family With "
                          "Cornelia de Lange Syndrome'; four affected relatives "
                          "across three generations, mild phenotype"),
        "classification_raw": "novel pathogenic variant (author assertion)",
        "classification": "LP",
        "classification_method": ("author assertion, segregating with a mild CdLS "
                                  "phenotype across three generations; no ACMG "
                                  "criteria stated, so mapped conservatively to LP"),
        "table": "case description and Table 1",
    },
    "arefeshghi2020_episig_curated.tsv": {
        "source": "ArefEshghi_2020",
        "citation": "Aref-Eshghi 2020",
        "pmid": "",
        "design": "episignature cohort study (42 Mendelian NDDs)",
        "disease_reported": "Cornelia de Lange syndrome (by episignature)",
        "disease_harmonised": "CdLS",
        "disease_basis": ("attribution by DNA methylation episignature: the "
                          "subject clustered with other CdLS cases on "
                          "genome-wide methylation analysis, correcting an "
                          "initial misclassification. This is an orthogonal, "
                          "assay-based attribution rather than a clinical "
                          "impression"),
        "classification_raw": "rare variant in a CdLS-related gene, episignature-positive",
        "classification": "LP",
        "classification_method": ("episignature-positive for CdLS; no ACMG "
                                  "criteria stated, so mapped conservatively to LP"),
        "table": "Figure 5 / text",
    },
    "astorino2025_table2_curated.tsv": {
        "source": "Astorino_2025",
        "citation": "Astorino 2025",
        "pmid": "",
        "design": "REVIEW -- overlay source only, never a source of variants",
        "disease_reported": "SMC1A-related DEE",
        "disease_harmonised": "DEE85",
        "disease_basis": ("Table 2 compiles seizure onset, XCI, speech, ID, "
                          "walking and MRI from THREE case series, and they do "
                          "not share a diagnosis. The table's own footnote "
                          "assigns rows 1-41 to Bozarth 2023 and rows 42-58 to "
                          "Gibellato 2024, both DEE series, but rows 59-64 to "
                          "Borck 2007, a CdLS series. Disease is therefore taken "
                          "PER ROW from `disease_harmonised_row` in the "
                          "extraction file, not from this paper-level default, "
                          "which applies only to rows with no per-row value "
                          "(D37, D61, D138). Attribution is accepted as an "
                          "OVERLAY on variants established elsewhere; the review "
                          "is not a source of variants or coordinates"),
        "classification_raw": "compiled from published case series",
        "classification": "LP",
        "classification_method": ("compiled review; no independent "
                                  "classification. Overlay only"),
        "table": "Table 2",
        "overlay_only": True,
    },
    "yuan2019_tableS2_curated.tsv": {
        "source": "Yuan_2019",
        "citation": "Yuan 2015",
        "pmid": "30158690",
        "design": "genetically-defined CdLS cohort (supplementary Table S2)",
        "disease_reported": "Cornelia de Lange syndrome",
        "disease_harmonised": "CdLS",
        "disease_basis": ("genetically-ascertained CdLS cohort; the authors note "
                          "the phenotypes skew towards the mild end of the CdLS "
                          "spectrum compared with phenotype-driven cohorts. "
                          "Table S2 gives transcript, hg19 coordinate, zygosity "
                          "and inheritance per case"),
        "classification_raw": "pathogenic or likely pathogenic (stated in the abstract)",
        "classification": "LP",
        "classification_method": ("stated as pathogenic or likely pathogenic in a "
                                  "diagnostic cohort; no ACMG criteria itemised, "
                                  "so mapped conservatively to LP"),
        "table": "Table S2",
    },
    "liu2009_table1_curated.tsv": {
        "source": "Liu_2009",
        "citation": "Liu 2009",
        "pmid": "19701948",
        "design": "primary cohort, 29 unrelated CdLS probands",
        "disease_reported": "Cornelia de Lange syndrome",
        "disease_harmonised": "CdLS",
        "disease_basis": ("Table 1 is the in-frame mutation spectrum across 29 "
                          "unrelated CdLS probands, with a per-variant proband "
                          "count from which affected family members are "
                          "explicitly excluded"),
        "classification_raw": "reported as CdLS-causing mutations",
        "classification": "LP",
        "classification_method": ("reported in a primary CdLS mutation-spectrum "
                                  "cohort; no ACMG criteria stated, so mapped "
                                  "conservatively to LP"),
        "table": "Table 1",
    },
    "musio2006_curated.tsv": {
        "source": "Musio_2006",
        "citation": "Musio 2006",
        "pmid": "16604071",
        "design": "primary report, the original SMC1A/SMC1L1 CdLS description",
        "disease_reported": "X-linked Cornelia de Lange syndrome",
        "disease_harmonised": "CdLS",
        "disease_basis": ("the original report establishing SMC1A (then SMC1L1) "
                          "as an X-linked CdLS gene: two families, NIPBL "
                          "excluded by sequencing, FISH and linkage analysis"),
        "classification_raw": "disease-causing mutation (original gene-disease report)",
        "classification": "P",
        "classification_method": ("the original gene-disease establishing report, "
                                  "with segregation in a multiplex family and a "
                                  "confirmed de novo occurrence in a second "
                                  "family; classified P rather than LP on that "
                                  "basis"),
        "table": "Figure 2 and text",
    },
    "borck2007_curated.tsv": {
        "source": "Borck_2007",
        "citation": "Borck 2007",
        "pmid": "",
        "design": "primary report, boys with unexplained CdLS",
        "disease_reported": "Cornelia de Lange syndrome",
        "disease_harmonised": "CdLS",
        "disease_basis": ("primary report of two de novo SMC1L1 (=SMC1A) missense variants in boys with unexplained CdLS; NIPBL variants in the same paper are stated on NM_133433.2 and excluded"),
        "classification_raw": "reported as CdLS-causing",
        "classification": "LP",
        "classification_method": ("reported in a primary CdLS cohort; no ACMG "
                                  "criteria stated, so mapped conservatively to LP"),
        "table": "text",
    },
    "pie2010_curated.tsv": {
        "source": "Pie_2010",
        "citation": "Pie 2010",
        "pmid": "",
        "design": "primary cohort, 30 CdLS patients",
        "disease_reported": "Cornelia de Lange syndrome",
        "disease_harmonised": "CdLS",
        "disease_basis": ("primary CdLS cohort reporting three SMC1A mutations among 30 patients (10% prevalence)"),
        "classification_raw": "reported as CdLS-causing",
        "classification": "LP",
        "classification_method": ("reported in a primary CdLS cohort; no ACMG "
                                  "criteria stated, so mapped conservatively to LP"),
        "table": "text",
    },
    "gervasini2013_curated.tsv": {
        "source": "Gervasini_2013",
        "citation": "Gervasini 2013",
        "pmid": "24124034",
        "design": "primary cohort, 8 CdLS patients with SMC1A mutations",
        "disease_reported": "Cornelia de Lange syndrome",
        "disease_harmonised": "CdLS",
        "disease_basis": ("primary report titled 'Cornelia de Lange Individuals "
                          "With New and Recurrent SMC1A Mutations Enhance "
                          "Delineation of...'; Tables I and II give growth, "
                          "facial features, limb reduction, milestones, "
                          "neurological features, seizures, intellectual "
                          "disability and a clinical severity score per patient"),
        "classification_raw": "reported as CdLS-causing mutations",
        "classification": "LP",
        "classification_method": ("reported in a primary CdLS cohort with "
                                  "per-patient phenotyping and in-silico "
                                  "prediction; no ACMG criteria stated, so "
                                  "mapped conservatively to LP"),
        "table": "Tables I and II",
    },
    "yuan2019_table1_curated.tsv": {
        "source": "Yuan_2019",
        "citation": "Yuan 2015",
        "pmid": "30158690",
        "design": "clinical exome sequencing cohort (Baylor Genetics), Table 1",
        "disease_reported": "cohesinopathy-gene variant on clinical exome",
        # Attribution is decided PER CASE from Yuan's own "CdLS as a
        # differential diagnosis?" column, not assumed from the journal's
        # CdLS framing. Handled in the row loop below.
        "disease_harmonised": "",
        "disease_basis": ("Yuan records per case whether CdLS was a differential "
                          "diagnosis. It was NOT in 13 of 14, so no CdLS "
                          "attribution is taken from this source for those cases "
                          "(D69). The paper is titled 'Clinical exome sequencing "
                          "reveals locus heterogeneity and phenotypic variability "
                          "of cohesinopathies'"),
        "classification_raw": "P or LP as formally classified by the authors",
        "classification": "LP",
        "classification_method": ("the authors give a formal P/LP classification "
                                  "per variant, which is used directly via "
                                  "classification_override rather than "
                                  "conservatively downgraded"),
        "table": "Table 1",
    },
    "kruszka2019_table1_curated.tsv": {
        "source": "Kruszka_2019",
        "citation": "Kruszka 2019",
        "pmid": "",
        "design": "primary cohort, holoprosencephaly with cohesin-complex variants",
        "disease_reported": "cohesin-associated holoprosencephaly",
        "disease_harmonised": "DEE85",
        "disease_basis": ("holoprosencephaly is WITHIN the DEE85 phenotypic spectrum, not a separate disease. OMIM 301044 is titled \"Developmental and epileptic encephalopathy 85 WITH OR WITHOUT MIDLINE BRAIN DEFECTS\", and holoprosencephaly is the canonical midline brain defect. OMIM additionally cites Kruszka as a primary source for DEE85 itself, alongside Goldstein 2015, Hansen 2013, Jansen 2016, Lebrun 2015 and Symonds 2017. Kruszka's own cohort reports intellectual disability in 13/13 and seizures in 3/12, and discusses SMC1A drug-resistant epilepsy citing Symonds 2017. Mapped to DEE85 (D108, reverting D105)."),
        "classification_raw": "reported as causative in an HPE cohort",
        "classification": "LP",
        "classification_method": ("reported with CADD scores and inheritance in a "
                                  "primary HPE cohort; no ACMG criteria stated, "
                                  "so mapped conservatively to LP"),
        "table": "Table 1",
    },
    "chuan2022_curated.tsv": {
        "source": "Chuan_2022",
        "citation": "Chuan 2022",
        "pmid": "",
        "design": "epilepsy gene-panel / exome diagnostic cohort",
        "disease_reported": "Cornelia de Lange syndrome-2",
        "disease_harmonised": "CdLS2",
        "disease_basis": ("supplementary Table S1 assigns a specific disease "
                          "(Cornelia de Lange-2) and inheritance mode (XLD) to this "
                          "individual, and Table S3 gives an explicit ACMG "
                          "classification for the variant"),
        "classification_raw": "Likely pathogenic (PS2, PM2, PM6, PP3)",
        "classification": "LP",
        "classification_method": ("authors' own ACMG classification with criteria "
                                  "stated explicitly; adopted as given"),
        "table": "Supplementary Tables S1 and S3",
    },
    "deardorff2007_curated.tsv": {
        "source": "Deardorff_2007",
        "citation": "Deardorff 2007",
        "pmid": "17273969",
        "design": "NIPBL-mutation-negative CdLS proband series",
        "disease_reported": "Cornelia de Lange syndrome",
        "disease_harmonised": "CdLS2",
        "disease_basis": ("clinically diagnosed CdLS probands ascertained as "
                          "NIPBL-mutation-negative and then found to carry SMC1A "
                          "variants; the paper characterises the SMC1A-positive "
                          "patients as a milder CdLS phenotype"),
        "classification_raw": "causative variants in clinically diagnosed CdLS probands",
        "classification": "LP",
        "classification_method": ("cohort-membership framing, applied consistently "
                                  "with the other cohort papers (METHODS 4A.16); no "
                                  "ACMG criteria stated, so mapped conservatively to LP"),
        "table": "Tables 2 and 3",
    },
    "dinardo2026_cdls_curated.tsv": {
        "source": "DiNardo_2026_CdLS",
        "citation": "DiNardo 2026",
        "pmid": "",
        "design": "patient-derived cell lines (CdLS arm of the transcriptome cohort)",
        "disease_reported": "Cornelia de Lange syndrome",
        "disease_harmonised": "CdLS2",
        "disease_basis": ("Supplementary Table 1 labels these samples CdL* and states "
                          "they are CdLS cell lines, distinguishing them from the EP "
                          "(epilepsy) samples in the same table"),
        "classification_raw": "pathogenic (patient cell lines in a characterised cohort)",
        "classification": "LP",
        "classification_method": ("cohort-membership framing, consistent with the other "
                                  "cohort papers (METHODS 4A.16); no ACMG criteria stated"),
        "table": "Supplementary Table 1",
    },
    "dinardo2026_curated.tsv": {
        "source": "DiNardo_2026",
        "citation": "DiNardo 2026",
        "pmid": "",
        "design": "primary functional/therapeutic study on patient LCLs",
        "disease_reported": "SMC1A-related developmental and epileptic encephalopathy",
        "disease_harmonised": "DEE85",
        "disease_basis": ("primary study of SMC1A-related DEE using "
                          "patient-derived lymphoblastoid cell lines, with "
                          "mutation-type-specific transcriptomic signatures and "
                          "readthrough therapy response"),
        "classification_raw": "pathogenic (patient-derived cell lines studied for therapy response)",
        "classification": "LP",
        "classification_method": ("patients recruited on an established SMC1A-DEE "
                                  "diagnosis with functional confirmation of "
                                  "reduced protein; no ACMG criteria stated, so "
                                  "mapped conservatively to LP"),
        "table": "figures and text",
    },
    "limongelli2009_curated.tsv": {
        "source": 'Limongelli_2009',
        "citation": 'Limongelli 2009',
        "pmid": '19449418',
        "design": 'single case report',
        "disease_reported": 'Cornelia de Lange syndrome with cardiac involvement',
        "disease_harmonised": 'CdLS2',
        "disease_basis": 'primary case report of a CdLS patient with a cardiac phenotype, variant identified by direct sequencing of SMC1A',
        "classification_raw": 'pathogenic (novel causative missense reported in a CdLS patient)',
        "classification": 'LP',
        "classification_method": 'reported as the novel causative variant in a clinically diagnosed CdLS patient; no ACMG criteria stated, mapped conservatively to LP',
        "table": 'text and Fig. 2A',
    },
    "hansen2013_curated.tsv": {
        "source": 'Hansen_2013',
        "citation": 'Hansen 2013',
        "pmid": '',
        "design": 'exome sequencing cohort',
        "disease_reported": 'Cornelia de Lange syndrome-2',
        "disease_harmonised": 'CdLS2',
        "disease_basis": 'variant table entry explicitly annotating the SMC1A change as de novo and assigning Cornelia de Lange syndrome-2 with X-linked dominant inheritance',
        "classification_raw": "pathogenic (de novo, novel, disease assigned in the authors' own table)",
        "classification": 'LP',
        "classification_method": "authors' table assigns a specific disease and de novo status but state no ACMG criteria; mapped conservatively to LP",
        "table": 'candidate variant table',
    },
    "cetinkaya2025_curated.tsv": {
        "source": "Cetinkaya_2025",
        "citation": "Cetinkaya 2025",
        "pmid": "41230206",
        "design": "17-patient CdLS case series with multi-gene panel and WES",
        "disease_reported": "Cornelia de Lange syndrome",
        "disease_harmonised": "CdLS2",
        "disease_basis": ("clinically diagnosed CdLS series; 2 of 17 patients carried "
                          "SMC1A variants, with CdLS clinical scores recorded for the "
                          "cohort"),
        "classification_raw": "pathogenic (one previously reported, one novel and pathogenic in silico)",
        "classification": "LP",
        "classification_method": ("authors assert pathogenicity without itemising ACMG "
                                  "criteria per variant; mapped conservatively to LP "
                                  "(METHODS 4A.16)"),
        "table": "Results text",
    },
    "amllal2025_curated.tsv": {
        "source": "Amllal_2025",
        "citation": "Amllal 2025",
        "pmid": "39513746",
        "design": "clinical vignette, two patients with truncating SMC1A variants",
        "disease_reported": ("intellectual disability with therapy-resistant epilepsy "
                             "(ID-TRE), the entity described by Jansen 2016"),
        "disease_harmonised": "DEE85",
        "disease_basis": ("both patients have intellectual disability with "
                          "therapy-resistant epilepsy, which the paper frames as the "
                          "ID-TRE entity of Jansen 2016 -- a source OMIM cites as a "
                          "DEE85 primary. NOTE patient 2 additionally has CdLS-type "
                          "facial dysmorphology and hand oligodactyly, so that "
                          "attribution is the less clean of the two and is flagged in "
                          "the variant note"),
        "classification_raw": "likely pathogenic (both); p.Gln222* pathogenic under ACMG on confirmed de novo status",
        "classification": "LP",
        "classification_method": ("authors state ACMG classification; both novel "
                                  "variants likely pathogenic, with c.664C>T reaching "
                                  "pathogenic on confirmed de novo status. Mapped to LP "
                                  "as the conservative common value"),
        "table": "text and Figure 1",
    },
    "yang2025_curated.tsv": {
        "source": "Yang_2025",
        "citation": "Yang 2025",
        "pmid": "39831465",
        "design": "single case report with a literature comparison",
        "disease_reported": "non-classic Cornelia de Lange syndrome with epilepsy",
        "disease_harmonised": "CdLS2",
        "disease_basis": ("the paper's own patient, diagnosed as non-classic CdLS with "
                          "epilepsy and reported as such in the title"),
        "classification_raw": "de novo frameshift reported as causative",
        "classification": "LP",
        "classification_method": ("de novo by trio WES and reported as causative; no "
                                  "ACMG criteria itemised, so mapped conservatively to LP"),
        "table": "Table 1 row 27 and Results text",
    },
    "yang2025_table1_curated.tsv": {
        "source": "Yang_2025_table1",
        "citation": "Yang 2025",
        "pmid": "39831465",
        "design": "literature compilation of published SMC1A epilepsy cases",
        "disease_reported": "classic and non-classic Cornelia de Lange syndrome with epilepsy",
        "disease_harmonised": "unknown",
        "disease_basis": ("a compilation, not a primary report, and its 'classic or "
                          "non-classic CdLS with epilepsy' framing overlaps "
                          "SMC1A-related DEE; taking it as a disease attribution would "
                          "inject CdLS labels for variants whose primaries call them "
                          "DEE85, so no disease is asserted"),
        "classification_raw": "not stated per row",
        "classification": "not_classified",
        "classification_method": "no classification asserted by the compilation",
        "table": "Table 1, rows 1-26",
        "overlay_only": True,
    },
    "chinen2019_curated.tsv": {
        "source": 'Chinen_2019',
        "citation": 'Chinen 2019',
        "pmid": '31123597',
        "design": 'single case report',
        "disease_reported": 'SMC1A-associated Cornelia de Lange syndrome',
        "disease_harmonised": 'CdLS2',
        "disease_basis": ("primary case report; the authors diagnose "
                          "SMC1A-associated Cornelia de Lange syndrome, recorded "
                          "here as the source states it. The variant is however "
                          "protein-truncating, a class otherwise associated with "
                          "DEE85, and other sources attribute the same variant to "
                          "DEE85 -- so `conflicting` is DERIVED by the consensus "
                          "rules from the disagreement between sources rather than "
                          "asserted by this entry. An earlier registration set "
                          "disease_harmonised directly to 'conflicting', which is a "
                          "consensus outcome and not a disease token, and therefore "
                          "never resolved (D112)."),
        "classification_raw": 'pathogenic (de novo nonsense in a clinically diagnosed patient)',
        "classification": 'LP',
        "classification_method": 'de novo confirmed by parental testing and a truncating consequence, but no ACMG criteria stated; mapped conservatively to LP',
        "table": 'text (Results) and Fig. 2',
    },
    "trajkova2024_curated.tsv": {
        "source": "Trajkova_2024",
        "citation": "Trajkova 2024",
        "pmid": "",
        "design": "DNA methylation episignature validation cohort",
        "disease_reported": "Cornelia de Lange syndrome-2",
        "disease_harmonised": "CdLS2",
        "disease_basis": ("the expected clinical condition (CDLS2) was tested "
                          "against a genome-wide DNA methylation profile and the "
                          "episignature returned CDLS, so the attribution is "
                          "supported by an orthogonal assay rather than by "
                          "phenotype alone"),
        "classification_raw": "LP (as stated in Table 1)",
        "classification": "LP",
        "classification_method": ("classification stated directly by the authors "
                                  "in Table 1; adopted as given"),
        "table": "Table 1, validation cohort row 32",
    },
    "tszchach2015_curated.tsv": {
        "source": 'Tzschach_2015',
        "citation": 'Tzschach 2015',
        "pmid": '26235285',
        "design": 'X-linked intellectual disability gene-panel cohort',
        "disease_reported": 'X-linked intellectual disability (no syndrome-level diagnosis)',
        "disease_harmonised": 'unknown',
        "disease_basis": 'XLID panel cohort; the authors state only that clinical features are compatible with an SMC1A phenotype and make no CdLS or DEE diagnosis, so no disease can be attributed',
        "classification_raw": "pathogenic per authors' panel reporting, supported by in silico predictions only",
        "classification": 'VUS',
        "classification_method": 'no functional or segregation evidence and no syndrome-level diagnosis; in silico predictions alone do not support a confident pathogenic call, so mapped to VUS',
        "table": 'Table 1, family 17',
    },
    "revenkova2009_curated.tsv": {
        "source": "Revenkova_2009",
        "citation": "Revenkova 2009",
        "pmid": "19218272",
        "design": "functional study of CdLS-associated SMC1A missense variants",
        "disease_reported": "Cornelia de Lange syndrome",
        "disease_harmonised": "CdLS2",
        "disease_basis": ("functional study of variants already established as "
                          "CdLS-associated; contributes functional evidence rather "
                          "than new patients, so it asserts no new disease attribution"),
        "classification_raw": "CdLS-associated missense variants selected for functional study",
        "classification": "not_classified",
        "classification_method": ("no per-variant classification is given; the source "
                                  "contributes functional evidence only and the "
                                  "classification comes from the primary reports"),
        "table": "text (legacy protein-only notation)",
    },
    "fateh2024_curated.tsv": {
        "source": "Fateh_2024",
        "citation": "Fateh 2024",
        "pmid": "",
        "design": "single case report (whole exome sequencing)",
        "disease_reported": "Cornelia de Lange syndrome (CdLSp2)",
        "disease_harmonised": "CdLS2",
        "disease_basis": ("clinically diagnosed non-classic CdLS, de novo with both "
                          "parents confirmed normal, and the SMC1A variant reported as "
                          "confirming the diagnosis. CORRECTED 2026-09-25 (D149): this "
                          "entry previously carried a caveat that the same patient also "
                          "carried a NIPBL variant. It does not. Fateh reports TWO "
                          "cases -- case 1 is classic CdLS whose blood exome was "
                          "negative and whose mosaic NIPBL variant "
                          "(NM_133433.4:c.6534_6535del) was found only in skin-derived "
                          "DNA, and case 2 is this non-classic CdLS patient with the "
                          "SMC1A variant and no NIPBL finding. The caveat came from "
                          "misreading 'Similar to the NIPBL variant, the SMC1A variant "
                          "was classified as likely pathogenic', which compares the "
                          "classification of the two separate cases. The attribution is "
                          "therefore cleaner than previously recorded, not weaker"),
        "classification_raw": "likely pathogenic (authors' own classification)",
        "classification": "LP",
        "classification_method": ("authors state the variant was classified as likely "
                                  "pathogenic; no ACMG criteria itemised"),
        "table": "text and Fig. 2",
    },
    "gibellato2024_curated.tsv": {
        "source": "Gibellato_2024",
        "citation": "Gibellato 2024",
        "pmid": "",
        "design": "large international SMC1A epilepsy syndrome (DEE85) cohort",
        "disease_reported": ("SMC1A epilepsy syndrome / developmental and epileptic "
                             "encephalopathy-85 with or without midline brain defects"),
        "disease_harmonised": "DEE85",
        "disease_basis": ("dedicated DEE85 cohort; the paper defines its population as "
                          "SMC1A epilepsy syndrome (DEE85, OMIM 301044) and states the "
                          "syndrome is described only in female patients"),
        "classification_raw": "pathogenic (cohort inclusion as characterised SMC1A-DEE variants)",
        "classification": "LP",
        "classification_method": ("cohort-membership framing; the paper gives a variant "
                                  "class per patient but no per-variant ACMG tier, so "
                                  "mapped conservatively to LP"),
        "table": "Table 4",
    },
    "huisman2017_table5_curated.tsv": {
        "source": "Huisman_2017",
        "citation": "Huisman 2017",
        "pmid": "",
        "design": ("cohort of 51 individuals plus a literature-wide "
                   "concordance table (Table 5)"),
        "disease_reported": "SMC1A-related disorder (CdLS and/or epilepsy)",
        # Huisman spans BOTH phenotypes and does not assign a single condition
        # per variant in Table 5, so no disease is attributed from this source.
        # Its value is provenance: primary patient identifiers for the older
        # CdLS literature, which the cross-source individual keys require.
        "disease_harmonised": "",
        "disease_basis": ("Table 5 is a genotype concordance table spanning "
                          "both CdLS and epilepsy presentations and does not "
                          "assign a condition per variant; no disease "
                          "attribution is taken from it. Used for provenance "
                          "and for cross-referencing primary patient ids."),
        "classification_raw": "reported in an SMC1A cohort; no classification stated",
        "classification": "LP",
        "classification_method": (
            "cohort-membership framing, applied consistently with the other cohort "
            "papers (Bozarth 2023, Symonds 2017, Gibellato 2024). Huisman's stated "
            "inclusion criterion is explicit: both cohorts were assembled by asking "
            "centres to identify individuals with *pathological* SMC1A variants "
            "(sections 2.2 and 2.3), so cohort membership is itself a pathogenicity "
            "assertion by the submitting molecular laboratories. No ACMG criteria are "
            "stated, so mapped conservatively to LP rather than P. NOTE this licenses "
            "the CLASSIFICATION only and deliberately does NOT license a disease "
            "attribution: Huisman states its cohort includes individuals in whom CdLS "
            "was not clinically suspected and whose main manifestation is an epileptic "
            "encephalopathy, so it is a mixed CdLS/DEE series. Disease is left to the "
            "normal consensus rules (D92)."),
        "table": "Table 5",
    },
    "symonds2017_table2_curated.tsv": {
        "source": "Symonds_2017",
        "citation": "Symonds 2017",
        "pmid": "",
        "design": "primary case series (10 new cases), DEE85-defining report",
        "disease_reported": "severe early onset epilepsy with cluster seizures",
        "disease_harmonised": "DEE85",
        "disease_basis": ("primary series titled 'Heterozygous truncation "
                          "mutations of SMC1A cause a severe early onset "
                          "epilepsy with cluster seizures in females: Detailed "
                          "phenotyping of 10 new cases'; all cases de novo and "
                          "heterozygous, with per-case phenotyping and "
                          "X-inactivation studies"),
        "classification_raw": "pathogenic (author assertion; no ACMG criteria stated)",
        "classification": "LP",
        "classification_method": ("author assertion in the primary DEE85-defining "
                                  "series; no ACMG criteria stated, so mapped "
                                  "conservatively to LP rather than P"),
        "table": "Table 2",
    },
    "baranano2022_table1_curated.tsv": {
        "source": "Baranano_2022",
        "citation": "Baranano 2022",
        "pmid": "",
        "design": "primary case series (13 patients)",
        "disease_reported": "SMC1A loss-of-function epilepsy",
        "disease_harmonised": "DEE85",
        "disease_basis": ("primary series titled 'Further Characterization of "
                          "Loss of Function Epilepsy Distinct From Cornelia de "
                          "Lange Syndrome'; eligibility required a documented "
                          "SMC1A variant AND epilepsy, with per-patient de novo "
                          "status and X-inactivation studies"),
        "classification_raw": "pathogenic (author assertion; no ACMG criteria stated)",
        "classification": "LP",
        "classification_method": ("author assertion in a peer-reviewed primary "
                                  "epilepsy series; no ACMG criteria stated, so "
                                  "mapped conservatively to LP rather than P"),
        "table": "Table 1",
    },
}

# Citation strings seen in extraction files -> compact form used in keys, and
# the PMID where independently known (from LOVD's citation markup).
CITATION_PMID = {
    "Deardorff 2007": "17273969",
    "Musio 2006": "16604071",
    "Liu 2009": "19701948",
    "Borck 2007": "17221863",
    "Goldstein 2015": "",
    "Symonds 2017": "",
    "Baranano 2022": "",
    "Huisman 2017": "",
    "Gorman 2017": "",
    "Kruszka 2019": "",
    "Chinen 2019": "",
    "Fang 2021": "",
    "Jansen 2016": "",
    "Lebrun 2015": "",
    "Elwan 2022": "",
    "Oguni 2019": "",
    "Wenger 2016": "",
    "Bozarth 2023": "",
}


def compact(citation: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", citation)


REF_SPLIT = re.compile(r"\s*;\s*")
REF_PARSE = re.compile(r"^(.*?)\s*(?:\((P\d+[^)]*)\))?$")


def parse_refs(field: str):
    """'Baranano 2022(P6); Symonds 2017' -> [('Baranano 2022','P6'), ...]"""
    out = []
    for part in REF_SPLIT.split(L.clean(field)):
        if not part:
            continue
        m = REF_PARSE.match(part.strip())
        cite = L.clean(m.group(1))
        pt = L.clean(m.group(2))
        if cite:
            out.append((cite, pt))
    # collapse duplicates that arise from a citation appearing with and
    # without a patient id in the same cell
    seen, dedup = set(), []
    for cite, pt in out:
        if pt:
            seen.add(cite)
            dedup.append((cite, pt))
    for cite, pt in out:
        if not pt and cite not in seen:
            dedup.append((cite, pt))
    return dedup


# XCI wording is free text; normalise to a comparable value while keeping the
# original, because skewing is mechanistically important for X-linked female
# DEE85 but is reported inconsistently.
def normalise_xci(s: str):
    t = L.clean(s).lower().replace(" ", "")
    if not t or t in ("n/a", "na"):
        return "", ""
    ratio = re.search(r"(\d{1,3}):(\d{1,3})", t)
    r = f"{ratio.group(1)}:{ratio.group(2)}" if ratio else ""
    if "random" in t:
        return "random", r
    if "highly" in t:
        return "highly_skewed", r
    if "moderately" in t:
        return "moderately_skewed", r
    if "skewed" in t:
        return "skewed", r
    return "other", r


def main():
    files = sorted(glob.glob(os.path.join(RAW, "*_curated.tsv")))
    if not files:
        sys.exit(f"FATAL: no *_curated.tsv extraction files in {RAW}")

    say("# Literature ingest report")
    say()
    say("Extraction (transcribing a paper) is kept separate from ingestion "
        "(reshaping onto the ledger schema). This script only does the latter; "
        "every hand correction made during extraction is recorded in "
        "`docs/DECISIONS_LOG.md`.")
    say()

    obs = []
    for f in files:
        base = os.path.basename(f)
        meta = PAPERS.get(base)
        if meta is None:
            say(f"**SKIPPED** `{base}`: no entry in the PAPERS registry. Add one "
                f"so its provenance and disease basis are explicit.")
            continue
        rows = list(csv.DictReader(open(f, encoding="utf-8"), delimiter="\t"))
        say(f"## `{base}`")
        say()
        say(f"- source: **{meta['citation']}** ({meta['design']}), {meta['table']}")
        say(f"- rows: {len(rows)}")
        per_row = [L.clean(r.get("disease_harmonised_row")) for r in rows]
        if any(per_row):
            from collections import Counter as _C
            brk = _C(x or f"(paper default: {meta['disease_harmonised']})"
                     for x in per_row)
            say(f"- disease attributed **per row**: "
                + ", ".join(f"`{k}` x{v}" for k, v in sorted(brk.items()))
                + f" -- {meta['disease_basis']}")
        else:
            say(f"- disease attributed: `{meta['disease_harmonised']}` "
                f"-- {meta['disease_basis']}")
        say(f"- classification: `{meta['classification']}` "
            f"from \"{meta['classification_raw']}\"")
        say()

        n_noref = 0
        for r in rows:
            variant = L.clean(r.get("variant")).replace("C.", "c.")
            refs = parse_refs(r.get("refs", ""))
            if not refs:
                refs = [(meta["citation"], "")]
                n_noref += 1
            xci_class, xci_ratio = normalise_xci(r.get("xci"))

            notes = []
            if L.clean(r.get("de_novo_reported")):
                notes.append(f"de novo (per {meta['citation']}): "
                             f"{L.clean(r['de_novo_reported'])}")
            if L.clean(r.get("predicted_effect_reported")):
                notes.append(f"predicted effect (per {meta['citation']}): "
                             f"{L.clean(r['predicted_effect_reported'])}")
            if L.clean(r.get("clinical_detail")):
                notes.append(L.clean(r["clinical_detail"]))
            if L.clean(r.get("domain")):
                notes.append(f"protein domain (per {meta['citation']}): "
                             f"{L.clean(r['domain'])}")
            if xci_class:
                notes.append(f"XCI: {xci_class}"
                             + (f" ({xci_ratio})" if xci_ratio else ""))
            if L.clean(r.get("note")):
                notes.append(f"extraction note: {L.clean(r['note'])}")

            # optional explicit family map, e.g. "7P=7; 7S=7; 8P=8; 8S=8"
            fam_of = {}
            for part in L.clean(r.get("family_by_patient")).split(";"):
                if "=" in part:
                    k, val = part.split("=", 1)
                    if k.strip() and val.strip():
                        fam_of[k.strip()] = val.strip()

            def fam_hit(cite_s, pt):
                """Family for this patient, or None.

                The patient identifier may arrive either as the parsed patient
                field or folded into the citation string, depending on how the
                source wrote it, so both are checked. Longest identifier first
                so `10P` is not matched by `0P`.
                """
                if pt in fam_of:
                    return fam_of[pt]
                for k in sorted(fam_of, key=len, reverse=True):
                    if k in str(cite_s):
                        return fam_of[k]
                return None

            for cite, patient in refs:
                is_rereport = compact(cite) != compact(meta["citation"])
                key_pt = patient or "unlabelled"
                ind_key = f"PUB:{compact(cite)}:{key_pt}"
                obs.append({
                    "source": meta["source"],
                    "source_detail": f"{meta['citation']} {meta['table']}",
                    "source_record_id": f"{meta['source']}:row{r.get('row')}",
                    "individual_key": ind_key,
                    # Family grouping is taken ONLY from an explicit `family_by
                    # _patient` mapping in the extraction file, never inferred
                    # from the shape of a patient identifier. Sibling suffixes
                    # are a convention, not a guarantee, and a regex that merged
                    # unrelated patients would UNDERCOUNT probands silently --
                    # the mirror of over-attributing from siblings. Where a
                    # source states the relationship, it is recorded; otherwise
                    # each patient is its own family.
                    "family_key": (
                        f"PUB:{compact(meta['citation'])}:FAM:{fam_hit(cite, key_pt)}"
                        if fam_hit(cite, key_pt) else ind_key),
                    "family_basis": (
                        f"source states this patient belongs to family "
                        f"{fam_hit(cite, key_pt)}" if fam_hit(cite, key_pt)
                        else "no relationship stated; treated as unrelated"),
                    "is_independent_proband": "yes" if patient else "assumed",
                    "proband_independence_basis": (
                        "patient individually numbered in the primary report"
                        if patient else
                        "no patient identifier given in the citing table; "
                        "independence assumed, may be a duplicate of another row"),
                    "kinship_needs_adjudication": str(not patient),

                    "hgvs_c_input": variant,
                    "hgvs_c_input_transcript": L.TRANSCRIPT,
                    "hgvs_p_input": L.clean(r.get("protein")),
                    "variant_class_supplied": L.clean(r.get("vtype")),

                    "pathogenicity_reported": (L.clean(r.get("classification_override"))
                                               or meta["classification"]),
                    "pathogenicity_reported_raw": meta["classification_raw"],
                    "classification_method_raw": meta["classification_method"],

                    # A per-row disease in the extraction file overrides the
                    # paper-level default. Needed where one table compiles
                    # patients from several primary series that do not share a
                    # diagnosis: Astorino 2025 Table 2 states in its own
                    # footnote that rows 1-41 come from Bozarth 2023 and rows
                    # 42-58 from Gibellato 2024, both DEE series, but rows
                    # 59-64 from Borck 2007, a CdLS series. Attributing the
                    # review's own headline disease to every row labelled four
                    # Borck CdLS patients as DEE85, and for c.587G>A and
                    # c.3254A>G that put Borck 2007 in contradiction with
                    # itself -- read directly it asserts CdLS, read through the
                    # review it asserted DEE85. Three variants were excluded
                    # from the CdLS arm on that manufactured conflict.
                    "disease_reported_raw": (L.clean(r.get("disease_reported_row"))
                                             or meta["disease_reported"]),
                    # Per-case CdLS qualifier overrides both. A source that
                    # records "CdLS was NOT the differential diagnosis" has made
                    # a determination, and what it leaves behind is a disease
                    # that cannot be mapped to either arm -- not the absence of
                    # any disease record. Emitting `unknown` rather than an
                    # empty token keeps that distinction: the variant lands in
                    # `Exclude_disease_uncertain` (a disease was recorded and is
                    # unresolvable) rather than `Exclude_pending_classification`,
                    # which promises that more patient detail could still
                    # upgrade it. For a source that has already considered the
                    # diagnosis and declined it, nothing is pending. Upgradeable
                    # is the right reading for the cohort records that carry HPO
                    # terms but no syndrome label; it is the wrong reading here
                    # (D140).
                    "disease_harmonised": ("CdLS"
                        if L.clean(r.get("cdls_differential_reported")) == "yes"
                        else "unknown"
                        if L.clean(r.get("cdls_differential_reported")) == "no"
                        else (L.clean(r.get("disease_harmonised_row"))
                              or meta["disease_harmonised"])),
                    "cdls_differential_reported": L.clean(r.get("cdls_differential_reported")),
                    "dual_molecular_diagnosis": L.clean(r.get("dual_molecular_diagnosis")),
                    "zygosity_reported": L.clean(r.get("zygosity_reported")),
                    "attribution_basis": meta["disease_basis"],

                    "xci_skewing": xci_class,
                    "xci_ratio": xci_ratio,
                    "xci_raw": L.clean(r.get("xci")),
                    "protein_domain_reported": L.clean(r.get("domain")),

                    "is_rereport": str(is_rereport),
                    "rereport_of_citation": cite if is_rereport else "",
                    "rereport_of_patient": patient if is_rereport else "",
                    "rereport_of_pmid": CITATION_PMID.get(cite, ""),
                    "primary_citation": cite,
                    "primary_patient_id": patient,
                    "pmids": CITATION_PMID.get(cite, ""),
                    "citations": cite,

                    "functional_evidence": L.clean(r.get("functional_evidence")),
                    # Sex is a property of the PATIENT, so a per-variant row
                    # carries the union across that variant's patients, and
                    # `sex_by_patient` preserves the per-patient assignment so
                    # the union can always be decomposed.
                    "sex": L.clean(r.get("sex")),
                    "sex_by_patient": L.clean(r.get("sex_by_patient")),
                    "hgvs_c_derived_from_protein": L.clean(
                        r.get("hgvs_c_derived_from_protein")),
                    "episignature_result": L.clean(r.get("episignature_result")),
                    "episignature_source": L.clean(r.get("episignature_source")),
                    "is_overlay_source": str(bool(meta.get("overlay_only"))),
                    "competing_attribution": L.clean(r.get("competing_attribution")),
                    "additional_notes": "; ".join(notes),
                    "phenotype_detail": "; ".join(
                        x for x in [L.clean(r.get("clinical_detail")),
                                    (f"X-inactivation {xci_class}"
                                     + (f" {xci_ratio}" if xci_ratio else ""))
                                    if xci_class else "",
                                    (f"de novo: {L.clean(r.get('de_novo_reported'))}"
                                     if L.clean(r.get("de_novo_reported")) else "")]
                        if x),
                    "inheritance": ("de_novo" if L.clean(r.get("de_novo_reported","")).upper().startswith("Y")
                                    else ""),
                    "eligible_for_curated_set": "True",
                })
        say(f"- rows with no primary attribution parsed (attributed to the "
            f"citing paper itself): {n_noref}")
        say()

    # ---- summary ----
    say("## Observations built")
    say()
    say(f"- observation rows: {len(obs)}")
    say(f"- distinct variants: {len({o['hgvs_c_input'] for o in obs})}")
    say(f"- distinct individuals (keyed on the PRIMARY report): "
        f"{len({o['individual_key'] for o in obs})}")
    say(f"- rows that are re-reports: "
        f"{len([o for o in obs if o['is_rereport'] == 'True'])}")
    say()

    say("### Individuals per primary report")
    say()
    say("These are the papers to read directly next; reading them will attach "
        "to the SAME individual keys, so the patients collapse rather than "
        "duplicating.")
    say()
    say("| primary report | individuals | variants |")
    say("|---|---|---|")
    per = defaultdict(lambda: [set(), set()])
    for o in obs:
        per[o["primary_citation"]][0].add(o["individual_key"])
        per[o["primary_citation"]][1].add(o["hgvs_c_input"])
    for cite, (inds, vars_) in sorted(per.items(), key=lambda kv: -len(kv[1][0])):
        say(f"| {cite} | {len(inds)} | {len(vars_)} |")
    say()

    say("### XCI skewing captured")
    say()
    say("| XCI | observations |")
    say("|---|---|")
    for k, n in Counter(o["xci_skewing"] or "(not reported)" for o in obs).most_common():
        say(f"| {k} | {n} |")
    say()
    say(f"Observations with a numeric skewing ratio: "
        f"{len([o for o in obs if o['xci_ratio']])}")
    say()

    say("### Variants with no patient identifier")
    say()
    unl = sorted({o["hgvs_c_input"] for o in obs
                  if not o["primary_patient_id"]})
    say(f"{len(unl)} variants are cited without a patient number, so their "
        f"individual identity cannot be established and they are flagged "
        f"`kinship_needs_adjudication`: "
        + ", ".join(f"`{v}`" for v in unl))
    say()

    cols = sorted({k for o in obs for k in o})
    out = os.path.join(INTERIM, "literature_observations.tsv")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for o in obs:
            fh.write("\t".join(str(o.get(c, "")).replace("\t", " ").replace("\n", " ")
                               for c in cols) + "\n")
    say(f"Wrote `data/interim/literature_observations.tsv` "
        f"({len(obs)} rows, {len(cols)} columns).")

    os.makedirs(DOCS, exist_ok=True)
    with open(os.path.join(DOCS, "literature_ingest.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(report) + "\n")
    print("Wrote docs/verification/literature_ingest.md")


if __name__ == "__main__":
    main()
