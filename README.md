# SMC1A SGE

## HDR Oligo Design

### Background:
To design HDR libraries for SMC1A, all coding exons present in the MANE transcript were first extracted from a whole genome GTF file. Exon coordinates and variant design parameters were then entered into a manifest and submitted to the Wellcome Sanger Gene Editing informatics (MAVE-SGE) design pipeline. Once the manifest had been processed, [VaLiAnT](https://github.com/cancerit/VaLiAnT) (Variant Library Annotation Tool) output files containing all HDR oligo sequences underwent QC and approval. This was done to ensure appropriate sgRNA selection, PAM protection variant selection and the presence of the correct variant profile. Following design approval, HDR libraries were ordered by the gene editing team from TWIST Bioscience and preparatory cloning was initiated.

*Oligo design (manifest, VaLiAnT output, QC/approval records) was carried out via the
[Step1_valiant](https://gitlab.internal.sanger.ac.uk/team302/sge/sge_production_screen_qc/-/tree/feature/run46503_Jan2023/Step1_valiant?ref_type=heads)
stage of Sanger's internal SGE production screen QC pipeline (internal link — not
accessible outside the Sanger network). Wet-lab cloning and library assembly were
executed via Sanger's internal SGE pipeline by the Wellcome Sanger gene editing team.
Neither is part of this repo.*

---

## Data Analysis

### Background:
Downstream analysis of SGE screening data is performed using the established Sanger **sge-metapipeline**, which wraps sample retrieval (iRODS to lustre) and quantification (QUANTS) into a single pipeline.

* QUANTS output docs: [QUANTS output.md](https://github.com/cancerit/QUANTS/blob/develop/docs/output.md)
* `irods_to_lustre` is run via the internal `sge-metapipeline` module — no public repo currently exists for this component.


---

### STEP ONE: Load the sge-metapipeline

#### Requirements:
* Access to the `HGI/common/sge-metapipeline` module

#### Running the script:
Run:
```bash
module load HGI/common/sge-metapipeline/v0.1.1
```
Once loaded, there are two basic commands, `fetch` and `execute`, each with a `--help` option:
```bash
sge-metapipeline fetch --help
sge-metapipeline execute --help
```

#### Notes:
* `fetch` retrieves sample information from the MLWH and creates the initial Sample Manifest
* `execute` runs the `irods_to_lustre` and `QUANTS` sub-commands

---

### STEP TWO: Fetch sample manifest and run iRODS to lustre

#### Requirements:
* Study ID and run ID for the samples of interest
* `config.ini` file with warehouse credentials

#### Running the script:
First, load the iRODS module and initiate a session:
```bash
module load ISG/IRODS/1.0
iinit
```
*`iinit` will prompt you to enter your Sanger password.*

Fetch sample information from the MLWH (example):
```bash
sge-metapipeline fetch --study-id <study_id> --run-id <run_id> --warehouse-creds config.ini
```
Modify the resulting `manifest_fetched.tsv` if needed, to include only the desired samples.

Then run iRODS to lustre (example):
```bash
module load badger/samtools/1.20
sge-metapipeline execute irods_to_lustre --manifest-file manifest_fetched.tsv --output irods_to_lustre --pipeline-config config.ini
```

#### Output:
`manifest_fetched.irods_to_lustre.tsv`, along with merged cram/fastq files on lustre, to be used as input for QUANTS.

#### Notes:
* You may need to load `samtools` (`badger/samtools/1.20`) before running `irods_to_lustre`

---

### STEP THREE: Run QUANTS

#### Requirements:
* `meta.csv` file (oligo design file, from VaLiAnT output)
* `manifest_fetched.irods_to_lustre.tsv`, populated with the path to the `meta.csv` file
* pipeline `config.ini` file

#### Running the script:
Run (example):
```bash
sge-metapipeline execute QUANTS --manifest-file irods_to_lustre/manifest_fetched.irods_to_lustre.tsv --output SMC1A_QUANTS --pipeline-config sge_metapipeline.config.ini
```

#### Output:
See expected output structure: [QUANTS output.md](https://github.com/cancerit/QUANTS/blob/develop/docs/output.md)

---

### STEP FOUR: Run MAVEQC

#### Background:
QC is performed using an adapted version of [MAVE-QC](https://github.com/wtsi-hgi/MAVEQC), with the adapted script [`run_maveqc_VB.R`](Code/pipeline/run_maveqc_VB.R).

**Lineage**: this script is adapted from `run_maveqc_VB.R` in the
[`vb9Sanger/5-UTR-SGE`](https://github.com/vb9Sanger/5-UTR-SGE) pilot repo. One
SMC1A-specific change from that original: the 5'UTR project targeted non-coding
sequence, so its variants were annotated generically as `SNV`, and that class was
included in the background-fitting pool (`fit_consequences`) used to build the
PTV/no-impact anchors (see the terminology note at the top of STEP FIVE). SMC1A is coding, so variants carry `Missense_Variant`/
`Synonymous_Variant`/etc. labels instead — `SNV` never applies here and has been
removed from that list. Everything else differing from the 5'UTR original is a bug fix,
documented in this script's own header comment (report generation was previously never
triggered, two QC output TSVs were silently never written, an off-by-one loop crashed on
empty file globs, etc.) rather than a change to the analysis itself.

#### Requirements:
* `HGI/common/sge-metapipeline` module
* `manifest_fetched.irods_to_lustre.quants.tsv` manifest generated by QUANTS
* pipeline `config.ini` file
* [`run_maveqc_VB.R`](Code/pipeline/run_maveqc_VB.R)

#### Running the script:
First, load the sge-metapipeline module (if not already loaded):
```bash
module load HGI/common/sge-metapipeline/v0.1.1
```
Then run MAVEQC (example):
```bash
sge-metapipeline execute MAVEQC --manifest-file SMC1A_QUANTS/manifest_fetched.irods_to_lustre.quants.tsv --output SMC1A_MAVEQC --pipeline-config sge_metapipeline.config.ini
```
This will generate a folder called **'input'**. Once this has been generated, the job can be killed. Delete the *contents* of the output folder (but not the output folder itself).

Then run the adapted script:
```bash
Rscript Code/pipeline/run_maveqc_VB.R <path/to/input_dir> <path/to/output_dir>
```
*where `<path/to/input_dir>` is the folder one level above the 'input' folder generated above, and `<path/to/output_dir>` is the desired output location.*

Run once per targeton, once against a Day4 reference and once against a Plasmid
reference — see STEP SIX for why both are needed.

#### Output:
`plasmid_qc/`, `screen_qc/`, `experiment_qc/` subfolders, DESeq2 results, normalized counts, `MAVEQC_report.html`.

---

### STEP FIVE: Classify variants: shrinkage-based Gaussian anchor model

#### Background:

**Terminology note**: this step's anchors are built from an annotation-based consequence
class — nonsense/frameshift calls, i.e. **protein-truncating variants (PTVs)** — not from
a confirmed functional outcome. The scripts, their CLI flags, and their output columns
all name this class `LOF` (e.g. `anchor_mu_lof`, `--lof_label`), for historical reasons;
this README uses **PTV** in prose instead, reserving "loss-of-function" for the actual
functional/mechanistic claim (which is what the assay's `anchor_tier` result speaks to,
not what goes into building the anchor). The distinction matters concretely: not every
PTV causes loss of function (NMD escape being the clearest exception, handled explicitly
below), so treating the annotation label and the functional outcome as interchangeable
would beg the question this step's classification is meant to test.

Run on the Day15-vs-reference DESeq2 output files from MAVEQC (e.g. `APDY_exon2_all_deseq2_results_condition_Day15_vs_Day4.tsv`), [`gaussian_shrinkage_classifier.R`](Code/pipeline/gaussian_shrinkage_classifier.R) classifies every variant into `enriched` / `no impact` / `weakly depleting` / `strongly depleting`.

It sits between two extremes:
* a **gene-wide anchor** (one shared PTV/no-impact distribution pooled across all targetons) — maximum statistical power, but can't adapt to a targeton that's genuinely different
* a **per-targeton anchor** (each targeton fit in total isolation) — adapts fully to each exon, but small-n exons can be yanked around by a handful of noisy points

For each targeton, this script computes both the gene-wide anchor and that targeton's own local anchor, then **shrinks** the local anchor toward the gene-wide one, weighted by how much data that targeton has (a pseudo-count rule, the same logic as DESeq2's own dispersion shrinkage):

```
mu_shrunk_i  = (n_i * mu_local_i  + k * mu_global)  / (n_i + k)
var_shrunk_i = (n_i * var_local_i + k * var_global) / (n_i + k)
```

applied separately to the PTV anchor and the no-impact anchor. A well-powered targeton (n_i >> k) stays close to its own local estimate; a thin targeton (n_i << k) is pulled toward the gene-wide value. `k` defaults to `"auto"` (median per-targeton n across the gene, separately for PTV and no-impact). Its sensitivity to that choice is checked in `Code/investigations/compare_shrinkage_k_runs.py` (see STEP NINE).

Every row in every targeton is then classified using *that targeton's own shrunk anchors* and shrunk 95th-percentile PTV threshold (posterior-probability depleted/no-impact call, with weak/strong tiering).

An `exon_map` (`SMC1A_exon_map.tsv`) was supplied for this run, enabling NMD-escape-aware PTV filtering — PTVs falling in the final exon, or within `--nmd_escape_distance` (default 50bp) of the final exon-exon junction in the penultimate exon, are treated as plausible NMD-escape variants and excluded from the PTV fitting pool (both gene-wide and local), since they may not behave like true loss-of-function alleles.

#### Requirements:
* R package `data.table` (required); `ggplot2` + `ragg` (only if `--plot_dir` is set)
* Day15-vs-reference DESeq2 output TSVs from MAVEQC (e.g. `*_all_deseq2_results_condition_Day15_vs_Day4.tsv`)
* `SMC1A_exon_map.tsv` — columns `Exon_position` (chrN:start-end), `EXON` (integer, transcript order), `Targeton_ID`
* [`gaussian_shrinkage_classifier.R`](Code/pipeline/gaussian_shrinkage_classifier.R)

#### Running the script:
Run (example):
```bash
Rscript Code/pipeline/gaussian_shrinkage_classifier.R \
  --input "deseq2_results/*_all_deseq2_results_condition_Day15_vs_Day4.tsv" \
  --out_dir GMM_shrinkage_results \
  --plot_dir GMM_shrinkage_results/plots \
  --exon_map SMC1A_exon_map.tsv
```
Full flag list: run with `--help`. Run once per reference condition (Day4, Plasmid) — see STEP SIX.

#### Output:
* Per-targeton output TSVs (same rows/columns as input, plus `anchor_mu_lof`, `anchor_sd_lof`, `anchor_mu_lof_local`, `anchor_mu_lof_global`, `anchor_weight_lof_local`, `anchor_mu_noimpact`, `anchor_sd_noimpact`, `anchor_mu_noimpact_local`, `anchor_mu_noimpact_global`, `anchor_weight_noimpact_local`, `anchor_direction`, `anchor_lof_threshold`, `anchor_post_lof`, `anchor_call`, and final `anchor_tier`)
* `gmm_shrinkage_anchor_summary.tsv` — one row per targeton, with local/global/shrunk fit parameters and tier counts
* If `--plot_dir` is set: 2 diagnostic PNGs per targeton (`<targeton>_anchor_fit.png` — all variants; `<targeton>_missense_anchor_fit.png` — missense only), showing the shrunk fit (solid) vs gene-wide fit (dotted) and the shrunk 95th-percentile PTV threshold

#### Notes:
* Rows with a pre-existing `"enriched"` status (from `stat_pos_raw` by default) are passed through as `enriched` and excluded from the anchor-fitting pool.

---

### STEP SIX: Choose a reference condition per targeton

#### Background:
[`qc_signal_noise_metrics_with_composite.R`](Code/pipeline/qc_signal_noise_metrics_with_composite.R) computes per-targeton signal-to-noise QC metrics (median |Z|, control neutrality, positional-bias LOESS residual diagnostic, Cohen's d for PTV vs controls) for the Day4- and Plasmid-referenced runs of STEPS FOUR/FIVE, and ranks the two references per targeton via a hierarchical/cascading composite score, to decide which one's output should actually be used downstream.

#### Requirements:
* R package `optparse`
* Day15-vs-reference DESeq2/anchor-tier TSVs for **both** reference conditions
* [`qc_signal_noise_metrics_with_composite.R`](Code/pipeline/qc_signal_noise_metrics_with_composite.R)

#### Running the script:
Run (example):
```bash
Rscript Code/pipeline/qc_signal_noise_metrics_with_composite.R \
  --input "GMM_shrinkage_results_D4/*_all_deseq2_results_condition_*.tsv,GMM_shrinkage_results_Plasmid/*_all_deseq2_results_condition_*.tsv" \
  --out_prefix reference_qc
```

#### Output:
* `<out_prefix>.metrics.tsv` — per-targeton/per-reference metrics and composite rank
* `<out_prefix>.median_abs_z.png`, `.median_abs_z_functional.png`, `.composite_score.png`


---

### STEP SEVEN: Variant curation and overlap with the assay

The primary variant reference set for this project is a curated SMC1A truth set —
assembled from LOVD, GeneDx, DDD/DECIPHER, and per-publication extraction — specifically
because ClinVar's own disease-condition attribution turned out to be unreliable for this
gene (see STEP TEN). It separates the two clinically distinct SMC1A disorders, CdLS
(missense/altered-function) and DEE85 (protein-truncating/loss-of-function), on a
per-case rather than gene-level basis. [`Code/variant_curation/`](Code/variant_curation)
contains that curation pipeline and its join against this screen's results, run in order:

| script | purpose |
|---|---|
| `01_build_transcript_reference.py` | builds the transcript model, exon table and reference slice; verifies the model reconstructs the MANE protein |
| `02_ingest_lovd.py` | LOVD relational export to the observation ledger |
| `02b_verify_lovd_capture.py` | proves the public export is a superset of the authenticated web view |
| `03_resolve_variants.py` | resolves every variant description to GRCh38 coordinates through six guards |
| `04_ingest_genedx.py` | GeneDx de novo cohort |
| `04b_ingest_literature.py` | per-publication extractions |
| `04d_ingest_ddd.py` | DDD/DECIPHER cohort |
| `05_build_registry.py` | collapses observations to variants, assigns CdLS/DEE85 groups, applies inclusion criteria |
| `06_join_assay.py` | joins the curated set to the SGE screen results; decodes each oligo from its own sequence, detects the library's fixed background edits, and reports per variant whether it was screened |
| `07_make_data_dictionary.py` | regenerates a data dictionary from `smc1a_schema.py` |
| `08_audit_curated_set.py` | independent audit of the curated set; exits non-zero on failure |

`smc1a_lib.py` / `smc1a_schema.py` are shared library code, not run directly.

#### Data note

None of this code's inputs or outputs are in this repo: source publications are
copyrighted, the GeneDx/DDD-DECIPHER cohorts are access-controlled and carry real patient
identifiers, the LOVD export is a bulk relational download best re-fetched than
redistributed, and the SGE screen results/library design used by `06_join_assay.py` are
unpublished, pre-publication experimental data. This code is provided for methodological
transparency and reuse; see each script's own docstring for exactly what it expects as
input.

#### Requirements:
```bash
python -m venv .venv && .venv/bin/pip install -r Code/variant_curation/requirements.txt
export SMC1A_CONTACT="you@example.org"     # sent to the Ensembl REST API, as it requests
```
Network access is needed on first run for the Ensembl REST API; responses are cached so
later runs are offline and deterministic.

**Planned, not yet done:** formal clinical validation/calibration of the assay against this
curated set — sensitivity, specificity, and OddsPath, following the Brnich et al. 2019 /
ClinGen SVI framework — is intended but has not been carried out yet.

---

### STEP EIGHT: Intersect with gnomAD and plot

#### Background:
[`sge_gnomad_intersect_plot.py`](Code/pipeline/sge_gnomad_intersect_plot.py) intersects the SGE DESeq2/GMM-anchor results with gnomAD (exomes + genomes) for every assay variant across all targetons, and tests whether the depleted/enriched calls are under-represented in the population relative to "no impact" — a purifying-selection sanity check.

It mirrors the ClinVar intersect script's structure (`sge_clinvar_intersect.py`, STEP TEN — same `find_deseq2_files`/`load_meta_consequences`/`get_depletion_column` helpers, same merge key: `vcf_pos` + `vcf_ref` + `vcf_alt`), with two differences:
1. **Two** gnomAD sources per targeton (exomes + genomes), pooled per variant, rather than one ClinVar VCF.
2. The cross-targeton summary TSV includes **every** assayed variant (matched or not) — testing for absence requires the full denominator, not just the hits.

**Important:** this script assumes gnomAD v4.1 INFO field names (`AC`, `AN`, `AF`, `grpmax`, `AF_grpmax`), which have **not** been verified against the actual downloaded VCFs. Before running in earnest, check with `bcftools view -h <TID>_gnomad_exomes.vcf.gz | grep '^##INFO'` and adjust `--ac_field`/`--an_field`/`--af_field`/`--grpmax_af_field`/`--grpmax_pop_field` if the names differ (the script warns once per file if a configured field is never found).

#### Requirements:
* `pip install pandas numpy scipy matplotlib`
* Directory of per-targeton DESeq2/anchor TSVs (filenames starting with `Targeton_ID`)
* `*_meta_consequences.tsv` files
* `<Targeton_ID>_gnomad_exomes.vcf.gz` / `<Targeton_ID>_gnomad_genomes.vcf.gz` — generated by [`extract_gnomad_vcfs.sh`](Code/pipeline/extract_gnomad_vcfs.sh) (per-targeton, queried remotely from gnomAD v4.1 sites VCFs via `bcftools` + HTTPS range requests, no local download needed)
* `targeton_regions.tsv` — generated by [`extract_targeton_regions.py`](Code/pipeline/extract_targeton_regions.py) from a directory of `*_meta_consequences.tsv`/`*_meta.csv` pairs
* [`sge_gnomad_intersect_plot.py`](Code/pipeline/sge_gnomad_intersect_plot.py)

#### Running the script:
Run (example):
```bash
python Code/pipeline/sge_gnomad_intersect_plot.py \
  --deseq2_dir     all_deseq2_results/ \
  --meta_dir       /path/to/meta_consequence/files \
  --gnomad_vcf_dir gnomad_vcfs/ \
  --regions        targeton_regions.tsv \
  --outdir         gnomad_intersect_results/
```
By default only `FILTER=PASS` gnomAD records are kept; add `--include_non_pass` to include non-PASS records too.

#### Output:
* `<TID>_gnomad_annotated.tsv` — every DESeq2/anchor row for that targeton, plus gnomAD exomes/genomes/pooled columns
* `sge_gnomad_summary.tsv` — **all** assay variants across all targetons (the full denominator table), with gnomAD + `anchor_tier` columns
* `gnomad_absence_stats.tsv` — per-tier gnomAD-match rate, plus the depleted/enriched-vs-no-impact contingency test (Fisher's exact) results
* `gnomad_match_rate_by_tier.png` — bar chart of % gnomAD-matched per tier
* `gnomad_af_by_tier.png` — boxplot of log10(pooled AF) per tier, gnomAD-matched variants only

#### Notes:
* Verify gnomAD INFO field names against your actual VCFs before trusting the output (see Background above).
* The `TARGETON_EXON` dict (targeton → exon number) used for the cross-targeton summary is
  hardcoded here and **separately duplicated** in `sge_clinvar_intersect.py` (STEP TEN,
  same dict, same values) — there's no shared source of truth, so if the exon map ever
  changes, both copies need updating by hand.

---

### STEP NINE: In silico predictor concordance and other supporting analyses

A set of standalone analyses in [`Code/investigations/`](Code/investigations), each answering
one specific supporting question rather than being pipeline infrastructure — none of them feed
into a later step:

* [`insilico_concordance.R`](Code/investigations/insilico_concordance.R) — tests concordance
  between in-silico variant-effect predictors (SIFT, PolyPhen, CADD, REVEL, SpliceAI) and the
  assay's functional calls, gene-wide, stratified by consequence class (Spearman correlation,
  Kruskal-Wallis across `anchor_tier`, and AUC + sensitivity-at-Youden's-J treating `anchor_tier`
  as ground truth).
* [`cdls_vs_dee85_missense_predictors.R`](Code/investigations/cdls_vs_dee85_missense_predictors.R) —
  follow-up to `insilico_concordance.R`, restricted to the curated missense variants from
  STEP SEVEN: for `CdLS_pathogenic` vs `DEE85_pathogenic` missense variants specifically, do
  the predictor scores (SIFT, PolyPhen, CADD, REVEL) separate the two disease groups the way
  the assay's own depletion tier does? Answer: no — none of the four predictors distinguish
  the groups (Mann-Whitney p = 0.32-0.89; CADD nominally p=0.029 but on near-identical
  medians), while the assay itself sharply does (0% of CdLS missense strongly depleting vs
  64% of DEE85 missense). Predictors score "is this damaging", not "damaging by which
  mechanism" — this is a concrete demonstration of that distinction, using the disease
  groups this project's own curation established. Joins `06_join_assay.py`'s missense
  analysis table (STEP SEVEN) with `insilico_concordance.R`'s merged predictor table on
  oligo name + targeton.
* [`cassette_exon_consequence_analysis.py`](Code/investigations/cassette_exon_consequence_analysis.py) —
  tests whether cassette-exon-annotated exons show different depletion behaviour for synonymous
  and splice-region variants than non-cassette exons (Mann-Whitney, Fisher's exact).
* [`compare_shrinkage_k_runs.py`](Code/investigations/compare_shrinkage_k_runs.py) — sensitivity
  of the STEP FIVE classifier's tier calls to the shrinkage strength parameter `k`, which is
  auto-selected heuristically rather than fit from the data.

---

### STEP TEN: ClinVar

**Framing:** ClinVar was the original reference variant set for this project. This chain of
analyses is what surfaced that ClinVar's disease-condition field (CdLS vs DEE85) is frequently
unreliable for SMC1A — not a validation of the assay in its own right. It's presented here as
the reason a separately curated truth set exists (STEP SEVEN), rather than as this project's
primary evidence.

#### Extract per-targeton ClinVar VCFs

[`extract_clinvar_vcfs.sh`](Code/pipeline/extract_clinvar_vcfs.sh) generates the per-targeton ClinVar VCFs used below, from a whole-genome ClinVar VCF via `bcftools`, using the same `targeton_regions.tsv` as the gnomAD step.

```bash
bash Code/pipeline/extract_clinvar_vcfs.sh \
  --regions  targeton_regions.tsv \
  --clinvar  clinvar.vcf.gz \
  --outdir   clinvar_vcfs/
```

Requirements: `bcftools` (>=1.9), `tabix`. Note the 'chr' prefix is stripped from `chrom` before querying (ClinVar's GRCh38 VCF uses plain contig names) — the opposite convention from the gnomAD script.

#### Intersect with ClinVar and plot

[`sge_clinvar_intersect.py`](Code/pipeline/sge_clinvar_intersect.py) intersects the SGE anchor-tier results with ClinVar variants and produces per-targeton lollipop plots (shape = depletion tier, colour = disease condition, border = P/LP status) plus a condition-split plot (DEE-only vs CdLS-only, coloured by consequence).

```bash
python Code/pipeline/sge_clinvar_intersect.py \
  --deseq2_dir  all_deseq2_results.tsv \
  --meta_dir    /path/to/meta_consequence/files \
  --vcf_dir     /path/to/per_targeton_vcfs/ \
  --regions     targeton_regions.tsv \
  --outdir      clinvar_intersect_plots/
```
Add `--de_novo_only` to restrict the **plots** to ClinVar variants with the de-novo bit set (`ORIGIN`); this never affects the TSVs. Produces `clinvar_variants_summary.tsv`, the cross-targeton table everything below is built from.

Portability caveat, fine for this single-gene, chrX-only screen but worth knowing before
adapting the code elsewhere: the oligo-name-to-variant-key parser (`chrX:POS_REF>ALT`)
hardcodes `chrX` in its regex — it won't match oligo names for a gene on another
chromosome. Also carries its own copy of the `TARGETON_EXON` dict (targeton → exon
number) — see the STEP EIGHT note on why that's duplicated rather than shared.

#### Plot domain maps

[`sge_domain_map_clinvar.py`](Code/pipeline/sge_domain_map_clinvar.py) fetches SMC1A domain annotations from UniProt, maps them to genomic coordinates using a hardcoded SMC1A/hg38 exon table, and produces per-targeton lollipop + heatmap plots, optionally overlaying ClinVar P/LP variants.

The UniProt ID, Ensembl transcript, protein length, and the genomic↔AA exon coordinate
table are all hardcoded SMC1A/hg38 defaults (`--uniprot`/`--transcript` can override the
first two, but the exon table itself is SMC1A-specific code, not a flag) — expect to edit
the script, not just its arguments, before pointing it at a different gene.

```bash
python Code/pipeline/sge_domain_map_clinvar.py \
  --input SMC1A_maveqc/all_screens/thesis/GMM_shrinkage/D4_ref/results/*.tsv \
  --outdir SMC1A_maveqc/all_screens/thesis/GMM_shrinkage/D4_ref/plots/domain \
  --clinvar SMC1A_maveqc/all_screens/clin_var/clinvar.vcf.gz
```

#### Investigate the discordant cases

Discordant cases surfaced by the intersection above (functional score vs ClinVar
classification disagreeing) are followed up with a three-script chain, in order:

1. [`investigate_clinvar_variants.py`](Code/investigations/investigate_clinvar_variants.py) —
   pulls submission-level detail (review status, dates, submitters) from ClinVar's
   `submission_summary.txt.gz` for a filtered or explicit subset of variants from
   `clinvar_variants_summary.tsv`. Run with **no filters** to process every ClinVar-matched
   variant gene-wide, which the next two scripts require.
2. [`check_omim_coding_patterns.py`](Code/investigations/check_omim_coding_patterns.py) — checks,
   gene-wide, how often individual submissions actually cite a DEE85-specific vs CdLS-specific
   OMIM code/term, rather than relying on ClinVar's own aggregated `condition` column.
3. [`classify_both_mechanism.py`](Code/investigations/classify_both_mechanism.py) — for every
   variant labelled `condition == "Both"`, determines whether that's one submission citing both
   conditions together, or genuine cross-submitter disagreement.

Also: [`check_de_novo_agreement.py`](Code/investigations/check_de_novo_agreement.py), a minor
diagnostic checking that the VCF-level `ORIGIN` de-novo flag and the submission-level de-novo
flag never silently disagree.

#### ClinVar concordance

[`sge_clinvar_concordance.py`](Code/pipeline/sge_clinvar_concordance.py) quantitatively compares
SGE classification approaches against ClinVar P/LP vs B/LB calls (sensitivity/specificity/AUC),
stratified by CdLS vs DEE. Given the condition-labelling issues above, read this as a
demonstration of ClinVar's limits as ground truth for this gene — not as this project's primary
validation, which is the (currently in-progress) clinical calibration against the curated truth
set in STEP SEVEN.

#### Assay sensitivity check

[`nondepleting_plp_vs_synonymous.py`](Code/investigations/nondepleting_plp_vs_synonymous.py) —
are ClinVar P/LP variants the assay calls "no impact" truly indistinguishable from synonymous
controls, or is there a subtle residual signal below the calling threshold? Uses
`clinvar_variants_summary.tsv` from the intersection above.

---
