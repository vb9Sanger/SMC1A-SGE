# SMC1A SGE

Saturation genome editing of the entire *SMC1A* coding sequence in HAP1 cells: 27 targetons
across 25 coding exons, of which 25 targetons returned analysable data. This repository holds
the analysis code, in the order it runs.

---

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
| `09_calibrate_sensitivity_oddspath.py` | formal sensitivity/specificity/OddsPath calibration against the curated set (see below) |
| `10_calibrate_by_tier.py` | the same calibration per depletion tier rather than as a binary call, across four benign reference compositions, reporting both zero-cell corrections (see below) |

`smc1a_lib.py` / `smc1a_schema.py` are shared library code, not run directly.

#### One consequence vocabulary (`smc1a_lib.py`, pipeline v0.55)

The source meta-consequence files carry three consequence columns, and mixing them
was a real source of error, so the rule is fixed in one place:

| column | what it is | used? |
|---|---|---|
| `Consequence` | raw VEP string, most specific, **does not prioritise** when a variant has several terms | no — selecting on it over-includes (e.g. `startswith("missense")` matches `missense_variant,splice_region_variant`) |
| `Summary_Consequence` | the prioritised call, specific | **yes — canonical** |
| `Summary_Plot` | the same prioritisation, but collapses `Nonsense_Variant` + `Frameshift_Variant` into `LOF` | no |

`Summary_Plot` is not used because **`LOF` is a misnomer** — loss of function is a
mechanism and a missense variant can be LOF; what that label actually denotes is
*PTV*. It also breaks concretely: `10_calibrate_by_tier.py`'s `SUMMARY_TO_CLASS`
map keys on `Nonsense_Variant` and `Frameshift_Variant`, which exist only in
`Summary_Consequence`.

`CONSEQUENCE_ORDER` in `smc1a_lib.py` reproduces `Summary_Consequence`'s
prioritisation (splice above every protein-level consequence) while staying
finer-grained — `Others` resolves to seven specific classes and `Splice_Variant`
to three. It is validated 32/32 against the source column.

**PTV status is a separate field, not derived from the class.** With splice
outranking everything, a `stop_gained,splice_region_variant` variant classes as
`splice_region`, which is not a PTV class — deriving PTV status from the class
would silently drop 81 library variants. `is_predicted_ptv(so_terms)` reads the
SO terms directly instead:

```python
PTV_SO_TERMS = {"stop_gained", "frameshift_variant", "splice_acceptor_variant",
                "splice_donor_variant", "start_lost", "transcript_ablation",
                "feature_truncation", "feature_elongation"}
```

The class answers *what category is this*; the flag answers *is this a predicted
PTV*. Different questions, different fields. Any new analysis that selects on
consequence must use `consequence_class` (curated/join tables) or
`Summary_Consequence` (source meta and ClinVar/gnomAD summary tables) — never
`Consequence` or `Summary_Plot`.

Four further scripts are verification rather than pipeline stages, run on demand rather
than in sequence. They were written during a full audit of the curated set and are kept
because each one caught something:

| script | what it checks |
|---|---|
| `audit_attribution_evidence.py` | what independent, non-suppressed, per-case disease claim actually stands behind each curated variant's arm assignment, and whether any source claims the other arm. Found that 22 attributions rested on a single genotype-ascertained or rule-inferred source, which led to the D142/D144 exclusions |
| `verify_extractions_against_sources.py` | checks every extracted cDNA and protein description against its source PDF, XLSX or DOCX, main text and supplementary, notation-insensitively. Found D143 — a 15 bp in-frame deletion recorded as a 1 bp frameshift — and confirmed the curation silently corrects six errors in its own sources |
| `audit_all_variants_rederivation.py` | re-derives protein change, protein position and group assignment for **all** variants rather than the curated subset, independently of what the pipeline recorded. The other 153 feed the benign reference and the ClinVar comparisons, so an error there reaches the headline statistics without touching the curated set |
| `check_deposit_individual_matching.py` | whether each database deposit (LOVD, GeneDx, DDD) can be matched one-to-one onto a named individual in the publication it came from, which is what decides whether two records are one case or two. Open question Q38 |

#### Two things the step numbering doesn't tell you

**`annotate_domain_assignment.py` runs between 05 and 06, not after
everything.** Step 05 does not write the `domain_assignment` / `in_atpase_head` /
`in_hinge` / `in_reported_cluster` columns, nor the `literature_domain_hotspot*` ones
retained alongside them; that investigation script (STEP NINE) adds them afterwards, to *both*
`smc1a_variants_curated.tsv` and `smc1a_variants_all.tsv`. A registry rebuild drops them
silently — nothing errors, the columns are simply gone, and `06_join_assay.py` then carries
a narrower table through to the join. Re-run it after any rebuild of step 05, before step
06. The full order is `01` → `05` → annotate → `06` → `07` → `08`.

**Which targetons are screened is determined by the classifier output, not by a flag.** A
targeton counts as screened iff it has a results file in the classifier's output directory:
25 of the 27 designed targetons, i.e. all except **CQEJ** (exon 6, second tiling window) and
**NLVE** (exon 23), which are designed and present in the library but not screened. Those
two carry `screening_status = not_screened` in `01_build_transcript_reference.py`, and
variants whose only design window is one of them get `assay_testable_now = False`. Note that
`targeton_screening_status` is pipe-joined over a variant's targetons, so a variant spanning
a screened and an unscreened targeton reads `not_screened|screened` — split on `|` rather
than substring-matching, since `not_screened` contains `screened`. `08_audit_curated_set.py`
reports such variants without failing: an unscreened targeton is a true state of the world
that curation cannot fix, unlike a missing design window, which remains a hard failure.

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

#### Sensitivity/specificity/OddsPath calibration

`09_calibrate_sensitivity_oddspath.py` computes, from `06_join_assay.py`'s output: per-group
depletion rate with Clopper-Pearson 95% CIs, and, per pair of groups (group_a = pathogenic/case
reference, group_b = benign/control reference), Brnich et al. 2019's own OddsPath quantities:
LR+ = sensitivity/(1-specificity) for an abnormal (depleting) result, tiered against the PS3
thresholds, and LR- = (1-sensitivity)/specificity for a normal result, tiered against the BS3
thresholds (exact reciprocals of the PS3 ones) — each with a log-scale 95% CI and Fisher's exact
p on the underlying 2x2 table. The pooled case-control odds ratio (DOR = LR+/LR-) is also
reported, as a supplementary omnibus statistic, but is *not* itself OddsPath and is not tiered
against Table 3 (an earlier version of this script tiered the DOR directly, which is only a good
approximation of OddsPath when specificity is high; see the script's own docstring for the
Bayes-factor derivation showing why LR+/LR- are the paper's actual quantities). Computed
separately per `curation_group`
(default: `CdLS_pathogenic`, `DEE85_pathogenic` vs `Exclude_Benign`) rather than pooling all
pathogenic variants together, since CdLS and DEE85 differ in mechanism (dominant-negative vs
loss-of-function) and a single pooled figure would blend two very different evidence
strengths into one misleading average.

```bash
python Code/variant_curation/09_calibrate_sensitivity_oddspath.py \
  --join_tsv assay_join_all.tsv \
  --outdir sensitivity_results/
```

De-duplicates to one row per `variant_key` before computing anything (`--on_duplicate`
controls how a variant matched by more than one oligo with discordant calls is resolved —
see the script's own docstring, in particular the note on `position_only_verify_manually`
rows needing manual verification before use, not automatic resolution).

Three further options widen what a single run can compare, each reported alongside —
never instead of — the base per-group numbers:

* `--extra_group_tsv NAME=path.tsv` (repeatable) — adds a group read from any plain
  `(variant_key, anchor_call)` TSV instead of `--join_tsv`'s `curation_group` column, e.g.
  additional benign controls from ClinVar (see
  [`extract_clinvar_benign_controls.py`](Code/investigations/extract_clinvar_benign_controls.py),
  STEP TEN) rather than only the curated set's own small `Exclude_Benign` group. A group
  is treated as playing the benign/control role in reporting if `"benign"` appears in its
  name (case-insensitive).
* `--pool NAME=GROUP1,GROUP2,...` (repeatable) — adds a virtual group whose count is the
  sum of already-computed groups' counts, e.g. `Combined_PLP=CdLS_pathogenic,DEE85_pathogenic`
  or a pooled curated+ClinVar benign reference — the mechanism for "what if you combine
  the two mechanisms/sources," done transparently rather than by hand.
* `--consequence_filter CLASS [CLASS ...]` — restricts every group (pathogenic, benign,
  and any `--extra_group_tsv`, which must then also carry a `consequence_class` column) to
  the given consequence classes before counting, e.g.
  `--consequence_filter missense inframe_deletion inframe_insertion` to calibrate
  specifically on variants where annotation alone can't already call pathogenicity —
  excluding PTVs, whose classification often doesn't need functional evidence in the first
  place (ACMG PVS1), so a calibration that includes them risks non-independent,
  double-counted evidence if later stacked with a PVS1 call for the same variant.

Actual sensitivity/specificity/OddsPath figures belong in the thesis write-up, not this
README, per the "no results in this repo" policy elsewhere in this document.

#### Calibration by tier, and the benign reference

[`10_calibrate_by_tier.py`](Code/variant_curation/10_calibrate_by_tier.py) extends the above
in two ways that turned out to interact, and reports both corrections for a zero cell.

**Per-tier likelihood ratios instead of a binary call.** `09_...py` collapses `anchor_tier`
to depleted / not depleted, which pools "weakly depleting" with "strongly depleting". Those
are not equivalent evidence, and pooling lets the confident calls carry the uncertain ones.
This computes LR(tier) = P(tier | pathogenic) / P(tier | benign) for each observable result,
which is what OddsPath reduces to per result — the binary LR+ being the special case where
the result is "depleted". In this dataset the difference is not academic: for one arm,
binarising promotes the evidence a whole tier, because the weakly-depleting variants enter
with an empty benign comparator.

**A benign reference that carries no clinical assertion.**
[`extract_gnomad_benign_controls.py`](Code/investigations/extract_gnomad_benign_controls.py)
builds population-frequency controls from the assay's own gnomAD intersection. For a severe,
early-onset, essentially always de novo disorder, a pathogenic allele has close to zero
reproductive fitness and should not persist in a population database, so a variant seen
repeatedly in gnomAD is evidence against pathogenicity that is independent of whether anyone
classified it correctly. Running the calibration with and without the ClinVar-derived controls
then answers directly how far the result depends on trusting ClinVar.

The frequency threshold is derived rather than conventional, using the maximum credible
population allele frequency of Whiffin et al. 2017 (*Genet Med* 19:1151–1158) — ACMG's BA1 (5%)
is far too permissive for a disorder of this severity and BS1 has no default. For SMC1A/CdLS
that value is 5.0e-7, which turns out to sit **below the resolution of gnomAD**: at the observed
AN of ~1.1M a single allele is already 9.2e-7, so every variant observed at all exceeds the
maximum credible pathogenic frequency.

The script offers two rules on that basis, and the choice between them is the substance of the
decision rather than a parameter:

- `--proxy_clinical` — **the default for this work.** Every gnomAD-observed variant becomes a
  benign control, per CanVIG-UK Recommendation 4, which extends the maximum-tolerated-allele-
  frequency principle that far for very rare, highly penetrant, early-onset disease.
- `--max_credible_af 5.0e-7` — Whiffin's rigorous form, filtering on the *filtering* allele
  frequency (the CI lower bound) rather than the point estimate. The script derives the
  necessary allele count rather than taking one: it finds the smallest AC whose one-sided 95%
  lower bound clears maxAF, which here is AC ≥ 3. Retained as a robustness comparator.

Two further filters apply in both modes and each was added after seeing what the data does
without it: deduplicate by genomic variant (the gnomAD table is per oligo and SMC1A's targeton
windows overlap), and exclude anything asserted pathogenic elsewhere. Note that under the
AC ≥ 3 rule that last filter removes nothing — no variant asserted pathogenic in the curated
arms or in ClinVar reaches AC ≥ 3 in gnomAD.

**Zero benign cells.** Several comparisons have no depleting benign control at all, and the
gene is the reason — SMC1A is missense-constrained, so benign missense variants barely exist
and no amount of further *clinical* curation produces a larger matched reference (the curated
benign group contributes 2 of them, and all of ClinVar adds 9). Where that happens the
maximum likelihood estimate of LR+ is infinite and any finite figure comes from a correction,
not from data. The script therefore reports **both** the Haldane-Anscombe convention `09_...py`
uses (+0.5 to every cell) and the correction Brnich et al. themselves use (+1 misclassified
variant per set), because the two can differ by a factor of two and straddle a tier boundary.

It also reports a **`fisher_p`** column, and that is what settles significance: Fisher's exact
test on the raw 2×2 table is the only quantity here needing no zero-cell correction at all. It
is one-sided in the direction under test, so the `no impact` rows are tested in the benign
direction rather than the depleting one. The reporting convention built on this — gate on
`fisher_p`, tier the **Brnich point estimate** per that paper's own Table 3, quote the 95% CI
alongside, and give no evidence tier to a bound — is recorded as **D159**. It supersedes the
earlier position that the confidence bound was the defensible statement to tier; tiering a
bound turned out to be this project's own extension rather than anything the source frameworks
do. The argument is set out in full in the thesis write-up (not in this repo, per the results
policy above) and recorded as decisions **D152** (which correction), **D153** (which benign
reference, and why) and **D159** (how all of it is reported). Enlarging the reference turned out to matter more than the choice of
correction: on the matched cut the two corrections disagreed across a tier boundary against
the smallest reference and agree against the largest.

One related fix lives in the code rather than the write-up. Where *both* cells are zero — no
pathogenic variant occupies a tier and no benign one does either — a Haldane correction returns
a large ratio purely from n₂ ≫ n₁, and the table would label a tier nothing occupies. `lr_ci`
and `lr_brnich` now return `None` whenever the pathogenic cell is zero, and the tables print a
dash.

```bash
# the benign controls, in the mode used for the primary analysis
python Code/investigations/extract_gnomad_benign_controls.py \
  --gnomad_summary sge_gnomad_summary.tsv \
  --curated_tsv smc1a_variants_all.tsv \
  --clinvar_summary clinvar_variants_summary.tsv \
  --proxy_clinical --output gnomad_proxy_clinical_controls.tsv

python Code/variant_curation/10_calibrate_by_tier.py \
  --join_tsv assay_join_all.tsv \
  --clinvar_benign clinvar_benign_controls.tsv \
  --gnomad_benign gnomad_proxy_clinical_controls.tsv \
  --consequence_map clinvar_variants_summary.tsv sge_gnomad_summary.tsv \
  --outdir calibration/ --label incl_PTV
```

`--consequence_filter` behaves as in `09_...py` but applies to the external control files too,
which is why `--consequence_map` is needed: most controls are not in the curated join and
carry no `consequence_class`, so without a class for them a matched run silently drops them
rather than matching on nothing.

#### Pilot screen

[`pilot_screen_analysis.py`](Code/investigations/pilot_screen_analysis.py) tabulates the
pilot SGE screen of two exons, which was run ahead of the full screen to establish that the
assay discriminates variant classes.

**Input:** the pilot workbook of per-library, per-oligo DESeq2 results (`--pilot_xlsx`), plus
the exon map (`--exon_map`) for exon spans. It takes no full-screen input: the pilot LFCs are
not adjusted for positional bias, which `pos_adj_log2FoldChange_raw` does account for, so
pilot and full-screen scores are not combined anywhere in this pipeline.

**Output:** a per-variant table with LFC at each timepoint and the Day 15 call; depletion rate
and median LFC per consequence class within each HDR library; the per-timepoint trajectory of
protein-truncating against synonymous variants; and a per-library QC flag for libraries
yielding no significance calls. Consequence classes follow `Summary_Consequence`, as
everywhere else here.

#### Per-residue mutational constraint

[`smc1a_mii.py`](Code/investigations/smc1a_mii.py) summarises the per-variant functional map
as a per-residue constraint index, stratified by protein domain, and
[`smc1a_mii_curated_overlay.py`](Code/investigations/smc1a_mii_curated_overlay.py) overlays the
curated CdLS/DEE85/benign variants onto it.

**Two indices per residue**, both collapsed to distinct amino acid substitutions first:
**mean LFC** (the mean assay score across substitutions at that residue — the standard
"mutational tolerance" of the deep-mutational-scanning literature, and the primary index here)
and **MII** (the fraction of substitutions classified as functionally abnormal, following the
definition in Jaramillo Calle's 2025 Cambridge thesis §5.3.4.6, applied there to CTCF). The two
are strongly correlated and agree on essentially every conclusion; mean LFC leads because it
retains the magnitude of each substitution's effect, where MII thresholds each substitution and
discards it.

Collapsing to substitutions rather than variants matters for this library specifically: it
includes multi-nucleotide codon replacements, so one amino acid substitution is often reached by
several genomic variants, and how many reach it depends on the wild-type codon rather than on
the residue. Counting variants would weight a substitution by how many codon paths reach it.
The collapse also yields a free reproducibility statistic, written to the run summary: the rate
at which independently encoded routes to the same substitution agree on depleting vs not.

Domain boundaries are fetched live from UniProt and InterPro and cached, never hard-coded. The
ATPase head comes from InterPro IPR028468, which annotates it as **two** fragments because the
domain is discontinuous in sequence — UniProt has no equivalent feature, since a UniProt
`Domain` carries a single contiguous range and the SMC head is formed by the N- and C-termini
folding together. The script exits rather than guessing if InterPro does not return exactly two
fragments.

The overlay uses **leave-one-out** indices throughout: each curated variant contributes ~1/19 of
its own residue's index, and the analysis tests exactly that direction, so the variant's own
substitution is removed before comparing. The effect turns out to be negligible (≤0.02 on
either index) but it is removed rather than argued away.

#### A second calibration method, for contrast

[`build_excalibr_input.py`](Code/investigations/build_excalibr_input.py) supports a calibration of
the same assay by **ExCALIBR** (Zeiberg, Stewart et al. 2026, bioRxiv 2025.04.29.651326;
`github.com/rosstewart/exCALIBR`), which is run from its own repository rather than from here.
Where the Brnich scripts tier a categorical call, ExCALIBR models the four score distributions
(P/LP, B/LB, gnomAD, synonymous) as skew-normal mixtures and returns a per-variant posterior in
ACMG evidence points.

The builder writes one input per disease arm, because SMC1A's two disorders act by different
mechanisms and a single pathogenic sample would be bimodal — which is also the extension the
paper flags as unevaluated. Three things in it are judgement calls rather than format
requirements, and the docstring records each: splice-region synonymous variants are excluded
from the functionally-normal reference (they deplete at 11.4% against 0.9% elsewhere, because
SGE reports splicing), no SpliceAI filter is applied (that filter exists for cDNA assays), and
curated-pathogenic variants seen in gnomAD keep both labels.

It also uses the VCF triple rather than parsing `oligo_name`. The `chrX:pos_REF>ALT` pattern
used elsewhere in this repo matches SNVs only — deletions are named `chrX:53380079_1del` — so
keying on it drops 25,279 of 39,135 assayed rows. Anything else in this repo that parses
variants out of oligo names inherits that limitation.

[`excalibr_pathomechanism.py`](Code/investigations/excalibr_pathomechanism.py) recovers
P(M=1|Y=1) — the fraction of a disease arm's pathogenic variants whose mechanism the assay
measures — which ExCALIBR exposes as `--pathomechanism-prior` but which cannot be obtained that
way for a gene whose population prior falls below its trust floor: the reported quantity is
P(Y=1) × P(M=1|Y=1), so a missing first factor makes it NaN, and the `--manual-prior` escape
raises a NameError when combined with it. P(M=1|Y=1) itself does not depend on the prior, so
this calls their own boundary estimator directly against fits the pipeline already wrote.

**Reconciling the two calibration scripts.** Both normalise `variant_key` before counting
(the three benign sources disagree on whether to write the contig, so without it no ClinVar key
ever equals a curated or gnomAD one), and both take `--pool` as a union by variant_key rather
than a sum of counts — the benign sources share 254 variants, so summing inflated the
denominator. Both of those were defects in `09_...py` until 2026-09-28.

One difference remains, by design: `10_...py` restricts to the three depletion tiers and so
drops `enriched` controls, while `09_...py` keeps them in the denominator on the grounds that
an enriched call is not a depleted call and therefore contributes to specificity. For the
curated + ClinVar + gnomAD reference that is n=1408 against n=1402, a difference of exactly the
6 enriched controls, with the same 6 depleted in both.

---

### STEP EIGHT: Intersect with gnomAD and plot

#### Background:
[`sge_gnomad_intersect_plot.py`](Code/pipeline/sge_gnomad_intersect_plot.py) intersects the SGE DESeq2/GMM-anchor results with gnomAD (exomes + genomes) for every assay variant across all targetons, and tests whether the depleted/enriched calls are under-represented in the population relative to "no impact" — a purifying-selection sanity check.

It mirrors the ClinVar intersect script's structure (`sge_clinvar_intersect.py`, STEP TEN — same `find_deseq2_files`/`load_meta_consequences`/`get_depletion_column` helpers, same merge key: `vcf_pos` + `vcf_ref` + `vcf_alt`), with two differences:
1. **Two** gnomAD sources per targeton (exomes + genomes), pooled per variant, rather than one ClinVar VCF.
2. The cross-targeton summary TSV includes **every** assayed variant (matched or not) — testing for absence requires the full denominator, not just the hits.

The script reads gnomAD v4.1 INFO field names (`AC`, `AN`, `AF`, `grpmax`, `AF_grpmax`). All five were checked against the headers of the downloaded per-targeton VCFs and are present in both the exomes and the genomes files. If adapting this to another gnomAD release, re-check with `bcftools view -h <TID>_gnomad_exomes.vcf.gz | grep '^##INFO'` and adjust `--ac_field`/`--an_field`/`--af_field`/`--grpmax_af_field`/`--grpmax_pop_field` — the script warns once per file if a configured field is never found, but a renamed field would otherwise read as absent data rather than as an error.

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

#### Representation by consequence class

The pooled depleted-vs-no-impact test above is confounded, and
[`gnomad_representation_by_consequence.py`](Code/investigations/gnomad_representation_by_consequence.py)
is the comparison to read instead. Baseline gnomAD representation varies by more than
twentyfold across consequence classes — synonymous and intronic variants are well
represented, missense variants barely at all — and the assay's depleting calls are
concentrated in precisely the classes where that baseline is lowest. A pooled odds ratio
therefore partly measures the consequence composition of the depleting set rather than the
effect of depletion, and overstates the contrast. Stratifying holds the baseline fixed.

It runs from `sge_gnomad_summary.tsv` alone and emits both, so the size of the confound is
visible rather than hidden. Classes where almost every variant is classified as depleting
leave no within-class comparison pool; those rows are computed and flagged
`no_comparison_pool` rather than dropped, since their absence from a reported table is
itself something to be able to justify.

```bash
python Code/investigations/gnomad_representation_by_consequence.py \
  --gnomad_summary gnomad_intersect_results/sge_gnomad_summary.tsv \
  --outdir gnomad_intersect_results/
```

`--systematic_only` repeats the analysis over variants reached by the systematic mutators
alone. Variants already reported in gnomAD v2.1/v3.1 were deliberately incorporated at the
HDR oligo design stage, so a fraction of the matches arrived by that route rather than by
the assay independently sampling population variation; those variants are overwhelmingly
intronic, and the intronic row is the only one that moves materially. The flag exits rather
than returning a null result if it matches no oligo-name tag, since a silently empty filter
would read as a clean sensitivity check.

#### Notes:
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
  the assay's own depletion tier does? Predictors score "is this damaging", not "damaging by
  which mechanism", and this tests that distinction directly, using the disease groups this
  project's own curation established rather than a database's condition field. Joins
  `06_join_assay.py`'s missense analysis table (STEP SEVEN) with `insilico_concordance.R`'s
  merged predictor table on oligo name + targeton. Figures belong in the write-up, not here,
  per the results policy in STEP SEVEN.

  It also records **predictor coverage** across the missense library, which is a fact about
  the library rather than a comparison: SIFT and PolyPhen score every assayed missense
  variant, while CADD and REVEL score single-nucleotide substitutions only. This library is
  largely multi-nucleotide codon replacements, so their coverage is a minority of it.
* [`build_predictor_context.py`](Code/investigations/build_predictor_context.py) — assembles
  per-variant **orthogonal** annotation to set against the assay's calls, where the
  predictors above are not orthogonal at all. Four sources chosen for independence rather
  than performance: CADD (a trained pathogenicity predictor, and so partly circular with
  clinical labels), phyloP per-base conservation from the UCSC 100-way vertebrate alignment,
  relative solvent accessibility retrieved per residue from AlphaSync, and predicted change
  in folding free energy from a ThermoMPNN site-saturation scan of the AlphaFold model of
  Q14683. phyloP, RSA and ddG know nothing about pathogenicity, which is the point: CADD
  agreeing with the assay is weak evidence, the other three agreeing is not.

  Three things the script handles rather than assumes. RSA and pLDDT come from AlphaSync's
  stored residue features, which for this entry were computed from an earlier AlphaFold DB
  release than the model the ddG scan runs on, so the two are **not guaranteed to describe
  identical coordinates** -- the docstring records where they diverge and by how much. ddG is
  joined by parsing protein HGVS with the wild-type residue carried into the join key, so a
  mismatch arrives as missing rather than silently picking up another substitution's score.
  pLDDT is carried through per residue so that low-confidence regions of a predicted model
  can be excluded downstream rather than quietly weighted in.
* [`annotate_domain_assignment.py`](Code/investigations/annotate_domain_assignment.py) —
  a review of the CdLS/DEE85 mechanism literature turned up a specific,
  independently-replicated claim (Baranano et al. 2022; Bozarth et al. 2023; Di Nardo
  et al. 2026): the small non-truncating (missense/in-frame) subset of DEE85 variants
  clusters in SMC1A's N-/C-terminal ATPase head and central hinge, rather than being
  scattered like the truncating majority. This script assigns every curated variant to
  a domain and flags whether it falls in that region, adding `domain_assignment`,
  `in_atpase_head`, `in_hinge` and `in_reported_cluster` to
  `smc1a_variants_curated.tsv`/`smc1a_variants_all.tsv`.

  Boundaries come from `smc1a_domain_annotation.json` — the same UniProt/InterPro
  annotation `smc1a_mii.py` uses, so the two analyses speak one structural language.
  **This replaces an earlier publication-based approach**, which marked a variant
  as being in the cluster only where a named publication explicitly placed it there.
  That could only ever mark True a variant some paper had already discussed, so a
  variant inside the head that nobody had written about was marked False; testing the
  literature's claim with a method that can only recover what the literature named is
  close to circular (**D161**). The publication-based columns are retained as
  provenance, not as the test. `VARIANT_COLS` in `06_join_assay.py` (STEP SEVEN)
  carries all of these through to `assay_join_all.tsv`, which is why this script runs
  *before* step 06 and must be re-run after any rebuild of step 05 (see STEP SEVEN's
  run-order note).
* [`cassette_exon_consequence_analysis.py`](Code/investigations/cassette_exon_consequence_analysis.py) —
  tests whether cassette-exon-annotated exons show different depletion behaviour for synonymous
  and splice-region variants than non-cassette exons (Mann-Whitney, Fisher's exact). Its
  `--exon_events_tsv` input (which exons are UCSC-annotated cassette exons) is built by
  [`build_smc1a_exon_splice_events_table.py`](Code/investigations/build_smc1a_exon_splice_events_table.py)
  from UCSC's "Alt Events" track (`knownAlt`), fetched live via the UCSC REST API and
  intersected against each exon's genomic boundaries — no manual UCSC Table Browser export
  needed. Synonymous variants in the same exons are compared in parallel as a specificity
  control: they do not perturb a splice signal, so a difference among them would mean these
  exons simply tolerate more or less variation overall, rather than being specifically
  sensitive to splice disruption.

  Three tests on the splice-region side, reported as Tables 4–6: Mann-Whitney on the
  continuous z-scores, Fisher's exact on the proportion depleting, and — added 2026-10-08 —
  the same comparison with the **exon as the unit of analysis**. The first two pool
  thousands of variants and so assume an independence that variants sharing an exon's splice
  context do not have; Table 6 repeats the contrast on one figure per targeton, which has far
  less power and is the conservative reading to quote.
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

#### Investigate the discordant cases

Discordant cases surfaced by the intersection above (functional score vs ClinVar
classification disagreeing) are followed up with a two-script chain, in order:

1. [`investigate_clinvar_variants.py`](Code/investigations/investigate_clinvar_variants.py) —
   pulls submission-level detail (review status, dates, submitters) from ClinVar's
   `submission_summary.txt.gz` for a filtered or explicit subset of variants from
   `clinvar_variants_summary.tsv`. Run with **no filters** to process every ClinVar-matched
   variant gene-wide, which the next script requires.
2. [`check_omim_coding_patterns.py`](Code/investigations/check_omim_coding_patterns.py) — checks,
   gene-wide, how often individual submissions actually cite a DEE85-specific vs CdLS-specific
   OMIM code/term, rather than relying on ClinVar's own aggregated `condition` column.

#### Assay sensitivity check

[`nondepleting_plp_vs_synonymous.py`](Code/investigations/nondepleting_plp_vs_synonymous.py) —
are ClinVar P/LP variants the assay calls "no impact" truly indistinguishable from synonymous
controls, or is there a subtle residual signal below the calling threshold? Uses
`clinvar_variants_summary.tsv` from the intersection above, so it lives here rather than under
STEP NINE's general supporting analyses, despite testing the assay's own sensitivity rather
than ClinVar's reliability.

#### Widening the STEP SEVEN calibration with ClinVar

Two further scripts connect this step back to STEP SEVEN's calibration, using ClinVar data for
purposes that don't depend on its disease-condition attribution being reliable:

* [`extract_clinvar_benign_controls.py`](Code/investigations/extract_clinvar_benign_controls.py) —
  pulls every ClinVar Benign/Likely benign call for the gene region, joined to the assay's own
  `anchor_tier` (via `clinvar_variants_summary.tsv`) and to ClinVar's `CLNREVSTAT` review-status
  field (fetched directly from the VCF, since the intersect script above never captures it), for
  use as an `--extra_group_tsv` benign reference in `09_calibrate_sensitivity_oddspath.py` —
  widening it well beyond the curated set's own small `Exclude_Benign` group. Defensible because
  disease-condition attribution (what this step found unreliable) doesn't apply to a benign call;
  `--min_review_tier multi_submitter` is available for a stricter, smaller, higher-confidence cut.
* [`reclassify_clinvar_vus.py`](Code/investigations/reclassify_clinvar_vus.py) — applies the
  calibrated OddsPath tiers to the gene's actual ClinVar Uncertain-significance calls. Mapping
  is **by depletion tier**, reading the per-tier calibration (`--ptv_splice_tsv`,
  `--missense_tsv`) rather than a single binary depleted-vs-not odds ratio, so a strongly
  depleting and a weakly depleting call do not collapse to the same evidence strength.

  Three outcomes, and the distinctions between them are the point of the script:

  - **PTV and splice-class VUS that deplete** take the tier calibrated on the arm given by
    `--ptv_splice_arm`.
  - **Missense VUS that deplete** are *not* assigned a disease-specific tier. Missense is the
    predominant class in one disorder but variants of the other are also missense, so the arm
    cannot be inferred from consequence class alone — the routing the earlier version did.
    They are emitted as `observation_depleting_missense` instead: an observation of functional
    abnormality, carrying no ACMG/AMP code. `--missense_policy route` restores the old
    behaviour for comparison, but `observe` is the default and the documented choice.
  - **Everything else** — in-frame indels, enriched calls, non-depleting calls of any class —
    is left unmapped and flagged for manual review rather than forced into a tier.

  Two gates apply before any tier is borrowed: `--alpha` on the cell's Fisher exact p, and
  `--min_k_path`, a minimum number of pathogenic variants in the calibrating cell. The latter
  is a judgement made here rather than a published rule, and the script's header says so.

---
