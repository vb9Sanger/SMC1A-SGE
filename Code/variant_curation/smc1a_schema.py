"""
smc1a_schema.py -- single source of truth for output column order and meaning.

Both the column ordering of the delivered tables and `docs/DATA_DICTIONARY.md`
are generated from the structures below, so the documentation cannot drift out
of step with the files.

Each table is a list of sections; each section is
    (section_title, section_note, [(column, type, description), ...])

Columns produced by the pipeline but absent from this schema are still written
(appended at the end, under "Unclassified") and are listed in the dictionary as
undocumented -- so a new column is surfaced rather than silently dropped.
"""

from __future__ import annotations

# ==========================================================================
# Variant-level tables: smc1a_variants_all.tsv and smc1a_variants_curated.tsv
# ==========================================================================

VARIANTS = [

    ("Variant identity (GRCh38, plus strand)",
     "The join key and the canonical description. `pos`/`ref`/`alt` are "
     "**plus-strand** GRCh38, matching VCF and the SGE oligo names, so "
     "`variant_key` joins directly to the screen results. SMC1A is on the "
     "minus strand, so these alleles are the complement of the bases named in "
     "the c. description.",
     [
      ("variant_key", "string", "Canonical identifier `chrX:<pos>:<REF>:<ALT>`. Primary key. `UNRESOLVED::<desc>` where no coordinate could be assigned."),
      ("chrom", "string", "Chromosome, always `chrX`."),
      ("pos_grch38", "int", "1-based plus-strand GRCh38 position."),
      ("ref", "string", "Plus-strand reference allele."),
      ("alt", "string", "Plus-strand alternate allele."),
      ("hgvs_g", "string", "HGVS genomic description on NC_000023.11."),
      ("hgvs_c_mane", "string", "HGVS coding description on ENST00000322213.9 (MANE Select)."),
      ("hgvs_p_mane", "string", "HGVS protein description on ENSP00000323421.3."),
      ("hgvs_c_as_submitted", "list", "Every distinct c. description any source used, `|`-separated. Retained so the original wording is never lost."),
      ("transcript", "string", "Reference transcript, always ENST00000322213.9."),
      ("protein", "string", "Reference protein, always ENSP00000323421.3."),
      ("assembly", "string", "Genome build, always GRCh38."),
     ]),

    ("Molecular consequence",
     "Recomputed on the MANE transcript by Ensembl VEP; never taken from the "
     "source. `ref_aa` is derived independently from the cached reference "
     "sequence and is what guard G5 checks against.",
     [
      ("protein_position", "int", "1-based codon number."),
      ("ref_aa", "string", "Reference amino acid at that codon, recomputed from the cached GRCh38 sequence."),
      ("edit_type", "string", "HGVS edit class: `sub`, `del`, `dup`, `ins`, `delins`, `inv`, `ident`."),
      ("consequence_class", "string", "Harmonised class: `missense`, `synonymous`, `nonsense`, `frameshift`, `inframe_deletion`, `inframe_insertion`, `splice_donor`, `splice_acceptor`, `splice_region`, `intronic`, `utr5`, `utr3`, `start_lost`, `stop_lost`, `cnv_*`."),
      ("consequence_terms", "list", "Raw VEP SO terms, `;`-separated."),
      ("is_predicted_ptv", "bool", "True for nonsense, frameshift, canonical splice, start-loss and CNV classes -- the classes SGE depletion reports on most directly."),
      ("exon", "int", "Exon number (1-25) if exonic."),
      ("intron", "int", "Intron number if intronic (intron *n* lies between exons *n* and *n*+1)."),
     ]),

    ("SGE assay context",
     "Three nested targeton regions are distinguished because conflating them "
     "gives the wrong answer to 'is this variant in the library?'. The "
     "**design window** is the operative one. `in_sge_library` is NOT inferred "
     "from these flags -- it is determined empirically by joining to the oligo "
     "list in step 06.",
     [
      ("targetons_design", "list", "Targeton(s) whose **design window** contains the variant -- i.e. where oligos were systematically designed. The operative field."),
      ("targetons_amplicon", "list", "Targeton(s) whose sequenced amplicon contains the variant. Wider than the design window; amplicon-only means sequenced but no oligo designed."),
      ("in_design_window", "bool", "True if the variant falls in at least one design window."),
      ("targeton_screening_status", "string", "`screened`, `not_screened` (CQEJ and NLVE -- designed and in the library, but not yet screened, so no results exist), or `not_in_design_window`. Pipe-joined over a variant's design-window targetons, so a variant spanning a screened and an unscreened targeton reads `not_screened|screened`; split on `|` rather than substring-matching, since `not_screened` contains `screened`."),
      ("assay_testable_now", "bool", "True only if at least one containing targeton has completed screening."),
     ]),

    ("Evidence counts",
     "Recurrence is evidence of a hotspot only when observations are "
     "independent, so the counts are deliberately separate. Use "
     "`n_independent_probands` for any hotspot inference; it counts ONE PER "
     "FAMILY. Check `proband_count_is_upper_bound` before quoting it.",
     [
      ("n_observations", "int", "Rows in the observation ledger for this variant. Inflated by re-reporting and by multiple submitters depositing the same case."),
      ("n_individuals", "int", "Distinct named individuals across all sources."),
      ("n_independent_probands", "int", "**One per family.** The hotspot-safe count."),
      ("n_families", "int", "Distinct families (equal to the above by construction; both emitted for clarity)."),
      ("n_siblings_collapsed", "int", "`n_individuals - n_families`: how much family grouping reduced the naive count."),
      ("is_recurrent", "bool", "True if seen in more than one independent proband."),
      ("proband_count_is_upper_bound", "bool", "True when independence is not established -- because a source supplies no relatedness field (GeneDx) or because the variant appears in more than one source and individuals cannot be linked across them."),
      ("proband_independence_caveat", "string", "Why the count is an upper bound."),
      ("kinship_needs_adjudication", "bool", "True if any contributing record states kinship only in free text. Such records are flagged for a human decision, not parsed by regex."),
      ("n_observations_without_individual", "int", "Contributing rows with no individual attached (LOVD submitter classification records)."),
     ]),

    ("Provenance",
     "Every claim traces back to a source record. PMIDs are extracted from "
     "LOVD's structured citation markup.",
     [
      ("n_sources", "int", "Number of distinct sources reporting this variant."),
      ("sources", "list", "Source identifiers, e.g. `lovd`, `genedx`, `lovd|genedx`."),
      ("source_records", "list", "Source-native record identifiers."),
      ("individual_keys", "list", ("Individual identifiers, namespaced by source: `LOVD:<8-digit LOVD id>`, "
      "`GDX:<6-digit anonymised sample id>`, `DDD:<DECIPHER patient id>`, "
      "`PUB:<CitationCompact>:<patient label>`.")),
      ("lovd_dbids", "list", "LOVD DB-IDs, e.g. `SMC1A_000034`. Note the prefix names the hosting gene database, not the gene affected."),
      ("pmids", "list", "PubMed IDs of the primary publications, `;`-separated. Parsed from the source's citation markup, not the raw text."),
      ("citations", "list", "Author/year labels accompanying each PMID."),
      ("dois", "list", "DOIs where the source records one."),
      ("omim_allelic_variants", "list", "OMIM allelic-variant ids (OMIM 300040 is the SMC1A gene entry) -- an independent attribution source."),
      ("reference_raw", "string", "Citation field(s) verbatim, retained so the source wording is recoverable."),
      ("clinvar_ids", "list", "ClinVar identifiers. Cross-reference only -- ClinVar is not used for disease attribution."),
      ("dbsnp_id", "list", "dbSNP rsIDs."),
      ("published_as", "list", "The description as it appeared in the original publication, verbatim."),
     ]),

    ("Pathogenicity classification",
     "Consensus across sources. A variant reported both pathogenic and benign "
     "is `conflicting` -- a finding, not something to average away.",
     [
      ("pathogenicity_consensus", "string", "`P`, `LP`, `VUS`, `LB`, `B`, `conflicting`, or `not_classified`."),
      ("pathogenicity_all_reported", "list", "Every classification as reported, verbatim."),
      ("pathogenicity_note", "string", "Explanation where sources differ."),
      ("classification_methods", "list", "Classification method where the source records one."),
     ]),

    ("Disease attribution",
     "The binding constraint of this curation. Resolved only to `CdLS` or "
     "`DEE85`; anything else routes the variant to an `Exclude_*` group. "
     "ClinVar disease labels are deliberately not used.",
     [
      ("disease_attribution", "string", "`CdLS`, `DEE85`, `conflicting` (sources disagree), `unresolved` (a condition was recorded but does not map to either), or `none` (no condition recorded)."),
      ("disease_attribution_note", "string", "Basis for the value, including any additional codes carried."),
      ("disease_attribution_declined", "bool",
       ("True where a source explicitly considered a CdLS/DEE85 attribution for "
        "this variant's patient and recorded that it did not apply, as distinct "
        "from never addressing it. Read from the source's own per-case field, "
        "never inferred. Such a variant is excluded as "
        "`disease_attribution_declined_by_source` rather than "
        "`no_disease_attribution_yet`, because the latter group is upgradeable "
        "by further clinical detail and this one is not.")),
      ("disease_attribution_declined_by", "list",
       "Sources that explicitly declined the attribution."),
      ("disease_all_reported", "list", "Every disease code as reported, verbatim."),
      ("has_clinical_detail", "bool", "True if at least one observation carries phenotype detail, a recorded diagnosis, or an episignature result -- i.e. something beyond a bare disease label. Inclusion criterion 4."),
     ]),

    ("Individual-level detail",
     "Sex is recorded because it is mechanistically load-bearing for SMC1A: "
     "DEE85 PTVs are almost exclusively female and hemizygous male PTVs are "
     "thought largely non-viable.",
     [
      ("sex_M", "int", "Independent probands recorded male."),
      ("sex_F", "int", "Independent probands recorded female."),
      ("sex_unknown", "int", "Probands with no sex recorded."),
      ("inheritance_summary", "string", "Harmonised inheritance with counts, e.g. `de_novo:2|inherited_maternal:1`."),
      ("n_de_novo", "int", "Observations recorded as de novo."),
      ("phenotype_detail", "string", "Concatenated phenotype fields from all sources (LOVD Phenotypes section, incl. definite/initial diagnosis)."),
      ("hpo_ids", "list", "HPO term ids, carried verbatim. No syndrome inference is drawn from them."),
      ("individual_remarks", "string", "Free-text individual remarks from the source."),
     ]),

    ("Episignature and functional evidence",
     "Orthogonal evidence. Populated at the literature step; empty until then.",
     [
      ("episignature_available", "bool", "True if any source reports an episignature result."),
      ("episignature_result", "list", "Result as reported, e.g. CdLS-episignature positive/negative."),
      ("episignature_source", "list", "Source of the episignature result."),
      ("functional_evidence", "list", "Functional assay evidence reported by a source."),
     ]),

    ("Quality flags",
     "Each flag marks a specific hazard identified during curation. None is "
     "cosmetic; see DECISIONS_LOG for the incident behind each.",
     [
      ("not_smc1a_only", "bool", "Event spans more than one gene (assessed on genomic span, not on the hosting database), so a phenotype cannot be attributed to SMC1A alone."),
      ("alt_isoform_numbering", "bool", "A source described this variant on NM_001281463.1, a different SMC1A isoform whose c. numbering is offset by +66 nt / +22 codons. Requires conversion, not face-value use."),
      ("possible_cross_source_duplicate", "bool", "Reported by more than one source with no way to link individuals across them, so the same person may be counted twice."),
      ("resolution_status", "string", "`resolved`, `resolved_with_flag`, `quarantined`, or `direct_coordinates` (source supplied verified GRCh38 coordinates)."),
      ("guard_failed", "string", "Which coordinate-resolution guard failed (G1-G6). See METHODS 3.2."),
      ("resolution_notes", "string", "Detail of any resolution problem."),
      ("notes", "string", "Free-text variant remarks from sources."),
     ]),

    ("Curation outcome",
     "The assignment and its justification. `curation_group` is exactly one of "
     "the six specified groups; `exclusion_reason` carries the specifics.",
     [
      ("curation_group", "string", "`CdLS_pathogenic`, `DEE85_pathogenic`, `Exclude_VUS`, `Exclude_Benign`, `Exclude_disease_uncertain`, or `Exclude_pending_classification`."),
      ("curation_group_reason", "string", "Prose justification for the assignment."),
      ("include_high_confidence", "bool", "True for rows in smc1a_variants_curated.tsv, i.e. meeting all five inclusion criteria."),
      ("exclusion_reason", "string", "`no_disease_attribution_yet` (outstanding, expected to be upgradable), `disease_recorded_but_unresolvable` (irreducibly ambiguous), `conflicting_disease_attribution`, `conflicting_classification`, `vus`, `benign`, `insufficient_clinical_detail`, `not_smc1a_only`, `unresolved_coordinates`, `not_classified`."),
      ("curation_date", "date", "Date this row was generated."),
      ("pipeline_version", "string", "Pipeline version, cross-referenced to docs/CHANGELOG.md."),
     ]),
    ("Source-reported annotations carried through to variant level",
     "Fields contributed by individual sources and retained on the collapsed variant row. Each records what a source stated; none is inferred.",
     [
         ("cdls_differential_reported", "str", "Source explicitly considered CdLS in the differential."),
         ("dual_molecular_diagnosis", "str", "Individual carries a second molecular diagnosis."),
         ("zygosity_reported", "str", "Zygosity or allelic state as stated by the source."),
         ("cdls_presentation_absent_in_a_source", "str", "A source states CdLS features were absent."),
     ]),
    ("Residue-level mutational hotspots",
     "The same amino acid residue altered by different variants in different "
     "individuals -- neither a variant-identity duplicate nor a patient-identity "
     "duplicate. Weak independent evidence that the residue matters to protein "
     "function. Computed by systematic scan, not taken from any source's own "
     "cross-referencing.",
     [
         ("residue_hotspot", "bool", "This residue is altered by more than one distinct variant in the registry."),
         ("residue_hotspot_residue", "str", "The residue, e.g. Arg1049."),
         ("residue_hotspot_n_variants", "int", "Number of distinct variants altering this residue."),
         ("residue_hotspot_groups", "str", "Analysis arms the co-hitting variants fall into; a residue spanning both arms is directly relevant to the CdLS/DEE85 mechanism contrast."),
         ("residue_hotspot_partners", "str", "The other variant(s) altering this residue."),
     ]),
    ("Literature-reported domain hotspot",
     "Whether a DEE85_pathogenic missense variant sits in the N-/C-terminal "
     "ATPase head domain that the literature reports as where non-loss-of-"
     "function (missense/in-frame) DEE85 variants specifically cluster "
     "(Baranano et al. 2022; Bozarth et al. 2023; Di Nardo et al. 2026). "
     "Unlike every other column here these are **not** written by step 05: "
     "they are added afterwards by "
     "`Code/investigations/annotate_dee85_missense_domain_hotspot.py`, which "
     "must be re-run against both processed tables after any registry rebuild "
     "or they are silently lost. Only the 11 DEE85_pathogenic missense "
     "variants are evaluated; every other row is False with the domain and "
     "source columns blank, since the claim was not assessed for them.",
     [
         ("literature_domain_hotspot", "bool", "A specific paper explicitly places this variant in the N- or C-terminal ATPase head domain. Variants the literature discusses but places in the coiled-coil arm, outside the reported cluster, are False."),
         ("literature_domain_hotspot_domain", "str", "The domain as the citing paper describes it, e.g. `N-terminal ATPase head (aa 4-148)`."),
         ("literature_domain_hotspot_source", "str", "The citation placing it there, including the basis (cohort observation or molecular-dynamics modelling)."),
     ]),
    ("Coordinate provenance",
     "Whether the cDNA description was read from the source or derived from a "
     "protein change under the unique-single-nucleotide-route rule. Kept explicit "
     "so a derived coordinate is never mistaken for an extracted one.",
     [
         ("hgvs_c_derived_from_protein", "bool", "True if at least one contributing observation's cDNA was back-derived from a protein change (unique route only; ambiguous routes and frameshifts are refused)."),
     ]),
    ("Contradicted CdLS attribution",
     "A CdLS label can arrive either from a publication that clinically "
     "diagnosed the patient or from a bare gene-level tag on a database record. "
     "Where a source states CdLS was NOT in the differential and no publication "
     "supports the label, the attribution is not carried.",
     [
         ("cdls_attribution_contradicted", "bool", "True if the CdLS consensus is contradicted by a source and no publication-derived observation supports it; such variants are excluded from the CdLS arm."),
         ("cdls_attribution_note", "str", "Why the contradiction was or was not acted on, naming the contradicting source and any publication that independently supports CdLS."),
     ]),
    ("Protein description provenance",
     "How the protein description was chosen where sources disagree or write it "
     "malformed, and why a column is empty when it is.",
     [
         ("hgvs_p_note", "str", "Records a normalisation (a missing `p.` prefix supplied, naming the source) or why no usable protein description exists. Empty when the description was taken unchanged."),
     ]),
]

# ==========================================================================
# Observation ledger: smc1a_observations.tsv
# ==========================================================================

OBSERVATIONS = [

    ("Keys",
     "One row per (variant x source x individual/report). This is the "
     "lossless evidence ledger from which both variant tables are derived; "
     "use it to audit any variant-level value.",
     [
      ("variant_key_resolved", "string", "Foreign key to `variant_key` in the variant tables."),
      ("source", "string", "`LOVD` or `GeneDx_180k_DNM`."),
      ("source_detail", "string", "Full source description."),
      ("source_record_id", "string", "Source-native record identifier."),
      ("individual_key", "string", ("Individual identifier, namespaced by source -- `LOVD:<8-digit LOVD id>`, "
      "`GDX:<6-digit anonymised sample id>`, `DDD:<DECIPHER patient id>`, "
      "`PUB:<CitationCompact>:<patient label>`. Empty for submitter "
      "classification records, which name no individual.")),
     ]),

    ("Family and independence",
     "Kinship is taken only from structured fields. Free-text kinship is "
     "flagged for human adjudication, never parsed by regex, because most "
     "remarks describe a de novo singleton (implying independence) and only "
     "some imply an affected relative.",
     [
      ("family_key", "string", "Family identifier. Assigned by union-find over LOVD `panelid` links and the `<PMID>.<family><P|S>` submitter-label convention."),
      ("family_key_basis", "string", "Which rule assigned the family, or `singleton`."),
      ("is_independent_proband", "string", "`yes`, `assumed` (no relatedness field supplied by the source), or `no_individual`."),
      ("proband_independence_basis", "string", "Why independence is or is not established."),
      ("kinship_needs_adjudication", "bool", "True where kinship is stated only in free text, or the entry stands for more than one person."),
      ("kinship_note", "string", "The free-text kinship wording, verbatim, for adjudication."),
      ("lovd_panel_id", "string", "LOVD panel this individual belongs to."),
      ("lovd_panel_size", "int", "Number of people the LOVD entry represents. >1 means one row stands for several individuals."),
      ("individual_label_raw", "string", "Submitter's own patient label, e.g. `16604071.7P` (PMID 16604071, family 7, proband)."),
     ]),

    ("Variant description as submitted",
     "Retained verbatim. Nothing here is overwritten by the recalibrated "
     "values, so the original claim is always recoverable.",
     [
      ("hgvs_c_input", "string", "c. description as submitted by the source."),
      ("hgvs_c_input_transcript", "string", "Transcript the source described it on."),
      ("hgvs_p_input", "string", "Protein description as submitted."),
      ("hgvs_r_input", "string", "RNA description as submitted."),
      ("exon_input", "string", "Exon as stated by the source."),
      ("published_as", "string", "Description as it appeared in the original publication."),
      ("published_as_transcript", "string", "Transcript accession cited by the original publication."),
      ("conversion_audit", "string", "Result of auditing the source's transcript conversion against the derived isoform offset."),
      ("alt_isoform_numbering", "bool", "Source used NM_001281463.1 (+66 nt / +22 codons offset)."),
      ("hgvs_g_hg19_input", "string", "LOVD genomic description on hg19."),
      ("hgvs_g_hg38_input", "string", "LOVD genomic description on hg38, used as an independent cross-check."),
      ("lovd_pos_g_start_hg19", "int", "LOVD hg19 start."),
      ("lovd_pos_g_end_hg19", "int", "LOVD hg19 end."),
      ("genomic_span_bp", "int", "Event span in bp; >1000 cannot be SMC1A-only at this locus."),
      ("lovd_variant_type", "string", "LOVD variant type."),
      ("iscn", "string", "ISCN description for cytogenetic events."),
      ("not_smc1a_only", "bool", "Event spans more than one gene."),
     ]),

    ("Resolved coordinates",
     "Plus-strand GRCh38, recalibrated onto the MANE transcript.",
     [
      ("variant_key", "string", "Canonical plus-strand key, where the source supplied coordinates directly."),
      ("chrom", "string", "Chromosome."),
      ("pos", "int", "Plus-strand GRCh38 position."),
      ("ref", "string", "Plus-strand reference allele."),
      ("alt", "string", "Plus-strand alternate allele."),
      ("ref_allele_verified", "bool", "Reference allele confirmed against the cached GRCh38 sequence."),
      ("ref_allele_in_grch38", "string", "What the reference genome actually has at that position."),
      ("consequence_class", "string", "Harmonised consequence class on MANE."),
      ("consequence_terms", "string", "Raw VEP SO terms."),
      ("protein_position", "int", "Codon number."),
      ("exon", "int", "Exon number (1-25), empty if intronic."),
      ("intron", "int", "Intron number, empty if exonic."),
     ]),

    ("SGE assay context",
     "",
     [
      ("targetons_design", "list", "Targeton(s) whose design window contains the variant."),
      ("targetons_exon_core", "list", "Targeton(s) whose exonic core (r2) contains it."),
      ("targetons_amplicon", "list", "Targeton(s) whose amplicon contains it."),
      ("in_design_window", "bool", "True if in at least one design window."),
      ("exon_source_supplied", "string", "Exon as supplied by the source, before recomputation."),
      ("targeton_source_supplied", "string", "Targeton as supplied by the source, before recomputation."),
     ]),

    ("Classification as reported",
     "",
     [
      ("pathogenicity_reported", "string", "Mapped to the 5-tier ACMG vocabulary."),
      ("pathogenicity_reported_raw", "string", "Classification exactly as the source stated it."),
      ("classification_method_raw", "string", "Method as stated."),
      ("classification_note", "string", "Interpretation note, e.g. what a laboratory diagnostic status means."),
      ("lovd_effect_code", "string", "LOVD two-digit effect code."),
      ("lovd_effect_decoded", "string", "Decoded as reported/concluded effect on function."),
      ("eligible_for_curated_set", "bool", "Source-specific inclusion gate. For GeneDx, true only for `diagnostic_status = positive`."),
      ("genedx_diagnostic_positive", "bool", "GeneDx laboratory judged the SMC1A variant to be the diagnostic finding."),
     ]),

    ("Disease and phenotype",
     "",
     [
      ("disease_reported_raw", "string", "Disease code exactly as the source stated it."),
      ("disease_reported_name", "string", "Full disease name from the source."),
      ("disease_harmonised", "string", "Mapped to `CdLS`, `DEE85`, `DEE_unspecified`, `epilepsy_only`, `unspecified_NDD`, `unknown`, or `other`."),
      ("attribution_basis", "string", "How the attribution was arrived at."),
      ("phenotype_detail", "string", "Phenotype fields from the source, incl. LOVD definite/initial diagnosis."),
      ("hpo_ids", "list", "HPO term ids, verbatim."),
      ("hpo_names", "list", "HPO term labels, verbatim."),
      ("n_hpo_terms", "int", "Number of distinct HPO terms."),
      ("methylation_raw", "string", "LOVD methylation field."),
      ("individual_remarks", "string", "Free-text individual remarks."),
      ("variant_remarks", "string", "Free-text variant remarks."),
     ]),

    ("Individual detail",
     "",
     [
      ("sex", "string", "`M`, `F`, or empty. Harmonised."),
      ("sex_raw", "string", "Source value verbatim (LOVD labels this field 'Gender')."),
      ("country", "string", "Geographic origin."),
      ("population", "string", "Population origin."),
      ("consanguinity", "string", "Consanguinity as stated."),
      ("age_of_death", "string", "Age at death where recorded."),
      ("data_availability", "string", "LOVD Data_av field, e.g. a DECIPHER link."),
     ]),

    ("Inheritance",
     "LOVD records inheritance in two fields that frequently disagree; "
     "disagreements are surfaced as `conflicting` rather than resolved by "
     "precedence.",
     [
      ("inheritance", "string", "`de_novo`, `inherited_maternal`, `inherited_paternal`, `inherited_unspecified_parent`, `somatic`, `conflicting`, `unknown`, `not_applicable`."),
      ("inheritance_detail", "string", "Basis, including any conflict between source fields."),
      ("genetic_origin_raw", "string", "LOVD Genetic_origin verbatim."),
      ("allele_raw", "string", "LOVD allele assignment, decoded from its numeric code."),
      ("allele_code_raw", "string", "The raw LOVD numeric allele code."),
      ("segregation_raw", "string", "Segregation as stated."),
     ]),

    ("References",
     "",
     [
      ("pmids", "list", "PubMed IDs extracted from the source's citation markup."),
      ("citations", "list", "Author/year labels accompanying each PMID."),
      ("dois", "list", "DOIs."),
      ("omim_allelic_variants", "list", "OMIM allelic-variant ids (OMIM 300040 is the SMC1A gene entry)."),
      ("reference_is_unpublished", "bool", "Source cites unpublished data or a personal communication."),
      ("reference_raw", "string", "Citation field verbatim."),
      ("clinvar_id", "string", "ClinVar identifier."),
      ("dbsnp_id", "string", "dbSNP rsID."),
      ("lovd_dbid", "string", "LOVD DB-ID."),
      ("lovd_individual_id", "string", "LOVD internal individual id."),
      ("lovd_owner", "string", "LOVD submitter id."),
      ("frequency_raw", "string", "Allele frequency as stated by the source."),
      ("screening_technique", "string", "Assay used to detect the variant."),
      ("screening_template", "string", "Template used (DNA/RNA)."),
      ("variant_class_supplied", "string", "Variant class as supplied by the source."),
      ("consequence_supplied", "string", "Consequence as supplied, before recomputation."),
      ("cds_change_supplied", "string", "c. change as supplied, before recomputation."),
      ("protein_change_supplied", "string", "Protein change as supplied, before recomputation."),
      ("net_length_change", "int", "Net length change in bp."),
     ]),

    ("Curation outcome",
     "",
     [
      ("curation_group", "string", "Group assigned to the parent variant."),
      ("exclusion_reason", "string", "Why excluded, if it was."),
      ("notes", "string", "Free text."),
     ]),
    ("Source-reported annotations, provenance and re-report tracking",
     "Observation-level fields recording what each source stated verbatim, plus the re-report resolution that prevents one individual being counted twice.",
     [
         ("xci_raw", "str", "X-inactivation status exactly as printed by the source, before harmonisation."),
         ("xci_skewing", "str", "Harmonised X-inactivation skewing category."),
         ("xci_ratio", "str", "Numeric X-inactivation ratio where the source reports one."),
         ("episignature_source", "str", "Source of the DNA methylation episignature result."),
         ("episignature_result", "str", "Episignature classification as reported."),
         ("functional_evidence", "str", "Functional or experimental evidence reported, verbatim."),
         ("protein_domain_reported", "str", "Protein domain assignment as stated by the source."),
         ("zygosity_reported", "str", "Zygosity or allelic state as stated by the source."),
         ("additional_notes", "str", "Free-text notes that do not fit an existing column."),
         ("cdls_differential_reported", "str", "Source explicitly considered CdLS in the differential."),
         ("dual_molecular_diagnosis", "str", "Individual carries a second molecular diagnosis."),
         ("competing_attribution", "str", "A second, conflicting disease attribution for this variant."),
         ("is_rereport", "str", "Observation re-reports a previously published individual."),
         ("rereport_of_pmid", "str", "PMID of the original report."),
         ("disease_attribution_suppressed", "str",
          ("Set where this observation's disease claim was deliberately not "
           "counted, with the reason. Applies to a database record that is an "
           "identified deposit of a publication ingested separately: the "
           "publication governs the disease attribution and the deposit adds "
           "no independent claim, so `disease_harmonised` is cleared while "
           "classification, sex, inheritance and phenotype are retained. "
           "Empty for every other observation.")),
         ("rereport_of_citation", "str", "Citation of the original report."),
         ("rereport_of_patient", "str", "Patient identifier in the original report."),
         ("primary_patient_id", "str", "Patient identifier resolved to the primary report."),
         ("primary_citation", "str", "Citation resolved to the primary report."),
         ("is_overlay_source", "str", "Observation comes from an overlay-only source."),
     ]),
    ("Coordinate provenance",
     "Whether the cDNA description was read from the source or derived from a "
     "protein change under the unique-single-nucleotide-route rule. Kept explicit "
     "so a derived coordinate is never mistaken for an extracted one.",
     [
         ("hgvs_c_derived_from_protein", "bool", "True if at least one contributing observation's cDNA was back-derived from a protein change (unique route only; ambiguous routes and frameshifts are refused)."),
     ]),
    ("Cohort observation fields",
     "Fields emitted by the unpublished cohort ingests (GeneDx, DDD/DECIPHER), "
     "which supply coordinates and structured phenotype directly rather than an "
     "HGVS description and free text.",
     [
         ("hgvs_c_mane", "str", "MANE cDNA description carried on the observation where the source supplies coordinates rather than an HGVS string."),
         ("hgvs_p_mane", "str", "MANE protein description carried on the observation."),
         ("pmid", "str", "PubMed identifier of the source, where published."),
         ("citation", "str", "Citation string for the source."),
         ("coordinate_basis", "str", "How the coordinate was obtained -- supplied directly by the source, or derived."),
         ("source_transcript", "str", "Transcript the source reported against."),
         ("source_transcript_is_mane", "str", "Whether the source transcript is the MANE Select transcript."),
         ("source_consequence", "str", "Consequence term as stated by the source, before recomputation."),
         ("vep_exon", "str", "Exon number returned by VEP for this observation."),
         ("vep_intron", "str", "Intron number returned by VEP for this observation."),
         ("has_clinical_detail", "str", "Whether the observation carries phenotype detail beyond a bare disease label."),
         ("associated_condition_raw", "str", "Disease label as entered in the cohort file, verbatim."),
         ("associated_condition_notes", "str", "Qualifier accompanying the disease label (e.g. a Rett-like presentation)."),
         ("patient_id_source", "str", "Patient identifier as used by the source cohort."),
         ("decipher_variant_id", "str", "DECIPHER internal variant identifier."),
     ]),
    ("Per-patient attributes and kinship",
     "Sex and family membership are properties of the PATIENT, not the variant. "
     "A per-variant row carries the union across that variant's patients, with "
     "the per-patient assignment preserved so the union can be decomposed.",
     [
         ("sex_by_patient", "str", "Per-patient sex assignment as stated by the source, e.g. \"7P=female; 7S=female\"."),
         ("family_basis", "str", "Why this individual was assigned to its family: an explicit statement in the source, or the default that no relationship was stated."),
     ]),
]

TABLES = {
    "smc1a_variants_all.tsv": VARIANTS,
    "smc1a_variants_curated.tsv": VARIANTS,
    "smc1a_observations.tsv": OBSERVATIONS,
}


def order_columns(present, schema):
    """Return `present` ordered per `schema`, with any unlisted column
    appended (never dropped), plus the list of unlisted columns."""
    wanted = [c for _t, _n, cols in schema for c, _ty, _d in cols]
    ordered = [c for c in wanted if c in present]
    extra = [c for c in present if c not in wanted]
    return ordered + extra, extra
