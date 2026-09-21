#!/usr/bin/env Rscript
#
# cdls_vs_dee85_missense_predictors.R
#
# Follow-up to the assay's own CdLS-vs-DEE85 missense split (06_join_assay.py's
# default missense-only comparison; OR=0.0163, p=3.83e-05 on the promoted
# curated set): for the same curated missense variants, are in-silico
# predictors (SIFT, PolyPhen, CADD, REVEL) tracking the same CdLS-vs-DEE85
# split the assay's own depletion data shows, or missing it?
#
# Joins two already-computed outputs rather than raw data:
#   - the curated-set-vs-assay join's missense analysis table
#     (06_join_assay.py --class_filter missense; one row per curated missense
#     variant, with curation_group and anchor_tier)
#   - the in-silico concordance merge (insilico_concordance.R; one row per
#     assayed oligo, with predictor scores)
# joined on oligo_name + targeton.
#
# ---------------------------------------------------------------------
# INPUTS
# ---------------------------------------------------------------------
# --assay_join   assay_join_analysis.tsv from 06_join_assay.py, missense
#                class_filter. Needs: oligo_name, matched_targeton,
#                curation_group, hgvs_p_mane, anchor_tier,
#                pos_adj_log2FoldChange_raw.
# --insilico     merged_variant_data.tsv from insilico_concordance.R. Needs:
#                oligo_name, Targeton_ID, sift_damaging, polyphen_score,
#                CADD_num, REVEL_num.
# --outdir       Output directory (default: current directory)
#
# ---------------------------------------------------------------------
# OUTPUTS
# ---------------------------------------------------------------------
#   cdls_dee85_missense_predictors.tsv   one row per matched variant: group,
#                                        depletion status, all four predictor
#                                        scores
#   cdls_dee85_missense_summary.tsv      per-group, per-predictor: n, median,
#                                        Mann-Whitney U p-value (CdLS vs DEE85)
#   cdls_dee85_missense_boxplots.png     one panel per predictor, boxplot by
#                                        curation_group
#
# ---------------------------------------------------------------------
# USAGE
# ---------------------------------------------------------------------
#   Rscript cdls_vs_dee85_missense_predictors.R \
#       --assay_join  data/assay_join/assay_join_analysis.tsv \
#       --insilico    insilico_concordance/merged_variant_data.tsv \
#       --outdir      results/
#
# ---------------------------------------------------------------------

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) stop("Package 'optparse' is required.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
  library(optparse)
  library(data.table)
})

have_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--assay_join", type = "character"),
  make_option("--insilico", type = "character"),
  make_option("--outdir", type = "character", default = ".")
)))

if (is.null(opt$assay_join) || is.null(opt$insilico)) {
  stop("Both --assay_join and --insilico are required. See header comment for usage.")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

assay <- fread(opt$assay_join)
insilico <- fread(opt$insilico)

required_assay <- c("oligo_name", "matched_targeton", "curation_group", "hgvs_p_mane",
                     "anchor_tier", "pos_adj_log2FoldChange_raw")
missing_assay <- setdiff(required_assay, names(assay))
if (length(missing_assay) > 0) stop(sprintf("%s missing column(s): %s", opt$assay_join, paste(missing_assay, collapse = ", ")))

required_insilico <- c("oligo_name", "Targeton_ID", "sift_damaging", "polyphen_score", "CADD_num", "REVEL_num")
missing_insilico <- setdiff(required_insilico, names(insilico))
if (length(missing_insilico) > 0) stop(sprintf("%s missing column(s): %s", opt$insilico, paste(missing_insilico, collapse = ", ")))

# insilico's Targeton_ID carries a "_exonN" suffix (e.g. "APDY_exon2") that
# assay_join's matched_targeton (e.g. "APDY") doesn't -- strip it before
# joining.
insilico$targeton_code <- sub("_exon[0-9]+$", "", insilico$Targeton_ID)

merged <- merge(
  assay[, c("oligo_name", "matched_targeton", "curation_group", "hgvs_p_mane",
            "anchor_tier", "pos_adj_log2FoldChange_raw")],
  insilico[, c("oligo_name", "targeton_code", "sift_damaging", "polyphen_score", "CADD_num", "REVEL_num")],
  by.x = c("oligo_name", "matched_targeton"), by.y = c("oligo_name", "targeton_code"),
  all.x = TRUE
)

n_unmatched <- sum(is.na(merged$sift_damaging) & is.na(merged$polyphen_score) &
                    is.na(merged$CADD_num) & is.na(merged$REVEL_num))
if (n_unmatched > 0) {
  cat(sprintf("NOTE: %d / %d variant(s) had no matching row in --insilico (all predictor columns NA).\n",
              n_unmatched, nrow(merged)))
}

setnames(merged, c("sift_damaging", "polyphen_score", "CADD_num", "REVEL_num"),
         c("SIFT_damaging", "PolyPhen", "CADD_PHRED", "REVEL"))

out_path <- file.path(opt$outdir, "cdls_dee85_missense_predictors.tsv")
fwrite(merged, out_path, sep = "\t")
cat(sprintf("Wrote %s (%d rows)\n", out_path, nrow(merged)))

cat("\n=== Depletion status by group ===\n")
print(table(merged$curation_group, merged$anchor_tier))

predictors <- c("SIFT_damaging", "PolyPhen", "CADD_PHRED", "REVEL")
summary_rows <- list()
for (p in predictors) {
  cdls_vals <- merged[[p]][merged$curation_group == "CdLS_pathogenic"]
  dee85_vals <- merged[[p]][merged$curation_group == "DEE85_pathogenic"]
  cdls_vals <- cdls_vals[!is.na(cdls_vals)]
  dee85_vals <- dee85_vals[!is.na(dee85_vals)]
  wt <- tryCatch(wilcox.test(dee85_vals, cdls_vals), error = function(e) NULL)
  summary_rows[[p]] <- data.frame(
    predictor = p,
    n_CdLS = length(cdls_vals), median_CdLS = if (length(cdls_vals) > 0) median(cdls_vals) else NA,
    n_DEE85 = length(dee85_vals), median_DEE85 = if (length(dee85_vals) > 0) median(dee85_vals) else NA,
    mann_whitney_p = if (!is.null(wt)) wt$p.value else NA
  )
}
summary_df <- do.call(rbind, summary_rows)
summary_path <- file.path(opt$outdir, "cdls_dee85_missense_summary.tsv")
fwrite(summary_df, summary_path, sep = "\t")
cat("\n=== Predictor score by group (CdLS vs DEE85, missense only) ===\n")
print(summary_df, row.names = FALSE)
cat(sprintf("\nWrote %s\n", summary_path))

if (have_ggplot2) {
  suppressPackageStartupMessages(library(ggplot2))
  long <- do.call(rbind, lapply(predictors, function(p) {
    data.frame(predictor = p, curation_group = merged$curation_group, score = merged[[p]])
  }))
  long <- long[!is.na(long$score) & !is.na(long$curation_group), ]

  p <- ggplot(long, aes(x = curation_group, y = score, fill = curation_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
    facet_wrap(~predictor, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(x = NULL, y = "predictor score (higher = more damaging)",
         title = "Missense variants: predictor scores by curated disease group")

  plot_path <- file.path(opt$outdir, "cdls_dee85_missense_boxplots.png")
  ggsave(plot_path, p, width = 9, height = 7, dpi = 300)
  cat(sprintf("Wrote %s\n", plot_path))
} else {
  cat("NOTE: ggplot2 not installed, skipping boxplot.\n")
}
