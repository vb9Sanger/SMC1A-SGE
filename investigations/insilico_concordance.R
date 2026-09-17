#!/usr/bin/env Rscript
#
# insilico_concordance.R
#
# Tests concordance between in-silico variant-effect predictors (SIFT,
# PolyPhen, CADD, REVEL, SpliceAI) and the SMC1A SGE assay's functional
# calls, gene-wide, stratified by consequence class (predictors don't all
# apply to all classes -- SIFT/PolyPhen/REVEL are missense-only; CADD and
# SpliceAI are genome-wide).
#
# For each (consequence class, applicable predictor) pair with enough data:
#   1. Spearman correlation: predictor score vs z_noimpact (continuous).
#      Spearman rather than Pearson because these predictor scores are
#      typically skewed/bounded, not normally distributed.
#   2. Kruskal-Wallis test: predictor score across anchor_tier groups
#      (no impact / weakly depleting / strongly depleting). "enriched"
#      calls are excluded throughout -- none of these predictors are
#      designed to predict gain-of-function-like behaviour, so including
#      them would just add noise to a comparison they were never meant to
#      inform.
#   3. AUC: treating anchor_tier ("depleted" = weakly+strongly depleting,
#      vs "no impact") as ground truth and the predictor score as the
#      classifier being evaluated. This is the standard "how well does this
#      predictor recover what my assay already tells me" framing. Uses
#      pROC if installed; falls back to a dependency-free Mann-Whitney-U-
#      based AUC calculation otherwise (numerically identical result).
#      Alongside AUC, also reports the data-derived (Youden's J-optimal)
#      threshold and the sensitivity/specificity AT that threshold -- AUC
#      is threshold-free and rigorous, but "sensitivity at the best cutoff"
#      answers the more directly interpretable question: what fraction of
#      our assay-depleted variants would this predictor actually catch.
#
# z_noimpact = (LFC - anchor_mu_noimpact) / anchor_sd_noimpact
#   i.e. how many no-impact-control SDs a variant's LFC sits from its own
#   targeton's no-impact anchor. More negative = more depletion-like.
#
# ---------------------------------------------------------------------
# INPUTS
# ---------------------------------------------------------------------
# --anchor_dir   Directory of per-targeton anchor-tier/DESeq2 output TSVs
#                (gaussian_shrinkage_classifier.R output). Needs columns:
#                oligo_name, consequence, anchor_tier, pos_adj_log2FoldChange_raw,
#                anchor_mu_noimpact, anchor_sd_noimpact.
#                Targeton_ID is inferred from each filename (the part
#                before "_all_deseq2_results_condition_").
#
# --meta_dir     Directory of per-targeton meta_consequences TSVs (VEP-style
#                annotation with predictor scores). Needs columns:
#                unique_oligo_name, SIFT, PolyPhen, CADD_PHRED, REVEL,
#                SpliceAI_pred.
#                Targeton_ID is inferred from each filename (the part
#                before "_chrX_"/"_chr<N>_").
#
#                The two directories are joined per targeton by matching
#                inferred Targeton_ID, then per-variant by stripping the
#                trailing "_<region_start>_<region_end>" suffix from
#                unique_oligo_name so it matches anchor_dir's oligo_name
#                exactly (verified against real data: 100% match rate).
#
# --outdir       Output directory (default: current directory)
#
# ---------------------------------------------------------------------
# PREDICTOR FORMAT NOTES (handled automatically, documented here for
# transparency)
# ---------------------------------------------------------------------
# - SIFT: e.g. "tolerated(0.83)", "deleterious(0)", "-" for missing.
#   Lower raw score = more damaging (SIFT convention). We report
#   sift_damaging = 1 - raw_score, so higher = more damaging, consistent
#   with every other predictor here.
# - PolyPhen: e.g. "probably_damaging(0.989)", "benign(0.009)", "-".
#   Already higher = more damaging.
# - CADD_PHRED: numeric, "-" for missing. Higher = more damaging.
# - REVEL: numeric 0-1, "-" for missing. Higher = more damaging.
# - SpliceAI_pred: pipe-delimited, e.g. "SMC1A|0.00|0.00|0.01|0.00|-29|-36|44|36"
#   = GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL (the four DS_*
#   fields are delta scores for acceptor-gain/acceptor-loss/donor-gain/
#   donor-loss; the four DP_* fields are delta positions, not used here).
#   We take the max of the four DS_* delta scores as a single
#   "spliceai_max" summary, standard practice for this kind of analysis.
#   "-" for missing.
#
# ---------------------------------------------------------------------
# OUTPUTS
# ---------------------------------------------------------------------
#   merged_variant_data.tsv:
#     Full merged table, one row per variant, ALL consequence classes and
#     enriched calls INCLUDED (this file is not filtered for the tier
#     analyses below) -- use this to look up any specific variant's
#     predictor scores, e.g. checking the SpliceAI score of a known
#     outlier synonymous variant.
#   insilico_concordance_tables.tsv (one file, four sections):
#     Table 0: tier definitions and per-class notes (read this first)
#     Table 1: Spearman correlation, predictor vs z_noimpact
#     Table 2: Kruskal-Wallis test, predictor across anchor_tier
#     Table 3: AUC, predictor discriminating depleted vs no-impact
#   auc_summary_chart.png:
#     Headline horizontal bar chart of AUC per (consequence, predictor),
#     VALIDATION classes on top, DIAGNOSTIC below (EXPLORATORY omitted --
#     see Table 3 for those numbers instead)
#   insilico_concordance_boxplots.pdf:
#     One boxplot per (consequence, predictor) pair: score by anchor_tier
#
# ---------------------------------------------------------------------
# USAGE
# ---------------------------------------------------------------------
#   Rscript insilico_concordance.R \
#       --anchor_dir per_targeton_results/ \
#       --meta_dir meta_consequences/ \
#       --outdir results/
#
# ---------------------------------------------------------------------

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("Package 'optparse' is required. Install with install.packages('optparse').")
  }
  library(optparse)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with install.packages('ggplot2').")
  }
  library(ggplot2)
})

HAS_PROC <- requireNamespace("pROC", quietly = TRUE)
if (HAS_PROC) {
  suppressPackageStartupMessages(library(pROC))
} else {
  message("NOTE: package 'pROC' not found -- using a dependency-free Mann-Whitney-based AUC instead (numerically equivalent).")
}

# --------------------------- CLI args ---------------------------------

option_list <- list(
  make_option("--anchor_dir", type = "character", help = "Directory of per-targeton anchor-tier/DESeq2 output TSVs"),
  make_option("--meta_dir", type = "character", help = "Directory of per-targeton meta_consequences TSVs"),
  make_option("--outdir", type = "character", default = ".", help = "Output directory [default: current directory]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$anchor_dir) || is.null(opt$meta_dir)) {
  stop("Both --anchor_dir and --meta_dir are required. See header comment for usage.")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# --------------------------- Helpers -----------------------------------

infer_targeton_anchor <- function(fname) {
  sub("_all_deseq2_results_condition_.*", "", basename(fname))
}

infer_targeton_meta <- function(fname) {
  sub("_chr[0-9XYM]+_.*", "", basename(fname))
}

parse_bracket_score <- function(x) {
  # e.g. "tolerated(0.83)" -> 0.83 ; "-" -> NA
  x <- ifelse(x == "-", NA, x)
  as.numeric(sub(".*\\(([-0-9.]+)\\).*", "\\1", x))
}

parse_spliceai_max <- function(x) {
  # e.g. "SMC1A|0.00|0.00|0.01|0.00|-29|-36|44|36" -> max(0.00,0.00,0.01,0.00) = 0.01 ; "-" -> NA
  sapply(x, function(s) {
    if (is.na(s) || s == "-") return(NA_real_)
    parts <- strsplit(s, "\\|")[[1]]
    if (length(parts) < 5) return(NA_real_)
    vals <- suppressWarnings(as.numeric(parts[2:5]))
    if (all(is.na(vals))) return(NA_real_)
    max(vals, na.rm = TRUE)
  }, USE.NAMES = FALSE)
}

manual_auc <- function(scores, labels) {
  # AUC via the Mann-Whitney U relationship: AUC = U / (n_pos * n_neg).
  # labels: logical, TRUE = positive class (depleted).
  pos <- scores[labels]
  neg <- scores[!labels]
  pos <- pos[!is.na(pos)]
  neg <- neg[!is.na(neg)]
  if (length(pos) < 2 || length(neg) < 2) return(NA_real_)
  u <- suppressWarnings(wilcox.test(pos, neg, exact = FALSE)$statistic)
  as.numeric(u) / (length(pos) * length(neg))
}

compute_auc <- function(scores, labels) {
  if (HAS_PROC) {
    ok <- !is.na(scores) & !is.na(labels)
    if (sum(ok) < 4 || length(unique(labels[ok])) < 2) return(NA_real_)
    roc_obj <- tryCatch(
      pROC::roc(response = labels[ok], predictor = scores[ok], quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(roc_obj)) return(NA_real_)
    as.numeric(pROC::auc(roc_obj))
  } else {
    manual_auc(scores, labels)
  }
}

manual_youden <- function(scores, labels) {
  # Dependency-free Youden's J optimal threshold: the point on the ROC curve
  # that maximises sensitivity + specificity - 1. Tries every observed score
  # as a candidate threshold ("predicted positive" = score >= threshold).
  ok <- !is.na(scores) & !is.na(labels)
  scores <- scores[ok]; labels <- labels[ok]
  if (length(unique(labels)) < 2) return(list(threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_))

  thresholds <- sort(unique(scores))
  best_j <- -Inf
  best <- list(threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_)
  for (t in thresholds) {
    pred_pos <- scores >= t
    tp <- sum(pred_pos & labels); fn <- sum(!pred_pos & labels)
    tn <- sum(!pred_pos & !labels); fp <- sum(pred_pos & !labels)
    sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    if (is.na(sens) || is.na(spec)) next
    j <- sens + spec - 1
    if (j > best_j) {
      best_j <- j
      best <- list(threshold = t, sensitivity = sens, specificity = spec)
    }
  }
  best
}

compute_youden <- function(scores, labels) {
  # Returns the data-derived (Youden's J-optimal) threshold, and the
  # sensitivity/specificity of the predictor AT that threshold -- i.e. "if we
  # used this predictor's own best cutoff, what fraction of assay-depleted
  # variants would it catch (sensitivity), and what fraction of no-impact
  # variants would it correctly leave alone (specificity)".
  if (HAS_PROC) {
    ok <- !is.na(scores) & !is.na(labels)
    if (sum(ok) < 4 || length(unique(labels[ok])) < 2) {
      return(list(threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_))
    }
    roc_obj <- tryCatch(
      pROC::roc(response = labels[ok], predictor = scores[ok], quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(roc_obj)) return(list(threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_))
    best <- pROC::coords(roc_obj, "best", best.method = "youden",
                          ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
    if (is.data.frame(best) && nrow(best) > 1) best <- best[1, ]  # tie-break: take first
    list(threshold = as.numeric(best["threshold"]), sensitivity = as.numeric(best["sensitivity"]), specificity = as.numeric(best["specificity"]))
  } else {
    manual_youden(scores, labels)
  }
}

# --------------------------- Load + join --------------------------------

anchor_files <- list.files(opt$anchor_dir, pattern = "\\.tsv$", full.names = TRUE)
anchor_files <- anchor_files[!grepl("^gmm_shrinkage_anchor_summary", basename(anchor_files))]
meta_files <- list.files(opt$meta_dir, pattern = "\\.tsv$", full.names = TRUE)

if (length(anchor_files) == 0) stop(sprintf("No .tsv files found in --anchor_dir %s", opt$anchor_dir))
if (length(meta_files) == 0) stop(sprintf("No .tsv files found in --meta_dir %s", opt$meta_dir))

anchor_map <- setNames(anchor_files, sapply(anchor_files, infer_targeton_anchor))
meta_map <- setNames(meta_files, sapply(meta_files, infer_targeton_meta))

common_targetons <- intersect(names(anchor_map), names(meta_map))
if (length(common_targetons) == 0) {
  stop("No matching targetons found between --anchor_dir and --meta_dir filenames. Check naming conventions match the assumptions in the header comment.")
}
cat(sprintf("Found %d matching targetons: %s\n", length(common_targetons), paste(common_targetons, collapse = ", ")))

missing_anchor <- setdiff(names(anchor_map), names(meta_map))
missing_meta <- setdiff(names(meta_map), names(anchor_map))
if (length(missing_anchor) > 0) cat(sprintf("NOTE: %d anchor_dir targeton(s) have no matching meta_dir file, excluded: %s\n", length(missing_anchor), paste(missing_anchor, collapse = ", ")))
if (length(missing_meta) > 0) cat(sprintf("NOTE: %d meta_dir targeton(s) have no matching anchor_dir file, excluded: %s\n", length(missing_meta), paste(missing_meta, collapse = ", ")))

load_and_merge_one <- function(targeton) {
  anchor <- read.delim(anchor_map[[targeton]], stringsAsFactors = FALSE)
  meta <- read.delim(meta_map[[targeton]], stringsAsFactors = FALSE)

  required_anchor <- c("oligo_name", "consequence", "anchor_tier", "pos_adj_log2FoldChange_raw", "anchor_mu_noimpact", "anchor_sd_noimpact")
  missing_cols <- setdiff(required_anchor, names(anchor))
  if (length(missing_cols) > 0) stop(sprintf("%s is missing required column(s): %s", anchor_map[[targeton]], paste(missing_cols, collapse = ", ")))

  required_meta <- c("unique_oligo_name", "SIFT", "PolyPhen", "CADD_PHRED", "REVEL", "SpliceAI_pred")
  missing_cols <- setdiff(required_meta, names(meta))
  if (length(missing_cols) > 0) stop(sprintf("%s is missing required column(s): %s", meta_map[[targeton]], paste(missing_cols, collapse = ", ")))

  meta$oligo_name_stripped <- sub("_[0-9]+_[0-9]+$", "", meta$unique_oligo_name)

  merged <- merge(anchor, meta, by.x = "oligo_name", by.y = "oligo_name_stripped", all.x = TRUE)
  merged$Targeton_ID <- targeton

  merged$z_noimpact <- (merged$pos_adj_log2FoldChange_raw - merged$anchor_mu_noimpact) / merged$anchor_sd_noimpact

  merged$sift_damaging <- 1 - parse_bracket_score(merged$SIFT)
  merged$polyphen_score <- parse_bracket_score(merged$PolyPhen)
  merged$CADD_num <- suppressWarnings(as.numeric(ifelse(merged$CADD_PHRED == "-", NA, merged$CADD_PHRED)))
  merged$REVEL_num <- suppressWarnings(as.numeric(ifelse(merged$REVEL == "-", NA, merged$REVEL)))
  merged$spliceai_max <- parse_spliceai_max(merged$SpliceAI_pred)

  merged
}

all_data <- do.call(rbind, lapply(common_targetons, load_and_merge_one))
cat(sprintf("Loaded and merged %d total variants across %d targetons\n", nrow(all_data), length(common_targetons)))

n_unmatched <- sum(is.na(all_data$SpliceAI_pred) & is.na(all_data$SIFT))
if (n_unmatched > 0) {
  cat(sprintf("NOTE: %d variant(s) had no matching row in the meta_consequences file (all predictor columns NA for these) -- check join keys if this seems high.\n", n_unmatched))
}

# --------------------------- Write merged per-variant data ---------------
#
# The full merged table (assay results + cleaned predictor scores), one row
# per variant, BEFORE the "enriched" exclusion below -- so it's a complete
# reference file for looking up any specific variant (e.g. checking the
# SpliceAI score of a known outlier synonymous variant), not just the subset
# used in the tier analyses. Includes both the raw predictor strings (SIFT,
# PolyPhen, CADD_PHRED, REVEL, SpliceAI_pred) and the cleaned numeric/parsed
# versions used throughout this script, plus HGVSc/HGVSp for readability.
merged_out_cols <- c(
  "Targeton_ID", "oligo_name", "position", "consequence", "HGVSc", "HGVSp",
  "anchor_tier", "pos_adj_log2FoldChange_raw", "z_noimpact",
  "SIFT", "sift_damaging", "PolyPhen", "polyphen_score",
  "CADD_PHRED", "CADD_num", "REVEL", "REVEL_num", "SpliceAI_pred", "spliceai_max"
)
merged_out_cols <- intersect(merged_out_cols, names(all_data))  # tolerate a missing optional column (e.g. no HGVSc/HGVSp/position in some anchor files)
missing_optional <- setdiff(c("position", "HGVSc", "HGVSp"), names(all_data))
if (length(missing_optional) > 0) {
  cat(sprintf("NOTE: optional column(s) not found, omitted from merged output: %s\n", paste(missing_optional, collapse = ", ")))
}

merged_out_path <- file.path(opt$outdir, "merged_variant_data.tsv")
write.table(all_data[, merged_out_cols], merged_out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote merged per-variant data (%d rows, all consequence classes, enriched INCLUDED) to:\n  %s\n", nrow(all_data), merged_out_path))

# Exclude "enriched" calls throughout -- these predictors are not designed
# to predict gain-of-function-like assay behaviour.
n_before <- nrow(all_data)
all_data <- all_data[is.na(all_data$anchor_tier) | all_data$anchor_tier != "enriched", ]
cat(sprintf("Excluded %d 'enriched' variant(s) from all predictor comparisons.\n", n_before - nrow(all_data)))

all_data$is_depleted <- all_data$anchor_tier %in% c("weakly depleting", "strongly depleting")

# --------------------------- Analysis grid -----------------------------

# Which predictors are meaningful for which consequence class.
# SIFT/PolyPhen/REVEL are missense-only. CADD and SpliceAI are genome-wide.
# LOF is included with CADD/SpliceAI only -- useful for checking whether
# "non-depleting" LOF calls (e.g. the known exon-15 case) show elevated
# SpliceAI scores, consistent with an exon-skipping / cassette-exon story.
predictor_map <- list(
  Missense_Variant = c("sift_damaging", "polyphen_score", "REVEL_num", "CADD_num"),
  Synonymous_Variant = c("CADD_num", "spliceai_max"),
  Splice_Variant = c("CADD_num", "spliceai_max"),
  Splice_Polypyrimidine_Tract_Variant = c("CADD_num", "spliceai_max"),
  Intronic_Variant = c("CADD_num", "spliceai_max"),
  LOF = c("CADD_num", "spliceai_max"),
  Others = c("CADD_num", "spliceai_max")  # typically 3'UTR variants -- CADD applies genome-wide;
                                          # SpliceAI kept too in case any sit near the terminal exon boundary
)

pretty_names <- c(
  sift_damaging = "SIFT (1-score, higher=damaging)",
  polyphen_score = "PolyPhen",
  REVEL_num = "REVEL",
  CADD_num = "CADD_PHRED",
  spliceai_max = "SpliceAI (max delta)"
)

MIN_N <- 10  # minimum non-NA variants required to run a comparison
MIN_CLASS_N <- 30  # minimum size for BOTH the depleted and no-impact groups to trust an AUC/sensitivity result

# Which consequence classes represent genuine predictor VALIDATION (real
# heterogeneity of effect within the class -- missense ranges benign to
# severe, splice-region ranges benign to junction-abolishing) vs DIAGNOSTIC
# use (the class is expected to behave a certain way almost by definition,
# so a strong result says something about the ANCHOR/exceptions, not about
# predictor quality). This distinction is standard in the SGE literature:
# e.g. the DDX3X SGE paper (Radford et al.) tests SpliceAI specifically
# against canonical splice donor/acceptor variants (a validation use), while
# LOF/nonsense variants are universally treated as an anchor/reference class
# rather than a predictor test set, precisely because they are ~uniformly
# damaging by construction -- there's very little "no impact" LOF to
# discriminate against, so an AUC there mostly measures class imbalance,
# not predictor quality.
VALIDATION_CLASSES <- c("Missense_Variant", "Splice_Variant", "Splice_Polypyrimidine_Tract_Variant")
DIAGNOSTIC_CLASSES <- c("Synonymous_Variant", "LOF")
# Anything else (Intronic_Variant, Others/3'UTR) is exploratory by default --
# not enough prior literature precedent either way, and typically low n here.

classify_tier <- function(conseq, n_dep, n_ni) {
  base_tier <- if (conseq %in% VALIDATION_CLASSES) {
    "VALIDATION (standard predictor benchmark)"
  } else if (conseq %in% DIAGNOSTIC_CLASSES) {
    "DIAGNOSTIC (not a predictor-quality test -- see note)"
  } else {
    "EXPLORATORY (limited precedent for this class)"
  }
  imbalanced <- !is.na(n_dep) && !is.na(n_ni) && (min(n_dep, n_ni) < MIN_CLASS_N)
  if (imbalanced) {
    paste0(base_tier, " -- CAUTION: imbalanced (min class n=", min(n_dep, n_ni), " < ", MIN_CLASS_N, "), AUC/sensitivity unstable")
  } else {
    base_tier
  }
}

tier_note <- function(conseq) {
  if (conseq %in% VALIDATION_CLASSES) {
    "Real heterogeneity of effect expected within this class -- a strong AUC here is genuine evidence the predictor recovers assay-measured function."
  } else if (conseq == "Synonymous_Variant") {
    "Most synonymous variants are expected to be neutral -- a strong AUC does NOT validate the predictor generally. It means the predictor independently flags the same rare outlier variants the assay calls depleted, which is useful as a contamination/QC check (cf. splice-hotspot finding), not as a predictor-performance benchmark."
  } else if (conseq == "LOF") {
    "LOF variants are ~uniformly damaging by construction, so there is very little 'no impact' class to discriminate against -- AUC here mostly reflects class imbalance, not predictor quality. More useful read: do the RARE tolerated LOF exceptions get low predictor scores (consistent with e.g. an exon-skipping/NMD-escape mechanism)?"
  } else {
    "Limited precedent in the SGE literature for benchmarking predictors against this class specifically -- treat as exploratory."
  }
}

corr_rows <- list()
kw_rows <- list()
auc_rows <- list()

for (conseq in names(predictor_map)) {
  sub_df <- all_data[all_data$consequence == conseq & !is.na(all_data$consequence), ]
  if (nrow(sub_df) < MIN_N) next

  for (pred in predictor_map[[conseq]]) {
    vals <- sub_df[[pred]]
    z <- sub_df$z_noimpact
    ok <- !is.na(vals) & !is.na(z)
    if (sum(ok) < MIN_N) next

    # 1. Spearman correlation vs continuous z-score
    sp <- suppressWarnings(cor.test(vals[ok], z[ok], method = "spearman"))
    corr_rows[[length(corr_rows) + 1]] <- data.frame(
      consequence = conseq, predictor = pretty_names[[pred]], n = sum(ok),
      spearman_rho = round(unname(sp$estimate), 3), p_value = sp$p.value,
      tier = classify_tier(conseq, NA, NA),
      stringsAsFactors = FALSE
    )

    # 2. Kruskal-Wallis across anchor_tier (enriched already excluded above)
    tier_ok <- ok & !is.na(sub_df$anchor_tier)
    if (length(unique(sub_df$anchor_tier[tier_ok])) >= 2) {
      kw <- kruskal.test(vals[tier_ok], sub_df$anchor_tier[tier_ok])
      kw_rows[[length(kw_rows) + 1]] <- data.frame(
        consequence = conseq, predictor = pretty_names[[pred]], n = sum(tier_ok),
        kruskal_stat = round(unname(kw$statistic), 3), p_value = kw$p.value,
        tier = classify_tier(conseq, NA, NA),
        stringsAsFactors = FALSE
      )
    }

    # 3. AUC: depleted vs no-impact (anchor_tier as ground truth), plus the
    #    data-derived (Youden-optimal) threshold and sensitivity/specificity
    #    at that threshold -- i.e. "what fraction of our assay-depleted
    #    variants would this predictor's own best cutoff actually catch".
    auc_mask <- ok & sub_df$anchor_tier %in% c("no impact", "weakly depleting", "strongly depleting")
    n_dep <- sum(sub_df$is_depleted[auc_mask], na.rm = TRUE)
    n_ni <- sum(auc_mask) - n_dep
    if (sum(auc_mask) >= MIN_N && n_dep >= 2 && n_ni >= 2) {
      auc_val <- compute_auc(vals[auc_mask], sub_df$is_depleted[auc_mask])
      yj <- compute_youden(vals[auc_mask], sub_df$is_depleted[auc_mask])
      auc_rows[[length(auc_rows) + 1]] <- data.frame(
        consequence = conseq, predictor = pretty_names[[pred]],
        n_depleted = n_dep, n_no_impact = n_ni, AUC = round(auc_val, 3),
        youden_threshold = round(yj$threshold, 3),
        sensitivity_at_threshold = round(yj$sensitivity, 3),
        specificity_at_threshold = round(yj$specificity, 3),
        pct_depleted_variants_caught = round(100 * yj$sensitivity, 1),
        tier = classify_tier(conseq, n_dep, n_ni),
        stringsAsFactors = FALSE
      )
    }
  }
}

corr_df <- if (length(corr_rows) > 0) do.call(rbind, corr_rows) else data.frame()
kw_df <- if (length(kw_rows) > 0) do.call(rbind, kw_rows) else data.frame()
auc_df <- if (length(auc_rows) > 0) do.call(rbind, auc_rows) else data.frame()

# Sort so VALIDATION rows appear first, DIAGNOSTIC/EXPLORATORY last -- makes
# it obvious at a glance which results are standard predictor benchmarks vs
# which need the caveats in tier_notes_df below.
tier_sort_key <- function(df) {
  order(factor(substr(df$tier, 1, 10), levels = c("VALIDATION", "DIAGNOSTIC", "EXPLORATOR")))
}
if (nrow(corr_df) > 0) corr_df <- corr_df[tier_sort_key(corr_df), ]
if (nrow(kw_df) > 0) kw_df <- kw_df[tier_sort_key(kw_df), ]
if (nrow(auc_df) > 0) auc_df <- auc_df[tier_sort_key(auc_df), ]

# One-paragraph explanation per consequence class, for anyone reading the
# output cold without this conversation's context.
tier_notes_df <- data.frame(
  consequence = names(predictor_map),
  tier = sapply(names(predictor_map), function(c) classify_tier(c, NA, NA)),
  note = sapply(names(predictor_map), tier_note),
  stringsAsFactors = FALSE
)

# --------------------------- Write combined TSV -------------------------

out_tsv <- file.path(opt$outdir, "insilico_concordance_tables.tsv")
con <- file(out_tsv, "w")
writeLines("# Table 0: how to read the 'tier' column in Tables 1-3 -- VALIDATION classes are standard predictor benchmarks (real effect heterogeneity expected); DIAGNOSTIC classes are not predictor-quality tests (the class is expected to behave a certain way almost by construction -- a strong result says something about the assay's exceptions/anchor purity, not predictor quality); EXPLORATORY classes have limited literature precedent either way, usually low n. See per-class notes below.", con)
write.table(tier_notes_df, con, sep = "\t", row.names = FALSE, quote = FALSE)
writeLines("", con)
writeLines("# Table 1: Spearman correlation, predictor score vs z_noimpact (continuous), by consequence class -- sorted VALIDATION first", con)
write.table(corr_df, con, sep = "\t", row.names = FALSE, quote = FALSE)
writeLines("", con)
writeLines("# Table 2: Kruskal-Wallis test, predictor score across anchor_tier groups (enriched excluded) -- sorted VALIDATION first", con)
write.table(kw_df, con, sep = "\t", row.names = FALSE, quote = FALSE)
writeLines("", con)
writeLines("# Table 3: AUC, predictor score discriminating depleted vs no-impact (anchor_tier as ground truth), plus the data-derived (Youden-optimal) threshold and sensitivity/specificity at that threshold -- pct_depleted_variants_caught is the plain-language number: what fraction of assay-depleted variants this predictor's own best cutoff would flag. Sorted VALIDATION first.", con)
write.table(auc_df, con, sep = "\t", row.names = FALSE, quote = FALSE)
close(con)

cat("\n=== How to read the tiers below (see Table 0 in the output file for full notes) ===\n")
print(tier_notes_df[, c("consequence", "tier")], row.names = FALSE)
cat("\n=== Spearman correlation (predictor vs z_noimpact) -- VALIDATION first ===\n")
print(corr_df, row.names = FALSE)
cat("\n=== Kruskal-Wallis (predictor across anchor_tier) -- VALIDATION first ===\n")
print(kw_df, row.names = FALSE)
cat("\n=== AUC (depleted vs no-impact) -- VALIDATION first ===\n")
print(auc_df, row.names = FALSE)

# --------------------------- AUC summary chart ---------------------------
#
# Headline bar chart: AUC per (consequence, predictor), VALIDATION and
# DIAGNOSTIC only (EXPLORATORY dropped -- too few variants per class to
# present as a headline result; see insilico_concordance_tables.tsv for
# those numbers instead). VALIDATION sits on top, DIAGNOSTIC below, each
# sorted so its strongest result is nearest the group boundary.

summary_chart_path <- file.path(opt$outdir, "auc_summary_chart.png")
chart_df <- auc_df[grepl("^(VALIDATION|DIAGNOSTIC)", auc_df$tier), ]

if (nrow(chart_df) > 0) {
  chart_df$tier_short <- ifelse(grepl("^VALIDATION", chart_df$tier), "VALIDATION", "DIAGNOSTIC")
  chart_df$label <- paste0(
    gsub("_", " ", gsub("_Variant$", "", chart_df$consequence)), "  |  ",
    gsub(" \\(.*\\)", "", chart_df$predictor)
  )

  # DIAGNOSTIC plotted first (bottom), VALIDATION last (top); ascending AUC
  # within each group so the strongest bar of each group sits at that
  # group's outer edge (mirrors the chart used on the results slide).
  chart_df$tier_rank <- ifelse(chart_df$tier_short == "DIAGNOSTIC", 0, 1)
  chart_df <- chart_df[order(chart_df$tier_rank, chart_df$AUC), ]
  chart_df$label <- factor(chart_df$label, levels = chart_df$label)  # lock in plot order

  chart_colors <- c(VALIDATION = "#34A853", DIAGNOSTIC = "#9AA0A6")

  p_summary <- ggplot(chart_df, aes(x = label, y = AUC, fill = tier_short)) +
    geom_col(width = 0.65) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "#888888") +
    geom_text(aes(label = sprintf("%.2f", AUC)), hjust = -0.15, size = 3.3) +
    scale_fill_manual(values = chart_colors, name = NULL) +
    coord_flip(ylim = c(0.3, 1.03)) +
    labs(x = NULL, y = "AUC (predictor discriminating assay-depleted vs no-impact)") +
    theme_minimal() +
    theme(legend.position = "top", panel.grid.minor = element_blank())

  ggsave(summary_chart_path, p_summary, width = 10.5, height = max(4, 0.5 * nrow(chart_df) + 1), dpi = 200)
  cat(sprintf("Wrote AUC summary chart (VALIDATION + DIAGNOSTIC only) to:\n  %s\n", summary_chart_path))
} else {
  cat("NOTE: no VALIDATION/DIAGNOSTIC rows available -- skipped auc_summary_chart.png.\n")
}

# --------------------------- Boxplots ------------------------------------

plot_df <- all_data[all_data$anchor_tier %in% c("no impact", "weakly depleting", "strongly depleting"), ]
plot_df$anchor_tier <- factor(plot_df$anchor_tier, levels = c("no impact", "weakly depleting", "strongly depleting"))

boxplot_path <- file.path(opt$outdir, "insilico_concordance_boxplots.pdf")
pdf(boxplot_path, width = 8, height = 6)
n_plots <- 0
plot_order <- c(VALIDATION_CLASSES, DIAGNOSTIC_CLASSES, setdiff(names(predictor_map), c(VALIDATION_CLASSES, DIAGNOSTIC_CLASSES)))
for (conseq in plot_order) {
  if (!(conseq %in% names(predictor_map))) next
  sub_df <- plot_df[plot_df$consequence == conseq & !is.na(plot_df$consequence), ]
  if (nrow(sub_df) < MIN_N) next
  short_tier <- if (conseq %in% VALIDATION_CLASSES) "VALIDATION" else if (conseq %in% DIAGNOSTIC_CLASSES) "DIAGNOSTIC" else "EXPLORATORY"
  for (pred in predictor_map[[conseq]]) {
    d <- sub_df[!is.na(sub_df[[pred]]), c(pred, "anchor_tier")]
    if (nrow(d) < MIN_N) next
    names(d) <- c("score", "anchor_tier")
    p <- ggplot(d, aes(x = anchor_tier, y = score, fill = anchor_tier)) +
      geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
      geom_jitter(width = 0.15, size = 0.6, alpha = 0.3) +
      scale_fill_manual(values = c("no impact" = "#4285F4", "weakly depleting" = "#FBBC04", "strongly depleting" = "#EA4335")) +
      labs(title = paste0("[", short_tier, "] ", conseq, ": ", pretty_names[[pred]]),
           subtitle = sprintf("n=%d (enriched excluded)", nrow(d)),
           x = "anchor_tier", y = pretty_names[[pred]]) +
      theme_minimal() +
      theme(legend.position = "none")
    print(p)
    n_plots <- n_plots + 1
  }
}
dev.off()

cat(sprintf("\nWrote %d boxplot page(s).\n", n_plots))
cat("\nWrote:\n")
cat(sprintf("  %s\n", merged_out_path))
cat(sprintf("  %s\n", out_tsv))
cat(sprintf("  %s\n", summary_chart_path))
cat(sprintf("  %s\n", boxplot_path))
