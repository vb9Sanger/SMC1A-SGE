#!/usr/bin/env Rscript
# gmm_shrinkage_anchor_pipeline.R
#
# Partial-pooling ("empirical Bayes shrinkage") variant, sitting between
# gmm_gene_wide_anchor_pipeline.R (one shared anchor for the whole gene)
# and gmm_per_targeton_anchor_pipeline.R (each targeton fit in total
# isolation). Both of those are the two extremes of the same tradeoff:
# gene-wide borrows maximum statistical power but can't adapt to a
# targeton that's genuinely different; per-targeton adapts fully to each
# exon but small-n exons can be yanked around by a handful of noisy
# points. This script blends the two:
#
#   1. Compute the GENE-WIDE anchor exactly as in
#      gmm_gene_wide_anchor_pipeline.R (pooled LOF/no-impact controls
#      across all files, same NMD-escape-aware filtering, same exclusion
#      of pre-called-"enriched" rows from the fitting pool).
#
#   2. Compute each targeton's own LOCAL anchor exactly as in
#      gmm_per_targeton_anchor_pipeline.R (that targeton's own LOF/no-
#      impact rows only, same NMD-escape and enriched-exclusion logic,
#      scoped to that file).
#
#   3. SHRINK each targeton's local anchor toward the gene-wide anchor,
#      weighted by how much data that targeton actually has, using a
#      simple pseudo-count rule (the same logic as DESeq2's own
#      dispersion shrinkage -- borrowing strength across features rather
#      than trusting each one in total isolation):
#
#        mu_shrunk_i  = (n_i * mu_local_i  + k * mu_global)  / (n_i + k)
#        var_shrunk_i = (n_i * var_local_i + k * var_global) / (n_i + k)
#        sd_shrunk_i  = sqrt(var_shrunk_i)
#
#      applied separately to the LOF anchor (using each targeton's n_lof
#      and --shrinkage_k_lof) and the no-impact anchor (n_ctrl and
#      --shrinkage_k_noimpact). A well-powered targeton (n_i >> k) stays
#      close to its own local estimate; a thin targeton (n_i << k) is
#      pulled close to the gene-wide value; as n_i -> 0 the formula
#      gracefully degrades to the pure gene-wide estimate.
#
#      k defaults to "auto", which sets it to the MEDIAN per-targeton n
#      across the gene (separately for LOF and no-impact) -- i.e. a
#      "typical" targeton gets roughly equal weight between its own data
#      and the gene-wide prior. Override with a fixed number via
#      --shrinkage_k_lof / --shrinkage_k_noimpact if you want more or
#      less aggressive shrinkage.
#
#   4. Classify every row in every targeton (any consequence class) using
#      THAT TARGETON'S OWN SHRUNK anchors and shrunk 95th-percentile LOF
#      threshold -- same posterior-probability depleted/no-impact call
#      and weak/strong tiering rule as both other scripts.
#
# Final column: anchor_tier, same 4 levels as the other two scripts
# ("enriched" / "no impact" / "weakly depleting" / "strongly depleting").
#
# Per-row output also includes the pre-shrinkage local and gene-wide
# values side by side with the final shrunk ones, plus the shrinkage
# weight actually applied, so it's fully auditable how much each
# targeton's own data vs. the gene-wide prior contributed.
#
# Dependencies: data.table (required). ggplot2 + ragg (only if --plot_dir set).
#
# Usage:
#   Rscript gmm_shrinkage_anchor_pipeline.R \
#     --input "deseq2_results/*_all_deseq2_results_condition_*.tsv" \
#     --out_dir GMM_shrinkage_results \
#     --plot_dir GMM_shrinkage_results/plots \
#     --exon_map SMC1A_exon_map.tsv
#
# Full flag list: run with --help

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Requires data.table")
})
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

# ---- Helpers ----------------------------------------------------------------

print_help <- function() {
  cat(
"gmm_shrinkage_anchor_pipeline.R

Required:
  --input <glob>              e.g. \"deseq2_results/*_all_deseq2_results_condition_*.tsv\"
  --out_dir <dir>

Optional:
  --lfc_col <name>             default: pos_adj_log2FoldChange_raw
  --position_col <name>         default: position (only needed if --exon_map is used)
  --consequence_col <name>     default: consequence
  --lof_label <name>           default: LOF
  --no_impact_labels <csv>     default: Synonymous_Variant,Intronic_Variant
  --missense_label <name>      default: Missense_Variant (used only for the missense-specific diagnostic plot)
  --status_col <name>          default: stat_pos_raw (set to \"\" to disable enriched exclusion/passthrough)
  --lof_percentile <0-1>       default: 0.95
  --post_threshold <0-1>       default: 0.5
  --min_controls <int>         default: 20 (minimum n per anchor class, GENE-WIDE pool)
  --min_local_n <int>          default: 5 (minimum per-targeton n to compute a local estimate
                                 at all; below this, that targeton's anchor is pure gene-wide,
                                 i.e. local weight = 0, rather than fitting a variance off a
                                 handful of points)
  --shrinkage_k_lof <num|auto>       default: auto (median per-targeton n_lof across the gene)
  --shrinkage_k_noimpact <num|auto> default: auto (median per-targeton n_no_impact across the gene)
  --exon_map <path>             optional TSV/CSV with columns Exon_position (chrN:start-end),
                                 EXON (integer, transcript order), Targeton_ID.
                                 If given, enables NMD-escape-aware LOF filtering, applied to
                                 both the gene-wide pool and each targeton's local pool.
  --nmd_escape_distance <int>   default: 50 (bp into the penultimate exon, from the final
                                 exon-exon junction, also treated as plausible escape zone)
  --plot_dir <dir>             writes 2 diagnostic PNGs per exon: one all-variants, one missense-only
  --suffix <string>            appended before .tsv in output filenames (default: \"\")

Example:
  Rscript gmm_shrinkage_anchor_pipeline.R \\
    --input \"deseq2_results/*_all_deseq2_results_condition_*.tsv\" \\
    --out_dir GMM_shrinkage_results --plot_dir GMM_shrinkage_results/plots

", sep = "")
}

get_flag <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx[1] == length(args)) stop("Flag provided without value: ", flag)
  args[idx[1] + 1]
}
has_flag <- function(flag) flag %in% args

if (has_flag("--help") || has_flag("-h")) { print_help(); quit(save = "no", status = 0) }

# ---- Parse args -------------------------------------------------------------

input_glob <- get_flag("--input")
out_dir    <- get_flag("--out_dir")
if (is.null(input_glob) || is.null(out_dir)) {
  print_help()
  stop("Missing required arg(s): --input and --out_dir are both required.")
}

lfc_col              <- get_flag("--lfc_col", "pos_adj_log2FoldChange_raw")
position_col_for_nmd <- get_flag("--position_col", "position")
consequence_col      <- get_flag("--consequence_col", "consequence")
lof_label            <- get_flag("--lof_label", "LOF")
no_impact_labels     <- strsplit(get_flag("--no_impact_labels", "Synonymous_Variant,Intronic_Variant"), ",")[[1]]
missense_label       <- get_flag("--missense_label", "Missense_Variant")
status_col           <- get_flag("--status_col", "stat_pos_raw")
lof_percentile       <- as.numeric(get_flag("--lof_percentile", "0.95"))
post_threshold       <- as.numeric(get_flag("--post_threshold", "0.5"))
min_controls         <- as.integer(get_flag("--min_controls", "20"))
min_local_n          <- as.integer(get_flag("--min_local_n", "5"))
shrinkage_k_lof_arg       <- get_flag("--shrinkage_k_lof", "auto")
shrinkage_k_noimpact_arg  <- get_flag("--shrinkage_k_noimpact", "auto")
exon_map_path        <- get_flag("--exon_map", NULL)
nmd_escape_distance  <- as.integer(get_flag("--nmd_escape_distance", "50"))
plot_dir             <- get_flag("--plot_dir", NULL)
suffix               <- get_flag("--suffix", "")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat("Created output directory: ", out_dir, "\n", sep = "")
}
plotting_enabled <- FALSE
if (!is.null(plot_dir)) {
  have_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
  have_ragg    <- requireNamespace("ragg", quietly = TRUE)
  if (have_ggplot2 && have_ragg) {
    plotting_enabled <- TRUE
    if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  } else {
    cat("--plot_dir given but ggplot2/ragg not both available; skipping all plots.\n")
  }
}

files <- Sys.glob(input_glob)
if (length(files) == 0) stop("No files matched: ", input_glob)
cat("Matched ", length(files), " file(s).\n\n", sep = "")

infer_targeton_name <- function(path) {
  bn <- basename(path)
  bn <- sub("\\.tsv$", "", bn, ignore.case = TRUE)
  bn <- sub("_all_deseq2_results_condition_.*$", "", bn, ignore.case = TRUE)
  bn
}

not_enriched_mask_for <- function(dt) {
  if (!is.null(status_col) && nzchar(status_col) && status_col %in% names(dt)) {
    s <- tolower(trimws(as.character(dt[[status_col]])))
    is.na(s) | s != "enriched"
  } else {
    rep(TRUE, nrow(dt))
  }
}

# ---- Derive NMD-escape zone boundaries ONCE from --exon_map (if given) -----

exon_map_available <- FALSE
if (!is.null(exon_map_path)) {
  if (!file.exists(exon_map_path)) stop("--exon_map file not found: ", exon_map_path)
  emap <- fread(exon_map_path, header = TRUE, data.table = TRUE)  # sep auto-detected: works for csv or tsv
  req_emap_cols <- c("Exon_position", "EXON")
  missing_emap <- setdiff(req_emap_cols, names(emap))
  if (length(missing_emap) > 0) stop("--exon_map is missing required column(s): ", paste(missing_emap, collapse = ", "))

  emap[, c("emap_chrom", "emap_range") := tstrsplit(Exon_position, ":", fixed = TRUE)]
  emap[, c("emap_start", "emap_end") := tstrsplit(emap_range, "-", fixed = TRUE)]
  emap[, emap_start := as.numeric(emap_start)]
  emap[, emap_end   := as.numeric(emap_end)]
  emap_u <- unique(emap, by = "EXON")
  setorder(emap_u, EXON)

  if (nrow(emap_u) < 2) stop("--exon_map needs at least 2 distinct EXON values to infer strand/terminal exon.")

  strand_is_minus <- cor(emap_u$EXON, emap_u$emap_start) < 0
  strand_str <- if (strand_is_minus) "minus" else "plus"

  terminal_exon_num <- max(emap_u$EXON)
  penult_exon_num    <- terminal_exon_num - 1
  terminal_row <- emap_u[EXON == terminal_exon_num]
  penult_row   <- emap_u[EXON == penult_exon_num]

  if (nrow(terminal_row) == 0 || nrow(penult_row) == 0) {
    stop("Could not find terminal (EXON=", terminal_exon_num, ") and/or penultimate (EXON=",
         penult_exon_num, ") exon rows in --exon_map. Check EXON numbering is consecutive.")
  }

  if (strand_is_minus) {
    pen_zone_start <- penult_row$emap_start
    pen_zone_end   <- penult_row$emap_start + nmd_escape_distance - 1
  } else {
    pen_zone_start <- penult_row$emap_end - nmd_escape_distance + 1
    pen_zone_end   <- penult_row$emap_end
  }

  cat(sprintf("NMD-escape zone from --exon_map: strand=%s | terminal exon %d [%d-%d] | penultimate exon %d junction-proximal [%d-%d]\n\n",
              strand_str, terminal_exon_num, terminal_row$emap_start, terminal_row$emap_end,
              penult_exon_num, pen_zone_start, pen_zone_end))

  nmd_zone_start1 <- terminal_row$emap_start; nmd_zone_end1 <- terminal_row$emap_end
  nmd_zone_start2 <- pen_zone_start;          nmd_zone_end2 <- pen_zone_end
  exon_map_available <- TRUE
}

# Returns a logical vector: is this position within the plausible
# NMD-escape zone (terminal exon or penultimate-exon junction buffer)?
# Always length(pos), FALSE everywhere if no --exon_map was given.
nmd_zone_mask <- function(pos) {
  if (!exon_map_available) return(rep(FALSE, length(pos)))
  (pos >= nmd_zone_start1 & pos <= nmd_zone_end1) |
    (pos >= nmd_zone_start2 & pos <= nmd_zone_end2)
}

# Applies the self-correcting NMD-escape filter to a set of LOF values/
# positions: a LOF variant is "true escape" only if it's in the plausible
# zone AND already looks non-depleting under an initial, uncorrected fit
# on this same set of values. Returns a logical "true escape" mask (same
# length as lof_vals).
apply_nmd_filter <- function(lof_vals, lof_pos, ctl_vals) {
  if (!exon_map_available || length(lof_vals) == 0) return(rep(FALSE, length(lof_vals)))
  in_zone <- nmd_zone_mask(lof_pos)
  if (sum(in_zone, na.rm = TRUE) == 0 || length(ctl_vals) == 0) return(in_zone & FALSE)
  mu0 <- mean(lof_vals); sd0 <- sd(lof_vals)
  muc0 <- mean(ctl_vals); sdc0 <- sd(ctl_vals)
  if (!is.finite(sd0) || sd0 == 0 || !is.finite(sdc0) || sdc0 == 0) return(in_zone & FALSE)
  d0 <- dnorm(lof_vals, mu0, sd0); dc0 <- dnorm(lof_vals, muc0, sdc0)
  post0 <- d0 / (d0 + dc0)
  in_zone & (post0 <= post_threshold)
}

# ---- PASS 1: read every file; compute gene-wide pool AND each targeton's ----
# ---- own local (pre-shrinkage) anchor values --------------------------------

cat("PASS 1: reading files, computing gene-wide pool and each targeton's local anchor...\n")

dt_list <- vector("list", length(files))
targeton_names <- character(length(files))
local_stats <- vector("list", length(files))  # one row per targeton: n/mu/sd for lof & ctl (pre-shrinkage)

pooled_lof_vals <- list(); pooled_lof_pos <- list(); pooled_ctl_vals <- list()

for (i in seq_along(files)) {
  f <- files[[i]]
  targeton_names[i] <- infer_targeton_name(f)
  dt <- tryCatch(
    fread(f, sep = "\t", header = TRUE, data.table = TRUE),
    error = function(e) stop("Failed reading TSV: ", f, "\n", conditionMessage(e))
  )
  setDT(dt)
  required_cols <- c(lfc_col, consequence_col)
  if (exon_map_available) required_cols <- c(required_cols, position_col_for_nmd)
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("File ", basename(f), " is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  lfc_val <- suppressWarnings(as.numeric(dt[[lfc_col]]))
  conseq  <- as.character(dt[[consequence_col]])
  pos_val <- if (exon_map_available) suppressWarnings(as.numeric(dt[[position_col_for_nmd]])) else rep(NA_real_, length(lfc_val))
  finite_mask <- is.finite(lfc_val)
  not_enr <- not_enriched_mask_for(dt)

  lof_mask <- !is.na(conseq) & conseq == lof_label & finite_mask & not_enr
  ctl_mask <- !is.na(conseq) & conseq %in% no_impact_labels & finite_mask & not_enr

  lof_vals_i <- lfc_val[lof_mask]
  lof_pos_i  <- pos_val[lof_mask]
  ctl_vals_i <- lfc_val[ctl_mask]

  # NMD-escape filtering, scoped to this targeton's own LOF pool (used both
  # for this targeton's local estimate AND contributes the filtered values
  # to the gene-wide pool, so the two stay consistent with each other).
  true_escape_i <- apply_nmd_filter(lof_vals_i, lof_pos_i, ctl_vals_i)
  lof_vals_i_clean <- lof_vals_i[!true_escape_i]

  n_lof_i <- length(lof_vals_i_clean)
  n_ctl_i <- length(ctl_vals_i)

  local_stats[[i]] <- data.table(
    targeton = targeton_names[i],
    n_lof_before_nmd = length(lof_vals_i), n_lof = n_lof_i,
    n_nmd_zone = sum(nmd_zone_mask(lof_pos_i), na.rm = TRUE),
    n_nmd_true_escape = sum(true_escape_i, na.rm = TRUE),
    mu_lof_local = if (n_lof_i >= min_local_n) mean(lof_vals_i_clean) else NA_real_,
    var_lof_local = if (n_lof_i >= min_local_n) var(lof_vals_i_clean) else NA_real_,
    n_ctl = n_ctl_i,
    mu_ctl_local = if (n_ctl_i >= min_local_n) mean(ctl_vals_i) else NA_real_,
    var_ctl_local = if (n_ctl_i >= min_local_n) var(ctl_vals_i) else NA_real_
  )

  pooled_lof_vals[[i]] <- lof_vals_i_clean
  pooled_ctl_vals[[i]] <- ctl_vals_i

  dt_list[[i]] <- dt
}

local_dt <- rbindlist(local_stats)

# ---- Gene-wide pooled anchor (same definition as gmm_gene_wide_anchor_pipeline.R) --

lof_vals_global <- unlist(pooled_lof_vals)
ctl_vals_global <- unlist(pooled_ctl_vals)
n_lof_global <- length(lof_vals_global)
n_ctl_global <- length(ctl_vals_global)

if (n_lof_global < min_controls || n_ctl_global < min_controls) {
  stop("Insufficient gene-wide controls to fit anchors (need >= ", min_controls, " each). ",
       "LOF n=", n_lof_global, ", no-impact n=", n_ctl_global, ".")
}

mu_lof_global <- mean(lof_vals_global); var_lof_global <- var(lof_vals_global)
mu_ctl_global <- mean(ctl_vals_global); var_ctl_global <- var(ctl_vals_global)

cat(sprintf("Gene-wide pooled anchor: mu_lof=%.4f sd_lof=%.4f (n=%d) | mu_no_impact=%.4f sd_no_impact=%.4f (n=%d)\n\n",
            mu_lof_global, sqrt(var_lof_global), n_lof_global,
            mu_ctl_global, sqrt(var_ctl_global), n_ctl_global))

# ---- Determine shrinkage strength k (auto = median per-targeton n) ---------

k_lof <- if (identical(shrinkage_k_lof_arg, "auto")) {
  median(local_dt$n_lof, na.rm = TRUE)
} else {
  as.numeric(shrinkage_k_lof_arg)
}
k_ctl <- if (identical(shrinkage_k_noimpact_arg, "auto")) {
  median(local_dt$n_ctl, na.rm = TRUE)
} else {
  as.numeric(shrinkage_k_noimpact_arg)
}

cat(sprintf("Shrinkage strength: k_lof=%.1f (%s) | k_no_impact=%.1f (%s)\n\n",
            k_lof, if (identical(shrinkage_k_lof_arg, "auto")) "auto=median n_lof" else "user-specified",
            k_ctl, if (identical(shrinkage_k_noimpact_arg, "auto")) "auto=median n_no_impact" else "user-specified"))

# ---- Compute each targeton's SHRUNK anchor ----------------------------------

local_dt[, `:=`(
  n_lof_for_shrink = fifelse(is.na(mu_lof_local), 0, as.numeric(n_lof)),
  n_ctl_for_shrink = fifelse(is.na(mu_ctl_local), 0, as.numeric(n_ctl))
)]
local_dt[, mu_lof_shrunk := (n_lof_for_shrink * fifelse(is.na(mu_lof_local), mu_lof_global, mu_lof_local) + k_lof * mu_lof_global) /
                            (n_lof_for_shrink + k_lof)]
local_dt[, var_lof_shrunk := (n_lof_for_shrink * fifelse(is.na(var_lof_local), var_lof_global, var_lof_local) + k_lof * var_lof_global) /
                             (n_lof_for_shrink + k_lof)]
local_dt[, mu_ctl_shrunk := (n_ctl_for_shrink * fifelse(is.na(mu_ctl_local), mu_ctl_global, mu_ctl_local) + k_ctl * mu_ctl_global) /
                            (n_ctl_for_shrink + k_ctl)]
local_dt[, var_ctl_shrunk := (n_ctl_for_shrink * fifelse(is.na(var_ctl_local), var_ctl_global, var_ctl_local) + k_ctl * var_ctl_global) /
                             (n_ctl_for_shrink + k_ctl)]
local_dt[, sd_lof_shrunk := sqrt(var_lof_shrunk)]
local_dt[, sd_ctl_shrunk := sqrt(var_ctl_shrunk)]
local_dt[, weight_lof_local := n_lof_for_shrink / (n_lof_for_shrink + k_lof)]
local_dt[, weight_ctl_local := n_ctl_for_shrink / (n_ctl_for_shrink + k_ctl)]

local_dt[, direction_lof_lower := mu_lof_shrunk < mu_ctl_shrunk]
local_dt[, lof_threshold := fifelse(
  direction_lof_lower,
  qnorm(lof_percentile, mean = mu_lof_shrunk, sd = sd_lof_shrunk),
  qnorm(1 - lof_percentile, mean = mu_lof_shrunk, sd = sd_lof_shrunk)
)]

cat("Per-targeton shrinkage summary:\n")
print(local_dt[, .(targeton, n_lof, weight_lof_local = round(weight_lof_local, 2),
                    mu_lof_local = round(mu_lof_local, 3), mu_lof_shrunk = round(mu_lof_shrunk, 3),
                    n_ctl, weight_ctl_local = round(weight_ctl_local, 2),
                    mu_ctl_local = round(mu_ctl_local, 3), mu_ctl_shrunk = round(mu_ctl_shrunk, 3))])
cat("\n")

# ---- PASS 2: classify every row in every file against ITS OWN SHRUNK anchor --

cat("PASS 2: classifying every row per exon against its shrunk anchor...\n\n")

summary_rows <- vector("list", length(files))

for (i in seq_along(files)) {
  f  <- files[[i]]
  bn <- basename(f)
  targeton <- targeton_names[i]
  dt <- dt_list[[i]]
  ts <- local_dt[targeton == targeton_names[i]]

  mu_lof <- ts$mu_lof_shrunk; sd_lof <- ts$sd_lof_shrunk
  mu_ctl <- ts$mu_ctl_shrunk; sd_ctl <- ts$sd_ctl_shrunk
  direction_lof_lower <- ts$direction_lof_lower
  lof_threshold <- ts$lof_threshold

  lfc_val <- suppressWarnings(as.numeric(dt[[lfc_col]]))
  finite_mask <- is.finite(lfc_val)

  d_lof <- dnorm(lfc_val, mean = mu_lof, sd = sd_lof)
  d_ctl <- dnorm(lfc_val, mean = mu_ctl, sd = sd_ctl)
  post_lof <- d_lof / (d_lof + d_ctl)

  anchor_call <- ifelse(!finite_mask, NA_character_,
                         ifelse(post_lof > post_threshold, "depleted", "no impact"))

  # Scalar direction_lof_lower -> plain if/else (NOT nested inside ifelse(),
  # which would shape its result to length(test) == 1 and silently recycle
  # a single value across every row).
  if (direction_lof_lower) {
    weak_or_strong <- ifelse(lfc_val > lof_threshold, "weakly depleting", "strongly depleting")
  } else {
    weak_or_strong <- ifelse(lfc_val < lof_threshold, "weakly depleting", "strongly depleting")
  }

  anchor_tier <- ifelse(!finite_mask, NA_character_,
                         ifelse(anchor_call == "depleted", weak_or_strong, "no impact"))

  if (!is.null(status_col) && nzchar(status_col) && status_col %in% names(dt)) {
    orig_status <- tolower(trimws(as.character(dt[[status_col]])))
    is_enriched <- !is.na(orig_status) & orig_status == "enriched"
    anchor_call[is_enriched] <- "enriched"
    anchor_tier[is_enriched] <- "enriched"
  }

  dt[, anchor_mu_lof := mu_lof]
  dt[, anchor_sd_lof := sd_lof]
  dt[, anchor_mu_lof_local := ts$mu_lof_local]
  dt[, anchor_mu_lof_global := mu_lof_global]
  dt[, anchor_weight_lof_local := ts$weight_lof_local]
  dt[, anchor_mu_noimpact := mu_ctl]
  dt[, anchor_sd_noimpact := sd_ctl]
  dt[, anchor_mu_noimpact_local := ts$mu_ctl_local]
  dt[, anchor_mu_noimpact_global := mu_ctl_global]
  dt[, anchor_weight_noimpact_local := ts$weight_ctl_local]
  dt[, anchor_direction := ifelse(direction_lof_lower, "lof_is_lower", "lof_is_higher")]
  dt[, anchor_lof_threshold := lof_threshold]
  dt[, anchor_post_lof := post_lof]
  dt[, anchor_call := anchor_call]
  dt[, anchor_tier := anchor_tier]

  if (plotting_enabled) {
    conseq <- as.character(dt[[consequence_col]])
    lof_x <- lfc_val[conseq == lof_label]
    ctl_x <- lfc_val[conseq %in% no_impact_labels]

    tryCatch({
      ragg::agg_png(file.path(plot_dir, paste0(targeton, "_anchor_fit.png")),
                    width = 900, height = 600, res = 120)
      other_x <- lfc_val[finite_mask]
      xr <- range(c(lof_x, ctl_x, other_x), na.rm = TRUE)
      hist(other_x, breaks = 40, freq = FALSE, col = rgb(0.6, 0.6, 0.6, 0.5), border = NA,
           xlim = xr, main = paste0(targeton, ": shrunk anchor fit applied (all variants)"), xlab = lfc_col)
      if (length(ctl_x) > 0) hist(ctl_x, breaks = 15, freq = FALSE, col = rgb(0.2, 0.5, 0.9, 0.5), border = NA, add = TRUE)
      if (length(lof_x) > 0) hist(lof_x, breaks = 15, freq = FALSE, col = rgb(0.9, 0.2, 0.2, 0.5), border = NA, add = TRUE)
      xs <- seq(xr[1], xr[2], length.out = 400)
      lines(xs, dnorm(xs, mu_lof, sd_lof), col = "red", lwd = 2)
      lines(xs, dnorm(xs, mu_ctl, sd_ctl), col = "blue", lwd = 2)
      lines(xs, dnorm(xs, mu_lof_global, sqrt(var_lof_global)), col = "red", lwd = 1, lty = 3)
      lines(xs, dnorm(xs, mu_ctl_global, sqrt(var_ctl_global)), col = "blue", lwd = 1, lty = 3)
      abline(v = lof_threshold, lty = 2, lwd = 2)
      legend("topleft", bty = "n",
             legend = c("all variants (grey)", "no-impact (blue)", "LOF (red)",
                        "shrunk fit (solid)", "gene-wide fit (dotted)", "shrunk 95th pct LOF threshold"),
             fill = c(rgb(0.6,0.6,0.6,0.5), rgb(0.2,0.5,0.9,0.5), rgb(0.9,0.2,0.2,0.5), NA, NA, NA),
             border = NA, lty = c(NA, NA, NA, 1, 3, 2), col = c(NA, NA, NA, "black", "black", "black"),
             cex = 0.8)
    }, error = function(e) cat("  (all-variants plot failed for ", targeton, ": ", conditionMessage(e), ")\n", sep = ""),
    finally = { if (!is.null(dev.list())) dev.off() })

    if (missense_label %in% conseq) {
      tryCatch({
        ragg::agg_png(file.path(plot_dir, paste0(targeton, "_missense_anchor_fit.png")),
                      width = 900, height = 600, res = 120)
        mis_x <- lfc_val[conseq == missense_label & finite_mask]
        xr <- range(c(lof_x, ctl_x, mis_x), na.rm = TRUE)
        hist(mis_x, breaks = 40, freq = FALSE, col = rgb(0.6, 0.6, 0.6, 0.5), border = NA,
             xlim = xr, main = paste0(targeton, ": shrunk anchor fit (missense only)"), xlab = lfc_col)
        if (length(ctl_x) > 0) hist(ctl_x, breaks = 15, freq = FALSE, col = rgb(0.2, 0.5, 0.9, 0.5), border = NA, add = TRUE)
        if (length(lof_x) > 0) hist(lof_x, breaks = 15, freq = FALSE, col = rgb(0.9, 0.2, 0.2, 0.5), border = NA, add = TRUE)
        xs <- seq(xr[1], xr[2], length.out = 400)
        lines(xs, dnorm(xs, mu_lof, sd_lof), col = "red", lwd = 2)
        lines(xs, dnorm(xs, mu_ctl, sd_ctl), col = "blue", lwd = 2)
        lines(xs, dnorm(xs, mu_lof_global, sqrt(var_lof_global)), col = "red", lwd = 1, lty = 3)
        lines(xs, dnorm(xs, mu_ctl_global, sqrt(var_ctl_global)), col = "blue", lwd = 1, lty = 3)
        abline(v = lof_threshold, lty = 2, lwd = 2)
        legend("topleft", bty = "n",
               legend = c("missense (grey)", "no-impact (blue)", "LOF (red)",
                          "shrunk fit (solid)", "gene-wide fit (dotted)", "shrunk 95th pct LOF threshold"),
               fill = c(rgb(0.6,0.6,0.6,0.5), rgb(0.2,0.5,0.9,0.5), rgb(0.9,0.2,0.2,0.5), NA, NA, NA),
               border = NA, lty = c(NA, NA, NA, 1, 3, 2), col = c(NA, NA, NA, "black", "black", "black"),
               cex = 0.8)
      }, error = function(e) cat("  (missense plot failed for ", targeton, ": ", conditionMessage(e), ")\n", sep = ""),
      finally = { if (!is.null(dev.list())) dev.off() })
    }
  }

  out_name <- if (nzchar(suffix)) sub("\\.tsv$", paste0(suffix, ".tsv"), bn, ignore.case = TRUE) else bn
  out_path <- file.path(out_dir, out_name)
  fwrite(dt, file = out_path, sep = "\t", quote = FALSE)

  tier_counts <- table(factor(anchor_tier, levels = c("no impact", "weakly depleting", "strongly depleting", "enriched")))
  cat(sprintf("%-10s  no_impact=%-4d weak=%-4d strong=%-4d enriched=%-4d  -> %s\n",
              targeton, tier_counts["no impact"], tier_counts["weakly depleting"],
              tier_counts["strongly depleting"], tier_counts["enriched"], basename(out_path)))

  summary_rows[[i]] <- data.table(
    targeton = targeton, input_file = f, output_file = out_path,
    n_lof_local = ts$n_lof, weight_lof_local = ts$weight_lof_local,
    mu_lof_local = ts$mu_lof_local, mu_lof_global = mu_lof_global, mu_lof_shrunk = mu_lof,
    sd_lof_shrunk = sd_lof,
    n_ctl_local = ts$n_ctl, weight_ctl_local = ts$weight_ctl_local,
    mu_ctl_local = ts$mu_ctl_local, mu_ctl_global = mu_ctl_global, mu_ctl_shrunk = mu_ctl,
    sd_ctl_shrunk = sd_ctl,
    lof_threshold = lof_threshold,
    n_lof_in_nmd_zone = ts$n_nmd_zone, n_lof_excluded_true_nmd_escape = ts$n_nmd_true_escape,
    n_no_impact = as.integer(tier_counts["no impact"]),
    n_weakly_depleting = as.integer(tier_counts["weakly depleting"]),
    n_strongly_depleting = as.integer(tier_counts["strongly depleting"]),
    n_enriched = as.integer(tier_counts["enriched"])
  )
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
summary_path <- file.path(out_dir, "gmm_shrinkage_anchor_summary.tsv")
fwrite(summary_dt, summary_path, sep = "\t")

cat("\n---\n")
cat("Wrote per-targeton summary (local/global/shrunk fit params + tier counts): ", summary_path, "\n", sep = "")
if (plotting_enabled) cat("Wrote diagnostic plots to: ", plot_dir, "\n", sep = "")
