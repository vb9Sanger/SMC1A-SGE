#!/usr/bin/env Rscript

# qc_signal_noise_metrics_with_composite_V3_anchor.R
#
# For each all_deseq2_results_condition_*.tsv file, compute QC metrics for
# comparing reference and LOESS choices across screens.
#
# Recognises the depletion-status column from either the anchor pipelines
# (anchor_tier, falling back to anchor_call) or the older un-anchored GMM
# (GMM_status), whichever is present -- so this works unchanged whether
# --input points at gmm_shrinkage_anchor_pipeline.R output, the gene-wide/
# per-targeton anchor variants, or the original GMM_status files.
#
# --input accepts multiple comma-separated glob patterns, so results from
# separate reference-condition runs (e.g. a D4-referenced shrinkage-anchor
# output directory and a separately-run Plasmid-referenced one -- the
# anchor pipelines refuse to run both references in one pass, since they
# are different statistical comparisons, not replicate data for the same
# exon) can be combined into a single comparison here.
#
# Outputs:
#   <out_prefix>.metrics.tsv
#   <out_prefix>.median_abs_z.png
#   <out_prefix>.median_abs_z_functional.png
#   <out_prefix>.composite_score.png
#
# By default, the metrics TSV is compact and contains only the main
# decision-making columns. Use --verbose_output to include all intermediate
# metrics and ranking components.
#
# LOESS NOTE:
#   The upstream pipeline fits LOESS #1 to remove positional bias from raw
#   counts/LFCs. This script fits LOESS #2 purely as a diagnostic: after
#   LOESS #1 corrected the data, how much positional structure remains in
#   neutral/control variants? A flat residual (low residual_positional_bias_*)
#   means LOESS #1 did its job well. The span for this diagnostic LOESS is
#   adaptive by default: wider when few control variants are available
#   (reducing overfitting of noise), narrower when many are available
#   (allowing detection of genuine residual waves).
#
# RANKING:
#   Within each targeton, candidate analysis conditions are ranked
#   hierarchically:
#     1. Passing basic filters (controls present; positional bias estimable)
#     2. Residual positional correction (lowest residual positional structure)
#     3. Control neutrality (controls behave as expected nulls)
#     4. Functional recovery (LOF variants are separated from controls)
#
#   pct_controls_significant is reported as an informational metric only;
#   it does not gate any filter. Domain scores use rank-based scoring within
#   each targeton rather than z-score standardisation. With only 2 candidates
#   per targeton (D4 reference vs. Plasmid reference), z-scores are unstable;
#   rank scores are robust regardless of the value distribution. Note the
#   positional-correction tier uses a raw-value tie margin
#   (--tie_margin_positional_bias) rather than pure rank scoring, since a
#   rank-based score with only 2 candidates can't distinguish a genuinely
#   tiny gap from the clearest gap in the dataset -- see that flag's help
#   text for why this matters.
#
#   composite_score is retained only as an audit/plotting summary;
#   rank_within_targeton is the decision column.
#
# Usage example (combining the separately-run D4-referenced and Plasmid-
# referenced shrinkage-anchor output directories):
# Rscript qc_signal_noise_metrics_with_composite_V3_anchor.R \
#   --input "GMM_shrinkage_D4/*all_deseq2_results_condition_*.tsv,GMM_shrinkage_Plasmid/*all_deseq2_results_condition_*.tsv" \
#   --out_prefix out/all_targetons_signal_noise \
#   --metric positional \
#   --control "Synonymous_Variant,Intronic_Variant" \
#   --signal "LOF"

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

# ----------------------------
# DEFAULT SETS
# ----------------------------

CONTROL_CONSEQS_DEFAULT  <- c("Synonymous_Variant", "Intronic_Variant")
FUNCTIONAL_CONSEQS_DEFAULT <- c("LOF")

# ----------------------------
# HELPERS
# ----------------------------

parse_csv_arg <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(character(0))
  trimws(unlist(strsplit(x, ",", fixed = TRUE)))
}

guess_ref_from_contrast <- function(contrast_name) {
  if (grepl("_vs_Plasmid", contrast_name)) return("Plasmid")
  if (grepl("_vs_Day4",    contrast_name)) return("Day4")
  return("Other")
}

infer_targeton_name <- function(path) {
  bn <- basename(path)
  bn <- sub("\\.tsv$", "", bn)
  bn <- sub("_missense_all_deseq2_results_condition_.*$", "", bn, ignore.case = TRUE)
  bn <- sub("_all_deseq2_results_condition_.*$",          "", bn, ignore.case = TRUE)
  bn
}

infer_contrast_name <- function(path, dt) {
  if ("contrast" %in% names(dt)) {
    u <- unique(dt$contrast)
    u <- u[!is.na(u)]
    if (length(u) == 1 && nzchar(u)) return(u)
  }
  bn <- basename(path)
  bn <- sub("\\.tsv$", "", bn)
  bn <- sub("^.*_all_deseq2_results_condition_", "", bn)
  bn
}

pick_first_existing <- function(dt, candidates) {
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

# ----------------------------
# COLUMN SELECTION
# ----------------------------

choose_cols <- function(dt, metric) {
  metric <- tolower(metric)

  cfg <- switch(metric,
    pos = , positional = list(
      pfx      = "pos",
      score_fb = "pos_adj_log2FoldChange_raw",
      se_fb    = "pos_total_se_raw",
      label    = "positional"
    ),
    adj = , adjusted = list(
      pfx      = "adj",
      score_fb = "adj_log2FoldChange_raw",
      se_fb    = "lfcSE_raw",
      label    = "adjusted"
    ),
    stop("Unknown --metric. Use 'positional' or 'adjusted'.")
  )

  p <- cfg$pfx
  score_col <- pick_first_existing(dt, c(
    paste0(p, "_combined_LFC"), paste0("combined_LFC_", p), "combined_LFC"
  ))
  se_col <- pick_first_existing(dt, c(
    paste0(p, "_combined_SE"), paste0("combined_SE_", p), "combined_SE"
  ))
  fdr_col <- pick_first_existing(dt, c(
    paste0(p, "_combined_FDR"), paste0("combined_FDR_", p), "combined_FDR",
    paste0(p, "_adj_fdr_raw"), "padj_raw"
  ))
  stat_col <- pick_first_existing(dt, c(
    "anchor_tier", "anchor_call", "GMM_status",
    paste0(p, "_combined_status"), paste0("combined_status_", p), "combined_status",
    paste0("stat_", p, "_raw")
  ))

  if (is.na(score_col)) score_col <- cfg$score_fb
  if (is.na(se_col))    se_col    <- cfg$se_fb

  if (!all(c(score_col, se_col) %in% names(dt))) {
    stop(
      "Metric '", cfg$label, "' requires score/se columns. Expected either ",
      "combined LFC/SE columns or ", cfg$score_fb, " and ", cfg$se_fb, "."
    )
  }
  list(score = score_col, se = se_col, fdr = fdr_col, stat = stat_col)
}

choose_position_col <- function(dt, position_col = NA_character_) {
  if (!is.na(position_col) && nzchar(position_col)) {
    if (!(position_col %in% names(dt))) {
      stop(paste0("Requested --position_col not found: ", position_col))
    }
    return(position_col)
  }
  pick_first_existing(dt, c(
    "genomic_coordinate", "Genomic_Coordinate", "genomic_coord",
    "coordinate", "coord", "position", "pos", "start",
    "start_position", "variant_position"
  ))
}

# ----------------------------
# SAFE SUMMARY STATS
# ----------------------------

safe_median <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) return(NA_real_)
  stats::median(v, na.rm = TRUE)
}

safe_mean <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) return(NA_real_)
  mean(v, na.rm = TRUE)
}

safe_sd <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 2) return(NA_real_)
  stats::sd(v, na.rm = TRUE)
}

safe_mad <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 2) return(NA_real_)
  stats::mad(v, constant = 1.4826, na.rm = TRUE)
}


# ----------------------------
# ADAPTIVE LOESS SPAN
# ----------------------------

# The diagnostic LOESS (#2) estimates residual positional structure that
# remains after the upstream correction LOESS (#1). Its span needs to balance:
#   - Too wide  -> misses genuine residual waves (under-penalises bad correction)
#   - Too narrow -> overfits noise in the residuals (over-penalises good correction)
#
# The adaptive rule widens the span as n decreases, protecting against
# overfitting when few control variants are available. Bounds are clamped to
# [min_span, max_span]. Override with --positional_bias_loess_span to fix the
# span at a constant value instead.

adaptive_loess_span <- function(n,
                                target_effective_points = 30,
                                min_span = 0.30,
                                max_span = 0.90) {
  if (!is.finite(n) || n <= 0) return(max_span)
  raw <- target_effective_points / n
  pmin(pmax(raw, min_span), max_span)
}

# ----------------------------
# RESIDUAL POSITIONAL BIAS
# ----------------------------

calc_residual_positional_bias <- function(x, score_col, position_col, stat_col,
                                          control_set,
                                          loess_span = NA_real_,
                                          min_n = 20) {
  na_result <- function(used = NA_character_, n = 0L) {
    list(
      position_col                             = position_col,
      positional_bias_set_used                 = used,
      n_positional_bias_variants               = n,
      diagnostic_loess_span_used               = NA_real_,
      residual_positional_bias_mean_abs_fitted = NA_real_,
      residual_positional_bias_r2              = NA_real_
    )
  }

  if (is.na(position_col) || !(position_col %in% names(x))) {
    return(na_result())
  }

  y <- copy(x)
  y[, .pos_for_bias   := suppressWarnings(as.numeric(get(position_col)))]
  y[, .score_for_bias := as.numeric(get(score_col))]
  y <- y[is.finite(.pos_for_bias) & is.finite(.score_for_bias)]

  # Prefer "no impact" variants if sufficient, else fall back to controls.
  used <- "controls"
  if (!is.na(stat_col) && stat_col %in% names(y)) {
    y_noimpact <- y[as.character(get(stat_col)) == "no impact"]
    if (nrow(y_noimpact) >= min_n && length(unique(y_noimpact$.pos_for_bias)) >= 5) {
      y    <- y_noimpact
      used <- "noimpact"
    } else {
      y <- y[consequence %in% control_set]
    }
  } else {
    y <- y[consequence %in% control_set]
  }

  n <- nrow(y)
  if (n < min_n || length(unique(y$.pos_for_bias)) < 5) {
    return(na_result(used, n))
  }

  span_used <- if (is.finite(loess_span)) loess_span else adaptive_loess_span(n)

  fit <- try(
    stats::loess(
      .score_for_bias ~ .pos_for_bias,
      data      = y,
      span      = span_used,
      degree    = 1,
      surface   = "direct",
      na.action = stats::na.exclude
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) return(na_result(used, n))

  fitted_vals <- suppressWarnings(as.numeric(stats::predict(fit, newdata = y)))
  fitted_vals <- fitted_vals[is.finite(fitted_vals)]
  if (length(fitted_vals) < min_n) return(na_result(used, n))

  fitted_centered <- fitted_vals - stats::median(fitted_vals, na.rm = TRUE)

  rss <- sum(stats::residuals(fit)^2, na.rm = TRUE)
  tss <- sum((y$.score_for_bias - mean(y$.score_for_bias, na.rm = TRUE))^2,
             na.rm = TRUE)
  r2  <- ifelse(is.finite(tss) && tss > 0, 1 - rss / tss, NA_real_)

  list(
    position_col                             = position_col,
    positional_bias_set_used                 = used,
    n_positional_bias_variants               = n,
    diagnostic_loess_span_used               = span_used,
    residual_positional_bias_mean_abs_fitted = safe_mean(abs(fitted_centered)),
    residual_positional_bias_r2              = r2
  )
}

# ----------------------------
# COHEN'S D
# ----------------------------

cohens_d <- function(x_signal, x_control) {
  x_signal  <- x_signal[is.finite(x_signal)]
  x_control <- x_control[is.finite(x_control)]
  if (length(x_signal) < 2 || length(x_control) < 2) return(NA_real_)

  s1 <- stats::var(x_signal)
  s0 <- stats::var(x_control)
  sp <- sqrt(
    ((length(x_signal)  - 1) * s1 +
     (length(x_control) - 1) * s0) /
    (length(x_signal) + length(x_control) - 2)
  )
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(x_signal) - mean(x_control)) / sp
}

# ----------------------------
# PER-FILE SUMMARY
# ----------------------------

summarise_contrast <- function(dt, contrast_name, control_set, signal_set,
                               fdr_cutoff                = 0.05,
                               metric                    = "positional",
                               min_n_for_d               = 10,
                               cohens_d_on_lfc           = FALSE,
                               position_col              = NA_character_,
                               positional_bias_loess_span = NA_real_,
                               min_n_for_positional_bias = 20) {

  if (!("consequence" %in% names(dt))) {
    stop("Input must contain a 'consequence' column.")
  }

  ref  <- guess_ref_from_contrast(contrast_name)
  cols <- choose_cols(dt, metric)

  score_col <- cols$score
  se_col    <- cols$se
  fdr_col   <- cols$fdr
  stat_col  <- cols$stat
  pos_col   <- choose_position_col(dt, position_col)

  x <- as.data.table(copy(dt))
  x <- x[!is.na(get(score_col)) & !is.na(get(se_col)) & get(se_col) > 0]
  x[, z     := get(score_col) / get(se_col)]
  x[, abs_z := abs(z)]

  median_abs_z <- safe_median(x$abs_z)

  xc0 <- x[consequence %in% CONTROL_CONSEQS_DEFAULT]
  xf0 <- x[consequence %in% FUNCTIONAL_CONSEQS_DEFAULT]

  median_abs_z_controls   <- safe_median(xc0$abs_z)
  median_abs_z_functional <- safe_median(xf0$abs_z)
  n_controls_median       <- nrow(xc0)
  n_functional_median     <- nrow(xf0)

  # --- Control neutrality ---
  n_controls               <- 0L
  pct_controls_significant <- NA_real_
  control_sig_basis        <- NA_character_
  control_scores           <- numeric(0)
  control_z                <- numeric(0)

  if (length(control_set) > 0) {
    xc             <- x[consequence %in% control_set]
    n_controls     <- nrow(xc)
    control_scores <- xc[[score_col]]
    control_z      <- xc$z

    if (n_controls > 0) {
      if (!is.na(stat_col) && stat_col %in% names(xc)) {
        pct_controls_significant <-
          100 * mean(as.character(xc[[stat_col]]) != "no impact", na.rm = TRUE)
        control_sig_basis <- paste0("stat:", stat_col)
      } else if (!is.na(fdr_col) && fdr_col %in% names(xc)) {
        pct_controls_significant <-
          100 * mean(xc[[fdr_col]] < fdr_cutoff, na.rm = TRUE)
        control_sig_basis <- paste0("fdr:", fdr_col, "<", fdr_cutoff)
      }
    }
  }

  # --- Cohen's d (signal vs control) ---
  # By default computed on z-scores, which rewards both large effect size and
  # precision. Set cohens_d_on_lfc = TRUE to use raw LFCs instead, which
  # measures pure effect-size separation independent of SE.
  d        <- NA_real_
  n_signal <- 0L
  n_ctrl_d <- 0L

  if (length(signal_set) > 0 && length(control_set) > 0) {
    xs       <- x[consequence %in% signal_set]
    xc       <- x[consequence %in% control_set]
    n_signal <- nrow(xs)
    n_ctrl_d <- nrow(xc)

    if (n_signal >= min_n_for_d && n_ctrl_d >= min_n_for_d) {
      if (cohens_d_on_lfc) {
        d <- cohens_d(xs[[score_col]], xc[[score_col]])
      } else {
        d <- cohens_d(xs$z, xc$z)
      }
    }
  }

  # --- Residual positional bias (diagnostic LOESS #2) ---
  pos_bias <- calc_residual_positional_bias(
    x           = x,
    score_col   = score_col,
    position_col = pos_col,
    stat_col    = stat_col,
    control_set = control_set,
    loess_span  = positional_bias_loess_span,
    min_n       = min_n_for_positional_bias
  )

  data.table(
    contrast        = contrast_name,
    ref             = ref,
    input_n         = nrow(dt),
    n_variants_used = nrow(x),
    metric          = metric,
    score_col       = score_col,
    se_col          = se_col,
    fdr_col         = ifelse(is.na(fdr_col),  NA_character_, fdr_col),
    stat_col        = ifelse(is.na(stat_col), NA_character_, stat_col),
    position_col    = ifelse(is.na(pos_bias$position_col), NA_character_,
                             pos_bias$position_col),

    median_abs_z = median_abs_z,

    control_consequence_set_for_median    = paste(CONTROL_CONSEQS_DEFAULT, collapse = ","),
    n_controls_for_median                 = n_controls_median,
    median_abs_z_controls                 = median_abs_z_controls,

    functional_consequence_set_for_median = paste(FUNCTIONAL_CONSEQS_DEFAULT, collapse = ","),
    n_functional_for_median               = n_functional_median,
    median_abs_z_functional               = median_abs_z_functional,

    control_set              = paste(control_set, collapse = ","),
    n_controls               = n_controls,
    mad_control_score        = safe_mad(control_scores),
    median_abs_control_score = safe_median(abs(control_scores)),
    mean_control_z           = safe_mean(control_z),
    sd_control_z             = safe_sd(control_z),

    positional_bias_set_used                 = pos_bias$positional_bias_set_used,
    n_positional_bias_variants               = pos_bias$n_positional_bias_variants,
    diagnostic_loess_span_used               = pos_bias$diagnostic_loess_span_used,
    residual_positional_bias_mean_abs_fitted = pos_bias$residual_positional_bias_mean_abs_fitted,
    residual_positional_bias_r2              = pos_bias$residual_positional_bias_r2,

    pct_controls_significant = pct_controls_significant,
    control_sig_basis        = control_sig_basis,

    signal_set               = paste(signal_set, collapse = ","),
    cohens_d_computed_on     = ifelse(cohens_d_on_lfc, "lfc", "z_score"),
    n_signal                 = n_signal,
    n_control_for_d          = n_ctrl_d,
    cohens_d_z_signal_vs_control = d
  )
}

# ----------------------------
# SANITY CHECK: CANDIDATE COUNT
# ----------------------------

# Warns when a targeton does not have the expected number of candidates.
# Mismatches almost always indicate a filename parsing failure silently
# corrupting the targeton grouping.

check_candidate_counts <- function(summary_dt, expected_n = 4L) {
  counts <- summary_dt[, .N, by = targeton]
  bad    <- counts[N != expected_n]

  if (nrow(bad) > 0) {
    warning(
      "The following targetons do not have exactly ", expected_n,
      " candidates. This may indicate a filename parsing error:\n",
      paste0("  ", bad$targeton, " (found ", bad$N, ")", collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(counts)
}

# ----------------------------
# MIXED SIGNIFICANCE METHOD FLAG
# ----------------------------

# Flags targetons where pct_controls_significant is computed using different
# methods (stat column vs FDR fallback) across candidates. Since
# pct_controls_significant is now informational only, this flag is a
# data-quality note rather than a filter concern.

flag_mixed_sig_method <- function(summary_dt) {
  x <- copy(summary_dt)

  x[, control_sig_method_mixed := {
    methods <- unique(na.omit(control_sig_basis))
    length(methods) > 1
  }, by = targeton]

  mixed_targetons <- unique(x[control_sig_method_mixed == TRUE, targeton])
  if (length(mixed_targetons) > 0) {
    warning(
      "The following targetons have pct_controls_significant computed using ",
      "mixed methods (stat column for some candidates, FDR for others). ",
      "Comparisons of pct_controls_significant across candidates may be ",
      "inconsistent:\n",
      paste0("  ", mixed_targetons, collapse = "\n"),
      call. = FALSE
    )
  }

  x[]
}

# ----------------------------
# RANK-BASED DOMAIN SCORING
# ----------------------------

# Replaces z-score standardisation with rank-based scoring within each
# targeton. With typically only 4 candidates per targeton, z-score
# standardisation is unstable: one outlier compresses all other scores toward
# zero. Rank scores are robust regardless of the value distribution.
#
# Returns values in [0, 1]: 1 = best candidate in the group, 0 = worst.
# Ties are averaged (ties.method = "average"). NAs receive NA.

rank_score_by_group <- function(x, higher_is_better = TRUE) {
  x   <- as.numeric(x)
  ok  <- is.finite(x)
  out <- rep(NA_real_, length(x))
  n   <- sum(ok)

  if (n == 0L) return(out)
  if (n == 1L) { out[ok] <- 1; return(out) }

  r <- rank(x[ok], ties.method = "average")
  if (higher_is_better) {
    out[ok] <- (r - 1) / (n - 1)
  } else {
    out[ok] <- (n - r) / (n - 1)
  }
  out
}

# ----------------------------
# CASCADING, MAGNITUDE-AWARE TIE-BREAKING
# ----------------------------

# The plain hierarchical setorder() previously used has a blind spot with
# few candidates per targeton (as few as 2): rank-based domain scores
# convert ANY nonzero raw-metric gap into the same maximal 0-vs-1 score
# difference, so a genuinely tiny (likely noise-level) top-priority gap
# gets treated as equally decisive as the clearest gap seen anywhere in
# the dataset -- and because it's checked first, it can override a much
# larger, clearly meaningful gap on a lower-priority tier.
#
# This function fixes that: at each tier, candidates are only separated if
# their gap on that tier's VALUE COLUMN exceeds that tier's margin; within
# a margin, they're treated as tied and the decision cascades to the next
# tier. margin = 0 preserves the original fully-decisive behaviour for that
# tier (only exact ties defer further) -- so tiers without an explicit
# magnitude gate behave exactly as before; only the tier(s) given a
# nonzero margin gain magnitude-awareness.
#
# criteria: list of list(value_col = <column, already oriented so higher
# is better>, margin = <numeric, 0 disables magnitude-gating>)
# Returns an integer rank (1 = best) aligned to dt_group's row order.
cascading_tier_rank <- function(dt_group, criteria) {
  n <- nrow(dt_group)
  if (n == 0L) return(integer(0))
  if (n == 1L) return(1L)

  assign_order <- function(idx, tier_i) {
    if (length(idx) <= 1 || tier_i > length(criteria)) return(idx)
    crit   <- criteria[[tier_i]]
    vals   <- dt_group[[crit$value_col]][idx]
    margin <- if (is.na(crit$margin)) 0 else crit$margin

    ord        <- order(vals, decreasing = TRUE, na.last = TRUE)
    idx_sorted <- idx[ord]
    vals_sorted <- vals[ord]

    groups   <- list()
    cur      <- idx_sorted[1]
    cur_best <- vals_sorted[1]
    if (length(idx_sorted) > 1) {
      for (k in 2:length(idx_sorted)) {
        gap <- cur_best - vals_sorted[k]
        if (is.na(gap) || gap <= margin) {
          cur <- c(cur, idx_sorted[k])
        } else {
          groups[[length(groups) + 1]] <- cur
          cur      <- idx_sorted[k]
          cur_best <- vals_sorted[k]
        }
      }
    }
    groups[[length(groups) + 1]] <- cur

    unlist(lapply(groups, assign_order, tier_i = tier_i + 1))
  }

  final_order <- assign_order(seq_len(n), 1)
  rank_out <- integer(n)
  rank_out[final_order] <- seq_len(n)
  rank_out
}

# ----------------------------
# HIERARCHICAL SCORING & RANKING
# ----------------------------

add_hierarchical_scores <- function(summary_dt,
                                    require_controls        = TRUE,
                                    require_positional_bias = TRUE,
                                    near_tie_margin         = 0.10,
                                    tie_margin_positional_bias = 0.02) {
  x <- copy(summary_dt)

  # --- Hard/pass-fail filters ---
  # pct_controls_significant is informational only; not a filter.
  x[, pass_control_n       := !require_controls | n_controls > 0]
  x[, pass_positional_bias := !require_positional_bias |
        is.finite(residual_positional_bias_mean_abs_fitted)]

  x[, passes_basic_filters              := pass_control_n & pass_positional_bias]
  x[, all_conditions_fail_basic_filters := !any(passes_basic_filters), by = targeton]

  # --- Rank-based domain components (within each targeton) ---
  # Each component is in [0, 1]; higher = better after direction is applied.
  x[, rc_control_mad :=
      rank_score_by_group(mad_control_score, higher_is_better = FALSE),
    by = targeton]

  x[, rc_control_median_abs :=
      rank_score_by_group(median_abs_control_score, higher_is_better = FALSE),
    by = targeton]

  x[, rc_residual_positional_bias :=
      rank_score_by_group(residual_positional_bias_mean_abs_fitted,
                          higher_is_better = FALSE),
    by = targeton]

  x[, rc_functional_abs_z :=
      rank_score_by_group(median_abs_z_functional, higher_is_better = TRUE),
    by = targeton]

  x[, rc_cohens_d :=
      rank_score_by_group(abs(cohens_d_z_signal_vs_control), higher_is_better = TRUE),
    by = targeton]

  # --- Domain scores (means of rank components) ---
  control_cols    <- c("rc_control_mad", "rc_control_median_abs")
  positional_cols <- c("rc_residual_positional_bias")
  functional_cols <- c("rc_functional_abs_z", "rc_cohens_d")

  x[, control_neutrality_score    := rowMeans(.SD, na.rm = TRUE), .SDcols = control_cols]
  x[, positional_correction_score := rowMeans(.SD, na.rm = TRUE), .SDcols = positional_cols]
  x[, functional_recovery_score   := rowMeans(.SD, na.rm = TRUE), .SDcols = functional_cols]

  x[rowSums(!is.na(x[, ..control_cols]))    == 0, control_neutrality_score    := NA_real_]
  x[rowSums(!is.na(x[, ..positional_cols])) == 0, positional_correction_score := NA_real_]
  x[rowSums(!is.na(x[, ..functional_cols])) == 0, functional_recovery_score   := NA_real_]

  # --- Hierarchical ranking ---
  # Priority: pass filters > positional correction > control neutrality >
  #           functional recovery.
  # Rationale: once controls are below the false-positive threshold, the main
  # purpose of comparing reference conditions is minimising residual positional
  # structure. Control neutrality and LOF recovery break ties between
  # similarly flat fits.
  #
  # Positional correction gets a MAGNITUDE-AWARE tie margin
  # (tie_margin_positional_bias, on the RAW residual_positional_bias_
  # mean_abs_fitted scale, not the rank score) because with few candidates
  # per targeton, a rank-based score cannot distinguish "the clearest gap
  # in the dataset" from "a hair's-breadth, likely-noise gap" -- both
  # become the same maximal score difference. Without this gate, a tiny
  # top-tier win can override a large, clearly meaningful control-
  # neutrality difference (this is exactly what happened for one targeton
  # with only a 0.008 raw positional-bias gap, the smallest in the
  # dataset, while every other targeton's gap was >= 0.0096 and typically
  # 0.03-0.09). Lower tiers keep margin = 0 (fully decisive, as before)
  # since there's no evidence they need magnitude-gating too -- extend if
  # a similar edge case turns up there.
  x[, ranking_system := paste0(
    "hierarchical (magnitude-aware cascade): pass_filters > ",
    "positional_correction (raw-value tie margin=", tie_margin_positional_bias,
    ") > control_neutrality > functional_recovery | ",
    "scoring: rank_based [0,1] within targeton"
  )]

  x[, rank_within_targeton := NA_integer_]

  criteria <- list(
    list(value_col = "neg_residual_positional_bias", margin = tie_margin_positional_bias),
    list(value_col = "control_neutrality_score",      margin = 0),
    list(value_col = "functional_recovery_score",      margin = 0)
  )
  # Oriented so higher = better, matching cascading_tier_rank's convention
  # (residual bias itself is lower-is-better in raw units).
  x[, neg_residual_positional_bias := -residual_positional_bias_mean_abs_fitted]

  for (t in unique(x$targeton)) {
    pass_idx <- which(x$targeton == t & x$passes_basic_filters)
    fail_idx <- which(x$targeton == t & !x$passes_basic_filters)
    if (length(pass_idx) > 0) {
      ranks <- cascading_tier_rank(x[pass_idx], criteria)
      x[pass_idx, rank_within_targeton := ranks]
    }
    if (length(fail_idx) > 0) {
      # Failing basic filters always ranks after every passing candidate;
      # order among themselves by composite_score (informational only,
      # since none of them should be selected anyway).
      fail_order <- order(-x$composite_score[fail_idx], na.last = TRUE)
      x[fail_idx[fail_order], rank_within_targeton := length(pass_idx) + seq_along(fail_idx)]
    }
  }
  x[, neg_residual_positional_bias := NULL]

  # --- Audit summary score (for plotting/inspection only; not used for ranking) ---
  domain_cols <- c("control_neutrality_score", "positional_correction_score",
                   "functional_recovery_score")
  x[, composite_score := rowMeans(.SD, na.rm = TRUE), .SDcols = domain_cols]
  x[rowSums(!is.na(x[, ..domain_cols])) == 0, composite_score := NA_real_]

  # Near-tie flag: audit-score gap to the top-ranked candidate. Interpreted as
  # a prompt for manual inspection, not as a ranking rule.
  x[, best_composite_score        := composite_score[rank_within_targeton == 1], by = targeton]
  x[, composite_score_gap_to_best := best_composite_score - composite_score]
  x[, near_tie_with_best          := is.finite(composite_score_gap_to_best) &
        composite_score_gap_to_best <= near_tie_margin]
  x[rank_within_targeton == 1, near_tie_with_best := FALSE]

  x[]
}

# ----------------------------
# OUTPUT TABLE
# ----------------------------

make_output_table <- function(summary_dt, verbose_output = FALSE) {
  if (isTRUE(verbose_output)) return(summary_dt[])

  compact_cols <- c(
    "targeton",
    "contrast",
    "ref",
    "n_controls",
    "n_signal",
    "control_sig_method_mixed",
    "pct_controls_significant",   # informational; not used as a filter
    "mad_control_score",
    "median_abs_control_score",
    "diagnostic_loess_span_used",
    "residual_positional_bias_mean_abs_fitted",
    "median_abs_z_functional",
    "cohens_d_computed_on",
    "cohens_d_z_signal_vs_control",
    "control_neutrality_score",
    "positional_correction_score",
    "functional_recovery_score",
    "passes_basic_filters",
    "all_conditions_fail_basic_filters",
    "ranking_system",
    "composite_score",
    "composite_score_gap_to_best",
    "near_tie_with_best",
    "rank_within_targeton"
  )

  compact_cols <- compact_cols[compact_cols %in% names(summary_dt)]
  summary_dt[, ..compact_cols]
}

# ----------------------------
# CLI
# ----------------------------

opt_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = paste0(
      "Input TSV(s). Can be a glob like 'path/*all_deseq2_results_condition_*.tsv'. ",
      "Accepts multiple comma-separated globs to combine files from separate ",
      "directories in one run, e.g. results living in separate per-reference-",
      "condition output folders: 'D4_results/*.tsv,Plasmid_results/*.tsv'"
    )
  ),
  make_option(
    c("-o", "--out_prefix"),
    type = "character",
    default = "signal_noise",
    help = "Output prefix [default %default]"
  ),
  make_option(
    c("--metric"),
    type = "character",
    default = "positional",
    help = "Metric set to use: positional or adjusted [default %default]"
  ),
  make_option(
    c("--control"),
    type = "character",
    default = "Synonymous_Variant,Intronic_Variant",
    help = "Comma-separated neutral/control consequence labels [default %default]"
  ),
  make_option(
    c("--signal"),
    type = "character",
    default = "LOF",
    help = "Comma-separated signal consequence labels used for Cohen's d [default %default]"
  ),
  make_option(
    c("--fdr_cutoff"),
    type = "double",
    default = 0.05,
    help = "FDR cutoff for fallback significance calling [default %default]"
  ),
  make_option(
    c("--min_n_for_d"),
    type = "integer",
    default = 10,
    help = "Minimum n per group for Cohen's d [default %default]"
  ),
  make_option(
    c("--cohens_d_on_lfc"),
    action = "store_true",
    default = FALSE,
    help = paste0(
      "Compute Cohen's d on raw LFCs instead of z-scores. ",
      "LFC-based d measures pure effect-size separation independent of SE; ",
      "z-score-based d (default) rewards both effect size and precision."
    )
  ),
  make_option(
    c("--position_col"),
    type = "character",
    default = "",
    help = "Genomic coordinate column. If blank, common names are auto-detected [default: auto]"
  ),
  make_option(
    c("--positional_bias_loess_span"),
    type = "double",
    default = NA_real_,
    help = paste0(
      "Span for the diagnostic LOESS used to measure residual positional bias. ",
      "If omitted (default), an adaptive span is used: wider when few control ",
      "variants are available (reducing noise overfitting), narrower when many ",
      "are available (resolving genuine residual waves). ",
      "Provide a fixed value in (0, 1] to override."
    )
  ),
  make_option(
    c("--min_n_for_positional_bias"),
    type = "integer",
    default = 20,
    help = "Minimum variants required to estimate residual positional bias [default %default]"
  ),
  make_option(
    c("--expected_candidates_per_targeton"),
    type = "integer",
    default = 2L,
    help = paste0(
      "Expected number of candidate conditions per targeton ",
      "(default 2 = D4 reference vs. Plasmid reference). ",
      "A warning is emitted for any targeton with a different count, ",
      "which usually indicates a filename parsing error [default %default]"
    )
  ),
  make_option(
    c("--allow_no_controls"),
    action = "store_true",
    default = FALSE,
    help = "Allow candidates with zero control variants to pass basic filters [default %default]"
  ),
  make_option(
    c("--allow_no_positional_bias"),
    action = "store_true",
    default = FALSE,
    help = paste0(
      "Allow candidates with missing residual positional-bias metrics to ",
      "pass basic filters [default %default]"
    )
  ),
  make_option(
    c("--near_tie_margin"),
    type = "double",
    default = 0.10,
    help = paste0(
      "Audit-score gap to best below which a non-winning candidate is ",
      "flagged as a near tie warranting manual inspection [default %default]"
    )
  ),
  make_option(
    c("--tie_margin_positional_bias"),
    type = "double",
    default = 0.02,
    help = paste0(
      "Raw-value tie margin (units of residual_positional_bias_mean_abs_fitted) ",
      "below which two candidates are treated as tied on positional correction, ",
      "deferring the ranking decision to control neutrality instead of letting ",
      "a possibly noise-level gap decide outright. Default chosen from observed ",
      "gap sizes across a real 22-targeton run: the smallest genuine (correctly ",
      "decisive) gap was 0.0096, the one problem case was 0.008 -- 0.02 sits ",
      "safely between noise-level and real gaps for that dataset, but check ",
      "your own gap distribution and adjust if needed [default %default]"
    )
  ),
  make_option(
    c("--verbose_output"),
    action = "store_true",
    default = FALSE,
    help = paste0(
      "Write all intermediate QC metrics and rank components instead of ",
      "the compact default table [default %default]"
    )
  )
)

opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$input) || !nzchar(opt$input)) stop("Please provide --input")
if (is.null(opt$out_prefix) || !nzchar(opt$out_prefix)) stop("Please provide --out_prefix")

# out_prefix may include a directory path (e.g. "thesis/GMM_shrinkage/
# signal_to_noise/all_targetons_signal_noise") that doesn't exist yet --
# fwrite/ggsave won't create it, so create it up front rather than failing
# after all the analysis work is done.
out_dir_for_prefix <- dirname(opt$out_prefix)
if (nzchar(out_dir_for_prefix) && out_dir_for_prefix != "." && !dir.exists(out_dir_for_prefix)) {
  dir.create(out_dir_for_prefix, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir_for_prefix)) {
    stop("Could not create output directory: ", out_dir_for_prefix,
         " -- check permissions and disk space.")
  }
  message("Created output directory: ", out_dir_for_prefix)
}

# --input accepts one or more comma-separated glob patterns, so results
# living in separate directories (e.g. separate D4-referenced and
# Plasmid-referenced shrinkage-anchor output folders, run separately per
# the shrinkage pipeline's own duplicate-targeton-name guard) can be
# combined into one comparison run.
input_patterns <- parse_csv_arg(opt$input)
files <- unique(unlist(lapply(input_patterns, Sys.glob)))
if (length(files) == 0) stop(paste0("No files matched: ", opt$input))

control_set <- parse_csv_arg(opt$control)
signal_set  <- parse_csv_arg(opt$signal)

# NA_real_ is the sentinel for "use adaptive span"
loess_span_arg <- if (is.na(opt$positional_bias_loess_span)) {
  NA_real_
} else {
  as.double(opt$positional_bias_loess_span)
}

all_summ <- vector("list", length(files))

for (i in seq_along(files)) {
  f  <- files[[i]]
  dt <- fread(f)
  cn <- infer_contrast_name(f, dt)

  s <- summarise_contrast(
    dt                        = dt,
    contrast_name             = cn,
    control_set               = control_set,
    signal_set                = signal_set,
    fdr_cutoff                = opt$fdr_cutoff,
    metric                    = opt$metric,
    min_n_for_d               = opt$min_n_for_d,
    cohens_d_on_lfc           = opt$cohens_d_on_lfc,
    position_col              = opt$position_col,
    positional_bias_loess_span = loess_span_arg,
    min_n_for_positional_bias = opt$min_n_for_positional_bias
  )

  s[, targeton              := infer_targeton_name(f)]
  s[, input_file            := f]
  all_summ[[i]] <- s
}

summary_dt <- rbindlist(all_summ, fill = TRUE)

setcolorder(
  summary_dt,
  c("targeton", "contrast", "ref",
    setdiff(names(summary_dt),
            c("targeton", "contrast", "ref")))
)
setorder(summary_dt, targeton, ref, contrast)

# Sanity check: warn if any targeton has an unexpected number of candidates.
check_candidate_counts(summary_dt,
                       expected_n = opt$expected_candidates_per_targeton)

# Flag targetons where significance method differs across candidates.
summary_dt <- flag_mixed_sig_method(summary_dt)

summary_dt <- add_hierarchical_scores(
  summary_dt,
  require_controls        = !opt$allow_no_controls,
  require_positional_bias = !opt$allow_no_positional_bias,
  near_tie_margin         = opt$near_tie_margin,
  tie_margin_positional_bias = opt$tie_margin_positional_bias
)

output_dt <- make_output_table(summary_dt, verbose_output = opt$verbose_output)

out_metrics <- paste0(opt$out_prefix, ".metrics.tsv")
fwrite(output_dt, out_metrics, sep = "\t")

if (isTRUE(opt$verbose_output)) {
  message("Wrote verbose metrics table: ", out_metrics)
} else {
  message("Wrote compact metrics table: ", out_metrics)
}

# ----------------------------
# PLOTS
# ----------------------------

# Both plots use a fixed-column facet_wrap grid rather than facet_grid with
# one row per targeton -- height then scales with the number of GRID ROWS
# (ceiling(n_targetons / n_facet_cols)), not the raw targeton count, so it
# stays well within ggsave's 50-inch limit regardless of how many targetons
# are in a run (facet_grid with 22 targetons at 2.5in/row hit 55in and
# crashed ggsave; this bounds height sensibly instead).
n_facet_cols <- 5

# facet_wrap orders panels by the faceting variable's factor levels, which
# default to alphabetical for a character column -- and since the exon
# number is embedded as a string suffix (e.g. "APDY_exon2" vs
# "AXFG_exon19"), alphabetical order puts exon2 after exon19. Parse the
# exon number out of each targeton name and use it to set explicit factor
# levels instead, so panels appear in exon/genomic order. Falls back to
# alphabetical order (after all exon-numbered ones) for any targeton name
# without a parseable exon number, with a warning.
extract_exon_number <- function(targeton_name) {
  has_match <- grepl("exon[_]?[0-9]+", targeton_name, ignore.case = TRUE)
  num <- rep(NA_integer_, length(targeton_name))
  num[has_match] <- as.integer(
    sub(".*exon[_]?([0-9]+).*", "\\1", targeton_name[has_match], ignore.case = TRUE)
  )
  num
}

targeton_levels <- sort(unique(summary_dt$targeton))
exon_nums <- extract_exon_number(targeton_levels)
if (any(is.na(exon_nums))) {
  warning(
    "Could not parse an exon number from the following targeton name(s); ",
    "these will be sorted alphabetically after the exon-numbered ones in ",
    "plot facet order: ",
    paste(targeton_levels[is.na(exon_nums)], collapse = ", "),
    call. = FALSE
  )
}
targeton_levels_by_exon <- targeton_levels[order(is.na(exon_nums), exon_nums, targeton_levels)]

plot_dt <- summary_dt[ref %in% c("Plasmid", "Day4") & is.finite(median_abs_z)]

if (nrow(plot_dt) > 0) {
  plot_dt[, targeton := factor(targeton, levels = targeton_levels_by_exon)]
  n_rows_1 <- ceiling(length(unique(plot_dt$targeton)) / n_facet_cols)

  p1 <- ggplot(
    plot_dt,
    aes(x = ref, y = median_abs_z, shape = ref)
  ) +
    geom_point(size = 2.5) +
    facet_wrap(~ targeton, scales = "free_x", ncol = n_facet_cols) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste0("Signal-to-noise summary: median |Z| (metric=", opt$metric, ")"),
      x     = "Reference",
      y     = "median(|Z|)",
      shape = "Reference"
    )

  out_png_1 <- paste0(opt$out_prefix, ".median_abs_z.png")
  ggsave(
    filename = out_png_1, plot = p1,
    width  = min(45, max(10, 3.0 * n_facet_cols)),
    height = min(45, max(5,  2.5 * n_rows_1)),
    units  = "in", dpi = 300
  )
  message("Wrote: ", out_png_1)
}

# median_abs_z_functional (LOF-restricted median |Z|) -- unlike the
# all-variant median_abs_z above, this isolates whether the true-signal
# class specifically shows more separation, rather than reflecting a
# generalised noise shift across every variant type. This is the more
# appropriate metric for judging "did an alternative reference pull out
# more true depletion signal", and should be read alongside
# cohens_d_z_signal_vs_control (which additionally accounts for whether
# any apparent gain here is just noise inflation rather than real added
# separation) rather than in isolation.
plot_dt_func <- summary_dt[ref %in% c("Plasmid", "Day4") & is.finite(median_abs_z_functional)]

if (nrow(plot_dt_func) > 0) {
  plot_dt_func[, targeton := factor(targeton, levels = targeton_levels_by_exon)]
  n_rows_1b <- ceiling(length(unique(plot_dt_func$targeton)) / n_facet_cols)

  p1b <- ggplot(
    plot_dt_func,
    aes(x = ref, y = median_abs_z_functional, shape = ref)
  ) +
    geom_point(size = 2.5) +
    facet_wrap(~ targeton, scales = "free_x", ncol = n_facet_cols) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste0("Signal-to-noise summary: median |Z|, LOF variants only (metric=", opt$metric, ")"),
      x     = "Reference",
      y     = "median(|Z|), LOF variants only",
      shape = "Reference"
    )

  out_png_1b <- paste0(opt$out_prefix, ".median_abs_z_functional.png")
  ggsave(
    filename = out_png_1b, plot = p1b,
    width  = min(45, max(10, 3.0 * n_facet_cols)),
    height = min(45, max(5,  2.5 * n_rows_1b)),
    units  = "in", dpi = 300
  )
  message("Wrote: ", out_png_1b)
}

plot_ranked <- summary_dt[ref %in% c("Plasmid", "Day4") & is.finite(composite_score)]

if (nrow(plot_ranked) > 0) {
  plot_ranked[, targeton := factor(targeton, levels = targeton_levels_by_exon)]
  n_rows_2 <- ceiling(length(unique(plot_ranked$targeton)) / n_facet_cols)

  p2 <- ggplot(
    plot_ranked,
    aes(x = ref, y = composite_score, shape = passes_basic_filters)
  ) +
    geom_point(size = 2.5) +
    facet_wrap(~ targeton, scales = "free_x", ncol = n_facet_cols) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Audit summary score by candidate analysis condition",
      x     = "Reference",
      y     = "Audit summary score; ranking uses rank_within_targeton",
      shape = "Passed basic filters"
    )

  out_png_2 <- paste0(opt$out_prefix, ".composite_score.png")
  ggsave(
    filename = out_png_2, plot = p2,
    width  = min(45, max(10, 3.0 * n_facet_cols)),
    height = min(45, max(5,  2.5 * n_rows_2)),
    units  = "in", dpi = 300
  )
  message("Wrote: ", out_png_2)
}
