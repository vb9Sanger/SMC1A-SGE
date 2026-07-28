#!/usr/bin/env Rscript

# =============================================================================
# run_maveqc_VB_fixed.R
#
# Author: Vanessa Burns 
# Version: 1.1
# Date: 2026-07-17
#
# Description:
# End-to-end analysis pipeline for saturation genome editing (SGE) screens.
# Performs read processing, positional correction, normalisation, log fold-
# change calculation, quality control, and generation of analysis outputs,
# INCLUDING the final MAVEQC_report.html for plasmid and screen QC.
#
# =============================================================================

# Parse command line arguments (only enforced when this script is Rscript'd
# directly with two args, e.g. the original full pipeline run; safe to
# source() this file for just the function definitions, e.g. to call
# create_qc_reports() on existing outputs, without supplying these).
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 2) {
  input_dir  <- args[1]
  output_dir <- args[2]
} else {
  input_dir  <- NULL
  output_dir <- NULL
}

# Create user library directory if it doesn't exist
user_lib <- Sys.getenv("R_LIBS_USER", unset = "~/R/library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(user_lib)

# Global variables
maveqc_config <- NULL
maveqc_ref_time_point <- NULL
maveqc_ref_time_point_samples <- NULL
maveqc_deseq_coldata <- NULL
qc_samplesheet <- NULL

# Function to safely install packages
install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        try(install.packages(pkg, 
                           lib = user_lib,
                           repos = "https://cloud.r-project.org",
                           quiet = TRUE))
    }
}

# Install required packages if not already installed
required_packages <- c(
  "jsonlite", "configr", "vroom", "data.table", "Ckmeans.1d.dp", 
  "gplots", "ggplot2", "plotly", "ggcorrplot", "corrplot", "see",
  "ggbeeswarm", "reactable", "reshape2", "htmltools", "sparkline",
  "dendextend", "gtools", "BiocManager"
)

# Try to install each package
for (pkg in required_packages) {
    install_if_missing(pkg)
}

# Install Bioconductor packages if not already installed
if (requireNamespace("BiocManager", quietly = TRUE)) {
    BiocManager::install(c("DESeq2", "DEGreport", "apeglm"),
                        lib = user_lib,
                        quiet = TRUE,
                        ask = FALSE)
}

# Load required libraries
library(configr)
library(vroom)
library(data.table)
library(Ckmeans.1d.dp)
library(gplots)
library(ggplot2)
library(plotly)
library(ggcorrplot)
library(corrplot)
library(see)
library(ggbeeswarm)
library(reactable)
library(htmltools)
library(sparkline)
library(dendextend)
library(reshape2)
library(gtools)
library(DESeq2)
library(DEGreport)
library(apeglm)

# - classes.R
#' A class representing a SGE object
#'
#' @export
#' @name SGE
#' @slot sample             the sample name
#' @slot libname            library name
#' @slot libtype            library type
#' @slot adapt5             adaptor sequence at 5 prime end
#' @slot adapt3             adaptor sequence at 3 prime end
#' @slot per_r1_adaptor     percentage of r1 adatpor in the sample sheet
#' @slot per_r2_adaptor     percentage of r2 adatpor in the sample sheet
#' @slot refseq             reference sequence
#' @slot pamseq             sequence with pam variants
#' @slot libcounts          QUANTS library-dependent counts, per sequence per count
#' @slot allcounts          QUANTS library-independent counts, per sequence per count
#' @slot valiant_meta       VaLiAnT meta file
#' @slot vep_anno           vep consequence annotation file
#' @slot meta_mseqs         non-redundant mseq in VaLiAnT meta file
#' @slot missing_meta_seqs  missing sequenced in library compared to VaLiAnT meta file
#' @slot libstats           summaries of library dependent counts
#' @slot allstats           summaries of library independent counts
#' @slot libstats_qc        qc stats of library dependent counts
#' @slot allstats_qc        qc stats of library independent counts

setClass("SGE",
    slots = list(
        sample = "character",
        libname = "character",
        libtype = "character",
        adapt5 = "character",
        adapt3 = "character",
        per_r1_adaptor = "numeric",
        per_r2_adaptor = "numeric",
        refseq = "character",
        pamseq = "character",
        libcounts = "data.frame",
        allcounts = "data.frame",
        valiant_meta = "data.frame",
        vep_anno = "data.frame",
        meta_mseqs = "character",
        missing_meta_seqs = "character",
        libstats = "data.frame",
        allstats = "data.frame",
        libstats_qc = "data.frame",
        allstats_qc = "data.frame"
    ),
    prototype = list(
        sample = character(),
        libname = character(),
        libtype = character(),
        adapt5 = character(),
        adapt3 = character(),
        per_r1_adaptor = numeric(),
        per_r2_adaptor = numeric(),
        refseq = character(),
        pamseq = character(),
        libcounts = data.frame(),
        allcounts = data.frame(),
        valiant_meta = data.frame(),
        vep_anno = data.frame(),
        meta_mseqs = character(),
        missing_meta_seqs = character(),
        libstats = data.frame(),
        allstats = data.frame(),
        libstats_qc = data.frame(),
        allstats_qc = data.frame()
    )
)

#' Create a new SGE object
#'
#' @export
#' @name create_sge_object
#' @param file_libcount            QUANTS library-dependent count file, per sequence per count
#' @param file_allcount            QUANTS library-independent count file, per sequence per count
#' @param file_valiant_meta        VaLiAnT meta file
#' @param file_vep_anno            vep annotation file
#' @param file_libcount_hline      line number of header in library-dependent count file
#' @param file_allcount_hline      line number of header in library-independent count file
#' @param file_valiant_meta_hline  line number of header in VaLiAnT meta file
#' @param file_vep_anno_hline      line number of header in vep annotation file
#' @param file_libcount_cols       a vector of numbers of selected columns in library-dependent count file, default is none
#' @param file_allcount_cols       a vector of numbers of selected columns in library-independent count file, default is none
#' @param file_valiant_meta_cols   a vector of numbers of selected columns in VaLiAnT meta file, default is none
#' @param file_vep_anno_cols       a vector of numbers of selected columns in vep annotation file, default is none
#' @return An object of class SGE
create_sge_object <- function(file_libcount,
                              file_allcount,
                              file_valiant_meta,
                              file_vep_anno = NULL,
                              file_libcount_hline = 3,
                              file_allcount_hline = 3,
                              file_valiant_meta_hline = 1,
                              file_vep_anno_hline = 1,
                              file_libcount_cols = vector(),
                              file_allcount_cols = vector(),
                              file_valiant_meta_cols = vector(),
                              file_vep_anno_cols = vector()) {
    # Read files
    libcounts <- read_sge_file(file_libcount, file_libcount_hline, file_libcount_cols)
    allcounts <- read_sge_file(file_allcount, file_allcount_hline, file_allcount_cols)
    valiant_meta <- read_sge_file(file_valiant_meta, file_valiant_meta_hline, file_valiant_meta_cols)

    # vep is only required for screen qc
    if (is.null(file_vep_anno)) {
        vep_anno <- data.frame()
    } else {
        vep_anno <- read_sge_file(file_vep_anno, file_vep_anno_hline, file_vep_anno_cols)

        if (vep_anno[1, ]$unique_oligo_name %nin% valiant_meta$oligo_name) {
            vep_anno$unique_oligo_name <- sapply(vep_anno$unique_oligo_name, function (x) paste(head(unlist(strsplit(x, "_")), -2), collapse = "_"))
        }

        if (vep_anno[1, ]$seq %nin% libcounts$sequence) {
            if (revcomp(vep_anno[1, ]$seq) %in% libcounts$sequence) {
                vep_anno$seq <- sapply(vep_anno$seq, function (x) revcomp(x))
            }
        }
    }

    # initializing
    cols <- c("total_num_oligos",
              "total_num_unique_oligos",
              "total_counts",
              "max_counts",
              "min_counts",
              "median_counts",
              "mean_counts",
              "num_oligos_nocount",
              "num_oligos_lowcount",
              "max_len_oligos",
              "min_len_oligos")
    df_libstats <- data.frame(matrix(NA, 1, length(cols)))
    colnames(df_libstats) <- cols

    cols <- c("total_num_oligos",
              "total_num_unique_oligos",
              "total_counts",
              "max_counts",
              "min_counts",
              "median_counts",
              "mean_counts",
              "num_oligos_nocount",
              "num_oligos_lowcount",
              "max_len_oligos",
              "min_len_oligos")
    df_allstats <- data.frame(matrix(NA, 1, length(cols)))
    colnames(df_allstats) <- cols

    cols <- c("num_ref_reads",
              "per_ref_reads",
              "num_pam_reads",
              "per_pam_reads",
              "num_eff_reads",
              "per_eff_reads",
              "num_unmapped_reads",
              "per_unmapped_reads",
              "num_missing_var",
              "per_missing_var",
              "gini_coeff")
    df_libstats_qc <- data.frame(matrix(NA, 1, length(cols)))
    colnames(df_libstats_qc) <- cols

    cols <- c("num_ref_reads",
              "per_ref_reads",
              "num_pam_reads",
              "per_pam_reads",
              "num_eff_reads",
              "per_eff_reads",
              "num_unmapped_reads",
              "per_unmapped_reads",
              "num_missing_var",
              "per_missing_var",
              "gini_coeff")
    df_allstats_qc <- data.frame(matrix(NA, 1, length(cols)))
    colnames(df_allstats_qc) <- cols

    # Create the object
    sge_object <- new("SGE",
        libcounts = libcounts,
        allcounts = allcounts,
        valiant_meta = valiant_meta,
        vep_anno = vep_anno,
        libstats = df_libstats,
        allstats = df_allstats,
        libstats_qc = df_libstats_qc,
        allstats_qc = df_allstats_qc)

    sge_object@libcounts <- as.data.table(sge_object@libcounts)
    sge_object@allcounts <- as.data.table(sge_object@allcounts)
    sge_object@valiant_meta <- as.data.table(sge_object@valiant_meta)
    sge_object@vep_anno <- as.data.table(sge_object@vep_anno)

    # Return the object
    return(sge_object)
}

#' A class representing a sample QC object
#'
#' @export
#' @name sampleQC
#' @slot cutoffs                  a data frame of cutoffs using in sample QC
#' @slot samples                  a list of SGE objects
#' @slot samples_ref              a list of SGE objects which are the references for screen QC
#' @slot counts                   a list of sample libraray-independent counts
#' @slot lengths                  a list of sequence lengths
#' @slot seq_clusters             a list of dataframes of sequences and cluster IDs
#' @slot accepted_counts          a list of filtered counts of all the samples
#' @slot library_counts           a list of library counts of all the samples
#' @slot unmapped_counts          a list of unmapped counts against meta library sequences of all the samples
#' @slot library_counts_chr       a list of chromosomes for library counts
#' @slot library_counts_pos       a list of library counts of all the samples sorted by position in meta
#' @slot library_counts_anno      a data frame of library counts of all the samples, annotated with consequences
#' @slot library_counts_pos_anno  a data frame of library counts of all the samples, annotated with consequences, sorted by position
#' @slot stats                    a data frame of samples and stats, eg. total no, filtered no.
#' @slot bad_seqs_bycluster       a list of filter-out sequences by cluster
#' @slot bad_seqs_bydepth         a list of filter-out sequences by depth
#' @slot bad_seqs_bylib           a list of filter-out sequences by library mapping
#' @slot filtered_samples         a vector of filtered sample names
setClass("sampleQC",
    slots = list(
        cutoffs = "data.frame",
        samples = "list",
        samples_ref = "list",
        counts = "list",
        lengths = "list",
        seq_clusters = "list",
        accepted_counts = "list",
        library_counts = "list",
        unmapped_counts = "list",
        library_counts_chr = "list",
        library_counts_pos = "list",
        library_counts_anno = "data.frame",
        library_counts_pos_anno = "data.frame",
        stats = "data.frame",
        bad_seqs_bycluster = "list",
        bad_seqs_bydepth = "list",
        bad_seqs_bylib = "list",
        filtered_samples = "character"
    ),
    prototype = list(
        cutoffs = data.frame(),
        samples = list(),
        samples_ref = list(),
        counts = list(),
        lengths = list(),
        seq_clusters = list(),
        accepted_counts = list(),
        library_counts = list(),
        unmapped_counts = list(),
        library_counts_chr = list(),
        library_counts_pos = list(),
        library_counts_anno = data.frame(),
        library_counts_pos_anno = data.frame(),
        stats = data.frame(),
        bad_seqs_bycluster = list(),
        bad_seqs_bydepth = list(),
        bad_seqs_bylib = list(),
        filtered_samples = character()
    )
)

#' Create a new sample QC object
#'
#' @export
#' @name create_sampleqc_object
#' @param samples a list of SGE objects
#' @return An object of class sampleQC
create_sampleqc_object <- function(samples) {
    # checking
    if (length(samples) == 0) {
         stop(paste0("====> Error: no sample found in the input!"))
    }

    # initializing
    num_samples <- length(samples)
    sample_names <- character()
    for (s in samples) {
        sample_names <- append(sample_names, s@sample)
    }
    if (anyDuplicated(sample_names) != 0) {
        dup_names <- paste(unique(sample_names[duplicated(sample_names)]), collapse = ",")
        stop(paste0("====> Error: duplicated sample names:", " ", dup_names))
    }

    # get reference sampels if ref_time_point is not null in the sample sheet
    if (length(maveqc_ref_time_point_samples) >= 1) {
        ref_samples <- select_objects(samples, maveqc_ref_time_point_samples)
    } else {
        ref_samples <- list()
    }

    list_counts <- list()
    list_lengths <- list()
    for (s in samples) {
        counts <- s@allcounts[, c("sequence", "count")]

        lengths <- s@allcounts[, "sequence", drop = FALSE]
        lengths$length <- nchar(lengths$sequence) - nchar(s@adapt5) - nchar(s@adapt3)

        # in case negative length after deduction
        lengths$length[lengths$length < 0] <- 0

        list_counts[[s@sample]] <- counts
        list_lengths[[s@sample]] <- lengths
    }

    cols <- c("per_r1_adaptor",
              "per_r2_adaptor",
              "total_reads",
              "excluded_reads",
              "accepted_reads",
              "library_seqs",
              "missing_meta_seqs",
              "per_missing_meta_seqs",
              "library_reads",
              "per_library_reads",
              "unmapped_reads",
              "per_unmapped_reads",
              "ref_reads",
              "per_ref_reads",
              "pam_reads",
              "per_pam_reads",
              "median_cov",
              "library_cov",
              "gini_coeff_before_qc",
              "gini_coeff_after_qc",
              "qcpass_total_reads",
              "qcpass_missing_per",
              "qcpass_accepted_reads",
              "qcpass_mapping_per",
              "qcpass_ref_per",
              "qcpass_library_per",
              "qcpass_library_cov",
              "qcpass")
    df_stats <- data.frame(matrix(NA, num_samples, length(cols)))
    rownames(df_stats) <- sample_names
    colnames(df_stats) <- cols

    # Create the object
    sampleqc_object <- new("sampleQC",
        samples = samples,
        samples_ref = ref_samples,
        counts = list_counts,
        lengths = list_lengths,
        stats = df_stats)

    # Return the object
    return(sampleqc_object)
}

setClass("hclust")
setClass("prcomp")
#' A class representing a experiment QC object
#'
#' @export
#' @name experimentQC
#' @slot samples                    a list of SGE objects
#' @slot coldata                    a data frame of coldata for DESeq2
#' @slot ref_condition              the reference condition, like D4, others VS D4 in DESeq2
#' @slot vep_anno                   a data frame of consequence annotations (should be the same in all the samples for screen qc)
#' @slot accepted_counts            a data frame of accepted counts of all the samples
#' @slot library_counts_anno        a data frame of library counts of all the samples, annotated with consequences
#' @slot library_counts_pos_anno    a data frame of library counts of all the samples, annotated with consequences, sorted by position
#' @slot comparisons                a list of comparisons for degComps
#' @slot lib_deseq_rlog             a data frame of deseq rlog counts of all the samples using library counts
#' @slot lib_hclust_res             a hclust object for all the samples using library counts
#' @slot lib_corr_res               the correlation results for all the samples using library counts
#' @slot lib_pca_res                the pca results for all the samples using library counts
#' @slot lib_deseq_res              a list of deseq results of all the comparison against reference using library counts
#' @slot lib_deseq_res_anno         a list of deseq results with consequence annotations using library counts
#' @slot all_deseq_rlog             a data frame of deseq rlog counts of all the samples using all counts
#' @slot all_deseq_res              a list of deseq results of all the comparison against reference using all counts
#' @slot all_deseq_res_anno         a list of deseq results with consequence annotations using all counts
#' @slot all_deseq_res_anno_adj     a list of deseq results with consequence annotations using all counts and adjusted lfc and p value
setClass("experimentQC",
    slots = list(
        samples = "list",
        coldata = "data.frame",
        ref_condition = "character",
        vep_anno = "data.frame",
        accepted_counts = "data.frame",
        library_counts_anno = "data.frame",
        library_counts_pos_anno = "data.frame",
        comparisons = "list",
        lib_deseq_rlog = "data.frame",
        lib_hclust_res = "hclust",
        lib_corr_res = "matrix",
        lib_pca_res = "prcomp",
        lib_deseq_res = "list",
        lib_deseq_res_anno = "list",
        all_deseq_rlog = "data.frame",
        all_deseq_res = "list",
        all_deseq_res_anno = "list",
        all_deseq_res_anno_adj = "list"
    ),
    prototype = list(
        samples = list(),
        coldata = data.frame(),
        ref_condition = character(),
        vep_anno = data.frame(),
        accepted_counts = data.frame(),
        library_counts_anno = data.frame(),
        library_counts_pos_anno = data.frame(),
        comparisons = list(),
        lib_deseq_rlog = data.frame(),
        lib_hclust_res = hclust(dist(matrix(seq(1:9), nrow = 3))),
        lib_corr_res = matrix(),
        lib_pca_res = prcomp(as.data.frame(matrix(round(runif(n = 25, min = 1, max = 20), 0), nrow = 5))),
        lib_deseq_res = list(),
        lib_deseq_res_anno = list(),
        all_deseq_rlog = data.frame(),
        all_deseq_res =  list(),
        all_deseq_res_anno =  list(),
        all_deseq_res_anno_adj = list()
    )
)

#' Create a new experiment QC object
#'
#' @export
#' @name create_experimentqc_object
#' @param samqc_obj a sampleQC object
#' @param coldata   a data frame of coldata for DESeq2
#' @param refcond   the reference condition, eg. D4
#' @return An object of class sampleQC
create_experimentqc_object <- function(samqc_obj,
                                       coldata = maveqc_deseq_coldata,
                                       refcond = maveqc_ref_time_point) {
    # checking
    if (is.null(coldata)) {
         stop(paste0("====> Error: no coldata found in the input!"))
    }

    if (refcond %nin% coldata$condition) {
        stop(paste0("====> Error: reference condition is not in the coldata!"))
    }

    # initializing
    if ("condition" %nin% colnames(coldata)) {
        stop(paste0("====> Error: coldata must have condition values!"))
    } else {
        coldata <- as.data.frame(coldata)

        coldata$condition <- factor(coldata$condition)
        coldata$condition <- factor(coldata$condition, levels = mixedsort(levels(coldata$condition)))

        coldata$replicate <- factor(coldata$replicate)
        coldata$replicate <- factor(coldata$replicate, levels = mixedsort(levels(coldata$replicate)))
    }

    conds <- levels(coldata$condition)
    ds_contrast <- list()
    for (i in 1:length(conds)) {
        if (conds[i] != refcond) {
            ds_contrast <- append(ds_contrast, paste0("condition_", conds[i], "_vs_", refcond))
        }
    }

    # Create the object
    experimentqc_object <- new("experimentQC",
        samples = samqc_obj@samples,
        coldata = coldata,
        ref_condition = refcond,
        vep_anno = samqc_obj@samples[[1]]@vep_anno,
        accepted_counts = merge_list_to_dt(samqc_obj@accepted_counts, "sequence", "count"),
        library_counts_anno = samqc_obj@library_counts_anno,
        library_counts_pos_anno = samqc_obj@library_counts_pos_anno,
        comparisons = ds_contrast)

    experimentqc_object@vep_anno <- as.data.table(experimentqc_object@vep_anno)
    experimentqc_object@accepted_counts <- as.data.table(experimentqc_object@accepted_counts)
    experimentqc_object@library_counts_anno <- as.data.table(experimentqc_object@library_counts_anno)
    experimentqc_object@library_counts_pos_anno <- as.data.table(experimentqc_object@library_counts_pos_anno)

    # Return the object
    return(experimentqc_object)
}

#' initialize function
setGeneric("format_count", function(object, ...) {
  standardGeneric("format_count")
})

#' format library dependent and independent counts with extra info and remove duplicate/useless info
#'
#' @export
#' @name format_count
#' @param object SGE object
#' @return object
setMethod(
    "format_count",
    signature = "SGE",
    definition = function(object) {
        #----------#
        # checking #
        #----------#
        if (length(object@adapt5) == 0 | length(object@adapt3) == 0) {
            if ((length(object@refseq) == 0)) {
                stop(paste0("====> Error: no reference sequence, please provide adaptor sequences instead!"))
            }

            if ((length(object@pamseq) == 0)) {
                stop(paste0("====> Error: no pam sequence, please provide adaptor sequences instead!"))
            }
        }

        #----------------------------#
        # 1. valiant ref and pam seq #
        #----------------------------#
        if ((length(object@refseq) == 0)) {
            tmp_refseq <- unique(object@valiant_meta$ref_seq)

            if (tmp_refseq %in% object@allcounts$sequence) {
                object@refseq <- tmp_refseq
            } else {
                tmp_refseq <- revcomp(tmp_refseq)
                if (tmp_refseq %in% object@allcounts$sequence) {
                    object@refseq <- tmp_refseq
                } else {
                    object@refseq <- unique(object@valiant_meta$ref_seq)
                    message("         Warning: reference sequence cannot be found in library-independent counts.")
                }
            }
        }

        if ((length(object@pamseq) == 0)) {
            tmp_pamseq <- unique(object@valiant_meta$pam_seq)

            if (tmp_pamseq %in% object@allcounts$sequence) {
                object@pamseq <- tmp_pamseq
            } else {
                tmp_pamseq <- revcomp(tmp_pamseq)
                if (tmp_pamseq %in% object@allcounts$sequence) {
                    object@pamseq <- tmp_pamseq
                } else {
                    object@pamseq <- unique(object@valiant_meta$pam_seq)
                    message("         Warning: pam sequence cannot be found in library-independent counts.")
                }
            }
        }

        #----------------------------#
        # 2. library dependent count #
        #----------------------------#
        object@libcounts[, is_ref := fifelse(sequence == object@refseq, 1, 0)]
        object@libcounts[, is_pam := fifelse(sequence == object@pamseq, 1, 0)]

        #------------------------------#
        # 3. library independent count #
        #------------------------------#
        object@allcounts[, is_ref := fifelse(sequence == object@refseq, 1, 0)]
        object@allcounts[, is_pam := fifelse(sequence == object@pamseq, 1, 0)]

        #--------------------------#
        # 4. mseq in valiant meta  #
        #--------------------------#
        # use library dependent sequences to get meta mseqs
        # may change as meta and counts will have the same seqs with adaptors
        # library dependent sequences are not unique
        # library dependent sequences have no ref and pam, but independent do

        tmp_mseqs <- unique(object@libcounts$sequence)
        object@meta_mseqs <- tmp_mseqs[tmp_mseqs %nin% c(object@refseq, object@pamseq)]

        #object@missing_meta_seqs <- object@meta_mseqs[object@meta_mseqs %nin% object@allcounts$sequence]
        object@missing_meta_seqs <- unique(object@libcounts[count == 0]$sequence)

        return(object)
    }
)

#' initialize function
setGeneric("sge_stats", function(object, ...) {
  standardGeneric("sge_stats")
})

#' format library dependent and independent counts with extra info and remove duplicate/useless info
#'
#' @export
#' @name sge_stats
#' @param object SGE object
#' @param lowcut cutoff which determines the oligo count is low, user's definition
#' @return object
setMethod(
    "sge_stats",
    signature = "SGE",
    definition = function(object,
                          lowcut = 10) {
        # library dependent counts
        unique_libcounts <- unique(object@libcounts[, c("sequence", "count")])

        object@libstats$total_num_oligos <- nrow(object@libcounts)
        object@libstats$total_num_unique_oligos <- nrow(unique_libcounts)
        object@libstats$total_counts <- sum(unique_libcounts$count)
        object@libstats$max_counts <- max(unique_libcounts$count)
        object@libstats$min_counts <- min(unique_libcounts$count)
        object@libstats$median_counts <- median(unique_libcounts$count)
        object@libstats$mean_counts <- mean(unique_libcounts$count)
        object@libstats$num_oligos_nocount <- nrow(unique_libcounts[unique_libcounts$count == 0, ])
        object@libstats$num_oligos_lowcount <- nrow(unique_libcounts[unique_libcounts$count <= lowcut, ])
        object@libstats$max_len_oligos <- max(nchar(unique_libcounts$sequence))
        object@libstats$min_len_oligos <- min(nchar(unique_libcounts$sequence))

        # library independent counts
        object@allstats$total_num_oligos <- nrow(object@allcounts)
        object@allstats$total_num_unique_oligos <- nrow(object@allcounts[object@allcounts$unique == 1, ])
        object@allstats$total_counts <- sum(object@allcounts$count)
        object@allstats$max_counts <- max(object@allcounts$count)
        object@allstats$min_counts <- min(object@allcounts$count)
        object@allstats$median_counts <- median(object@allcounts$count)
        object@allstats$mean_counts <- mean(object@allcounts$count)
        object@allstats$num_oligos_nocount <- nrow(object@allcounts[object@allcounts$count == 0, ])
        object@allstats$num_oligos_lowcount <- nrow(object@allcounts[object@allcounts$count <= lowcut, ])
        object@allstats$max_len_oligos <- max(nchar(object@allcounts$sequence))
        object@allstats$min_len_oligos <- min(nchar(object@allcounts$sequence))

        return(object)
    }
)

#' initialize function
setGeneric("sge_qc_stats", function(object, ...) {
  standardGeneric("sge_qc_stats")
})

#' format library dependent and independent counts with extra info and remove duplicate/useless info
#'
#' @export
#' @name sge_qc_stats
#' @param object SGE object
#' @param lowcut cutoff which determines the oligo count is low, user's definition
#' @return object
setMethod(
    "sge_qc_stats",
    signature = "SGE",
    definition = function(object) {
        # issue: total_counts is counts, not no. of seqeunced reads, need from qc report
        # now assume total counts of library independent is total no of reads
        total_num_sequenced_reads <- object@allstats$total_counts

        # library dependent counts -- notice lib-dependent sequences are not unique
        unique_libcounts <- unique(object@libcounts[, c("sequence", "count", "is_ref", "is_pam")])

        qc_count <- as.numeric(unique_libcounts[unique_libcounts$is_ref == 1, "count"])
        object@libstats_qc$num_ref_reads <- ifelse(is.na(qc_count), 0, qc_count)
        object@libstats_qc$per_ref_reads <- object@libstats_qc$num_ref_reads / total_num_sequenced_reads * 100
        object@libstats_qc$per_ref_reads <- round(object@libstats_qc$per_ref_reads, 2)

        qc_count <- as.numeric(unique_libcounts[unique_libcounts$is_pam == 1, "count"])
        object@libstats_qc$num_pam_reads <- ifelse(is.na(qc_count), 0, qc_count)
        object@libstats_qc$per_pam_reads <- object@libstats_qc$num_pam_reads / total_num_sequenced_reads * 100
        object@libstats_qc$per_pam_reads <- round(object@libstats_qc$per_pam_reads, 2)

        qc_count <- sum(unique_libcounts[unique_libcounts$is_ref == 0 & unique_libcounts$is_pam == 0, "count"])
        object@libstats_qc$num_eff_reads <- ifelse(length(qc_count) == 0, 0, qc_count)
        object@libstats_qc$per_eff_reads <- object@libstats_qc$num_eff_reads / total_num_sequenced_reads * 100
        object@libstats_qc$per_eff_reads <- round(object@libstats_qc$per_eff_reads, 2)

        qc_count <- total_num_sequenced_reads - object@libstats_qc$num_ref_reads - object@libstats_qc$num_pam_reads - object@libstats_qc$num_eff_reads
        object@libstats_qc$num_unmapped_reads <- qc_count
        object@libstats_qc$per_unmapped_reads <- object@libstats_qc$num_unmapped_reads / total_num_sequenced_reads * 100
        object@libstats_qc$per_unmapped_reads <- round(object@libstats_qc$per_unmapped_reads, 2)

        object@libstats_qc$num_missing_var <- length(object@missing_meta_seqs)
        object@libstats_qc$per_missing_var <- object@libstats_qc$num_missing_var / length(object@meta_mseqs) * 100
        object@libstats_qc$per_missing_var <- round(object@libstats_qc$per_missing_var, 2)

        object@libstats_qc$gini_coeff <- cal_gini(unique_libcounts$count, corr = FALSE, na.rm = TRUE)
        object@libstats_qc$gini_coeff <- round(object@libstats_qc$gini_coeff, 3)

        # library independent counts
        qc_count <- as.numeric(object@allcounts[object@allcounts$is_ref == 1, "count"])
        object@allstats_qc$num_ref_reads <- ifelse(is.na(qc_count), 0, qc_count)
        object@allstats_qc$per_ref_reads <- object@allstats_qc$num_ref_reads / total_num_sequenced_reads * 100
        object@allstats_qc$per_ref_reads <- round(object@allstats_qc$per_ref_reads, 2)

        qc_count <- as.numeric(object@allcounts[object@allcounts$is_pam == 1, "count"])
        object@allstats_qc$num_pam_reads <- ifelse(is.na(qc_count), 0, qc_count)
        object@allstats_qc$per_pam_reads <- object@allstats_qc$num_pam_reads / total_num_sequenced_reads * 100
        object@allstats_qc$per_pam_reads <- round(object@allstats_qc$per_pam_reads, 2)

        qc_count <- sum(object@allcounts[object@allcounts$is_ref == 0 & object@allcounts$is_pam == 0, "count"])
        object@allstats_qc$num_eff_reads <- ifelse(length(qc_count) == 0, 0, qc_count)
        object@allstats_qc$per_eff_reads <- object@allstats_qc$num_eff_reads / total_num_sequenced_reads * 100
        object@allstats_qc$per_eff_reads <- round(object@allstats_qc$per_eff_reads, 2)

        qc_count <- total_num_sequenced_reads - object@allstats_qc$num_ref_reads - object@allstats_qc$num_pam_reads - object@allstats_qc$num_eff_reads
        object@allstats_qc$num_unmapped_reads <- qc_count
        object@allstats_qc$per_unmapped_reads <- object@allstats_qc$num_unmapped_reads / total_num_sequenced_reads * 100
        object@allstats_qc$per_unmapped_reads <- round(object@allstats_qc$per_unmapped_reads, 2)

        object@allstats_qc$gini_coeff <- cal_gini(object@allcounts$count, corr = FALSE, na.rm = TRUE)
        object@allstats_qc$gini_coeff <- round(object@allstats_qc$gini_coeff, 3)

        return(object)
    }
)

# - process.R
#' transparent color function
#'
#' creating a transparent color using color name,
#' alpha rate is 0 to 1
#'
#' @name t_col
#' @param col  color name
#' @param rate alpha rate
#' @return transparent color
t_col <- function(col, rate) {
    newcol <- rgb(col2rgb(col)["red", ],
                  col2rgb(col)["green", ],
                  col2rgb(col)["blue", ],
                  as.integer(rate * 255),
                  maxColorValue = 255)
    return(newcol)
}

#' capitalise the first character in the string
#'
#' creating the capitalised names for data frame
#'
#' @name capital_names
#' @param x  a vector of strings
#' @return capitalised strings
capital_names <- function(x) {
    y <- x
    for (i in 1:length(x)) {
        y[i] <- paste(toupper(substring(x[i], 1, 1)), substring(x[i], 2), sep = "", collapse = " ")
    }
    return(y)
}

#' not in function
#'
#' @name %nin%
#' @param x X
#' @param y Y
#' @return True or False
`%nin%` <- function(x, y) !(x %in% y)

#' reverse complement
#'
#' @name revcomp
#' @param seq sequence
#' @return string
revcomp <- function(seq) {
    seq <- toupper(seq)
    splits <- strsplit(seq, "")[[1]]
    reversed <- rev(splits)
    seq_rev <- paste(reversed, collapse = "")
    seq_rev_comp <- chartr("ATCG", "TAGC", seq_rev)
    return(seq_rev_comp)
}

#' trim adaptor sequences
#'
#' @name trim_adaptor
#' @param seq    sequence
#' @param adapt5 5 prime adaptor sequence
#' @param adapt3 3 prime adaptor sequence
#' @return string
trim_adaptor <- function(seq, adapt5, adapt3) {
    adapt5_pos <- regexpr(adapt5, seq, fixed = TRUE)[1]
    adapt3_pos <- regexpr(adapt3, seq, fixed = TRUE)[1]

    is_revcomp <- FALSE
    # ? could adaptor revcomp ?
    if (adapt5_pos < 0 & adapt3_pos < 0) {
        adapt5_revcomp <- revcomp(adapt5)
        adapt3_revcomp <- revcomp(adapt3)

        adapt5_revcomp_pos <- regexpr(adapt5_revcomp, seq, fixed = TRUE)[1]
        adapt3_revcomp_pos <- regexpr(adapt3_revcomp, seq, fixed = TRUE)[1]

        if (adapt3_revcomp_pos < 0 & adapt5_revcomp_pos < 0) {
            return(seq)
        } else {
            is_revcomp <- TRUE
        }
    }

    if (is_revcomp == FALSE) {
        if (adapt3_pos > adapt5_pos) {
            if (adapt5_pos > 0 & adapt3_pos > 0) {
                return(substr(seq, adapt5_pos + nchar(adapt5), adapt3_pos - 1))
            } else if (adapt5_pos > 0 & adapt3_pos < 0) {
                return(substr(seq, adapt5_pos + nchar(adapt5), nchar(seq)))
            } else if (adapt5_pos < 0 & adapt3_pos > 0) {
                return(substr(seq, 1, adapt3_pos - 1))
            }
        } else {
            stop(paste0("====> Error: 3 prime adaptor found before 5 prime adaptor in the sequence: ", seq))
        }
    } else {
        if (adapt5_revcomp_pos > adapt3_revcomp_pos) {
            if (adapt3_revcomp_pos > 0 & adapt5_revcomp_pos > 0) {
                return(substr(seq, adapt3_revcomp_pos + nchar(adapt3_revcomp), adapt5_revcomp_pos - 1))
            } else if (adapt3_revcomp_pos > 0 & adapt5_revcomp_pos < 0) {
                return(substr(seq, adapt3_revcomp_pos + nchar(adapt3_revcomp), nchar(seq)))
            } else if (adapt3_revcomp_pos < 0 & adapt5_revcomp_pos > 0) {
                return(substr(seq, 1, adapt3_revcomp_pos - 1))
            }
        } else {
            stop(paste0("====> Error: 5 prime adaptor (RC) found before 3 prime adaptor (RC) in the sequence: ", seq))
        }
    }
}

#' column binding with filling NAs
#'
#' @name cbind_fill
#' @return matrix
cbind_fill <- function(...) {
    nm <- list(...)
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow))
    do.call(cbind, lapply(nm, function (x) rbind(x, matrix(, n - nrow(x), ncol(x)))))
}

#' calculate gini coefficiency for a sample
#'
#' @name cal_gini
#' @param x a vector
#' @return a value
cal_gini <- function(x, corr = FALSE, na.rm = TRUE) {
    if (!na.rm && any(is.na(x))) return(NA_real_)
    x <- as.numeric(na.omit(x))
    n <- length(x)
    x <- sort(x)
    G <- sum(x * 1L:n)
    G <- 2 * G/sum(x) - (n + 1L)

    if (corr) {
        return(G / (n - 1L))
    } else {
        return(G / n)
    }
}

#' merge a list of data tables into a data table
#'
#' @name merge_list_to_dt
#' @param list_dt   a list of data tables
#' @param by_val    join data tables by which column
#' @param join_val  join which column in the data tables
#' @return a data table
merge_list_to_dt <- function(list_dt, by_val, join_val) {
    dt_out <- data.table()

    for (i in 1:length(list_dt)) {
        cols <- c(by_val, join_val)
        dt_tmp <- list_dt[[i]][, ..cols]

        if (nrow(dt_out) == 0) {
            dt_out <- dt_tmp
            colnames(dt_out) <- c(by_val, names(list_dt)[i])
        } else {
            coln <- colnames(dt_out)
            dt_out <- merge(dt_out, dt_tmp, by = by_val, all = TRUE)
            colnames(dt_out) <- c(coln, names(list_dt)[i])
        }
    }

    return(dt_out)
}

#' color blind friendly
#'
#' @name select_colorblind
#' @param col_id a character to select colors
#' @return a vector of colors
select_colorblind <- function(col_id) {
    col8 <- c("#D55E00", "#56B4E9", "#E69F00",
              "#009E73", "#F0E442", "#0072B2",
              "#CC79A7", "#000000")

    col12 <- c("#88CCEE", "#CC6677", "#DDCC77",
               "#117733", "#332288", "#AA4499",
               "#44AA99", "#999933", "#882255",
               "#661100", "#6699CC", "#888888")

    col15 <- c("red",       "royalblue", "olivedrab",
               "purple",    "violet",    "maroon1",
               "seagreen1", "navy",      "pink",
               "coral",     "steelblue", "turquoise1",
               "red4",      "skyblue",   "yellowgreen")

    col21 <- c("#F60239", "#009503", "#FFDC3D",
               "#9900E6", "#009FFA", "#FF92FD",
               "#65019F", "#FF6E3A", "#005A01",
               "#00E5F8", "#DA00FD", "#AFFF2A",
               "#00F407", "#00489E", "#0079FA",
               "#560133", "#EF0096", "#000000",
               "#005745", "#00AF8E", "#00EBC1")

    if (col_id == "col8") {
        return(col8)
    } else if (col_id == "col12") {
        return(col12)
    } else if (col_id == "col15") {
        return(col15)
    } else if (col_id == "col21") {
        return(col21)
    } else {
        stop(paste0("====> Error: wrong col_id"))
    }
}

#' fetch objects from list by the names or indexes
#'
#' @export
#' @name select_objects
#' @param objects a list of objects
#' @param tags    a vector of names or indexes
#' @return a list of objects
select_objects <- function(objects, tags) {
    if (length(objects) == 0) {
        stop(paste0("====> Error: no object found in the list"))
    }

    if (length(tags) == 0) {
        stop(paste0("====> Error: please provide tags to fetch objects"))
    } else {
        if (class(tags) %in% c("character", "numeric")) {
            if (class(tags) == "numeric") {
                tags <- as.integer(tags)
            }
        } else {
            stop(paste0("====> Error: wrong tag type, must be integer or character"))
        }
    }

    list_select <- list()
    if (class(tags) == "character") {
        for (i in 1:length(tags)) {
            for (j in 1:length(objects)) {
                if (tags[i] == objects[[j]]@sample) {
                    list_select <- append(list_select, objects[[j]])
                    break
                }
            }
        }
    } else {
        for (i in 1:length(tags)) {
            list_select <- append(list_select, objects[[tags[i]]])
        }
    }

    return(list_select)
}


# - ingest.R
#' check file format and read it
#' file has header but beginning with #
#'
#' @export
#' @name read_sge_file
#' @param file_path the file path
#' @param hline     header line number, default is 0
#' @param colnums   a vector of selected colummns, default is none
#' @return a dataframe
read_sge_file <- function(file_path,
                          hline = 0,
                          colnums = vector()) {
    # check if file path is a symbolic link
    if (!file.exists(file_path)) {
        stop(paste0("====> Error: ", file_path, " doesn't exist"))
    } else {
        tmp_path <- Sys.readlink(file_path)
        file_path <- ifelse(nzchar(tmp_path), tmp_path, file_path)
    }

    if (hline < 0) {
        stop(paste0("====> Error: ", hline, " must be >= 0"))
    }

    # read data and check file is csv or tsv
    # speed: vroom > fread > read.table
    csv_pattern <- "\\.csv(\\.gz)?$"
    tsv_pattern <- "\\.tsv(\\.gz)?$"
    if (grepl(csv_pattern, file_path)) {
        #filedata <- read.table(file_path, sep = ",", comment.char = "#", skip = hline)
        #filedata <- fread(file_path, sep = ",", skip = hline)        
        suppressWarnings(filedata <- vroom(file_path, delim = ",", comment = "#", skip = hline, col_names = FALSE, show_col_types = FALSE))
    } else if (grepl(tsv_pattern, file_path)) {
        #filedata <- read.table(file_path, sep = "\t", comment.char = "#", skip = hline)
        #filedata <- fread(file_path, sep = "\t", skip = hline)
        suppressWarnings(filedata <- vroom(file_path, delim = "\t", comment = "#", skip = hline, col_names = FALSE, show_col_types = FALSE))
    } else {
        stop(paste0("====> Error: wrong format, ", file_path, " is not .csv(.gz) or .tsv(.gz)!"))
    }

    # examine header
    if (hline > 0) {
        headers <- list()
        conn <- file(file_path, "r")
        lines <- readLines(conn, n = hline)
        for (l in lines) {
            if (length(l) > 0) {
                headers <- append(headers, lines)
            }
        }
        close(conn)

        if (grepl(csv_pattern, file_path)) {
            header <- strsplit(sub("#", "", headers[hline]), ",")[[1]]
        } else {
            header <- strsplit(sub("#", "", headers[hline]), "\t")[[1]]
        }
    }

    filedata <- as.data.frame(filedata)
    # add header to dataframe, and transfer to lower characters
    if (hline > 0) {
        colnames(filedata) <- tolower(header)
    }
    # select columns
    if (length(colnums) > 0) {
        filedata <- filedata[, colnums]
    }

    return(filedata)
}

#' import files from the sample sheet
#' create objects based on the sample names
#'
#' @export
#' @name import_sge_files
#' @param dir_path                the directory path
#' @param sample_sheet            the file name of the sample sheet which is in the directory
#' @param file_libcount_hline     line number of header in library-dependent count file
#' @param file_allcount_hline     line number of header in library-independent count file
#' @param file_valiant_meta_hline line number of header in VaLiAnT meta file
#' @param file_vep_anno_hline     line number of header in vep annotation file
#' @param file_libcount_cols      a vector of numbers of selected columns in library-dependent count file, default is none
#' @param file_allcount_cols      a vector of numbers of selected columns in library-independent count file, default is none
#' @param file_valiant_meta_cols  a vector of numbers of selected columns in VaLiAnT meta file, default is none
#' @param file_vep_anno_cols      a vector of numbers of selected columns in vep annotation file, default is none
#' @return a list of objects
import_sge_files <- function(dir_path = NULL,
                             sample_sheet = NULL,
                             file_libcount_hline = 3,
                             file_allcount_hline = 3,
                             file_valiant_meta_hline = 1,
                             file_vep_anno_hline = 1,
                             file_libcount_cols = vector(),
                             file_allcount_cols = vector(),
                             file_valiant_meta_cols = vector(),
                             file_vep_anno_cols = vector()) {
    # check input format
    if (is.null(dir_path)) {
        stop(paste0("====> Error: please provide the directory of input files!"))
    } else {
        if (!dir.exists(dir_path)) {
            stop(paste0("====> Error: ", dir_path, " doesn't exist"))
        }
    }

    if (is.null(sample_sheet)) {
        stop(paste0("====> Error: please provide the file name of sample sheet!"))
    } else {
        if (!file.exists(paste0(dir_path, "/", sample_sheet))) {
            stop(paste0("====> Error: ", sample_sheet, " doesn't exist in the directory!"))
        }
    }

    # read sample sheet and check format
    qc_samplesheet <<- read.table(paste0(dir_path, "/", sample_sheet), sep = "\t", comment.char = "#", header = TRUE, fill = TRUE)
    targeton_name  <<- unique(qc_samplesheet$targeton_id)[1]
    message(sprintf("|--> Targeton: %s", targeton_name))
    require_cols <- c("sample_name",
                      "replicate",
                      "condition",
                      "ref_time_point",
                      "library_independent_count",
                      "library_dependent_count",
                      "valiant_meta",
                      "vep_anno",
                      "adapt5",
                      "adapt3",
                      "per_r1_adaptor",
                      "per_r2_adaptor",
                      "library_name",
                      "library_type")
    for (s in require_cols) {
        if (s %nin% colnames(qc_samplesheet)) {
            stop(paste0("====> Error: ", s, " must be in the sample sheet as the header"))
        }
    }

    if (length(unique(qc_samplesheet$sample_name)) != nrow(qc_samplesheet)) {
        stop(paste0("====> Error: ", sample_sheet, " has duplicated sample names!"))
    }

    maveqc_ref_time_point <<- unique(qc_samplesheet$ref_time_point)
    if (is.null(maveqc_ref_time_point) || is.na(maveqc_ref_time_point) || maveqc_ref_time_point == "") {
        maveqc_ref_time_point_samples <<- "NoRef"
    } else {
        if (maveqc_ref_time_point %in% qc_samplesheet$condition) {
            if (length(maveqc_ref_time_point) == 1) {
                maveqc_ref_time_point_samples <<- qc_samplesheet[qc_samplesheet$condition == maveqc_ref_time_point, ]$sample_name
            } else {
                stop(paste0("====> Error: ", sample_sheet, " has duplicated ref_time_point! It must be only one time point!"))
            }
        } else {
            stop(paste0("====> Error: ref_time_point cannot be found in condition of ", sample_sheet))
        }
    }

    maveqc_deseq_coldata <<- qc_samplesheet[, c("replicate", "condition")]
    rownames(maveqc_deseq_coldata) <<- qc_samplesheet$sample_name

    list_objects <- list()
    cat("Importing files for samples:", "\n", sep = "")
    for (i in 1:nrow(qc_samplesheet)) {
        cat("    |--> ", qc_samplesheet[i, ]$sample_name, "\n", sep = "")

        # leave the access in case user provides vep anno
        if (is.null(qc_samplesheet[i, ]$vep_anno)) {
            file_vep_anno <- NULL
        } else {
            if (is.na(qc_samplesheet[i, ]$vep_anno)) {
                file_vep_anno <- NULL
            } else {
                file_vep_anno <- paste0(dir_path, "/", qc_samplesheet[i, ]$vep_anno)
            }
        }

        tmp_obj <- create_sge_object(file_libcount = paste0(dir_path, "/", qc_samplesheet[i, ]$library_dependent_count),
                                     file_allcount = paste0(dir_path, "/", qc_samplesheet[i, ]$library_independent_count),
                                     file_valiant_meta = paste0(dir_path, "/", qc_samplesheet[i, ]$valiant_meta),
                                     file_vep_anno = file_vep_anno,
                                     file_libcount_hline = file_libcount_hline,
                                     file_allcount_hline = file_allcount_hline,
                                     file_valiant_meta_hline = file_valiant_meta_hline,
                                     file_vep_anno_hline = file_vep_anno_hline,
                                     file_libcount_cols = file_libcount_cols,
                                     file_allcount_cols = file_allcount_cols,
                                     file_valiant_meta_cols = file_valiant_meta_cols,
                                     file_vep_anno_cols = file_vep_anno_cols)
        tmp_obj@sample <- qc_samplesheet[i, ]$sample_name
        tmp_obj@libname <- qc_samplesheet[i, ]$library_name
        tmp_obj@libtype <- qc_samplesheet[i, ]$library_type

        tmp_obj@adapt5 <- ifelse(is.na(qc_samplesheet[i, ]$adapt5), "", qc_samplesheet[i, ]$adapt5)
        tmp_obj@adapt3 <- ifelse(is.na(qc_samplesheet[i, ]$adapt3), "", qc_samplesheet[i, ]$adapt3)
        tmp_obj@per_r1_adaptor <- ifelse(is.na(qc_samplesheet[i, ]$per_r1_adaptor), 0, qc_samplesheet[i, ]$per_r1_adaptor)
        tmp_obj@per_r2_adaptor <- ifelse(is.na(qc_samplesheet[i, ]$per_r2_adaptor), 0, qc_samplesheet[i, ]$per_r2_adaptor)

        tmp_obj <- format_count(tmp_obj)
        tmp_obj <- sge_stats(tmp_obj)
        tmp_obj <- sge_qc_stats(tmp_obj)

        list_objects <- append(list_objects, tmp_obj)
    }

    return(list_objects)
}

# - display.R
#' show basic info of a SGE object
#'
#' @export
#' @name show
#' @param object SGE object
setMethod(
    "show",
    signature = "SGE",
    definition = function(object) {
        cat("An object of class ", class(object), "\n", sep = "")
        cat("|--> sample name: ", object@sample, "\n", sep = "")
        cat("|--> library type: ", object@libtype, "\n", sep = "")
        cat("|--> library name: ", object@libname, "\n", sep = "")
        cat("    |--> 5' adaptor: ", object@adapt5, "\n", sep = "")
        cat("    |--> 3' adaptor: ", object@adapt3, "\n", sep = "")
        cat("    |--> ref seq: ", object@refseq, "\n", sep = "")
        cat("    |--> pam seq: ", object@pamseq, "\n", sep = "")
        cat("    |--> No. of library-dependent sequences: ", nrow(object@libcounts), "\n", sep = "")
        cat("    |--> No. of library-independent sequences: ", nrow(object@allcounts), "\n", sep = "")
        cat("|--> valiant meta: ", nrow(object@valiant_meta), " records and ", ncol(object@valiant_meta), " fields", "\n", sep = "")
        cat("    |--> ", sum(object@libcounts$name %in% object@valiant_meta$oligo_name), " library-dependent sequence IDs matched in valiant meta oligo names", "\n", sep = "")
    }
)

#' initialize function
setGeneric("show_stats", function(object, ...) {
  standardGeneric("show_stats")
})

#' show basic stats of a SGE object
#'
#' @export
#' @name show_stats
#' @param object SGE object
setMethod(
    "show_stats",
    signature = "SGE",
    definition = function(object) {
        colstrs <- colnames(object@libstats)
        dash_line <- paste0("|", strrep("-", 20 + nchar("type")), "-|-")
        dash_line <- paste0(dash_line, strrep("-", nchar("library dependent counts")), "-|-")
        dash_line <- paste0(dash_line, strrep("-", nchar("library independent counts")), "-|")
        cat("Basic stats of sample: ", object@sample, "\n", sep = "")
        cat(dash_line, "\n", sep = "")
        header_line <- paste0("|", strrep(" ", 20), "type | library dependent counts | library independent counts |")
        cat(header_line, "\n", sep = "")
        cat(dash_line, "\n", sep = "")
        for (i in 1:length(colstrs)) {
            spacestr <- strrep(" ", 20 - nchar(colstrs[i]) + nchar("type"))
            info_line <- paste0("|", spacestr, colstrs[i], " | ")
            spacestr <- strrep(" ", nchar("library dependent counts") - nchar(toString(object@libstats[, i])))
            info_line <- paste0(info_line, spacestr, object@libstats[, i], " | ")
            spacestr <- strrep(" ", nchar("library independent counts") - nchar(toString(object@allstats[, i])))
            info_line <- paste0(info_line, spacestr, object@allstats[, i], " | ")
            cat(info_line, "\n", sep = "")
        }
        cat(dash_line, "\n", sep = "")
    }
)

#' initialize function
setGeneric("show_qc_stats", function(object, ...) {
  standardGeneric("show_qc_stats")
})

#' show qc stats of a SGE object
#'
#' @export
#' @name show_qc_stats
#' @param object SGE object
setMethod(
    "show_qc_stats",
    signature = "SGE",
    definition = function(object) {
        colstrs <- colnames(object@libstats_qc)
        dash_line <- paste0("|", strrep("-", 20 + nchar("type")), "-|-")
        dash_line <- paste0(dash_line, strrep("-", nchar("library dependent counts")), "-|-")
        dash_line <- paste0(dash_line, strrep("-", nchar("library independent counts")), "-|")
        cat("QC stats of sample: ", object@sample, "\n", sep = "")
        cat(dash_line, "\n", sep = "")
        header_line <- paste0("|", strrep(" ", 20), "type | library dependent counts | library independent counts |")
        cat(header_line, "\n", sep = "")
        cat(dash_line, "\n", sep = "")
        for (i in 1:length(colstrs)) {
            spacestr <- strrep(" ", 20 - nchar(colstrs[i]) + nchar("type"))
            info_line <- paste0("|", spacestr, colstrs[i], " | ")
            spacestr <- strrep(" ", nchar("library dependent counts") - nchar(toString(object@libstats_qc[, i])))
            info_line <- paste0(info_line, spacestr, object@libstats_qc[, i], " | ")
            spacestr <- strrep(" ", nchar("library independent counts") - nchar(toString(object@allstats_qc[, i])))
            info_line <- paste0(info_line, spacestr, object@allstats_qc[, i], " | ")
            cat(info_line, "\n", sep = "")
        }
        cat(dash_line, "\n", sep = "")
    }
)

#' show basic info of a sample QC object
#'
#' @export
#' @name show
#' @param object sampleQC object
setMethod(
    "show",
    signature = "sampleQC",
    definition = function(object) {
        cat("An object of class ", class(object), "\n", sep = "")
        cat("|--> samples: ", "\n", sep = "")
        for (s in object@samples) {
            cat("    |--> ", s@sample, "\n", sep = "")
        }

        cat("|--> reference samples: ", "\n", sep = "")
        if (length(object@samples_ref) == 0) {
            cat("    |--> no sample found", "\n", sep = "")
        } else {
            for (s in object@samples_ref) {
                cat("    |--> ", s@sample, "\n", sep = "")
            }
        }

        cat("|--> QC results: ", "\n", sep = "")
        for (i in 1:nrow(object@stats)) {
            cat("    |--> ", rownames(object@stats)[i], ": ", object@stats[i, ]$qcpass, "\n", sep = "")
        }
    }
)

#' show basic info of a experiment QC object
#'
#' @export
#' @name show
#' @param object experimentQC object
setMethod(
    "show",
    signature = "experimentQC",
    definition = function(object) {
        cat("An object of class ", class(object), "\n", sep = "")
        cat("|--> samples: ", "\n", sep = "")
        for (s in object@samples) {
            cat("    |--> ", s@sample, "\n", sep = "")
        }

        cat("|--> reference condition: ", object@ref_condition, "\n", sep = "")

        cat("|--> DESeq coldata: ", "\n", sep = "")
        tmp_coldata <- as.matrix(object@coldata)
        for (i in 1:nrow(object@coldata)) {
            cat("    |--> ", rownames(tmp_coldata)[i], "\t", tmp_coldata[i, 1], "\t", tmp_coldata[i, 2], "\n", sep = "")
        }
    }
)

#' transparent color function
#'
#' creating a transparent color using color name,
#' alpha rate is 0 to 1
#'
#' @name t_col
#' @param col  color name
#' @param rate alpha rate
#' @return transparent color
t_col <- function(col, rate) {
    newcol <- rgb(col2rgb(col)["red", ],
                  col2rgb(col)["green", ],
                  col2rgb(col)["blue", ],
                  as.integer(rate * 255),
                  maxColorValue = 255)
    return(newcol)
}

#' capitalise the first character in the string
#'
#' creating the capitalised names for data frame
#'
#' @name capital_names
#' @param x  a vector of strings
#' @return capitalised strings
capital_names <- function(x) {
    y <- x
    for (i in 1:length(x)) {
        y[i] <- paste(toupper(substring(x[i], 1, 1)), substring(x[i], 2), sep = "", collapse = " ")
    }
    return(y)
}

#' not in function
#'
#' @name %nin%
#' @param x X
#' @param y Y
#' @return True or False
`%nin%` <- function(x, y) !(x %in% y)

#' reverse complement
#'
#' @name revcomp
#' @param seq sequence
#' @return string
revcomp <- function(seq) {
    seq <- toupper(seq)
    splits <- strsplit(seq, "")[[1]]
    reversed <- rev(splits)
    seq_rev <- paste(reversed, collapse = "")
    seq_rev_comp <- chartr("ATCG", "TAGC", seq_rev)
    return(seq_rev_comp)
}

#' trim adaptor sequences
#'
#' @name trim_adaptor
#' @param seq    sequence
#' @param adapt5 5 prime adaptor sequence
#' @param adapt3 3 prime adaptor sequence
#' @return string
trim_adaptor <- function(seq, adapt5, adapt3) {
    adapt5_pos <- regexpr(adapt5, seq, fixed = TRUE)[1]
    adapt3_pos <- regexpr(adapt3, seq, fixed = TRUE)[1]

    is_revcomp <- FALSE
    # ? could adaptor revcomp ?
    if (adapt5_pos < 0 & adapt3_pos < 0) {
        adapt5_revcomp <- revcomp(adapt5)
        adapt3_revcomp <- revcomp(adapt3)

        adapt5_revcomp_pos <- regexpr(adapt5_revcomp, seq, fixed = TRUE)[1]
        adapt3_revcomp_pos <- regexpr(adapt3_revcomp, seq, fixed = TRUE)[1]

        if (adapt3_revcomp_pos < 0 & adapt5_revcomp_pos < 0) {
            return(seq)
        } else {
            is_revcomp <- TRUE
        }
    }

    if (is_revcomp == FALSE) {
        if (adapt3_pos > adapt5_pos) {
            if (adapt5_pos > 0 & adapt3_pos > 0) {
                return(substr(seq, adapt5_pos + nchar(adapt5), adapt3_pos - 1))
            } else if (adapt5_pos > 0 & adapt3_pos < 0) {
                return(substr(seq, adapt5_pos + nchar(adapt5), nchar(seq)))
            } else if (adapt5_pos < 0 & adapt3_pos > 0) {
                return(substr(seq, 1, adapt3_pos - 1))
            }
        } else {
            stop(paste0("====> Error: 3 prime adaptor found before 5 prime adaptor in the sequence: ", seq))
        }
    } else {
        if (adapt5_revcomp_pos > adapt3_revcomp_pos) {
            if (adapt3_revcomp_pos > 0 & adapt5_revcomp_pos > 0) {
                return(substr(seq, adapt3_revcomp_pos + nchar(adapt3_revcomp), adapt5_revcomp_pos - 1))
            } else if (adapt3_revcomp_pos > 0 & adapt5_revcomp_pos < 0) {
                return(substr(seq, adapt3_revcomp_pos + nchar(adapt3_revcomp), nchar(seq)))
            } else if (adapt3_revcomp_pos < 0 & adapt5_revcomp_pos > 0) {
                return(substr(seq, 1, adapt3_revcomp_pos - 1))
            }
        } else {
            stop(paste0("====> Error: 5 prime adaptor (RC) found before 3 prime adaptor (RC) in the sequence: ", seq))
        }
    }
}

#' column binding with filling NAs
#'
#' @name cbind_fill
#' @return matrix
cbind_fill <- function(...) {
    nm <- list(...)
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow))
    do.call(cbind, lapply(nm, function (x) rbind(x, matrix(, n - nrow(x), ncol(x)))))
}

#' calculate gini coefficiency for a sample
#'
#' @name cal_gini
#' @param x a vector
#' @return a value
cal_gini <- function(x, corr = FALSE, na.rm = TRUE) {
    if (!na.rm && any(is.na(x))) return(NA_real_)
    x <- as.numeric(na.omit(x))
    n <- length(x)
    x <- sort(x)
    G <- sum(x * 1L:n)
    G <- 2 * G/sum(x) - (n + 1L)

    if (corr) {
        return(G / (n - 1L))
    } else {
        return(G / n)
    }
}

#' merge a list of data tables into a data table
#'
#' @name merge_list_to_dt
#' @param list_dt   a list of data tables
#' @param by_val    join data tables by which column
#' @param join_val  join which column in the data tables
#' @return a data table
merge_list_to_dt <- function(list_dt, by_val, join_val) {
    dt_out <- data.table()

    for (i in 1:length(list_dt)) {
        cols <- c(by_val, join_val)
        dt_tmp <- list_dt[[i]][, ..cols]

        if (nrow(dt_out) == 0) {
            dt_out <- dt_tmp
            colnames(dt_out) <- c(by_val, names(list_dt)[i])
        } else {
            coln <- colnames(dt_out)
            dt_out <- merge(dt_out, dt_tmp, by = by_val, all = TRUE)
            colnames(dt_out) <- c(coln, names(list_dt)[i])
        }
    }

    return(dt_out)
}

#' color blind friendly
#'
#' @name select_colorblind
#' @param col_id a character to select colors
#' @return a vector of colors
select_colorblind <- function(col_id) {
    col8 <- c("#D55E00", "#56B4E9", "#E69F00",
              "#009E73", "#F0E442", "#0072B2",
              "#CC79A7", "#000000")

    col12 <- c("#88CCEE", "#CC6677", "#DDCC77",
               "#117733", "#332288", "#AA4499",
               "#44AA99", "#999933", "#882255",
               "#661100", "#6699CC", "#888888")

    col15 <- c("red",       "royalblue", "olivedrab",
               "purple",    "violet",    "maroon1",
               "seagreen1", "navy",      "pink",
               "coral",     "steelblue", "turquoise1",
               "red4",      "skyblue",   "yellowgreen")

    col21 <- c("#F60239", "#009503", "#FFDC3D",
               "#9900E6", "#009FFA", "#FF92FD",
               "#65019F", "#FF6E3A", "#005A01",
               "#00E5F8", "#DA00FD", "#AFFF2A",
               "#00F407", "#00489E", "#0079FA",
               "#560133", "#EF0096", "#000000",
               "#005745", "#00AF8E", "#00EBC1")

    if (col_id == "col8") {
        return(col8)
    } else if (col_id == "col12") {
        return(col12)
    } else if (col_id == "col15") {
        return(col15)
    } else if (col_id == "col21") {
        return(col21)
    } else {
        stop(paste0("====> Error: wrong col_id"))
    }
}

#' fetch objects from list by the names or indexes
#'
#' @export
#' @name select_objects
#' @param objects a list of objects
#' @param tags    a vector of names or indexes
#' @return a list of objects
select_objects <- function(objects, tags) {
    if (length(objects) == 0) {
        stop(paste0("====> Error: no object found in the list"))
    }

    if (length(tags) == 0) {
        stop(paste0("====> Error: please provide tags to fetch objects"))
    } else {
        if (class(tags) %in% c("character", "numeric")) {
            if (class(tags) == "numeric") {
                tags <- as.integer(tags)
            }
        } else {
            stop(paste0("====> Error: wrong tag type, must be integer or character"))
        }
    }

    list_select <- list()
    if (class(tags) == "character") {
        for (i in 1:length(tags)) {
            for (j in 1:length(objects)) {
                if (tags[i] == objects[[j]]@sample) {
                    list_select <- append(list_select, objects[[j]])
                    break
                }
            }
        }
    } else {
        for (i in 1:length(tags)) {
            list_select <- append(list_select, objects[[tags[i]]])
        }
    }

    return(list_select)
}

#' initialize function
setGeneric("run_experiment_qc", function(object, ...) {
  standardGeneric("run_experiment_qc")
})

#' run DESeq2 for the list of samples
#'
#' @export
#' @name run_experiment_qc
#' @param object  experimentQC object
#' @param pcut    the padj cutoff
#' @param dcut    the depleted log2 fold change cutoff
#' @param ecut    the enriched log2 fold change cutoff
#' @param ntop    the number of top variances
#' @return object
setMethod(
    "run_experiment_qc",
    signature = "experimentQC",
    definition = function(object,
                          pcut = maveqc_config$expqc_padj,
                          dcut = maveqc_config$expqc_lfc_depleted,
                          ecut = maveqc_config$expqc_lfc_enriched,
                          ntop = maveqc_config$expqc_top_variants) {
        cat("Running DESeq2 on library counts for PCA...", "\n", sep = "")
        object <- run_experiment_qc_lib_lfc(object, pcut = pcut, dcut = dcut, ecut = ecut, ntop = ntop)

        cat("Running DESeq2 on all counts after filtering...", "\n", sep = "")
        object <- run_experiment_qc_all_lfc(object, pcut = pcut, dcut = dcut, ecut = ecut)
        
        return(object)
    }
)


#' initialize function
setGeneric("run_experiment_qc_lib_lfc", function(object, ...) {
  standardGeneric("run_experiment_qc_lib_lfc")
})

#' run DESeq2 for the list of samples
#'
#' @export
#' @name run_experiment_qc_lib_lfc
#' @param object  experimentQC object
#' @param pcut    the padj cutoff
#' @param dcut    the depleted log2 fold change cutoff
#' @param ecut    the enriched log2 fold change cutoff
#' @param ntop    the number of top variances
#' @return object
setMethod(
    "run_experiment_qc_lib_lfc",
    signature = "experimentQC",
    definition = function(object,
                          pcut,
                          dcut,
                          ecut,
                          ntop) {
        #----------------------------#
        # 1. calculating size factor #
        #----------------------------#
        # run control
        cat("    |--> Running control DESeq2 to get size factor...", "\n", sep = "")

        library_counts_anno <- as.data.frame(object@library_counts_anno)
        library_counts_anno[is.na(library_counts_anno)] <- 0
        ds_coldata <- object@coldata

		# rownames are necessary for DESeq2, otherwise error happens
		rownames(library_counts_anno) <- library_counts_anno$sequence

		# --------------------------
		# Choose control set
		# --------------------------
		syn_counts <- library_counts_anno[library_counts_anno$consequence == "Synonymous_Variant", rownames(ds_coldata)]

		if (nrow(syn_counts) > 0) {
			message("    |--> Using synonymous variants for normalization")
			control_counts <- syn_counts
		} else {
			message("    |--> No synonymous variants found, using SNVs for normalization")
			control_counts <- library_counts_anno[library_counts_anno$consequence == "SNV", rownames(ds_coldata)]
		}

		# --------------------------
		# Build DESeq2 object for control counts
		# --------------------------
		suppressMessages(control_ds_obj <- DESeqDataSetFromMatrix(
			countData = control_counts,
			colData   = ds_coldata,
			design    = ~condition
		))
		control_ds_obj <- control_ds_obj[rowSums(counts(control_ds_obj)) > 0, ]
		control_ds_obj$condition <- factor(control_ds_obj$condition, levels = mixedsort(levels(control_ds_obj$condition)))
		control_ds_obj$condition <- relevel(control_ds_obj$condition, ref = object@ref_condition)
		control_ds_obj <- estimateSizeFactors(control_ds_obj)

		# --------------------------
		# Apply size factors to all counts
		# --------------------------
		cat("    |--> Applying size factor to get DESeq2 normalised counts...", "\n", sep = "")
		deseq_counts <- library_counts_anno[, rownames(ds_coldata)]

		suppressMessages(ds_obj <- DESeqDataSetFromMatrix(
			countData = deseq_counts,
			colData   = ds_coldata,
			design    = ~condition
		))
		ds_obj <- ds_obj[rowSums(counts(ds_obj)) > 0, ]
		ds_obj$condition <- factor(ds_obj$condition, levels = mixedsort(levels(ds_obj$condition)))
		ds_obj$condition <- relevel(ds_obj$condition, ref = object@ref_condition)
		sizeFactors(ds_obj) <- sizeFactors(control_ds_obj)

		suppressMessages(ds_obj <- DESeq(ds_obj, quiet = TRUE))
		ds_rlog <- rlog(ds_obj)

		object@lib_deseq_rlog <- as.data.frame(assay(ds_rlog))


        #-----------------------#
        # 2. clustering and PCA #
        #-----------------------#
        cat("    |--> Clustering and PCA...", "\n", sep = "")

        sample_dist <- dist(t(object@lib_deseq_rlog), method = "euclidean")
        sample_hclust <- hclust(d = sample_dist, method = "ward.D2")

        object@lib_hclust_res <- sample_hclust
        object@lib_corr_res <- cor(scale(as.matrix(object@lib_deseq_rlog)))

        pca_input <- as.matrix(object@lib_deseq_rlog)
        rv <- rowVars(pca_input)
        select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        pca <- prcomp(t(pca_input[select, ]), center = TRUE, scale = TRUE)

        object@lib_pca_res <- pca

        #-------------------#
        # 3. DESeq2 results #
        #-------------------#
        cat("    |--> Calculating DESeq2 LFC...", "\n", sep = "")

        suppressMessages(ds_res <- degComps(ds_obj,
                                            combs = "condition",
                                            contrast = object@comparisons,
                                            alpha = 0.05,
                                            skip = FALSE,
                                            type = "apeglm",
                                            pairs = FALSE,
                                            fdr = "default"))

        object@lib_deseq_res <- ds_res

        library_counts_pos_anno <- object@library_counts_pos_anno

        comparisons <- names(object@lib_deseq_res)
        for (i in seq_along(object@lib_deseq_res)) {
            # get raw & shrunk results (if available)
            res_raw <- object@lib_deseq_res[[i]]$raw[, c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"), drop = FALSE]
            colnames(res_raw) <- paste0(colnames(res_raw), "_raw")
            res_raw$sequence <- rownames(object@lib_deseq_res[[i]]$raw)

            res_shr <- object@lib_deseq_res[[i]]$shrunken[, c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"), drop = FALSE]
            colnames(res_shr) <- paste0(colnames(res_shr), "_shrunk")
            res_shr$sequence <- rownames(object@lib_deseq_res[[i]]$shrunken)

            # merge raw and shrunk (inner join to avoid rows only present on one side)
            res <- merge(as.data.table(res_raw), as.data.table(res_shr), by = "sequence", all = FALSE)

            # annotate
            res[library_counts_pos_anno, c("oligo_name", "position", "consequence") := .(oligo_name, position, consequence), on = .(sequence)]

            # Keep only library-designed sequences
            n_before <- nrow(res)
            res <- res[sequence %in% library_counts_pos_anno$sequence]
            message(sprintf("    |--> lib: filtered to library variants: %d -> %d rows", n_before, nrow(res)))

            # compute stat based on the user's thresholds for raw and shrunk separately
            res[, stat_raw := "no impact"]
            res[(padj_raw < pcut) & (log2FoldChange_raw > ecut), stat_raw := "enriched"]
            res[(padj_raw < pcut) & (log2FoldChange_raw < dcut), stat_raw := "depleted"]
            res[, stat_raw := factor(stat_raw, levels = c("no impact", "enriched", "depleted"))]

            res[, stat_shrunk := "no impact"]
            res[(padj_raw < pcut) & (log2FoldChange_shrunk > ecut), stat_shrunk := "enriched"]
            res[(padj_raw < pcut) & (log2FoldChange_shrunk < dcut), stat_shrunk := "depleted"]
            res[, stat_shrunk := factor(stat_shrunk, levels = c("no impact", "enriched", "depleted"))]

            # drop rows lacking consequence annotation 
            res <- na.omit(res, cols = "consequence")

            # desired column order for exp_lib (robust to missing columns)
            desired_order_lib <- c(
              "sequence", "oligo_name", "position", "consequence",
              "baseMean_raw", "log2FoldChange_raw", "lfcSE_raw", "pvalue_raw", "padj_raw", "stat_raw",
              "baseMean_shrunk", "log2FoldChange_shrunk", "lfcSE_shrunk", "pvalue_shrunk", "padj_shrunk", "stat_shrunk"
            )
            present_lib <- intersect(desired_order_lib, names(res))
            setcolorder(res, c(present_lib, setdiff(names(res), present_lib)))

            object@lib_deseq_res_anno[[comparisons[i]]] <- res
        }


        return(object)
    }
)

# -------------------------------
# Positional bias correction helpers
# -------------------------------

fit_positional_loess_fullrange <- function(df_fit,
                                           all_positions,
                                           span = 0.2) {
  # df_fit: data.frame with columns position (numeric/integer), lfc (numeric)
  # all_positions: integer/numeric vector for the full library coordinate range

  df_fit <- df_fit[!is.na(df_fit$position) & !is.na(df_fit$lfc), ]
  df_fit$position <- as.numeric(df_fit$position)
  df_fit$lfc <- as.numeric(df_fit$lfc)

  if (nrow(df_fit) < 10) {
    stop("Too few variants to fit LOESS positional model.")
  }

  lo <- loess(lfc ~ position, data = df_fit, span = span)

  # Predict on full coordinate range (unique sorted)
  pred_pos <- sort(unique(as.numeric(all_positions)))

  pred <- predict(
    lo,
    newdata = data.frame(position = pred_pos),
    se = TRUE
  )

  fit_vals <- as.numeric(pred$fit)
  fit_se   <- as.numeric(pred$se.fit)

  # Flat extension for any NA predictions: fill from nearest non-NA
  ok <- which(!is.na(fit_vals))
  if (length(ok) == 0) stop("LOESS returned only NA predictions.")
  first_ok <- ok[1]
  last_ok  <- ok[length(ok)]

  if (first_ok > 1) {
    fit_vals[1:(first_ok - 1)] <- fit_vals[first_ok]
    fit_se[1:(first_ok - 1)]   <- fit_se[first_ok]
  }
  if (last_ok < length(fit_vals)) {
    fit_vals[(last_ok + 1):length(fit_vals)] <- fit_vals[last_ok]
    fit_se[(last_ok + 1):length(fit_vals)]   <- fit_se[last_ok]
  }

  data.table(position = pred_pos,
             pos_fit = fit_vals,
             pos_fit_se = fit_se)
}

apply_positional_correction <- function(res_dt, fit_dt, lfc_cols) {
  # res_dt: data.table with 'position' column and LFC columns
  # fit_dt: data.table(position, pos_fit, pos_fit_se)
  # lfc_cols: character vector of columns in res_dt to adjust

  x <- as.data.table(copy(res_dt))

  # Ensure numeric position match
  x[, position := as.numeric(position)]
  fit_dt <- fit_dt[, .(
    position = as.numeric(position),
    pos_fit,
    pos_fit_se
  )]

  # Join fit info onto result rows
  x[fit_dt, `:=`(
    pos_fit    = i.pos_fit,
    pos_fit_se = i.pos_fit_se
  ), on = .(position)]

  # safeguard for any missing
  x[is.na(pos_fit), `:=`(pos_fit = 0, pos_fit_se = 0)]

  # Apply correction to requested LFC columns
  for (col in lfc_cols) {
    new_col <- paste0("pos_adj_", col)
    x[, (new_col) := get(col) - pos_fit]
  }

  x
}

recenter_by_controls <- function(dt_res,
                                 lfc_col,
                                 control_consequences = c("Synonymous_Variant", "Intronic_Variant"),
                                 fallback = "SNV") {
  ctrl <- dt_res[consequence %in% control_consequences & !is.na(get(lfc_col))]
  if (nrow(ctrl) == 0) {
    message("    |--> No synonymous/intronic controls for re-centering; using SNVs instead")
    ctrl <- dt_res[consequence == fallback & !is.na(get(lfc_col))]
  }
  if (nrow(ctrl) == 0) {
    warning(paste0("No controls available for re-centering column ", lfc_col, "; leaving unchanged"))
    return(dt_res)
  }

  med <- median(ctrl[[lfc_col]], na.rm = TRUE)
  dt_res[, (lfc_col) := get(lfc_col) - med]
  dt_res
}

#' initialize function
setGeneric("run_experiment_qc_all_lfc", function(object, ...) {
  standardGeneric("run_experiment_qc_all_lfc")
})
#' run DESeq2 for the list of samples
#'
#' @export
#' @name run_experiment_qc_all_lfc
#' @param object  experimentQC object
#' @param pcut    the padj cutoff
#' @param dcut    the depleted log2 fold change cutoff
#' @param ecut    the enriched log2 fold change cutoff
#' @return object
setMethod(
    "run_experiment_qc_all_lfc",
    signature = "experimentQC",
    definition = function(object,
                          pcut,
                          dcut,
                          ecut) {
        #----------------------------#
        # 1. calculating size factor #
        #----------------------------#
        # run control
        cat("    |--> Running normalisation using total number of counts...", "\n", sep = "")

        deseq_counts <- as.data.frame(object@accepted_counts)
        deseq_counts[is.na(deseq_counts)] <- 0
        ds_coldata <- object@coldata

        rownames(deseq_counts) <- deseq_counts$sequence
        deseq_counts <- deseq_counts[, rownames(ds_coldata)]

        sizeFactor <- colSums(deseq_counts) / 1000000
        normFactor <- t(replicate(nrow(deseq_counts), sizeFactor))

        suppressMessages(ds_obj <- DESeqDataSetFromMatrix(countData = deseq_counts, colData = ds_coldata, design = ~condition))
        ds_obj$condition <- factor(ds_obj$condition, levels = mixedsort(levels(ds_obj$condition)))
        ds_obj$condition <- relevel(ds_obj$condition, ref = object@ref_condition)
        normalizationFactors(ds_obj) <- normFactor[, rownames(ds_coldata)]
        suppressMessages(ds_obj <- DESeq(ds_obj, quiet = TRUE))
        ds_rlog <- rlog(ds_obj)

        object@all_deseq_rlog <- as.data.frame(assay(ds_rlog))

        #---------------------------#
        # 2. calculating DESeq2 LFC #
        #---------------------------#
        cat("    |--> Calculating DESeq2 LFC...", "\n", sep = "")

        # Get both raw and shrunk results
        suppressMessages(ds_res <- degComps(ds_obj,
                                          combs = "condition",
                                          contrast = object@comparisons,
                                          alpha = 0.05,
                                          skip = FALSE,
                                          type = "apeglm",
                                          pairs = FALSE,
                                          fdr = "default"))

        object@all_deseq_res <- ds_res

		vep_anno <- as.data.table(object@vep_anno)[, .(
		  oligo_name = unique_oligo_name,
		  sequence = seq,
		  position = vcf_pos,
		  consequence = summary_plot
		)]

        # ----- build annotated merged tables (raw + shrunk) -----
        contrast_names_raw <- names(object@all_deseq_res)
        for (i in seq_along(object@all_deseq_res)) {
            # Get raw results
            res_raw <- object@all_deseq_res[[i]]$raw[, c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"), drop = FALSE]
            colnames(res_raw) <- paste0(colnames(res_raw), "_raw")
            res_raw$sequence <- rownames(object@all_deseq_res[[i]]$raw)

            # Get shrunk results 
            res_shrunk <- object@all_deseq_res[[i]]$shrunken[, c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"), drop = FALSE]
            colnames(res_shrunk) <- paste0(colnames(res_shrunk), "_shrunk")
            res_shrunk$sequence <- rownames(object@all_deseq_res[[i]]$shrunken)

            # Merge raw and shrunk results (inner join)
            res <- merge(as.data.table(res_raw), as.data.table(res_shrunk), by = "sequence", all = FALSE)

            # annotate with vep info
            res[vep_anno, c("oligo_name", "position", "consequence") := .(oligo_name, position, consequence), on = .(sequence)]

            # Keep only library-designed sequences
            n_before <- nrow(res)
            res <- res[sequence %in% vep_anno$sequence]
            message(sprintf("    |--> all: filtered to designed variants: %d -> %d rows", n_before, nrow(res)))

            object@all_deseq_res_anno[[contrast_names_raw[i]]] <- res
        }

		#-----------------------------------#
		# 3. adjusting DESeq2 LFC & p value #
		#-----------------------------------#
        cat("    |--> Adjusting DESeq2 LFC & p value (computing adj_* for positional stats)...\n")

        control_consequences <- c("Synonymous_Variant", "Intronic_Variant")
        contrast_names_anno <- names(object@all_deseq_res_anno)

        for (nm in contrast_names_anno) {
            res <- copy(object@all_deseq_res_anno[[nm]])  # work on a copy

            control_res <- res[consequence %in% control_consequences]
            if (nrow(control_res) == 0) {
                message("No synonymous/intronic controls found, using SNVs instead")
                control_res <- res[consequence == "SNV"]
            }

            # compute medians on the raw and shrunk log2FCs
            control_median_lfc_raw <- median(control_res$log2FoldChange_raw, na.rm = TRUE)
            control_median_lfc_shrunk <- median(control_res$log2FoldChange_shrunk, na.rm = TRUE)

            # adjust
            res[, adj_log2FoldChange_raw := log2FoldChange_raw - control_median_lfc_raw]
            res[, adj_log2FoldChange_shrunk := log2FoldChange_shrunk - control_median_lfc_shrunk]

            # adjusted scores / p-values (both raw and shrunk)
            res[, adj_score_raw := adj_log2FoldChange_raw / lfcSE_raw]
            res[, adj_score_shrunk := adj_log2FoldChange_shrunk / lfcSE_shrunk]

            res[, adj_pval_raw := pnorm(abs(adj_score_raw), lower.tail = FALSE) * 2]
            res[, adj_pval_shrunk := pnorm(abs(adj_score_shrunk), lower.tail = FALSE) * 2]

            res[, adj_fdr_raw := p.adjust(adj_pval_raw, method = "fdr")]
            res[, adj_fdr_shrunk := p.adjust(adj_pval_shrunk, method = "fdr")]
            
            res[, stat_adj_raw := "no impact"]
			res[(adj_fdr_raw < pcut) & (adj_log2FoldChange_raw > ecut), stat_adj_raw := "enriched"]
			res[(adj_fdr_raw < pcut) & (adj_log2FoldChange_raw < dcut), stat_adj_raw := "depleted"]
			res[, stat_adj_raw := factor(stat_adj_raw, levels = c("no impact", "enriched", "depleted"))]

			res[, stat_adj_shrunk := "no impact"]
			res[(adj_fdr_raw < pcut) & (adj_log2FoldChange_shrunk > ecut), stat_adj_shrunk := "enriched"]
			res[(adj_fdr_raw < pcut) & (adj_log2FoldChange_shrunk < dcut), stat_adj_shrunk := "depleted"]
			res[, stat_adj_shrunk := factor(stat_adj_shrunk, levels = c("no impact", "enriched", "depleted"))]
            
            # drop rows that lack consequence annotation
            res <- na.omit(res, cols = "consequence")

            object@all_deseq_res_anno_adj[[nm]] <- res
        }

		#-----------------------------------#
		# 3.5 Positional bias correction     #
		#-----------------------------------#
		cat("    |--> Positional bias correction (LOESS) using time series vs reference...\n")

		ref <- as.character(object@ref_condition)

		# Full library coordinate range
		all_positions <- as.numeric(vep_anno$position)
		all_positions <- all_positions[!is.na(all_positions)]
		if (length(all_positions) < 2) {
		  stop("Not enough vcf_pos values in vep_anno to define library coordinate range.")
		}

		# Fit from unadjusted tables...
		fit_list <- object@all_deseq_res_anno
		# ...apply correction onto adjusted tables for output consistency
		res_list <- object@all_deseq_res_anno_adj

		# Helper: require a contrast to exist
		get_contrast_fit <- function(name) {
		  if (!name %in% names(fit_list)) stop(paste0("Missing contrast in all_deseq_res_anno: ", name))
		  fit_list[[name]]
		}

		# Dynamically discover contrasts for this reference, sorted by numeric day value.
		# This handles any set of timepoints (e.g. Day4/Day5/Day15 or Day4/Day7/Day15).
		all_contrast_names <- names(fit_list)
		pattern <- paste0("^condition_Day(\\d+)_vs_", ref, "$")
		matched <- grep(pattern, all_contrast_names, value = TRUE)
		if (length(matched) == 0) {
		  stop(paste0("Positional correction: no contrasts matching 'condition_Day*_vs_", ref,
		              "' found in all_deseq_res_anno. Available contrasts: ",
		              paste(all_contrast_names, collapse = ", ")))
		}
		# Sort contrasts by the numeric day number
		day_nums <- as.integer(sub(pattern, "\\1", matched))
		matched  <- matched[order(day_nums)]
		cat(paste0("    |--> Detected contrasts for positional correction (ref=", ref, "): ",
		           paste(matched, collapse = ", "), "\n"))

		# Assign to c1..c5 (NULL if fewer than 5 timepoints)
		c1 <- if (length(matched) >= 1) matched[1] else NULL
		c2 <- if (length(matched) >= 2) matched[2] else NULL
		c3 <- if (length(matched) >= 3) matched[3] else NULL
		c4 <- if (length(matched) >= 4) matched[4] else NULL
		c5 <- if (length(matched) >= 5) matched[5] else NULL

		# Pull contrasts (only the ones that exist)
		d1 <- if (!is.null(c1)) get_contrast_fit(c1) else NULL
		d2 <- if (!is.null(c2)) get_contrast_fit(c2) else NULL
		d3 <- if (!is.null(c3)) get_contrast_fit(c3) else NULL
		d4 <- if (!is.null(c4)) get_contrast_fit(c4) else NULL
		d5 <- if (!is.null(c5)) get_contrast_fit(c5) else NULL

		fit_consequences <- c("Synonymous_Variant", "Intronic_Variant", "Missense_Variant")
		message(sprintf("    |--> fit_consequences: %s", paste(fit_consequences, collapse = ", ")))
		lfc_raw_col    <- "log2FoldChange_raw"
		lfc_shrunk_col <- "log2FoldChange_shrunk"

		# Build joined table for filtering / fallback -- dynamically over however many
		# timepoints were detected (d1 always exists; d2..d5 may be NULL).
		active_d <- Filter(Negate(is.null), list(d1, d2, d3, d4, d5))
		dt_parts <- lapply(seq_along(active_d), function(i) {
		  raw_col <- paste0("lfc_t", i, "_raw")
		  shr_col <- paste0("lfc_t", i, "_shr")
		  dt_i <- as.data.table(active_d[[i]])
		  if (i == 1) {
		    out <- dt_i[, .(sequence, position, consequence,
		                    lfc_raw_tmp = get(lfc_raw_col),
		                    lfc_shr_tmp = get(lfc_shrunk_col))]
		  } else {
		    out <- dt_i[, .(sequence,
		                    lfc_raw_tmp = get(lfc_raw_col),
		                    lfc_shr_tmp = get(lfc_shrunk_col))]
		  }
		  setnames(out, "lfc_raw_tmp", raw_col)
		  setnames(out, "lfc_shr_tmp", shr_col)
		  out
		})
		dt0 <- Reduce(function(x, y) merge(x, y, by = "sequence", all = FALSE), dt_parts)

		dt0 <- as.data.table(dt0)
		dt0 <- dt0[!is.na(position) & consequence %in% fit_consequences]
		dt0[, position := as.numeric(position)]
		dt0 <- dt0[!is.na(position)]

		# -----------------------------
		# Strict fit-point selection (Findlay-style)
		# -----------------------------
		log2_08 <- log2(0.8)

		# Compute t1 threshold from the median of synonymous/intronic variants only,
		# regardless of whether missense is included in fit_consequences. This makes
		# the threshold distribution-relative (centre - 0.5) rather than fixed at
		# log2(0.5) = -1, so that clearly depleting variants are excluded even when
		# the bulk of the distribution sits near zero.
		syn_intr_only <- c("Synonymous_Variant", "Intronic_Variant")
		t1_centre <- median(
		  dt0[consequence %in% syn_intr_only & !is.na(lfc_t1_raw), lfc_t1_raw],
		  na.rm = TRUE
		)
		if (!is.finite(t1_centre)) {
		  message("    |--> Could not compute t1 centre from synonymous/intronic variants; falling back to log2(0.5) threshold.")
		  t1_centre <- log2(0.5) + 0.5  # preserves the -1 fallback
		}
		t1_threshold <- t1_centre - 0.5
		message(sprintf("    |--> t1 LFC centre (syn/intronic median): %.3f; t1 threshold: %.3f", t1_centre, t1_threshold))

		n_timepoints <- length(Filter(Negate(is.null), list(c1, c2, c3, c4, c5)))
		need_cols <- paste0("lfc_t", seq_len(n_timepoints), "_raw")

		dt_strict <- dt0[Reduce(`&`, lapply(need_cols, function(cc) !is.na(dt0[[cc]])))]

		# Exclude variants substantially depleted at t1 (Findlay-style), using a
		# distribution-relative threshold rather than a fixed value.
		if (n_timepoints >= 1) {
		  dt_strict <- dt_strict[get("lfc_t1_raw") > t1_threshold]
		}
		if (n_timepoints >= 2) {
		  for (ti in seq_len(n_timepoints - 1)) {
		    col_curr <- paste0("lfc_t", ti + 1, "_raw")
		    col_prev <- paste0("lfc_t", ti,     "_raw")
		    dt_strict <- dt_strict[(get(col_curr) - get(col_prev)) > log2_08]
		  }
		}

		# -----------------------------
		# Fallback trigger: do fit points span the library?
		# -----------------------------
		span_ok <- function(pos_vec, all_pos, min_frac = 0.70) {
		  pos_vec <- as.numeric(pos_vec)
		  pos_vec <- pos_vec[!is.na(pos_vec)]
		  if (length(pos_vec) < 10) return(FALSE)

		  rng_fit <- diff(range(pos_vec))
		  rng_all <- diff(range(all_pos))
		  if (!is.finite(rng_fit) || !is.finite(rng_all) || rng_all <= 0) return(FALSE)

		  (rng_fit / rng_all) >= min_frac
		}

		use_fallback <- FALSE

		# If strict set is too small OR doesn’t span positions, fallback to broader set
		if (nrow(dt_strict) < 10 || !span_ok(dt_strict$position, all_positions, min_frac = 0.70)) {
		  use_fallback <- TRUE
		  message("    |--> Strict fit points insufficient or not spanning positions; falling back to all variants in fit_consequences for LOESS fit.")
		}

		# Broader fit set: all variants in fit_consequences with a usable t1 LFC (and position)
		dt_fitpts <- if (!use_fallback) dt_strict else dt0[!is.na(lfc_t1_raw)]


		if (nrow(dt_fitpts) < 10) {
		  stop("Too few variants available to fit LOESS even after fallback.")
		}

		# Fit LOESS using the earliest post-ref contrast (t1) as function of position
		df_fit <- data.frame(position = dt_fitpts$position, lfc = dt_fitpts$lfc_t1_raw)
		fit_dt <- fit_positional_loess_fullrange(df_fit = df_fit, all_positions = all_positions, span = 0.2)

		# Store diagnostics
		res_list[["positional_loess_fit"]] <- fit_dt
		res_list[["positional_loess_fit_points"]] <- dt_fitpts
		res_list[["positional_loess_meta"]] <- data.table::data.table(
		  span = 0.2,
		  ref_condition = ref,
		  fit_consequences = paste(fit_consequences, collapse = ","),
		  used_fallback = use_fallback,
		  fallback_reason = if (use_fallback) "insufficient_or_low_span_strict_points" else NA_character_,
		  t1_centre = t1_centre,
		  t1_min = t1_threshold,
		  step_min = log2_08,
		  c1 = if (is.null(c1)) NA_character_ else c1,
		  c2 = if (is.null(c2)) NA_character_ else c2,
		  c3 = if (is.null(c3)) NA_character_ else c3,
		  c4 = if (is.null(c4)) NA_character_ else c4,
		  c5 = if (is.null(c5)) NA_character_ else c5
		)

		# Apply positional correction onto the adjusted tables
		apply_one <- function(cn) {
		  if (is.null(cn)) return(NULL)
		  if (!cn %in% names(res_list)) stop(paste0("Missing contrast in all_deseq_res_anno_adj: ", cn))
		  apply_positional_correction(res_list[[cn]], fit_dt, c(lfc_raw_col, lfc_shrunk_col))
		}

		                         res_list[[c1]] <- apply_one(c1)
		if (!is.null(c2)) res_list[[c2]] <- apply_one(c2)
		if (!is.null(c3)) res_list[[c3]] <- apply_one(c3)
		if (!is.null(c4)) res_list[[c4]] <- apply_one(c4)
		if (!is.null(c5)) res_list[[c5]] <- apply_one(c5)

		# Compute which contrasts we will process from here on
		pos_contrasts <- Filter(Negate(is.null), c(c1, c2, c3, c4, c5))

		# ---- Recenter positional-adjusted LFCs so controls are at ~0 ----
		for (nm in pos_contrasts) {
		  dt_res <- res_list[[nm]]

		  if ("pos_adj_log2FoldChange_raw" %in% names(dt_res)) {
			dt_res <- recenter_by_controls(
			  dt_res, "pos_adj_log2FoldChange_raw",
			  control_consequences = control_consequences,
			  fallback = "SNV"
			)
		  }

		  if ("pos_adj_log2FoldChange_shrunk" %in% names(dt_res)) {
			dt_res <- recenter_by_controls(
			  dt_res, "pos_adj_log2FoldChange_shrunk",
			  control_consequences = control_consequences,
			  fallback = "SNV"
			)
		  }

		  res_list[[nm]] <- dt_res
		}

		desired_order_all <- c(
		  "sequence", "oligo_name", "position", "consequence",
		  "baseMean_raw", "log2FoldChange_raw", "lfcSE_raw", "pvalue_raw", "padj_raw",
		  "baseMean_shrunk", "log2FoldChange_shrunk", "lfcSE_shrunk", "pvalue_shrunk", "padj_shrunk",

		  "adj_log2FoldChange_raw", "adj_score_raw", "adj_pval_raw", "adj_fdr_raw", "stat_adj_raw",
		  "adj_log2FoldChange_shrunk", "adj_score_shrunk", "adj_pval_shrunk", "adj_fdr_shrunk", "stat_adj_shrunk",

		  "pos_fit", "pos_fit_se", "pos_total_se_raw", "pos_total_se_shrunk",
		  "pos_adj_log2FoldChange_raw", "pos_adj_score_raw", "pos_adj_pval_raw", "pos_adj_fdr_raw", "stat_pos_raw",
		  "pos_adj_log2FoldChange_shrunk", "pos_adj_score_shrunk", "pos_adj_pval_shrunk", "pos_adj_fdr_shrunk", "stat_pos_shrunk"
		)


		# Compute pos-based stats for contrasts 
		pos_contrasts <- Filter(Negate(is.null), c(c1, c2, c3, c4, c5))
		for (nm in pos_contrasts) {
		  dt_res <- res_list[[nm]]

		  # Total SE including LOESS uncertainty
		  dt_res[, pos_total_se_raw := sqrt(lfcSE_raw^2 + pos_fit_se^2)]
		  dt_res[, pos_total_se_shrunk := sqrt(lfcSE_shrunk^2 + pos_fit_se^2)]

		  dt_res[, pos_adj_score_raw := pos_adj_log2FoldChange_raw / pos_total_se_raw]
		  dt_res[, pos_adj_score_shrunk := pos_adj_log2FoldChange_shrunk / pos_total_se_shrunk]

		  dt_res[, pos_adj_pval_raw := pnorm(abs(pos_adj_score_raw), lower.tail = FALSE) * 2]
		  dt_res[, pos_adj_pval_shrunk := pnorm(abs(pos_adj_score_shrunk), lower.tail = FALSE) * 2]

		  dt_res[, pos_adj_fdr_raw := p.adjust(pos_adj_pval_raw, method = "fdr")]
		  dt_res[, pos_adj_fdr_shrunk := p.adjust(pos_adj_pval_shrunk, method = "fdr")]

		  dt_res[, stat_pos_raw := "no impact"]
		  dt_res[(pos_adj_fdr_raw < pcut) & (pos_adj_log2FoldChange_raw > ecut), stat_pos_raw := "enriched"]
		  dt_res[(pos_adj_fdr_raw < pcut) & (pos_adj_log2FoldChange_raw < dcut), stat_pos_raw := "depleted"]
		  dt_res[, stat_pos_raw := factor(stat_pos_raw, levels = c("no impact", "enriched", "depleted"))]

		  dt_res[, stat_pos_shrunk := "no impact"]
		  dt_res[(pos_adj_fdr_raw < pcut) & (pos_adj_log2FoldChange_shrunk > ecut), stat_pos_shrunk := "enriched"]
		  dt_res[(pos_adj_fdr_raw < pcut) & (pos_adj_log2FoldChange_shrunk < dcut), stat_pos_shrunk := "depleted"]
		  dt_res[, stat_pos_shrunk := factor(stat_pos_shrunk, levels = c("no impact", "enriched", "depleted"))]

		  present_all <- intersect(desired_order_all, names(dt_res))
		  setcolorder(dt_res, c(present_all, setdiff(names(dt_res), present_all)))

		  res_list[[nm]] <- dt_res
		}

		object@all_deseq_res_anno_adj <- res_list
		cat("    |--> Positional correction done. Final status columns: stat_pos_raw, stat_pos_shrunk\n")
		
        return(object)
    }
)

#' scale_override for facet_custom
#'
#' @name scale_override
#' @param which  which grid in the facet
#' @param scale  x or y scale
#' @return a structure
scale_override <- function(which, scale) {
    if (!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
        stop("which must be an integer of length 1")
    }

    if (is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
        stop("scale must be an x or y position scale")
    }

    structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
    "CustomFacetWrap", FacetWrap,
    init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
        # make the initial x, y scales list
        scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)

        if(is.null(params$scale_overrides)) return(scales)

        max_scale_x <- length(scales$x)
        max_scale_y <- length(scales$y)

        # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
        for (scale_override in params$scale_overrides) {
            which <- scale_override$which
            scale <- scale_override$scale

            if ("x" %in% scale$aesthetics) {
                if (!is.null(scales$x)) {
                    if (which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
                    scales$x[[which]] <- scale$clone()
                }
            } else if ("y" %in% scale$aesthetics) {
                if (!is.null(scales$y)) {
                    if (which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
                    scales$y[[which]] <- scale$clone()
                }
            } else {
                stop("Invalid scale")
            }
        }

        # return scales
        scales
    }
)

#' facet_wrap_custom
#'
#' @name facet_wrap_custom
#' @param scale_overrides  scale_overrides
facet_wrap_custom <- function(..., scale_overrides = NULL) {
    # take advantage of the sanitizing that happens in facet_wrap
    facet_super <- facet_wrap(...)

    # sanitize scale overrides
    if (inherits(scale_overrides, "scale_override")) {
        scale_overrides <- list(scale_overrides)
    } else if (!is.list(scale_overrides) || !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
        stop("scale_overrides must be a scale_override object or a list of scale_override objects")
    }

    facet_super$params$scale_overrides <- scale_overrides
    ggproto(NULL, CustomFacetWrap, shrink = facet_super$shrink, params = facet_super$params)
}

# Initialize global config variable
maveqc_config <- list(
    gini_coeff = 0.5,
    sqc_total = 1000000,
    sqc_missing = 0.01,
    sqc_low_count = 5,
    sqc_low_sample_per = 0.25,
    sqc_accepted = 1000000,
    sqc_mapping_per = 0.6,
    sqc_ref_per = 0.1,
    sqc_library_per = 0.4,
    sqc_library_cov = 100,
    sqc_low_per = 0.00005,
    sqc_low_lib_per = 0.7,
    expqc_padj = 0.05,
    expqc_lfc_depleted = 0,
    expqc_lfc_enriched = 0,
    expqc_top_variants = 500,
    expqc_lfc_min = -6,
    expqc_lfc_max = 2
)

#' Create the template of config file
#' @param config_dir  the directory of config file
create_config <- function(config_dir) {
    config_path <- file.path(config_dir, "config.yaml")
    
    sink(config_path)
    
    cat("# user defined thresholds for QC\n\n")
    
    cat("#-----------#\n")
    cat("# Sample QC #\n")
    cat("#-----------#\n\n")
    
    cat("# gini coefficient must be lower than 0.5\n")
    cat(sprintf("gini_coeff: %g\n\n", maveqc_config$gini_coeff))
    
    cat("# the sample must have more than 1000000 total reads\n")
    cat(sprintf("sqc_total: %d\n\n", maveqc_config$sqc_total))
    
    cat("# the missing variants in the library must be less than 10%\n")
    cat(sprintf("sqc_missing: %g\n\n", maveqc_config$sqc_missing))
    
    cat("# the sequence must have at least 5 counts in at least 25% of the samples\n")
    cat(sprintf("sqc_low_count: %d\n", maveqc_config$sqc_low_count))
    cat(sprintf("sqc_low_sample_per: %g\n\n", maveqc_config$sqc_low_sample_per))
    
    cat("# the sample must have more than 1000000 reads after the low count filtering\n")
    cat(sprintf("sqc_accepted: %d\n\n", maveqc_config$sqc_accepted))
    
    cat("# the sample must have more than 60% of reads aligned to the library including reference and PAM reads\n")
    cat(sprintf("sqc_mapping_per: %g\n\n", maveqc_config$sqc_mapping_per))
    
    cat("# the sample must have less than 10% of reads aligned to reference sequence\n")
    cat(sprintf("sqc_ref_per: %g\n\n", maveqc_config$sqc_ref_per))
    
    cat("# the sample must have more than 40% of reads aligned to the library\n")
    cat(sprintf("sqc_library_per: %g\n\n", maveqc_config$sqc_library_per))
    
    cat("# the sample must have more than 100x average coverage\n")
    cat("# the number of library reads divided by the number of sequences\n")
    cat(sprintf("sqc_library_cov: %d\n\n", maveqc_config$sqc_library_cov))
    
    cat("# the majority of the variants (>70%) distributed above the 0.005% cutoff for the reference samples\n")
    cat(sprintf("sqc_low_per: %g\n", maveqc_config$sqc_low_per))
    cat(sprintf("sqc_low_lib_per: %g\n\n", maveqc_config$sqc_low_lib_per))
    
    cat("#---------------#\n")
    cat("# Experiment QC #\n")
    cat("#---------------#\n\n")
    
    cat("# DESeq2 relevant cutoffs\n")
    cat(sprintf("expqc_padj: %g\n", maveqc_config$expqc_padj))
    cat(sprintf("expqc_lfc_depleted: %d\n", maveqc_config$expqc_lfc_depleted))
    cat(sprintf("expqc_lfc_enriched: %d\n", maveqc_config$expqc_lfc_enriched))
    cat(sprintf("expqc_top_variants: %d\n\n", maveqc_config$expqc_top_variants))
    
    cat("# Log2 Fold Change cutoffs\n")
    cat(sprintf("expqc_lfc_min: %d\n", maveqc_config$expqc_lfc_min))
    cat(sprintf("expqc_lfc_max: %d\n", maveqc_config$expqc_lfc_max))
    
    sink()
}

#' Load the user's config file
#' @param config_path  the path of config file
load_config <- function(config_path) {
    if (file.exists(config_path)) {
        user_config <- read.config(file = config_path)
        
        # Update global config with user settings while preserving defaults
        config_names <- names(maveqc_config)
        for (name in config_names) {
            if (!is.null(user_config[[name]])) {
                maveqc_config[[name]] <<- user_config[[name]]
            }
        }
        
        # Make config available globally
        assign("maveqc_config", maveqc_config, envir = .GlobalEnv)
    } else {
        stop(paste0("====> Error: ", config_path, " is not found!"))
    }
}

#' initialize function
setGeneric("qcout_samqc_all", function(object, ...) {
  standardGeneric("qcout_samqc_all")
})

#' create all the output files
#'
#' @export
#' @name qcout_samqc_all
#' @param object   sampleQC object
#' @param qc_type  qc type
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_all",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          out_dir = NULL) {
        if (is.null(out_dir)) {
            stop(paste0("====> Error: out_dir is not provided, no output directory."))
        }

        qc_type <- match.arg(qc_type)

        qcout_samqc_cutoffs(object = object, out_dir = out_dir)
        qcout_samqc_readlens(object = object, out_dir = out_dir)
        qcout_samqc_total(object = object, out_dir = out_dir)
        qcout_samqc_missing(object = object, out_dir = out_dir)
        qcout_samqc_accepted(object = object, out_dir = out_dir)
        qcout_samqc_libcov(object = object, out_dir = out_dir)
        qcout_samqc_pos_cov(object = object, qc_type = qc_type, out_dir = out_dir)
        qcout_samqc_results(object = object, qc_type = qc_type, out_dir = out_dir)

        if (qc_type == "screen") {
            qcout_samqc_pos_anno(object = object, out_dir = out_dir)
        }

        qcout_samqc_badseqs(object = object, qc_type = qc_type, out_dir = out_dir)
    }
)

#' initialize function
setGeneric("qcout_samqc_cutoffs", function(object, ...) {
  standardGeneric("qcout_samqc_cutoffs")
})

#' create output file of bad seqs which fail filtering
#'
#' @export
#' @name qcout_samqc_cutoffs
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_cutoffs",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        if (is.null(out_dir)) {
            stop(paste0("====> Error: out_dir is not provided, no output directory."))
        }

        write.table(object@cutoffs,
                    file = paste0(out_dir, "/", "sample_qc_cutoffs.tsv"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)
    }
)

#' initialize function
setGeneric("qcout_samqc_badseqs", function(object, ...) {
  standardGeneric("qcout_samqc_badseqs")
})

#' create output file of bad seqs which fail filtering for screen qc
#'
#' @export
#' @name qcout_samqc_badseqs
#' @param object   sampleQC object
#' @param qc_type  plasmid or screen
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_badseqs",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          out_dir = NULL) {
        if (is.null(out_dir)) {
            stop(paste0("====> Error: out_dir is not provided, no output directory."))
        }

        qc_type <- match.arg(qc_type)

        if (qc_type == "screen") {
            cat("Outputing bad sequences filtered out by clustering...", "\n", sep = "")
            bad_seqs <- data.table()
            for (i in 1:length(object@bad_seqs_bycluster)) {
                bad_tmp <- object@bad_seqs_bycluster[[i]]
                colnames(bad_tmp)[1] <- "seq"
                lib_tmp <- object@samples[[i]]@vep_anno
                res_tmp <- bad_tmp[lib_tmp, on = .(seq), nomatch = 0][, c("unique_oligo_name", "seq", "count")]
                colnames(res_tmp)[3] <- names(object@bad_seqs_bycluster)[i]

                if (nrow(bad_seqs) == 0) {
                    bad_seqs <- res_tmp
                } else {
                    bad_seqs <- merge(bad_seqs, res_tmp, by = c("unique_oligo_name", "seq"), all = TRUE)
                }
            }

            write.table(bad_seqs,
                        file = paste0(out_dir, "/", "failed_variants_by_cluster.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)

            cat("Outputing bad sequences filtered out by sequencing depth...", "\n", sep = "")
            bad_seqs <- data.table()
            for (i in 1:length(object@bad_seqs_bydepth)) {
                bad_tmp <- object@bad_seqs_bydepth[[i]]
                colnames(bad_tmp)[1] <- "seq"
                lib_tmp <- object@samples[[i]]@vep_anno
                res_tmp <- bad_tmp[lib_tmp, on = .(seq), nomatch = 0][, c("unique_oligo_name", "seq", "count")]
                colnames(res_tmp)[3] <- names(object@bad_seqs_bydepth)[i]

                if (nrow(bad_seqs) == 0) {
                    bad_seqs <- res_tmp
                } else {
                    bad_seqs <- merge(bad_seqs, res_tmp, by = c("unique_oligo_name", "seq"), all = TRUE)
                }
            }

            write.table(bad_seqs,
                        file = paste0(out_dir, "/", "failed_variants_by_depth.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)

            cat("Outputing bad sequences filtered out by library mapping...", "\n", sep = "")
            bad_seqs <- data.table()
            for (i in 1:length(object@bad_seqs_bylib)) {
                bad_tmp <- object@bad_seqs_bylib[[i]]
                colnames(bad_tmp)[1] <- "seq"
                lib_tmp <- object@samples[[i]]@vep_anno
                res_tmp <- bad_tmp[lib_tmp, on = .(seq), nomatch = 0][, c("unique_oligo_name", "seq", "count")]
                colnames(res_tmp)[3] <- names(object@bad_seqs_bylib)[i]

                if (nrow(bad_seqs) == 0) {
                    bad_seqs <- res_tmp
                } else {
                    bad_seqs <- merge(bad_seqs, res_tmp, by = c("unique_oligo_name", "seq"), all = TRUE)
                }
            }

            write.table(bad_seqs,
                        file = paste0(out_dir, "/", "failed_variants_by_mapping.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }

        cat("Outputing missing variants in the library...", "\n", sep = "")
        bad_seqs <- data.table()
        for (i in 1:length(object@samples)) {
            bad_tmp <- object@samples[[i]]@libcounts[count == 0]

            if (nrow(bad_seqs) == 0) {
                bad_seqs <- bad_tmp
            } else {
                bad_seqs <- merge(bad_seqs, bad_tmp, all = TRUE)
            }
        }

        write.table(bad_seqs,
                    file = paste0(out_dir, "/", "missing_variants_in_library.tsv"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)
    }
)

#' initialize function
setGeneric("qcout_samqc_readlens", function(object, ...) {
  standardGeneric("qcout_samqc_readlens")
})

#' create output file of total reads stats
#'
#' @export
#' @name qcout_samqc_readlens
#' @param object    sampleQC object
#' @param len_bins  the bins of length distribution
#' @param out_dir   the output directory
setMethod(
    "qcout_samqc_readlens",
    signature = "sampleQC",
    definition = function(object,
                          len_bins = seq(0, 300, 50),
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Total Reads",
                  "% 0 ~ 50",
                  "% 50 ~ 100",
                  "% 100 ~ 150",
                  "% 150 ~ 200",
                  "% 200 ~ 250",
                  "% 250 ~ 300",
                  "Pass Threshold (%)",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@stats$total_reads

        bin_per <- data.frame()
        for (i in 1:length(object@lengths)) {
            tmp_lens <- object@lengths[[i]]$length
            h <- hist(tmp_lens, breaks = len_bins, plot = FALSE)
            bin_per <- rbind(bin_per, round(h$counts / nrow(object@samples[[i]]@allcounts) * 100, 1))
        }

        df_outs[, 4] <- bin_per[, 1]
        df_outs[, 5] <- bin_per[, 2]
        df_outs[, 6] <- bin_per[, 3]
        df_outs[, 7] <- bin_per[, 4]
        df_outs[, 8] <- bin_per[, 5]
        df_outs[, 9] <- bin_per[, 6]
        df_outs[, 10] <- 90
        df_outs[, 11] <- (df_outs[, 8] + df_outs[, 9]) > df_outs[, 10]

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Total Reads" = colDef(format = colFormat(separators = TRUE)),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_read_length.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_missing", function(object, ...) {
  standardGeneric("qcout_samqc_missing")
})

#' create output file of bad seqs which fail filtering
#'
#' @export
#' @name qcout_samqc_missing
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_missing",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Library Sequences",
                  "Missing Library Sequences",
                  "% Missing Library Sequences",
                  "Pass Threshold (%)",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@stats$library_seqs
        df_outs[, 4] <- object@stats$missing_meta_seqs
        tmp_out <- object@stats$per_missing_meta_seqs * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 5] <- tmp_out
        tmp_out <- object@cutoffs$per_missing_variants * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 6] <- tmp_out
        df_outs[, 7] <- object@stats$qcpass_missing_per

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Library Sequences" = colDef(format = colFormat(separators = TRUE)),
                          "Missing Library Sequences" = colDef(format = colFormat(separators = TRUE)),
                          "% Missing Library Sequences" = colDef(format = colFormat(separators = TRUE),
                                                                 style = function(value) {
                                                                             if (value > object@cutoffs$per_missing_variants * 100) {
                                                                                 color <- "red"
                                                                                 fweight <- "bold"
                                                                             } else {
                                                                                 color <- "forestgreen"
                                                                                 fweight <- "plain"
                                                                             }
                                                                             list(color = color, fontWeight = fweight)}),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_stats_missing.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_total", function(object, ...) {
  standardGeneric("qcout_samqc_total")
})

#' create output file of total reads stats
#'
#' @export
#' @name qcout_samqc_total
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_total",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Accepted Reads",
                  "% Accepted Reads",
                  "Excluded Reads",
                  "% Excluded Reads",
                  "Total Reads",
                  "Pass Threshold",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@stats$accepted_reads
        tmp_out <- object@stats$accepted_reads / object@stats$total_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 4] <- tmp_out
        df_outs[, 5] <- object@stats$excluded_reads
        tmp_out <- object@stats$excluded_reads / object@stats$total_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 6] <- tmp_out
        df_outs[, 7] <- object@stats$total_reads
        df_outs[, 8] <- object@cutoffs$num_total_reads
        df_outs[, 9] <- object@stats$qcpass_accepted_reads

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Accepted Reads" = colDef(format = colFormat(separators = TRUE)),
                          "Excluded Reads" = colDef(format = colFormat(separators = TRUE)),
                          "Total Reads" = colDef(format = colFormat(separators = TRUE),
                                                 style = function(value) {
                                                             if (value < object@cutoffs$num_total_reads) {
                                                                 color <- "red"
                                                                 fweight <- "bold"
                                                             } else {
                                                                 color <- "forestgreen"
                                                                 fweight <- "plain"
                                                             }
                                                             list(color = color, fontWeight = fweight)}),
                          "Pass Threshold" = colDef(format = colFormat(separators = TRUE)),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_stats_total.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_accepted", function(object, ...) {
  standardGeneric("qcout_samqc_accepted")
})

#' create output file of library reads stats
#'
#' @export
#' @name qcout_samqc_accepted
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_accepted",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "% Library Reads",
                  "% Reference Reads",
                  "% PAM Reads",
                  "% Unmapped Reads",
                  "Pass Threshold (%)",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        tmp_out <- object@stats$per_library_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 3] <- tmp_out
        tmp_out <- object@stats$per_ref_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 4] <- tmp_out
        tmp_out <- object@stats$per_pam_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 5] <- tmp_out
        tmp_out <- object@stats$per_unmapped_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 6] <- tmp_out
        df_outs[, 7] <- object@cutoffs$per_library_reads * 100
        df_outs[, 8] <- object@stats$qcpass_library_per

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "% Library Reads" = colDef(style = function(value) {
                                                                 if (value < object@cutoffs$per_library_reads * 100) {
                                                                    color <- "red"
                                                                    fweight <- "bold"
                                                                } else {
                                                                    color <- "forestgreen"
                                                                    fweight <- "plain"
                                                                }
                                                                list(color = color, fontWeight = fweight)}),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_stats_accepted.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_libcov", function(object, ...) {
  standardGeneric("qcout_samqc_libcov")
})

#' create output file of library coverage
#'
#' @export
#' @name qcout_samqc_libcov
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_libcov",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Total Library Reads",
                  "Total Library Sequences",
                  "Library Coverage",
                  "Median Coverage",
                  "Pass Threshold",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@stats$library_reads
        df_outs[, 4] <- object@stats$library_seqs
        tmp_out <- object@stats$library_cov
        tmp_out <- sapply(tmp_out, function(x) round(x, 0))
        df_outs[, 5] <- tmp_out
        df_outs[, 6] <- object@stats$median_cov
        df_outs[, 7] <- object@cutoffs$library_cov
        df_outs[, 8] <- object@stats$qcpass_library_cov

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Total Library Reads" = colDef(format = colFormat(separators = TRUE)),
                          "Total Library Sequences" = colDef(format = colFormat(separators = TRUE)),
                          "Library Coverage" = colDef(format = colFormat(separators = TRUE),
                                                      style = function(value) {
                                                                  if (value < object@cutoffs$library_cov) {
                                                                      color <- "red"
                                                                      fweight <- "bold"
                                                                  } else {
                                                                      color <- "forestgreen"
                                                                      fweight <- "plain"
                                                                  }
                                                                  list(color = color, fontWeight = fweight)}),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_stats_coverage.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_pos_cov", function(object, ...) {
  standardGeneric("qcout_samqc_pos_cov")
})

#' create output file of lof percentages
#'
#' @export
#' @name qcout_samqc_pos_cov
#' @param object   sampleQC object
#' @param qc_type  screen or plasmid
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_pos_cov",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          out_dir = NULL) {
        qc_type <- match.arg(qc_type)

        library_dependent_counts <- data.table()
        if (qc_type == "screen") {
            library_dependent_counts <- merge_list_to_dt(object@library_counts, "sequence", "count")
            library_dependent_counts[, sequence := NULL]
        } else {
            list_counts <- list()
            for (i in 1:length(object@samples)) {
                sample <- object@samples[[i]]@sample
                list_counts[[sample]] <- object@samples[[i]]@libcounts[, c("id", "count")]
            }

            library_dependent_counts <- merge_list_to_dt(list_counts, "id", "count")
            library_dependent_counts[, id := NULL]
        }

        write.table(library_dependent_counts,
                    file = paste0(out_dir, "/", "sample_qc_stats_pos_counts.tsv"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)

        cols <- c("Group",
                  "Sample",
                  "Chromosome",
                  "Strand",
                  "Genomic Start",
                  "Genomic End",
                  "% Low Abundance",
                  "Low Abundance cutoff",
                  "Pass Threshold (%)",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- sapply(object@library_counts_chr, function (x) x[[1]])
        df_outs[, 4] <- sapply(object@library_counts_chr, function (x) x[[2]])
        df_outs[, 5] <- sapply(object@library_counts_chr, function (x) x[[3]])
        df_outs[, 6] <- sapply(object@library_counts_chr, function (x) x[[4]])

        low_per <- vector()
        if (qc_type == "screen") {
            for (s in object@samples) {
                tmp_num <- sum(object@library_counts_pos[[s@sample]]$count < object@cutoffs$seq_low_count)
                low_per <- append(low_per, round(tmp_num / nrow(object@library_counts_pos[[s@sample]]) * 100, 2))
            }
        } else {
            for (s in object@samples) {
                tmp_num <- sum(s@libcounts$count < object@cutoffs$seq_low_count)
                low_per <- append(low_per, round(tmp_num / nrow(object@library_counts_pos[[s@sample]]) * 100, 2))
            }
        }
        df_outs[, 7] <- low_per

        df_outs[, 8] <- object@cutoffs$seq_low_count
        df_outs[, 9] <- (1 - object@cutoffs$low_abundance_lib_per) * 100
        df_outs[, 10] <- df_outs[, 7] < (1 - object@cutoffs$low_abundance_lib_per) * 100

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Genomic Start" = colDef(format = colFormat(separators = TRUE)),
                          "Genomic End" = colDef(format = colFormat(separators = TRUE)),
                          "% Low Abundance" = colDef(minWidth = 200,
                                                           style = function(value) {
                                                                       if (value > (1 - object@cutoffs$low_abundance_lib_per) * 100) {
                                                                           color <- "red"
                                                                           fweight <- "bold"
                                                                       } else {
                                                                           color <- "forestgreen"
                                                                           fweight <- "plain"
                                                                       }
                                                                       list(color = color, fontWeight = fweight)}),
                          "Low Abundance cutoff" = colDef(minWidth = 200),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_stats_pos_coverage.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_pos_anno", function(object, ...) {
  standardGeneric("qcout_samqc_pos_anno")
})

#' create output file of lof percentages
#'
#' @export
#' @name qcout_samqc_pos_anno
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_pos_anno",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Chromosome",
                  "Strand",
                  "Genomic Start",
                  "Genomic End",
                  "% Low Abundance (LOF)",
                  "% Low Abundance (Others)",
                  "% Low Abundance (ALL)",
                  "% Low Abundance cutoff",
                  "Pass Threshold",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- sapply(object@library_counts_chr, function (x) x[[1]])
        df_outs[, 4] <- sapply(object@library_counts_chr, function (x) x[[2]])
        df_outs[, 5] <- sapply(object@library_counts_chr, function (x) x[[3]])
        df_outs[, 6] <- sapply(object@library_counts_chr, function (x) x[[4]])

        libcounts_pos <- as.data.frame(object@library_counts_pos_anno)
        libcounts_pos <- libcounts_pos[, c(rownames(object@stats), "consequence")]
        libcounts_pos$consequence <- ifelse(libcounts_pos$consequence == "LOF", "LOF", "Others")
        libcounts_pos[, rownames(object@stats)] <- t(t(libcounts_pos[, rownames(object@stats)]) / object@stats$accepted_reads * 100)

        # what about NA?
        #libcounts_pos[is.na(libcounts_pos)] <- 0

        lof_counts <- libcounts_pos[libcounts_pos$consequence == "LOF", rownames(object@stats)]
        # the number of seqs with low abundance
        lof_low_num <- colSums(lof_counts < object@cutoffs$low_abundance_per * 100, na.rm = TRUE)
        # the percentage of seqs with low abundance
        lof_low_per <- lof_low_num / nrow(libcounts_pos) * 100
        lof_low_per <- round(lof_low_per, 1)

        others_counts <- libcounts_pos[libcounts_pos$consequence == "Others", rownames(object@stats)]
        # the number of seqs with low abundance
        others_low_num <- colSums(others_counts < object@cutoffs$low_abundance_per * 100, na.rm = TRUE)
        # the percentage of seqs with low abundance
        others_low_per <- others_low_num / nrow(libcounts_pos) * 100
        others_low_per <- round(others_low_per, 1)

        df_outs[, 7] <- lof_low_per
        df_outs[, 8] <- others_low_per
        df_outs[, 9] <- lof_low_per + others_low_per
        df_outs[, 10] <- object@cutoffs$low_abundance_per * 100
        df_outs[, 11] <- (1 - object@cutoffs$low_abundance_lib_per) * 100
        df_outs[, 12] <- df_outs[, 9] < (1 - object@cutoffs$low_abundance_lib_per) * 100

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Genomic Start" = colDef(format = colFormat(separators = TRUE)),
                          "Genomic End" = colDef(format = colFormat(separators = TRUE)),
                          "% Low Abundance (LOF)" = colDef(minWidth = 200),
                          "% Low Abundance (Others)" = colDef(minWidth = 200),
                          "% Low Abundance (ALL)" = colDef(minWidth = 200,
                                                           style = function(value) {
                                                                       if (value > (1 - object@cutoffs$low_abundance_lib_per) * 100) {
                                                                           color <- "red"
                                                                           fweight <- "bold"
                                                                       } else {
                                                                           color <- "forestgreen"
                                                                           fweight <- "plain"
                                                                       }
                                                                       list(color = color, fontWeight = fweight)}),
                          "% Low Abundance cutoff" = colDef(minWidth = 200),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_stats_pos_percentage.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_results", function(object, ...) {
  standardGeneric("qcout_samqc_results")
})

#' create all the output files
#'
#' @export
#' @name qcout_samqc_results
#' @param object   sampleQC object
#' @param qc_type  qc type
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_results",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          out_dir = NULL) {
        qc_type <- match.arg(qc_type)

        cols <- c("Group",
                  "Sample",
                  "Gini Coefficient",
                  "Gini Coefficient Pass",
                  "Total Reads",
                  "% Missing Variants",
                  "Accepted Reads",
                  "% Mapping Reads",
                  "% Reference Reads",
                  "% Library Reads",
                  "Library Coverage",
                  "% R1 Adaptor",
                  "% R2 Adaptor")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@stats$gini_coeff_before_qc
        df_outs[, 4] <- ifelse(df_outs[, 3] < maveqc_config$gini_coeff, TRUE, FALSE)
        df_outs[, 5] <- object@stats$qcpass_total_reads
        df_outs[, 6] <- object@stats$qcpass_missing_per
        df_outs[, 7] <- object@stats$qcpass_accepted_reads
        df_outs[, 8] <- object@stats$qcpass_mapping_per
        df_outs[, 9] <- object@stats$qcpass_ref_per
        df_outs[, 10] <- object@stats$qcpass_library_per
        df_outs[, 11] <- object@stats$qcpass_library_cov
        df_outs[, 12] <- object@stats$per_r1_adaptor
        df_outs[, 13] <- object@stats$per_r2_adaptor

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (qc_type == "screen") {
            df_outs <- df_outs
        }

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Total Reads" = colDef(cell = function(value) { if (value) "\u2705" else "\u274c" }),
                          "% Missing Variants" = colDef(cell = function(value) { if (value) "\u2705" else "\u274c" }),
                          "Accepted Reads" = colDef(cell = function(value) { if (value) "\u2705" else "\u274c" }),
                          "% Mapping Reads" = colDef(cell = function(value) { if (value) "\u2705" else "\u274c" }),
                          "% Reference Reads" = colDef(cell = function(value) { if (value) "\u2705" else "\u274c" }),
                          "% Library Reads" = colDef(cell = function(value) { if (value) "\u2705" else "\u274c" }),
                          "Library Coverage" = colDef(cell = function(value) { if (value) "\u2705" else "\u274c" }))
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_results.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)


#####################################################################################################################################################

#' initialize function
setGeneric("qcout_expqc_all", function(object, ...) {
  standardGeneric("qcout_expqc_all")
})

#' create all the output files
#'
#' @export
#' @name qcout_expqc_all
#' @param object   experimentQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_expqc_all",
    signature = "experimentQC",
    definition = function(object,
                          out_dir = NULL) {
        if (is.null(out_dir)) {
            stop(paste0("====> Error: out_dir is not provided, no output directory."))
        }

        qcout_expqc_corr(object = object, out_dir = out_dir)
        qcout_expqc_deseq(object = object, eqc_type = "all", out_dir = out_dir)
    }
)

#' initialize function
setGeneric("qcout_expqc_corr", function(object, ...) {
  standardGeneric("qcout_expqc_corr")
})

#' create output file of clustering and correlation results
#'
#' @export
#' @name qcout_expqc_corr
#' @param object   experimentQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_expqc_corr",
    signature = "experimentQC",
    definition = function(object,
                          out_dir = NULL) {
        df_outs <- as.data.frame(object@lib_corr_res)

        num_clusters <- length(unique(object@coldata$condition))
        name_clusters <- as.vector(unique(object@coldata$condition))
        sample_clusters <- cutree(object@lib_hclust_res, num_clusters)

        df_outs <- cbind(object@coldata[rownames(df_outs), ], df_outs)
        colnames(df_outs)[1] <- "Replicate"
        colnames(df_outs)[2] <- "Condition"
        df_outs <- cbind(sample_clusters[rownames(df_outs)], df_outs)
        colnames(df_outs)[1] <- "Cluster"
        df_outs <- cbind(rownames(df_outs), df_outs)
        colnames(df_outs)[1] <- "Sample"

        df_outs$Pass <- NA
        for (i in 1:length(name_clusters)) {
            tmp_clusters <- df_outs[df_outs$Condition == name_clusters[i], ]$Cluster
            if (length(unique(tmp_clusters)) == 1) {
                df_outs[df_outs$Condition == name_clusters[i], ]$Pass <- TRUE
            } else {
                df_outs[df_outs$Condition == name_clusters[i], ]$Pass <- FALSE
            }
        }

        df_outs_tmp <- unique(df_outs[df_outs$Pass, c("Condition", "Cluster")])
        for (i in 1:nrow(df_outs_tmp)) {
            pass_check <- TRUE
            for (j in 1:nrow(df_outs_tmp)) {
                if (i != j) {
                    if (df_outs_tmp[i, ]$Cluster == df_outs_tmp[j, ]$Cluster) {
                        pass_check <- FALSE
                    }
                }
            }
            df_outs[df_outs$Condition == df_outs_tmp[i, ]$Condition, ]$Pass <- pass_check
        }

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        is.num <- sapply(df_outs, is.numeric)
        df_outs[is.num] <- lapply(df_outs[is.num], round, 2)

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Sample" = colDef(minWidth = 100),
                          "Cluster" = colDef(minWidth = 100),
                          "Replicate" = colDef(minWidth = 100),
                          "Condition" = colDef(minWidth = 100),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "experiment_qc_corr.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_expqc_deseq", function(object, ...) {
  standardGeneric("qcout_expqc_deseq")
})

#' create all the output files
#'
#' @export
#' @name qcout_expqc_deseq
#' @param object   experimentQC object
#' @param eqc_type library counts or all counts
#' @param out_dir  the output directory
setMethod(
    "qcout_expqc_deseq",
    signature = "experimentQC",
    definition = function(object,
                          eqc_type = c("lib", "all"),
                          out_dir  = NULL,
                          stat_col = NULL,
                          lfc_col  = NULL) {
        if (is.null(out_dir)) {
            stop(paste0("====> Error: out_dir is not provided, no output directory."))
        }

        eqc_type <- match.arg(eqc_type)

        # Current schema splits everything into _raw/_shrunk (and, for "all",
        # further into adj_*/pos_adj_* variants) -- there is no longer a bare
        # "stat" or "adj_log2FoldChange" column on either table, so pick the
        # matching current column names (same defaults as qcplot_expqc_deseq_fc).
        if (eqc_type == "lib") {
            if (is.null(stat_col)) stat_col <- "stat_raw"
            if (is.null(lfc_col))  lfc_col  <- "log2FoldChange_raw"
        } else {
            if (is.null(stat_col)) stat_col <- "stat_pos_raw"
            if (is.null(lfc_col))  lfc_col  <- "pos_adj_log2FoldChange_raw"
        }

        for (i in 1:length(object@comparisons)) {
            if (eqc_type == "lib") {
                df_outs <- object@lib_deseq_res_anno[[i]]
            } else {
                df_outs <- object@all_deseq_res_anno_adj[[i]]
            }

            stopifnot(stat_col %in% colnames(df_outs))
            stopifnot(lfc_col  %in% colnames(df_outs))

            tmpcons <- sort(unique(df_outs$consequence))
            cols <- c("consequence",
                      "number of depleted",
                      "number of no impact",
                      "number of enriched",
                      "lfc range",
                      "total number",
                      "total number of outside range")
            df_outs_summary <- data.frame(matrix(NA, length(tmpcons), length(cols)))
            colnames(df_outs_summary) <- cols

            df_outs_summary[, 1] <- tmpcons
            df_outs_summary[, 5] <- paste0(maveqc_config$expqc_lfc_min, " ~ ", maveqc_config$expqc_lfc_max)
            for (j in 1:nrow(df_outs_summary)) {
                df_tmp <- df_outs[df_outs$consequence == df_outs_summary$consequence[j], ]

                df_outs_summary[j, 2] <- sum(df_tmp[[stat_col]] == "depleted", na.rm = TRUE)
                df_outs_summary[j, 3] <- sum(df_tmp[[stat_col]] == "no impact", na.rm = TRUE)
                df_outs_summary[j, 4] <- sum(df_tmp[[stat_col]] == "enriched", na.rm = TRUE)

                df_outs_summary[j, 6] <- nrow(df_tmp)
                df_outs_summary[j, 7] <- sum(df_tmp[[lfc_col]] < maveqc_config$expqc_lfc_min | df_tmp[[lfc_col]] > maveqc_config$expqc_lfc_max, na.rm = TRUE)
            }

            dt_outs_summary <- as.data.table(df_outs_summary)
            dt1 <- setorder(dt_outs_summary[consequence %in% c("Synonymous_Variant", "LOF")], -consequence)
            dt2 <- dt_outs_summary[consequence %nin% c("Synonymous_Variant", "LOF")]
            dt_outs_summary <- rbind(dt1, dt2)

            write.table(df_outs,
                        file = paste0(out_dir, "/", "experiment_qc_deseq_fc.", object@comparisons[[i]], ".", eqc_type, ".tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)

            write.table(dt_outs_summary,
                        file = paste0(out_dir, "/", "experiment_qc_deseq_fc.", object@comparisons[[i]], ".", eqc_type, "_sum.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcplot_samqc_all", function(object, ...) {
  standardGeneric("qcplot_samqc_all")
})

#' create all the plot figures
#'
#' @export
#' @name qcplot_samqc_all
#' @param object   sampleQC object
#' @param qc_type  qc type
#' @param samples  samples for LOF annotation plot
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_samqc_all",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          samples = maveqc_ref_time_point_samples,
                          plot_dir = NULL) {
        if (is.null(plot_dir)) {
            stop(paste0("====> Error: plot_dir is not provided, no output directory."))
        }

        qc_type <- match.arg(qc_type)

        if (qc_type == "plasmid") {
            qcplot_samqc_readlens(object = object, plot_dir = plot_dir)
            qcplot_samqc_total(object = object, plot_dir = plot_dir)
            qcplot_samqc_accepted(object = object, plot_dir = plot_dir)
            qcplot_samqc_pos_cov(object = object, qc_type = qc_type, plot_dir = plot_dir)
        } else {
            if (is.null(samples)) {
                stop(paste0("====> Error: please provide samples, a vector."))
            }
            qcplot_samqc_readlens(object = object, plot_dir = plot_dir)
            qcplot_samqc_total(object = object, plot_dir = plot_dir)
            qcplot_samqc_accepted(object = object, plot_dir = plot_dir)
            qcplot_samqc_pos_cov(object = object, qc_type = qc_type, plot_dir = plot_dir)
            qcplot_samqc_pos_anno(object = object, samples = samples, plot_dir = plot_dir)
        }
    }
)

#' initialize function
setGeneric("qcplot_samqc_readlens", function(object, ...) {
  standardGeneric("qcplot_samqc_readlens")
})

#' create the read length plot
#'
#' @export
#' @name qcplot_samqc_readlens
#' @param object   sampleQC object
#' @param len_bins the bins of length distribution
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_samqc_readlens",
    signature = "sampleQC",
    definition = function(object,
                          len_bins = seq(0, 300, 50),
                          plot_dir = NULL) {
        read_lens <- data.table()
        for (i in 1:length(object@lengths)) {
            tmp_lens <- object@lengths[[i]][, "length", drop = FALSE]
            tmp_lens$sample <- names(object@lengths)[i]
            tmp_lens <- as.data.table(tmp_lens)

            if (nrow(read_lens) == 0) {
                read_lens <- tmp_lens
            } else {
                read_lens <- rbind(read_lens, tmp_lens)
            }
        }

        sample_names <- vector()
        for (s in object@samples) {
            sample_names <- append(sample_names, s@sample)
        }

        read_lens <- as.data.frame(read_lens)
        read_lens$sample <- factor(read_lens$sample, levels = sample_names)

        p1 <- ggplot(read_lens, aes(x = length)) +
                geom_histogram(aes(y = after_stat(width * density)), breaks = len_bins, color = "black", fill = "grey") +
                geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", color = "yellowgreen", linewidth = 0.4) +
                scale_y_continuous(labels = scales::percent) +
                coord_transform(y = "sqrt") +
                labs(x = "Length Distribution", y = "Composition Percentage", title = "Sample QC read lengths") +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 12, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 12, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 8, face = "bold")) +
                facet_wrap(~sample, scales = "free", dir = "h", ncol = 3)

        p2 <- ggplot(read_lens, aes(x = length)) +
                geom_histogram(aes(y = after_stat(width * density)), breaks = len_bins, color = "black", fill = "grey") +
                geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", color = "yellowgreen", linewidth = 0.4) +
                scale_y_continuous(labels = scales::percent) +
                labs(x = "Length Distribution", y = "Composition Percentage", title = "Sample QC read lengths") +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 12, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 12, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 8, face = "bold")) +
                facet_wrap(~sample, scales = "free", dir = "h", ncol = 3)

        pheight <- 400 * ceiling((length(sample_names) / 3))

        if (is.null(plot_dir)) {
            ggplotly(p2)
        } else {
            png(paste0(plot_dir, "/", "sample_qc_read_length.png"), width = 1200, height = pheight, res = 200)
            print(p1)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_samqc_clusters", function(object, ...) {
  standardGeneric("qcplot_samqc_clusters")
})

#' create the sequence counts and clusters plot
#'
#' @export
#' @name qcplot_samqc_clusters
#' @param object    sampleQC object
#' @param qc_type   qc type for plot
#' @param plot_dir  the output plot directory
setMethod(
    "qcplot_samqc_clusters",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          plot_dir = NULL) {
        qc_type <- match.arg(qc_type)

        if (qc_type == "screen") {
            seq_clusters <- object@seq_clusters[[1]]

            seq_breaks <- seq(0, round(max(seq_clusters$count_log2)), 2)
            select_colors <- select_colorblind("col8")[1:2]
            fill_colors <- sapply(select_colors, function(x) t_col(x, 0.5), USE.NAMES = FALSE)

            p1 <- ggplot(seq_clusters, aes(x = 1:dim(seq_clusters_new)[1], y = count_log2, color = factor(cluster))) +
                    geom_point(shape = 21, size = 0.3, aes(fill = factor(cluster), color = factor(cluster))) +
                    coord_polar() +
                    scale_fill_manual(values = fill_colors) +
                    scale_color_manual(values = select_colors) +
                    labs(x = "sequence index", y = "log2(count+1)", title = "Sample QC clusters") +
                    annotate("text", x = 0, y = seq_breaks, label = seq_breaks, size = 3) +
                    scale_y_continuous(breaks = seq_breaks) +
                    theme(panel.grid.major.x = element_blank()) +
                    theme(panel.grid.major.y = element_line(color = "darkgrey", linewidth = 0.1)) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title.x = element_blank()) +
                    theme(axis.title.y = element_text(size = 12, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 12, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_blank(), axis.ticks = element_blank())

            p2 <- ggplot(seq_clusters_new, aes(x = count_log2, color = factor(cluster))) +
                    geom_density(aes(fill = factor(cluster), color = factor(cluster))) +
                    scale_fill_manual(values = c(t_col("tomato", 0.5), t_col("royalblue", 0.5))) +
                    scale_color_manual(values = c("tomato", "royalblue")) +
                    labs(x = "log2(count+1)", y = "frequency", title = "Sample QC clusters") +
                    theme(legend.position = "none", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 12, face = "bold"))
        } else {
            seq_clusters <- data.table()
            for (i in 1:length(object@seq_clusters)) {
                tmp_cluster <- object@seq_clusters[[i]][, c("count_log2", "cluster")]
                tmp_cluster$samples <- names(object@seq_clusters)[i]
                tmp_cluster[cluster == 1, group := "low-count cluster"]
                tmp_cluster[cluster == 2, group := "high-count cluster"]

                if (nrow(seq_clusters) == 0) {
                    seq_clusters <- tmp_cluster
                } else {
                    seq_clusters <- rbind(seq_clusters, tmp_cluster)
                }
            }

            p1 <- ggplot(seq_clusters, aes(x = count_log2, color = samples)) +
                    geom_density() +
                    labs(x = "log2(count+1)", y = "frequency", title = "Sample QC clusters") +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 12, face = "bold")) +
                    facet_wrap(~group, scales = "free", dir = "h")

            p2 <- p1
        }

        if (is.null(plot_dir)) {
            ggplotly(p2)
        } else {
            png(paste0(plot_dir, "/", "sample_qc_seq_clusters.png"), width = 1200, height = 1200, res = 200)
            print(p1)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_samqc_total", function(object, ...) {
  standardGeneric("qcplot_samqc_total")
})

#' create the stats plot
#'
#' @export
#' @name qcplot_samqc_total
#' @param object   sampleQC object
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_samqc_total",
    signature = "sampleQC",
    definition = function(object,
                          plot_dir = NULL) {
        df_total <- object@stats[, c("excluded_reads", "accepted_reads")]
        df_total$samples <- rownames(df_total)
        dt_total <- reshape2::melt(as.data.table(df_total), id.vars = "samples", variable.name = "types", value.name = "counts")

        dt_total$samples <- factor(dt_total$samples, levels = mixedsort(levels(factor(dt_total$samples))))

        select_colors <- select_colorblind("col8")[1:2]
        fill_colors <- sapply(select_colors, function(x) t_col(x, 0.5), USE.NAMES = FALSE)

        p1 <- ggplot(dt_total,  aes(x = samples, y = counts, fill = types)) +
                geom_bar(stat = "identity") +
                scale_fill_manual(values = fill_colors) +
                scale_color_manual(values = select_colors) +
                labs(x = "samples", y = "counts", title = "Sample QC Stats") +
                scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
                theme(legend.position = "right", legend.title = element_blank()) +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 8, face = "bold")) +
                theme(axis.text.x = element_text(angle = 90))

        pwidth <- 150 * nrow(df_total)

        if (is.null(plot_dir)) {
            ggplotly(p1)
        } else {
            png(paste0(plot_dir, "/", "sample_qc_stats_total.png"), width = pwidth, height = 1200, res = 200)
            print(p1)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_samqc_accepted", function(object, ...) {
  standardGeneric("qcplot_samqc_accepted")
})

#' create the stats plot
#'
#' @export
#' @name qcplot_samqc_accepted
#' @param object   sampleQC object
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_samqc_accepted",
    signature = "sampleQC",
    definition = function(object,
                          plot_dir = NULL) {
        df_accepted <- object@stats[, c("per_unmapped_reads", "per_ref_reads", "per_pam_reads", "per_library_reads")]
        colnames(df_accepted) <- c("unmapped_reads", "ref_reads", "pam_reads", "library_reads")
        df_accepted <- round(df_accepted * 100, 1)
        df_accepted$samples <- rownames(df_accepted)
        dt_filtered <- reshape2::melt(as.data.table(df_accepted), id.vars = "samples", variable.name = "types", value.name = "percent")

        dt_filtered$samples <- factor(dt_filtered$samples, levels = mixedsort(levels(factor(dt_filtered$samples))))

        df_cov <- object@stats[, c("total_reads", "library_reads", "library_cov")]
        colnames(df_cov) <- c("num_total_reads", "num_library_reads", "library_cov")
        df_cov$samples <- rownames(df_cov)
        df_cov$type <- "coverage"

        select_colors <- select_colorblind("col8")[1:4]
        fill_colors <- sapply(select_colors, function(x) t_col(x, 0.5), USE.NAMES = FALSE)

        y_scale <- max(df_cov$library_cov) * 2

        p1 <- ggplot(dt_filtered,  aes(x = samples, y = percent, fill = types)) +
                geom_bar(stat = "identity", position = "fill") +
                geom_line(data = df_cov, aes(x = samples, y = library_cov / y_scale, group = 1), linetype = "dashed", color = "red", inherit.aes = FALSE) +
                geom_point(data = df_cov, aes(x = samples, y = library_cov / y_scale, color = type), shape = 18, size = 3, inherit.aes = FALSE) +
                scale_y_continuous(labels = scales::percent, sec.axis = sec_axis(~. * y_scale, name = "library coverage")) +
                scale_fill_manual(values = fill_colors) +
                scale_color_manual(values = "red") +
                labs(x = "samples", y = "percent", title = "Sample QC Stats") +
                theme(legend.position = "right", legend.title = element_blank()) +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 8, face = "bold")) +
                theme(axis.text.x = element_text(angle = 90)) +
                geom_text(aes(label = percent), position = position_fill(vjust = 0.5), size = 3)

        pwidth <- 150 * nrow(df_accepted)

        if (is.null(plot_dir)) {
            dt_filtered$types <- factor(dt_filtered$types, levels = rev(levels(dt_filtered$types)))

            ay <- list(overlaying = "y",
                       side = "right",
                       title = "Library Coverage")

            mk <- list(size = 12,
                       symbol = "diamond",
                       color = "red")

            plot_ly(data = dt_filtered, x = ~samples, y = ~percent, color = ~types, type = "bar", colors = rev(fill_colors)) %>%
                layout(barmode = "stack") %>%
                add_markers(data = df_cov, x = ~samples, y = ~library_cov, inherit = FALSE, yaxis = "y2", marker = mk, name = "library") %>%
                layout(yaxis2 = ay)
        } else {
            png(paste0(plot_dir, "/", "sample_qc_stats_accepted.png"), width = pwidth, height = 1200, res = 200)
            print(p1)
            dev.off()
        }

        # bubble plot, may be useful, leave it here

        # p2 <- ggplot(df_cov,  aes(x = total_reads, y = library_reads, color = samples)) +
        #         geom_point(alpha = 0.7, aes(size = library_cov)) +
        #         geom_text(size = 2, color = "black", aes(label = library_cov)) +
        #         geom_text(size = 2, color = "black", vjust = -1, aes(label = samples)) +
        #         labs(x = "total reads", y = "library reads", title = "Sample QC Stats") +
        #         scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
        #         scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
        #         scale_size_continuous(range = c(6, 12)) +
        #         theme(legend.position = "right") +
        #         theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
        #         theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
        #         theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
        #         theme(axis.text = element_text(size = 12, face = "bold"))

        # png(paste0(plot_dir, "/", "sample_qc_stats_cov.png"), width = 1200, height = 1200, res = 200)
        # print(p2)
        # dev.off()
    }
)

#' initialize function
setGeneric("qcplot_samqc_gini", function(object, ...) {
  standardGeneric("qcplot_samqc_gini")
})

#' create the gini plot
#'
#' @export
#' @name qcplot_samqc_gini
#' @param object   sampleQC object
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_samqc_gini",
    signature = "sampleQC",
    definition = function(object,
                          plot_dir = NULL) {
        sample_names <- character()
        all_gini <- character()
        for (s in object@samples) {
            sample_names <- append(sample_names, s@sample)
            all_gini <- append(all_gini, s@allstats_qc$gini_coeff)
        }
        names(all_gini) <- sample_names

        lib_gini <- object@stats$gini_coeff_before_qc
        names(lib_gini) <- rownames(object@stats)
        qc_gini <- object@stats$gini_coeff_after_qc
        names(qc_gini) <- rownames(object@stats)

        num_samples <- length(sample_names)
        df_gini <- data.frame(matrix(NA, num_samples * 3, 3))
        colnames(df_gini) <- c("gini", "sample", "type")
        df_gini$gini <- c(all_gini, lib_gini, qc_gini)
        df_gini$sample <- c(names(all_gini), names(lib_gini), names(qc_gini))
        df_gini$type <- c(rep("independent", num_samples), rep("dependent", num_samples), rep("after_qc", num_samples))

        df_gini$gini <- as.numeric(df_gini$gini)
        df_gini$sample <- factor(df_gini$sample, levels = sample_names)
        df_gini$type <- factor(df_gini$type, levels = c("independent", "dependent", "after_qc"))

        gg_colors_fill <- c(t_col("tomato", 0.5), t_col("royalblue", 0.5), t_col("yellowgreen", 0.5))
        gg_colors <- c(c("tomato", "royalblue", "yellowgreen"))
        p1 <- ggplot(df_gini,  aes(x = sample, y = gini, fill = type)) +
                geom_bar(position = "dodge", stat = "identity") +
                scale_fill_manual(values = gg_colors_fill) +
                scale_color_manual(values = gg_colors) +
                labs(x = "samples", y = "score", title = "Sample QC Gini Efficiency") +
                theme(legend.position = "right", legend.title = element_blank()) +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 12, face = "bold")) +
                theme(axis.text.x = element_text(angle = 90)) +
                scale_y_continuous(limits = c(0, 1))

        pwidth <- 150 * num_samples

        if (is.null(plot_dir)) {
            ggplotly(p1)
        } else {
            png(paste0(plot_dir, "/", "sample_qc_gini.png"), width = pwidth, height = 1200, res = 200)
            print(p1)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_samqc_pos_cov", function(object, ...) {
  standardGeneric("qcplot_samqc_pos_cov")
})

#' create the position plot
#'
#' @export
#' @name qcplot_samqc_pos_cov
#' @param object   sampleQC object
#' @param qc_type  plot type, screen or plasmid
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_samqc_pos_cov",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          plot_dir = NULL) {
        if (is.null(plot_dir)) {
            stop(paste0("====> Error: plot_dir is not provided, no output directory."))
        }

        qc_type <- match.arg(qc_type)

        sample_names <- vector()
        libcounts_pos <- data.table()
        for (s in object@samples) {
            sample_names <- append(sample_names, s@sample)

            tmp_counts <- object@library_counts_pos[[s@sample]][, c("sequence", "position", "count")]
            tmp_counts[, sample := s@sample]

            if (nrow(libcounts_pos) == 0) {
                libcounts_pos <- tmp_counts
            } else {
                libcounts_pos <- rbind(libcounts_pos, tmp_counts)
            }
        }
        libcounts_pos[, log2p1 := log2(count + 1)]

        if (qc_type == "plasmid") {
            libcounts_dependent_pos <- data.table()

            for (s in object@samples) {
                tmp_counts <- s@libcounts[, c("name", "count")]
                colnames(tmp_counts) <- c("oligo_name", "count")
                tmp_counts <- as.data.table(tmp_counts)
                tmp_counts[, sample := s@sample]

                tmp_meta <- s@valiant_meta[, c("oligo_name", "mut_position")]
                tmp_meta <- as.data.table(tmp_meta)

                tmp_counts[tmp_meta, position := i.mut_position, on = .(oligo_name)]
                setorder(tmp_counts, cols = "position")

                if (nrow(libcounts_dependent_pos) == 0) {
                libcounts_dependent_pos <- tmp_counts
                } else {
                    libcounts_dependent_pos <- rbind(libcounts_dependent_pos, tmp_counts)
                }
            }

            libcounts_pos <- libcounts_dependent_pos
            libcounts_pos[, log2p1 := log2(count+1)]
        }

        libcounts_pos$sample <- factor(libcounts_pos$sample, levels = mixedsort(levels(factor(libcounts_pos$sample))))
        libcounts_pos_range <- libcounts_pos[, .(min = min(position, na.rm = TRUE), max = max(position, na.rm = TRUE)), by = sample]
        list_scales <- list()
        for (i in 1:nrow(libcounts_pos_range)) {
            list_scales[[i]] <- scale_override(i, scale_x_continuous(breaks = c(libcounts_pos_range$min[i], libcounts_pos_range$max[i])))
        }

        p1 <- ggplot(libcounts_pos, aes(x = position, y = log2p1)) +
                geom_point(shape = 16, size = 0.5, color = "tomato", alpha = 0.8) +
                geom_hline(yintercept = log2(object@cutoffs$seq_low_count+1), linetype = "dashed", color = "springgreen4", linewidth = 0.4) +
                labs(x = "Genomic Coordinate", y = "log2(count+1)", title = "Sample QC position coverage") +
                ylim(0, as.integer(max(libcounts_pos$log2p1)) + 1) +
                theme(legend.position = "none", panel.grid.major = element_blank()) +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 12, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 12, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 8, face = "bold")) +
                theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
                facet_wrap_custom(~sample, scales = "free", scale_overrides = list_scales, ncol = 3)

        pheight <- 400 * ceiling((length(sample_names) / 3))

        if (is.null(plot_dir)) {
            stop(paste0("====> Error: plot_dir is not provided, no output directory."))
        } else {
            png(paste0(plot_dir, "/", "sample_qc_position_cov.dots.png"), width = 2400, height = pheight, res = 200)
            print(p1)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_samqc_pos_anno", function(object, ...) {
  standardGeneric("qcplot_samqc_pos_anno")
})

#' create the position plot
#'
#' @export
#' @name qcplot_samqc_pos_anno
#' @param object    sampleQC object
#' @param samples   a vector of sample names
#' @param type      plot type, lof or all
#' @param plot_dir  the output plot directory
setMethod(
    "qcplot_samqc_pos_anno",
    signature = "sampleQC",
    definition = function(object,
                          samples = NULL,
                          type = "lof",
                          plot_dir = NULL) {
        if (is.null(plot_dir)) {
            stop(paste0("====> Error: plot_dir is not provided, no output directory."))
        }

        if (is.null(samples)) {
            stop(paste0("====> Error: please provide samples, a vector."))
        }

        if (type %nin% c("lof", "all")) {
            stop(paste0("====> Error: wrong type, please use lof or all."))
        }

        libcounts_pos <- as.data.frame(object@library_counts_pos_anno)
        libcounts_pos <- libcounts_pos[, c(samples, "position", "consequence")]
        libcounts_pos_range <- c(min(libcounts_pos$position, na.rm = TRUE), max(libcounts_pos$position, na.rm = TRUE))

        if (type == "lof") {
            libcounts_pos$consequence <- ifelse(libcounts_pos$consequence == "LOF", "LOF", "Others")

            # be careful, df / vec is by row, not column
            libcounts_pos[, samples] <- t(t(libcounts_pos[, samples]) / object@stats[samples, ]$accepted_reads * 100)

            df_libcounts_pos <- reshape2::melt(libcounts_pos, id.vars = c("consequence", "position"), variable.name = "samples", value.name = "counts")
            df_libcounts_pos$samples <- factor(df_libcounts_pos$samples, levels = samples)

            tmp_cutoff <- object@cutoffs$low_abundance_per * 100

            p1 <- ggplot(df_libcounts_pos, aes(x = position, y = counts)) +
                    geom_point(shape = 19, size = 0.5, aes(color = factor(consequence))) +
                    geom_hline(yintercept = tmp_cutoff, linetype = "dashed", color = "springgreen4", linewidth = 0.4) +
                    scale_color_manual(values = c(t_col("red", 1), t_col("royalblue", 0.2)), labels = c("LOF", "Others")) +
                    labs(x = "Genomic Coordinate", y = "Percentage", title = "Sample QC position percentage", color = "Type") +
                    scale_x_continuous(limits = libcounts_pos_range, breaks = libcounts_pos_range) +
                    coord_transform(y = "log2") +
                    scale_y_continuous(breaks = c(0.005, 0.01, 0.05, 0.2, 0.5, 1)) +
                    theme(legend.position = "right", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 12, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 12, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 8, face = "bold")) +
                    facet_wrap(~samples, scales = "free_x",  dir = "v", ncol = 1)

            pheight <- 400 * length(samples)

            if (is.null(plot_dir)) {
                stop(paste0("====> Error: plot_dir is not provided, no output directory."))
            } else {
                png(paste0(plot_dir, "/", "sample_qc_position_anno.lof_dots.png"), width = 1200, height = pheight, res = 200)
                print(p1)
                dev.off()
            }
        } else {
            libcounts_pos[, samples] <- t(t(libcounts_pos[, samples]) / object@stats[samples, ]$accepted_reads * 100)

            df_libcounts_pos <- reshape2::melt(libcounts_pos, id.vars = c("consequence", "position"), variable.name = "samples", value.name = "counts")
            df_libcounts_pos$samples <- factor(df_libcounts_pos$samples, levels = samples)

            df_libcounts_pos[df_libcounts_pos == 0] <- NA

            num_colors <- length(unique(libcounts_pos$consequence))
            index_colors <- sample(seq(1, length(select_colorblind("col15"))), num_colors)
            select_colors <- select_colorblind("col15")[index_colors]

            freq_cons <- table(libcounts_pos$consequence)
            names(select_colors) <- names(freq_cons)

            freq_cons <- sort(freq_cons, decreasing = TRUE)
            freq_cons <- names(freq_cons)
            rate_cons <- seq(0.2, 0.1 + length(freq_cons)/10, 0.1)
            names(rate_cons) <- freq_cons

            for (i in 1:(length(select_colors) - 1)) {
                select_colors[i] <- t_col(select_colors[i], rate_cons[names(select_colors[i])])
            }
            select_colors <- as.vector(select_colors)

            tmp_cutoff <- object@cutoffs$low_abundance_per * 100

            p1 <- ggplot(df_libcounts_pos, aes(x = position, y = counts)) +
                    geom_point(shape = 19, size = 0.5, aes(color = factor(consequence))) +
                    geom_hline(yintercept = tmp_cutoff, linetype = "dashed", color = "springgreen4", linewidth = 0.4) +
                    scale_color_manual(values = select_colors) +
                    labs(x = "Genomic Coordinate", y = "Percentage", title = "Sample QC position percentage", color = "Type") +
                    scale_x_continuous(limits = libcounts_pos_range, breaks = libcounts_pos_range) +
                    coord_transform(y = "log2") +
                    scale_y_continuous(breaks = c(0.005, 0.01, 0.05, 0.2, 0.5, 1)) +
                    theme(legend.position = "right", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 12, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 12, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 8, face = "bold")) +
                    facet_wrap(~samples, scales = "free_x", dir = "v", ncol = 1)

            pheight <- 400 * length(samples)

            if (is.null(plot_dir)) {
                stop(paste0("====> Error: plot_dir is not provided, no output directory."))
            } else {
                png(paste0(plot_dir, "/", "sample_qc_position_anno.all_dots.png"), width = 1200, height = pheight, res = 200)
                print(p1)
                dev.off()
            }
        }
    }
)

#####################################################################################################################################################

#' initialize function
setGeneric("qcplot_expqc_all", function(object, ...) {
  standardGeneric("qcplot_expqc_all")
})

#' create all the plot figures
#'
#' @export
#' @name qcplot_expqc_all
#' @param object   sampleQC object
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_expqc_all",
    signature = "experimentQC",
    definition = function(object,
                          plot_dir = NULL) {
        if (is.null(plot_dir)) {
            stop(paste0("====> Error: plot_dir is not provided, no output directory."))
        }

        qcplot_expqc_sample_corr(object = object, plot_dir = plot_dir)
        qcplot_expqc_sample_pca(object = object, plot_dir = plot_dir)
        
	    # Use positional-corrected RAW values for "all"
        qcplot_expqc_deseq_fc(
            object    = object,
            eqc_type  = "all",
            plot_type = "beeswarm",
            plot_dir  = plot_dir,
            stat_col  = "stat_pos_raw",
            lfc_col   = "pos_adj_log2FoldChange_raw"
        )

        qcplot_expqc_deseq_fc_pos(
            object   = object,
            eqc_type = "all",
            plot_dir = plot_dir,
            stat_col = "stat_pos_raw",
            lfc_col  = "pos_adj_log2FoldChange_raw"
        )
        
        qcplot_expqc_positional_loess_diag(object = object, plot_dir = plot_dir, use = "raw")  # targeton label not available here; call from main() for titled output
    }
)

#' initialize function
setGeneric("qcplot_expqc_sample_corr", function(object, ...) {
  standardGeneric("qcplot_expqc_sample_corr")
})

#' create the heatmap of samples
#'
#' @export
#' @name qcplot_expqc_sample_corr
#' @param object   experimentQC object
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_expqc_sample_corr",
    signature = "experimentQC",
    definition = function(object,
                          plot_dir = NULL) {
        sample_dend <- as.dendrogram(object@lib_hclust_res)

        num_clusters <- length(unique(object@coldata$condition))
        name_clusters <- as.vector(unique(object@coldata$condition))

        sample_dend <- dendextend::set(sample_dend, "branches_lwd", 1)
        sample_dend <- dendextend::set(sample_dend, "branches_k_color", select_colorblind("col8")[1:num_clusters], k = num_clusters)
        sample_dend <- dendextend::set(sample_dend, "labels_cex", 0.6)
        sample_dend <- dendextend::set(sample_dend, "labels_colors", select_colorblind("col8")[1:num_clusters], k = num_clusters)

        pheight <- 50 * length(object@lib_hclust_res$labels)
        png(paste0(plot_dir, "/", "experiment_qc_samples_tree.png"), width = 800, height = pheight, res = 200)
        par(mar = c(1, 1, 1, 5))
        plot(sample_dend, axes = FALSE, horiz = TRUE)
        dev.off()

        sample_rlog <- as.matrix(object@lib_deseq_rlog)

        min_rlog <- round(min(sample_rlog))
        max_rlog <- round(max(sample_rlog))

        sample_corr <- cor(scale(sample_rlog))
        min_corr <- floor(min(sample_corr) * 10) / 10

        p <- ggcorrplot(sample_corr,
                        method = "square",
	                    hc.method = "ward.D2",
                        hc.order = TRUE,
  		                lab = TRUE,
  		                lab_col = "black",
  		                lab_size = 3,
		                p.mat = cor_pmat(sample_corr),
		                sig.level = 0.05,
                        tl.col = "black",
                        tl.cex = 12)
        p1 <- p + scale_fill_gradient2(limit = c(min_corr, 1),
                                       low = "royalblue",
                                       high =  "red",
                                       mid = "ivory",
                                       midpoint = (1 + min_corr) / 2,
                                       name = "correlation")

        if (is.null(plot_dir)) {
            ggplotly(p1)
        } else {
            png(paste0(plot_dir, "/", "experiment_qc_samples_corr.png"), width = 1200, height = 1200, res = 200)
            corrplot(sample_corr,
                     method = "color",
                     order = "hclust",
                     col = colorpanel(100, "royalblue", "ivory", "red"),
                     col.lim = c(min_corr, 1),
                     is.corr = FALSE,
                     addrect = 3,
                     rect.col = "black",
                     rect.lwd = 1.5,
                     addgrid.col = "white",
                     tl.col = "black",
                     tl.cex = 0.75,
                     addCoef.col = "black",
                     number.cex = 0.75)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_expqc_sample_pca", function(object, ...) {
  standardGeneric("qcplot_expqc_sample_pca")
})

#' create the pca of samples
#'
#' @export
#' @name qcplot_expqc_sample_pca
#' @param object     experimentQC object
#' @param plot_dir   the output plot directory
setMethod(
    "qcplot_expqc_sample_pca",
    signature = "experimentQC",
    definition = function(object,
                      plot_dir = NULL,
                      extra_pca_exclude_condition = NULL,
                      extra_pca_auto_when_ref_is  = "Plasmid",
                      extra_pca_ntop = 500) {
        if (is.null(plot_dir)) {
            stop(paste0("====> Error: plot_dir is not provided, no output directory."))
        }

        pca <- object@lib_pca_res
        percentVar <- pca$sdev^2 / sum(pca$sdev^2)
        percentVar <- round(percentVar, digits = 3) * 100

        pc1_set <- c((min(pca$x[, 1]) - sd(pca$x[, 1])), (max(pca$x[, 1]) + sd(pca$x[, 1])))
        pc2_set <- c((min(pca$x[, 2]) - sd(pca$x[, 2])), (max(pca$x[, 2]) + sd(pca$x[, 2])))
        pc3_set <- c((min(pca$x[, 3]) - sd(pca$x[, 3])), (max(pca$x[, 3]) + sd(pca$x[, 3])))

        ds_coldata <- object@coldata
        # mark conditions
        default_colors <- c("tomato", "royalblue", "yellowgreen", "orange", "pink", "purple", "coral", "cyan")
        select_colors <- default_colors[1:length(levels(ds_coldata$condition))]
        names(select_colors) <- levels(ds_coldata$condition)

        pca_colors <- 1:nrow(ds_coldata)
        for (i in 1:nrow(ds_coldata)) {
            pca_colors[i] <- select_colors[ds_coldata[i, ]$condition]
        }

        pca_bgs <- sapply(pca_colors, function(x) t_col(x, 0.5))

        # mark replicates
        default_pchs <- c(21, 22, 23, 24, 25)
        select_pchs <- default_pchs[1:length(levels(ds_coldata$replicate))]
        names(select_pchs) <- levels(ds_coldata$replicate)

        pca_pchs <- 1:nrow(ds_coldata)
        for (i in 1:nrow(ds_coldata)) {
            pca_pchs[i] <- select_pchs[ds_coldata[i, ]$replicate]
        }

        png(paste0(plot_dir, "/", "experiment_qc_pca_samples.png"), width = 1200, height = 1200, res = 200)
        par(mfrow = c(2, 2), mar = c(4, 4, 4, 1))
        plot(pca$x[, 1], pca$x[, 2], xlab = "PC1", ylab = "PC2", pch = pca_pchs, col = pca_colors, bg = pca_bgs, lwd = 1, cex = 2, xlim = pc1_set, ylim = pc2_set, main = "PC1 vs PC2")
        plot(pca$x[, 2], pca$x[, 3], xlab = "PC2", ylab = "PC3", pch = pca_pchs, col = pca_colors, bg = pca_bgs, lwd = 1, cex = 2, xlim = pc2_set, ylim = pc3_set, main = "PC2 vs PC3")
        plot(pca$x[, 1], pca$x[, 3], xlab = "PC1", ylab = "PC3", pch = pca_pchs, col = pca_colors, bg = pca_bgs, lwd = 1, cex = 2, xlim = pc1_set, ylim = pc3_set, main = "PC1 vs PC3")
        b <- barplot(percentVar, col = t_col("royalblue", 0.5), border = "royalblue", ylim = c(0, 105))
        text(b, percentVar + 5, paste0(percentVar, "%"), cex = 0.6)
        legend("topright", legend = levels(ds_coldata$replicate), pch = select_pchs, cex = 1, bty = "n")
        legend("top", legend = levels(ds_coldata$condition), pch = 19, col = select_colors, cex = 1, bty = "n")
        dev.off()
        
		# ----------------------------
		# Optional: extra PCA excluding a condition (e.g. Plasmid)
		# ----------------------------
		exclude_cond <- extra_pca_exclude_condition

		# auto-trigger when reference is plasmid (customizable)
		if (is.null(exclude_cond) && !is.null(extra_pca_auto_when_ref_is)) {
		  if (tolower(object@ref_condition) == tolower(extra_pca_auto_when_ref_is)) {
			exclude_cond <- object@ref_condition
		  }
		}

		if (!is.null(exclude_cond)) {

		  ds_coldata2 <- object@coldata
		  keep_samples <- rownames(ds_coldata2)[ds_coldata2$condition != exclude_cond]

		  if (length(keep_samples) >= 3) {  # need >=3 samples to make PCA meaningful
			pca_input2 <- as.matrix(object@lib_deseq_rlog)[, keep_samples, drop = FALSE]

			rv2 <- matrixStats::rowVars(pca_input2)
			select2 <- order(rv2, decreasing = TRUE)[seq_len(min(extra_pca_ntop, length(rv2)))]

			pca2 <- prcomp(t(pca_input2[select2, , drop = FALSE]), center = TRUE, scale. = TRUE)
			percentVar2 <- pca2$sdev^2 / sum(pca2$sdev^2)
			percentVar2 <- round(percentVar2, digits = 3) * 100

			# subset coldata for plotting
			ds_coldata2 <- ds_coldata2[keep_samples, , drop = FALSE]

			# colors by condition (subset)
			default_colors <- c("tomato", "royalblue", "yellowgreen", "orange", "pink", "purple", "coral", "cyan")
			select_colors2 <- default_colors[1:length(levels(ds_coldata2$condition))]
			names(select_colors2) <- levels(ds_coldata2$condition)

			pca_colors2 <- sapply(ds_coldata2$condition, function(x) select_colors2[[as.character(x)]])
			pca_bgs2 <- sapply(pca_colors2, function(x) t_col(x, 0.5))

			# shapes by replicate (subset)
			default_pchs <- c(21, 22, 23, 24, 25)
			select_pchs2 <- default_pchs[1:length(levels(ds_coldata2$replicate))]
			names(select_pchs2) <- levels(ds_coldata2$replicate)
			pca_pchs2 <- sapply(ds_coldata2$replicate, function(x) select_pchs2[[as.character(x)]])

			pc1_set2 <- c((min(pca2$x[, 1]) - sd(pca2$x[, 1])), (max(pca2$x[, 1]) + sd(pca2$x[, 1])))
			pc2_set2 <- c((min(pca2$x[, 2]) - sd(pca2$x[, 2])), (max(pca2$x[, 2]) + sd(pca2$x[, 2])))
			pc3_set2 <- c((min(pca2$x[, 3]) - sd(pca2$x[, 3])), (max(pca2$x[, 3]) + sd(pca2$x[, 3])))

			out <- paste0(plot_dir, "/", "experiment_qc_pca_samples_excluding_", exclude_cond, ".png")
			png(out, width = 1200, height = 1200, res = 200)
			par(mfrow = c(2, 2), mar = c(4, 4, 4, 1))
			plot(pca2$x[, 1], pca2$x[, 2], xlab = "PC1", ylab = "PC2",
				 pch = pca_pchs2, col = pca_colors2, bg = pca_bgs2, lwd = 1, cex = 2,
				 xlim = pc1_set2, ylim = pc2_set2, main = paste0("PC1 vs PC2 (no ", exclude_cond, ")"))
			plot(pca2$x[, 2], pca2$x[, 3], xlab = "PC2", ylab = "PC3",
				 pch = pca_pchs2, col = pca_colors2, bg = pca_bgs2, lwd = 1, cex = 2,
				 xlim = pc2_set2, ylim = pc3_set2, main = paste0("PC2 vs PC3 (no ", exclude_cond, ")"))
			plot(pca2$x[, 1], pca2$x[, 3], xlab = "PC1", ylab = "PC3",
				 pch = pca_pchs2, col = pca_colors2, bg = pca_bgs2, lwd = 1, cex = 2,
				 xlim = pc1_set2, ylim = pc3_set2, main = paste0("PC1 vs PC3 (no ", exclude_cond, ")"))
			b2 <- barplot(percentVar2, col = t_col("royalblue", 0.5), border = "royalblue", ylim = c(0, 105))
			text(b2, percentVar2 + 5, paste0(percentVar2, "%"), cex = 0.6)
			legend("topright", legend = levels(ds_coldata2$replicate), pch = select_pchs2, cex = 1, bty = "n")
			legend("top", legend = levels(ds_coldata2$condition), pch = 19, col = select_colors2, cex = 1, bty = "n")
			dev.off()

		  } else {
			warning(sprintf("Skipping extra PCA excluding '%s' (need at least 3 samples after exclusion).", exclude_cond))
		  }
		}
        
    }
)

#' initialize function
setGeneric("qcplot_expqc_deseq_fc", function(object, ...) {
  standardGeneric("qcplot_expqc_deseq_fc")
})

#' create fold change and consequence plot
#'
#' @export
#' @name qcplot_expqc_deseq_fc
#' @param object     experimentQC object
#' @param eqc_type   library counts or all counts 
#' @param cons       a vector of the selected consequences in the vep annotation file
#' @param plot_type  beeswarm or violin
#' @param plot_dir   the output plot directory
#' @param stat_col   name of the status column to use (e.g. "stat_pos_raw")
#' @param lfc_col    name of the log2FC column to use (e.g. "pos_adj_log2FoldChange_raw")
setMethod(
    "qcplot_expqc_deseq_fc",
    signature = "experimentQC",
    definition = function(object,
                          eqc_type  = c("lib", "all"),
                          cons      = c("Synonymous_Variant", "LOF", "Missense_Variant"),
                          plot_type = c("beeswarm", "violin"),
                          plot_dir  = NULL,
                          stat_col  = NULL,
                          lfc_col   = NULL) {

        if (length(plot_dir) == 0) {
            stop("====> Error: plot_dir is not provided, no output directory.")
        }

        eqc_type  <- match.arg(eqc_type)
        plot_type <- match.arg(plot_type)
        
        # pick correct columns automatically
		if (eqc_type == "lib") {
		  if (is.null(stat_col)) stat_col <- "stat_raw"
		  if (is.null(lfc_col))  lfc_col  <- "log2FoldChange_raw"
		} else {
		  # eqc_type == "all" -> positional-corrected RAW by default
		  if (is.null(stat_col)) stat_col <- "stat_pos_raw"
		  if (is.null(lfc_col))  lfc_col  <- "pos_adj_log2FoldChange_raw"
		}

        if (eqc_type == "lib") {
            comparisons <- names(object@lib_deseq_res_anno)
            df_list     <- object@lib_deseq_res_anno
        } else {
            comparisons <- names(object@all_deseq_res_anno_adj)
            df_list     <- object@all_deseq_res_anno_adj
            
            # drop non-contrast diagnostic entries added by positional correction
            keep <- !grepl("^positional_loess_", names(df_list))
            df_list <- df_list[keep]
            comparisons <- names(df_list)

        }

        #-----------------------------#
        # 1) Y-limits from chosen LFC #
        #-----------------------------#
        ylimits <- vector()
        
        for (i in seq_along(df_list)) {
            res      <- df_list[[i]]
            res_cons <- res[res$consequence %in% cons]

            if (!nrow(res_cons)) next

            # ensure the chosen columns exist
            stopifnot(stat_col %in% colnames(res_cons))
            stopifnot(lfc_col  %in% colnames(res_cons))

            # alias to canonical names
            res_cons$stat          <- res_cons[[stat_col]]
            res_cons$log2FoldChange <- res_cons[[lfc_col]]

            # factor with expected levels
            res_cons$stat <- factor(res_cons$stat,
                                    levels = c("no impact", "enriched", "depleted"))

            # drop NAs in the chosen LFC
            res_cons <- res_cons[!is.na(res_cons$log2FoldChange), ]

            if (!nrow(res_cons)) next

            ylimits <- append(ylimits, ceiling(max(res_cons$log2FoldChange)))
            ylimits <- append(ylimits, floor(min(res_cons$log2FoldChange)))
        }

        if (!length(ylimits)) {
            warning("No data available for DESeq FC QC plot with given cons/stat/lfc settings.")
            return(invisible(NULL))
        }

        ylimits <- sort(ylimits)
        ymin    <- head(ylimits, n = 1)
        ymax    <- tail(ylimits, n = 1)

        # user overrides from config
        ymin <- maveqc_config$expqc_lfc_min
        ymax <- maveqc_config$expqc_lfc_max

        #-----------------------------#
        # 2) Plot per comparison      #
        #-----------------------------#
        for (i in seq_along(df_list)) {
            res      <- df_list[[i]]
            res_cons <- res[res$consequence %in% cons]

            if (!nrow(res_cons)) next

            stopifnot(stat_col %in% colnames(res_cons))
            stopifnot(lfc_col  %in% colnames(res_cons))

            res_cons$stat           <- res_cons[[stat_col]]
            res_cons$log2FoldChange <- res_cons[[lfc_col]]
            res_cons$stat <- factor(res_cons$stat,
                                    levels = c("no impact", "enriched", "depleted"))
            res_cons <- res_cons[!is.na(res_cons$log2FoldChange), ]

            if (!nrow(res_cons)) next

            stat_unique <- unique(res_cons$stat)
            stat_level  <- levels(res_cons$stat)
            stat_size   <- c(0.5, 1, 1)
            stat_color  <- c(t_col("black", 0.4),
                             t_col("red", 0.8),
                             t_col("yellowgreen", 0.8))

            stat_size_plot  <- vector()
            stat_color_plot <- vector()
            for (j in seq_along(stat_level)) {
                if (stat_level[j] %in% stat_unique) {
                    stat_size_plot  <- append(stat_size_plot,  stat_size[j])
                    stat_color_plot <- append(stat_color_plot, stat_color[j])
                }
            }

            if (plot_type == "beeswarm") {
                p1 <- ggplot(res_cons, aes(x = consequence, y = log2FoldChange)) +
                    geom_violin(trim = FALSE, scale = "width",
                                fill = t_col("lightblue", 0.5),
                                color = "royalblue") +
                    geom_quasirandom(width = 0.4,
                                     aes(color = factor(stat),
                                         size  = factor(stat))) +
                    scale_color_manual(values = stat_color_plot) +
                    scale_size_manual(values  = stat_size_plot) +
                    ylim(ymin, ymax) +
                    coord_flip() +
                    labs(x = "log2FoldChange",
                         title = comparisons[i],
                         color = "Type") +
                    theme(legend.position = "right",
                          panel.grid.major = element_blank(),
                          panel.background = element_rect(fill = "ivory", colour = "white"),
                          axis.title.y = element_blank(),
                          axis.title.x = element_text(size = 12, face = "bold", family = "Arial"),
                          plot.title   = element_text(size = 12, face = "bold.italic", family = "Arial"),
                          axis.text    = element_text(size = 8, face = "bold")) +
                    guides(size = "none")
            } else {
                p1 <- ggplot(res_cons, aes(x = consequence, y = log2FoldChange)) +
                    geom_violinhalf(trim = FALSE, scale = "width",
                                    fill = t_col("lightblue", 0.5),
                                    color = "royalblue",
                                    position = position_nudge(x = .2, y = 0)) +
                    geom_jitter(width = 0.15,
                                aes(color = factor(stat),
                                    size  = factor(stat))) +
                    scale_color_manual(values = stat_color_plot) +
                    scale_size_manual(values  = stat_size_plot) +
                    ylim(ymin, ymax) +
                    coord_flip() +
                    labs(x = "log2FoldChange",
                         title = comparisons[i],
                         color = "Type") +
                    theme(legend.position = "right",
                          panel.grid.major = element_blank(),
                          panel.background = element_rect(fill = "ivory", colour = "white"),
                          axis.title.y = element_blank(),
                          axis.title.x = element_text(size = 12, face = "bold", family = "Arial"),
                          plot.title   = element_text(size = 12, face = "bold.italic", family = "Arial"),
                          axis.text    = element_text(size = 8, face = "bold")) +
                    guides(size = "none")
            }

            pheight <- 200 * length(cons)
            file_path <- paste0(plot_dir, "/",
                                "experiment_qc_deseq_fc.",
                                comparisons[i], ".",
                                eqc_type, "_", plot_type, ".png")

            png(file_path, width = 1500, height = pheight, res = 200)
            print(p1)
            dev.off()
        }
    }
)



#' initialize function
setGeneric("qcplot_expqc_deseq_fc_pos", function(object, ...) {
  standardGeneric("qcplot_expqc_deseq_fc_pos")
})

#' create fold change and consequence plot (by position)
#'
#' @export
#' @name qcplot_expqc_deseq_fc_pos
#' @param object    experimentQC object
#' @param eqc_type  library counts or all counts
#' @param cons      a vector of all the consequences in the vep annotation file
#' @param plot_dir  the output plot directory
#' @param stat_col  name of the status column to use (e.g. "stat_pos_raw")
#' @param lfc_col   name of the log2FC column to use (e.g. "pos_adj_log2FoldChange_raw")
setMethod(
    "qcplot_expqc_deseq_fc_pos",
    signature = "experimentQC",
    definition = function(object,
                          eqc_type = c("lib", "all"),
                          cons     = c("Synonymous_Variant",
                                       "LOF",
                                       "Missense_Variant",
                                       "Intronic_Variant",
                                       "Inframe_Deletion",
                                       "Splice_Variant",
                                       "Splice_Polypyrimidine_Tract_Variant",
                                       "Others"),
                          plot_dir = NULL,
                          stat_col = NULL,
                          lfc_col  = NULL) {

        if (length(plot_dir) == 0) {
            stop("====> Error: plot_dir is not provided, no output directory.")
        }

        eqc_type <- match.arg(eqc_type)

		# pick correct columns automatically
		if (eqc_type == "lib") {
		  if (is.null(stat_col)) stat_col <- "stat_raw"
		  if (is.null(lfc_col))  lfc_col  <- "log2FoldChange_raw"
		} else {
		  # eqc_type == "all" -> positional-corrected RAW by default
		  if (is.null(stat_col)) stat_col <- "stat_pos_raw"
		  if (is.null(lfc_col))  lfc_col  <- "pos_adj_log2FoldChange_raw"
		}

        if (eqc_type == "lib") {
            comparisons <- names(object@lib_deseq_res_anno)
            df_list     <- object@lib_deseq_res_anno
        } else {
            comparisons <- names(object@all_deseq_res_anno_adj)
            df_list     <- object@all_deseq_res_anno_adj
            
			# drop non-contrast diagnostic entries added by positional correction
			keep <- !grepl("^positional_loess_", names(df_list))
			df_list <- df_list[keep]
			comparisons <- names(df_list)
        }

        colors        <- select_colorblind("col15")[1:length(cons)]
        select_colors <- sapply(colors, function(x) t_col(x, 0.3), USE.NAMES = FALSE)
        fill_colors   <- sapply(colors, function(x) t_col(x, 0.8), USE.NAMES = FALSE)

        #-----------------------------#
        # 1) Y-limits from chosen LFC #
        #-----------------------------#
        ylimits <- vector()
        for (i in seq_along(df_list)) {
            dt_res <- df_list[[i]]
            if (!nrow(dt_res)) next

            stopifnot(stat_col %in% colnames(dt_res))
            stopifnot(lfc_col  %in% colnames(dt_res))

            dt_res$stat           <- dt_res[[stat_col]]
            dt_res$log2FoldChange <- dt_res[[lfc_col]]
            dt_res <- dt_res[!is.na(dt_res$log2FoldChange), ]

            if (!nrow(dt_res)) next

            ylimits <- append(ylimits, ceiling(max(dt_res$log2FoldChange)))
            ylimits <- append(ylimits, floor(min(dt_res$log2FoldChange)))
        }

        if (!length(ylimits)) {
            warning("No data available for DESeq FC position QC plot with given stat/lfc settings.")
            return(invisible(NULL))
        }

        ylimits <- sort(ylimits)
        ymin    <- head(ylimits, n = 1)
        ymax    <- tail(ylimits, n = 1)

        # user override
        ymin <- maveqc_config$expqc_lfc_min
        ymax <- maveqc_config$expqc_lfc_max

        #-----------------------------#
        # 2) Plot per comparison      #
        #-----------------------------#
        for (i in seq_along(df_list)) {
            dt_res <- df_list[[i]]
            if (!nrow(dt_res)) next

            stopifnot(stat_col %in% colnames(dt_res))
            stopifnot(lfc_col  %in% colnames(dt_res))

            dt_res$consequence <- factor(dt_res$consequence, levels = cons)

            dt_res$stat           <- dt_res[[stat_col]]
            dt_res$log2FoldChange <- dt_res[[lfc_col]]

            dt_res <- dt_res[!is.na(dt_res$log2FoldChange), ]

            if (!nrow(dt_res)) next

            pos_tmp <- unique(sort(dt_res$position))
            pos_min <- head(pos_tmp, n = 1)
            pos_max <- tail(pos_tmp, n = 1)
            pos_by  <- floor((pos_max - pos_min) / 5)

            dt_res$stat <- factor(dt_res$stat,
                                  levels = c("no impact", "enriched", "depleted"))
            stat_unique <- unique(dt_res$stat)
            stat_level  <- levels(dt_res$stat)
            stat_size   <- c(0.5, 2, 2)
            stat_shape  <- c(16, 24, 25)

            stat_size_plot  <- vector()
            stat_shape_plot <- vector()
            for (j in seq_along(stat_level)) {
                if (stat_level[j] %in% stat_unique) {
                    stat_size_plot  <- append(stat_size_plot,  stat_size[j])
                    stat_shape_plot <- append(stat_shape_plot, stat_shape[j])
                }
            }

            p1 <- ggplot(dt_res, aes(x = position, y = log2FoldChange)) +
                geom_point(aes(size  = factor(stat),
                               shape = factor(stat),
                               fill  = factor(consequence),
                               color = factor(consequence))) +
                scale_size_manual(values  = stat_size_plot) +
                scale_shape_manual(values = stat_shape_plot) +
                scale_color_manual(values = select_colors) +
                scale_fill_manual(values  = fill_colors) +
                labs(x = "Genomic Coordinate",
                     y = "Log2 Fold Change",
                     title = comparisons[i]) +
                theme(legend.position = "right",
                      legend.title    = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.background = element_rect(fill = "ivory", colour = "white"),
                      axis.title      = element_text(size = 12, face = "bold", family = "Arial"),
                      plot.title      = element_text(size = 12, face = "bold.italic", family = "Arial"),
                      axis.text       = element_text(size = 8, face = "bold")) +
                scale_x_continuous(limits = c(pos_min, pos_max),
                                   breaks = seq(pos_min, pos_max, pos_by)) +
                scale_y_continuous(limits = c(ymin, ymax),
                                   breaks = seq(ymin, ymax)) +
                guides(fill = guide_legend(override.aes = list(shape = 21)))

            file_path <- paste0(plot_dir, "/",
                                "experiment_qc_deseq_fc.",
                                comparisons[i], ".",
                                eqc_type, "_position.png")

            png(file_path, width = 1500, height = 1000, res = 200)
            print(p1)
            dev.off()
        }
    }
)


#' initialize function
setGeneric("qcplot_expqc_positional_loess_diag", function(object, ...) {
  standardGeneric("qcplot_expqc_positional_loess_diag")
})

#' Positional LOESS diagnostics: fit points + LOESS curve (ggplot)
#'
#' @export
#' @name qcplot_expqc_positional_loess_diag
#' @param object      experimentQC object
#' @param plot_dir    output plot directory
#' @param targeton    targeton name for plot title (parsed from input_dir in main)
#' @param use         "raw" or "shrunk" for the y-values plotted for fit points
setMethod(
  "qcplot_expqc_positional_loess_diag",
  signature = "experimentQC",
  definition = function(object,
                        plot_dir    = NULL,
                        targeton    = NULL,
                        use         = c("raw", "shrunk")) {

    if (is.null(plot_dir)) {
      stop("====> Error: plot_dir is not provided, no output directory.")
    }
    use <- match.arg(use)

    fit_dt  <- object@all_deseq_res_anno_adj[["positional_loess_fit"]]
    fit_pts <- object@all_deseq_res_anno_adj[["positional_loess_fit_points"]]
    meta    <- object@all_deseq_res_anno_adj[["positional_loess_meta"]]

    if (is.null(fit_dt) || is.null(fit_pts)) {
      warning("No stored LOESS diagnostics found (positional_loess_fit / fit_points). Skipping LOESS diagnostic PNG.")
      return(invisible(NULL))
    }

    fit_dt  <- as.data.table(fit_dt)
    fit_pts <- as.data.table(fit_pts)

    # Pick which fit-point LFC column to use as y
    ycol_fit <- if (use == "raw") {
      if ("lfc_t1_raw" %in% names(fit_pts)) "lfc_t1_raw"
      else if ("lfc_d4_raw" %in% names(fit_pts)) "lfc_d4_raw"
      else stop("fit_points missing expected column: 'lfc_t1_raw' or 'lfc_d4_raw'")
    } else {
      if ("lfc_t1_shr" %in% names(fit_pts)) "lfc_t1_shr"
      else if ("lfc_d4_shr" %in% names(fit_pts)) "lfc_d4_shr"
      else stop("fit_points missing expected column: 'lfc_t1_shr' or 'lfc_d4_shr'")
    }

    setorder(fit_dt, position)

    # Build plot title from targeton name and reference condition
    ref_condition <- if (!is.null(meta) && "ref_condition" %in% names(meta)) {
      meta$ref_condition[1]
    } else {
      "unknown"
    }
    targeton_label <- if (!is.null(targeton)) targeton else "unknown targeton"
    plot_title <- paste0(targeton_label, " LOESS positional fit - ", ref_condition, " reference")

    # Rename y column for clean aes mapping
    fit_pts[, lfc_plot := get(ycol_fit)]

    p <- ggplot() +
      geom_point(data = fit_pts,
                 aes(x = position, y = lfc_plot),
                 colour = "grey40", size = 0.8, alpha = 0.6) +
      geom_line(data = fit_dt,
                aes(x = position, y = pos_fit),
                colour = "black", linewidth = 1) +
      coord_cartesian(ylim = c(-3, 3)) +
      labs(
        title = plot_title,
        x     = "Genomic position",
        y     = "t1 raw LFC"
      ) +
      theme_bw() +
      theme(
        plot.title   = element_text(size = 12, face = "bold"),
        axis.title   = element_text(size = 10),
        axis.text    = element_text(size = 8),
        panel.grid.minor = element_blank()
      )

    out <- file.path(plot_dir, paste0("experiment_qc_positional_loess_fitpoints_",
                                      ref_condition, ".png"))

    ggsave(out, plot = p, width = 8, height = 5, dpi = 200)

    invisible(out)
  }
)

#' barplot in reactable
#'
#' @name bar_style
#' @param length   the length of bar
#' @param height   the height of bar
#' @param fill     the bar color
#' @param align    the alignment of bar
#' @param color    the font color
#' @param fweight  the font weight
bar_style <- function(length = 1,
                      height = "80%",
                      fill = "#00FFFF7F",
                      align = c("right", "left"),
                      color = "black",
                      fweight = "plain") {
    align <- match.arg(align)
    if (align == "left") {
        position <- paste0(length * 100, "%")
        image <- sprintf("linear-gradient(90deg, %1$s %2$s, transparent %2$s)", fill, position)
    } else {
        position <- paste0(100 - length * 100, "%")
        image <- sprintf("linear-gradient(90deg, transparent %1$s, %2$s %1$s)", position, fill)
    }
    list(
        backgroundImage = image,
        backgroundSize = paste("100%", height),
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center",
        color = color,
        fontWeight = fweight
    )
}

#' barplot in reactable
#'
#' @name bar_chart_pos_neg
#' @param lable   the length of bar
#' @param value   the height of bar
#' @param max_value     the bar color
#' @param height    the alignment of bar
#' @param pos_fill    the font color
#' @param neg_fill  the font weight
bar_chart_pos_neg <- function(label,
                              value,
                              max_value = 1,
                              height = "1rem",
                              pos_fill = t_col("tomato", 0.8),
                              neg_fill = t_col("yellowgreen", 0.8)) {
    neg_chart <- div(style = list(flex = "1 1 0"))
    pos_chart <- div(style = list(flex = "1 1 0"))
    width <- paste0(abs(value / max_value) * 100, "%")

    if (value < 0) {
        bar <- div(style = list(marginLeft = "0.5rem", background = neg_fill, width = width, height = height))
        chart <- div(style = list(display = "flex", alignItems = "center", justifyContent = "flex-end"),
                     label,
                     bar)
        neg_chart <- tagAppendChild(neg_chart, chart)
    } else {
        bar <- div(style = list(marginRight = "0.5rem", background = pos_fill, width = width, height = height))
        chart <- div(style = list(display = "flex", alignItems = "center"), bar, label)
        pos_chart <- tagAppendChild(pos_chart, chart)
    }

    div(style = list(display = "flex"), neg_chart, pos_chart)
}

#' create QC reports
#'
#' @export
#' @name create_qc_reports
#' @param samplesheet    the path of sample sheet file
#' @param qc_type         screen or plasmid
#' @param qc_dir          the directory of QC plots and outs
create_qc_reports <- function(samplesheet = NULL,
                              qc_type = c("plasmid", "screen"),
                              qc_dir = NULL,
                              experiment_qc_dir = NULL) {
        #----------#
        # checking #
        #----------#
        if (is.null(samplesheet)) {
            stop(paste0("====> Error: please provide the path of sample sheet file!"))
        }

        if (is.null(qc_dir)) {
            stop(paste0("====> Error: qc_dir is not provided, no output directory."))
        }

        # Resolve to absolute paths up front. rmarkdown::render() can change
        # the working directory while knitting, so any relative path we bake
        # into the generated .Rmd (qc_dir, experiment_qc_dir, samplesheet)
        # is fragile unless it's made absolute here, before we ever write it
        # into a cat() line. Fail fast with a clear message if something
        # doesn't actually exist, rather than erroring deep inside the knit.
        if (!file.exists(samplesheet)) {
            stop(paste0("====> Error: samplesheet not found: ", samplesheet))
        }
        samplesheet <- normalizePath(samplesheet, mustWork = TRUE)

        if (!dir.exists(qc_dir)) {
            stop(paste0("====> Error: qc_dir not found: ", qc_dir))
        }
        qc_dir <- normalizePath(qc_dir, mustWork = TRUE)

        if (!file.exists(paste0(qc_dir, "/sample_qc_cutoffs.tsv"))) {
            stop(paste0("====> Error: sample_qc_cutoffs.tsv is not in ", qc_dir, ". Please use qcout_samqc_cutoffs to create it."))
        }

        # The "Run Experiment QC" section (2.3) below reads experiment_qc_*
        # files, which main() writes to their own experiment_qc/ directory,
        # a sibling of screen_qc/ and plasmid_qc/ -- NOT inside qc_dir itself.
        # If the caller doesn't tell us where that is, guess the sibling
        # folder next to qc_dir; fall back to qc_dir itself (old behaviour)
        # if that guess doesn't exist.
        if (is.null(experiment_qc_dir)) {
            guessed_dir <- file.path(dirname(qc_dir), "experiment_qc")
            experiment_qc_dir <- if (dir.exists(guessed_dir)) guessed_dir else qc_dir
        }
        if (!dir.exists(experiment_qc_dir)) {
            stop(paste0("====> Error: experiment_qc_dir not found: ", experiment_qc_dir))
        }
        experiment_qc_dir <- normalizePath(experiment_qc_dir, mustWork = TRUE)

        qc_type <- match.arg(qc_type)

        #------------------#
        # creating reports #
        #------------------#
        package_version <- if (requireNamespace("MAVEQC", quietly = TRUE)) {
            paste0("MAVEQC", "-v", packageVersion("MAVEQC"))
        } else {
            "MAVEQC (standalone script, unversioned)"
        }
        report_path <- paste0(qc_dir, "/", "MAVEQC_report.Rmd")
        sink(report_path)

        cat("---", "\n", sep = "")
        cat("title: \"MAVE QC Report\"", "\n", sep = "")
        cat("author: \"", package_version, "\"", "\n", sep = "")
        cat("date: \"`r Sys.time()`\"", "\n", sep = "")
        cat("output:", "\n", sep = "")
        cat("    html_document:", "\n", sep = "")
        cat("        toc: true", "\n", sep = "")
        cat("        toc_depth: 4", "\n", sep = "")
        cat("        theme: united", "\n", sep = "")
        cat("        highlight: tango", "\n", sep = "")
        cat("---", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{r setup, include = FALSE}", "\n", sep = "")
        cat("knitr::opts_chunk$set(echo = TRUE, fig.align = \"center\")", "\n", sep = "")
        cat("library(reactable)", "\n", sep = "")
        cat("outdir <- \"", qc_dir, "\"", "\n", sep = "")
        cat("outdir_exp <- \"", experiment_qc_dir, "\"", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{js, echo = FALSE}", "\n", sep = "")
        cat("function formatNumber(num, precision = 1) {", "\n", sep = "")
        cat("    const map = [", "\n", sep = "")
        cat("        { suffix: 'T', threshold: 1e12 },", "\n", sep = "")
        cat("        { suffix: 'B', threshold: 1e9 },", "\n", sep = "")
        cat("        { suffix: 'M', threshold: 1e6 },", "\n", sep = "")
        cat("        { suffix: 'K', threshold: 1e3 },", "\n", sep = "")
        cat("        { suffix: '', threshold: 1 },", "\n", sep = "")
        cat("    ];", "\n", sep = "")
        cat("    const found = map.find((x) => Math.abs(num) >= x.threshold);", "\n", sep = "")
        cat("    if (found) {", "\n", sep = "")
        cat("        const formatted = (num / found.threshold).toFixed(precision) + found.suffix;", "\n", sep = "")
        cat("        return formatted;", "\n", sep = "")
        cat("    }", "\n", sep = "")
        cat("    return num;", "\n", sep = "")
        cat("}", "\n", sep = "")
        cat("", "\n", sep = "")
        cat("function rangeMore(column, state) {", "\n", sep = "")
        cat("    let min = Infinity", "\n", sep = "")
        cat("    let max = 0", "\n", sep = "")
        cat("    state.data.forEach(function(row) {", "\n", sep = "")
        cat("        const value = row[column.id]", "\n", sep = "")
        cat("        if (value < min) {", "\n", sep = "")
        cat("            min = Math.floor(value)", "\n", sep = "")
        cat("        }", "\n", sep = "")
        cat("        if (value > max) {", "\n", sep = "")
        cat("            max = Math.ceil(value)", "\n", sep = "")
        cat("        }", "\n", sep = "")
        cat("    })", "\n", sep = "")
        cat("", "\n", sep = "")
        cat("    const filterValue = column.filterValue || min", "\n", sep = "")
        cat("    const input = React.createElement('input', {", "\n", sep = "")
        cat("        type: 'range',", "\n", sep = "")
        cat("        value: filterValue,", "\n", sep = "")
        cat("        min: min,", "\n", sep = "")
        cat("        max: max,", "\n", sep = "")
        cat("        onChange: function(event) {", "\n", sep = "")
        cat("            column.setFilter(event.target.value || undefined)", "\n", sep = "")
        cat("        },", "\n", sep = "")
        cat("        style: { width: '100%', marginRight: '8px' },", "\n", sep = "")
        cat("        'aria-label': 'Filter ' + column.name", "\n", sep = "")
        cat("    })", "\n", sep = "")
        cat("", "\n", sep = "")
        cat("    return React.createElement(", "\n", sep = "")
        cat("        'div',", "\n", sep = "")
        cat("        { style: { display: 'flex', alignItems: 'center', height: '100%' } },", "\n", sep = "")
        cat("        [input, formatNumber(filterValue)]", "\n", sep = "")
        cat("    )", "\n", sep = "")
        cat("}", "\n", sep = "")
        cat("", "\n", sep = "")
        cat("function filterMinValue(rows, columnId, filterValue) {", "\n", sep = "")
        cat("    return rows.filter(function(row) {", "\n", sep = "")
        cat("        return row.values[columnId] >= filterValue", "\n", sep = "")
        cat("    })", "\n", sep = "")
        cat("}", "\n", sep = "")
        cat("", "\n", sep = "")
        cat("function rangeLess(column, state) {", "\n", sep = "")
        cat("    let min = Infinity", "\n", sep = "")
        cat("    let max = 0", "\n", sep = "")
        cat("    state.data.forEach(function(row) {", "\n", sep = "")
        cat("        const value = row[column.id]", "\n", sep = "")
        cat("        if (value < min) {", "\n", sep = "")
        cat("            min = Math.floor(value)", "\n", sep = "")
        cat("        }", "\n", sep = "")
        cat("        if (value > max) {", "\n", sep = "")
        cat("            max = Math.ceil(value)", "\n", sep = "")
        cat("        }", "\n", sep = "")
        cat("    })", "\n", sep = "")
        cat("", "\n", sep = "")
        cat("    const filterValue = column.filterValue || max", "\n", sep = "")
        cat("    const input = React.createElement('input', {", "\n", sep = "")
        cat("        type: 'range',", "\n", sep = "")
        cat("        value: filterValue,", "\n", sep = "")
        cat("        min: min,", "\n", sep = "")
        cat("        max: max,", "\n", sep = "")
        cat("        onChange: function(event) {", "\n", sep = "")
        cat("            column.setFilter(event.target.value || undefined)", "\n", sep = "")
        cat("        },", "\n", sep = "")
        cat("        style: { width: '100%', marginRight: '8px' },", "\n", sep = "")
        cat("        'aria-label': 'Filter ' + column.name", "\n", sep = "")
        cat("    })", "\n", sep = "")
        cat("", "\n", sep = "")
        cat("    return React.createElement(", "\n", sep = "")
        cat("        'div',", "\n", sep = "")
        cat("        { style: { display: 'flex', alignItems: 'center', height: '100%' } },", "\n", sep = "")
        cat("        [input, formatNumber(filterValue)]", "\n", sep = "")
        cat("    )", "\n", sep = "")
        cat("}", "\n", sep = "")
        cat("", "\n", sep = "")
        cat("function filterMaxValue(rows, columnId, filterValue) {", "\n", sep = "")
        cat("    return rows.filter(function(row) {", "\n", sep = "")
        cat("        return row.values[columnId] <= filterValue", "\n", sep = "")
        cat("    })", "\n", sep = "")
        cat("}", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("\n", sep = "")

        cat("---", "\n", sep = "")
        cat("\n", sep = "")

        cat("## 1. Introduction", "\n", sep = "")
        cat("MAVEQC is a flexible R-package that provides QC analysis of Saturation Genome Editing (SGE) experimental data. ",
            "Available under GPL 3.0 from https://github.com/wtsi-hgi/MAVEQC", "\n", sep = "")
        cat("\n", sep = "")

        cat("---", "\n", sep = "")
        cat("\n", sep = "")

        if (qc_type == "screen") {
            cat("## 2. Screen QC", "\n", sep = "")
        } else {
            cat("## 2. Plasmid QC", "\n", sep = "")
        }
        cat("Displays QC plots and statistics for all samples for QC.", "\n", sep = "")
        cat("\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("samqc_cutoffs <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_cutoffs.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("### 2.1. Sample Sheet", "\n", sep = "")
        cat("\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", samplesheet, "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = min_row, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("### 2.2. Run Sample QC", "\n", sep = "")
        cat("\n", sep = "")
        cat("#### 2.2.1. Read Length Distribution", "\n", sep = "")
        cat("Displays the percentage of reads for each sample, based on 50 nucleotide increments, using the total number of raw reads.", "\n", sep = "")
        cat("<p style=\"color:red\">Note: expected read length is 300.</p>", "\n", sep = "")
        cat("<p style=\"color:red\">Note: the lengths of primers are deducted from the read length based on the sample sheet information. ",
            "(see 2.1. Sample Sheet: quants_append_start and quants_append_end)</p>", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{r, echo = FALSE, out.height = \"50%\", out.width = \"50%\"}", "\n", sep = "")
        cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_read_length.png\"), rel_path = FALSE)", "\n", sep = "")
        cat("```", "\n", sep = "")

        cat("**Pass criterion:** more than 90% of reads are longer than 200 nucleotides", "\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_read_length.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = min_row, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"% 0 ~ 50\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                               style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"% 50 ~ 100\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                 style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"% 100 ~ 150\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                  style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"% 150 ~ 200\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                  style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"% 200 ~ 250\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                  style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"% 250 ~ 300\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                  style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"Total Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                  format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })),", "\n", sep = "")
        cat("          rowStyle = function(index) { if (!(df[index, \"Pass\"])) { list(background = t_col(\"tomato\", 0.2)) }})", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("#### 2.2.2. Missing Variants", "\n", sep = "")
        cat("Stats of missing variants in the library", "\n", sep = "")
        cat("\n", sep = "")
        cat("**Pass criterion:** less than 1% of expected variants are missing", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_stats_missing.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = min_row, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"Library Sequences\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Missing Library Sequences\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"% Missing Library Sequences\" = colDef(style = function(value) { ", "\n", sep = "")
        cat("                                                                                          if (value > samqc_cutoffs$per_missing_variants * 100) {", "\n", sep = "")
        cat("                                                                                              color <- \"red\"", "\n", sep = "")
        cat("                                                                                              fweight <- \"bold\"", "\n", sep = "")
        cat("                                                                                          } else {", "\n", sep = "")
        cat("                                                                                              color <- \"forestgreen\"", "\n", sep = "")
        cat("                                                                                              fweight <- \"bold\" }", "\n", sep = "")
        cat("                                                                                          list(color = color, fontWeight = fweight)}),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })),", "\n", sep = "")
        cat("          rowStyle = function(index) { if (!(df[index, \"Pass\"])) { list(background = t_col(\"tomato\", 0.2)) }})", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("Records of missing variants in the library", "\n", sep = "")
        cat("<p style=\"color:red\">Note: Unique indicates that a template sequence occurs only once in the VaLiAnT meta file. (1: Unique, 0: Not Unique)",
            "This is important as a template sequence can occur more than once depending on the mutation types applied in VaLiAnT.</p>", "\n", sep = "")
        cat("<p style=\"color:red\">Note: Table below shows all the missing variants in all the samples, so the variants may occur multiple times.</p>", "\n", sep = "")

        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qc_dir, "/missing_variants_in_library.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("df <- df[, c(\"name\", \"length\", \"count\", \"unique\", \"sample\")]", "\n", sep = "")
        cat("colnames(df) <- capital_names(colnames(df))", "\n", sep = "")
        cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = min_row, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"Name\" = colDef(minWidth = 300),", "\n", sep = "")
        cat("                         \"Sample\" = colDef(minWidth = 300)))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("#### 2.2.3. Total Reads (Counts)", "\n", sep = "")
        cat("Displays the total number of reads per sample. ",
            "Filtering based on 1-dimensional Kmean clustering that excludes unique sequences with low read counts.", "\n", sep = "")
        cat("\n", sep = "")
        cat("* **Accepted reads:** Total read count for all unique sequences with sufficient reads based 1D Kmean clustering.", "\n", sep = "")
        cat("* **Excluded reads:** Total read count for all unique sequences with insufficient reads based 1D Kmean clustering.", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{r, echo = FALSE, out.height = \"50%\", out.width = \"50%\"}", "\n", sep = "")
        cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_stats_total.png\"), rel_path = FALSE)", "\n", sep = "")
        cat("```", "\n", sep = "")

        cat("Total Reads: the total number of raw reads", "\n", sep = "")
        cat("\n", sep = "")
        cat("**Pass criterion:** more than 1,000,000 total reads", "\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_stats_total.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = min_row, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"Accepted Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                     format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"% Accepted Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                       style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"Excluded Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                     format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"% Excluded Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                       style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"Total Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                  format = colFormat(separators = TRUE),", "\n", sep = "")
        cat("                                                  style = function(value) { ", "\n", sep = "")
        cat("                                                              if (value < samqc_cutoffs$num_total_reads) {", "\n", sep = "")
        cat("                                                                  color <- \"red\"", "\n", sep = "")
        cat("                                                                  fweight <- \"bold\"", "\n", sep = "")
        cat("                                                              } else {", "\n", sep = "")
        cat("                                                                  color <- \"forestgreen\"", "\n", sep = "")
        cat("                                                                  fweight <- \"bold\" }", "\n", sep = "")
        cat("                                                              list(color = color, fontWeight = fweight)}),", "\n", sep = "")
        cat("                         \"Pass Threshold\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })),", "\n", sep = "")
        cat("          rowStyle = function(index) { if (!(df[index, \"Pass\"])) { list(background = t_col(\"tomato\", 0.2)) }})", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("#### 2.2.4. Accepted Reads (Percentage)", "\n", sep = "")
        cat("Displays the percentage of library reads vs non-library reads (ie. Reference, PAM and Unmapped) for Accepted Reads (see 2.2.3 explanation).", "\n", sep = "")
        cat("\n", sep = "")
        cat("* **Library Reads:** Percentage reads mapping to template oligo sequences, including intended variants.", "\n", sep = "")
        cat("* **Reference Reads:** Percentage reads mapping to Reference.", "\n", sep = "")
        cat("* **PAM Reads:** Percentage reads mapping to PAM/Protospacer Protection Edits (PPEs), without intended variant.", "\n", sep = "")
        cat("* **Unmapped Reads:** Percentage of Unmapped Reads (not mapped to library sequences, PAM sequence, and reference sequence).", "\n", sep = "")
        cat("* **Library Coverage:** Mean read count per template oligo sequence.", "\n", sep = "")
        cat("\n", sep = "")

        cat("\n", sep = "")
        df <- as.data.frame(read.table(paste0(qc_dir, "/sample_qc_stats_accepted.tsv"), header = TRUE, sep = "\t", check.names = FALSE))
        samples_ref0 <- df[df[, 4] == 0, 2]
        samples_pam0 <- df[df[, 5] == 0, 2]
        if (length(samples_ref0) != 0) {
            cat("<p style=\"color:red; font-weight: bold\">Warning: found samples with 0 counts of resequence sequence! Please check table below.</p>", "\n", sep = "")
        }
        if (length(samples_pam0) != 0) {
            cat("<p style=\"color:red; font-weight: bold\">Warning: found samples with 0 counts of pam sequence! Please check table below.</p>", "\n", sep = "")
        }
        cat("\n", sep = "")

        cat("```{r, echo = FALSE, out.height = \"50%\", out.width = \"50%\"}", "\n", sep = "")
        cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_stats_accepted.png\"), rel_path = FALSE)", "\n", sep = "")
        cat("```", "\n", sep = "")

        cat("**Pass criterion:** more than 40% of accepted reads are library reads", "\n", sep = "")
        cat("\n", sep = "")
        cat("<p style=\"color:red\">Note: Accepted reads are the filtered reads based on 2.2.3</p>", "\n", sep = "")

        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_stats_accepted.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = min_row, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"% Library Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                      style = function(value) {", "\n", sep = "")
        cat("                                                                  if (value < samqc_cutoffs$per_library_reads * 100) {", "\n", sep = "")
        cat("                                                                      bar_style(length = value/100, color = \"red\", fweight = \"bold\")", "\n", sep = "")
        cat("                                                                  } else {", "\n", sep = "")
        cat("                                                                      bar_style(length = value/100, color = \"forestgreen\", fweight = \"bold\") }}),", "\n", sep = "")
        cat("                         \"% Reference Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                        style = function(value) {", "\n", sep = "")
        cat("                                                                    if (value == 0) {", "\n", sep = "")
        cat("                                                                        bar_style(length = value/100, color = \"red\", fweight = \"bold\")", "\n", sep = "")
        cat("                                                                    } else {", "\n", sep = "")
        cat("                                                                        bar_style(length = value/100)}}),", "\n", sep = "")
        cat("                         \"% PAM Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                  style = function(value) {", "\n", sep = "")
        cat("                                                                    if (value == 0) {", "\n", sep = "")
        cat("                                                                        bar_style(length = value/100, color = \"red\", fweight = \"bold\")", "\n", sep = "")
        cat("                                                                    } else {", "\n", sep = "")
        cat("                                                                        bar_style(length = value/100)}}),", "\n", sep = "")
        cat("                         \"% Unmapped Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                       style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })),", "\n", sep = "")
        cat("          rowStyle = function(index) { if (!(df[index, \"Pass\"])) { list(background = t_col(\"tomato\", 0.2)) }})", "\n", sep = "")
        cat("```", "\n", sep = "")

        cat("\n", sep = "")

        cat("Defines the mean read count per template oligo sequence (dividing the total number of library reads by the total number of library sequences).", "\n", sep = "")
        cat("\n", sep = "")
        cat("**Pass criterion:** library coverage is more than 100 reads", "\n", sep = "")

        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_stats_coverage.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = min_row, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"Total Library Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                          format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Total Library Sequences\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                                     format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Library Coverage\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                       format = colFormat(separators = TRUE),", "\n", sep = "")
        cat("                                                       style = function(value) {", "\n", sep = "")
        cat("                                                                   if (value < samqc_cutoffs$library_cov) {", "\n", sep = "")
        cat("                                                                       color <- \"red\"", "\n", sep = "")
        cat("                                                                       fweight <- \"bold\"", "\n", sep = "")
        cat("                                                                   } else {", "\n", sep = "")
        cat("                                                                       color <- \"forestgreen\"", "\n", sep = "")
        cat("                                                                       fweight <- \"bold\" }", "\n", sep = "")
        cat("                                                                   list(color = color, fontWeight = fweight)}),", "\n", sep = "")
        cat("                         \"Median Coverage\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeMore\"),", "\n", sep = "")
        cat("                                                      format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })),", "\n", sep = "")
        cat("          rowStyle = function(index) { if (!(df[index, \"Pass\"])) { list(background = t_col(\"tomato\", 0.2)) }})", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("#### 2.2.5. Genomic Coverage", "\n", sep = "")
        cat("Distribution of variants across targeton region based on log2(count+1) values.", "\n", sep = "")
        if (qc_type == "screen") {
            cat("<p style=\"color:red\">Note: Does not show missing varaints (0 count in the libary).</p>", "\n", sep = "")
        }
        cat("\n", sep = "")

        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_position_cov.dots.png\"), rel_path = FALSE)", "\n", sep = "")
        cat("```", "\n", sep = "")

        cat("Low Abundance cutoff: the green dashed line indicates the threshold which is used to determine if the variant is low abundance (less than 5 reads)", "\n", sep = "")
        cat("\n", sep = "")
        cat("**Pass criterion:** the percentage of low-abundance variants is lower than 30%", "\n", sep = "")
        cat("\n", sep = "")
        cat("% Low Abundance: the percentage of variants below the low abundance cutoff", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_stats_pos_coverage.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("df_counts <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_stats_pos_counts.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("df_counts <- log2(df_counts + 1)", "\n", sep = "")
        cat("df$data <- NA", "\n", sep = "")
        cat("for (i in 1:nrow(df)) {", "\n", sep = "")
        cat("    tmp_data <- df_counts[, df[i, ]$Sample]", "\n", sep = "")
        cat("    df[i, ]$data <- list(tmp_data[!is.na(tmp_data)])", "\n", sep = "")
        cat("}", "\n", sep = "")
        cat("df$boxplot <- NA", "\n", sep = "")
        cat("boxplot_min <- min(df_counts, na.rm = TRUE)", "\n", sep = "")
        cat("boxplot_max <- max(df_counts, na.rm = TRUE)", "\n", sep = "")
        cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
        cat("reactable(df[, c(1, 2, 12, 7, 3:6, 8:10)], highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = min_row, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"Genomic Start\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Genomic End\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"boxplot\" = colDef(cell = function(value, index) {", "\n", sep = "")
        cat("                                                         if (length(df$data[[index]]) > 5) {", "\n", sep = "")
        cat("                                                             sparkline(df$data[[index]], type = \"box\", width = 120,", "\n", sep = "")
        cat("                                                                       chartRangeMin = boxplot_min, chartRangeMax = boxplot_max) }}),", "\n", sep = "")
        cat("                         \"% Low Abundance\" = colDef(style = function(value) {", "\n", sep = "")
        cat("                                                                  if (value > (1 - samqc_cutoffs$low_abundance_lib_per) * 100) {", "\n", sep = "")
        cat("                                                                      color <- \"red\"", "\n", sep = "")
        cat("                                                                      fweight <- \"bold\"", "\n", sep = "")
        cat("                                                                  } else {", "\n", sep = "")
        cat("                                                                      color <- \"forestgreen\"", "\n", sep = "")
        cat("                                                                      fweight <- \"bold\" }", "\n", sep = "")
        cat("                                                                  list(color = color, fontWeight = fweight)}),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })),", "\n", sep = "")
        cat("          rowStyle = function(index) { if (!(df[index, \"Pass\"])) { list(background = t_col(\"tomato\", 0.2)) }})", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        if (qc_type == "screen") {
            cat("#### 2.2.6. Genomic Position Percentage", "\n", sep = "")
            cat("Displays distribution of \"LOF\" (loss-of-function) vs all \"Other\" variants across the targeton region, ",
                "based on read percentages for reference timepoint. ",
                "Requires concordant distribution of LOF and Other variants.", "\n", sep = "")
            cat("<p style=\"color:red\">Note: Does not show missing varaints (0 count in the libary).</p>", "\n", sep = "")
            cat("\n", sep = "")

            cat("```{r, echo = FALSE, out.height = \"65%\", out.width = \"65%\"}", "\n", sep = "")
            cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_position_anno.lof_dots.png\"), rel_path = FALSE)", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("```{r, echo = FALSE}", "\n", sep = "")
            cat("df <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_stats_pos_percentage.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
            cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
            cat("reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
            cat("          filterable = TRUE, minRows = min_row, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
            cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
            cat("          columns = list(\"Genomic Start\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
            cat("                         \"Genomic End\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
            cat("                         \"% Low Abundance (ALL)\" = colDef(style = function(value) {", "\n", sep = "")
            cat("                                                                        if (value > (1 - samqc_cutoffs$low_abundance_lib_per) * 100) {", "\n", sep = "")
            cat("                                                                            color <- \"red\"", "\n", sep = "")
            cat("                                                                            fweight <- \"bold\"", "\n", sep = "")
            cat("                                                                        } else {", "\n", sep = "")
            cat("                                                                            color <- \"forestgreen\"", "\n", sep = "")
            cat("                                                                            fweight <- \"bold\" }", "\n", sep = "")
            cat("                                                                        list(color = color, fontWeight = fweight)}),", "\n", sep = "")
            cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })),", "\n", sep = "")
            cat("          rowStyle = function(index) { if (!(df[index, \"Pass\"])) { list(background = t_col(\"tomato\", 0.2)) }})", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")

            cat("### 2.3. Run Experiment QC", "\n", sep = "")
            cat("\n", sep = "")
            cat("#### 2.3.1. Sample Correlations", "\n", sep = "")
            cat("\n", sep = "")
            cat("```{r, echo = FALSE, out.height = \"50%\", out.width = \"50%\"}", "\n", sep = "")
            cat("knitr::include_graphics(paste0(outdir_exp, \"/experiment_qc_samples_tree.png\"), rel_path = FALSE)", "\n", sep = "")
            cat("knitr::include_graphics(paste0(outdir_exp, \"/experiment_qc_samples_corr.png\"), rel_path = FALSE)", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")
            cat("```{r, echo = FALSE}", "\n", sep = "")
            cat("df <- as.data.frame(read.table(\"", experiment_qc_dir, "/experiment_qc_corr.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
            cat("df$Pass <- NULL", "\n", sep = "")
            cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
            cat("reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
            cat("          filterable = TRUE, minRows = min_row, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
            cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
            cat("          columns = list(\"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })))", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")

            cat("#### 2.3.2. Sample PCA", "\n", sep = "")
            cat("\n", sep = "")
            cat("```{r, echo = FALSE, out.height = \"75%\", out.width = \"75%\"}", "\n", sep = "")
            cat("knitr::include_graphics(paste0(outdir_exp, \"/experiment_qc_pca_samples.png\"), rel_path = FALSE)", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")

            cat("#### 2.3.3. Fold Change (by category)", "\n", sep = "")
            cat("\n", sep = "")
            figs <- list.files(path = experiment_qc_dir, pattern = "experiment_qc_deseq_fc.*.all_beeswarm.png", full.names = TRUE)
            figs <- mixedsort(figs)
            tsvs <- list.files(path = experiment_qc_dir, pattern = "experiment_qc_deseq_fc.*.all.tsv", full.names = TRUE)
            tsvs <- mixedsort(tsvs)
            sum_tsvs <- list.files(path = experiment_qc_dir, pattern = "experiment_qc_deseq_fc.*.all_sum.tsv", full.names = TRUE)
            sum_tsvs <- mixedsort(sum_tsvs)
            if (length(figs) == 0) {
                warning(paste0("No experiment_qc_deseq_fc.*.all_beeswarm.png files found in ", experiment_qc_dir, "; skipping section 2.3.3."))
            }
            for (i in seq_along(figs)) {
                tmp_header <- strsplit(tail(strsplit(figs[i], "/", fixed = TRUE)[[1]], n = 1), ".", fixed = TRUE)[[1]][2]
                cat("**", tmp_header, "**", "\n", sep = "")

                cat("```{r, echo = FALSE, out.height = \"75%\", out.width = \"75%\"}", "\n", sep = "")
                cat("knitr::include_graphics(\"", figs[i], "\", rel_path = FALSE)", "\n", sep = "")
                cat("```", "\n", sep = "")
                cat("\n", sep = "")

                cat("```{r, echo = FALSE}", "\n", sep = "")
                cat("df <- as.data.frame(read.table(\"", tsvs[i], "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
                cat("df_sum <- as.data.frame(read.table(\"", sum_tsvs[i], "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
                cat("df_sum$data <- NA", "\n", sep = "")
                cat("for (i in 1:nrow(df_sum)) {", "\n", sep = "")
                cat("    df_sum[i, ]$data <- list(sort(df[df$consequence == df_sum[i, ]$consequence, ]$pos_adj_log2FoldChange_raw))", "\n", sep = "")
                cat("}", "\n", sep = "")
                cat("df_sum$boxplot <- NA", "\n", sep = "")
                cat("df_sum$barplot <- NA", "\n", sep = "")
                cat("min_row <- ifelse(nrow(df_sum) > 10, 10, nrow(df_sum))", "\n", sep = "")
                cat("boxplot_min <- min(df$pos_adj_log2FoldChange_raw, na.rm = TRUE)", "\n", sep = "")
                cat("boxplot_max <- max(df$pos_adj_log2FoldChange_raw, na.rm = TRUE)", "\n", sep = "")
                cat("reactable(df_sum[, c(1, 9, 10, 2:7)], highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
                cat("          filterable = FALSE, sortable = FALSE, minRows = min_row, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
                cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
                cat("          columns = list(\"number of depleted\"  = colDef(style = function(value) { list(color = \"forestgreen\", fontWeight = \"bold\") }),", "\n", sep = "")
                cat("                         \"number of no impact\"  = colDef(style = function(value) { list(color = \"black\", fontWeight = \"bold\") }),", "\n", sep = "")
                cat("                         \"number of enriched\"  = colDef(style = function(value) { list(color = \"red\", fontWeight = \"bold\") }),", "\n", sep = "")
                cat("                         \"total number\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
                cat("                         \"total number of outside range\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
                cat("                         \"boxplot\" = colDef(cell = function(value, index) {", "\n", sep = "")
                cat("                                                         if (length(df_sum$data[[index]]) > 5) {", "\n", sep = "")
                cat("                                                             sparkline(df_sum$data[[index]], type = \"box\", width = 120,", "\n", sep = "")
                cat("                                                                       chartRangeMin = boxplot_min, chartRangeMax = boxplot_max) }}),", "\n", sep = "")
                cat("                         \"barplot\" = colDef(cell = function(value, index) {", "\n", sep = "")
                cat("                                                         if (length(df_sum$data[[index]]) > 5) {", "\n", sep = "")
                cat("                                                             sparkline(df_sum$data[[index]], type = \"tristate\", width = 120,", "\n", sep = "")
                cat("                                                                       posBarColor = \"red\", negBarColor = \"yellowgreen\") }}) ))", "\n", sep = "")
                cat("```", "\n", sep = "")
                cat("\n", sep = "")
            }
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")

            cat("#### 2.3.4. Fold Change (by position)", "\n", sep = "")
            cat("\n", sep = "")
            figs <- list.files(path = experiment_qc_dir, pattern = "experiment_qc_deseq_fc.*.all_position.png", full.names = TRUE)
            figs <- mixedsort(figs)
            tsvs <- list.files(path = experiment_qc_dir, pattern = "experiment_qc_deseq_fc.*.all.tsv", full.names = TRUE)
            tsvs <- mixedsort(tsvs)
            if (length(figs) == 0) {
                warning(paste0("No experiment_qc_deseq_fc.*.all_position.png files found in ", experiment_qc_dir, "; skipping section 2.3.4."))
            }
            for (i in seq_along(figs)) {
                tmp_header <- strsplit(tail(strsplit(figs[i], "/", fixed = TRUE)[[1]], n = 1), ".", fixed = TRUE)[[1]][2]
                cat("**", tmp_header, "**", "\n", sep = "")

                cat("```{r, echo = FALSE, out.height = \"75%\", out.width = \"75%\"}", "\n", sep = "")
                cat("knitr::include_graphics(\"", figs[i], "\", rel_path = FALSE)", "\n", sep = "")
                cat("```", "\n", sep = "")
                cat("\n", sep = "")

                cat("```{r, echo = FALSE}", "\n", sep = "")
                cat("df <- as.data.frame(read.table(\"", tsvs[i], "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
                cat("df <- df[df$stat_pos_raw != \"no impact\", c(\"oligo_name\", \"consequence\", \"position\", \"pos_adj_log2FoldChange_raw\", \"pos_adj_fdr_raw\", \"stat_pos_raw\")]", "\n", sep = "")
                cat("max_log2fc <- ifelse(nrow(df) == 0, 0, max(abs(df$pos_adj_log2FoldChange_raw)))", "\n", sep = "")
                cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
                cat("reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
                cat("          filterable = TRUE, minRows = min_row, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
                cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
                cat("          columns = list(\"oligo_name\" = colDef(minWidth = 300),", "\n", sep = "")
                cat("                         \"pos_adj_log2FoldChange_raw\" = colDef(filterMethod = JS(\"filterMaxValue\"),", "\n", sep = "")
                cat("                                                     filterInput = JS(\"rangeLess\"),", "\n", sep = "")
                cat("                                                     format = colFormat(digits = 2),", "\n", sep = "")
                cat("                                                     cell = function(value) {", "\n", sep = "")
                cat("                                                                label <- round(value, 2)", "\n", sep = "")
                cat("                                                                bar_chart_pos_neg(label, value, max_value = max_log2fc)},", "\n", sep = "")
                cat("                                                     align = \"center\",", "\n", sep = "")
                cat("                                                     minWidth = 300),", "\n", sep = "")
                cat("                         \"pos_adj_fdr_raw\" = colDef(format = colFormat(digits = 3)),", "\n", sep = "")
                cat("                         \"stat_pos_raw\"  = colDef(style = function(value) {", "\n", sep = "")
                cat("                                                        if (value == \"enriched\") {", "\n", sep = "")
                cat("                                                            color <- \"red\"", "\n", sep = "")
                cat("                                                            fweight <- \"bold\"", "\n", sep = "")
                cat("                                                        } else {", "\n", sep = "")
                cat("                                                            color <- \"forestgreen\"", "\n", sep = "")
                cat("                                                            fweight <- \"bold\" }", "\n", sep = "")
                cat("                                                        list(color = color, fontWeight = fweight)})))", "\n", sep = "")
                cat("```", "\n", sep = "")
                cat("\n", sep = "")
            }
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")
        }

        cat("---", "\n", sep = "")
        cat("\n", sep = "")

        cat("## 3. QC Results", "\n", sep = "")
        cat("Summarising the final results, below are the cutoffs using for PASS/FAIL", "\n", sep = "")
        cat("\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(matrix(1:24, nrow = 8))", "\n", sep = "")
        cat("colnames(df) <- c(\"QC Tag\", \"Cutoff\", \"Description\")", "\n", sep = "")
        cat("df[, 1] <- c(\"Gini coefficient\",", "\n", sep = "")
        cat("             \"Total Reads\",", "\n", sep = "")
        cat("             \"% Missing Variants\",", "\n", sep = "")
        cat("             \"Accepted Reads\",", "\n", sep = "")
        cat("             \"% Mapping Reads\",", "\n", sep = "")
        cat("             \"% Reference Reads\",", "\n", sep = "")
        cat("             \"% Library Reads\",", "\n", sep = "")
        cat("             \"Library Coverage\")", "\n", sep = "")
        cat("df[, 2] <- c(maveqc_config$gini_coeff,", "\n", sep = "")
        cat("             samqc_cutoffs$num_total_reads,", "\n", sep = "")
        cat("             samqc_cutoffs$per_missing_variants,", "\n", sep = "")
        cat("             samqc_cutoffs$num_accepted_reads,", "\n", sep = "")
        cat("             samqc_cutoffs$per_mapping_reads,", "\n", sep = "")
        cat("             samqc_cutoffs$per_ref_reads,", "\n", sep = "")
        cat("             samqc_cutoffs$per_library_reads,", "\n", sep = "")
        cat("             samqc_cutoffs$library_cov)", "\n", sep = "")
        cat("df[, 3] <- c(\"Gini coefficient must be lower than [cutoff]\",", "\n", sep = "")
        cat("             \"Sample must have more than [cutoff] total reads\",", "\n", sep = "")
        cat("             \"Missing varaints in the library must be less than [cutoff]\",", "\n", sep = "")
        cat("             \"Sample must have more than [cutoff] reads after the low count filtering\",", "\n", sep = "")
        cat("             \"Sample must have more than [cutoff] of reads aligned to the library including reference and PAM reads\",", "\n", sep = "")
        cat("             \"Sample must have less than [cutoff] of reads aligned to reference sequence\",", "\n", sep = "")
        cat("             \"Sample must have more than [cutoff] of reads aligned to the library\",", "\n", sep = "")
        cat("             \"Sample must have more than [cutoff] average coverage\")", "\n", sep = "")
        cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          minRows = min_row, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"QC Tag\" = colDef(minWidth = 200),", "\n", sep = "")
        cat("                         \"Description\" = colDef(minWidth = 800)))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("### 3.1. Sample QC Results", "\n", sep = "")
        cat("\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_results.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("col_pal <- function(x) rgb(colorRamp(c(\"seagreen\", \"limegreen\", \"yellow3\", \"orange\", \"red\"))(x), maxColorValue = 255)", "\n", sep = "")
        cat("gini_min <- 0", "\n", sep = "")
        cat("gini_max <- 1", "\n", sep = "")
        cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = min_row, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"Gini Coefficient\" = colDef(style = function(value) {", "\n", sep = "")
        cat("                                                                   normalized <- (value - gini_min) / (gini_max - gini_min)", "\n", sep = "")
        cat("                                                                   list(color = col_pal(normalized))}),", "\n", sep = "")
        cat("                         \"Gini Coefficient Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" }),", "\n", sep = "")
        cat("                         \"Total Reads\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" }),", "\n", sep = "")
        cat("                         \"% Missing Variants\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" }),", "\n", sep = "")
        cat("                         \"Accepted Reads\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" }),", "\n", sep = "")
        cat("                         \"% Mapping Reads\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" }),", "\n", sep = "")
        cat("                         \"% Reference Reads\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" }),", "\n", sep = "")
        cat("                         \"% Library Reads\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" }),", "\n", sep = "")
        cat("                         \"Library Coverage\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("## 4. Methods and Glossary", "\n", sep = "")
        cat("### 4.1. Methods", "\n", sep = "")
        cat("#### 4.1.1. Methods of generating accepted reads (refer to 2.2.3)", "\n", sep = "")
        cat("1) Apply 1D Kmean clustering on each sequence (variant sequence) using log2 read count. Low read count sequences are removed in this step", "\n", sep = "")
        cat("2) A valid sequence must have at least 5 count in at least 25% of the samples in an experiment.", "\n", sep = "")
        cat("\n", sep = "")

        if (qc_type == "screen") {
            cat("#### 4.1.2. Methods of DESeq2 calculation (refer to 2.3.3 and 2.3.4)", "\n", sep = "")
            cat("1) using the total number of accepted reads to calculate the size factor which is applied in DESeq2 normalisation", "\n", sep = "")
            cat("2) run DESeq2 for each consequence", "\n", sep = "")
            cat("3) select synonymous variants and intronic variants as the control, then calculate the median log2 fold change of the control variants", "\n", sep = "")
            cat("4) re-calculate log2 fold change of other consequences by deducting the median log2 fold change of the control variants", "\n", sep = "")
            cat("5) re-calculate p value and adjusted p value", "\n", sep = "")
            cat("\n", sep = "")

            cat("### 4.2. Glossary", "\n", sep = "")
            cat("#### 4.2.1. Glossary of DESeq2 calculation (refer to 2.3.3 and 2.3.4)", "\n", sep = "")
            cat("| name | description |", "\n", sep = "")
            cat("| ----:|:---------------- |", "\n", sep = "")
            cat("| log2FoldChange | This is the initial log2FoldChange from DESeq2 using all the accepted reads |", "\n", sep = "")
            cat("| lfcSE | This is the initial lfcSE (log2FoldChange Standard Error) from DESeq2 using all the accepted reads |", "\n", sep = "")
            cat("| padj | This is the initial corrected p-value from DESeq2 using all the accepted reads |", "\n", sep = "")
            cat("| median control value | This is the median log2 fold change of the control variants (synonymous variants and intronic variants) |", "\n", sep = "")
            cat("| adj_log2FoldChange | This is the adjusted log2FoldChange calculated by deducting the median control value |", "\n", sep = "")
            cat("| adj_score | This is the adjusted score that is calculated from the adj_log2FoldChange divided by the lfcSE |", "\n", sep = "")
            cat("| adj_pval | This is the adjusted p-value derived from adj_score |", "\n", sep = "")
            cat("| adj_fdr | This is the adjusted FDR derived from adj_pval |", "\n", sep = "")
            cat("| stat | This indicates an enriched or a depleted status. adj_fdr < 0.05 & adj_log2FoldChange > 0 is enriched, adj_fdr < 0.05 & adj_log2FoldChange < 0 is depleted |", "\n", sep = "")
        }

        sink()

        rmarkdown::render(paste0(qc_dir, "/MAVEQC_report.Rmd"), clean = TRUE, quiet = TRUE)
        invisible(file.remove(paste0(qc_dir, "/MAVEQC_report.Rmd")))
}

#' initialize function
setGeneric("run_sample_qc", function(object, ...) {
  standardGeneric("run_sample_qc")
})

#' run sample QC for the list of samples
#'
#' @export
#' @name run_sample_qc
#' @param object                 sampleQC object
#' @param qc_type                plasmid or screen
#' @param cutoff_total           qc cutoff of the total reads
#' @param cutoff_missing_per     qc cutoff of the missing variant percentage
#' @param cutoff_low_count       count cutoff of library reads
#' @param cutoff_low_sample_per  sample percentage cutoff of library reads
#' @param cutoff_accepted        qc cutoff of the total accepted reads
#' @param cutoff_mapping_per     qc cutoff of mapping percentage (ref + pam + library)
#' @param cutoff_ref_per         qc cutoff of reference percentage
#' @param cutoff_library_per     qc cutoff of library reads percentage
#' @param cutoff_library_cov     qc cutoff of library coverage
#' @param cutoff_low_per         qc cutoff of low abundance percentage for LOF plot
#' @param cutoff_low_lib_per     qc cutoff of the percentage of library sequences with low abundance for LOF plot
#' @return object
setMethod(
    "run_sample_qc",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          cutoff_total = maveqc_config$sqc_total,
                          cutoff_missing_per = maveqc_config$sqc_missing,
                          cutoff_low_count = maveqc_config$sqc_low_count,
                          cutoff_low_sample_per = maveqc_config$sqc_low_sample_per,
                          cutoff_accepted = maveqc_config$sqc_accepted,
                          cutoff_mapping_per = maveqc_config$sqc_mapping_per,
                          cutoff_ref_per = maveqc_config$sqc_ref_per,
                          cutoff_library_per = maveqc_config$sqc_library_per,
                          cutoff_library_cov = maveqc_config$sqc_library_cov,
                          cutoff_low_per = maveqc_config$sqc_low_per,
                          cutoff_low_lib_per = maveqc_config$sqc_low_lib_per) {
        #----------#
        # checking #
        #----------#
        if (length(object@samples) == 0) {
            stop(paste0("====> Error: no sample found in the sampleQC object!"))
        }

        qc_type <- match.arg(qc_type)
        if (qc_type == "screen") {
            if (length(object@samples_ref) == 0) {
                stop(paste0("====> Error: samples_ref is empty! Screen QC must have reference samples."))
            }
        }

        cols <- c("num_total_reads",
                  "per_missing_variants",
                  "seq_low_count",
                  "seq_low_sample_per",
                  "num_accepted_reads",
                  "per_mapping_reads",
                  "per_ref_reads",
                  "per_library_reads",
                  "library_cov",
                  "low_abundance_per",
                  "low_abundance_lib_per")
        df_cutoffs <- data.frame(matrix(NA, 1, length(cols)))
        colnames(df_cutoffs) <- cols

        df_cutoffs$num_total_reads <- cutoff_total
        df_cutoffs$per_missing_variants <- cutoff_missing_per
        df_cutoffs$seq_low_count <- cutoff_low_count
        df_cutoffs$seq_low_sample_per <- cutoff_low_sample_per
        df_cutoffs$num_accepted_reads <- cutoff_accepted
        df_cutoffs$per_mapping_reads <- cutoff_mapping_per
        df_cutoffs$per_ref_reads <- cutoff_ref_per
        df_cutoffs$per_library_reads <- cutoff_library_per
        df_cutoffs$library_cov <- cutoff_library_cov
        df_cutoffs$low_abundance_per <- cutoff_low_per
        df_cutoffs$low_abundance_lib_per <- cutoff_low_lib_per

        object@cutoffs <- df_cutoffs

        #-------------------------------------------#
        # 1. Filtering by the total number of reads #
        #-------------------------------------------#
        cat("Filtering by the total number of reads...", "\n", sep = "")

        sample_names <- character()
        for (s in object@samples) {
            sample_names <- append(sample_names, s@sample)

            object@stats[s@sample, ]$per_r1_adaptor <- s@per_r1_adaptor
            object@stats[s@sample, ]$per_r2_adaptor <- s@per_r2_adaptor

            object@stats[s@sample, ]$total_reads <- s@allstats$total_counts
            object@stats[s@sample, ]$ref_reads <- s@allstats_qc$num_ref_reads
            object@stats[s@sample, ]$pam_reads <- s@allstats_qc$num_pam_reads

            object@stats[s@sample, ]$gini_coeff_before_qc <- s@libstats_qc$gini_coeff
        }

        #---------------------------------------#
        # 2. Filtering by low counts            #
        #    a) k-means clustering on screen QC #
        #    a) hard cutoff on Plasmid QC       #
        #---------------------------------------#
        cat("Filtering by low counts...", "\n", sep = "")

        if (qc_type == "screen") {
            cat("    |--> Creating k-means clusters...", "\n", sep = "")

            ref_counts <- data.table()
            for (s in object@samples_ref) {
                tmp_counts <- s@allcounts[, c("sequence", "count")]
                tmp_counts <- as.data.table(tmp_counts)

                if (nrow(ref_counts) == 0) {
                    ref_counts <- tmp_counts
                    colnames(ref_counts) <- c("sequence", s@sample)
                } else {
                    tmp_cols <- colnames(ref_counts)
                    ref_counts <- merge(ref_counts, tmp_counts, by = "sequence", all = TRUE)
                    colnames(ref_counts) <- c(tmp_cols, s@sample)
                }
            }
            ref_counts[, count := rowSums(.SD, na.rm = TRUE), .SDcols = 2:ncol(ref_counts)]
            ref_counts[, count_log2 := log2(count + 1)]

            kmeans_res <- Ckmeans.1d.dp(ref_counts$count_log2, k = 2, y = 1)
            ref_counts$cluster <- kmeans_res$cluster

            cat("    |--> Filtering using clusters...", "\n", sep = "")

            # filtering sequences on input samples by filtered set
            for (s in object@samples) {
                cat("        |--> Filtering on ", s@sample, "\n", sep = "")

                unfiltered_counts <- s@allcounts[, c("sequence", "count")]
                object@seq_clusters[[s@sample]] <- ref_counts

                # considering missing seqs
                object@accepted_counts[[s@sample]] <- unfiltered_counts[ref_counts[cluster == 2], on = .(sequence), nomatch = 0]

                object@bad_seqs_bycluster[[s@sample]] <- unfiltered_counts[!ref_counts[cluster == 2], on = .(sequence)]
            }
        } else {
            cat("    |--> Creating k-means clusters...", "\n", sep = "")

            for (s in object@samples) {
                cat("        |--> Filtering on ", s@sample, "\n", sep = "")

                tmp_counts <- s@allcounts[, c("sequence", "count")]
                tmp_counts <- as.data.table(tmp_counts)
                tmp_counts[, count_log2 := log2(count + 1)]

                kmeans_res <- Ckmeans.1d.dp(tmp_counts$count_log2, k = 2, y = 1)
                tmp_counts$cluster <- kmeans_res$cluster

                object@seq_clusters[[s@sample]] <- tmp_counts

                object@accepted_counts[[s@sample]] <- tmp_counts[cluster == 2]

                object@bad_seqs_bycluster[[s@sample]] <- tmp_counts[cluster == 1, c("sequence", "count")]
            }
        }

        #-------------------------------------#
        # 3. Filtering by depth and samples   #
        #    a) count >= X                    #
        #    b) in >= X% of samples           #
        #-------------------------------------#

        # if plasmid qc, don't apply percentage filtering as samples have different seqs
        if (qc_type == "screen") {
            cat("Filtering by depth and percentage in samples...", "\n", sep = "")

            # note library independent counts may have different sequences
            accepted_counts <- merge_list_to_dt(object@accepted_counts, "sequence", "count")
            accepted_counts[, sample_number := rowSums(.SD >= cutoff_low_count, na.rm = TRUE), .SDcols = 2:ncol(accepted_counts)]

            #accepted_counts[, sample_percentage := sample_number / length(sample_names)]

            # now use round the minimal number of samples, rather than percentage
            # like 0.25 * 9 = 2.25 minimal samples, 2.25 sounds not right, so use round(2.25) = 2
            min_samples <- round(cutoff_low_sample_per * length(sample_names))

            for (s in object@samples) {
                cols <- c("sequence", s@sample)
                #object@accepted_counts[[s@sample]] <- na.omit(accepted_counts[sample_percentage >= cutoff_low_sample_per, ..cols], cols = s@sample)
                object@accepted_counts[[s@sample]] <- na.omit(accepted_counts[sample_number >= min_samples, ..cols], cols = s@sample)
                colnames(object@accepted_counts[[s@sample]]) <- c("sequence", "count")

                cols <- c("sequence", s@sample)
                #object@bad_seqs_bydepth[[s@sample]] <- na.omit(accepted_counts[sample_percentage < cutoff_low_sample_per, ..cols], cols = s@sample)
                object@bad_seqs_bydepth[[s@sample]] <- na.omit(accepted_counts[sample_number < min_samples, ..cols], cols = s@sample)
                colnames(object@bad_seqs_bydepth[[s@sample]]) <- c("sequence", "count")
            }
        } else {
            cat("Filtering by depth...", "\n", sep = "")

            for (s in object@samples) {
                tmp_counts <- object@accepted_counts[[s@sample]]
                cols <- c("sequence", "count")

                object@accepted_counts[[s@sample]] <- tmp_counts[count >= cutoff_low_count, ..cols]

                object@bad_seqs_bydepth[[s@sample]] <- tmp_counts[count < cutoff_low_count, ..cols]
            }
        }

        #--------------------------------------#
        # 4. Filtering by library mapping      #
        #    a) reads mapped to VaLiAnT output #
        #--------------------------------------#
        cat("Filtering by library mapping...", "\n", sep = "")

        for (s in object@samples) {
            tmp_counts <- object@accepted_counts[[s@sample]]

            # meta_mseqs without ref and pam by format_count
            # and using library dependent sequences instead of meta sequences
            # accepted_counts have ref and pam
            object@library_counts[[s@sample]] <- object@accepted_counts[[s@sample]][sequence %in% s@meta_mseqs]
            object@unmapped_counts[[s@sample]] <- object@accepted_counts[[s@sample]][sequence %nin% c(s@meta_mseqs, s@refseq, s@pamseq)]

            object@bad_seqs_bylib[[s@sample]] <- object@unmapped_counts[[s@sample]]
        }

        #--------------------------------------#
        # 5. Filtering by library coverage     #
        #    a) library reads / oligos in meta #
        #--------------------------------------#
        cat("Filtering by library coverage...", "\n", sep = "")

        for (s in object@samples) {
            object@stats[s@sample, ]$accepted_reads <- sum(object@accepted_counts[[s@sample]]$count, na.rm = TRUE)
            object@stats[s@sample, ]$excluded_reads <- object@stats[s@sample, ]$total_reads - object@stats[s@sample, ]$accepted_reads
            object@stats[s@sample, ]$library_reads <- sum(object@library_counts[[s@sample]]$count, na.rm = TRUE)
            object@stats[s@sample, ]$unmapped_reads <- sum(object@unmapped_counts[[s@sample]]$count, na.rm = TRUE)

            object@stats[s@sample, ]$per_library_reads <- object@stats[s@sample, ]$library_reads / object@stats[s@sample, ]$accepted_reads
            object@stats[s@sample, ]$per_unmapped_reads <- object@stats[s@sample, ]$unmapped_reads / object@stats[s@sample, ]$accepted_reads
            object@stats[s@sample, ]$per_ref_reads <- object@stats[s@sample, ]$ref_reads / object@stats[s@sample, ]$accepted_reads
            object@stats[s@sample, ]$per_pam_reads <- object@stats[s@sample, ]$pam_reads / object@stats[s@sample, ]$accepted_reads

            object@stats[s@sample, ]$missing_meta_seqs <- length(s@missing_meta_seqs)
            object@stats[s@sample, ]$per_missing_meta_seqs <- length(s@missing_meta_seqs) / length(s@meta_mseqs)

            object@stats[s@sample, ]$library_seqs <- length(s@meta_mseqs)
            object@stats[s@sample, ]$median_cov <- median(object@library_counts[[s@sample]]$count)
            object@stats[s@sample, ]$library_cov <- as.integer(object@stats[s@sample, ]$library_reads / length(s@meta_mseqs))
        }

        #--------------------- -----------------#
        # 6. Sorting library counts by position #
        #---------------------------------------#
        cat("Sorting library counts by position...", "\n", sep = "")

        # main issue:
        # meta seqs have adaptors
        # oligo names are not unique
        # a seq has a consequence, but has many names, don't know which name is right

        for (s in object@samples) {
            cat("    |--> Sorting on ", s@sample, "\n", sep = "")

            tmp_meta <- s@valiant_meta[, c("oligo_name", "mut_position")]

            # fecth oligo name using library dependent counts
            tmp_map <- s@libcounts[, c("name", "sequence")]

            libcounts_pos <- object@library_counts[[s@sample]]
            libcounts_pos[tmp_map, oligo_name := i.name, on = .(sequence)]
            libcounts_pos[tmp_meta, position := i.mut_position, on = .(oligo_name)]
            setorder(libcounts_pos, cols = "position")

            # in case some sequences don't match in the library
            libcounts_pos <- na.omit(libcounts_pos, cols = "position")

            object@library_counts_pos[[s@sample]] <- libcounts_pos
            object@library_counts_chr[[s@sample]] <- c(unique(s@valiant_meta$ref_chr),
                                                       unique(s@valiant_meta$ref_strand),
                                                       unique(s@valiant_meta$ref_start),
                                                       unique(s@valiant_meta$ref_end))
        }

        #------------------------#
        # 7. Gini coeff after qc #
        #------------------------#
        cat("Calculating gini coefficiency...", "\n", sep = "")

        for (s in object@samples) {
            gini_coeff <- cal_gini(object@library_counts[[s@sample]]$count, corr = FALSE, na.rm = TRUE)
            object@stats[s@sample, ]$gini_coeff_after_qc <- round(gini_coeff, 3)
        }

        #------------------#
        # 8. QC results    #
        #------------------#
        object@stats$per_missing_meta_seqs <- unlist(lapply(object@stats$per_missing_meta_seqs, function(x) round(x, 4)))
        object@stats$per_library_reads <- unlist(lapply(object@stats$per_library_reads, function(x) round(x, 4)))
        object@stats$per_unmapped_reads <- unlist(lapply(object@stats$per_unmapped_reads, function(x) round(x, 4)))
        object@stats$per_ref_reads <- unlist(lapply(object@stats$per_ref_reads, function(x) round(x, 4)))
        object@stats$per_pam_reads <- unlist(lapply(object@stats$per_pam_reads, function(x) round(x, 4)))
        object@stats$per_r1_adaptor <- unlist(lapply(object@stats$per_r1_adaptor, function(x) round(x, 4)))
        object@stats$per_r2_adaptor <- unlist(lapply(object@stats$per_r2_adaptor, function(x) round(x, 4)))

        object@stats$qcpass_total_reads <- unlist(lapply(object@stats$total_reads, function(x) ifelse(x >= cutoff_total, TRUE, FALSE)))
        object@stats$qcpass_missing_per <- unlist(lapply(object@stats$per_missing_meta_seqs, function(x) ifelse(x < cutoff_missing_per, TRUE, FALSE)))
        object@stats$qcpass_accepted_reads <- unlist(lapply(object@stats$accepted_reads, function(x) ifelse(x >= cutoff_accepted, TRUE, FALSE)))
        object@stats$qcpass_mapping_per <- unlist(lapply(object@stats$per_unmapped_reads, function(x) ifelse(x < (1 - cutoff_mapping_per), TRUE, FALSE)))
        object@stats$qcpass_ref_per <- unlist(lapply(object@stats$per_ref_reads, function(x) ifelse(x < cutoff_ref_per, TRUE, FALSE)))
        object@stats$qcpass_library_per <- unlist(lapply(object@stats$per_library_reads, function(x) ifelse(x >= cutoff_library_per, TRUE, FALSE)))
        object@stats$qcpass_library_cov <- unlist(lapply(object@stats$library_cov, function(x) ifelse(x >= cutoff_library_cov, TRUE, FALSE)))

        qc_lables <- c("qcpass_total_reads",
                       "qcpass_missing_per",
                       "qcpass_accepted_reads",
                       "qcpass_mapping_per",
                       "qcpass_ref_per",
                       "qcpass_library_per",
                       "qcpass_library_cov")
        object@stats$qcpass <- apply(object@stats[, qc_lables], 1, function(x) all(x))

        #------------------------#
        # 9. Filtered samples    #
        #------------------------#

        object@filtered_samples <- rownames(object@stats[object@stats$qcpass == TRUE, ])

        #-------------------------#
        # 10. map vep consequence #
        #-------------------------#

        # main issue:
        # oligo names are not unique to sequence, cannot use as reference
        # a seq has a consequence, but has many names, don't which name is right
        # so, vep must have right seq, otherwise cannot determine the consequence

        # if plasmid qc, don't apply
        if (qc_type == "screen") {
            cat("Mapping consequencing annotation...", "\n", sep = "")

            # assuming all the samples have the sample library sequences and corresponding consequences
            # using sequences in vep annotation to identify consequences
            # assuming unique_oligo_name in vep_anno is unique
            vep_anno <- object@samples[[1]]@vep_anno[, c("unique_oligo_name", "seq", "summary_plot")]
            colnames(vep_anno) <- c("oligo_name", "sequence", "consequence")

            # merge counts, but use data table in case rownames of data frame is not unique
            library_counts_anno <- merge_list_to_dt(object@library_counts, "sequence", "count")
            cols <- colnames(library_counts_anno)
            library_counts_anno[vep_anno, c("oligo_name", "consequence") := .(oligo_name, consequence), on = .(sequence)]
            # Remove "oligo_name" and "consequence" from cols
            cols <- cols[!cols %in% c("oligo_name", "consequence")]

            # Now create the order_vector
            order_vector <- c("oligo_name", "consequence", cols)
            setcolorder(library_counts_anno, order_vector)
            object@library_counts_anno <- library_counts_anno

            # merge all the sorted library counts
            library_counts_pos_anno <- data.table()
            for (s in object@samples) {
                tmp_pos <- object@library_counts_pos[[s@sample]][, c("sequence", "position", "count")]
                colnames(tmp_pos) <- c("sequence", "position", s@sample)

                if (nrow(library_counts_pos_anno) == 0) {
                    library_counts_pos_anno <- tmp_pos
                } else {
                    library_counts_pos_anno <- merge(library_counts_pos_anno, tmp_pos, by = c("sequence", "position"), all = TRUE)
                }
            }

            cols <- colnames(library_counts_pos_anno)[-2]
            library_counts_pos_anno[vep_anno, c("oligo_name", "consequence") := .(oligo_name, consequence), on = .(sequence)]
            setcolorder(library_counts_pos_anno, c("oligo_name", "consequence", "position", cols))
            object@library_counts_pos_anno <- library_counts_pos_anno
        }

        return(object)
    }
)

#'
#' @param object experimentQC object
#' @param output_dir directory to save files
#' @return writes tsv files with complete results
extract_deseq2_results <- function(object, output_dir) {
    
    # Process library analysis results
	if (!length(object@lib_deseq_res_anno)) {
		warning("object@lib_deseq_res_anno is empty; did run_experiment_qc_lib_lfc() run?")
	  }

	  for (comparison in names(object@lib_deseq_res_anno)) {
		lib_res <- as.data.frame(object@lib_deseq_res_anno[[comparison]])

		# desired column order for exp_lib
		desired_order_lib <- c(
		  "sequence", "oligo_name", "position", "consequence",
		  "baseMean_raw", "log2FoldChange_raw", "lfcSE_raw", "pvalue_raw", "padj_raw", "stat_raw",
		  "baseMean_shrunk", "log2FoldChange_shrunk", "lfcSE_shrunk", "pvalue_shrunk", "padj_shrunk", "stat_shrunk"
		)
		present_lib <- intersect(desired_order_lib, colnames(lib_res))
		lib_res <- lib_res[, c(present_lib, setdiff(colnames(lib_res), present_lib)), drop = FALSE]

		# write out
		write.table(
		  lib_res,
		  file       = file.path(output_dir, paste0("library_deseq2_results_", comparison, ".tsv")),
		  sep        = "\t",
		  quote      = FALSE,
		  row.names  = FALSE
		)
	  }
    
	# Process all counts analysis results
	if (!length(object@all_deseq_res_anno_adj)) {
		warning("object@all_deseq_res_anno_adj is empty; did run_experiment_qc_all_lfc() run?")
	  }

	  for (comparison in names(object@all_deseq_res_anno_adj)) {
		all_res <- as.data.frame(object@all_deseq_res_anno_adj[[comparison]])

		# desired column order for exp_all
		desired_order_all <- c(
		  "sequence", "oligo_name", "position", "consequence",
		  "baseMean_raw", "log2FoldChange_raw", "lfcSE_raw", "pvalue_raw", "padj_raw",
		  "baseMean_shrunk", "log2FoldChange_shrunk", "lfcSE_shrunk", "pvalue_shrunk", "padj_shrunk",

		  "adj_log2FoldChange_raw", "adj_score_raw", "adj_pval_raw", "adj_fdr_raw", "stat_adj_raw",
		  "adj_log2FoldChange_shrunk", "adj_score_shrunk", "adj_pval_shrunk", "adj_fdr_shrunk", "stat_adj_shrunk",

		  "pos_fit", "pos_fit_se", "pos_total_se_raw", "pos_total_se_shrunk",
		  "pos_adj_log2FoldChange_raw", "pos_adj_score_raw", "pos_adj_pval_raw", "pos_adj_fdr_raw", "stat_pos_raw",
		  "pos_adj_log2FoldChange_shrunk", "pos_adj_score_shrunk", "pos_adj_pval_shrunk", "pos_adj_fdr_shrunk", "stat_pos_shrunk"
		)
		
		present_all <- intersect(desired_order_all, colnames(all_res))
		all_res <- all_res[, c(present_all, setdiff(colnames(all_res), present_all)), drop = FALSE]

		write.table(
		  all_res,
		  file       = file.path(output_dir, paste0("all_deseq2_results_", comparison, ".tsv")),
		  sep        = "\t",
		  quote      = FALSE,
		  row.names  = FALSE
		)
	  }
	}

#' Extract and save normalized read counts
#' 
#' @param object experimentQC object
#' @param output_dir directory to save files
#' @return writes tsv files with normalized counts
extract_normalized_counts <- function(object, output_dir) {
    
    # Get normalized counts from both lib and all analyses
    
    # Library counts
    lib_norm_counts <- as.data.frame(object@lib_deseq_rlog)
    lib_norm_counts$sequence <- rownames(lib_norm_counts)
    
    # Add annotations
    lib_norm_counts <- merge(lib_norm_counts,
                           object@library_counts_anno[, c("sequence", "oligo_name", "consequence")],
                           by = "sequence", 
                           all.x = TRUE)
                           
    # Reorder columns to put metadata first
    col_order <- c("sequence", "oligo_name", "consequence", 
                   colnames(lib_norm_counts)[!colnames(lib_norm_counts) %in% 
                                           c("sequence", "oligo_name", "consequence")])
    lib_norm_counts <- lib_norm_counts[, col_order]
    
    # All counts (post filtering)
    all_norm_counts <- as.data.frame(object@all_deseq_rlog)
    all_norm_counts$sequence <- rownames(all_norm_counts)
    
    # Add annotations 
    all_norm_counts <- merge(all_norm_counts,
                           object@library_counts_anno[, c("sequence", "oligo_name", "consequence")],
                           by = "sequence",
                           all.x = TRUE)
                           
    # Reorder columns
    all_norm_counts <- all_norm_counts[, col_order]
    
    # Add position information if available
    if(!is.null(object@library_counts_pos_anno)) {
        pos_info <- object@library_counts_pos_anno[, c("sequence", "position")]
        
        lib_norm_counts <- merge(lib_norm_counts, pos_info,
                               by = "sequence",
                               all.x = TRUE)
                               
        all_norm_counts <- merge(all_norm_counts, pos_info,
                               by = "sequence", 
                               all.x = TRUE)
    }
    
    # Save to files
    write.table(lib_norm_counts,
                file = file.path(output_dir, "library_normalized_counts.tsv"),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
                
    write.table(all_norm_counts, 
                file = file.path(output_dir, "all_normalized_counts.tsv"),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
                
    # Return the data frames invisibly
    invisible(list(library = lib_norm_counts,
                  all = all_norm_counts))
}

# Main execution flow
main <- function(input_dir, output_dir) {
  # Create output directories
  screen_qc_dir <- file.path(output_dir, "screen_qc")
  plasmid_qc_dir <- file.path(output_dir, "plasmid_qc")
  experiment_qc_dir <- file.path(output_dir, "experiment_qc")
  
  dir.create(screen_qc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plasmid_qc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(experiment_qc_dir, recursive = TRUE, showWarnings = FALSE)

  # Import dat

  sge_objs <- import_sge_files(file.path(input_dir, "input"), "sample_sheet.tsv")

  # Create and load config
  config_path <- file.path(output_dir, "config.yaml")
  create_config(output_dir)
  maveqc_config <- read.config(file = config_path)

  # Run Plasmid QC
  samqc <- create_sampleqc_object(sge_objs)
  samqc <- run_sample_qc(samqc, "plasmid")
  qcplot_samqc_all(samqc, qc_type = "plasmid", plot_dir = plasmid_qc_dir)
  qcout_samqc_all(samqc, qc_type = "plasmid", out_dir = plasmid_qc_dir)

  # Run Screen QC
  samqc <- create_sampleqc_object(sge_objs)
  samqc <- run_sample_qc(samqc, "screen")
  qcplot_samqc_all(samqc, qc_type = "screen", plot_dir = screen_qc_dir)
  qcout_samqc_all(samqc, qc_type = "screen", out_dir = screen_qc_dir)

  # Run Experiment QC
  expqc <- create_experimentqc_object(samqc)
  expqc <- run_experiment_qc(expqc)

  qcplot_expqc_sample_corr(expqc, plot_dir = experiment_qc_dir)
  qcplot_expqc_sample_pca(expqc, plot_dir = experiment_qc_dir)

  # Write out the TSV tables that go alongside the plots above
  # (experiment_qc_corr.tsv, experiment_qc_deseq_fc.*.all.tsv/_sum.tsv) --
  # these are required by create_qc_reports() / MAVEQC_report.html below.
  qcout_expqc_all(expqc, out_dir = experiment_qc_dir)

  # DESeq plots using positional raw metrics
  qcplot_expqc_deseq_fc(
    expqc,
    eqc_type = "all",
    plot_type = "beeswarm",
    plot_dir = experiment_qc_dir,
    stat_col = "stat_pos_raw",
    lfc_col  = "pos_adj_log2FoldChange_raw"
  )
  
  qcplot_expqc_deseq_fc_pos(
    expqc,
    eqc_type = "all",
    plot_dir = experiment_qc_dir,
    stat_col = "stat_pos_raw",
    lfc_col  = "pos_adj_log2FoldChange_raw"
  )

  qcplot_expqc_positional_loess_diag(
    expqc,
    plot_dir  = experiment_qc_dir,
    targeton  = targeton_name,
    use       = "raw"
  )

  extract_deseq2_results(expqc, experiment_qc_dir)

  extract_normalized_counts(expqc, experiment_qc_dir)

  # Generate the MAVEQC_report.html for both plasmid and screen QC
  samplesheet_path <- file.path(input_dir, "input", "sample_sheet.tsv")

  create_qc_reports(samplesheet       = samplesheet_path,
                     qc_type           = "plasmid",
                     qc_dir            = plasmid_qc_dir,
                     experiment_qc_dir = experiment_qc_dir)

  create_qc_reports(samplesheet       = samplesheet_path,
                     qc_type           = "screen",
                     qc_dir            = screen_qc_dir,
                     experiment_qc_dir = experiment_qc_dir)
}

# Run the main function
# (only auto-runs the FULL pipeline if this file is Rscript'd directly with
# two CLI args, exactly like the original script; sourcing this file to reuse
# the function definitions -- e.g. create_qc_reports() -- will NOT trigger
# the full pipeline)
if (!is.null(input_dir) && !is.null(output_dir)) {
  main(input_dir, output_dir)
} else {
  message("run_maveqc_functions.R sourced: function definitions loaded, ",
          "full pipeline NOT run (no input_dir/output_dir supplied).")
}
