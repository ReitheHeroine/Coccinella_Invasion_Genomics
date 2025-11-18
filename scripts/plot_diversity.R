#!/usr/bin/env Rscript

# title: 'plot_diversity.R'
# project: 'BIOL624 Final Project: Selection Detection in Lady Beetles'
# author: 'Reina Hastings'
# contact: 'reinahastings13@gmail.com'
# date created: 11/17/2025
# last modified: 11/17/2025
# purpose: 'Plot Tajima’s D and nucleotide diversity (pi) for one or more chromosomes.'
# usage examples:
#   Rscript scripts/plot_diversity.R NC_058189.1.filtered
#   Rscript scripts/plot_diversity.R NC_058189.1.filtered NC_058190.1.filtered
#   Rscript scripts/plot_diversity.R all

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_diversity.R <results_dir> <chr1|all> [<chr2> ...]")
}

results_dir <- args[1]
chrom_args  <- args[-1]

# extract population or dataset label from results directory name
results_label <- basename(normalizePath(results_dir))

# clean labels like "results_Europe" → "Europe"
results_label <- sub("^results[_-]?", "", results_label)

if (length(args) == 0) {
  cat("Usage:\n",
      "  Rscript plot_diversity.R <chr1> [<chr2> ...]\n",
      "  Rscript plot_diversity.R all\n\n",
      "Examples:\n",
      "  Rscript plot_diversity.R NC_058189.1.filtered\n",
      "  Rscript plot_diversity.R NC_058189.1.filtered NC_058190.1.filtered\n",
      "  Rscript plot_diversity.R all\n")
  quit(status = 1)
}

# directories where diversity.sh wrote the outputs
taj_dir <- file.path(results_dir, "tajima")
pi_dir  <- file.path(results_dir, "pi")

# output directory for plots
out_dir <- file.path(results_dir, "figs", "diversity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# figure out which chromosomes to plot
if (length(chrom_args) == 1 && chrom_args[1] == "all")
 {
  # List all TajimaD files and infer chromosome labels by stripping suffix
  taj_files <- list.files(taj_dir, pattern = "\\.Tajima\\.D$", full.names = FALSE)
  if (length(taj_files) == 0) {
    stop("No TajimaD files found in ", taj_dir)
  }
  chroms <- sub("\\.Tajima\\.D$", "", taj_files)
  message("Detected chromosomes from TajimaD files: ",
          paste(chroms, collapse = ", "))
} else {
  chroms <- chrom_args
  message("Using chromosomes specified on command line: ",
          paste(chroms, collapse = ", "))
}

# helper to safely read a file or warn
read_div_file <- function(path, chr, kind) {
  if (!file.exists(path)) {
    warning("Missing ", kind, " file for ", chr, ": ", path)
    return(NULL)
  }
  read.table(path, header = TRUE)
}

# loop over chromosomes and make plots

# to build combined faceted plots, store all per-window data
all_taj_list <- list()
all_pi_list  <- list()

for (chr in chroms) {
  message("Processing ", chr, " ...")

  taj_path <- file.path(taj_dir, paste0(chr, ".Tajima.D"))
  pi_path  <- file.path(pi_dir,  paste0(chr, ".windowed.pi"))

  taj <- read_div_file(taj_path, chr, "TajimaD")
  pi  <- read_div_file(pi_path,  chr, "pi")

  if (is.null(taj) || is.null(pi)) {
    message("  Skipping ", chr, " because one or both files are missing.")
    next
  }

  # drop windows with 0 SNPs/variants (can create weird values)
  if ("N_SNPS" %in% names(taj)) {
    taj <- taj %>% filter(N_SNPS > 0)
  }
  if ("N_VARIANTS" %in% names(pi)) {
    pi <- pi %>% filter(N_VARIANTS > 0)
  }

  # add chromosome ID column for later faceting
  taj$CHR_ID <- chr
  pi$CHR_ID  <- chr

  # compute medians (per chromosome)
  taj_median <- median(taj$TajimaD, na.rm = TRUE)
  pi_median  <- median(pi$PI,       na.rm = TRUE)

  message(sprintf("  Tajima's D median for %s: %.3f", chr, taj_median))
  message(sprintf("  pi median for %s: %.5f", chr, pi_median))

  # Tajima's D plot with:
  # - dashed line at 0
  # - solid line at median
  p_taj <- ggplot(taj, aes(x = BIN_START, y = TajimaD)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, alpha = 0.7) +
    geom_hline(yintercept = taj_median, color = "red", linewidth = 0.5, alpha = 0.7) +
    geom_point(alpha = 0.5, size = 0.7) +
    labs(
      title = paste("Tajima's D -", chr),
      subtitle = sprintf("Median Tajima's D = %.3f", taj_median),
      x = "Position (bp)",
      y = "Tajima's D"
    ) +
    theme_bw()

  taj_outfile <- file.path(out_dir, paste0("TajimaD_", chr, ".png"))
  ggsave(taj_outfile, p_taj, width = 7, height = 4, dpi = 300)
  message("  Saved: ", taj_outfile)

  # pi plot with:
  # - solid line at median pi
  p_pi <- ggplot(pi, aes(x = BIN_START, y = PI)) +
    geom_hline(yintercept = pi_median, color = "red", linewidth = 0.5, alpha = 0.7) +
    geom_point(alpha = 0.5, size = 0.7) +
    labs(
      title = paste("Nucleotide diversity (pi) -", chr),
      subtitle = sprintf("Median pi = %.5f", pi_median),
      x = "Position (bp)",
      y = "pi"
    ) +
    theme_bw()

  pi_outfile <- file.path(out_dir, paste0("pi_", chr, ".png"))
  ggsave(pi_outfile, p_pi, width = 7, height = 4, dpi = 300)
  message("  Saved: ", pi_outfile)

  # store for combined plots
  all_taj_list[[chr]] <- taj
  all_pi_list[[chr]]  <- pi
}

# build combined faceted plots if we have any data

if (length(all_taj_list) > 0) {
  all_taj <- dplyr::bind_rows(all_taj_list)

  # per-chromosome medians for hlines in facets
  taj_meds <- all_taj %>%
    group_by(CHR_ID) %>%
    summarise(taj_median = median(TajimaD, na.rm = TRUE), .groups = "drop")

  p_taj_all <- ggplot(all_taj, aes(x = BIN_START, y = TajimaD)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, alpha = 0.7) +
    geom_hline(
      data = taj_meds,
      aes(yintercept = taj_median),
      color = "red", linewidth = 0.4, alpha = 0.7,
      inherit.aes = FALSE
    ) +
    geom_point(alpha = 0.4, size = 0.5) +
    facet_wrap(~ CHR_ID, scales = "free_x") +
    labs(
      title = "Tajima's D across chromosomes",
      subtitle = "Red line = chromosome-wide median; dashed line = 0",
      x = "Position (bp)",
      y = "Tajima's D"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 8),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6)
    )

  taj_all_out <- file.path(out_dir, paste0("All_TajimaD_faceted_", results_label, ".png"))
  ggsave(taj_all_out, p_taj_all, width = 10, height = 8, dpi = 300)
  message("Saved combined Tajima's D plot: ", taj_all_out)
} else {
  message("No TajimaD data collected; skipping combined Tajima's D plot.")
}

if (length(all_pi_list) > 0) {
  all_pi <- dplyr::bind_rows(all_pi_list)

  # per-chromosome medians for hlines in facets
  pi_meds <- all_pi %>%
    group_by(CHR_ID) %>%
    summarise(pi_median = median(PI, na.rm = TRUE), .groups = "drop")

  p_pi_all <- ggplot(all_pi, aes(x = BIN_START, y = PI)) +
    geom_hline(
      data = pi_meds,
      aes(yintercept = pi_median),
      color = "red", linewidth = 0.4, alpha = 0.7,
      inherit.aes = FALSE
    ) +
    geom_point(alpha = 0.4, size = 0.5) +
    facet_wrap(~ CHR_ID, scales = "free_x") +
    labs(
      title = "Nucleotide diversity (pi) across chromosomes",
      subtitle = "Red line = chromosome-wide median",
      x = "Position (bp)",
      y = "pi"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 8),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6)
    )

  pi_all_out  <- file.path(out_dir, paste0("All_pi_faceted_", results_label, ".png"))
  ggsave(pi_all_out, p_pi_all, width = 10, height = 8, dpi = 300)
  message("Saved combined pi plot: ", pi_all_out)
} else {
  message("No pi data collected; skipping combined pi plot.")
}

message("Done. Plots written to: ", out_dir)