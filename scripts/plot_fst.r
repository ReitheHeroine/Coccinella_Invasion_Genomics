#!/usr/bin/env Rscript

# title: plot_fst.R
# project: BIOL624 Final Project — Selection Detection in Lady Beetles
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2025-11-17
# last modified: 2025-11-20
#
# purpose:
#   Plot windowed FST for pairwise comparisons between populations, using
#   vcftools --weir-fst-pop windowed outputs.
#
# inputs:
#   - results_dir: directory containing FST files named:
#       <chr>_FST.windowed.weir.fst
#     e.g., ../results/fst_CHI_vs_WEU_results
#   - chromosome arguments:
#       * one or more chromosome IDs (matching <chr> in the filenames), or
#       * the keyword 'all' to auto-detect all chromosomes from the FST files
#
# outputs:
#   - <results_dir>/figs/fst/FST_<chr>_<label>.png
#   - <results_dir>/figs/fst/All_FST_faceted_<label>.png
#
# required packages:
#   - tidyverse
#
# usage examples:
#   # Plot all chromosomes for CHI vs WEU comparison
#   Rscript plot_fst.R ../results/fst_CHI_vs_WEU_results all
#
#   # Plot a single chromosome
#   Rscript plot_fst.R ../results/fst_CHI_vs_WEU_results NC_058189.1.filtered

suppressPackageStartupMessages({
  library(tidyverse)
})

print_usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript plot_fst.R <results_dir> <chr1|all> [<chr2> ...]\n\n",
    "Arguments:\n",
    "  results_dir   Directory containing *_FST.windowed.weir.fst files\n",
    "  chr1|all      One chromosome ID (e.g., NC_058189.1.filtered) or 'all'\n",
    "  chr2 ...      Optional additional chromosome IDs\n\n",
    "Examples:\n",
    "  # CHI vs WEU, all chromosomes\n",
    "  Rscript plot_fst.R ../results/fst_CHI_vs_WEU_results all\n\n",
    "  # USA vs WEU, one chromosome\n",
    "  Rscript plot_fst.R ../results/fst_USA_vs_WEU_results NC_058189.1.filtered\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  print_usage()
  quit(save = "no", status = 1)
}

results_dir <- args[1]
chrom_args  <- args[-1]

fst_dir <- results_dir

# output directory
out_dir <- file.path(results_dir, "figs", "fst")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# derive comparison label from results_dir
# e.g., "../results/fst_CHI_vs_WEU_results" -> "CHI_vs_WEU"
results_label <- basename(normalizePath(results_dir))
results_label <- sub("^results[_-]?", "", results_label)
results_label <- sub("^fst[_-]?", "", results_label)
results_label <- sub("_results$", "", results_label)

message("Comparison label: ", results_label)

# identify chromosomes to plot
fst_files <- list.files(fst_dir, pattern = "_FST\\.windowed\\.weir\\.fst$", full.names = TRUE)

if (length(fst_files) == 0) {
  stop("No FST files found in: ", fst_dir)
}

all_chroms <- sub("_FST\\.windowed\\.weir\\.fst$", "", basename(fst_files))

if (length(chrom_args) == 1 && chrom_args[1] == "all") {
  chroms <- all_chroms
  message("Detected chromosomes: ", paste(chroms, collapse = ", "))
} else {
  chroms <- chrom_args
  missing <- setdiff(chroms, all_chroms)
  if (length(missing) > 0) {
    stop("These chromosomes were not found in results: ", paste(missing, collapse = ", "))
  }
  message("Using chromosomes: ", paste(chroms, collapse = ", "))
}

# helper to read FST file
read_fst_file <- function(chr) {
  fpath <- file.path(fst_dir, paste0(chr, "_FST.windowed.weir.fst"))
  if (!file.exists(fpath)) {
    warning("Missing FST file for chromosome: ", chr, " (", fpath, ")")
    return(NULL)
  }
  df <- read.table(fpath, header = TRUE)
  df$CHR_ID <- chr
  df
}

all_fst_list <- list()

# per-chromosome plots
for (chr in chroms) {

  fst <- read_fst_file(chr)
  if (is.null(fst)) next

  # drop NA values (windows with no variants / undefined FST)
  fst <- fst %>% filter(!is.na(WEIGHTED_FST))

  median_fst <- median(fst$WEIGHTED_FST, na.rm = TRUE)

  message(sprintf("Chrom %s median WEIGHTED_FST: %.4f", chr, median_fst))

  p_chr <- ggplot(fst, aes(x = BIN_START, y = WEIGHTED_FST)) +
    geom_hline(yintercept = median_fst, color = "red", alpha = 0.7) +
    geom_point(alpha = 0.5, size = 0.6) +
    labs(
      title    = paste("Windowed FST -", chr),
      subtitle = paste("Comparison:", results_label,
                       "| median WEIGHTED_FST =", round(median_fst, 4)),
      x = "Position (bp)",
      y = "Weighted FST"
    ) +
    theme_bw()

  outfile <- file.path(out_dir, paste0("FST_", chr, "_", results_label, ".png"))
  ggsave(outfile, p_chr, width = 7, height = 4, dpi = 300)
  message("Saved: ", outfile)

  all_fst_list[[chr]] <- fst
}

# combined faceted plot
if (length(all_fst_list) > 0) {

  all_fst <- bind_rows(all_fst_list)

  meds <- all_fst %>%
    group_by(CHR_ID) %>%
    summarize(median_fst = median(WEIGHTED_FST, na.rm = TRUE), .groups = "drop")

  p_all <- ggplot(all_fst, aes(x = BIN_START, y = WEIGHTED_FST)) +
    geom_hline(
      data = meds,
      aes(yintercept = median_fst),
      color = "red",
      alpha = 0.7
    ) +
    geom_point(alpha = 0.4, size = 0.4) +
    facet_wrap(~ CHR_ID, scales = "free_x") +
    labs(
      title    = "Windowed FST across chromosomes",
      subtitle = paste("Comparison:", results_label),
      x        = "Position (bp)",
      y        = "Weighted FST"
    ) +
    theme_bw() +
    theme(strip.text = element_text(size = 8))

  outfile_all <- file.path(out_dir, paste0("All_FST_faceted_", results_label, ".png"))
  ggsave(outfile_all, p_all, width = 10, height = 8, dpi = 300)
  message("Saved combined plot: ", outfile_all)

} else {
  message("No FST data found. No combined plot created.")
}

message("Done. All plots saved in: ", out_dir)