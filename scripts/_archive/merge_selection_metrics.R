#!/usr/bin/env Rscript

# title: 'merge_selection_metrics.R'
# project: 'BIOL624 Final Project: Selection Detection in Lady Beetles'
# author: 'Reina Hastings'
# contact: 'reinahastings13@gmail.com'
# date created: 11/17/2025
# last modified: 11/17/2025
# purpose: 'Merge π, Tajima’s D, and FST per window across chromosomes.'
# usage:
#   Rscript merge_selection_metrics.R -o ../results/merged_AllPops
#   Rscript merge_selection_metrics.R -o ../results/merged_AllPops NC_058189.1.filtered
#   Rscript merge_selection_metrics.R -o ../results/merged_AllPops all

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("Usage:\n",
      "  Rscript merge_selection_metrics.R -o <output_dir> [chr1 chr2 ...]\n",
      "  Rscript merge_selection_metrics.R -o <output_dir> all\n\n",
      "Assumes input structure under ../results:\n",
      "  results_Europe/{tajima,pi}\n",
      "  results_Asia/{tajima,pi}\n",
      "  results_North_America/{tajima,pi}\n",
      "  fst_Eur_vs_Asia, fst_Eur_vs_NA, fst_Asia_vs_NA\n")
  quit(status = 1)
}

# parse arguments
outdir <- NULL
chrom_args <- c()
i <- 1

while (i <= length(args)) {
  if (args[i] %in% c("-o", "--outdir")) {
    if (i + 1 > length(args)) {
      stop("Error: -o requires a directory argument.")
    }
    outdir <- args[i + 1]
    i <- i + 2
  } else {
    # any argument NOT a flag is interpreted as a chromosome name
    chrom_args <- c(chrom_args, args[i])
    i <- i + 1
  }
}

if (is.null(outdir)) {
  stop("Error: You must provide an output directory using -o <dir>")
}

# create output dir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(outdir, "figs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", outdir)

# define input directories (fixed layout)
base_results <- "../results"

pop_dirs <- list(
  Europe       = file.path(base_results, "results_Europe"),
  Asia         = file.path(base_results, "results_Asia"),
  NorthAmerica = file.path(base_results, "results_North_America")
)

fst_dirs <- list(
  Eur_Asia = file.path(base_results, "fst_Eur_vs_Asia"),
  Eur_NA   = file.path(base_results, "fst_Eur_vs_NA"),
  Asia_NA  = file.path(base_results, "fst_Asia_vs_NA")
)

# helper to check dir and warn
check_dir <- function(path, label) {
  if (!dir.exists(path)) {
    warning("Missing directory for ", label, ": ", path)
    return(FALSE)
  }
  TRUE
}

# discover chromosomes
chrom_sets <- list()

# use Europe Tajima as main reference if available
if (check_dir(file.path(pop_dirs$Europe, "tajima"), "Europe TajimaD")) {
  taj_files <- list.files(file.path(pop_dirs$Europe, "tajima"),
                          pattern = "\\.Tajima\\.D$", full.names = FALSE)
  if (length(taj_files) > 0) {
    chrom_sets[["Europe_Tajima"]] <- sub("\\.Tajima\\.D$", "", taj_files)
  }
}

# also include any FST chromosomes (in case a pop dir is missing)
for (comp in names(fst_dirs)) {
  if (check_dir(fst_dirs[[comp]], paste0("FST ", comp))) {
    fst_files <- list.files(fst_dirs[[comp]],
                            pattern = "_FST\\.windowed\\.weir\\.fst$",
                            full.names = FALSE)
    if (length(fst_files) > 0) {
      chrom_sets[[paste0("FST_", comp)]] <- sub("_FST\\.windowed\\.weir\\.fst$", "", fst_files)
    }
  }
}

if (length(chrom_sets) == 0) {
  stop("Could not find any TajimaD or FST files in expected locations under ../results")
}

all_chroms <- sort(unique(unlist(chrom_sets)))

if (length(chrom_args) == 0 || (length(chrom_args) == 1 && chrom_args[1] == "all")) {
  chroms <- all_chroms
  message("Merging all chromosomes: ", paste(chroms, collapse = ", "))
} else {
  chroms <- chrom_args
  missing <- setdiff(chroms, all_chroms)
  if (length(missing) > 0) {
    warning("These chromosomes were not found in inputs: ",
            paste(missing, collapse = ", "))
    chroms <- intersect(chroms, all_chroms)
  }
  message("Merging chromosomes: ", paste(chroms, collapse = ", "))
}

# helpers to read metrics

read_tajima <- function(pop_label, chr) {
  pdir <- pop_dirs[[pop_label]]
  tdir <- file.path(pdir, "tajima")
  f <- file.path(tdir, paste0(chr, ".Tajima.D"))
  if (!file.exists(f)) return(NULL)
  df <- read.table(f, header = TRUE)
  # expect CHROM, BIN_START, N_SNPS, TajimaD
  df %>%
    dplyr::select(CHROM, BIN_START, TajimaD) %>%
    dplyr::rename(!!paste0("TajD_", pop_label) := TajimaD)
}

read_pi <- function(pop_label, chr) {
  pdir <- pop_dirs[[pop_label]]
  pidir <- file.path(pdir, "pi")
  f <- file.path(pidir, paste0(chr, ".windowed.pi"))
  if (!file.exists(f)) return(NULL)
  df <- read.table(f, header = TRUE)
  # expect CHROM, BIN_START, N_VARIANTS, PI (vcftools default)
  col_pi <- if ("PI" %in% names(df)) "PI" else names(df)[grepl("PI", names(df), ignore.case = TRUE)][1]
  df %>%
    dplyr::select(CHROM, BIN_START, !!sym(col_pi)) %>%
    dplyr::rename(!!paste0("pi_", pop_label) := !!sym(col_pi))
}

read_fst <- function(comp_label, chr) {
  fdir <- fst_dirs[[comp_label]]
  f <- file.path(fdir, paste0(chr, "_FST.windowed.weir.fst"))
  if (!file.exists(f)) return(NULL)
  df <- read.table(f, header = TRUE)
  # expect CHROM, BIN_START, ..., WEIGHTED_FST
  df %>%
    dplyr::select(CHROM, BIN_START, WEIGHTED_FST) %>%
    dplyr::rename(!!paste0("FST_", comp_label) := WEIGHTED_FST)
}

# merge per chromosome
all_merged <- list()

for (chr in chroms) {
  message("Merging metrics for ", chr, " ...")

  merged_chr <- NULL

  # Tajima & pi per pop
  for (pop in names(pop_dirs)) {
    tdf <- read_tajima(pop, chr)
    pdf <- read_pi(pop, chr)

    for (df in list(tdf, pdf)) {
      if (!is.null(df)) {
        if (is.null(merged_chr)) {
          merged_chr <- df
        } else {
          merged_chr <- full_join(merged_chr, df, by = c("CHROM", "BIN_START"))
        }
      }
    }
  }

  # FST per comparison
  for (comp in names(fst_dirs)) {
    fdf <- read_fst(comp, chr)
    if (!is.null(fdf)) {
      if (is.null(merged_chr)) {
        merged_chr <- fdf
      } else {
        merged_chr <- full_join(merged_chr, fdf, by = c("CHROM", "BIN_START"))
      }
    }
  }

  if (is.null(merged_chr)) {
    warning("No data found for ", chr, "; skipping.")
    next
  }

  merged_chr <- merged_chr %>% arrange(BIN_START)

  # write per-chromosome CSV
  out_csv <- file.path(outdir, paste0("merged_metrics_", chr, ".csv"))
  write.csv(merged_chr, out_csv, row.names = FALSE)
  message("  Wrote: ", out_csv)

  # simple multi-metric plot (one facet per metric)
  long <- merged_chr %>%
    dplyr::select(-CHROM) %>%
    pivot_longer(-BIN_START, names_to = "metric", values_to = "value")

  p <- ggplot(long, aes(x = BIN_START, y = value)) +
    geom_point(alpha = 0.4, size = 0.4) +
    facet_wrap(~metric, ncol = 1, scales = "free_y") +
    labs(
      title = paste("Selection metrics -", chr),
      x = "Position (bp)",
      y = "Value"
    ) +
    theme_bw()

  out_png <- file.path(fig_dir, paste0("merged_tracks_", chr, ".png"))
  ggsave(out_png, p, width = 8, height = 10, dpi = 300)
  message("  Wrote plot: ", out_png)

  merged_chr$CHR_ID <- chr
  all_merged[[chr]] <- merged_chr
}

if (length(all_merged) > 0) {
  merged_all <- bind_rows(all_merged)
  out_all_csv <- file.path(outdir, "merged_metrics_all_chromosomes.csv")
  write.csv(merged_all, out_all_csv, row.names = FALSE)
  message("Wrote merged metrics for all chromosomes: ", out_all_csv)
} else {
  message("No merged data generated (no chromosomes had overlapping metrics).")
}

message("Done.")
