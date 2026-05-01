#!/usr/bin/env Rscript

# title: summarize_fst.R
# project: BIOL624 Final Project - Selection Detection in Coccinella septempunctata
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-01-13
#
# purpose:
#   Generate a single cumulative FST summary for all pairwise population comparisons.
#   Reads genome-wide and windowed FST outputs from fst.sh.
#
# inputs:
#   - results/fst/ directory containing <pop1>_vs_<pop2>/ subdirectories
#
# outputs:
#   - results/fst/fst_cumulative_summary.tsv
#   - results/fst/fst_cumulative_report.txt
#
# usage:
#   Rscript summarize_fst.R [results_dir]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "../results"

fst_base <- file.path(results_dir, "fst")

if (!dir.exists(fst_base)) {
  stop("FST directory not found: ", fst_base)
}

# Find all comparison directories (format: POP1_vs_POP2)
all_dirs <- list.dirs(fst_base, recursive = FALSE, full.names = TRUE)
fst_dirs <- all_dirs[grepl("^[A-Z]+_vs_[A-Z]+$", basename(all_dirs))]

if (length(fst_dirs) == 0) {
  stop("No comparison directories found in: ", fst_base)
}

cat("Found", length(fst_dirs), "FST comparisons\n")

# Process each comparison
all_results <- list()

for (fst_dir in fst_dirs) {
  comparison <- basename(fst_dir)
  cat("Processing:", comparison, "\n")

  # Read genome-wide FST table
  gw_file <- file.path(fst_dir, "fst_genomewide_table.tsv")

  if (!file.exists(gw_file)) {
    warning("  No genome-wide table found for ", comparison)
    next
  }

  gw_data <- read.delim(gw_file, stringsAsFactors = FALSE)

  gw_data <- gw_data %>%
    filter(!is.na(Weighted_FST) & Weighted_FST != "NA" & N_Sites != "NA") %>%
    mutate(
      Weighted_FST = as.numeric(Weighted_FST),
      Mean_FST = as.numeric(Mean_FST),
      N_Sites = as.numeric(N_Sites)
    )

  if (nrow(gw_data) == 0) {
    warning("  No valid data for ", comparison)
    next
  }

  # Calculate overall statistics
  overall_weighted <- sum(gw_data$Weighted_FST * gw_data$N_Sites) / sum(gw_data$N_Sites)
  overall_mean <- sum(gw_data$Mean_FST * gw_data$N_Sites) / sum(gw_data$N_Sites)
  total_sites <- sum(gw_data$N_Sites)
  n_chromosomes <- nrow(gw_data)

  chr_fst_values <- gw_data$Weighted_FST
  fst_sd <- sd(chr_fst_values)
  fst_min <- min(chr_fst_values)
  fst_max <- max(chr_fst_values)

  # Read windowed FST files
  fst_windowed_dir <- file.path(fst_dir, "fst")
  windowed_files <- list.files(fst_windowed_dir, pattern = "\\.windowed\\.weir\\.fst$", full.names = TRUE)

  windowed_stats <- NULL
  if (length(windowed_files) > 0) {
    all_windows <- lapply(windowed_files, function(f) {
      df <- read.delim(f, stringsAsFactors = FALSE)
      if ("WEIGHTED_FST" %in% names(df)) {
        return(df$WEIGHTED_FST)
      } else if ("MEAN_FST" %in% names(df)) {
        return(df$MEAN_FST)
      }
      return(NULL)
    })
    all_windows <- unlist(all_windows)
    all_windows <- all_windows[!is.na(all_windows) & all_windows != "nan"]

    if (length(all_windows) > 0) {
      windowed_stats <- list(
        n_windows = length(all_windows),
        median = median(all_windows),
        q1 = quantile(all_windows, 0.25),
        q3 = quantile(all_windows, 0.75),
        pct_negative = mean(all_windows < 0) * 100
      )
    }
  }

  # Determine comparison type
  pops <- strsplit(comparison, "_vs_")[[1]]
  comparison_type <- if (pops[1] == "CHI" || pops[2] == "CHI") {
    "Native_vs_Invasive"
  } else {
    "Among_Invasive"
  }

  # Store results
  result <- data.frame(
    Comparison = comparison,
    Type = comparison_type,
    Weighted_FST = round(overall_weighted, 6),
    Mean_FST = round(overall_mean, 6),
    SD = round(fst_sd, 6),
    Min = round(fst_min, 6),
    Max = round(fst_max, 6),
    Total_SNPs = total_sites,
    N_Chromosomes = n_chromosomes,
    stringsAsFactors = FALSE
  )

  if (!is.null(windowed_stats)) {
    result$N_Windows <- windowed_stats$n_windows
    result$Windowed_Median <- round(windowed_stats$median, 6)
    result$Windowed_Q1 <- round(windowed_stats$q1, 6)
    result$Windowed_Q3 <- round(windowed_stats$q3, 6)
    result$Pct_Negative <- round(windowed_stats$pct_negative, 1)
  }

  all_results[[comparison]] <- result
}

# Combine and sort
summary_df <- bind_rows(all_results)
summary_df <- summary_df %>% arrange(desc(Type), Comparison)

# Write outputs to results/fst/
tsv_file <- file.path(fst_base, "fst_cumulative_summary.tsv")
write.table(summary_df, tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote:", tsv_file, "\n")

# Write text report
report_file <- file.path(fst_base, "fst_cumulative_report.txt")
sink(report_file)

cat("================================================================================\n")
cat("CUMULATIVE FST SUMMARY - ALL POPULATION COMPARISONS\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================================\n\n")

cat("METHOD: Weir & Cockerham (1984) weighted FST estimator\n")
cat("SOFTWARE: vcftools\n\n")

cat("================================================================================\n")
cat("GENOME-WIDE FST\n")
cat("================================================================================\n\n")

cat(sprintf("%-15s  %10s  %10s  %10s  %8s  %8s  %10s\n",
            "Comparison", "Weighted", "Mean", "SD", "Min", "Max", "SNPs"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (i in seq_len(nrow(summary_df))) {
  row <- summary_df[i, ]
  cat(sprintf("%-15s  %10.4f  %10.4f  %10.4f  %8.4f  %8.4f  %10d\n",
              row$Comparison, row$Weighted_FST, row$Mean_FST, row$SD,
              row$Min, row$Max, row$Total_SNPs))
}

cat("\n")
cat("================================================================================\n")
cat("WINDOWED FST (5kb windows)\n")
cat("================================================================================\n\n")

if ("N_Windows" %in% names(summary_df)) {
  cat(sprintf("%-15s  %10s  %10s  %10s  %8s  %10s\n",
              "Comparison", "Median", "Q1", "Q3", "% Neg", "N_Windows"))
  cat(paste(rep("-", 70), collapse = ""), "\n")

  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, ]
    if (!is.na(row$N_Windows)) {
      cat(sprintf("%-15s  %10.4f  %10.4f  %10.4f  %7.1f%%  %10d\n",
                  row$Comparison, row$Windowed_Median, row$Windowed_Q1,
                  row$Windowed_Q3, row$Pct_Negative, row$N_Windows))
    }
  }
}

cat("\n")
cat("================================================================================\n")
cat("FILES\n")
cat("================================================================================\n")
cat("Summary table:", tsv_file, "\n")
cat("This report:  ", report_file, "\n")

sink()

cat("Wrote:", report_file, "\n")
cat("\nDone!\n")
