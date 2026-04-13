#!/usr/bin/env Rscript

# title: create_population_pi_table.R
# purpose: Generate publication-quality Pi table for all populations
#          Excludes sex chromosome (NC_058198.1) from genome-wide calculations
# author: Generated for BIOL624 Project

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

# Set paths
results_dir <- "../results"
output_dir <- file.path(results_dir, "whole_pop_results/summary")

# Define populations and their directories
populations <- c("CHI", "EEU", "WEU", "USA", "Whole Pop")
pop_dirs <- c(
  CHI = "CHI_results",
  EEU = "EEU_results",
  WEU = "WEU_results",
  USA = "USA_results",
  `Whole Pop` = "whole_pop_results"
)

# Sex chromosome to exclude from genome-wide calculations
sex_chrom <- "NC_058198.1"

# Create chromosome labels (shorter names for display)
chr_labels <- data.frame(
  CHROM = paste0("NC_0581", 89:98, ".1"),
  Chr = paste0("Chr ", 1:10),
  chr_order = 1:10,
  is_sex_chrom = c(rep(FALSE, 9), TRUE)  # Chr 10 is sex chromosome
)

# Function to read and combine pi files for a population
read_pop_pi <- function(pop_name, pop_dir) {
  pi_dir <- file.path(results_dir, pop_dir, "pi")

  if (!dir.exists(pi_dir)) {
    warning(paste("Directory not found:", pi_dir))
    return(NULL)
  }

  # Find pi files
  if (pop_name == "Whole Pop") {
    pi_files <- list.files(pi_dir, pattern = "\\.windowed\\.pi$", full.names = TRUE)
  } else {
    pi_files <- list.files(pi_dir, pattern = paste0("_", pop_name, "\\.windowed\\.pi$"), full.names = TRUE)
  }

  if (length(pi_files) == 0) {
    warning(paste("No pi files found in:", pi_dir))
    return(NULL)
  }

  # Read and combine all files
  all_data <- lapply(pi_files, function(f) {
    read.delim(f, header = TRUE)
  })

  combined <- bind_rows(all_data)
  combined$Population <- pop_name
  return(combined)
}

# Read data for all populations
cat("Reading Pi data for all populations...\n")
all_pi_data <- list()

for (pop in names(pop_dirs)) {
  cat("  Reading", pop, "...\n")
  data <- read_pop_pi(pop, pop_dirs[pop])
  if (!is.null(data)) {
    all_pi_data[[pop]] <- data
  }
}

# Combine all population data
combined_pi <- bind_rows(all_pi_data)

# Calculate summary statistics for each population and chromosome
pi_summary <- combined_pi %>%
  filter(!is.na(PI)) %>%
  left_join(chr_labels, by = "CHROM") %>%
  group_by(Population, Chr, chr_order, is_sex_chrom) %>%
  summarise(
    N_Windows = n(),
    Mean = mean(PI, na.rm = TRUE),
    Median = median(PI, na.rm = TRUE),
    SD = sd(PI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Population, chr_order)

# Calculate genome-wide summary EXCLUDING sex chromosome
cat("\nCalculating genome-wide statistics (excluding sex chromosome)...\n")

genome_wide_summary <- combined_pi %>%
  filter(!is.na(PI)) %>%
  left_join(chr_labels, by = "CHROM") %>%
  filter(!is_sex_chrom) %>%  # Exclude sex chromosome
  group_by(Population) %>%
  summarise(
    Chr = "Autosomal",
    chr_order = 11,
    is_sex_chrom = FALSE,
    N_Windows = n(),
    Mean = mean(PI, na.rm = TRUE),
    Median = median(PI, na.rm = TRUE),
    SD = sd(PI, na.rm = TRUE),
    .groups = "drop"
  )

# Combine per-chromosome and genome-wide
pi_summary_full <- bind_rows(pi_summary, genome_wide_summary) %>%
  arrange(Population, chr_order)

# Create wide format table for publication (populations as columns)
pi_wide <- pi_summary_full %>%
  mutate(
    Mean_fmt = sprintf("%.2e", Mean),
    Median_fmt = sprintf("%.2e", Median)
  ) %>%
  select(Chr, Population, Mean_fmt, Median_fmt, N_Windows) %>%
  pivot_wider(
    names_from = Population,
    values_from = c(Mean_fmt, Median_fmt, N_Windows),
    names_glue = "{Population}_{.value}"
  )

# Reorder columns for better presentation
# Format: Chr, then for each pop: Mean, Median
pop_order <- c("CHI", "EEU", "WEU", "USA", "Whole Pop")

# Create a cleaner table with just means for each population
pi_means_table <- pi_summary_full %>%
  mutate(Mean_fmt = sprintf("%.2e", Mean)) %>%
  select(Chr, Population, Mean_fmt) %>%
  pivot_wider(names_from = Population, values_from = Mean_fmt) %>%
  left_join(chr_labels %>% select(Chr, chr_order), by = "Chr") %>%
  mutate(chr_order = ifelse(Chr == "Autosomal", 11, chr_order)) %>%
  distinct(Chr, .keep_all = TRUE) %>%
  arrange(chr_order) %>%
  select(Chr, all_of(pop_order))

# Create table with means and sample sizes
pi_detailed_table <- pi_summary_full %>%
  mutate(
    Value = sprintf("%.2e (%d)", Mean, N_Windows)
  ) %>%
  select(Chr, Population, Value) %>%
  pivot_wider(names_from = Population, values_from = Value) %>%
  left_join(chr_labels %>% select(Chr, chr_order), by = "Chr") %>%
  mutate(chr_order = ifelse(Chr == "Autosomal", 11, chr_order)) %>%
  distinct(Chr, .keep_all = TRUE) %>%
  arrange(chr_order) %>%
  select(Chr, all_of(pop_order))

# Function to create publication-quality table
create_pop_table_figure <- function(data, title, filename, footnote = NULL) {

  n_rows <- nrow(data)
  n_cols <- ncol(data)

  # Create background colors
  core_bg <- matrix("white", nrow = n_rows, ncol = n_cols)
  core_bg[seq(1, n_rows, 2), ] <- "grey95"

  # Highlight sex chromosome row (Chr 10) in light red
  chr10_row <- which(data$Chr == "Chr 10")
  if (length(chr10_row) > 0) {
    core_bg[chr10_row, ] <- "#FFE4E1"  # Misty rose for sex chromosome
  }

  # Highlight autosomal row in light blue
  autosomal_row <- which(data$Chr == "Autosomal")
  if (length(autosomal_row) > 0) {
    core_bg[autosomal_row, ] <- "lightblue"
  }

  tt_custom <- ttheme_minimal(
    core = list(
      fg_params = list(fontsize = 9),
      bg_params = list(fill = core_bg, col = NA),
      padding = unit(c(3, 3), "mm")
    ),
    colhead = list(
      fg_params = list(fontsize = 10, fontface = "bold"),
      bg_params = list(fill = "grey70", col = NA),
      padding = unit(c(3, 3), "mm")
    )
  )

  # Create table grob
  table_grob <- tableGrob(data, rows = NULL, theme = tt_custom)

  # Add title
  title_grob <- textGrob(title, gp = gpar(fontsize = 12, fontface = "bold"))

  # Add footnote
  if (!is.null(footnote)) {
    footnote_grob <- textGrob(footnote, x = 0, hjust = 0,
                               gp = gpar(fontsize = 8, fontface = "italic"))

    final_grob <- arrangeGrob(
      title_grob,
      table_grob,
      footnote_grob,
      heights = unit(c(0.8, 5, 0.8), "null"),
      padding = unit(1, "line")
    )
  } else {
    final_grob <- arrangeGrob(
      title_grob,
      table_grob,
      heights = unit(c(0.8, 5), "null"),
      padding = unit(1, "line")
    )
  }

  # Calculate dimensions
  width <- 10
  height <- 0.35 * nrow(data) + 1.8

  # Save as PDF
  pdf_file <- file.path(output_dir, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  grid.draw(final_grob)
  dev.off()
  cat("Saved:", pdf_file, "\n")

  # Save as PNG
  png_file <- file.path(output_dir, paste0(filename, ".png"))
  png(png_file, width = width, height = height, units = "in", res = 300)
  grid.draw(final_grob)
  dev.off()
  cat("Saved:", png_file, "\n")
}

# Create the mean Pi table
cat("\nCreating population Pi table (means only)...\n")
create_pop_table_figure(
  pi_means_table,
  "Table. Mean nucleotide diversity (pi) across populations and chromosomes",
  "pi_population_means_table",
  "Note: Values are mean pi per 5 kb window. Chr 10 (sex chromosome, pink) excluded from Autosomal summary (blue). Populations: CHI=China, EEU=Eastern Europe, WEU=Western Europe, USA=North America."
)

# Create detailed table with sample sizes
cat("\nCreating detailed population Pi table...\n")
create_pop_table_figure(
  pi_detailed_table,
  "Table. Mean nucleotide diversity (pi) with window counts across populations",
  "pi_population_detailed_table",
  "Note: Values shown as mean pi (N windows). Chr 10 (sex chromosome, pink) excluded from Autosomal summary (blue). Populations: CHI=China, EEU=Eastern Europe, WEU=Western Europe, USA=North America."
)

# Save as TSV
write.table(pi_means_table, file.path(output_dir, "pi_population_means.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pi_detailed_table, file.path(output_dir, "pi_population_detailed.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary
cat("\n=== Population Pi Summary ===\n")
cat("Sex chromosome (NC_058198.1 = Chr 10) EXCLUDED from Autosomal calculations\n\n")
print(pi_means_table)

cat("\n=== Autosomal Mean Pi by Population ===\n")
genome_wide_summary %>%
  select(Population, Mean, N_Windows) %>%
  mutate(Mean = sprintf("%.4e", Mean)) %>%
  print()
