#!/usr/bin/env Rscript

# title: create_population_tajima_table.R
# purpose: Generate publication-quality Tajima's D table for all populations
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

# Function to read and combine Tajima's D files for a population
read_pop_tajima <- function(pop_name, pop_dir) {
  tajima_dir <- file.path(results_dir, pop_dir, "tajima")

  if (!dir.exists(tajima_dir)) {
    warning(paste("Directory not found:", tajima_dir))
    return(NULL)
  }

  # Find Tajima's D files
  if (pop_name == "Whole Pop") {
    tajima_files <- list.files(tajima_dir, pattern = "\\.Tajima\\.D$", full.names = TRUE)
    # Exclude population-specific files
    tajima_files <- tajima_files[!grepl("_CHI|_EEU|_WEU|_USA", tajima_files)]
  } else {
    tajima_files <- list.files(tajima_dir, pattern = paste0("_", pop_name, "\\.Tajima\\.D$"), full.names = TRUE)
  }

  if (length(tajima_files) == 0) {
    warning(paste("No Tajima's D files found in:", tajima_dir))
    return(NULL)
  }

  # Read and combine all files
  all_data <- lapply(tajima_files, function(f) {
    read.delim(f, header = TRUE)
  })

  combined <- bind_rows(all_data)
  combined$Population <- pop_name
  return(combined)
}

# Read data for all populations
cat("Reading Tajima's D data for all populations...\n")
all_tajima_data <- list()

for (pop in names(pop_dirs)) {
  cat("  Reading", pop, "...\n")
  data <- read_pop_tajima(pop, pop_dirs[pop])
  if (!is.null(data)) {
    all_tajima_data[[pop]] <- data
  }
}

# Combine all population data
combined_tajima <- bind_rows(all_tajima_data)

# Calculate summary statistics for each population and chromosome
tajima_summary <- combined_tajima %>%
  filter(!is.na(TajimaD) & TajimaD != "nan") %>%
  mutate(TajimaD = as.numeric(TajimaD)) %>%
  left_join(chr_labels, by = "CHROM") %>%
  group_by(Population, Chr, chr_order, is_sex_chrom) %>%
  summarise(
    N_Windows = n(),
    Mean = mean(TajimaD, na.rm = TRUE),
    Median = median(TajimaD, na.rm = TRUE),
    SD = sd(TajimaD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Population, chr_order)

# Calculate genome-wide summary EXCLUDING sex chromosome
cat("\nCalculating autosomal statistics (excluding sex chromosome)...\n")

autosomal_summary <- combined_tajima %>%
  filter(!is.na(TajimaD) & TajimaD != "nan") %>%
  mutate(TajimaD = as.numeric(TajimaD)) %>%
  left_join(chr_labels, by = "CHROM") %>%
  filter(!is_sex_chrom) %>%  # Exclude sex chromosome
  group_by(Population) %>%
  summarise(
    Chr = "Autosomal",
    chr_order = 11,
    is_sex_chrom = FALSE,
    N_Windows = n(),
    Mean = mean(TajimaD, na.rm = TRUE),
    Median = median(TajimaD, na.rm = TRUE),
    SD = sd(TajimaD, na.rm = TRUE),
    .groups = "drop"
  )

# Combine per-chromosome and autosomal
tajima_summary_full <- bind_rows(tajima_summary, autosomal_summary) %>%
  arrange(Population, chr_order)

# Define population order
pop_order <- c("CHI", "EEU", "WEU", "USA", "Whole Pop")

# Create a cleaner table with just means for each population
tajima_means_table <- tajima_summary_full %>%
  mutate(Mean_fmt = sprintf("%.4f", Mean)) %>%
  select(Chr, Population, Mean_fmt) %>%
  pivot_wider(names_from = Population, values_from = Mean_fmt) %>%
  left_join(chr_labels %>% select(Chr, chr_order), by = "Chr") %>%
  mutate(chr_order = ifelse(Chr == "Autosomal", 11, chr_order)) %>%
  distinct(Chr, .keep_all = TRUE) %>%
  arrange(chr_order) %>%
  select(Chr, all_of(pop_order))

# Create table with means and sample sizes
tajima_detailed_table <- tajima_summary_full %>%
  mutate(
    Value = sprintf("%.4f (%d)", Mean, N_Windows)
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

# Create the mean Tajima's D table
cat("\nCreating population Tajima's D table (means only)...\n")
create_pop_table_figure(
  tajima_means_table,
  "Table. Mean Tajima's D across populations and chromosomes",
  "tajima_d_population_means_table",
  "Note: Values are mean Tajima's D per 5 kb window. Chr 10 (sex chromosome, pink) excluded from Autosomal summary (blue). Populations: CHI=China, EEU=Eastern Europe, WEU=Western Europe, USA=North America."
)

# Create detailed table with sample sizes
cat("\nCreating detailed population Tajima's D table...\n")
create_pop_table_figure(
  tajima_detailed_table,
  "Table. Mean Tajima's D with window counts across populations",
  "tajima_d_population_detailed_table",
  "Note: Values shown as mean Tajima's D (N windows). Chr 10 (sex chromosome, pink) excluded from Autosomal summary (blue). Populations: CHI=China, EEU=Eastern Europe, WEU=Western Europe, USA=North America."
)

# Save as TSV
write.table(tajima_means_table, file.path(output_dir, "tajima_d_population_means.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tajima_detailed_table, file.path(output_dir, "tajima_d_population_detailed.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary
cat("\n=== Population Tajima's D Summary ===\n")
cat("Sex chromosome (NC_058198.1 = Chr 10) EXCLUDED from Autosomal calculations\n\n")
print(tajima_means_table)

cat("\n=== Autosomal Mean Tajima's D by Population ===\n")
autosomal_summary %>%
  select(Population, Mean, N_Windows) %>%
  mutate(Mean = sprintf("%.4f", Mean)) %>%
  arrange(match(Population, pop_order)) %>%
  print()
