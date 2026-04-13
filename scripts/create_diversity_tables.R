#!/usr/bin/env Rscript

# title: create_diversity_tables.R
# purpose: Generate publication-quality table figures for Tajima's D and Pi
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
base_dir <- "../results/whole_pop_results"
tajima_file <- file.path(base_dir, "summary/tajima_d_table.tsv")
pi_file <- file.path(base_dir, "summary/pi_table.tsv")
output_dir <- file.path(base_dir, "summary")

# Read data
cat("Reading Tajima's D data...\n")
tajima_data <- read.delim(tajima_file, header = TRUE)

cat("Reading Pi data...\n")
pi_data <- read.delim(pi_file, header = TRUE)

# Create chromosome labels (shorter names for display)
chr_labels <- data.frame(
  CHROM = paste0("NC_0581", 89:98, ".1"),
  Chr = paste0("Chr ", 1:10),
  chr_order = 1:10
)

# Calculate summary statistics for Tajima's D
tajima_summary <- tajima_data %>%
  filter(!is.na(TajimaD) & TajimaD != "nan") %>%
  mutate(TajimaD = as.numeric(TajimaD)) %>%
  group_by(CHROM) %>%
  summarise(
    `N Windows` = n(),
    Mean = round(mean(TajimaD, na.rm = TRUE), 4),
    Median = round(median(TajimaD, na.rm = TRUE), 4),
    SD = round(sd(TajimaD, na.rm = TRUE), 4),
    Min = round(min(TajimaD, na.rm = TRUE), 4),
    Max = round(max(TajimaD, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  left_join(chr_labels, by = "CHROM") %>%
  arrange(chr_order) %>%
  select(Chr, `N Windows`, Mean, Median, SD, Min, Max)

# Calculate summary statistics for Pi
pi_summary <- pi_data %>%
  filter(!is.na(PI)) %>%
  group_by(CHROM) %>%
  summarise(
    `N Windows` = n(),
    Mean = sprintf("%.2e", mean(PI, na.rm = TRUE)),
    Median = sprintf("%.2e", median(PI, na.rm = TRUE)),
    SD = sprintf("%.2e", sd(PI, na.rm = TRUE)),
    Min = sprintf("%.2e", min(PI, na.rm = TRUE)),
    Max = sprintf("%.2e", max(PI, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(chr_labels, by = "CHROM") %>%
  arrange(chr_order) %>%
  select(Chr, `N Windows`, Mean, Median, SD, Min, Max)

# Add genome-wide summary row for Tajima's D
tajima_genomewide <- tajima_data %>%
  filter(!is.na(TajimaD) & TajimaD != "nan") %>%
  mutate(TajimaD = as.numeric(TajimaD)) %>%
  summarise(
    Chr = "Genome-wide",
    `N Windows` = n(),
    Mean = round(mean(TajimaD, na.rm = TRUE), 4),
    Median = round(median(TajimaD, na.rm = TRUE), 4),
    SD = round(sd(TajimaD, na.rm = TRUE), 4),
    Min = round(min(TajimaD, na.rm = TRUE), 4),
    Max = round(max(TajimaD, na.rm = TRUE), 4)
  )

tajima_summary <- bind_rows(tajima_summary, tajima_genomewide)

# Add genome-wide summary row for Pi
pi_genomewide <- pi_data %>%
  filter(!is.na(PI)) %>%
  summarise(
    Chr = "Genome-wide",
    `N Windows` = n(),
    Mean = sprintf("%.2e", mean(PI, na.rm = TRUE)),
    Median = sprintf("%.2e", median(PI, na.rm = TRUE)),
    SD = sprintf("%.2e", sd(PI, na.rm = TRUE)),
    Min = sprintf("%.2e", min(PI, na.rm = TRUE)),
    Max = sprintf("%.2e", max(PI, na.rm = TRUE))
  )

pi_summary <- bind_rows(pi_summary, pi_genomewide)

# Function to create publication-quality table as a figure
create_table_figure <- function(data, title, filename, footnote = NULL) {

  # Theme for the table
  tt <- ttheme_minimal(
    core = list(
      fg_params = list(fontsize = 10, fontface = "plain"),
      bg_params = list(fill = c("grey95", "white"), col = NA),
      padding = unit(c(4, 4), "mm")
    ),
    colhead = list(
      fg_params = list(fontsize = 11, fontface = "bold"),
      bg_params = list(fill = "grey80", col = NA),
      padding = unit(c(4, 4), "mm")
    )
  )

  # Highlight genome-wide row
  n_rows <- nrow(data)
  core_bg <- matrix("white", nrow = n_rows, ncol = ncol(data))
  core_bg[seq(1, n_rows, 2), ] <- "grey95"
  core_bg[n_rows, ] <- "lightblue"  # Highlight genome-wide row

  tt_custom <- ttheme_minimal(
    core = list(
      fg_params = list(fontsize = 10),
      bg_params = list(fill = core_bg, col = NA),
      padding = unit(c(4, 4), "mm")
    ),
    colhead = list(
      fg_params = list(fontsize = 11, fontface = "bold"),
      bg_params = list(fill = "grey70", col = NA),
      padding = unit(c(4, 4), "mm")
    )
  )

  # Create table grob
  table_grob <- tableGrob(data, rows = NULL, theme = tt_custom)

  # Add title
  title_grob <- textGrob(title, gp = gpar(fontsize = 14, fontface = "bold"))

  # Add footnote if provided
  if (!is.null(footnote)) {
    footnote_grob <- textGrob(footnote, x = 0, hjust = 0,
                               gp = gpar(fontsize = 9, fontface = "italic"))

    # Combine elements
    final_grob <- arrangeGrob(
      title_grob,
      table_grob,
      footnote_grob,
      heights = unit(c(1, 0.8, 0.5), "null"),
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

  # Calculate appropriate dimensions
  width <- 8
  height <- 0.4 * nrow(data) + 1.5

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

# Create Tajima's D table figure
cat("\nCreating Tajima's D table figure...\n")
create_table_figure(
  tajima_summary,
  expression(bold("Table 1. ") * "Summary statistics of Tajima's D across chromosomes"),
  "tajima_d_summary_table",
  "Note: Statistics calculated from 5 kb non-overlapping windows. Genome-wide row highlighted in blue."
)

# Create Pi table figure
cat("\nCreating Pi table figure...\n")
create_table_figure(
  pi_summary,
  "Table 2. Summary statistics of nucleotide diversity (pi) across chromosomes",
  "pi_summary_table",
  "Note: Statistics calculated from 5 kb non-overlapping windows. Values in scientific notation. Genome-wide row highlighted in blue."
)

# Also save as TSV for reference
write.table(tajima_summary, file.path(output_dir, "tajima_d_summary_stats.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pi_summary, file.path(output_dir, "pi_summary_stats.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nDone! Summary statistics also saved as TSV files.\n")
cat("\nTajima's D Summary:\n")
print(tajima_summary, row.names = FALSE)
cat("\nPi Summary:\n
")
print(pi_summary, row.names = FALSE)
