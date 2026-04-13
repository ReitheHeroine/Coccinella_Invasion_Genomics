#!/usr/bin/env Rscript

# title: go_enrichment_analysis.R
# purpose: Perform GO enrichment analysis on outlier genes using clusterProfiler
#          Uses Tribolium GO annotations transferred via BLAST orthology
# author: Generated for BIOL624 Project

# Load required libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(enrichplot)
})

# Set paths
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  gene2go_file <- args[1]
} else {
  gene2go_file <- "../results/whole_pop_results/combined_outliers/outlier_gene2go.tsv"
}

results_dir <- dirname(gene2go_file)
go_data_dir <- "../data/go_annotations"

cat("GO Enrichment Analysis for Outlier Genes\n")
cat("=========================================\n\n")

# Read the outlier gene-to-GO mapping
cat("Loading outlier gene-to-GO mapping...\n")
outlier_genes <- read.delim(gene2go_file, stringsAsFactors = FALSE)
cat(sprintf("  Loaded %d genes with GO annotations\n\n", nrow(outlier_genes)))

# Read the full Tribolium gene2go database
cat("Loading Tribolium GO database...\n")
tribolium_gene2go <- read.delim(file.path(go_data_dir, "tribolium_gene2go.tsv"),
                                 stringsAsFactors = FALSE)
cat(sprintf("  Loaded %d GO annotations for %d Tribolium genes\n\n",
            nrow(tribolium_gene2go), length(unique(tribolium_gene2go$GeneID))))

# Create TERM2GENE mapping (GO term -> Gene)
term2gene <- tribolium_gene2go %>%
  select(GO_ID, GeneID) %>%
  distinct() %>%
  rename(term = GO_ID, gene = GeneID)

# Create TERM2NAME mapping (GO term -> Description)
term2name <- tribolium_gene2go %>%
  select(GO_ID, GO_term) %>%
  distinct() %>%
  rename(term = GO_ID, name = GO_term)

cat(sprintf("Created GO database with %d unique terms and %d gene-term associations\n\n",
            length(unique(term2gene$term)), nrow(term2gene)))

# Get the Tribolium GeneIDs from our outliers
outlier_tribolium_genes <- unique(outlier_genes$Tribolium_GeneID)
cat(sprintf("Outlier genes mapped to Tribolium: %d\n\n", length(outlier_tribolium_genes)))

# Background: all Tribolium genes with GO annotations
background_genes <- unique(tribolium_gene2go$GeneID)
cat(sprintf("Background universe: %d Tribolium genes\n\n", length(background_genes)))

# Perform enrichment analysis
cat("Performing GO enrichment analysis...\n")

# All GO terms (no category filter)
enrich_all <- enricher(
  gene = outlier_tribolium_genes,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  universe = background_genes,
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2,
  minGSSize = 2,
  maxGSSize = 500
)

cat("\n=== GO Enrichment Results ===\n")
if (!is.null(enrich_all) && nrow(as.data.frame(enrich_all)) > 0) {
  results_df <- as.data.frame(enrich_all)
  cat(sprintf("Significantly enriched GO terms: %d\n\n", nrow(results_df)))

  # Print top results
  print(head(results_df[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust")], 20))

  # Save results
  output_file <- file.path(results_dir, "go_enrichment_results.tsv")
  write.table(results_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("\nResults saved to: %s\n", output_file))

  # Create visualizations if we have results
  if (nrow(results_df) >= 2) {
    cat("\nCreating visualizations...\n")

    # Dot plot
    p1 <- dotplot(enrich_all, showCategory = min(20, nrow(results_df))) +
      ggtitle("GO Enrichment - Outlier Genes") +
      theme(axis.text.y = element_text(size = 8))

    ggsave(file.path(results_dir, "go_enrichment_dotplot.png"), p1,
           width = 10, height = 8, dpi = 300)
    ggsave(file.path(results_dir, "go_enrichment_dotplot.pdf"), p1,
           width = 10, height = 8)
    cat("  Saved: go_enrichment_dotplot.png/pdf\n")

    # Bar plot
    p2 <- barplot(enrich_all, showCategory = min(15, nrow(results_df))) +
      ggtitle("GO Enrichment - Top Terms") +
      theme(axis.text.y = element_text(size = 9))

    ggsave(file.path(results_dir, "go_enrichment_barplot.png"), p2,
           width = 10, height = 7, dpi = 300)
    ggsave(file.path(results_dir, "go_enrichment_barplot.pdf"), p2,
           width = 10, height = 7)
    cat("  Saved: go_enrichment_barplot.png/pdf\n")

    # Enrichment map (if enough terms and pairwise similarities can be computed)
    if (nrow(results_df) >= 3) {
      tryCatch({
        enrich_all_pt <- pairwise_termsim(enrich_all)
        p3 <- emapplot(enrich_all_pt, showCategory = min(30, nrow(results_df))) +
          ggtitle("GO Term Network")

        ggsave(file.path(results_dir, "go_enrichment_network.png"), p3,
               width = 12, height = 10, dpi = 300)
        cat("  Saved: go_enrichment_network.png\n")
      }, error = function(e) {
        cat("  Note: Could not create network plot (not enough term similarity)\n")
      })
    }
  }

} else {
  cat("No significantly enriched GO terms found.\n")
  cat("This could be due to:\n")
  cat("  - Small sample size of outlier genes\n")
  cat("  - Outliers spread across many functional categories\n")
  cat("  - Strict p-value/q-value cutoffs\n")

  # Try with more relaxed cutoffs
  cat("\nTrying with relaxed cutoffs (p < 0.25)...\n")
  enrich_relaxed <- enricher(
    gene = outlier_tribolium_genes,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    universe = background_genes,
    pvalueCutoff = 0.25,
    qvalueCutoff = 0.5,
    minGSSize = 1,
    maxGSSize = 1000
  )

  if (!is.null(enrich_relaxed) && nrow(as.data.frame(enrich_relaxed)) > 0) {
    results_relaxed <- as.data.frame(enrich_relaxed)
    cat(sprintf("Found %d terms with relaxed cutoffs:\n\n", nrow(results_relaxed)))
    print(head(results_relaxed[, c("ID", "Description", "GeneRatio", "pvalue")], 15))

    output_file <- file.path(results_dir, "go_enrichment_results_relaxed.tsv")
    write.table(results_relaxed, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("\nRelaxed results saved to: %s\n", output_file))
  }
}

# Summarize by GO category
cat("\n=== GO Categories Summary ===\n")
category_summary <- tribolium_gene2go %>%
  filter(GeneID %in% outlier_tribolium_genes) %>%
  group_by(Category) %>%
  summarise(
    N_Terms = n_distinct(GO_ID),
    N_Genes = n_distinct(GeneID),
    .groups = "drop"
  )
print(category_summary)

# Create a summary of functional categories
cat("\n=== Functional Summary of Outlier Genes ===\n")
outlier_functions <- tribolium_gene2go %>%
  filter(GeneID %in% outlier_tribolium_genes) %>%
  select(GeneID, GO_term, Category) %>%
  distinct()

cat("\nBiological Process terms:\n")
bp_terms <- outlier_functions %>% filter(Category == "Process")
if (nrow(bp_terms) > 0) {
  print(head(bp_terms$GO_term, 15))
}

cat("\nMolecular Function terms:\n")
mf_terms <- outlier_functions %>% filter(Category == "Function")
if (nrow(mf_terms) > 0) {
  print(head(mf_terms$GO_term, 15))
}

cat("\nCellular Component terms:\n")
cc_terms <- outlier_functions %>% filter(Category == "Component")
if (nrow(cc_terms) > 0) {
  print(head(cc_terms$GO_term, 10))
}

cat("\n\nDone!\n")
