# title: snp_to_gene_mapping.R
# project: BIOL624 Final Project
# author: Reina Hastings (reinahastings13@gmail.com)
# date created: 2026-04-13
# last modified: 2026-04-13
#
# purpose:
#   Map outlier SNPs to genes using GFF annotation with asymmetric flanking
#   (5 kb upstream, 2 kb downstream). Build foreground gene list (high-confidence
#   outliers, n_methods >= 2) and background gene list (all genes with at least
#   one autosomal SNP). Produces snp-to-gene mapping, foreground/background gene
#   tables, and a summary report. Implements Task 4.1 of the enrichment pipeline.
#
# inputs:
#   - results/cross_method_concordance/master_candidate_table.tsv (54,965 outlier SNPs)
#   - data/plink/ladybug_snps.bim (803,821 SNPs, full universe)
#   - data/coccinella_septempunctata_/ncbi_dataset/data/GCF_907165205.1/genomic.gff
#
# outputs (all in results/enrichment/):
#   - foreground_genes.tsv
#   - background_genes.tsv
#   - snp_to_gene_mapping.tsv
#   - snp_gene_summary.txt
#
# usage:
#   Rscript scripts/snp_to_gene_mapping.R

library(data.table)
library(dplyr)
library(stringr)

# --- Configuration ---

BASE_DIR   <- "."
GFF_PATH   <- file.path(BASE_DIR, "data/coccinella_septempunctata_/ncbi_dataset/data/GCF_907165205.1/genomic.gff")
MASTER_TSV <- file.path(BASE_DIR, "results/cross_method_concordance/master_candidate_table.tsv")
BIM_PATH   <- file.path(BASE_DIR, "data/plink/ladybug_snps.bim")
OUT_DIR    <- file.path(BASE_DIR, "results/enrichment")

# NC_058198.1 is chromosome 10 (sex); the other 9 are autosomes.
# NOTE (2026-04-23): previously set to NC_058189.1 (a typo). Corrected as part of Task 4.7.
SEX_CHROM       <- "NC_058198.1"
UPSTREAM_FLANK  <- 5000   # 5 kb upstream
DOWNSTREAM_FLANK <- 2000  # 2 kb downstream
MIN_METHODS     <- 2      # foreground threshold

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== SNP-to-Gene Mapping (Task 4.1) ===\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Parse GFF Gene Features ---

cat("--- Parsing GFF gene features ---\n")

gff_lines <- fread(
  cmd = paste("grep -v '^#'", shQuote(GFF_PATH)),
  sep = "\t",
  header = FALSE,
  select = c(1, 3, 4, 5, 7, 9),
  col.names = c("chrom", "feature", "start", "end", "strand", "attributes")
)

# Keep only 'gene' features on autosomes
genes <- gff_lines[feature == "gene" & chrom != SEX_CHROM]
rm(gff_lines)

cat("  Total gene features (autosomes):", nrow(genes), "\n")

# Extract gene_id (LOC*) and gene_name from attributes
genes[, gene_id := str_extract(attributes, "(?<=gene=)LOC[0-9]+")]
genes[, gene_name := str_extract(attributes, "(?<=Name=)[^;]+")]

# Some genes may be non-protein-coding (e.g., lncRNA, tRNA) with different ID patterns
# Also extract GeneID numeric for gene2go lookup later
genes[, ncbi_gene_id := str_extract(attributes, "(?<=GeneID:)[0-9]+")]

# For genes without LOC IDs, use gene-<ID> pattern from attributes
genes[is.na(gene_id), gene_id := str_extract(attributes, "(?<=ID=gene-)[^;]+")]

# Drop attributes column, no longer needed
genes[, attributes := NULL]

cat("  Genes with LOC IDs:", sum(str_detect(genes$gene_id, "^LOC"), na.rm = TRUE), "\n")
cat("  Genes without gene_id (will be dropped):", sum(is.na(genes$gene_id)), "\n")

genes <- genes[!is.na(gene_id)]

cat("  Retained genes:", nrow(genes), "\n\n")

# --- Build Flanked Gene Intervals ---
# Asymmetric flanking: 5 kb upstream, 2 kb downstream (strand-aware)

cat("--- Building flanked gene intervals ---\n")

genes[, `:=`(
  flank_start = ifelse(strand == "+",
                       pmax(1L, start - UPSTREAM_FLANK),
                       pmax(1L, start - DOWNSTREAM_FLANK)),
  flank_end   = ifelse(strand == "+",
                       end + DOWNSTREAM_FLANK,
                       end + UPSTREAM_FLANK)
)]

cat("  Flanking: ", UPSTREAM_FLANK/1000, "kb upstream, ",
    DOWNSTREAM_FLANK/1000, "kb downstream\n")
cat("  Median gene body size:", median(genes$end - genes$start + 1), "bp\n")
cat("  Median flanked interval:", median(genes$flank_end - genes$flank_start + 1), "bp\n\n")

# --- Load SNP Data ---

cat("--- Loading SNP data ---\n")

# Full SNP universe from BIM (autosomal only)
bim <- fread(BIM_PATH, header = FALSE,
             select = c(1, 2, 4),
             col.names = c("chrom", "snp_id", "pos"))
bim <- bim[chrom != SEX_CHROM]
cat("  Autosomal SNPs in BIM:", nrow(bim), "\n")

# Master candidate table (outlier SNPs)
master <- fread(MASTER_TSV)
cat("  Total outlier candidates:", nrow(master), "\n")

# Foreground: n_methods >= 2
foreground_snps <- master[n_methods >= MIN_METHODS]
cat("  High-confidence outliers (n_methods >= ", MIN_METHODS, "):", nrow(foreground_snps), "\n\n")

# --- SNP-to-Gene Overlap via data.table Non-Equi Join ---
# For each SNP, find all genes whose flanked interval contains it

cat("--- Computing SNP-to-gene overlaps ---\n")

# Function to map SNPs to genes using non-equi join
map_snps_to_genes <- function(snp_dt, gene_dt, label = "") {
  # Set keys for the join
  # Join condition: same chrom, snp pos >= flank_start, snp pos <= flank_end
  result <- gene_dt[snp_dt,
    on = .(chrom = chrom, flank_start <= pos, flank_end >= pos),
    .(snp_chrom = i.chrom, snp_pos = i.pos, snp_id = i.snp_id,
      gene_id = x.gene_id, gene_name = x.gene_name,
      gene_start = x.start, gene_end = x.end, gene_strand = x.strand,
      ncbi_gene_id = x.ncbi_gene_id),
    allow.cartesian = TRUE,
    nomatch = NA
  ]

  # Classify location: genic, upstream, or downstream
  result[, location_type := fcase(
    is.na(gene_id), "intergenic",
    snp_pos >= gene_start & snp_pos <= gene_end, "genic",
    gene_strand == "+" & snp_pos < gene_start, "upstream",
    gene_strand == "+" & snp_pos > gene_end, "downstream",
    gene_strand == "-" & snp_pos > gene_end, "upstream",
    gene_strand == "-" & snp_pos < gene_start, "downstream",
    default = "flanking"
  )]

  # Distance to gene body (0 if genic)
  result[, distance_to_gene := fcase(
    location_type == "intergenic", NA_integer_,
    location_type == "genic", 0L,
    snp_pos < gene_start, gene_start - snp_pos,
    snp_pos > gene_end, snp_pos - gene_end,
    default = 0L
  )]

  if (nchar(label) > 0) {
    n_mapped <- sum(!is.na(result$gene_id))
    n_total  <- uniqueN(result$snp_id)
    n_with_gene <- uniqueN(result[!is.na(gene_id)]$snp_id)
    cat("  ", label, ": ", n_with_gene, "/", n_total,
        " SNPs map to a gene (", round(100 * n_with_gene / n_total, 1), "%)\n", sep = "")
  }

  return(result)
}

# Map all autosomal SNPs (for background)
cat("  Mapping full autosomal SNP universe to genes...\n")
bg_mapping <- map_snps_to_genes(bim, genes, label = "Background")

# Map foreground SNPs
cat("  Mapping foreground outlier SNPs to genes...\n")
fg_mapping <- map_snps_to_genes(foreground_snps[, .(chrom, pos, snp_id)], genes, label = "Foreground")

# --- Build Foreground Gene Table ---

cat("\n--- Building foreground gene table ---\n")

# Get the consensus_pattern and detecting methods for each foreground SNP
fg_snp_info <- foreground_snps[, .(snp_id, n_methods, consensus_pattern,
                                    method_families)]

# Join pattern info onto foreground mapping
fg_mapped <- fg_mapping[!is.na(gene_id)]
fg_mapped <- merge(fg_mapped, fg_snp_info, by = "snp_id", all.x = TRUE)

# Per-gene summary for foreground
fg_genes <- fg_mapped[, .(
  chrom       = first(snp_chrom),
  gene_start  = first(gene_start),
  gene_end    = first(gene_end),
  strand      = first(gene_strand),
  gene_name   = first(gene_name),
  ncbi_gene_id = first(ncbi_gene_id),
  n_outlier_snps = .N,
  # Consensus pattern: majority rule, ties -> BOTH
  detecting_methods = paste(unique(na.omit(unlist(strsplit(method_families, ",")))), collapse = ","),
  location_types    = paste(unique(location_type), collapse = ","),
  patterns          = list(table(consensus_pattern))
), by = gene_id]

# Assign consensus pattern per gene (majority rule; tie -> BOTH)
fg_genes[, consensus_pattern := {
  sapply(patterns, function(tbl) {
    if (length(tbl) == 0) return(NA_character_)
    if (length(tbl) == 1) return(names(tbl))
    sorted <- sort(tbl, decreasing = TRUE)
    if (sorted[1] > sorted[2]) {
      return(names(sorted)[1])
    } else {
      return("BOTH")
    }
  })
}]
fg_genes[, patterns := NULL]

# Determine primary location type
fg_genes[, location_type := fifelse(
  str_detect(location_types, "genic"), "genic",
  fifelse(str_detect(location_types, "upstream"), "upstream", "downstream")
)]
fg_genes[, location_types := NULL]

cat("  Foreground genes:", nrow(fg_genes), "\n")
cat("  Pattern distribution:\n")
print(table(fg_genes$consensus_pattern))
cat("\n")

# --- Build Background Gene Table ---

cat("--- Building background gene table ---\n")

bg_mapped_genes <- bg_mapping[!is.na(gene_id)]

bg_genes <- bg_mapped_genes[, .(
  chrom       = first(snp_chrom),
  gene_start  = first(gene_start),
  gene_end    = first(gene_end),
  strand      = first(gene_strand),
  gene_name   = first(gene_name),
  ncbi_gene_id = first(ncbi_gene_id),
  n_total_snps = uniqueN(snp_id)
), by = gene_id]

cat("  Background genes (testable universe):", nrow(bg_genes), "\n")

# Verify foreground is subset of background
fg_in_bg <- sum(fg_genes$gene_id %in% bg_genes$gene_id)
cat("  Foreground genes in background:", fg_in_bg, "/", nrow(fg_genes), "\n")
if (fg_in_bg < nrow(fg_genes)) {
  warning("Some foreground genes are NOT in background! This should not happen.")
}

cat("\n")

# --- Build SNP-to-Gene Mapping Table (Foreground Only) ---

cat("--- Building SNP-to-gene mapping table ---\n")

snp_gene_map <- fg_mapping[, .(
  snp_id, chrom = snp_chrom, pos = snp_pos,
  gene_id, location_type, distance_to_gene
)]

cat("  Total mapping rows (including intergenic):", nrow(snp_gene_map), "\n")
cat("  Mapped to gene:", sum(!is.na(snp_gene_map$gene_id)), "\n")
cat("  Intergenic:", sum(is.na(snp_gene_map$gene_id)), "\n\n")

# --- Write Outputs ---

cat("--- Writing output files ---\n")

# Foreground genes
fg_out <- fg_genes[, .(gene_id, chrom, gene_start, gene_end, strand, gene_name,
                        ncbi_gene_id, n_outlier_snps, consensus_pattern,
                        detecting_methods, location_type)]
fwrite(fg_out, file.path(OUT_DIR, "foreground_genes.tsv"), sep = "\t")
cat("  foreground_genes.tsv:", nrow(fg_out), "genes\n")

# Background genes
bg_out <- bg_genes[, .(gene_id, chrom, gene_start, gene_end, strand, gene_name,
                        ncbi_gene_id, n_total_snps)]
fwrite(bg_out, file.path(OUT_DIR, "background_genes.tsv"), sep = "\t")
cat("  background_genes.tsv:", nrow(bg_out), "genes\n")

# SNP-to-gene mapping
fwrite(snp_gene_map, file.path(OUT_DIR, "snp_to_gene_mapping.tsv"), sep = "\t")
cat("  snp_to_gene_mapping.tsv:", nrow(snp_gene_map), "rows\n")

# --- Summary Report ---

cat("  Writing snp_gene_summary.txt...\n\n")

n_fg_snps       <- nrow(foreground_snps)
n_fg_mapped     <- uniqueN(fg_mapping[!is.na(gene_id)]$snp_id)
n_fg_intergenic <- n_fg_snps - n_fg_mapped
pct_mapped      <- round(100 * n_fg_mapped / n_fg_snps, 1)
pct_intergenic  <- round(100 * n_fg_intergenic / n_fg_snps, 1)

n_bg_snps       <- nrow(bim)
n_bg_mapped     <- uniqueN(bg_mapping[!is.na(gene_id)]$snp_id)
pct_bg_mapped   <- round(100 * n_bg_mapped / n_bg_snps, 1)

# Location type breakdown for foreground SNPs that map to genes
loc_table <- fg_mapping[!is.na(gene_id), .N, by = location_type]

summary_text <- paste0(
  "========================================================================\n",
  "SNP-TO-GENE MAPPING SUMMARY (Task 4.1)\n",
  "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "========================================================================\n\n",

  "--- Configuration ---\n",
  "Upstream flanking: ", UPSTREAM_FLANK/1000, " kb\n",
  "Downstream flanking: ", DOWNSTREAM_FLANK/1000, " kb\n",
  "Foreground threshold: n_methods >= ", MIN_METHODS, "\n",
  "Sex chromosome excluded: ", SEX_CHROM, "\n\n",

  "--- Gene Annotation ---\n",
  "Total genes in GFF (autosomes): ", nrow(genes), "\n",
  "Genes retained (with valid ID): ", nrow(genes), "\n",
  "Median gene body size: ", median(genes$end - genes$start + 1), " bp\n",
  "Median flanked interval: ", median(genes$flank_end - genes$flank_start + 1), " bp\n\n",

  "--- Foreground (High-Confidence Outlier SNPs) ---\n",
  "Total high-confidence outlier SNPs: ", n_fg_snps, "\n",
  "SNPs mapping to a gene (genic + flanking): ", n_fg_mapped, " (", pct_mapped, "%)\n",
  "SNPs intergenic: ", n_fg_intergenic, " (", pct_intergenic, "%)\n",
  "Unique foreground genes: ", nrow(fg_genes), "\n\n",

  "Location type breakdown (foreground SNP-gene hits):\n",
  paste(capture.output(print(loc_table, row.names = FALSE)), collapse = "\n"), "\n\n",

  "Foreground gene pattern distribution:\n",
  paste(capture.output(print(table(fg_genes$consensus_pattern))), collapse = "\n"), "\n\n",

  "--- Background (Testable Gene Universe) ---\n",
  "Total autosomal SNPs: ", n_bg_snps, "\n",
  "SNPs mapping to a gene: ", n_bg_mapped, " (", pct_bg_mapped, "%)\n",
  "Unique background genes: ", nrow(bg_genes), "\n",
  "Foreground genes in background: ", fg_in_bg, "/", nrow(fg_genes), "\n\n",

  "--- QC Flags ---\n",
  "Foreground gene count: ", nrow(fg_genes),
  ifelse(nrow(fg_genes) < 100, " ** RED FLAG: < 100 genes, consider expanding to full outlier set **",
         " (OK, above 100 threshold)"), "\n",
  "Intergenic rate: ", pct_intergenic, "%",
  ifelse(pct_intergenic > 80, " ** RED FLAG: > 80% intergenic, flanking may be too small **",
         " (OK, below 80% threshold)"), "\n",
  "Foreground genes > 2000: ",
  ifelse(nrow(fg_genes) > 2000, "YES -- flanking may be too wide", "No (OK)"), "\n\n",

  "========================================================================\n",
  "END OF SUMMARY\n",
  "========================================================================\n"
)

writeLines(summary_text, file.path(OUT_DIR, "snp_gene_summary.txt"))

cat(summary_text)
cat("\nDone. All outputs written to:", OUT_DIR, "\n")
