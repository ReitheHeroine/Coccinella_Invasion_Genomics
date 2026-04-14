# title: go_enrichment_ora.R
# project: BIOL624 Final Project
# author: Reina Hastings (reinahastings13@gmail.com)
# date created: 2026-04-14
# last modified: 2026-04-14
#
# purpose:
#   GO over-representation analysis (ORA) stratified by consensus_pattern.
#   Runs clusterProfiler::enricher() on outlier genes vs the testable background
#   (background genes with >=1 GO term). Tests INVASION_SPECIFIC, POST_INVASION,
#   BOTH, and ALL foreground sets against each ontology (BP, MF, CC). Applies
#   Benjamini-Hochberg FDR, reduces redundancy via ancestor pruning using GO.db,
#   and emits publication-quality dotplots + a cross-pattern comparison figure.
#   Implements Task 4.3 of the enrichment pipeline.
#
# inputs:
#   - results/enrichment/foreground_genes.tsv        (from Task 4.1)
#   - results/enrichment/background_genes.tsv        (from Task 4.1)
#   - results/enrichment/gene2go_foreground.tsv      (from Task 4.2)
#   - results/enrichment/gene2go_background.tsv      (from Task 4.2)
#
# outputs (all in results/enrichment/):
#   - enrichment_results_{pattern}_{ontology}.tsv
#   - enrichment_results_simplified_{pattern}_{ontology}.tsv
#   - enrichment_summary.txt
#   - enrichment_log.txt
#   - plots/dotplot_{pattern}_{ontology}.png
#   - plots/barplot_top_terms_{pattern}.png
#   - plots/comparison_dotplot.png
#
# usage:
#   Rscript scripts/go_enrichment_ora.R
#   Rscript scripts/go_enrichment_ora.R --qcut 0.05 --top 20

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(GO.db)
})

# --- CLI (simple key=value parser to avoid optparse dependency) ---

.defaults <- list(
  qcut         = 0.05,
  top          = 20L,
  min_genes    = 20L,
  min_gs_size  = 5L,
  max_gs_size  = 500L
)

.parse_args <- function(args, defaults) {
  opt <- defaults
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    key <- sub("^--", "", a)
    key <- gsub("-", "_", key)
    if (i + 1 > length(args) || startsWith(args[i + 1], "--")) {
      stop("Argument ", a, " requires a value (usage: --qcut 0.05 --top 20 --min-genes 20 --min-gs-size 5 --max-gs-size 500)")
    }
    val <- args[i + 1]
    if (!key %in% names(defaults)) stop("Unknown argument: ", a)
    opt[[key]] <- if (is.integer(defaults[[key]])) as.integer(val) else as.numeric(val)
    i <- i + 2
  }
  opt
}
opt <- .parse_args(commandArgs(trailingOnly = TRUE), .defaults)

# --- Configuration ---

BASE_DIR <- "."
ENR_DIR  <- file.path(BASE_DIR, "results/enrichment")
PLOT_DIR <- file.path(ENR_DIR, "plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

FG_GENES_PATH <- file.path(ENR_DIR, "foreground_genes.tsv")
BG_GENES_PATH <- file.path(ENR_DIR, "background_genes.tsv")
FG_GO_PATH    <- file.path(ENR_DIR, "gene2go_foreground.tsv")
BG_GO_PATH    <- file.path(ENR_DIR, "gene2go_background.tsv")

PATTERNS   <- c("ALL", "INVASION_SPECIFIC", "POST_INVASION", "BOTH")
ONTOLOGIES <- c("BP", "MF", "CC")

# Broad GO terms to flag if they dominate top results (potential background bias).
BROAD_TERMS <- c("metabolic process", "binding", "cellular process",
                 "biological process", "molecular function", "cellular component",
                 "organic substance metabolic process", "primary metabolic process",
                 "protein binding", "catalytic activity")

# --- Log helpers ---

LOG_PATH <- file.path(ENR_DIR, "enrichment_log.txt")
if (file.exists(LOG_PATH)) file.remove(LOG_PATH)

log_msg <- function(...) {
  msg <- paste0(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                paste0(..., collapse = ""))
  cat(msg, "\n", sep = "")
  cat(msg, "\n", file = LOG_PATH, append = TRUE, sep = "")
}

log_msg("=== GO Enrichment (ORA, Task 4.3) ===")
log_msg("qcut=", opt$qcut, " top=", opt$top, " min_genes=", opt$min_genes,
        " minGSSize=", opt$min_gs_size, " maxGSSize=", opt$max_gs_size)

# --- Input validation ---

required_inputs <- c(FG_GENES_PATH, BG_GENES_PATH, FG_GO_PATH, BG_GO_PATH)
missing_inputs  <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  log_msg("BLOCKER: missing required inputs:")
  for (f in missing_inputs) log_msg("  - ", f)
  stop("Task 4.2 outputs not found. Cannot run Task 4.3 until gene-to-GO mapping is complete.")
}

# --- Load data ---

fg_genes <- fread(FG_GENES_PATH)
bg_genes <- fread(BG_GENES_PATH)
fg_go    <- fread(FG_GO_PATH)
bg_go    <- fread(BG_GO_PATH)

log_msg("Loaded foreground genes: ",    nrow(fg_genes))
log_msg("Loaded background genes: ",    nrow(bg_genes))
log_msg("Loaded foreground gene-GO: ",  nrow(fg_go),  " rows, ", length(unique(fg_go$gene_id)), " genes")
log_msg("Loaded background gene-GO: ",  nrow(bg_go),  " rows, ", length(unique(bg_go$gene_id)), " genes")

# --- DISCORDANT decision ---

# DISCORDANT signals methodological disagreement across detection methods. Exclude
# from pattern-stratified runs (the pattern has no coherent biological meaning),
# but retain in ALL since ALL is the entire high-confidence foreground.
n_discordant <- sum(fg_genes$consensus_pattern == "DISCORDANT")
log_msg("DISCORDANT genes in foreground: ", n_discordant,
        " (retained in ALL; excluded from pattern-stratified runs)")

# --- TERM2GENE / TERM2NAME per ontology ---

# Universe = background genes with >=1 GO term. Filter enricher universe to this set.
testable_bg <- unique(bg_go$gene_id)
log_msg("Testable background (>=1 GO term): ", length(testable_bg),
        " of ", nrow(bg_genes), " (",
        sprintf("%.1f%%", 100 * length(testable_bg) / nrow(bg_genes)), ")")

# Build per-ontology TERM2GENE and TERM2NAME from bg_go (contains both fg and bg annotations
# as long as Task 4.2 annotated the full background; fg is subset of bg so bg_go is sufficient).
make_term_tables <- function(go_df, ont) {
  sub <- go_df[go_ontology == ont]
  t2g <- unique(sub[, .(term = go_id, gene = gene_id)])
  t2n <- unique(sub[, .(term = go_id, name = go_term)])
  list(term2gene = as.data.frame(t2g), term2name = as.data.frame(t2n))
}

# --- Ancestor-based simplification (OrgDb-free) ---

ancestor_map <- list(
  BP = as.list(GOBPANCESTOR),
  MF = as.list(GOMFANCESTOR),
  CC = as.list(GOCCANCESTOR)
)

simplify_by_ancestor <- function(enrich_df, ontology) {
  # Keep the most specific significant term in each lineage: drop any term that is
  # an ancestor of another significant term in the same result.
  if (nrow(enrich_df) == 0) return(enrich_df)
  amap <- ancestor_map[[ontology]]
  sig_ids <- enrich_df$ID
  drop <- rep(FALSE, length(sig_ids))
  for (i in seq_along(sig_ids)) {
    others <- setdiff(sig_ids, sig_ids[i])
    for (o in others) {
      anc <- amap[[o]]
      if (!is.null(anc) && sig_ids[i] %in% anc) { drop[i] <- TRUE; break }
    }
  }
  enrich_df[!drop, , drop = FALSE]
}

# --- Foreground gene sets by pattern ---

foreground_sets <- list(
  ALL               = fg_genes$gene_id,
  INVASION_SPECIFIC = fg_genes$gene_id[fg_genes$consensus_pattern == "INVASION_SPECIFIC"],
  POST_INVASION     = fg_genes$gene_id[fg_genes$consensus_pattern == "POST_INVASION"],
  BOTH              = fg_genes$gene_id[fg_genes$consensus_pattern == "BOTH"]
)

# Annotation coverage per pattern
ann_coverage <- sapply(foreground_sets, function(g) {
  n_total <- length(g)
  n_ann   <- sum(g %in% unique(fg_go$gene_id))
  c(n_total = n_total, n_annotated = n_ann,
    pct_annotated = if (n_total == 0) NA else 100 * n_ann / n_total)
})
log_msg("Foreground pattern sizes and annotation coverage:")
print(t(ann_coverage))

# --- Enrichment loop ---

all_results  <- list()
simp_results <- list()

for (patt in PATTERNS) {
  fg_set <- foreground_sets[[patt]]
  n_ann_fg <- sum(fg_set %in% unique(fg_go$gene_id))
  if (n_ann_fg < opt$min_genes) {
    log_msg("WARN: pattern ", patt, " has only ", n_ann_fg,
            " annotated foreground genes (< ", opt$min_genes, "); flagged as underpowered")
  }

  for (ont in ONTOLOGIES) {
    tt <- make_term_tables(bg_go, ont)
    if (nrow(tt$term2gene) == 0) {
      log_msg("SKIP: ", patt, " x ", ont, " (no GO terms in background for this ontology)")
      next
    }

    er <- tryCatch(
      enricher(
        gene          = fg_set,
        universe      = testable_bg,
        TERM2GENE     = tt$term2gene,
        TERM2NAME     = tt$term2name,
        pAdjustMethod = "BH",
        pvalueCutoff  = 1,   # retain all; filter after
        qvalueCutoff  = 1,
        minGSSize     = opt$min_gs_size,
        maxGSSize     = opt$max_gs_size
      ),
      error = function(e) { log_msg("ERROR in enricher for ", patt, " x ", ont, ": ", conditionMessage(e)); NULL }
    )

    if (is.null(er) || nrow(as.data.frame(er)) == 0) {
      log_msg("No enrichment output for ", patt, " x ", ont)
      next
    }

    res_df <- as.data.frame(er) %>%
      rename(go_id = ID, go_term = Description,
             gene_ratio = GeneRatio, bg_ratio = BgRatio,
             p_value = pvalue, q_value = qvalue, p_adjust = p.adjust,
             gene_list = geneID, count = Count) %>%
      arrange(q_value)

    out_path <- file.path(ENR_DIR, sprintf("enrichment_results_%s_%s.tsv", patt, ont))
    fwrite(res_df, out_path, sep = "\t")

    sig <- res_df[res_df$q_value <= opt$qcut, , drop = FALSE]
    sig_suggestive <- res_df[res_df$q_value <= 0.10, , drop = FALSE]
    broad_hits <- intersect(tolower(head(res_df$go_term, 10)), tolower(BROAD_TERMS))

    log_msg(sprintf("%s x %s: %d total, %d sig at q<=%.2f (%d suggestive at q<=0.10)%s",
                    patt, ont, nrow(res_df), nrow(sig), opt$qcut, nrow(sig_suggestive),
                    if (length(broad_hits) > 0)
                      paste0("; BROAD TERMS in top 10: ", paste(broad_hits, collapse = ", "))
                    else ""))

    all_results[[paste(patt, ont, sep = "_")]] <- list(
      pattern = patt, ontology = ont, result = res_df,
      n_sig = nrow(sig), n_suggestive = nrow(sig_suggestive),
      broad_hits = broad_hits, enrich_obj = er,
      n_annotated_fg = n_ann_fg
    )

    # Simplify via ancestor pruning on significant terms only (fall back to suggestive if none sig).
    base_for_simp <- if (nrow(sig) > 0) sig else sig_suggestive
    simp_df <- if (nrow(base_for_simp) > 0) simplify_by_ancestor(base_for_simp, ont) else base_for_simp
    simp_path <- file.path(ENR_DIR, sprintf("enrichment_results_simplified_%s_%s.tsv", patt, ont))
    fwrite(simp_df, simp_path, sep = "\t")
    simp_results[[paste(patt, ont, sep = "_")]] <- simp_df

    # Dotplot (top N by q-value among all terms; color by q; size by count).
    plot_df <- head(res_df, opt$top)
    if (nrow(plot_df) > 0) {
      plot_df$go_term <- factor(plot_df$go_term, levels = rev(plot_df$go_term))
      # Parse GeneRatio "k/n" into numeric
      plot_df$gene_ratio_num <- sapply(strsplit(as.character(plot_df$gene_ratio), "/"),
                                       function(x) as.numeric(x[1]) / as.numeric(x[2]))
      p <- ggplot(plot_df, aes(x = gene_ratio_num, y = go_term, size = count, color = q_value)) +
        geom_point() +
        scale_color_gradient(low = "red", high = "blue") +
        labs(title = sprintf("GO %s enrichment: %s (top %d)", ont, patt, opt$top),
             x = "Gene ratio", y = NULL, size = "Count", color = "q-value") +
        theme_bw(base_size = 10)
      ggsave(file.path(PLOT_DIR, sprintf("dotplot_%s_%s.png", patt, ont)),
             p, width = 8, height = 0.28 * nrow(plot_df) + 2, dpi = 300, limitsize = FALSE)
    }
  }

  # Combined barplot across ontologies: top 10 by q_value from significant terms.
  per_patt <- all_results[grep(paste0("^", patt, "_"), names(all_results))]
  if (length(per_patt) > 0) {
    combined <- do.call(rbind, lapply(per_patt, function(x) {
      df <- x$result
      if (nrow(df) == 0) return(NULL)
      df$ontology <- x$ontology
      df
    }))
    if (!is.null(combined) && nrow(combined) > 0) {
      combined_top <- combined %>%
        group_by(ontology) %>%
        slice_min(order_by = q_value, n = 10, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(label = paste0(go_term, " (", ontology, ")"),
               neg_log_q = -log10(q_value))
      combined_top$label <- factor(combined_top$label,
                                   levels = rev(combined_top$label[order(combined_top$neg_log_q)]))
      p_bar <- ggplot(combined_top, aes(x = neg_log_q, y = label, fill = ontology)) +
        geom_col() +
        geom_vline(xintercept = -log10(opt$qcut), linetype = "dashed", color = "grey40") +
        labs(title = sprintf("Top GO terms: %s (top 10 per ontology)", patt),
             x = expression(-log[10](q)), y = NULL, fill = "Ontology") +
        theme_bw(base_size = 10)
      ggsave(file.path(PLOT_DIR, sprintf("barplot_top_terms_%s.png", patt)),
             p_bar, width = 9, height = 0.25 * nrow(combined_top) + 2, dpi = 300, limitsize = FALSE)
    }
  }
}

# --- Comparison dotplot (BP across patterns) ---

bp_sets <- all_results[grep("_BP$", names(all_results))]
if (length(bp_sets) > 0) {
  cmp <- do.call(rbind, lapply(bp_sets, function(x) {
    df <- x$result[x$result$q_value <= opt$qcut, , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df$pattern <- x$pattern
    df
  }))
  # If nothing sig, fall back to suggestive so we still contrast patterns.
  if (is.null(cmp) || nrow(cmp) == 0) {
    cmp <- do.call(rbind, lapply(bp_sets, function(x) {
      df <- x$result[x$result$q_value <= 0.10, , drop = FALSE]
      if (nrow(df) == 0) return(NULL)
      df$pattern <- x$pattern
      df
    }))
    if (!is.null(cmp) && nrow(cmp) > 0) {
      log_msg("Comparison plot: no BP terms sig at q<=", opt$qcut, "; using q<=0.10 suggestive.")
    }
  }
  if (!is.null(cmp) && nrow(cmp) > 0) {
    cmp <- cmp %>%
      group_by(pattern) %>%
      slice_min(order_by = q_value, n = 10, with_ties = FALSE) %>%
      ungroup()
    cmp$gene_ratio_num <- sapply(strsplit(as.character(cmp$gene_ratio), "/"),
                                 function(x) as.numeric(x[1]) / as.numeric(x[2]))
    cmp$pattern <- factor(cmp$pattern, levels = PATTERNS)
    cmp$go_term <- factor(cmp$go_term, levels = unique(cmp$go_term[order(cmp$q_value)]))
    p_cmp <- ggplot(cmp, aes(x = pattern, y = go_term, size = gene_ratio_num, color = q_value)) +
      geom_point() +
      scale_color_gradient(low = "red", high = "blue") +
      labs(title = "GO BP enrichment across patterns",
           x = NULL, y = NULL, size = "Gene ratio", color = "q-value") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    ggsave(file.path(PLOT_DIR, "comparison_dotplot.png"),
           p_cmp, width = 8, height = 0.28 * length(unique(cmp$go_term)) + 2,
           dpi = 300, limitsize = FALSE)
    log_msg("Comparison dotplot written with ", nrow(cmp), " pattern x term rows.")
  } else {
    log_msg("Comparison plot skipped: no BP terms at q<=0.10 across any pattern.")
  }
}

# --- Summary ---

summary_path <- file.path(ENR_DIR, "enrichment_summary.txt")
sink(summary_path)
cat("GO Enrichment Summary (Task 4.3)\n")
cat("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("q cutoff: ", opt$qcut, "\n\n")

cat("Foreground pattern sizes and annotation coverage:\n")
print(as.data.frame(t(ann_coverage)))
cat("\nDISCORDANT genes excluded from pattern-stratified runs: ", n_discordant, "\n\n")

cat("Significant terms per pattern x ontology (q <=", opt$qcut, "):\n")
tab <- data.frame(
  pattern   = character(0),
  ontology  = character(0),
  n_total   = integer(0),
  n_sig     = integer(0),
  n_suggestive_q10 = integer(0),
  n_ann_fg  = integer(0),
  broad_hits_top10 = character(0),
  stringsAsFactors = FALSE
)
for (key in names(all_results)) {
  r <- all_results[[key]]
  tab <- rbind(tab, data.frame(
    pattern          = r$pattern,
    ontology         = r$ontology,
    n_total          = nrow(r$result),
    n_sig            = r$n_sig,
    n_suggestive_q10 = r$n_suggestive,
    n_ann_fg         = r$n_annotated_fg,
    broad_hits_top10 = paste(r$broad_hits, collapse = "; "),
    stringsAsFactors = FALSE
  ))
}
print(tab, row.names = FALSE)
sink()

log_msg("Summary written: ", summary_path)
log_msg("Done.")
