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
#   Benjamini-Hochberg FDR, reduces redundancy via rrvgo (REVIGO-like local
#   semantic similarity clustering, using org.Dm.eg.db as phylogenetic proxy
#   for information content since no Coleoptera OrgDb is available), and emits
#   publication-quality dotplots + rrvgo treemap/scatter + a cross-pattern
#   comparison figure. Implements Task 4.3 of the enrichment pipeline.
#
# inputs:
#   - results/enrichment/foreground_genes.tsv        (from Task 4.1)
#   - results/enrichment/background_genes.tsv        (from Task 4.1)
#   - results/enrichment/gene2go_foreground.tsv      (from Task 4.2)
#   - results/enrichment/gene2go_background.tsv      (from Task 4.2)
#
# outputs (all in results/enrichment/):
#   - enrichment_results_{pattern}_{ontology}.tsv
#   - enrichment_results_simplified_{pattern}_{ontology}.tsv   (rrvgo parent terms)
#   - rrvgo_sim_matrix_{pattern}_{ontology}.rds                (cached sim matrix)
#   - enrichment_summary.txt
#   - enrichment_log.txt
#   - plots/dotplot_{pattern}_{ontology}.png
#   - plots/barplot_top_terms_{pattern}.png
#   - plots/comparison_dotplot.png
#   - plots/rrvgo_treemap_{pattern}_{ontology}.png
#   - plots/rrvgo_scatter_{pattern}_{ontology}.png
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
  library(rrvgo)
  library(org.Dm.eg.db)  # phylogenetic proxy for IC (no Coleoptera OrgDb available)
})

# --- rrvgo configuration ---
# org.Dm.eg.db provides information-content weights for GO terms. D. melanogaster
# is the closest well-annotated insect to C. septempunctata on Bioconductor; T.
# castaneum has no OrgDb. IC values differ across species (D. melanogaster has
# deeper functional annotation than a non-model beetle), which is a known
# limitation and must be acknowledged in the manuscript.
RRVGO_ORGDB    <- "org.Dm.eg.db"
RRVGO_SIMCUT   <- 0.7   # "medium" in REVIGO terms: 0.4=tiny, 0.5=small, 0.7=medium, 0.9=large
RRVGO_METHOD   <- "Rel" # Relevance: Schlicker et al. similarity; works well with D. melanogaster IC

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

# --- rrvgo-based simplification ---
#
# rrvgo computes pairwise semantic similarity among significant GO terms via
# GOSemSim (using org.Dm.eg.db IC), then clusters terms and selects a parent
# term per cluster (highest-scoring member). Produces a reduced term table and
# two visualizations: treemap (parent-child hierarchy with cluster sizes) and
# MDS scatterplot (semantic space projection).

simplify_by_rrvgo <- function(enrich_df, ontology, patt, plot_dir, out_dir) {
  # Returns a list with: reduced (data.frame of parent terms), sim_matrix,
  # reduced_terms (full rrvgo reduceSimMatrix output), or NULL if <2 terms.
  if (nrow(enrich_df) < 2) {
    log_msg("rrvgo: ", patt, " x ", ontology, " has <2 terms; skipping simplification")
    return(list(reduced = enrich_df, sim_matrix = NULL, reduced_terms = NULL))
  }

  id_col    <- if ("go_id" %in% names(enrich_df)) "go_id" else "ID"
  score_col <- if ("q_value" %in% names(enrich_df)) "q_value" else "qvalue"
  go_ids    <- enrich_df[[id_col]]
  # rrvgo scores: higher = more important. Use -log10(q); clamp to avoid Inf.
  scores    <- -log10(pmax(enrich_df[[score_col]], 1e-300))
  names(scores) <- go_ids

  sim_mat <- tryCatch(
    rrvgo::calculateSimMatrix(
      go_ids,
      orgdb  = RRVGO_ORGDB,
      ont    = ontology,
      method = RRVGO_METHOD
    ),
    error = function(e) {
      log_msg("rrvgo calculateSimMatrix error (", patt, " x ", ontology, "): ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(sim_mat) || any(dim(sim_mat) < 2)) {
    log_msg("rrvgo: ", patt, " x ", ontology,
            " produced a degenerate similarity matrix; skipping")
    return(list(reduced = enrich_df, sim_matrix = sim_mat, reduced_terms = NULL))
  }

  saveRDS(sim_mat, file.path(out_dir,
          sprintf("rrvgo_sim_matrix_%s_%s.rds", patt, ontology)))

  reduced <- tryCatch(
    rrvgo::reduceSimMatrix(
      sim_mat,
      scores    = scores[rownames(sim_mat)],
      threshold = RRVGO_SIMCUT,
      orgdb     = RRVGO_ORGDB
    ),
    error = function(e) {
      log_msg("rrvgo reduceSimMatrix error (", patt, " x ", ontology, "): ",
              conditionMessage(e)); NULL
    }
  )
  if (is.null(reduced)) {
    return(list(reduced = enrich_df, sim_matrix = sim_mat, reduced_terms = NULL))
  }

  # Parent-term table: one row per unique parent (cluster representative).
  parent_df <- reduced %>%
    dplyr::distinct(parent, parentTerm) %>%
    dplyr::rename(parent_go_id = parent, parent_term = parentTerm) %>%
    dplyr::left_join(
      reduced %>%
        dplyr::group_by(parent) %>%
        dplyr::summarise(cluster_size = dplyr::n(),
                         cluster_score = max(score, na.rm = TRUE),
                         cluster_members = paste(unique(go), collapse = ";")),
      by = c("parent_go_id" = "parent")
    ) %>%
    dplyr::arrange(dplyr::desc(cluster_score))

  # Treemap
  tm_path <- file.path(plot_dir, sprintf("rrvgo_treemap_%s_%s.png", patt, ontology))
  tryCatch({
    png(tm_path, width = 1600, height = 1000, res = 150)
    rrvgo::treemapPlot(reduced,
                       title = sprintf("rrvgo treemap: %s x %s", patt, ontology))
    dev.off()
  }, error = function(e) {
    log_msg("rrvgo treemap failed (", patt, " x ", ontology, "): ", conditionMessage(e))
    if (dev.cur() > 1) dev.off()
  })

  # MDS scatter (needs >=3 terms to be meaningful)
  if (nrow(sim_mat) >= 3) {
    sc_path <- file.path(plot_dir, sprintf("rrvgo_scatter_%s_%s.png", patt, ontology))
    tryCatch({
      p_sc <- rrvgo::scatterPlot(sim_mat, reduced, size = "score",
                                 addLabel = TRUE, labelSize = 3)
      ggsave(sc_path, p_sc, width = 9, height = 7, dpi = 300)
    }, error = function(e) {
      log_msg("rrvgo scatter failed (", patt, " x ", ontology, "): ",
              conditionMessage(e))
    })
  }

  list(reduced = parent_df, sim_matrix = sim_mat, reduced_terms = reduced)
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

    res_df <- as.data.frame(er)
    res_df <- dplyr::rename(res_df,
      go_id      = "ID",
      go_term    = "Description",
      gene_ratio = "GeneRatio",
      bg_ratio   = "BgRatio",
      p_value    = "pvalue",
      q_value    = "qvalue",
      p_adjust   = "p.adjust",
      gene_list  = "geneID",
      count      = "Count"
    )
    res_df <- dplyr::arrange(res_df, q_value)

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

    # Simplify via rrvgo on significant terms only (fall back to suggestive if none sig).
    base_for_simp <- if (nrow(sig) > 0) sig else sig_suggestive
    if (nrow(base_for_simp) > 0) {
      rr <- simplify_by_rrvgo(base_for_simp, ont, patt, PLOT_DIR, ENR_DIR)
      simp_df <- rr$reduced
      log_msg(sprintf("rrvgo %s x %s: %d input terms -> %d parent clusters",
                      patt, ont, nrow(base_for_simp),
                      if (is.data.frame(simp_df)) nrow(simp_df) else NA_integer_))
    } else {
      simp_df <- base_for_simp
    }
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
