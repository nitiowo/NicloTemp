# GO and KEGG over-representation analysis for all DE contrasts

library(tidyverse)
library(clusterProfiler)

base_dir <- here::here()
annot_dir <- file.path(base_dir, "00_setup/03_GO_annotation/output")
de_dir <- file.path(base_dir, "02_de_analysis/output")
out_dir <- file.path(base_dir, "05_enrichment/output/ora")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Annotations ----

btru_g2go <- readRDS(file.path(annot_dir, "step05_go_objects/btru_gene2go.Rds"))
btru_go2name <- readRDS(file.path(annot_dir, "step05_go_objects/btru_go2name.Rds"))
btru_kegg <- readRDS(file.path(annot_dir, "step06_kegg_objects/btru_kegg_term2gene.Rds"))

shae_g2go <- readRDS(file.path(annot_dir, "step05_go_objects/shae_gene2go.Rds"))
shae_go2name <- readRDS(file.path(annot_dir, "step05_go_objects/shae_go2name.Rds"))
shae_kegg <- readRDS(file.path(annot_dir, "step06_kegg_objects/shae_kegg_term2gene.Rds"))

# Contrasts
btru_contrasts <- c("temp24v16", "temp32v16", "infection", "niclosamide")
shae_contrasts <- c("temp24v16", "temp32v16", "niclosamide")

# ---- ORA Function ----
# Runs GO + KEGG ORA for a set of DE genes against a universe
run_ora <- function(de_genes, universe, gene2go, go2name, kegg_t2g,
                    species, contrast, direction, out_dir) {
  
  prefix <- paste(species, contrast, direction, sep = "_")
  
  # GO ORA (all 3 ontologies for enricher with TERM2GENE)
  go_t2g <- gene2go %>% select(GO, gene_id) %>% rename(term = GO, gene = gene_id)
  go_t2n <- go2name %>% select(GO, name) %>% rename(term = GO, name = name)
  
  go_res <- enricher(
    gene = de_genes,
    universe = universe,
    TERM2GENE = go_t2g,
    TERM2NAME = go_t2n,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 5,
    maxGSSize = 500
  )
  
  if (!is.null(go_res) && nrow(go_res@result %>% filter(p.adjust < 0.05)) > 0) {
    write_tsv(as.data.frame(go_res), file.path(out_dir, paste0(prefix, "_go.tsv")))
    
    # Dot plot (top 20)
    n_sig <- sum(go_res@result$p.adjust < 0.05)
    if (n_sig > 0) {
      pdf(file.path(out_dir, paste0(prefix, "_go_dotplot.pdf")),
          width = 10, height = max(4, min(n_sig, 20) * 0.35 + 2))
      print(dotplot(go_res, showCategory = 20) +
              ggtitle(paste(species, contrast, direction, "- GO ORA")))
      dev.off()
    }
  }
  
  # KEGG ORA
  kegg_t2g_df <- kegg_t2g %>% select(KO, gene_id) %>% rename(term = KO, gene = gene_id)
  
  kegg_res <- enricher(
    gene = de_genes,
    universe = universe,
    TERM2GENE = kegg_t2g_df,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 3,
    maxGSSize = 500
  )
  
  if (!is.null(kegg_res) && nrow(kegg_res@result %>% filter(p.adjust < 0.05)) > 0) {
    write_tsv(as.data.frame(kegg_res), file.path(out_dir, paste0(prefix, "_kegg.tsv")))
  }
  
  # Return counts
  go_n <- ifelse(!is.null(go_res), sum(go_res@result$p.adjust < 0.05), 0)
  kegg_n <- ifelse(!is.null(kegg_res), sum(kegg_res@result$p.adjust < 0.05), 0)
  
  data.frame(species = species, contrast = contrast, direction = direction,
             go_sig = go_n, kegg_sig = kegg_n)
}

# ---- Run ORA for Each Contrast ----

results <- list()

# Btru
for (ct in btru_contrasts) {
  res_file <- file.path(de_dir, "btru", paste0("btru_res_", ct, ".tsv"))
  res <- read_tsv(res_file, show_col_types = FALSE) %>%
    filter(!is.na(padj))
  
  universe <- res$gene_id
  up <- res %>% filter(padj < 0.05, log2FoldChange > 0) %>% pull(gene_id)
  down <- res %>% filter(padj < 0.05, log2FoldChange < 0) %>% pull(gene_id)
  
  cat(sprintf("Btru %s: %d up, %d down, %d universe\n", ct, length(up), length(down), length(universe)))
  
  if (length(up) >= 5) {
    results[[paste0("btru_", ct, "_up")]] <- run_ora(
      up, universe, btru_g2go, btru_go2name, btru_kegg,
      "btru", ct, "up", out_dir)
  }
  if (length(down) >= 5) {
    results[[paste0("btru_", ct, "_down")]] <- run_ora(
      down, universe, btru_g2go, btru_go2name, btru_kegg,
      "btru", ct, "down", out_dir)
  }
}

# Shae
for (ct in shae_contrasts) {
  res_file <- file.path(de_dir, "shae", paste0("shae_res_", ct, ".tsv"))
  res <- read_tsv(res_file, show_col_types = FALSE) %>%
    filter(!is.na(padj))
  
  universe <- res$gene_id
  up <- res %>% filter(padj < 0.05, log2FoldChange > 0) %>% pull(gene_id)
  down <- res %>% filter(padj < 0.05, log2FoldChange < 0) %>% pull(gene_id)
  
  cat(sprintf("Shae %s: %d up, %d down, %d universe\n", ct, length(up), length(down), length(universe)))
  
  if (length(up) >= 5) {
    results[[paste0("shae_", ct, "_up")]] <- run_ora(
      up, universe, shae_g2go, shae_go2name, shae_kegg,
      "shae", ct, "up", out_dir)
  }
  if (length(down) >= 5) {
    results[[paste0("shae_", ct, "_down")]] <- run_ora(
      down, universe, shae_g2go, shae_go2name, shae_kegg,
      "shae", ct, "down", out_dir)
  }
}

summary_df <- bind_rows(results)
write_tsv(summary_df, file.path(out_dir, "ora_summary.tsv"))
