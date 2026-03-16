# GO and KEGG GSEA for all main contrasts

library(tidyverse)
library(clusterProfiler)
library(fgsea)
library(enrichplot)

base_dir <- here::here()
annot_dir <- file.path(base_dir, "00_setup/03_GO_annotation/output")
de_dir <- file.path(base_dir, "02_de_analysis/output")
out_dir <- file.path(base_dir, "05_enrichment/output/gsea")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Annotations ----

btru_g2go <- readRDS(file.path(annot_dir, "step05_go_objects/btru_gene2go.Rds"))
btru_go2name <- readRDS(file.path(annot_dir, "step05_go_objects/btru_go2name.Rds"))
btru_kegg <- readRDS(file.path(annot_dir, "step06_kegg_objects/btru_kegg_term2gene.Rds"))

shae_g2go <- readRDS(file.path(annot_dir, "step05_go_objects/shae_gene2go.Rds"))
shae_go2name <- readRDS(file.path(annot_dir, "step05_go_objects/shae_go2name.Rds"))
shae_kegg <- readRDS(file.path(annot_dir, "step06_kegg_objects/shae_kegg_term2gene.Rds"))

# ---- GSEA Function ----
run_gsea <- function(res_file, gene2go, go2name, kegg_t2g,
                     species, contrast, out_dir) {
  
  res <- read_tsv(res_file, show_col_types = FALSE) %>%
    filter(!is.na(stat))
  
  # Rank by Wald statistic (signed?)
  ranked <- setNames(res$stat, res$gene_id)
  ranked <- sort(ranked, decreasing = TRUE)
  
  prefix <- paste(species, contrast, sep = "_")
  
  go_t2g <- gene2go %>% select(GO, gene_id) %>% rename(term = GO, gene = gene_id)
  go_t2n <- go2name %>% select(GO, name) %>% rename(term = GO, name = name)
  
  set.seed(42)
  go_res <- GSEA(
    geneList = ranked,
    TERM2GENE = go_t2g,
    TERM2NAME = go_t2n,
    pvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    eps = 1e-10
  )
  
  go_n <- 0
  if (!is.null(go_res) && nrow(go_res@result %>% filter(p.adjust < 0.05)) > 0) {
    go_n <- sum(go_res@result$p.adjust < 0.05)
    write_tsv(as.data.frame(go_res), file.path(out_dir, paste0(prefix, "_go_gsea.tsv")))
    saveRDS(go_res, file.path(out_dir, paste0(prefix, "_go_gsea.rds")))
    
    # Dot plot
    pdf(file.path(out_dir, paste0(prefix, "_go_gsea_dotplot.pdf")),
        width = 11, height = max(4, min(go_n, 25) * 0.35 + 2))
    print(dotplot(go_res, showCategory = 25, split = ".sign") +
            facet_grid(~ .sign) +
            ggtitle(paste(species, contrast, "- GO GSEA")))
    dev.off()
    
    # Ridge plot
    pdf(file.path(out_dir, paste0(prefix, "_go_gsea_ridge.pdf")),
        width = 11, height = max(4, min(go_n, 20) * 0.4 + 2))
    print(ridgeplot(go_res, showCategory = 20) +
            ggtitle(paste(species, contrast, "- GO GSEA Ridges")))
    dev.off()
  }
  
  # KEGG GSEA
  kegg_t2g_df <- kegg_t2g %>%
    select(KO, gene_id) %>%
    rename(term = KO, gene = gene_id)
  
  set.seed(42)
  kegg_res <- GSEA(
    geneList = ranked,
    TERM2GENE = kegg_t2g_df,
    pvalueCutoff = 0.05,
    minGSSize = 5,
    maxGSSize = 500,
    eps = 1e-10
  )
  
  kegg_n <- 0
  if (!is.null(kegg_res) && nrow(kegg_res@result %>% filter(p.adjust < 0.05)) > 0) {
    kegg_n <- sum(kegg_res@result$p.adjust < 0.05)
    write_tsv(as.data.frame(kegg_res), file.path(out_dir, paste0(prefix, "_kegg_gsea.tsv")))
    saveRDS(kegg_res, file.path(out_dir, paste0(prefix, "_kegg_gsea.rds")))
  }
  
  cat(sprintf("  %s %s: GO=%d sig, KEGG=%d sig\n", species, contrast, go_n, kegg_n))
  
  data.frame(species = species, contrast = contrast,
             go_sig = go_n, kegg_sig = kegg_n)
}

# Run

results <- list()

btru_contrasts <- c("temp24v16", "temp32v16", "infection", "niclosamide")
shae_contrasts <- c("temp24v16", "temp32v16", "niclosamide")

for (ct in btru_contrasts) {
  res_file <- file.path(de_dir, "btru", paste0("btru_res_", ct, ".tsv"))
  results[[paste0("btru_", ct)]] <- run_gsea(
    res_file, btru_g2go, btru_go2name, btru_kegg,
    "btru", ct, out_dir)
}

for (ct in shae_contrasts) {
  res_file <- file.path(de_dir, "shae", paste0("shae_res_", ct, ".tsv"))
  results[[paste0("shae_", ct)]] <- run_gsea(
    res_file, shae_g2go, shae_go2name, shae_kegg,
    "shae", ct, out_dir)
}

summary_df <- bind_rows(results)
write_tsv(summary_df, file.path(out_dir, "gsea_summary.tsv"))
