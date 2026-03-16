# GO and KEGG ORA per WGCNA module for both species

library(tidyverse)
library(clusterProfiler)

base_dir <- here::here()
annot_dir <- file.path(base_dir, "00_setup/03_GO_annotation/output")
wgcna_dir <- file.path(base_dir, "03_wgcna/output")
out_dir <- file.path(base_dir, "05_enrichment/output/module_enrichment")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----Load Annotations ----

btru_g2go_temp <- readRDS(file.path(annot_dir, "step05_go_objects/btru_gene2go.Rds"))
btru_go2name_temp <- readRDS(file.path(annot_dir, "step05_go_objects/btru_go2name.Rds"))
btru_kegg_temp <- readRDS(file.path(annot_dir, "step06_kegg_objects/btru_kegg_term2gene.Rds"))

shae_g2go <- readRDS(file.path(annot_dir, "step05_go_objects/shae_gene2go.Rds"))
shae_go2name <- readRDS(file.path(annot_dir, "step05_go_objects/shae_go2name.Rds"))
shae_kegg <- readRDS(file.path(annot_dir, "step06_kegg_objects/shae_kegg_term2gene.Rds"))

# ---- Module Enrichment Function ----

run_module_ora <- function(gene_mods, gene2go, go2name, kegg_t2g,
                           species, out_dir) {
  
  universe <- gene_mods$gene_id
  colors <- unique(gene_mods$module_color)
  colors <- colors[colors != "grey"]
  
  go_t2g <- gene2go %>% select(GO, gene_id) %>% rename(term = GO, gene = gene_id)
  go_t2n <- go2name %>% select(GO, name) %>% rename(term = GO, name = name)
  kegg_t2g_df <- kegg_t2g %>% select(KO, gene_id) %>% rename(term = KO, gene = gene_id)
  
  all_results <- list()
  
  for (color in colors) {
    mod_genes <- gene_mods$gene_id[gene_mods$module_color == color]
    if (length(mod_genes) < 5) next
    
    # GO ORA
    go_res <- enricher(
      gene = mod_genes,
      universe = universe,
      TERM2GENE = go_t2g,
      TERM2NAME = go_t2n,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      minGSSize = 5,
      maxGSSize = 500
    )
    
    go_n <- 0
    if (!is.null(go_res)) {
      go_sig <- go_res@result %>% filter(p.adjust < 0.05)
      go_n <- nrow(go_sig)
      if (go_n > 0) {
        go_sig$module <- color
        go_sig$species <- species
        go_sig$db <- "GO"
        all_results[[paste0(color, "_go")]] <- go_sig
      }
    }
    
    # KEGG ORA
    kegg_res <- enricher(
      gene = mod_genes,
      universe = universe,
      TERM2GENE = kegg_t2g_df,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      minGSSize = 3,
      maxGSSize = 500
    )
    
    kegg_n <- 0
    if (!is.null(kegg_res)) {
      kegg_sig <- kegg_res@result %>% filter(p.adjust < 0.05)
      kegg_n <- nrow(kegg_sig)
      if (kegg_n > 0) {
        kegg_sig$module <- color
        kegg_sig$species <- species
        kegg_sig$db <- "KEGG"
        all_results[[paste0(color, "_kegg")]] <- kegg_sig
      }
    }
    
    cat(sprintf("  %s %-15s (%d genes): GO=%d, KEGG=%d\n",
                species, color, length(mod_genes), go_n, kegg_n))
  }
  
  if (length(all_results) > 0) {
    combined <- bind_rows(all_results) %>%
      select(species, module, db, ID, Description, GeneRatio, BgRatio,
             pvalue, p.adjust, qvalue, Count)
    return(combined)
  }
  return(NULL)
}

# ---- Run for Both Species ----

btru_mods <- read_tsv(file.path(wgcna_dir, "btru/btru_gene_modules.tsv"),
                      show_col_types = FALSE)
btru_enrich <- run_module_ora(btru_mods, btru_g2go, btru_go2name, btru_kegg,
                              "btru", out_dir)

shae_mods <- read_tsv(file.path(wgcna_dir, "shae/shae_gene_modules.tsv"),
                      show_col_types = FALSE)
shae_enrich <- run_module_ora(shae_mods, shae_g2go, shae_go2name, shae_kegg,
                              "shae", out_dir)


combined <- bind_rows(btru_enrich, shae_enrich)
write_tsv(combined, file.path(out_dir, "module_enrichment_all.tsv"))

if (!is.null(btru_enrich)) {
  write_tsv(btru_enrich, file.path(out_dir, "btru_module_enrichment.tsv"))
}
if (!is.null(shae_enrich)) {
  write_tsv(shae_enrich, file.path(out_dir, "shae_module_enrichment.tsv"))
}

# ---- Heatmap ----
for (sp in c("btru", "shae")) {
  sp_data <- combined %>% filter(species == sp, db == "GO")
  if (nrow(sp_data) == 0) next
  
  # Top n terms per module
  top_terms <- sp_data %>%
    group_by(module) %>%
    slice_min(p.adjust, n = 3) %>%
    ungroup()
  
  if (nrow(top_terms) == 0) next
  
  # Truncate long descriptions: TODO: fix this!! - check Neil code
  top_terms$label <- str_trunc(top_terms$Description, 50)
  
  p <- ggplot(top_terms, aes(x = module, y = label, fill = p.adjust))) +
    geom_tile(color = "white") +
    scale_fill_manual(low = "lightyellow", high = "red3", name = "-log10(padj)") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7)) +
    labs(title = paste(toupper(sp), "- Top GO Terms per Module"),
         x = "Module", y = "")
  
  ggsave(file.path(out_dir, paste0(sp, "_module_go_tile.pdf")),
         p, width = 12, height = max(6, nrow(top_terms) * 0.25 + 2))
}

# ---- Summarize ----

mod_summary <- combined %>%
  group_by(species, module, db) %>%
  summarize(n_sig = n(), .groups = "drop") %>%
  pivot_wider(names_from = db, values_from = n_sig, values_fill = 0)

write_tsv(mod_summary, file.path(out_dir, "module_enrichment_summary.tsv"))
