library(tidyverse)

base_dir <- here::here()
annot_dir <- file.path(base_dir, "00_setup/03_GO_annotation/output")
de_dir <- file.path(base_dir, "02_de_analysis/output")
wgcna_dir <- file.path(base_dir, "03_wgcna/output")
out_dir <- file.path(base_dir, "05_enrichment/output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

build_master_table <- function(species, contrasts, de_dir, wgcna_dir, annot_dir, out_dir) {
  
  base <- read_tsv(file.path(de_dir, species, paste0(species, "_res_temp32v24.tsv")),
                   show_col_types = FALSE) %>%
    select(gene_id, baseMean)
  
# join all DE contrasts
  contrast_tables <- map(contrasts, function(ct) {
    read_tsv(file.path(de_dir, species, paste0(species, "_res_", ct, "_shrunk.tsv")),
             show_col_types = FALSE) %>%
      select(gene_id,
             !!paste0("lfc_", ct) := log2FoldChange,
             !!paste0("padj_", ct) := padj)
  })
  base <- reduce(append(list(base), contrast_tables), left_join, by = "gene_id")
  
  mods <- read_tsv(file.path(wgcna_dir, species, paste0(species, "_gene_modules.tsv")),
                   show_col_types = FALSE) %>%
    select(gene_id, wgcna_module = module_color)
  base <- left_join(base, mods, by = "gene_id")
  
  g2go <- readRDS(file.path(annot_dir, "step05_go_objects", paste0(species, "_gene2go.Rds")))
  go2name <- readRDS(file.path(annot_dir, "step05_go_objects", paste0(species, "_go2name.Rds")))
  
  go_annot <- g2go %>%
    left_join(select(go2name, GO, name, ontology), by = "GO") %>%
    group_by(gene_id) %>%
    summarize(go_terms = paste(unique(GO), collapse = "; "),
              go_names = paste(unique(na.omit(name)), collapse = "; "),
              go_ontologies = paste(unique(na.omit(ontology)), collapse = "; "),
              n_go = n_distinct(GO),
              .groups = "drop")
  base <- left_join(base, go_annot, by = "gene_id")
  
  kegg <- readRDS(file.path(annot_dir, "step06_kegg_objects", paste0(species, "_kegg_term2gene.Rds")))
  kegg_annot <- kegg %>%
    group_by(gene_id) %>%
    summarize(kegg_kos = paste(unique(KO), collapse = "; "),
              n_kegg = n_distinct(KO),
              .groups = "drop")
  base <- left_join(base, kegg_annot, by = "gene_id")
  
# Calculate significance summaries across all contrasts
  padj_cols <- paste0("padj_", contrasts)
  
  base <- base %>%
    mutate(
      # Automatically sums rows across the padj columns
      n_sig_contrasts = rowSums(pick(all_of(padj_cols)) < 0.05, na.rm = TRUE),
      sig_any = n_sig_contrasts > 0
    )
  
  write_tsv(base, file.path(out_dir, paste0(species, "_gene_table.tsv")))
  return(base)
}

# Execute for both species

btru_base <- build_master_table(
  species   = "btru",
  contrasts = c("temp16v24", "temp32v24", "infection", "niclosamide"),
  de_dir = de_dir, wgcna_dir = wgcna_dir, annot_dir = annot_dir, out_dir = out_dir
)

shae_base <- build_master_table(
  species   = "shae",
  contrasts = c("temp16v24", "temp32v24", "niclosamide"),
  de_dir = de_dir, wgcna_dir = wgcna_dir, annot_dir = annot_dir, out_dir = out_dir
)
