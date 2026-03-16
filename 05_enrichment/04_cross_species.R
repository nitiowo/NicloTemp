# Compare enriched GO/KEGG terms between Btru and Shae for shared contrasts

library(tidyverse)

base_dir <- here::here()
ora_dir <- file.path(base_dir, "05_enrichment/output/ora")
gsea_dir <- file.path(base_dir, "05_enrichment/output/gsea")
out_dir <- file.path(base_dir, "05_enrichment/output/cross_species")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Shared contrasts between Btru and Shae
shared_contrasts <- c("temp24v16", "temp32v16", "niclosamide")


# ---- Compare ORA Results ----
# Wonky outputs - check. Skip ORA just use GSEA.
ora_comparison <- list()

for (ct in shared_contrasts) {
  for (dir in c("up", "down")) {
    btru_file <- file.path(ora_dir, paste0("btru_", ct, "_", dir, "_go.tsv"))
    shae_file <- file.path(ora_dir, paste0("shae_", ct, "_", dir, "_go.tsv"))
        
    btru_go <- read_tsv(btru_file, show_col_types = FALSE) %>%
      filter(p.adjust < 0.05)
    shae_go <- read_tsv(shae_file, show_col_types = FALSE) %>%
      filter(p.adjust < 0.05)
    
    shared <- setdiff(btru_go$ID, shae_go$ID)
    btru_only <- setdiff(btru_go$ID, shae_go$ID)
    shae_only <- setdiff(shae_go$ID, btru_go$ID)
        
    if (length(shared) > 0) {
      shared_df <- btru_go %>%
        filter(ID %in% shared) %>%
        select(ID, Description, btru_padj = p, btru_count = Count) %>%
        left_join(
          shae_go %>%
            filter(ID %in% shared) %>%
            select(ID, shae_padj = p, shae_count = Count),
          by = "ID"
        ) %>%
        mutate(contrast = ct, direction = dir)
      
      ora_comparison[[paste(ct, dir)]] <- shared_df
    }
  }
}

if (length(ora_comparison) > 0) {
  ora_shared <- bind_rows(ora_comparison)
  write_tsv(ora_shared, file.path(out_dir, "shared_go_terms_ora.tsv"))
} else {
  cat("\n  No shared GO terms found across species in ORA.\n")
}

# ---- Compare GSEA Results ----

gsea_comparison <- list()

for (ct in shared_contrasts) {
  btru_file <- file.path(gsea_dir, paste0("btru_", ct, "_go_gsea.tsv"))
  shae_file <- file.path(gsea_dir, paste0("shae_", ct, "_go_gsea.tsv"))
  
  if (!file.exists(btru_file) || !file.exists(shae_file)) {
    next
  }
  
  btru_gsea <- read_tsv(btru_file, show_col_types = FALSE) %>%
    filter(p.adjust < 0.05)
  shae_gsea <- read_tsv(shae_file, show_col_types = FALSE) %>%
    filter(p.adjust < 0.05)
  
  shared <- intersect(btru_gsea$ID, shae_gsea$ID)
  
  cat(sprintf("  %s: Btru=%d, Shae=%d, shared=%d\n",
              ct, nrow(btru_gsea), nrow(shae_gsea), length(shared)))
  
  if (length(shared) > 0) {
    shared_df <- btru_gsea %>%
      filter(ID %in% shared) %>%
      select(ID, Description, btru_NES = NES, btru_padj = p.adjust) %>%
      left_join(
        shae_gsea %>%
          filter(ID %in% shared) %>%
          select(ID, shae_NES = NES, shae_padj = p.adjust),
        by = "ID"
      ) %>%
      mutate(contrast = ct,
             same_direction = sign(btru_NES) == sign(shae_NES))
    
    gsea_comparison[[ct]] <- shared_df
  }
}

if (length(gsea_comparison) > 0) {
  gsea_shared <- bind_rows(gsea_comparison)
  write_tsv(gsea_shared, file.path(out_dir, "shared_go_terms_gsea.tsv"))
  
  # NES comparison scatter plot
  p <- ggplot(gsea_shared, aes(x = btru_NES, y = shae_NES, color = contrast)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey30") +
    labs(title = "Cross-Species GSEA: Shared GO Terms",
         x = "Btru NES", y = "Shae NES") +
    theme_bw(base_size = 12)
  
  ggsave(file.path(out_dir, "nes_comparison_scatter.pdf"),
         p, width = 8, height = 6)
} else {
  cat("\n  No shared GO terms found across species in GSEA.\n")
}

# ---- Compare Module Enrichment ----

mod_file <- file.path(base_dir, "05_enrichment/output/module_enrichment/module_enrichment_all.tsv")
if (file.exists(mod_file)) {
  mod_all <- read_tsv(mod_file, show_col_types = FALSE)
  
  btru_terms <- mod_all %>% filter(species == "btru") %>% pull(ID) %>% unique()
  shae_terms <- mod_all %>% filter(species == "shae") %>% pull(ID) %>% unique()
  shared_mod <- intersect(btru_terms, shae_terms)
  
  if (length(shared_mod) > 0) {
    shared_mod_df <- mod_all %>%
      filter(ID %in% shared_mod) %>%
      select(species, module, db, ID, Description, p.adjust) %>%
      pivot_wider(names_from = species,
                  values_from = c(module, p.adjust),
                  values_fn = list(module = ~ paste(., collapse = ","),
                                   p.adjust = ~ min(.)))
    write_tsv(shared_mod_df, file.path(out_dir, "shared_module_terms.tsv"))
  }
}
