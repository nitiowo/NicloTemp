# Cluster temperature-responsive genes by expression pattern across 16-24-32C

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(clusterProfiler)

base_dir <- here::here()
annot_dir <- file.path(base_dir, "00_setup/03_GO_annotation/output")
de_dir <- file.path(base_dir, "02_de_analysis/output")
out_dir <- file.path(base_dir, "05_enrichment/output/dose_response")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Process Each Species ----

for (species in c("btru", "shae")) {
  
  res_24 <- read_tsv(file.path(de_dir, species,
                                paste0(species, "_res_temp24v16.tsv")),
                     show_col_types = FALSE) %>%
    filter(!is.na(padj)) %>%
    select(gene_id, lfc_24 = log2FoldChange, padj_24 = padj)
  
  res_32 <- read_tsv(file.path(de_dir, species,
                                paste0(species, "_res_temp32v16.tsv")),
                     show_col_types = FALSE) %>%
    filter(!is.na(padj)) %>%
    select(gene_id, lfc_32 = log2FoldChange, padj_32 = padj)
  
  combined <- inner_join(res_24, res_32, by = "gene_id") %>%
    filter(padj_24 < 0.05 | padj_32 < 0.05)
  
  
  # Classify response patterns
  combined <- combined %>%
    mutate(
      pattern = case_when(
        lfc_24 > 0 & lfc_32 > 0 & lfc_32 > lfc_24 ~ "linear_up",
        lfc_24 < 0 & lfc_32 < 0 & lfc_32 < lfc_24 ~ "linear_down",
        lfc_24 > 0 & lfc_32 < 0 ~ "peak_24",
        lfc_24 < 0 & lfc_32 > 0 ~ "trough_24",
        lfc_24 > 0 & lfc_32 > 0 & lfc_32 <= lfc_24 ~ "plateau_up",
        lfc_24 < 0 & lfc_32 < 0 & lfc_32 >= lfc_24 ~ "plateau_down",
        abs(lfc_24) < 0.5 & abs(lfc_32) > 1 ~ "late_response",
        abs(lfc_32) < 0.5 & abs(lfc_24) > 1 ~ "early_response",
        TRUE ~ "other"
      )
    )
  
  pattern_counts <- combined %>% dplyr::count(pattern) %>% arrange(desc(n))
  for (i in seq_len(nrow(pattern_counts))) {
                pattern_counts$n[i]))
  }
  
  write_tsv(combined, file.path(out_dir, paste0(species, "_temp_patterns.tsv")))
  write_tsv(pattern_counts, file.path(out_dir, paste0(species, "_pattern_counts.tsv")))
  
  # ---- K-Means Clustering on LFC Profiles ----
  
  lfc_mat <- combined %>%
    select(gene_id, lfc_24, lfc_32) %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
  
  lfc_mat_full <- cbind(lfc_16 = 0, lfc_mat)
  
  lfc_scaled <- t(scale(t(lfc_mat_full)))
  lfc_scaled[is.nan(lfc_scaled)] <- 0
  
  k <- 6
  set.seed(42)
  km <- kmeans(lfc_scaled, centers = k, nstart = 25, iter.max = 50)
  
  combined$cluster <- km$cluster
    
  plot_data <- combined %>%
    select(gene_id, lfc_24, lfc_32, cluster) %>%
    mutate(lfc_16 = 0) %>%
    pivot_longer(cols = starts_with("lfc_"),
                 names_to = "temp", values_to = "lfc") %>%
    mutate(temp_C = as.numeric(gsub("lfc_", "", temp)),
           cluster = paste0("Cluster ", cluster))
  
  cluster_means <- plot_data %>%
    group_by(cluster, temp_C) %>%
    summarize(mean_lfc = mean(lfc), .groups = "drop")
  
  # Counts per cluster
  clust_n <- combined %>%
    dplyr::count(cluster) %>%
    mutate(label = paste0("Cluster ", cluster, " (n=", n, ")"))
  
  p <- ggplot(plot_data, aes(x = temp_C, y = lfc, group = gene_id)) +
    geom_line(alpha = 0.05, color = "grey50") +
    geom_line(data = cluster_means,
              aes(x = temp_C, y = mean_lfc, group = cluster),
              color = "red", linewidth = 1.2, inherit.aes = FALSE) +
    geom_point(data = cluster_means,
               aes(x = temp_C, y = mean_lfc),
               color = "red", size = 2, inherit.aes = FALSE) +
    facet_wrap(~ cluster, scales = "free_y") +
    labs(title = paste(toupper(species), "- Temperature Response Clusters"),
         x = "Temperature (C)", y = "log2 Fold Change vs 16C") +
    scale_x_continuous(breaks = c(16, 24, 32)) +
    theme_bw(base_size = 11)
  
  ggsave(file.path(out_dir, paste0(species, "_temp_clusters.pdf")),
         p, width = 12, height = 8)
  
  # ---- Enrichment per Cluster ----
  
  g2go <- readRDS(file.path(annot_dir, "step05_go_objects",
                            paste0(species, "_gene2go.Rds")))
  go2name <- readRDS(file.path(annot_dir, "step05_go_objects",
                               paste0(species, "_go2name.Rds")))
  
  go_t2g <- g2go %>% select(GO, gene_id) %>% rename(term = GO, gene = gene_id)
  go_t2n <- go2name %>% select(GO, name) %>% rename(term = GO, name = name)
  
  universe <- combined$gene_id
  
  cluster_enrich <- list()
  for (cl in sort(unique(combined$cluster))) {
    cl_genes <- combined$gene_id[combined$cluster == cl]
    if (length(cl_genes) < 10) next
    
    go_res <- enricher(
      gene = cl_genes,
      universe = universe,
      TERM2GENE = go_t2g,
      TERM2NAME = go_t2n,
      pvalueCutoff = 0.1,
      minGSSize = 5,
      maxGSSize = 500
    )
    
    if (!is.null(go_res)) {
      sig <- go_res@result %>% filter(p.adjust < 0.1)
      if (nrow(sig) > 0) {
        sig$cluster <- cl
        cluster_enrich[[as.character(cl)]] <- sig
      }
    }
  }
  
  if (length(cluster_enrich) > 0) {
    all_enrich <- bind_rows(cluster_enrich)
    write_tsv(all_enrich, file.path(out_dir, paste0(species, "_cluster_go.tsv")))
  }
  
  write_tsv(combined %>% select(gene_id, lfc_24, lfc_32, pattern, cluster),
            file.path(out_dir, paste0(species, "_temp_clusters_full.tsv")))
}

