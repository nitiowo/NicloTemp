# MA plots, heatmaps, dispersion plots, contrast overlap 

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(UpSetR)
library(patchwork)

# ---- Paths ----

base_dir <- here::here()
btru_dir <- file.path(base_dir, "02_de_analysis/output/btru")
shae_dir <- file.path(base_dir, "02_de_analysis/output/shae")
counts_dir <- file.path(base_dir, "02_de_analysis/output/counts")
viz_dir <- file.path(base_dir, "02_de_analysis/output/visualization")
dir.create(viz_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Data ----

btru_dds <- readRDS(file.path(btru_dir, "btru_dds.rds"))
btru_vsd <- readRDS(file.path(btru_dir, "btru_vsd.rds"))
shae_dds <- readRDS(file.path(shae_dir, "shae_dds.rds"))
shae_vsd <- readRDS(file.path(shae_dir, "shae_vsd.rds"))
meta <- readRDS(file.path(counts_dir, "metadata.rds"))

# ---- Color Palettes ----

col_inf <- c("uninfected" = "#008080", "infected" = "violet")
col_temp <- c("16" = "#2166ac", "24" = "#fdb863", "32" = "#b2182b")
col_niclo <- c("0" = "#d9d9d9", "0.05" = "#525252")

# ---- Dispersion Plots ----

pdf(file.path(viz_dir, "dispersion_plots.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plotDispEsts(btru_dds, main = "Btru dispersion estimates")
plotDispEsts(shae_dds, main = "Shae dispersion estimates")
par(mfrow = c(1, 1))
dev.off()

# ---- MA Plots ----

# Btru contrasts
btru_contrasts <- list(
  "Infection (S vs R)" = "infection_infected_vs_uninfected",
  "Temp 16C vs 24C" = "temp_C_16_vs_24",
  "Temp 32C vs 24C" = "temp_C_32_vs_24",
  "Niclosamide" = "niclo_ppm_0.05_vs_0"
)

pdf(file.path(viz_dir, "ma_plots_btru.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2))
for (nm in names(btru_contrasts)) {
  res <- results(btru_dds, name = btru_contrasts[[nm]], alpha = 0.05)
  plotMA(res, main = paste("Btru:", nm), ylim = c(-5, 5))
}
par(mfrow = c(1, 1))
dev.off()

# Shae contrasts
shae_contrasts <- list(
  "Temp 16C vs 24C" = "temp_C_16_vs_24",
  "Temp 32C vs 24C" = "temp_C_32_vs_24",
  "Niclosamide" = "niclo_ppm_0.05_vs_0"
)

pdf(file.path(viz_dir, "ma_plots_shae.pdf"), width = 10, height = 5)
par(mfrow = c(1, 3))
for (nm in names(shae_contrasts)) {
  res <- results(shae_dds, name = shae_contrasts[[nm]], alpha = 0.05)
  plotMA(res, main = paste("Shae:", nm), ylim = c(-5, 5))
}
par(mfrow = c(1, 1))
dev.off()

# ---- Top DE Gene Heatmaps ----

# Collect all DE genes across main contrasts for Btru
btru_res_list <- lapply(btru_contrasts, function(coef) {
  results(btru_dds, name = coef, alpha = 0.05)
})

btru_sig_genes <- unique(unlist(lapply(btru_res_list, function(res) {
  rownames(res)[which(res$padj < 0.05)]
})))

if (length(btru_sig_genes) >= 2) {
  # Z-score the VST values for heatmap coloring
  btru_vst_mat <- assay(btru_vsd)[btru_sig_genes, , drop = FALSE]
  btru_vst_scaled <- t(scale(t(btru_vst_mat)))

  btru_anno <- data.frame(
    infection = colData(btru_vsd)$infection,
    temp_C = colData(btru_vsd)$temp_C,
    niclo_ppm = colData(btru_vsd)$niclo_ppm,
    row.names = colnames(btru_vsd)
  )

  # Adjust height based on how many genes
  hmap_height <- max(4, length(btru_sig_genes) * 0.15 + 2)

  pdf(file.path(viz_dir, "btru_de_heatmap.pdf"), width = 12, height = hmap_height)
  pheatmap(btru_vst_scaled,
           annotation_col = btru_anno,
           annotation_colors = list(
             infection = col_inf, temp_C = col_temp, niclo_ppm = col_niclo
           ),
           show_rownames = length(btru_sig_genes) <= 80,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           main = paste("Btru DE genes (n =", length(btru_sig_genes), ") - z-scored VST"),
           fontsize_row = 6)
  dev.off()
}

# Shae: Same approach
shae_res_list <- lapply(shae_contrasts, function(coef) {
  results(shae_dds, name = coef, alpha = 0.05)
})

shae_sig_genes <- unique(unlist(lapply(shae_res_list, function(res) {
  rownames(res)[which(res$padj < 0.05)]
})))

if (length(shae_sig_genes) >= 2) {
  shae_vst_mat <- assay(shae_vsd)[shae_sig_genes, , drop = FALSE]
  shae_vst_scaled <- t(scale(t(shae_vst_mat)))

  shae_anno <- data.frame(
    temp_C = colData(shae_vsd)$temp_C,
    niclo_ppm = colData(shae_vsd)$niclo_ppm,
    row.names = colnames(shae_vsd)
  )

  hmap_height <- max(4, length(shae_sig_genes) * 0.15 + 2)

  pdf(file.path(viz_dir, "shae_de_heatmap.pdf"), width = 10, height = hmap_height)
  pheatmap(shae_vst_scaled,
           annotation_col = shae_anno,
           annotation_colors = list(temp_C = col_temp, niclo_ppm = col_niclo),
           show_rownames = length(shae_sig_genes) <= 80,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           main = paste("Shae DE genes (n =", length(shae_sig_genes), ") - z-scored VST"),
           fontsize_row = 6)
  dev.off()
}

# ---- UpSet Plots ----

btru_de_sets <- lapply(btru_res_list, function(res) {
  rownames(res)[which(res$padj < 0.05)]
})
names(btru_de_sets) <- names(btru_contrasts)

nonempty <- sapply(btru_de_sets, length) > 0
all_genes <- unique(unlist(btru_de_sets))

if (sum(nonempty) >= 2 && length(all_genes) >= 2) {
  btru_de_sets_filt <- btru_de_sets[nonempty]

  # Build binary matrix for UpSetR
  upset_mat <- sapply(btru_de_sets_filt, function(g) as.integer(all_genes %in% g))
  rownames(upset_mat) <- all_genes

  pdf(file.path(viz_dir, "btru_upset_de_overlap.pdf"), width = 8, height = 5)
  print(upset(as.data.frame(upset_mat),
              sets = names(btru_de_sets_filt),
              order.by = "freq",
              main.bar.color = "#2b8cbe",
              text.scale = 1.2))
  dev.off()
}

# ---- PCA with More Components ----

btru_pca <- prcomp(t(assay(btru_vsd)))
pct <- round(100 * btru_pca$sdev^2 / sum(btru_pca$sdev^2), 1)

btru_pca_df <- as.data.frame(btru_pca$x[, 1:4])
btru_pca_df$sample <- rownames(btru_pca_df)
btru_pca_df <- left_join(btru_pca_df,
                          as.data.frame(colData(btru_vsd)),
                          by = "sample")

p12 <- ggplot(btru_pca_df, aes(PC1, PC2, color = temp_C, shape = infection)) +
  geom_point(size = 3) +
  scale_color_manual(values = col_temp) +
  labs(x = paste0("PC1 (", pct[1], "%)"),
       y = paste0("PC2 (", pct[2], "%)"),
       title = "Btru PC1 vs PC2") +
  theme_minimal()

p34 <- ggplot(btru_pca_df, aes(PC3, PC4, color = temp_C, shape = infection)) +
  geom_point(size = 3) +
  scale_color_manual(values = col_temp) +
  labs(x = paste0("PC3 (", pct[3], "%)"),
       y = paste0("PC4 (", pct[4], "%)"),
       title = "Btru PC3 vs PC4") +
  theme_minimal()

# Scree plot
scree_df <- data.frame(
  PC = 1:min(10, length(pct)),
  pct_var = pct[1:min(10, length(pct))]
)
p_scree <- ggplot(scree_df, aes(PC, pct_var)) +
  geom_col(fill = "#2b8cbe") +
  geom_line(group = 1) +
  geom_point() +
  labs(x = "Principal Component", y = "% Variance", title = "Btru scree plot") +
  theme_minimal()

pdf(file.path(viz_dir, "btru_extended_pca.pdf"), width = 14, height = 5)
print(p12 + p34 + p_scree)
dev.off()

# ---- Per-Gene Count Plots for Top Hits ----
# Visualize top genes from the strongest contrast (temperature)

btru_top_res <- results(btru_dds, name = "temp_C_32_vs_24", alpha = 0.05)
top_genes <- head(rownames(btru_top_res[order(btru_top_res$padj), ]), 9)

if (length(top_genes) > 0) {
  btru_norm <- counts(btru_dds, normalized = TRUE)
  meta_btru <- as.data.frame(colData(btru_dds))

  plot_list <- lapply(top_genes, function(g) {
    df <- data.frame(
      count = btru_norm[g, ],
      temp_C = meta_btru$temp_C,
      infection = meta_btru$infection,
      niclo_ppm = meta_btru$niclo_ppm
    )
    ggplot(df, aes(temp_C, count, color = infection)) +
      geom_jitter(width = 0.15, size = 2) +
      scale_color_manual(values = col_inf) +
      labs(title = g, x = "Temperature", y = "Normalized count") +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 8))
  })

  pdf(file.path(viz_dir, "btru_top_gene_counts.pdf"), width = 12, height = 10)
  print(wrap_plots(plot_list, ncol = 3))
  dev.off()
}
