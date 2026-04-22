# ---- 03_deseq2_shae.R ----
# DESeq2 differential expression for Schistosoma haematobium
# Uses only the 18 infected (S) samples 

# Input: shae_counts.rds, metadata.rds from 01_qc_and_species_split.R
# Output: DESeq2 results, normalized counts, PCA, volcanos

library(tidyverse)
library(DESeq2)
library(pheatmap)

source("de_functions.R")

# ---- Paths ----

base_dir <- here::here()
counts_dir <- file.path(base_dir, "02_de_analysis/output/counts")
shae_dir <- file.path(base_dir, "02_de_analysis/output/shae")
dir.create(shae_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Data ----

shae_counts <- readRDS(file.path(counts_dir, "shae_counts.rds"))
meta <- readRDS(file.path(counts_dir, "metadata.rds"))

meta_inf <- meta[meta$infection == "infected", ]
shae_inf <- shae_counts[, meta_inf$sample]

# ---- Pre-filter Low-Count Genes ----
# Keep genes with >= 15 counts in >= 3 samples within at least one temperature group

temp_groups <- meta_inf$temp_C
keep <- apply(shae_inf, 1, function(x) {
  any(tapply(x, temp_groups, function(g) sum(g >= 15) >= 3))
})
shae_filt <- shae_inf[keep, ]

# ---- Create DESeqDataSet ----

dds <- DESeqDataSetFromMatrix(
  countData = round(shae_filt),
  colData = meta_inf,
  design = ~ temp_C + niclo_ppm + temp_C:niclo_ppm
)
dds$temp_C <- relevel(dds$temp_C, ref = "24")

# ---- Run DESeq2 ----

dds <- DESeq(dds)
print(resultsNames(dds))

# ---- LFC Shrinkage (apeglm) ----

res_temp_16v24_shr <- lfcShrink(dds, coef = "temp_C_16_vs_24", type = "apeglm")
res_temp_32v24_shr <- lfcShrink(dds, coef = "temp_C_32_vs_24", type = "apeglm")
res_niclo_shr <- lfcShrink(dds, coef = "niclo_ppm_0.05_vs_0", type = "apeglm")

# ---- Normalized Counts ----

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
norm_counts <- counts(dds, normalized = TRUE)

# ---- Exploratory PCA ----

pca_data <- plotPCA(vsd,
                    intgroup = c("temp_C", "niclo_ppm"),
                    returnData = TRUE,
                    ntop = 500)
pct_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = temp_C, shape = niclo_ppm)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("16" = "#2166ac", "24" = "#fdb863", "32" = "#b2182b")) +
  labs(
    x = paste0("PC1 (", pct_var[1], "%)"),
    y = paste0("PC2 (", pct_var[2], "%)"),
    title = "Shae PCA - by temperature and niclosamide (infected samples only)"
  ) +
  theme_minimal()

pdf(file.path(shae_dir, "shae_pca.pdf"), width = 10, height = 8)
print(p_pca)
dev.off()

# ---- Sample Distance Heatmap ----

sample_dists <- dist(t(assay(vsd)))
dist_mat <- as.matrix(sample_dists)

anno_df <- meta_inf %>%
  select(sample, temp_C, niclo_ppm) %>%
  column_to_rownames("sample")

anno_colors <- list(
  temp_C = c("16" = "#2166ac", "24" = "#fdb863", "32" = "#b2182b"),
  niclo_ppm = c("0" = "#d9d9d9", "0.05" = "#525252")
)

pdf(file.path(shae_dir, "shae_sample_distance.pdf"), width = 10, height = 8)
pheatmap(dist_mat,
         annotation_col = anno_df,
         annotation_colors = anno_colors,
         main = "Shae sample distance (VST, infected only)")
dev.off()

# ---- Extract Key Contrasts ----

res_temp_16v24 <- results(dds, name = "temp_C_16_vs_24", alpha = 0.05)
res_temp_32v24 <- results(dds, name = "temp_C_32_vs_24", alpha = 0.05)
res_niclo <- results(dds, name = "niclo_ppm_0.05_vs_0", alpha = 0.05)

summarize_res(res_temp_16v24, "Temperature (16C vs 24C)")
summarize_res(res_temp_32v24, "Temperature (32C vs 24C)")
summarize_res(res_niclo, "Niclosamide (0.05 vs 0)")

# ---- Volcano Plots ----

p_volc_t32 <- make_volcano(res_temp_32v24, "Shae: 32C vs 24C")
p_volc_niclo <- make_volcano(res_niclo, "Shae: Niclosamide 0.05 vs 0 ppm")

pdf(file.path(shae_dir, "shae_volcanos.pdf"), width = 10, height = 8)
print(p_volc_t32)
print(p_volc_niclo)
dev.off()

# Combined PDF
pdf(file.path(shae_dir, "shae_exploratory.pdf"), width = 10, height = 8)
print(p_pca)
pheatmap(dist_mat,
         annotation_col = anno_df,
         annotation_colors = anno_colors,
         main = "Shae sample distance (VST, infected only)")
print(p_volc_t32)
print(p_volc_niclo)
dev.off()

# ---- Save Data Objects ----

saveRDS(dds, file.path(shae_dir, "shae_dds.rds"))
saveRDS(vsd, file.path(shae_dir, "shae_vsd.rds"))
saveRDS(norm_counts, file.path(shae_dir, "shae_norm_counts.rds"))

write_results(res_temp_16v24, file.path(shae_dir, "shae_res_temp16v24.tsv"))
write_results(res_temp_32v24, file.path(shae_dir, "shae_res_temp32v24.tsv"))
write_results(res_niclo, file.path(shae_dir, "shae_res_niclosamide.tsv"))

write_results(res_temp_16v24_shr, file.path(shae_dir, "shae_res_temp16v24_shrunk.tsv"))
write_results(res_temp_32v24_shr, file.path(shae_dir, "shae_res_temp32v24_shrunk.tsv"))
write_results(res_niclo_shr, file.path(shae_dir, "shae_res_niclosamide_shrunk.tsv"))
