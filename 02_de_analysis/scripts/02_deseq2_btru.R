# ---- 02_deseq2_btru.R ----
# DESeq2 differential expression for Btru
# Uses all 36 samples. Three-factor design: temperature, infection, niclosamide

# Input: btru_counts.rds, metadata.rds from 01_qc_and_species_split.R
# Output: DESeq2 results, normalized counts, PCA, volcanos

library(tidyverse)
library(DESeq2)
library(pheatmap)

source("de_functions.R")

# ---- Paths ----

base_dir <- here::here()
counts_dir <- file.path(base_dir, "02_de_analysis/output/counts")
btru_dir <- file.path(base_dir, "02_de_analysis/output/btru")
dir.create(btru_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Data ----

btru_counts <- readRDS(file.path(counts_dir, "btru_counts.rds"))
meta <- readRDS(file.path(counts_dir, "metadata.rds"))

meta <- meta[match(colnames(btru_counts), meta$sample), ]
stopifnot(all(colnames(btru_counts) == meta$sample))

# ---- Pre-filter Low-Count Genes ----

keep <- rowSums(btru_counts >= 10) >= 3
btru_filt <- btru_counts[keep, ]

# ---- Create DESeqDataSet ----

dds <- DESeqDataSetFromMatrix(
  countData = round(btru_filt),
  colData = meta,
  design = ~ temp_C + infection + niclo_ppm +
    temp_C:infection + temp_C:niclo_ppm + infection:niclo_ppm
)
dds$temp_C <- relevel(dds$temp_C, ref = "24")

# ---- Run DESeq2 ----

dds <- DESeq(dds)
print(resultsNames(dds))

# ---- Custom Averaged Contrasts for Main Effects ----
# Grab the custom contrasts for each effect

mm <- model.matrix(design(dds), data = as.data.frame(colData(dds)))

cntr_infection <- colMeans(mm[dds$infection == "infected", ]) -
                  colMeans(mm[dds$infection == "uninfected", ])
res_infection <- results(dds, contrast = cntr_infection)

cntr_16v24 <- colMeans(mm[dds$temp_C == "16", ]) -
              colMeans(mm[dds$temp_C == "24", ])
res_temp_16v24 <- results(dds, contrast = cntr_16v24)

cntr_32v24 <- colMeans(mm[dds$temp_C == "32", ]) -
              colMeans(mm[dds$temp_C == "24", ])
res_temp_32v24 <- results(dds, contrast = cntr_32v24)

cntr_niclo <- colMeans(mm[dds$niclo_ppm == "0.05", ]) -
              colMeans(mm[dds$niclo_ppm == "0", ])
res_niclo <- results(dds, contrast = cntr_niclo)

# ---- LFC Shrinkage (apeglm) ----

res_infection_shr <- lfcShrink(dds, coef = "infection_infected_vs_uninfected", type = "apeglm")
res_temp_16v24_shr <- lfcShrink(dds, coef = "temp_C_16_vs_24", type = "apeglm")
res_temp_32v24_shr <- lfcShrink(dds, coef = "temp_C_32_vs_24", type = "apeglm")
res_niclo_shr <- lfcShrink(dds, coef = "niclo_ppm_0.05_vs_0", type = "apeglm")

# ---- Normalized Counts ----

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
norm_counts <- counts(dds, normalized = TRUE)

# ---- Exploratory PCA ----

topn <- 5000
pca_data <- plotPCA(vsd,
                    intgroup = c("infection", "temp_C", "niclo_ppm"),
                    returnData = TRUE,
                    ntop = topn)
pct_var <- round(100 * attr(pca_data, "percentVar"))

p_pca_inf <- ggplot(pca_data, aes(x = PC1, y = PC2, color = temp_C, shape = infection)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("16" = "#2166ac", "24" = "#fdb863", "32" = "#b2182b")) +
  labs(
    x = paste0("PC1 (", pct_var[1], "%)"),
    y = paste0("PC2 (", pct_var[2], "%)"),
    title = paste0("Btru PCA - by temperature and infection (top ", topn, " genes)")
  ) +
  theme_minimal()

pdf(file.path(btru_dir, "btru_pca_temp-inf.pdf"), width = 10, height = 8)
print(p_pca_inf)
dev.off()

p_pca_temp <- ggplot(pca_data, aes(x = PC1, y = PC2, color = temp_C, shape = niclo_ppm)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("16" = "#2166ac", "24" = "#fdb863", "32" = "#b2182b")) +
  labs(
    x = paste0("PC1 (", pct_var[1], "%)"),
    y = paste0("PC2 (", pct_var[2], "%)"),
    title = paste0("Btru PCA - by temperature and niclosamide (top ", topn, " genes)")
  ) +
  theme_minimal()

pdf(file.path(btru_dir, "btru_pca_temp-niclo.pdf"), width = 10, height = 8)
print(p_pca_temp)
dev.off()

p_pca_niclo <- ggplot(pca_data, aes(x = PC1, y = PC2, fill = infection, shape = niclo_ppm)) +
  geom_point(size = 3, stroke = 1) +
  scale_fill_manual(values = c("uninfected" = "white", "infected" = "black")) +
  scale_shape_manual(values = c("0" = 21, "0.05" = 24)) +
  labs(
    x = paste0("PC1 (", pct_var[1], "%)"),
    y = paste0("PC2 (", pct_var[2], "%)"),
    title = paste0("Btru PCA - by niclosamide and infection (top ", topn, " genes)")
  ) +
  theme_minimal()

pdf(file.path(btru_dir, "btru_pca_niclo-inf.pdf"), width = 10, height = 8)
print(p_pca_niclo)
dev.off()

p_pca_all <- ggplot(pca_data,
                    aes(x = PC1, y = PC2,
                        color = temp_C, shape = niclo_ppm, fill = infection)) +
  geom_point(size = 3, stroke = 2) +
  scale_color_manual(values = c("16" = "#2166ac", "24" = "#fdb863", "32" = "#b2182b")) +
  scale_fill_manual(values = c("uninfected" = "white", "infected" = "black")) +
  scale_shape_manual(values = c("0" = 21, "0.05" = 24)) +
  labs(
    x = paste0("PC1 (", pct_var[1], "%)"),
    y = paste0("PC2 (", pct_var[2], "%)"),
    title = paste0("Btru PCA - all variables (top ", topn, " genes)")
  ) +
  theme_minimal()

pdf(file.path(btru_dir, "btru_pca_allvars.pdf"), width = 10, height = 8)
print(p_pca_all)
dev.off()

# ---- Sample Distance Heatmap ----

sample_dists <- dist(t(assay(vsd)))
dist_mat <- as.matrix(sample_dists)

anno_df <- meta %>%
  select(sample, infection, temp_C, niclo_ppm) %>%
  column_to_rownames("sample")

anno_colors <- list(
  infection = c("uninfected" = "#008080", "infected" = "orange"),
  temp_C = c("16" = "#2166ac", "24" = "#fdb863", "32" = "#b2182b"),
  niclo_ppm = c("0" = "#d9d9d9", "0.05" = "#525252")
)

pdf(file.path(btru_dir, "btru_sample_distance.pdf"), width = 10, height = 8)
pheatmap(dist_mat,
         annotation_col = anno_df,
         annotation_colors = anno_colors,
         main = "Btru sample distance (VST)")
dev.off()

# ---- Summarize Main Effects ----

summarize_res(res_infection, "Infection (S vs R)")
summarize_res(res_temp_16v24, "Temperature (16C vs 24C)")
summarize_res(res_temp_32v24, "Temperature (32C vs 24C)")
summarize_res(res_niclo, "Niclosamide (0.05 vs 0)")

# ---- Volcano Plots ----

p_volc_inf <- make_volcano(res_infection, "Btru: Infected vs Uninfected")
p_volc_t32 <- make_volcano(res_temp_32v24, "Btru: 32C vs 24C")
p_volc_niclo <- make_volcano(res_niclo, "Btru: Niclosamide 0.05 vs 0 ppm")

pdf(file.path(btru_dir, "btru_volcanos.pdf"), width = 10, height = 8)
print(p_volc_inf)
print(p_volc_t32)
print(p_volc_niclo)
dev.off()

# Combined PDF
pdf(file.path(btru_dir, "btru_exploratory.pdf"), width = 10, height = 8)
print(p_pca_inf)
print(p_pca_temp)
pheatmap(dist_mat,
         annotation_col = anno_df,
         annotation_colors = anno_colors,
         main = "Btru sample distance (VST)")
print(p_volc_inf)
print(p_volc_t32)
print(p_volc_niclo)
dev.off()

# ---- Save Data Objects ----

saveRDS(dds, file.path(btru_dir, "btru_dds.rds"))
saveRDS(vsd, file.path(btru_dir, "btru_vsd.rds"))
saveRDS(norm_counts, file.path(btru_dir, "btru_norm_counts.rds"))

write_results(res_infection, file.path(btru_dir, "btru_res_infection.tsv"))
write_results(res_temp_16v24, file.path(btru_dir, "btru_res_temp16v24.tsv"))
write_results(res_temp_32v24, file.path(btru_dir, "btru_res_temp32v24.tsv"))
write_results(res_niclo, file.path(btru_dir, "btru_res_niclosamide.tsv"))

write_results(res_infection_shr, file.path(btru_dir, "btru_res_infection_shrunk.tsv"))
write_results(res_temp_16v24_shr, file.path(btru_dir, "btru_res_temp16v24_shrunk.tsv"))
write_results(res_temp_32v24_shr, file.path(btru_dir, "btru_res_temp32v24_shrunk.tsv"))
write_results(res_niclo_shr, file.path(btru_dir, "btru_res_niclosamide_shrunk.tsv"))
