# ---- 03_deseq2_shae.R ----
# DESeq2 differential expression for Schistosoma haematobium (parasite)
# Uses only the 18 infected (S) samples — Shae has zero counts in uninfected
# Two-factor design: temperature x niclosamide
# Reference temperature: 24C (baseline)
#
# Input: shae_counts.rds, metadata.rds from 01_qc_and_species_split.R
# Output: DESeq2 results, normalized counts, VST
# No plots yet — those come in batch 04

library(tidyverse)
library(DESeq2)

# ---- Paths ----

base_dir <- here::here()
counts_dir <- file.path(base_dir, "02_de_analysis/output/counts")
shae_dir <- file.path(base_dir, "02_de_analysis/output/shae")
dir.create(shae_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Data ----

shae_counts <- readRDS(file.path(counts_dir, "shae_counts.rds"))
meta <- readRDS(file.path(counts_dir, "metadata.rds"))

# Subset to infected samples only
meta_inf <- meta[meta$infection == "infected", ]
shae_inf <- shae_counts[, meta_inf$sample]

cat("Shae count matrix (infected only):", nrow(shae_inf), "genes x", ncol(shae_inf), "samples\n")

# ---- Pre-filter Low-Count Genes ----

# Lower threshold for Shae because read depth is much lower than Btru
keep <- rowSums(shae_inf >= 5) >= 3
shae_filt <- shae_inf[keep, ]
cat("Genes after filtering:", nrow(shae_filt), "/", nrow(shae_inf), "\n")

# ---- Create DESeqDataSet ----

dds <- DESeqDataSetFromMatrix(
  countData = round(shae_filt),
  colData = meta_inf,
  design = ~ temp_C + niclo_ppm + temp_C:niclo_ppm
)

# 24C as reference temperature (baseline condition)
dds$temp_C <- relevel(dds$temp_C, ref = "24")

cat("DESeqDataSet created\n")

# ---- Run DESeq2 ----

dds <- DESeq(dds)
cat("DESeq2 complete\n")
cat("Result names:\n")
print(resultsNames(dds))

# ---- LFC Shrinkage (apeglm) ----
# Shrinks log2FC estimates using an empirical Bayes prior
# Use these ranked lists for GSEA rather than the unshrunken results

res_temp_16v24_shr <- lfcShrink(dds, coef = "temp_C_16_vs_24", type = "apeglm")
res_temp_32v24_shr <- lfcShrink(dds, coef = "temp_C_32_vs_24", type = "apeglm")
res_niclo_shr <- lfcShrink(dds, coef = "niclo_ppm_0.05_vs_0", type = "apeglm")

cat("LFC shrinkage complete\n")

# ---- Normalized Counts ----

# varianceStabilizingTransformation() used instead of vst() — more robust when
# some genes have low counts, which is expected for Shae given lower read depth
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
norm_counts <- counts(dds, normalized = TRUE)

# ---- Extract Key Contrasts ----

res_temp_16v24 <- results(dds, name = "temp_C_16_vs_24", alpha = 0.05)
res_temp_32v24 <- results(dds, name = "temp_C_32_vs_24", alpha = 0.05)
res_niclo <- results(dds, name = "niclo_ppm_0.05_vs_0", alpha = 0.05)

summarize_res <- function(res, label) {
  n_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
  cat(sprintf("%-30s: %d up, %d down (padj < 0.05)\n", label, n_up, n_down))
}

cat("\n=== Shae DE Summary (padj < 0.05) ===\n")
summarize_res(res_temp_16v24, "Temperature (16C vs 24C)")
summarize_res(res_temp_32v24, "Temperature (32C vs 24C)")
summarize_res(res_niclo, "Niclosamide (0.05 vs 0)")

# ---- Shae Read Depth Assessment ----

cat("\n=== Shae Read Depth ===\n")
sample_totals <- colSums(shae_filt)
cat("Total filtered Shae counts per sample:\n")
cat("  Mean:", round(mean(sample_totals)), "\n")
cat("  Min:", min(sample_totals), "\n")
cat("  Max:", max(sample_totals), "\n")

# ---- Save Data Objects ----

saveRDS(dds, file.path(shae_dir, "shae_dds.rds"))
saveRDS(vsd, file.path(shae_dir, "shae_vsd.rds"))
saveRDS(norm_counts, file.path(shae_dir, "shae_norm_counts.rds"))

write_results <- function(res, filename) {
  res_df <- as.data.frame(res)
  res_df <- rownames_to_column(res_df, "gene_id")
  res_df <- arrange(res_df, padj)
  write_tsv(res_df, file.path(shae_dir, filename))
}

write_results(res_temp_16v24, "shae_res_temp16v24.tsv")
write_results(res_temp_32v24, "shae_res_temp32v24.tsv")
write_results(res_niclo, "shae_res_niclosamide.tsv")

# Shrunk results (for GSEA ranked gene lists)
write_results(res_temp_16v24_shr, "shae_res_temp16v24_shrunk.tsv")
write_results(res_temp_32v24_shr, "shae_res_temp32v24_shrunk.tsv")
write_results(res_niclo_shr, "shae_res_niclosamide_shrunk.tsv")

cat("\nAll Shae DESeq2 results saved to:", shae_dir, "\n")
cat("Done.\n")
