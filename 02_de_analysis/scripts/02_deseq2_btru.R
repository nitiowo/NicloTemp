# ---- 02_deseq2_btru.R ----
# DESeq2 differential expression for Bulinus truncatus (host)
# Uses all 36 samples. Three-factor design: temperature, infection, niclosamide
# Only fits 2-way interactions (no 3-way — underpowered with n=3)
# Reference temperature: 24C (baseline)
#
# Input: btru_counts.rds, metadata.rds from 01_qc_and_species_split.R
# Output: DESeq2 results, normalized counts, VST
# No plots yet — those come in batch 04

library(tidyverse)
library(DESeq2)

# ---- Paths ----

base_dir <- here::here()
counts_dir <- file.path(base_dir, "02_de_analysis/output/counts")
btru_dir <- file.path(base_dir, "02_de_analysis/output/btru")
dir.create(btru_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Data ----

btru_counts <- readRDS(file.path(counts_dir, "btru_counts.rds"))
meta <- readRDS(file.path(counts_dir, "metadata.rds"))

# Reorder metadata to match count matrix columns
meta <- meta[match(colnames(btru_counts), meta$sample), ]
stopifnot(all(colnames(btru_counts) == meta$sample))

cat("Btru count matrix:", nrow(btru_counts), "genes x", ncol(btru_counts), "samples\n")

# ---- Pre-filter Low-Count Genes ----

# Keep genes with >= 10 counts in at least 3 samples (smallest group size)
keep <- rowSums(btru_counts >= 10) >= 3
btru_filt <- btru_counts[keep, ]
cat("Genes after filtering:", nrow(btru_filt), "/", nrow(btru_counts), "\n")

# ---- Create DESeqDataSet ----

# Full model with all 2-way interactions
dds <- DESeqDataSetFromMatrix(
  countData = round(btru_filt),
  colData = meta,
  design = ~ temp_C + infection + niclo_ppm +
    temp_C:infection + temp_C:niclo_ppm + infection:niclo_ppm
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

res_infection_shr <- lfcShrink(dds, coef = "infection_infected_vs_uninfected", type = "apeglm")
res_temp_16v24_shr <- lfcShrink(dds, coef = "temp_C_16_vs_24", type = "apeglm")
res_temp_32v24_shr <- lfcShrink(dds, coef = "temp_C_32_vs_24", type = "apeglm")
res_niclo_shr <- lfcShrink(dds, coef = "niclo_ppm_0.05_vs_0", type = "apeglm")

cat("LFC shrinkage complete\n")

# ---- Normalized Counts ----

# varianceStabilizingTransformation() used instead of vst() because vst()
# can fail with unreliable size factors when low-count genes are present
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
norm_counts <- counts(dds, normalized = TRUE)

# ---- Extract Key Contrasts ----

res_infection <- results(dds, name = "infection_infected_vs_uninfected", alpha = 0.05)
res_temp_16v24 <- results(dds, name = "temp_C_16_vs_24", alpha = 0.05)
res_temp_32v24 <- results(dds, name = "temp_C_32_vs_24", alpha = 0.05)
res_niclo <- results(dds, name = "niclo_ppm_0.05_vs_0", alpha = 0.05)

# Summarize results
summarize_res <- function(res, label) {
  n_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
  cat(sprintf("%-30s: %d up, %d down (padj < 0.05)\n", label, n_up, n_down))
}

cat("\n=== Btru DE Summary (padj < 0.05) ===\n")
summarize_res(res_infection, "Infection (S vs R)")
summarize_res(res_temp_16v24, "Temperature (16C vs 24C)")
summarize_res(res_temp_32v24, "Temperature (32C vs 24C)")
summarize_res(res_niclo, "Niclosamide (0.05 vs 0)")

# ---- Save Data Objects ----

saveRDS(dds, file.path(btru_dir, "btru_dds.rds"))
saveRDS(vsd, file.path(btru_dir, "btru_vsd.rds"))
saveRDS(norm_counts, file.path(btru_dir, "btru_norm_counts.rds"))

# Write results as TSVs
write_results <- function(res, filename) {
  res_df <- as.data.frame(res)
  res_df <- rownames_to_column(res_df, "gene_id")
  res_df <- arrange(res_df, padj)
  write_tsv(res_df, file.path(btru_dir, filename))
}

write_results(res_infection, "btru_res_infection.tsv")
write_results(res_temp_16v24, "btru_res_temp16v24.tsv")
write_results(res_temp_32v24, "btru_res_temp32v24.tsv")
write_results(res_niclo, "btru_res_niclosamide.tsv")

# Shrunk results (for GSEA ranked gene lists)
write_results(res_infection_shr, "btru_res_infection_shrunk.tsv")
write_results(res_temp_16v24_shr, "btru_res_temp16v24_shrunk.tsv")
write_results(res_temp_32v24_shr, "btru_res_temp32v24_shrunk.tsv")
write_results(res_niclo_shr, "btru_res_niclosamide_shrunk.tsv")

cat("\nAll Btru DESeq2 results saved to:", btru_dir, "\n")
cat("Done.\n")
