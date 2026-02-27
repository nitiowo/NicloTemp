# ---- 01_qc_and_species_split.R ----
# Inspect nf-core rnaseq output, split counts by species (Btru vs Shae),
# assess alignment quality, and save split count matrices
#
# Input: salmon merged counts from nf-core, sample metadata
# Output: species-split count matrices, QC summary, metadata RDS
# No plots yet — those come in batch 04

library(tidyverse)

# ---- Paths ----

base_dir <- here::here()
nfcore_dir <- file.path(base_dir, "01_alignment/output/star_salmon")

counts_dir <- file.path(base_dir, "02_de_analysis/output/counts")
qc_dir <- file.path(base_dir, "02_de_analysis/output/qc")
dir.create(counts_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Metadata ----

meta <- read_csv(file.path(base_dir, "metadata/sample_metadata.csv"))
meta$infection <- factor(meta$infection, levels = c("uninfected", "infected"))
meta$temperature <- factor(meta$temperature, levels = c("A", "C", "E"))
meta$niclosamide <- factor(meta$niclosamide, levels = c("1", "2"))
meta$temp_C <- factor(meta$temp_C, levels = c(16, 24, 32))
meta$niclo_ppm <- factor(meta$niclo_ppm, levels = c(0, 0.05))

cat("Metadata loaded:", nrow(meta), "samples\n")

# ---- Load Merged Gene Counts ----

counts_raw <- read_tsv(file.path(nfcore_dir, "salmon.merged.gene_counts.tsv"))
tpm_raw <- read_tsv(file.path(nfcore_dir, "salmon.merged.gene_tpm.tsv"))

cat("Total genes in count matrix:", nrow(counts_raw), "\n")
cat("Total samples in count matrix:", ncol(counts_raw) - 2, "\n")

# Confirm all metadata samples are in the count matrix
stopifnot(all(meta$sample %in% colnames(counts_raw)))

# ---- Split by Species ----

# Btru genes start with "Btru_", Shae genes start with "MS3_"
btru_idx <- grepl("^Btru_", counts_raw$gene_id)
shae_idx <- grepl("^MS3_", counts_raw$gene_id)

cat("Btru genes:", sum(btru_idx), "\n")
cat("Shae genes:", sum(shae_idx), "\n")
cat("Other genes:", sum(!btru_idx & !shae_idx), "\n")

# Build a count matrix: gene_id as rownames, sample columns only
make_mat <- function(df, idx) {
  subset_df <- df[idx, ]
  subset_df <- subset_df %>%
    select(-gene_name) %>%
    column_to_rownames("gene_id")
  mat <- as.matrix(subset_df)
  mat[, meta$sample]
}

btru_counts <- make_mat(counts_raw, btru_idx)
shae_counts <- make_mat(counts_raw, shae_idx)
btru_tpm <- make_mat(tpm_raw, btru_idx)
shae_tpm <- make_mat(tpm_raw, shae_idx)

cat("\nBtru count matrix:", nrow(btru_counts), "genes x", ncol(btru_counts), "samples\n")
cat("Shae count matrix:", nrow(shae_counts), "genes x", ncol(shae_counts), "samples\n")

# ---- Species Read Proportions per Sample ----

btru_totals <- colSums(btru_counts)
shae_totals <- colSums(shae_counts)

species_summary <- meta
species_summary$btru_counts <- btru_totals[meta$sample]
species_summary$shae_counts <- shae_totals[meta$sample]
species_summary$total_counts <- species_summary$btru_counts + species_summary$shae_counts
species_summary$pct_btru <- round(species_summary$btru_counts / species_summary$total_counts * 100, 2)
species_summary$pct_shae <- round(species_summary$shae_counts / species_summary$total_counts * 100, 2)

cat("\n=== Species Read Proportions ===\n")
cat("\nUninfected (R) samples — Shae counts should be ~0:\n")
uninf <- species_summary[species_summary$infection == "uninfected", ]
print(uninf[, c("sample", "btru_counts", "shae_counts", "pct_shae")])

cat("\nInfected (S) samples — expect some Shae reads:\n")
inf <- species_summary[species_summary$infection == "infected", ]
print(inf[, c("sample", "btru_counts", "shae_counts", "pct_shae")])

# ---- Summary Stats ----

cat("\n=== Summary Statistics ===\n")
uninf_pct <- species_summary$pct_shae[species_summary$infection == "uninfected"]
inf_pct <- species_summary$pct_shae[species_summary$infection == "infected"]

cat("\nShae % in uninfected samples:\n")
cat("  Mean:", mean(uninf_pct), "\n")
cat("  Max:", max(uninf_pct), "\n")

cat("\nShae % in infected samples:\n")
cat("  Mean:", mean(inf_pct), "%\n")
cat("  Range:", range(inf_pct), "%\n")

# ---- Expressed Genes per Species ----

btru_expressed <- rowSums(btru_counts >= 1) > 0

inf_samples <- meta$sample[meta$infection == "infected"]
uninf_samples <- meta$sample[meta$infection == "uninfected"]
shae_expressed_inf <- rowSums(shae_counts[, inf_samples] >= 1) > 0
shae_expressed_uninf <- rowSums(shae_counts[, uninf_samples] >= 1) > 0

cat("\n=== Gene Detection ===\n")
cat("Btru genes detected (>= 1 count in any sample):", sum(btru_expressed), "/", nrow(btru_counts), "\n")
cat("Shae genes detected in infected samples:", sum(shae_expressed_inf), "/", nrow(shae_counts), "\n")
cat("Shae genes detected in uninfected samples:", sum(shae_expressed_uninf), "/", nrow(shae_counts), "\n")

# ---- STAR Alignment Stats ----

multiqc_dir <- file.path(base_dir, "01_alignment/output/multiqc/star_salmon/multiqc_report_data")

star_stats <- read_tsv(file.path(multiqc_dir, "multiqc_star.txt"))
star_stats <- rename(star_stats, sample = Sample)
star_stats <- left_join(star_stats, meta, by = "sample")

cat("\n=== STAR Mapping Rates ===\n")
cat("Overall uniquely mapped %: mean =", mean(star_stats$uniquely_mapped_percent),
    ", range =", min(star_stats$uniquely_mapped_percent), "-",
    max(star_stats$uniquely_mapped_percent), "\n")

# ---- Save Split Count Matrices ----

write.csv(species_summary, file.path(qc_dir, "species_summary.csv"), row.names = FALSE)

saveRDS(btru_counts, file.path(counts_dir, "btru_counts.rds"))
saveRDS(shae_counts, file.path(counts_dir, "shae_counts.rds"))
saveRDS(btru_tpm, file.path(counts_dir, "btru_tpm.rds"))
saveRDS(shae_tpm, file.path(counts_dir, "shae_tpm.rds"))
saveRDS(meta, file.path(counts_dir, "metadata.rds"))

cat("\nCount matrices and metadata saved to:", counts_dir, "\n")
cat("Done.\n")
