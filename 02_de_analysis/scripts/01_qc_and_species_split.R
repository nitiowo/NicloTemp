# ---- 01_qc_and_species_split.R ----
# Inspect nf-core rnaseq output, split counts by species (Btru vs Shae),
# assess alignment quality and make plots
#
# Input: salmon merged counts from nf-core, sample metadata
# Output: species-split count matrices, QC summary, alignment plots

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

# ---- Load Merged Gene Counts ----

counts_raw <- read_tsv(file.path(nfcore_dir, "salmon.merged.gene_counts.tsv"))
tpm_raw <- read_tsv(file.path(nfcore_dir, "salmon.merged.gene_tpm.tsv"))

stopifnot(all(meta$sample %in% colnames(counts_raw)))

# ---- Split By Species ----

btru_idx <- grepl("^Btru_", counts_raw$gene_id)
shae_idx <- grepl("^MS3_", counts_raw$gene_id)

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

# ---- Species read ratio summary ----

btru_totals <- colSums(btru_counts)
shae_totals <- colSums(shae_counts)

species_summary <- meta
species_summary$btru_counts <- btru_totals[meta$sample]
species_summary$shae_counts <- shae_totals[meta$sample]
species_summary$total_counts <- species_summary$btru_counts + species_summary$shae_counts
species_summary$pct_btru <- round(species_summary$btru_counts / species_summary$total_counts * 100, 2)
species_summary$pct_shae <- round(species_summary$shae_counts / species_summary$total_counts * 100, 2)

uninf <- species_summary[species_summary$infection == "uninfected", ]
print(uninf[, c("sample", "btru_counts", "shae_counts", "pct_shae")])

inf <- species_summary[species_summary$infection == "infected", ]
print(inf[, c("sample", "btru_counts", "shae_counts", "pct_shae")])

# ---- Summary Stats ----

uninf_pct <- species_summary$pct_shae[species_summary$infection == "uninfected"]
inf_pct <- species_summary$pct_shae[species_summary$infection == "infected"]

cat("\nShae % in uninfected samples:\n")
cat("  Mean:", mean(uninf_pct), "\n")
cat("  Max:", max(uninf_pct), "\n")

cat("\nShae % in infected samples:\n")
cat("  Mean:", mean(inf_pct), "%\n")
cat("  Range:", range(inf_pct), "%\n")

# ---- Genes per Species ----

btru_expressed <- rowSums(btru_counts >= 1) > 0

inf_samples <- meta$sample[meta$infection == "infected"]
uninf_samples <- meta$sample[meta$infection == "uninfected"]
shae_expressed_inf <- rowSums(shae_counts[, inf_samples] >= 1) > 0
shae_expressed_uninf <- rowSums(shae_counts[, uninf_samples] >= 1) > 0

# ---- STAR Alignment Stats ----

multiqc_dir <- file.path(base_dir, "01_alignment/output/multiqc/star_salmon/multiqc_report_data")

star_stats <- read_tsv(file.path(multiqc_dir, "multiqc_star.txt"))
star_stats <- rename(star_stats, sample = Sample)
star_stats <- left_join(star_stats, meta, by = "sample")

cat("Overall uniquely mapped %: mean =", mean(star_stats$uniquely_mapped_percent),
    ", range =", min(star_stats$uniquely_mapped_percent), "-",
    max(star_stats$uniquely_mapped_percent), "\n")

# ---- Alignment Plots ----

# Species proportion barplot
plot_data <- species_summary %>%
  select(sample, infection, temp_C, niclo_ppm, pct_btru, pct_shae) %>%
  pivot_longer(c(pct_btru, pct_shae), names_to = "species", values_to = "pct") %>%
  mutate(
    species = ifelse(species == "pct_btru", "B. truncatus", "S. haematobium"),
    label = paste(temp_C, "C /", niclo_ppm, "ppm")
  )

p_species <- ggplot(plot_data, aes(x = sample, y = pct, fill = species)) +
  geom_col() +
  facet_wrap(~ infection, scales = "free_x") +
  scale_fill_manual(values = c("B. truncatus" = "#2b8cbe", "S. haematobium" = "#e34a33")) +
  labs(y = "% of assigned counts", x = NULL, fill = "Species",
       title = "Species composition per sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# STAR mapping rate by infection
p_star <- ggplot(star_stats,
                 aes(x = reorder(sample, uniquely_mapped_percent),
                     y = uniquely_mapped_percent,
                     fill = infection)) +
  geom_col() +
  scale_fill_manual(values = c("uninfected" = "#008083", "infected" = "#F78104")) +
  labs(y = "Uniquely mapped (%)", x = NULL,
       title = "STAR unique mapping rate per sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# Shae counts by treatment
shae_treat <- species_summary[species_summary$infection == "infected", ]
p_shae_treat <- ggplot(shae_treat,
                       aes(x = temp_C, y = shae_counts, color = niclo_ppm)) +
  geom_jitter(width = 0.1, size = 3) +
  scale_color_manual(values = c("0" = "#b2b2b2", "0.05" = "#333333")) +
  labs(y = "Total Shae counts", x = "Temperature (C)", color = "Niclosamide",
       title = "Shae counts by treatment (infected only)") +
  theme_minimal()

pdf(file.path(qc_dir, "qc_species_split.pdf"), width = 12, height = 8)
print(p_species)
print(p_star)
print(p_shae_treat)
dev.off()

# ---- Save Split Count Matrices ----

write.csv(species_summary, file.path(qc_dir, "species_summary.csv"), row.names = FALSE)

saveRDS(btru_counts, file.path(counts_dir, "btru_counts.rds"))
saveRDS(shae_counts, file.path(counts_dir, "shae_counts.rds"))
saveRDS(btru_tpm, file.path(counts_dir, "btru_tpm.rds"))
saveRDS(shae_tpm, file.path(counts_dir, "shae_tpm.rds"))
saveRDS(meta, file.path(counts_dir, "metadata.rds"))
