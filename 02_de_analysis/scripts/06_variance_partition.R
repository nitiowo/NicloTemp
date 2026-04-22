library(tidyverse)
library(DESeq2)
library(variancePartition)
library(BiocParallel)

# ---- Paths ----

base_dir <- here::here()
btru_dir <- file.path(base_dir, "02_de_analysis/btru")
shae_dir <- file.path(base_dir, "02_de_analysis/shae")
counts_dir <- file.path(base_dir, "02_de_analysis/counts")
out_dir <- file.path(base_dir, "02_de_analysis/visualization/variance_partition")
dir.create(out_dir,recursive = TRUE)

# ---- Load Data ----

btru_vsd <- readRDS(file.path(btru_dir, "btru_vsd.rds"))
shae_vsd <- readRDS(file.path(shae_dir, "shae_vsd.rds"))
meta <- readRDS(file.path(counts_dir, "metadata.rds"))

meta_btru <- meta[match(colnames(btru_vsd), meta$sample), , drop = FALSE]
meta_shae <- meta[meta$infection == "infected", , drop = FALSE]
meta_shae <- meta_shae[match(colnames(shae_vsd), meta_shae$sample), , drop = FALSE]

# meta_btru <- meta_btru %>%
#   mutate(temp_C = factor(temp_C),
#          infection = factor(infection),
#          niclo_ppm = factor(niclo_ppm)) %>%
#   as.data.frame()
# rownames(meta_btru) <- meta_btru$sample

# meta_shae <- meta_shae %>%
#   mutate(temp_C = factor(temp_C),
#          niclo_ppm = factor(niclo_ppm)) %>%
#   as.data.frame()
# rownames(meta_shae) <- meta_shae$sample

btru_expr2 <- assay(btru_vsd)
shae_expr2 <- assay(shae_vsd)

# ---- Variance Partitioning ----

btru_form <- ~ temp_C + infection + niclo_ppm
shae_form <- ~ temp_C + niclo_ppm

btru_varpart <- fitExtractVarPartModel(btru_expr, btru_form, meta_btru,
                                       BPPARAM = SerialParam())
shae_varpart <- fitExtractVarPartModel(shae_expr, shae_form, meta_shae,
                                       BPPARAM = SerialParam())

saveRDS(btru_varpart, file.path(out_dir, "btru_varpart_full.rds"))
saveRDS(shae_varpart, file.path(out_dir, "shae_varpart_full.rds"))

summarize_varpart <- function(varpart, species) {
  tibble(
    species = species,
    factor = colnames(varpart),
    mean = apply(varpart, 2, mean, na.rm = TRUE),
    median = apply(varpart, 2, median, na.rm = TRUE),
    q25 = apply(varpart, 2, quantile, probs = 0.25, na.rm = TRUE),
    q75 = apply(varpart, 2, quantile, probs = 0.75, na.rm = TRUE)
  )
}

btru_summary <- summarize_varpart(btru_varpart, "Btru")
shae_summary <- summarize_varpart(shae_varpart, "Shae")

write_tsv(btru_summary, file.path(out, "btru_varpart_summary.tsv"))
write_tsv(shae_summary, file.path(out, "shae_varpart_summary.tsv"))

# ---- Figures ----

pdf(file.path(out_dir, "btru_varpart_violin.pdf"), width = 8, height = 5)
print(plotVarPart(btru_varpart) + ggtitle("Btru variance partitioning"))
dev.off()

pdf(file.path(out_dir, "shae_varpart_violin.pdf"), width = 8, height = 5)
print(plotVarPart(shae_varpart) + ggtitle("Shae variance partitioning"))
dev.off()
