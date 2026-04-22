library(tidyverse)
library(variancePartition)
library(BiocParallel)

# ---- Paths ----

base_dir <- here::here()
counts_dir <- file.path(base_dir, "02_de_analysis/output/counts")
btru_dir <- file.path(base_dir, "03_wgcna/output/btru")
shae_dir <- file.path(base_dir, "03_wgcna/output/shae")
out_dir <- file.path(base_dir, "03_wgcna/output/variance_partition")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Data ----

meta <- readRDS(file.path(counts_dir, "metadata.rds"))
btru_MEs <- readRDS(file.path(btru_dir, "btru_MEs.rds"))
shae_MEs <- readRDS(file.path(shae_dir, "shae_MEs.rds"))

meta_btru <- meta[match(rownames(btru_MEs), meta$sample), , drop = FALSE]
meta_shae <- meta[meta$infection == "infected", , drop = FALSE]
meta_shae <- meta_shae[match(rownames(shae_MEs), meta_shae$sample), , drop = FALSE]

meta_btru <- meta_btru %>%
  mutate(temp_C = factor(temp_C),
         infection = factor(infection),
         niclo_ppm = factor(niclo_ppm)) %>%
  as.data.frame()
rownames(meta_btru) <- meta_btru$sample

meta_shae <- meta_shae %>%
  mutate(temp_C = factor(temp_C),
         niclo_ppm = factor(niclo_ppm)) %>%
  as.data.frame()
rownames(meta_shae) <- meta_shae$sample

btru_expr <- t(as.matrix(btru_MEs))
shae_expr <- t(as.matrix(shae_MEs))

# ---- Variance Partitioning ----

btru_form <- ~ temp_C + infection + niclo_ppm
shae_form <- ~ temp_C + niclo_ppm

btru_varpart <- fitExtractVarPartModel(btru_expr, btru_form, meta_btru,
                                       BPPARAM = SerialParam())
shae_varpart <- fitExtractVarPartModel(shae_expr, shae_form, meta_shae,
                                       BPPARAM = SerialParam())

btru_tbl <- as_tibble(as.data.frame(btru_varpart), rownames = "module_num")
shae_tbl <- as_tibble(as.data.frame(shae_varpart), rownames = "module_num")

write_tsv(btru_tbl, file.path(out_dir, "btru_eigengene_varpart.tsv"))
write_tsv(shae_tbl, file.path(out_dir, "shae_eigengene_varpart.tsv"))

# ---- Figures ----

pdf(file.path(out_dir, "btru_eigengene_varpart.pdf"), width = 9, height = 5)
print(plotPercentBars(sortCols(btru_varpart)) +
        ggtitle("Btru eigengene variance partitioning"))
dev.off()

pdf(file.path(out_dir, "shae_eigengene_varpart.pdf"), width = 8, height = 5)
print(plotPercentBars(sortCols(shae_varpart)) +
        ggtitle("Shae eigengene variance partitioning"))
dev.off()
