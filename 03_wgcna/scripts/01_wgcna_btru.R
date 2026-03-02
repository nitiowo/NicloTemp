library(WGCNA)
library(DESeq2)
library(dplyr)
library(readr)
library(tibble)

allowWGCNAThreads()
cor <- WGCNA::cor

# ---- Paths ----
base_dir <- here::here()
btru_dir <- file.path(base_dir, "02_de_analysis/output/btru")
counts_dir <- file.path(base_dir, "02_de_analysis/output/counts")
out_dir <- file.path(base_dir, "03_wgcna/output/btru")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Data ----
btru_vsd <- readRDS(file.path(btru_dir, "btru_vsd.rds"))
meta <- readRDS(file.path(counts_dir, "metadata.rds"))
meta <- meta[match(colnames(btru_vsd), meta$sample), ]
expr_mat <- t(assay(btru_vsd))

# ---- Filter to n_top genes ----
gene_vars <- apply(expr_mat, 2, var)
n_top <- min(5000, ncol(expr_mat))
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_top]
expr_filt <- expr_mat[, top_genes]

gsg <- goodSamplesGenes(expr_filt, verbose = 3)
if (!gsg$allOK) {
    expr_filt <- expr_filt[gsg$goodSamples, gsg$goodGenes]
}

# ---- Soft Threshold ----
powers <- c(1:10, seq(12, 30, by = 2))
sft <- pickSoftThreshold(expr_filt, powerVector = powers, verbose = 5,
                         networkType = "signed",
                         corFnc = "cor", corOptions = list(use = "p"))

# Default 12
picked_power <- sft$powerEstimate
if (is.na(picked_power)) {
    picked_power <- 12
}

# ---- Build Network and Detect Modules ----
net <- blockwiseModules(
    expr_filt,
    power = picked_power,
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    verbose = 3,
    maxBlockSize = n_top,
    corType = "pearson"
)

module_colors <- labels2colors(net$colors)

# ---- Module Eigengenes ----
MEs <- net$MEs
MEs <- orderMEs(MEs)
colnames(MEs) <- gsub("ME", "", colnames(MEs))

# ---- Trait Correlation ----
traits <- data.frame(
    temp_16 = as.integer(meta$temp_C == "16"),
    temp_24 = as.integer(meta$temp_C == "24"),
    temp_32 = as.integer(meta$temp_C == "32"),
    infected = as.integer(meta$infection == "infected"),
    niclo_treated = as.integer(meta$niclo_ppm == "0.05"),
    row.names = meta$sample
)
traits <- traits[rownames(expr_filt), ]

n_samples <- nrow(expr_filt)
cor_mat <- cor(MEs, traits, use = "p")
pval_mat <- corPvalueStudent(cor_mat, n_samples)

# ---- Save Objects ----

# Gene modules
gene_modules <- data.frame(
    gene_id = colnames(expr_filt),
    module_num = net$colors,
    module_color = module_colors,
    stringsAsFactors = FALSE
)
write_tsv(gene_modules, file.path(out_dir, "btru_gene_modules.tsv"))

# Eigengenes with sample info
me_df <- MEs %>%
    rownames_to_column("sample") %>%
    left_join(meta, by = "sample")
write_tsv(me_df, file.path(out_dir, "btru_module_eigengenes.tsv"))

saveRDS(net, file.path(out_dir, "btru_wgcna_net.rds"))
saveRDS(MEs, file.path(out_dir, "btru_MEs.rds"))
saveRDS(list(cor = cor_mat, pval = pval_mat),
        file.path(out_dir, "btru_module_trait_cor.rds"))
