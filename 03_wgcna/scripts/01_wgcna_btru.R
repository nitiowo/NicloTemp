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

# Sample dendrogram
sample_tree <- hclust(dist(expr_filt), method = "average")

pdf(file.path(out_dir, "sample_dendrogram.pdf"), width = 12, height = 5)
plot(sample_tree, main = "Btru sample clustering",
     sub = "", xlab = "", cex = 0.7)
dev.off()

# ---- Soft Threshold ----
powers <- c(1:10, seq(12, 30, by = 2))
sft <- pickSoftThreshold(expr_filt, powerVector = powers, verbose = 5,
                         networkType = "signed",
                         corFnc = "cor", corOptions = list(use = "p"))

# Plot topology fit and  connectivity
pdf(file.path(out_dir, "soft_threshold.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R^2)",
     main = "Scale independence", type = "n")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.80, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, col = "red")

par(mfrow = c(1, 1))
dev.off()

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
print(table(module_colors))

# Module dendrogram
pdf(file.path(out_dir, "module_dendrogram.pdf"), width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]],
                    module_colors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Btru gene dendrogram and module colors")
dev.off()

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

# Pairwise interaction columns
traits$tempXinf_16   <- traits$temp_16 * traits$infected
traits$tempXinf_32   <- traits$temp_32 * traits$infected
traits$tempXniclo_16 <- traits$temp_16 * traits$niclo_treated
traits$tempXniclo_32 <- traits$temp_32 * traits$niclo_treated
traits$infXniclo     <- traits$infected * traits$niclo_treated

# Display order
trait_order <- c("temp_16", "temp_24", "temp_32", "infected", "niclo_treated",
                 "tempXinf_16", "tempXinf_32", "tempXniclo_16", "tempXniclo_32", "infXniclo")
traits <- traits[, trait_order]

n_samples <- nrow(expr_filt)
cor_mat <- cor(MEs, traits, use = "p")
pval_mat <- corPvalueStudent(cor_mat, n_samples)

text_mat <- matrix(
    paste0(signif(cor_mat, 2), "\n(", signif(pval_mat, 1), ")"),
    nrow = nrow(cor_mat)
)

pdf(file.path(out_dir, "module_trait_correlation.pdf"),
    width = 8, height = max(4, nrow(cor_mat) * 0.4 + 2))
par(mar = c(6, 8, 3, 3))
labeledHeatmap(
    Matrix = cor_mat,
    xLabels = colnames(traits),
    yLabels = colnames(MEs),
    ySymbols = colnames(MEs),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = text_mat,
    setStdMargins = FALSE,
    cex.text = 0.6,
    zlim = c(-1, 1),
    main = "Btru module-trait correlation"
)
dev.off()

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
saveRDS(top_genes, file.path(out_dir, "btru_top_genes.rds"))
saveRDS(picked_power, file.path(out_dir, "btru_soft_power.rds"))
saveRDS(list(cor = cor_mat, pval = pval_mat),
        file.path(out_dir, "btru_module_trait_cor.rds"))
