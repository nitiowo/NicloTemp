library(tidyverse)
library(WGCNA)
library(DESeq2)

allowWGCNAThreads()
cor <- WGCNA::cor

# ---- Paths ----

base_dir <- here::here()
btru_dir <- file.path(base_dir, "02_de_analysis/output/btru")
counts_dir <- file.path(base_dir, "02_de_analysis/output/counts")
full_out <- file.path(base_dir, "03_wgcna/output/btru")

# ---- Load ----

btru_vsd <- readRDS(file.path(btru_dir, "btru_vsd.rds"))
meta <- readRDS(file.path(counts_dir, "metadata.rds"))
meta <- meta[match(colnames(btru_vsd), meta$sample), ]

# Common genes and soft power from the full network
top_genes  <- readRDS(file.path(full_out, "btru_top_genes.rds"))
soft_power <- readRDS(file.path(full_out, "btru_soft_power.rds"))

# Full expression matrix subset to the top genes
expr_all <- t(assay(btru_vsd))
expr_all <- expr_all[, intersect(top_genes, colnames(expr_all))]

# Defining the subgourps
subgroups <- list(
  infected      = meta$sample[meta$infection == "infected"],
  uninfected    = meta$sample[meta$infection == "uninfected"],
  niclo_treated = meta$sample[meta$niclo_ppm == "0.05"],
  niclo_control = meta$sample[meta$niclo_ppm == "0"],
  temp_16       = meta$sample[meta$temp_C == "16"],
  temp_24       = meta$sample[meta$temp_C == "24"],
  temp_32       = meta$sample[meta$temp_C == "32"]
)

# ---- Building the subnetworks ----

build_subnet <- function(group_name, sample_ids) {
  out_dir <- file.path(base_dir, "03_wgcna/output", paste0("btru_", group_name))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  expr_sub <- expr_all[sample_ids, , drop = FALSE]
  meta_sub <- meta[meta$sample %in% sample_ids, ]

  gsg <- goodSamplesGenes(expr_sub, verbose = 0)
  if (!gsg$allOK) {
    expr_sub <- expr_sub[gsg$goodSamples, gsg$goodGenes]
  }

  net <- blockwiseModules(
    expr_sub,
    power          = soft_power,
    networkType    = "signed",
    TOMType        = "signed",
    minModuleSize  = 20,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels  = TRUE,
    pamRespectsDendro = FALSE,
    verbose        = 0,
    maxBlockSize   = ncol(expr_sub),
    corType        = "pearson"
  )

  module_colors <- labels2colors(net$colors)

  MEs <- orderMEs(net$MEs)
  colnames(MEs) <- gsub("ME", "", colnames(MEs))

  gene_modules <- data.frame(
    gene_id      = colnames(expr_sub),
    module_num   = net$colors,
    module_color = module_colors,
    stringsAsFactors = FALSE
  )
  write_tsv(gene_modules,
            file.path(out_dir, paste0("btru_", group_name, "_gene_modules.tsv")))

  me_df <- MEs %>%
    rownames_to_column("sample") %>%
    left_join(meta_sub, by = "sample")
  write_tsv(me_df,
            file.path(out_dir, paste0("btru_", group_name, "_module_eigengenes.tsv")))

  saveRDS(net, file.path(out_dir, paste0("btru_", group_name, "_wgcna_net.rds")))
  saveRDS(MEs, file.path(out_dir, paste0("btru_", group_name, "_MEs.rds")))

  data.frame(
    group     = group_name,
    n_samples = nrow(expr_sub),
    n_genes   = ncol(expr_sub),
    n_modules = length(unique(net$colors[net$colors != 0]))
  )
}

# ---- Run All Subgroups ----

summary_rows <- lapply(names(subgroups), function(g) build_subnet(g, subgroups[[g]]))

subnet_summary <- bind_rows(summary_rows)
print(subnet_summary)

write_tsv(subnet_summary,
          file.path(base_dir, "03_wgcna/output/btru_subnet_summary.tsv"))
