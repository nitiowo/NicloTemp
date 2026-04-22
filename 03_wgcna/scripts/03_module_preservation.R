library(tidyverse)
library(WGCNA)
library(DESeq2)

allowWGCNAThreads()
cor <- WGCNA::cor

# ---- Paths ----

base_dir <- here::here()
btru_dir <- file.path(base_dir, "02_de_analysis/output/btru")
counts_dir <- file.path(base_dir, "02_de_analysis/output/counts")
wgcna_dir <- file.path(base_dir, "03_wgcna/output/btru")
out_dir <- file.path(base_dir, "03_wgcna/output/preservation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Data ----

btru_vsd <- readRDS(file.path(btru_dir, "btru_vsd.rds"))
meta <- readRDS(file.path(counts_dir, "metadata.rds"))
meta <- meta[match(colnames(btru_vsd), meta$sample), ]

net <- readRDS(file.path(wgcna_dir, "btru_wgcna_net.rds"))
gene_mods <- read_tsv(file.path(wgcna_dir, "btru_gene_modules.tsv"),
                      show_col_types = FALSE)

expr_full <- t(assay(btru_vsd))
net_genes <- gene_mods$gene_id
expr_full <- expr_full[, net_genes]

module_colors <- setNames(gene_mods$module_color, gene_mods$gene_id)

# ----Contrasts ----

contrasts <- list(
  infection = list(
    ref_label = "uninfected",
    test_label = "infected",
    ref_idx = which(meta$infection == "uninfected"),
    test_idx = which(meta$infection == "infected")
  ),
  temp_32v24 = list(
    ref_label = "24C",
    test_label = "32C",
    ref_idx = which(meta$temp_C == "24"),
    test_idx = which(meta$temp_C == "32")
  ),
  temp_16v24 = list(
    ref_label = "24C",
    test_label = "16C",
    ref_idx = which(meta$temp_C == "24"),
    test_idx = which(meta$temp_C == "16")
  ),
  niclosamide = list(
    ref_label = "0ppm",
    test_label = "0.05ppm",
    ref_idx = which(meta$niclo_ppm == "0"),
    test_idx = which(meta$niclo_ppm == "0.05")
  )
)

# ----Run Module Preservation ----
# Check Neil's vignette for comments
all_results <- list()

for (cname in names(contrasts)) {
  ct <- contrasts[[cname]]

  if (length(ct$ref_idx) < 6 || length(ct$test_idx) < 6) next

  ref_expr <- expr_full[ct$ref_idx, ]
  test_expr <- expr_full[ct$test_idx, ]

  gsg_ref <- goodSamplesGenes(ref_expr, verbose = 0)
  gsg_test <- goodSamplesGenes(test_expr, verbose = 0)
  keep_genes <- gsg_ref$goodGenes & gsg_test$goodGenes

  if (sum(keep_genes) < 100) next

  ref_expr <- ref_expr[, keep_genes]
  test_expr <- test_expr[, keep_genes]

  multiExpr <- list(
    reference = list(data = ref_expr),
    test = list(data = test_expr)
  )
  multiColor <- list(reference = module_colors[colnames(ref_expr)])

  mp <- modulePreservation(
    multiExpr,
    multiColor,
    referenceNetworks = 1,
    testNetworks = 2,
    nPermutations = 200,
    randomSeed = 42,
    quickCor = 0,
    verbose = 3
  )

  saveRDS(mp, file.path(out_dir, paste0("mp_", cname, ".rds")))

  ref_stat <- mp$preservation$Z[[1]][[2]]
  obs_stat <- mp$preservation$observed[[1]][[2]]

  z_df <- data.frame(
    contrast = cname,
    module = rownames(ref_stat),
    size = ref_stat[, "moduleSize"],
    z_summary = ref_stat[, "Zsummary.pres"],
    z_density = ref_stat[, "Zdensity.pres"],
    z_connectivity = ref_stat[, "Zconnectivity.pres"],
    preserved = ifelse(ref_stat[, "Zsummary.pres"] > 10, "strongly",
                       ifelse(ref_stat[, "Zsummary.pres"] > 2, "moderately",
                              "not preserved"))
  )

  all_results[[cname]] <- z_df

  p <- ggplot(z_df[z_df$module != "gold", ], # exclude grey instead of gold here
              aes(x = reorder(module, z_summary), y = z_summary, fill = preserved)) +
    geom_col() +
    geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "blue") +
    coord_flip() +
    scale_fill_manual(values = c("strongly" = "#2166ac",
                                 "moderately" = "#fdb863",
                                 "not preserved" = "#b2182b")) +
    labs(title = paste("Module preservation:", ct$ref_label, "\u2192", ct$test_label),
         x = "Module", y = "Z-summary") +
    theme_minimal()

  ggsave(file.path(out_dir, paste0("mp_zsummary_", cname, ".pdf")),
         p, width = 8, height = 6)
}

# ---- Combined Results ----

combined <- bind_rows(all_results)
write_tsv(combined, file.path(out_dir, "preservation_summary.tsv"))

combined_plot <- combined[combined$module != "gold", ]
combined_plot$module <- factor(combined_plot$module,
                               levels = unique(combined_plot$module[order(combined_plot$size,
                                                                          decreasing = TRUE)]))

p_all <- ggplot(combined_plot, aes(x = module, y = z_summary, color = preserved)) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", linewidth = 0.3) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "blue", linewidth = 0.3) +
  geom_point(size = 1.8) +
  facet_wrap(~ contrast, ncol = 2) +
  scale_color_manual(values = c("strongly" = "#2166ac",
                                "moderately" = "#fdb863",
                                "not preserved" = "#b2182b")) +
  labs(title = "Module preservation Z-summary across comparisons",
       x = "Module", y = "Z-summary") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "mp_zsummary_all_contrasts.pdf"), p_all,
       width = 12, height = 8)
