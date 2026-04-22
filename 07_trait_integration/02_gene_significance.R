library(dplyr)
library(readr)
library(ggplot2)
library(SummarizedExperiment)

# ---- Paths ----

out_gs <- "07_trait_integration/output/gene_significance"
dir.create(out_gs, showWarnings = FALSE, recursive = TRUE)

# ---- Load  Data ----

btru_vsd <- readRDS("02_de_analysis/output/btru/btru_vsd.rds")
shae_vsd <- readRDS("02_de_analysis/output/shae/shae_vsd.rds")

btru_expr <- assay(btru_vsd)
shae_expr <- assay(shae_vsd)

meta <- readRDS("02_de_analysis/output/counts/metadata.rds")
meta_uninf <- meta %>% filter(infection == "uninfected")
meta_inf <- meta %>% filter(infection == "infected")

trait_btru_uninf <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_btru_uninf.rds")
trait_btru_inf <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_btru_inf.rds")
trait_shae <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_shae.rds")

btru_mods <- read_tsv("03_wgcna/output/btru/btru_gene_modules.tsv",
                       show_col_types = FALSE)
shae_mods <- read_tsv("03_wgcna/output/shae/shae_gene_modules.tsv",
                       show_col_types = FALSE)

me_btru <- read_tsv("03_wgcna/output/btru/btru_module_eigengenes.tsv",
                     show_col_types = FALSE)
me_shae <- read_tsv("03_wgcna/output/shae/shae_module_eigengenes.tsv",
                     show_col_types = FALSE)

me_cols_btru <- grep("^X?[0-9]+$", names(me_btru), value = TRUE)
me_cols_shae <- grep("^X?[0-9]+$", names(me_shae), value = TRUE)

# Gene significance function 
# Spearman correlation per gene vs each trait across ONLY shared samples

compute_gene_sig <- function(expr_mat, trait_mat) {
  shared <- intersect(colnames(expr_mat), rownames(trait_mat))
  expr_sub <- expr_mat[, shared, drop = FALSE]
  trait_sub <- as.matrix(trait_mat[shared, , drop = FALSE])

  n_genes <- nrow(expr_sub)
  n_traits <- ncol(trait_sub)
  trait_names <- colnames(trait_sub)

  gs_cor <- matrix(NA, nrow = n_genes, ncol = n_traits)
  gs_pval <- matrix(NA, nrow = n_genes, ncol = n_traits)
  rownames(gs_cor) <- rownames(expr_sub)
  colnames(gs_cor) <- paste0("GS.", trait_names)
  rownames(gs_pval) <- rownames(expr_sub)
  colnames(gs_pval) <- paste0("GS.", trait_names)

  for (j in seq_len(n_traits)) {
    if (sum(!is.na(trait_sub[, j])) < 5) next
    for (i in seq_len(n_genes)) {
      ct <- tryCatch(
        cor.test(expr_sub[i, ], trait_sub[, j], method = "spearman", exact = FALSE),
        error = function(e) NULL
      )
      if (!is.null(ct)) {
        gs_cor[i, j] <- ct$estimate
        gs_pval[i, j] <- ct$p.value
      }
    }
  }

  list(cor = gs_cor, pval = gs_pval)
}

#  Module membership function 
# Pearson correlation per gene vs each module eigengene

compute_module_membership <- function(expr_mat, me_df, me_cols) {
  shared <- intersect(colnames(expr_mat), me_df$sample)
  expr_sub <- expr_mat[, shared, drop = FALSE]
  me_sub <- me_df %>%
    filter(sample %in% shared) %>%
    arrange(match(sample, shared))
  me_mat <- as.matrix(me_sub[, me_cols])

  n_genes <- nrow(expr_sub)
  n_mods <- length(me_cols)

  mm_cor <- matrix(NA, nrow = n_genes, ncol = n_mods)
  rownames(mm_cor) <- rownames(expr_sub)
  colnames(mm_cor) <- paste0("MM.", me_cols)

  for (j in seq_len(n_mods)) {
    mm_cor[, j] <- cor(t(expr_sub), me_mat[, j],
                       use = "pairwise.complete.obs")
  }

  mm_cor
}

# ---- Btru Uninfected Gene Significance ----

gs_bu <- compute_gene_sig(btru_expr, trait_btru_uninf)

gs_bu_df <- data.frame(gene_id = rownames(gs_bu$cor), gs_bu$cor,
                        check.names = FALSE)
gs_bu_pval_df <- data.frame(gene_id = rownames(gs_bu$pval), gs_bu$pval,
                             check.names = FALSE)
names(gs_bu_pval_df)[-1] <- paste0(names(gs_bu_pval_df)[-1], ".pval")

gs_bu_out <- left_join(gs_bu_df, gs_bu_pval_df, by = "gene_id") %>%
  left_join(btru_mods, by = "gene_id")

write_tsv(gs_bu_out, file.path(out_gs, "btru_uninf_gene_significance.tsv"))

# ---- Btru Uninfected Module Membership ----

btru_me_uninf <- me_btru %>% filter(sample %in% meta_uninf$sample)

mm_bu <- compute_module_membership(btru_expr, btru_me_uninf, me_cols_btru)

mm_bu_df <- data.frame(gene_id = rownames(mm_bu), mm_bu, check.names = FALSE)
mm_bu_df <- left_join(mm_bu_df, btru_mods, by = "gene_id")

write_tsv(mm_bu_df, file.path(out_gs, "btru_uninf_module_membership.tsv"))

# ---- Btru Infected Gene Significance ----

gs_bi <- compute_gene_sig(btru_expr, trait_btru_inf)

gs_bi_df <- data.frame(gene_id = rownames(gs_bi$cor), gs_bi$cor,
                        check.names = FALSE)
gs_bi_pval_df <- data.frame(gene_id = rownames(gs_bi$pval), gs_bi$pval,
                             check.names = FALSE)
names(gs_bi_pval_df)[-1] <- paste0(names(gs_bi_pval_df)[-1], ".pval")

gs_bi_out <- left_join(gs_bi_df, gs_bi_pval_df, by = "gene_id") %>%
  left_join(btru_mods, by = "gene_id")

write_tsv(gs_bi_out, file.path(out_gs, "btru_inf_gene_significance.tsv"))

# ---- Btru Infected Module Membership ----

btru_me_inf <- me_btru %>% filter(sample %in% meta_inf$sample)

mm_bi <- compute_module_membership(btru_expr, btru_me_inf, me_cols_btru)

mm_bi_df <- data.frame(gene_id = rownames(mm_bi), mm_bi, check.names = FALSE)
mm_bi_df <- left_join(mm_bi_df, btru_mods, by = "gene_id")

write_tsv(mm_bi_df, file.path(out_gs, "btru_inf_module_membership.tsv"))

# ---- Shae Gene Significance ----

gs_sh <- compute_gene_sig(shae_expr, trait_shae)

gs_sh_df <- data.frame(gene_id = rownames(gs_sh$cor), gs_sh$cor,
                        check.names = FALSE)
gs_sh_pval_df <- data.frame(gene_id = rownames(gs_sh$pval), gs_sh$pval,
                             check.names = FALSE)
names(gs_sh_pval_df)[-1] <- paste0(names(gs_sh_pval_df)[-1], ".pval")

gs_sh_out <- left_join(gs_sh_df, gs_sh_pval_df, by = "gene_id") %>%
  left_join(shae_mods, by = "gene_id")

write_tsv(gs_sh_out, file.path(out_gs, "shae_gene_significance.tsv"))

# ---- Shae Module Membership ----

shae_me <- me_shae %>% filter(sample %in% meta_inf$sample)

mm_sh <- compute_module_membership(shae_expr, shae_me, me_cols_shae)

mm_sh_df <- data.frame(gene_id = rownames(mm_sh), mm_sh, check.names = FALSE)
mm_sh_df <- left_join(mm_sh_df, shae_mods, by = "gene_id")

write_tsv(mm_sh_df, file.path(out_gs, "shae_module_membership.tsv"))
