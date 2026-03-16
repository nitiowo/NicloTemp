# RDA / CCA — constrained ordination linking traits to expression

library(dplyr)
library(readr)
library(vegan)
library(ggplot2)
library(patchwork)
library(SummarizedExperiment)

# ---- Load stuff ----
out_mv <- "08_multivariate/output"
out_fig <- "08_multivariate/output/figures"
dir.create(out_mv, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

btru_vsd <- readRDS("02_de_analysis/output/btru/btru_vsd.rds")
shae_vsd <- readRDS("02_de_analysis/output/shae/shae_vsd.rds")

meta <- readRDS("02_de_analysis/output/counts/metadata.rds")
meta_uninf <- meta %>% filter(infection == "uninfected")
meta_inf <- meta %>% filter(infection == "infected")

trait_uninf <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_btru_uninf.rds")
trait_btru_inf <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_btru_inf.rds")
trait_shae <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_shae.rds")

# Make Expression Matrices 

# Extract assay, subset to samples, pick top variable genes, transpose
prep_expr <- function(vsd_obj, samples, n_top = 2000) {
  mat <- assay(vsd_obj)
  mat <- mat[, colnames(mat) %in% samples, drop = FALSE]
  vars <- apply(mat, 1, var)
  top_genes <- names(sort(vars, decreasing = TRUE))[1:min(n_top, length(vars))]
  t(mat[top_genes, ])
}

expr_btru_uninf <- prep_expr(btru_vsd, meta_uninf$sample)
expr_btru_inf <- prep_expr(btru_vsd, meta_inf$sample)
expr_shae <- prep_expr(shae_vsd, meta_inf$sample)

# Align Expression and Traits

# Intersect rownames, remove zero-variance or all-NA traits, impute NAs
align_for_rda <- function(expr_mat, trait_df) {
  shared <- intersect(rownames(expr_mat), rownames(trait_df))
  expr_mat <- expr_mat[shared, , drop = FALSE]
  trait_df <- trait_df[shared, , drop = FALSE]

  # Remove all-NA columns
  all_na <- sapply(trait_df, function(x) all(is.na(x)))
  trait_df <- trait_df[, !all_na, drop = FALSE]

  # Remove zero-variance columns
  zero_var <- sapply(trait_df, function(x) var(x, na.rm = TRUE) == 0)
  trait_df <- trait_df[, !zero_var, drop = FALSE]

  # Impute remaining NAs with column means
  for (j in seq_len(ncol(trait_df))) {
    na_idx <- is.na(trait_df[, j])
    if (any(na_idx)) {
      trait_df[na_idx, j] <- mean(trait_df[, j], na.rm = TRUE)
    }
  }

  list(expr = expr_mat, traits = trait_df)
}

dat_btru_uninf <- align_for_rda(expr_btru_uninf, trait_uninf)
dat_btru_inf <- align_for_rda(expr_btru_inf, trait_btru_inf)
dat_shae <- align_for_rda(expr_shae, trait_shae)

# ---- Run RDA ----

# Fit RDA, permutation test, compute constrained variance
run_rda <- function(dat, label) {
  rda_fit <- rda(dat$expr ~ ., data = dat$traits)

  perm_global <- anova(rda_fit, permutations = 999)
  perm_axes <- anova(rda_fit, by = "axis", permutations = 999)

  tot_var <- rda_fit$tot.chi
  constr_var <- rda_fit$CCA$tot.chi
  prop_constr <- constr_var / tot_var

  list(
    fit = rda_fit,
    perm_global = perm_global,
    perm_axes = perm_axes,
    constrained_var = prop_constr,
    pval = perm_global[["Pr(>F)"]][1],
    label = label
  )
}

rda_btru_uninf <- run_rda(dat_btru_uninf, "btru_uninf")
rda_btru_inf <- run_rda(dat_btru_inf, "btru_inf")
rda_shae <- run_rda(dat_shae, "shae")

saveRDS(rda_btru_uninf, file.path(out_mv, "rda_btru_uninf.rds"))
saveRDS(rda_btru_inf, file.path(out_mv, "rda_btru_inf.rds"))
saveRDS(rda_shae, file.path(out_mv, "rda_shae.rds"))

# ---- CCA ----

# Fit CCA with permutation test
run_cca <- function(dat, label) {
  fit <- tryCatch(
    cca(dat$expr ~ ., data = dat$traits),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  perm <- anova(fit, permutations = 999)
  pval <- perm[["Pr(>F)"]][1]

  list(fit = fit, pval = pval, label = label)
}

cca_btru_uninf <- run_cca(dat_btru_uninf, "btru_uninf")
cca_btru_inf <- run_cca(dat_btru_inf, "btru_inf")
cca_shae <- run_cca(dat_shae, "shae")

saveRDS(cca_btru_uninf, file.path(out_mv, "cca_btru_uninf.rds"))
saveRDS(cca_btru_inf, file.path(out_mv, "cca_btru_inf.rds"))
saveRDS(cca_shae, file.path(out_mv, "cca_shae.rds"))
