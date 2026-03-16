# sPLS linking traits to expression

library(dplyr)
library(readr)
library(mixOmics)
library(ggplot2)
library(patchwork)
library(SummarizedExperiment)

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

## Make master functions source script???
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

# ---- Clean Traits ----

# Intersect samples, remove near-zero-variance columns, impute NAs
clean_traits <- function(trait_df, samples) {
  shared <- intersect(samples, rownames(trait_df))
  trait_df <- trait_df[shared, , drop = FALSE]

  # Remove all-NA columns
  all_na <- sapply(trait_df, function(x) all(is.na(x)))
  trait_df <- trait_df[, !all_na, drop = FALSE]

  # Remove near-zero-variance columns
  nzv <- sapply(trait_df, function(x) {
    v <- var(x, na.rm = TRUE)
    is.na(v) || v < 1e-10
  })
  trait_df <- trait_df[, !nzv, drop = FALSE]

  # Impute remaining NAs with column means
  for (j in seq_len(ncol(trait_df))) {
    na_idx <- is.na(trait_df[, j])
    if (any(na_idx)) {
      trait_df[na_idx, j] <- mean(trait_df[, j], na.rm = TRUE)
    }
  }

  as.matrix(trait_df)
}

trait_clean_uninf <- clean_traits(trait_uninf, rownames(expr_btru_uninf))
trait_clean_btru_inf <- clean_traits(trait_btru_inf, rownames(expr_btru_inf))
trait_clean_shae <- clean_traits(trait_shae, rownames(expr_shae))

# ---- Fit sPLS ----

fit_spls <- function(X, Y, label, keepX = 100) {
  ncomp <- min(nrow(X) - 1, ncol(Y), 3)
  fit <- spls(X, Y, ncomp = ncomp,
              keepX = rep(keepX, ncomp),
              mode = "regression")
  fit
}

spls_btru_uninf <- fit_spls(expr_btru_uninf, trait_clean_uninf, "btru_uninf")
spls_btru_inf <- fit_spls(expr_btru_inf, trait_clean_btru_inf, "btru_inf")
spls_shae <- fit_spls(expr_shae, trait_clean_shae, "shae")

saveRDS(spls_btru_uninf, file.path(out_mv, "spls_btru_uninf.rds"))
saveRDS(spls_btru_inf, file.path(out_mv, "spls_btru_inf.rds"))
saveRDS(spls_shae, file.path(out_mv, "spls_shae.rds"))

# ---- Extract Selected Genes ----

# Get selected variables per component into a data.frame
extract_selected <- function(spls_fit, label) {
  ncomp <- spls_fit$ncomp
  res_list <- lapply(seq_len(ncomp), function(comp) {
    sv <- selectVar(spls_fit, comp = comp)
    df <- data.frame(
      gene_id = rownames(sv$value),
      loading = sv$value[, 1],
      component = comp,
      organism = label,
      stringsAsFactors = FALSE
    )
    df
  })
  do.call(rbind, res_list)
}

sel_btru_uninf <- extract_selected(spls_btru_uninf, "btru_uninf")
sel_btru_inf <- extract_selected(spls_btru_inf, "btru_inf")
sel_shae <- extract_selected(spls_shae, "shae")

# ---- sPLS Plots ----

plot_spls_samples <- function(spls_fit, meta_sub, title) {
  scores <- spls_fit$variates$X[, 1:2]
  df <- data.frame(
    comp1 = scores[, 1],
    comp2 = scores[, 2],
    sample = rownames(scores)
  )
  df <- merge(df, meta_sub, by = "sample")
  df$temp_C <- as.numeric(as.character(df$temp_C))
  df$niclo_ppm <- as.numeric(as.character(df$niclo_ppm))

  ggplot(df, aes(x = comp1, y = comp2)) +
    geom_point(aes(color = factor(temp_C),
                   shape = factor(niclo_ppm)), size = 3) +
    scale_color_manual(
      values = c("16" = "#377EB8", "24" = "#4DAF4A", "32" = "#E41A1C"),
      name = "Temp (C)"
    ) +
    scale_shape_manual(
      values = c("0" = 16, "0.05" = 17),
      name = "Niclo (ppm)"
    ) +
    labs(title = title, x = "sPLS Comp 1", y = "sPLS Comp 2") +
    theme_minimal(base_size = 10)
}

p_samp_uninf <- plot_spls_samples(spls_btru_uninf, meta_uninf, "B. truncatus uninf")
p_samp_inf <- plot_spls_samples(spls_btru_inf, meta_inf, "B. truncatus inf")
p_samp_shae <- plot_spls_samples(spls_shae, meta_inf, "S. haematobium")

pdf(file.path(out_fig, "spls_sample_plot.pdf"), width = 18, height = 6)
print(p_samp_uninf | p_samp_inf | p_samp_shae)
dev.off()

