library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ---- Paths ----

out_pt <- "07_trait_integration/output/pathway_trait"
out_fig <- "07_trait_integration/output/figures"
dir.create(out_pt, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

# ---- Load Pathway Scores ----

btru_pw <- read.csv("05_enrichment/output/pathway_activity/btru_pathway_scores.csv",
                     row.names = 1, check.names = FALSE)
shae_pw <- read.csv("05_enrichment/output/pathway_activity/shae_pathway_scores.csv",
                     row.names = 1, check.names = FALSE)

btru_pw_mat <- as.matrix(btru_pw)
shae_pw_mat <- as.matrix(shae_pw)

# ---- Load Traits ----

trait_bu <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_btru_uninf.rds")
trait_bi <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_btru_inf.rds")
trait_sh <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_shae.rds")

meta <- readRDS("02_de_analysis/output/counts/metadata.rds")
meta_uninf <- meta %>% filter(infection == "uninfected")
meta_inf <- meta %>% filter(infection == "infected")

# ---- Correlate Pathways vs Traits ----

cor_pw_trait <- function(pw_mat, trait_df) {
  # Samples are columns in pw_mat but rows in trait_df, so align on the shared sample IDs
  shared <- intersect(colnames(pw_mat), rownames(trait_df))
  pw_sub <- pw_mat[, shared]
  tr_sub <- trait_df[shared, , drop = FALSE]

  # Drop traits that are entirely NA or effectively constant (cor.test fails)
  keep <- sapply(tr_sub, function(x) !all(is.na(x)) && var(x, na.rm = TRUE) > 1e-10)
  tr_sub <- tr_sub[, keep, drop = FALSE]

  n_pw <- nrow(pw_sub)
  n_tr <- ncol(tr_sub)
  cor_mat <- matrix(NA, n_pw, n_tr)
  p_mat <- matrix(NA, n_pw, n_tr)
  rownames(cor_mat) <- rownames(p_mat) <- rownames(pw_sub)
  colnames(cor_mat) <- colnames(p_mat) <- names(tr_sub)

  for (j in seq_len(n_tr)) {
    y <- tr_sub[, j]
    ok <- !is.na(y)
    if (sum(ok) < 5) next
    for (i in seq_len(n_pw)) {
      ct <- cor.test(pw_sub[i, ok], y[ok], method = "pearson")
      cor_mat[i, j] <- ct$estimate
      p_mat[i, j] <- ct$p.value
    }
  }

  list(cor = cor_mat, pval = p_mat)
}

res_bu <- cor_pw_trait(btru_pw_mat, trait_bu)
res_bi <- cor_pw_trait(btru_pw_mat, trait_bi)
res_sh <- cor_pw_trait(shae_pw_mat, trait_sh)

write.csv(res_bu$cor, file.path(out_pt, "btru_uninf_pathway_trait_cor.csv"))
write.csv(res_bi$cor, file.path(out_pt, "btru_inf_pathway_trait_cor.csv"))
write.csv(res_sh$cor, file.path(out_pt, "shae_pathway_trait_cor.csv"))

# ---- Helper: Melt Matrix ----

melt_matrix <- function(mat, row_name = "row", col_name = "col", val_name = "value") {
  # Long-format conversion that matches as.vector(mat)'s column order
  # row names cycle fastest (times = ncol), col names repeat in blocks (each = nrow)
  data.frame(
    row = rep(rownames(mat), times = ncol(mat)),
    col = rep(colnames(mat), each = nrow(mat)),
    value = as.vector(mat),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |> setNames(c(row_name, col_name, val_name))
}

# ---- Top Pathway-Trait Pairs ----

top_pairs <- function(cor_mat, p_mat, n_top = 30) {
  cor_long <- melt_matrix(cor_mat, "pathway", "trait", "r")
  p_long <- melt_matrix(p_mat, "pathway", "trait", "pval")
  left_join(cor_long, p_long, by = c("pathway", "trait")) %>%
    filter(!is.na(pval)) %>%
    arrange(pval) %>%
    head(n_top)
}

top_bu <- top_pairs(res_bu$cor, res_bu$pval)
top_bi <- top_pairs(res_bi$cor, res_bi$pval)
top_sh <- top_pairs(res_sh$cor, res_sh$pval)

write_tsv(bind_rows(mutate(top_bu, context = "btru_uninf"),
                    mutate(top_bi, context = "btru_inf"),
                    mutate(top_sh, context = "shae")),
          file.path(out_pt, "top_pathway_trait_pairs.tsv"))

# ---- Heatmap Plot ----

plot_pw_trait_heatmap <- function(cor_mat, p_mat, n_top = 25, title) {
  # Rank pathways by their single best (smallest) p-value across any trait,
  # Keep the top n_top for plot
  min_p <- apply(p_mat, 1, function(x) min(x, na.rm = TRUE))
  top_pw <- names(sort(min_p))[1:min(n_top, length(min_p))]

  cor_sub <- cor_mat[top_pw, , drop = FALSE]
  p_sub <- p_mat[top_pw, , drop = FALSE]

  cor_long <- melt_matrix(cor_sub, "pathway", "trait", "r")
  p_long <- melt_matrix(p_sub, "pathway", "trait", "pval")

  df <- left_join(cor_long, p_long, by = c("pathway", "trait")) %>%
    mutate(sig = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE ~ ""
    ))

  df$pathway <- substr(as.character(df$pathway), 1, 50)  # Truncate long pathway names so y labels are readable

  ggplot(df, aes(x = trait, y = pathway, fill = r)) +
    geom_tile(color = "grey90") +
    geom_text(aes(label = sig), size = 2.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          axis.text.y = element_text(size = 6))
}

p_bu <- plot_pw_trait_heatmap(res_bu$cor, res_bu$pval, 25,
                               "Btru Uninf: Top Pathway-Trait Correlations")
p_bi <- plot_pw_trait_heatmap(res_bi$cor, res_bi$pval, 25,
                               "Btru Inf: Top Pathway-Trait Correlations")

pdf(file.path(out_fig, "pathway_trait_heatmap_btru.pdf"), width = 14, height = 12)
print(p_bu / p_bi)
dev.off()

p_sh <- plot_pw_trait_heatmap(res_sh$cor, res_sh$pval, 25,
                               "Shae: Top Pathway-Trait Correlations")

pdf(file.path(out_fig, "pathway_trait_heatmap_shae.pdf"), width = 12, height = 8)
print(p_sh)
dev.off()
