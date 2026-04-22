# Module-trait correlation — correlate WGCNA eigengenes with phenotype traits

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)


out_mt <- "07_trait_integration/output/module_trait"
out_fig <- "07_trait_integration/output/figures"
dir.create(out_mt, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)


me_btru <- read_tsv("03_wgcna/output/btru/btru_module_eigengenes.tsv",
                     show_col_types = FALSE)
me_shae <- read_tsv("03_wgcna/output/shae/shae_module_eigengenes.tsv",
                     show_col_types = FALSE)


me_cols_btru <- grep("^X?[0-9]+$", names(me_btru), value = TRUE)
me_cols_shae <- grep("^X?[0-9]+$", names(me_shae), value = TRUE)

# ---- Split Btru by Infection Status ----

meta <- readRDS("02_de_analysis/output/counts/metadata.rds")
meta_uninf <- meta %>% filter(infection == "uninfected")
meta_inf <- meta %>% filter(infection == "infected")

me_btru_uninf <- me_btru %>% filter(sample %in% meta_uninf$sample)
me_btru_inf <- me_btru %>% filter(sample %in% meta_inf$sample)

# ---- Load Trait Matrices ----

trait_uninf <- read_csv("06_phenotype/output/trait_summaries/trait_summaries_uninfected.csv",
                         show_col_types = FALSE)
trait_inf <- read_csv("06_phenotype/output/trait_summaries/trait_summaries_infected.csv",
                       show_col_types = FALSE)

trait_mat_uninf <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_btru_uninf.rds")
trait_mat_inf <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_btru_inf.rds")
trait_mat_shae <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_shae.rds")

# ---- Align Eigengenes and Traits ----

# Match samples between ME table and trait matrix, return aligned numeric matrices
align_me_trait <- function(me_df, trait_df, me_cols) {
  shared <- intersect(me_df$sample, rownames(trait_df))
  me_sub <- me_df %>%
    filter(sample %in% shared) %>%
    arrange(sample)
  trait_sub <- trait_df[shared[order(shared)], , drop = FALSE]

  me_mat <- as.matrix(me_sub[, me_cols])
  rownames(me_mat) <- me_sub$sample
  trait_mat <- as.matrix(trait_sub)

  list(me = me_mat, trait = trait_mat)
}

aligned_uninf <- align_me_trait(me_btru_uninf, trait_mat_uninf, me_cols_btru)
aligned_inf <- align_me_trait(me_btru_inf, trait_mat_inf, me_cols_btru)
aligned_shae <- align_me_trait(me_shae, trait_mat_shae, me_cols_shae)

# ---- Correlation and P-value Computation ----

# Pearson cor.test for each module-trait pair
cor_and_pval <- function(me_mat, trait_mat) {
  n_mod <- ncol(me_mat)
  n_trait <- ncol(trait_mat)
  cor_mat <- matrix(NA, nrow = n_mod, ncol = n_trait)
  pval_mat <- matrix(NA, nrow = n_mod, ncol = n_trait)
  rownames(cor_mat) <- colnames(me_mat)
  colnames(cor_mat) <- colnames(trait_mat)
  rownames(pval_mat) <- colnames(me_mat)
  colnames(pval_mat) <- colnames(trait_mat)

  for (j in seq_len(n_trait)) {
    for (i in seq_len(n_mod)) {
      ct <- tryCatch(
        cor.test(me_mat[, i], trait_mat[, j], method = "pearson"),
        error = function(e) NULL
      )
      if (!is.null(ct)) {
        cor_mat[i, j] <- ct$estimate
        pval_mat[i, j] <- ct$p.value
      }
    }
  }

  list(cor = cor_mat, pval = pval_mat)
}

res_uninf <- cor_and_pval(aligned_uninf$me, aligned_uninf$trait)
res_inf <- cor_and_pval(aligned_inf$me, aligned_inf$trait)
res_shae <- cor_and_pval(aligned_shae$me, aligned_shae$trait)

# BH adjustment for all module-trait pairs per analysis
bh_adjust_matrix <- function(pval_mat) {
  adj <- matrix(p.adjust(pval_mat, method = "BH"),
                nrow = nrow(pval_mat),
                dimnames = dimnames(pval_mat))
  adj
}

res_uninf$pval_bh <- bh_adjust_matrix(res_uninf$pval)
res_inf$pval_bh <- bh_adjust_matrix(res_inf$pval)
res_shae$pval_bh <- bh_adjust_matrix(res_shae$pval)

# ---- Save Results ----

write.csv(res_uninf$cor, file.path(out_mt, "btru_uninf_module_trait_cor.csv"))
write.csv(res_uninf$pval, file.path(out_mt, "btru_uninf_module_trait_pval_raw.csv"))
write.csv(res_uninf$pval_bh, file.path(out_mt, "btru_uninf_module_trait_pval_bh.csv"))
write.csv(res_inf$cor, file.path(out_mt, "btru_inf_module_trait_cor.csv"))
write.csv(res_inf$pval, file.path(out_mt, "btru_inf_module_trait_pval_raw.csv"))
write.csv(res_inf$pval_bh, file.path(out_mt, "btru_inf_module_trait_pval_bh.csv"))
write.csv(res_shae$cor, file.path(out_mt, "shae_module_trait_cor.csv"))
write.csv(res_shae$pval, file.path(out_mt, "shae_module_trait_pval_raw.csv"))
write.csv(res_shae$pval_bh, file.path(out_mt, "shae_module_trait_pval_bh.csv"))

# ---- Module-Trait Heatmap ----

# Heatmap of module-trait correlations
plot_module_trait_heatmap <- function(cor_mat, pval_mat, title) {
  cor_long <- as.data.frame(as.table(as.matrix(cor_mat)))
  names(cor_long) <- c("module", "trait", "r")

  pval_long <- as.data.frame(as.table(as.matrix(pval_mat)))
  names(pval_long) <- c("module", "trait", "pval")

  df <- left_join(cor_long, pval_long, by = c("module", "trait"))

  df$star <- ifelse(df$pval < 0.001, "***",
               ifelse(df$pval < 0.01, "**",
                 ifelse(df$pval < 0.05, "*", "")))

  ggplot(df, aes(x = trait, y = module, fill = r)) +
    geom_tile(color = "white") +
    geom_text(aes(label = star), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1))
}

# Make sure everything is BH-corrected
p_heat_uninf <- plot_module_trait_heatmap(res_uninf$cor, res_uninf$pval_bh,
                                          "Module-Trait Cor: B. truncatus Uninfected")
p_heat_inf <- plot_module_trait_heatmap(res_inf$cor, res_inf$pval_bh,
                                        "Module-Trait Cor: B. truncatus Infected")
p_heat_shae <- plot_module_trait_heatmap(res_shae$cor, res_shae$pval_bh,
                                         "Module-Trait Cor: S. haematobium")

pdf(file.path(out_fig, "module_trait_heatmap_btru_uninf.pdf"), width = 12, height = 8)
print(p_heat_uninf)
dev.off()

pdf(file.path(out_fig, "module_trait_heatmap_btru_inf.pdf"), width = 12, height = 8)
print(p_heat_inf)
dev.off()

pdf(file.path(out_fig, "module_trait_heatmap_shae.pdf"), width = 12, height = 8)
print(p_heat_shae)
dev.off()

# ---- Top Hits ----

top_hits_df <- function(cor_mat, pval_mat, pval_bh, label, n = 20) {
  cor_long <- as.data.frame(as.table(as.matrix(cor_mat)))
  names(cor_long) <- c("module", "trait", "r")
  pval_long <- as.data.frame(as.table(as.matrix(pval_mat)))
  names(pval_long) <- c("module", "trait", "pval")
  bh_long <- as.data.frame(as.table(as.matrix(pval_bh)))
  names(bh_long) <- c("module", "trait", "pval_bh")

  cor_long %>%
    left_join(pval_long, by = c("module", "trait")) %>%
    left_join(bh_long, by = c("module", "trait")) %>%
    mutate(label = label) %>%
    arrange(pval_bh, pval) %>%
    head(n)
}

top_hits <- bind_rows(
  top_hits_df(res_uninf$cor, res_uninf$pval, res_uninf$pval_bh, "btru_uninf"),
  top_hits_df(res_inf$cor, res_inf$pval, res_inf$pval_bh, "btru_inf"),
  top_hits_df(res_shae$cor, res_shae$pval, res_shae$pval_bh, "shae")
)
write_tsv(top_hits, file.path(out_mt, "module_trait_top_hits.tsv"))

# ---- Parasite as Continuous Variable ----
# Correlate module eigengenes with % Shae reads in infected Btru samples
# Infection proxy - check with Emily - use housekeeping gene??

species_sum <- read_csv("02_de_analysis/output/qc/species_summary.csv",
                         show_col_types = FALSE)
shae_pct <- species_sum %>%
  filter(sample %in% meta_inf$sample) %>%
  select(sample, pct_shae) %>%
  arrange(match(sample, pull(me_btru_inf, sample)))

me_inf_aligned <- me_btru_inf %>%
  filter(sample %in% shae_pct$sample) %>%
  arrange(match(sample, shae_pct$sample))

para_cor <- sapply(me_cols_btru, function(mod) {
  ct <- cor.test(me_inf_aligned[[mod]], shae_pct$pct_shae,
                 method = "spearman", exact = FALSE)
  c(r = unname(ct$estimate), p = ct$p.value)
})
para_df <- data.frame(
  module = me_cols_btru,
  rho = para_cor["r", ],
  pval = para_cor["p", ]
) %>% arrange(pval)

write_tsv(para_df, file.path(out_mt, "btru_parasite_burden_module_cor.tsv"))
