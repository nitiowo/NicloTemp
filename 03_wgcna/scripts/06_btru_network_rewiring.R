library(tidyverse)
library(WGCNA)

# ---- Paths ----

base_dir <- here::here()
wgcna_base <- file.path(base_dir, "03_wgcna/output")
out_dir <- file.path(base_dir, "03_wgcna/output/rewiring")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Eigengene Matrices ----

load_mes <- function(group_name) {
  me_file <- file.path(wgcna_base, paste0("btru_", group_name),
                       paste0("btru_", group_name, "_MEs.rds"))
  if (!file.exists(me_file)) return(NULL)
  MEs <- readRDS(me_file)
  MEs[, colnames(MEs) != "0", drop = FALSE]  # Drop grey (unassigned)
}

groups <- c("infected", "uninfected", "niclo_treated", "niclo_control",
            "temp_16", "temp_24", "temp_32")

me_list <- lapply(groups, load_mes)
names(me_list) <- groups

# ---- Compute Eigengene-Network Correlation Matrices ----

me_net <- lapply(me_list, function(MEs) {
  if (is.null(MEs) || ncol(MEs) < 2) return(NULL)
  cor(MEs, use = "pairwise.complete.obs")
})

# ---- Pairwise Delta-Correlation ----

pairwise_comps <- list(
  infection     = c("uninfected", "infected"),
  niclosamide   = c("niclo_control", "niclo_treated"),
  temp_16_vs_24 = c("temp_24", "temp_16"),
  temp_32_vs_24 = c("temp_24", "temp_32"),
  temp_16_vs_32 = c("temp_32", "temp_16")
)

delta_rows <- lapply(names(pairwise_comps), function(comp) {
  pair <- pairwise_comps[[comp]]
  mat_ref  <- me_net[[pair[1]]]
  mat_test <- me_net[[pair[2]]]
  if (is.null(mat_ref) || is.null(mat_test)) return(NULL)

  common_mods <- intersect(colnames(mat_ref), colnames(mat_test))
  if (length(common_mods) < 2) return(NULL)

  mat_ref  <- mat_ref[common_mods, common_mods]
  mat_test <- mat_test[common_mods, common_mods]
  diff_mat <- mat_test - mat_ref

  idx <- which(upper.tri(diff_mat), arr.ind = TRUE)
  tibble(
    comparison = comp,
    module_a   = rownames(diff_mat)[idx[, 1]],
    module_b   = colnames(diff_mat)[idx[, 2]],
    cor_ref    = mat_ref[idx],
    cor_test   = mat_test[idx],
    delta_r    = diff_mat[idx]
  ) %>%
    arrange(desc(abs(delta_r)))
})

delta_tbl <- bind_rows(delta_rows)
write_tsv(delta_tbl, file.path(out_dir, "btru_me_pairwise_delta.tsv"))

# ---- Rewiring Summary Table ----

delta_summary <- delta_tbl %>%
  group_by(comparison) %>%
  summarise(
    mean_abs_delta  = round(mean(abs(delta_r), na.rm = TRUE), 3),
    n_pairs         = n(),
    n_rewired_gt02  = sum(abs(delta_r) > 0.2, na.rm = TRUE),
    .groups = "drop"
  )
write_tsv(delta_summary, file.path(out_dir, "btru_me_rewiring_summary.tsv"))

# ---- Summary Barplot ----

comp_order <- delta_summary %>%
  arrange(desc(mean_abs_delta)) %>%
  pull(comparison)

p_bar <- ggplot(delta_summary,
                aes(x = factor(comparison, levels = comp_order),
                    y = mean_abs_delta)) +
  geom_col(fill = "#4575b4", width = 0.6) +
  geom_text(aes(label = paste0(n_rewired_gt02, " pairs\n> 0.2")),
            vjust = -0.3, size = 2.5) +
  labs(x = NULL, y = "Mean |\u0394r|",
       title = "Eigengene network rewiring magnitude by condition") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

invisible(ggsave(file.path(out_dir, "btru_rewiring_barplot.pdf"),
                 p_bar, width = 6, height = 4))

# ---- Heatmap ----

top_delta <- delta_tbl %>%
  group_by(comparison) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(pair = paste0("M", module_a, "-M", module_b))

p_heat <- ggplot(top_delta, aes(x = comparison, y = pair, fill = delta_r)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, name = "\u0394r") +
  labs(title = "Top module-pair rewiring across Btru subgroup networks",
       x = NULL, y = "Module pair") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

invisible(ggsave(file.path(out_dir, "btru_rewiring_heatmap.pdf"),
                 p_heat, width = 8, height = 6))
