# Trait overview — visualize phenotype data across treatment groups

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)

out_fig <- "06_phenotype/output/figures"
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

uninf <- read_csv("06_phenotype/output/trait_summaries/trait_summaries_uninfected.csv",
                   show_col_types = FALSE)
inf <- read_csv("06_phenotype/output/trait_summaries/trait_summaries_infected.csv",
                 show_col_types = FALSE)
indiv <- read_csv("06_phenotype/output/trait_summaries/trait_summaries_individual.csv",
                   show_col_types = FALSE)

# ---- Treatment Group Label ----

make_label <- function(df) {
  df %>% mutate(group = paste0(temp_c_nominal, "C_", conc_ppm, "ppm"))
}

uninf <- make_label(uninf)
inf <- make_label(inf)
indiv <- make_label(indiv)

# ----  Group-Level Heatmaps ----

# excluding sample-count columns
scale_traits <- function(df, id_cols = c("temp_c_nominal", "conc_ppm", "group")) {
  trait_cols <- setdiff(names(df), id_cols)
  trait_cols <- trait_cols[!grepl("^(lht_n_|snail_egg_n)", trait_cols)]
  scaled <- df[, trait_cols] %>%
    mutate(across(everything(), ~ as.numeric(scale(.x))))
  bind_cols(df[, id_cols], scaled)
}

uninf_sc <- scale_traits(uninf)
inf_sc <- scale_traits(inf)

plot_trait_heatmap <- function(df_sc, title) {
  long <- df_sc %>%
    pivot_longer(-c(temp_c_nominal, conc_ppm, group),
                 names_to = "trait", values_to = "z_score")
  ggplot(long, aes(x = group, y = trait, fill = z_score)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, name = "Z-score") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_heat_uninf <- plot_trait_heatmap(uninf_sc, "Uninfected Snail Traits (scaled)")
p_heat_inf <- plot_trait_heatmap(inf_sc, "Infected Snail / Parasite Traits (scaled)")

pdf(file.path(out_fig, "trait_heatmap_combined.pdf"), width = 12, height = 8)
print(p_heat_uninf / p_heat_inf)
dev.off()

# ---- Individual Boxplots ----

plot_boxplots <- function(indiv_df, traits, title) {
  df <- indiv_df %>%
    filter(trait %in% traits, rnaseq_group == TRUE) %>%
    mutate(conc_ppm = factor(conc_ppm))
  ggplot(df, aes(x = factor(temp_c_nominal), y = value, fill = conc_ppm)) +
    geom_boxplot(outlier.size = 0.8) +
    facet_wrap(~ trait, scales = "free_y") +
    scale_fill_manual(values = c("0" = "forestgreen", "0.05" = "firebrick")) +
    labs(title = title, x = "Temperature (C)", y = "Value", fill = "Niclo (ppm)") +
    theme_minimal(base_size = 10)
}

uninf_traits_list <- c("survived", "end_size_mm", "day_since_start",
                        "size_mm", "hatched_yn")
inf_traits_list <- c("total_cercariae", "cercariae_per_day",
                      "prop_dead_24hr", "inf_survived")

p_box_uninf <- plot_boxplots(indiv, uninf_traits_list,
                              "Uninfected Snail Traits by Treatment")
p_box_inf <- plot_boxplots(indiv, inf_traits_list,
                            "Infected / Parasite Traits by Treatment")

pdf(file.path(out_fig, "trait_overview_uninfected.pdf"), width = 10, height = 7)
print(p_box_uninf)
dev.off()

pdf(file.path(out_fig, "trait_overview_infected.pdf"), width = 10, height = 7)
print(p_box_inf)
dev.off()

# ---- Trait Correlation Panels ----




cor_uninf <- uninf %>%
  select(where(is.numeric)) %>%
  select(-matches("^(lht_n_|snail_egg_n)")) %>%
  cor(use = "pairwise.complete.obs")

cor_inf <- inf %>%
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs")

plot_cormat <- function(cormat, title) {
  df <- as.data.frame(as.table(cormat)) %>%
    rename(trait1 = Var1, trait2 = Var2, r = Freq)
  ggplot(df, aes(x = trait1, y = trait2, fill = r)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_cor_uninf <- plot_cormat(cor_uninf, "Uninfected Trait Correlations")
p_cor_inf <- plot_cormat(cor_inf, "Infected Trait Correlations")

pdf(file.path(out_fig, "trait_correlation_panels.pdf"), width = 16, height = 7)
print(p_cor_uninf | p_cor_inf)
dev.off()

# ---- Temperature Response  ----

plot_temp_profile <- function(df, traits, title) {
  long <- df %>%
    select(temp_c_nominal, conc_ppm, all_of(traits)) %>%
    pivot_longer(-c(temp_c_nominal, conc_ppm),
                 names_to = "trait", values_to = "value") %>%
    mutate(conc_ppm = factor(conc_ppm))
  ggplot(long, aes(x = temp_c_nominal, y = value,
                   color = conc_ppm, group = conc_ppm)) +
    geom_point(size = 2) +
    geom_line() +
    facet_wrap(~ trait, scales = "free_y") +
    scale_color_manual(values = c("0" = "forestgreen", "0.05" = "firebrick")) +
    labs(title = title, x = "Temperature (C)", y = "Value",
         color = "Niclo (ppm)") +
    theme_minimal(base_size = 10)
}

key_uninf <- c("lht_survival_rate", "lht_cum_eggs",
               "lht_mean_end_size", "snail_egg_hatch_rate")
key_inf <- c("cerc_total", "cerc_prop_dead_24hr",
             "inf_survival_rate", "para_prop_hatched")

p_temp_uninf <- plot_temp_profile(uninf, key_uninf,
                                   "Uninfected Traits vs Temperature")
p_temp_inf <- plot_temp_profile(inf, key_inf,
                                 "Infected Traits vs Temperature")

pdf(file.path(out_fig, "trait_temp_profiles.pdf"), width = 10, height = 8)
print(p_temp_uninf / p_temp_inf)
dev.off()
