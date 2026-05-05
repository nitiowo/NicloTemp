library(mediation)

# ---- Paths ----

base_dir <- here::here()
out_dir <- file.path(base_dir, "08_multivariate/output/mediation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Data ----

meta <- readRDS(file.path(base_dir, "02_de_analysis/output/counts/metadata.rds"))

me_btru <- read_tsv(file.path(base_dir, "03_wgcna/output/btru/btru_module_eigengenes.tsv"),
                    show_col_types = FALSE)

traits_inf   <- readRDS(file.path(base_dir, "06_phenotype/output/trait_summaries/trait_matrix_btru_inf.rds"))
traits_uninf <- readRDS(file.path(base_dir, "06_phenotype/output/trait_summaries/trait_matrix_btru_uninf.rds"))

# ---- Define Treatment-Module-Trait Combos ----
# Pick combos with strong module-trait correlations

trios <- tribble(
  ~treatment,  ~module, ~trait,                   ~subset,
  "temp_C",    "6",     "para_prop_hatched",       "infected",
  "temp_C",    "6",     "para_miracidia_per_egg",  "infected",
  "temp_C",    "8",     "cerc_total",              "infected",
  "temp_C",    "8",     "cerc_per_day",            "infected",
  "niclo_ppm", "5",     "inf_survival_rate",       "infected",
  "niclo_ppm", "5",     "inf_shed_cercariae",      "infected",
  "temp_C",    "1",     "inf_n_eggs_total",        "infected",
  "temp_C",    "3",     "cerc_prop_dead_24hr",     "infected",
  "temp_C",    "15",    "cerc_shed_days",          "infected",
  "niclo_ppm", "12",    "inf_n_egg_masses",        "infected"
)

# ---- Mediation ----

# Mediator model: module ~ treatment; outcome model: trait ~ treatment + module
run_mediation <- function(dat, n_sims = 1000) {
  mod_med <- lm(mod_eig ~ trt, data = dat)
  mod_out <- lm(trait_val ~ trt + mod_eig, data = dat)
  med_res <- mediate(mod_med, mod_out,
                     treat = "trt",
                     mediator = "mod_eig",
                     boot = TRUE,
                     sims = n_sims)
  tibble(
    acme       = med_res$d0,
    acme_ci_lo = med_res$d0.ci[1],
    acme_ci_hi = med_res$d0.ci[2],
    acme_p     = med_res$d0.p,
    ade        = med_res$z0,
    ade_ci_lo  = med_res$z0.ci[1],
    ade_ci_hi  = med_res$z0.ci[2],
    ade_p      = med_res$z0.p,
    total          = med_res$tau.coef,
    prop_mediated  = med_res$n0,
    prop_med_p     = med_res$n0.p
  )
}

results <- list()

for (i in seq_len(nrow(trios))) {
  tr <- trios[i, ]

  if (tr$subset == "infected") {
    samps     <- meta$sample[meta$infection == "infected"]
    trait_mat <- traits_inf
  } else {
    samps     <- meta$sample[meta$infection == "uninfected"]
    trait_mat <- traits_uninf
  }

  me_sub <- me_btru %>%
    dplyr::filter(sample %in% samps) %>%
    dplyr::select(sample, all_of(tr$module))

  trait_vals <- trait_mat[me_sub$sample, tr$trait, drop = TRUE]

  dat <- me_sub %>%
    mutate(
      trait_val = trait_vals,
      trt = if (tr$treatment == "temp_C") {
        as.numeric(as.character(meta$temp_C[match(sample, meta$sample)]))
      } else {
        as.numeric(as.character(meta$niclo_ppm[match(sample, meta$sample)]))
      }
    ) %>%
    dplyr::select(sample, mod_eig = !!tr$module, trait_val, trt) %>%
    drop_na()

  if (nrow(dat) < 8) next

  res <- tryCatch(
    run_mediation(dat),
    error = function(e) NULL
  )

  if (!is.null(res)) results[[i]] <- bind_cols(tr, res)
}

med_df <- bind_rows(results)
write_tsv(med_df, file.path(out_dir, "mediation_results.tsv"))

# ---- Forest Plot of ACME ----

if (nrow(med_df) > 0) {
  med_df <- med_df %>%
    mutate(label = paste(treatment, "\u2192", module, "\u2192", trait))

  p <- ggplot(med_df, aes(x = acme, y = reorder(label, acme))) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = acme_ci_lo, xmax = acme_ci_hi), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = "Average Causal Mediation Effect (ACME)",
         y = NULL,
         title = "Mediation: Treatment \u2192 Module \u2192 Trait") +
    theme_minimal(base_size = 11)

  invisible(ggsave(file.path(out_dir, "mediation_forest.pdf"), p, width = 10, height = 6))
}
