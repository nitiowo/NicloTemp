# Summarize cleaned phenotype data to treatment-group level

#### Check cleaning logic

library(dplyr)
library(readr)
library(tidyr)

pheno_dir <- "Data/Phenotype_data/cleaned"
out_dir <- "06_phenotype/output/trait_summaries"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
meta <- readRDS("02_de_analysis/output/counts/metadata.rds")

meta_uninf <- meta %>% filter(infection == "uninfected")
meta_inf <- meta %>% filter(infection == "infected")

lht_ms <- read_csv(file.path(pheno_dir, "lht_mortality_sizing.csv"), show_col_types = FALSE)
lht_egg <- read_csv(file.path(pheno_dir, "lht_egglaying.csv"), show_col_types = FALSE)
lht_gr <- read_csv(file.path(pheno_dir, "lht_growthrate.csv"), show_col_types = FALSE)
snail_hatch <- read_csv(file.path(pheno_dir, "snail_egghatching.csv"), show_col_types = FALSE)
cerc_em_sum <- read_csv(file.path(pheno_dir, "cerc_emergence_summary.csv"), show_col_types = FALSE)
cerc_fec <- read_csv(file.path(pheno_dir, "cerc_fecundity.csv"), show_col_types = FALSE)
cerc_mort <- read_csv(file.path(pheno_dir, "cerc_mortality.csv"), show_col_types = FALSE)
cerc_mort_inf <- read_csv(file.path(pheno_dir, "cerc_mortality_infsnail.csv"), show_col_types = FALSE)
para_hatch <- read_csv(file.path(pheno_dir, "parasite_egghatching.csv"), show_col_types = FALSE)
susc_mort <- read_csv(file.path(pheno_dir, "snail_susceptibility_mortality.csv"), show_col_types = FALSE)

# ---- Group Summarization Functions ----

# Group key: temp_c_nominal x conc_ppm, only rnaseq_group == TRUE
group_mean <- function(df, ..., .filter_rnaseq = TRUE) {
  if (.filter_rnaseq) df <- df %>% filter(rnaseq_group == TRUE)
  df %>%
    group_by(temp_c_nominal, conc_ppm, ...) %>%
    summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
              n = n(), .groups = "drop")
}

group_se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# ---- Uninfected BBtru Traits ----

uninf_mortality <- lht_ms %>%
  filter(rnaseq_group == TRUE) %>%
  group_by(temp_c_nominal, conc_ppm) %>%
  summarize(
    lht_survival_rate = mean(survived, na.rm = TRUE),
    lht_median_day_death = median(day_since_start[!survived], na.rm = TRUE),
    lht_mean_start_size = mean(start_size_mm, na.rm = TRUE),
    lht_mean_end_size = mean(end_size_mm, na.rm = TRUE),
    lht_size_change = mean(end_size_mm - start_size_mm, na.rm = TRUE),
    lht_n_snails = n(),
    .groups = "drop"
  )

uninf_egglaying <- lht_egg %>%
  filter(rnaseq_group == TRUE) %>%
  group_by(temp_c_nominal, conc_ppm, replicate) %>%
  summarize(
    cum_eggs = sum(eggs_total, na.rm = TRUE),
    cum_clutches = sum(n_clutches, na.rm = TRUE),
    mean_eggs_per_snail = mean(eggs_per_snail, na.rm = TRUE),
    n_dates = n(),
    .groups = "drop"
  ) %>%
  group_by(temp_c_nominal, conc_ppm) %>%
  summarize(
    lht_cum_eggs = mean(cum_eggs, na.rm = TRUE),
    lht_cum_clutches = mean(cum_clutches, na.rm = TRUE),
    lht_eggs_per_snail = mean(mean_eggs_per_snail, na.rm = TRUE),
    .groups = "drop"
  )

uninf_growth <- lht_gr %>%
  filter(rnaseq_group == TRUE) %>%
  group_by(temp_c_nominal, conc_ppm) %>%
  summarize(
    lht_mean_size_mm = mean(size_mm, na.rm = TRUE),
    lht_size_sd = sd(size_mm, na.rm = TRUE),
    lht_n_measured = n(),
    .groups = "drop"
  )

uninf_egghatching <- snail_hatch %>%
  filter(rnaseq_group == TRUE) %>%
  group_by(temp_c_nominal, conc_ppm) %>%
  summarize(
    snail_egg_hatch_rate = mean(hatched_yn, na.rm = TRUE),
    snail_egg_time_to_hatch = mean(time_to_hatch_days, na.rm = TRUE),
    snail_egg_n = n(),
    .groups = "drop"
  )

uninf_traits <- uninf_mortality %>%
  left_join(uninf_egglaying, by = c("temp_c_nominal", "conc_ppm")) %>%
  left_join(uninf_growth, by = c("temp_c_nominal", "conc_ppm")) %>%
  left_join(uninf_egghatching, by = c("temp_c_nominal", "conc_ppm"))

write_csv(uninf_traits, file.path(out_dir, "trait_summaries_uninfected.csv"))

# ---- Infected Btru / Shae Traits ----

inf_cerc_em <- cerc_em_sum %>%
  filter(rnaseq_group == TRUE) %>%
  group_by(temp_c_nominal, conc_ppm) %>%
  summarize(
    cerc_total = mean(total_cercariae, na.rm = TRUE),
    cerc_per_day = mean(mean_cercariae_per_day, na.rm = TRUE),
    cerc_shed_days = mean(n_shed_days, na.rm = TRUE),
    .groups = "drop"
  )

inf_cerc_mort <- cerc_mort %>%
  filter(rnaseq_group == TRUE) %>%
  mutate(prop_dead_24hr = n_dead_24hr / 6) %>%
  group_by(temp_c_nominal, conc_ppm) %>%
  summarize(
    cerc_prop_dead_24hr = mean(prop_dead_24hr, na.rm = TRUE),
    cerc_n_dead_24hr = mean(n_dead_24hr, na.rm = TRUE),
    .groups = "drop"
  )

inf_fec <- cerc_fec %>%
  filter(rnaseq_group == TRUE) %>%
  group_by(temp_c_nominal, conc_ppm) %>%
  summarize(
    inf_n_egg_masses = mean(n_egg_masses, na.rm = TRUE),
    inf_n_eggs_total = mean(n_eggs_total, na.rm = TRUE),
    .groups = "drop"
  )

inf_mort <- cerc_mort_inf %>%
  filter(rnaseq_group == TRUE) %>%
  group_by(temp_c_nominal, conc_ppm) %>%
  summarize(
    inf_survival_rate = mean(survived, na.rm = TRUE),
    inf_shed_cercariae = mean(shed_cerc, na.rm = TRUE),
    inf_laid_eggs = mean(laid_eggs, na.rm = TRUE),
    inf_sporocyst_rate = mean(has_sporocysts, na.rm = TRUE),
    inf_mean_cerc_per_hr = mean(mean_cerc_per_hour, na.rm = TRUE),
    inf_start_size = mean(start_size_mm, na.rm = TRUE),
    .groups = "drop"
  )

inf_susc <- susc_mort %>%
  filter(rnaseq_group == TRUE) %>%
  group_by(temp_c_nominal, conc_ppm) %>%
  summarize(
    susc_survival_rate = mean(survived, na.rm = TRUE),
    susc_mean_day_death = mean(day_since_start[!survived], na.rm = TRUE),
    .groups = "drop"
  )

inf_para_hatch <- para_hatch %>%
  filter(rnaseq_group == TRUE, timepoint == "4hr") %>%
  group_by(temp_c_nominal, conc_ppm) %>%
  summarize(
    para_prop_hatched = mean(prop_hatched, na.rm = TRUE),
    para_miracidia_per_egg = mean(miracidia_per_egg, na.rm = TRUE),
    .groups = "drop"
  )

inf_traits <- inf_cerc_em %>%
  left_join(inf_cerc_mort, by = c("temp_c_nominal", "conc_ppm")) %>%
  left_join(inf_fec, by = c("temp_c_nominal", "conc_ppm")) %>%
  left_join(inf_mort, by = c("temp_c_nominal", "conc_ppm")) %>%
  left_join(inf_susc, by = c("temp_c_nominal", "conc_ppm")) %>%
  left_join(inf_para_hatch, by = c("temp_c_nominal", "conc_ppm"))

write_csv(inf_traits, file.path(out_dir, "trait_summaries_infected.csv"))

# ---- Individual-Level  ----

indiv_list <- list(
  lht_ms %>%
    transmute(source = "lht_mortality", temp_c_nominal, conc_ppm,
              replicate, rnaseq_group, trait = "survived", value = as.numeric(survived)),
  lht_ms %>%
    transmute(source = "lht_mortality", temp_c_nominal, conc_ppm,
              replicate, rnaseq_group, trait = "end_size_mm", value = end_size_mm),
  lht_ms %>%
    filter(!survived) %>%
    transmute(source = "lht_mortality", temp_c_nominal, conc_ppm,
              replicate, rnaseq_group, trait = "day_since_start", value = day_since_start),
  lht_gr %>%
    transmute(source = "lht_growth", temp_c_nominal, conc_ppm,
              replicate, rnaseq_group, trait = "size_mm", value = size_mm),
  snail_hatch %>%
    transmute(source = "snail_egghatching", temp_c_nominal, conc_ppm,
              replicate, rnaseq_group, trait = "hatched_yn", value = as.numeric(hatched_yn)),
  cerc_em_sum %>%
    transmute(source = "cerc_emergence", temp_c_nominal, conc_ppm,
              replicate, rnaseq_group, trait = "total_cercariae", value = total_cercariae),
  cerc_em_sum %>%
    transmute(source = "cerc_emergence", temp_c_nominal, conc_ppm,
              replicate, rnaseq_group, trait = "cercariae_per_day", value = mean_cercariae_per_day),
  cerc_mort %>%
    transmute(source = "cerc_mortality", temp_c_nominal, conc_ppm,
              replicate, rnaseq_group, trait = "prop_dead_24hr", value = n_dead_24hr / 6),
  cerc_mort_inf %>%
    transmute(source = "cerc_infsnail_mort", temp_c_nominal, conc_ppm,
              replicate, rnaseq_group, trait = "inf_survived", value = as.numeric(survived))
)

indiv_all <- bind_rows(indiv_list)
write_csv(indiv_all, file.path(out_dir, "trait_summaries_individual.csv"))

# ---- Sample-Aligned Trait Matrices ----

# Each sample gets trait values based on its temp_C x niclo_ppm group
align_traits <- function(meta_sub, trait_df) {
  joined <- meta_sub %>%
    mutate(temp_C = as.numeric(as.character(temp_C)),
           niclo_ppm = as.numeric(as.character(niclo_ppm))) %>%
    left_join(trait_df,
              by = c("temp_C" = "temp_c_nominal", "niclo_ppm" = "conc_ppm")) %>%
    arrange(sample)
  trait_cols <- setdiff(names(trait_df), c("temp_c_nominal", "conc_ppm"))
  mat <- as.data.frame(joined[, trait_cols, drop = FALSE])
  rownames(mat) <- joined$sample
  mat
}

trait_mat_btru_uninf <- align_traits(meta_uninf, uninf_traits)
trait_mat_btru_inf <- align_traits(meta_inf, inf_traits)
trait_mat_shae <- align_traits(meta_inf, inf_traits)

saveRDS(trait_mat_btru_uninf, file.path(out_dir, "trait_matrix_btru_uninf.rds"))
saveRDS(trait_mat_btru_inf, file.path(out_dir, "trait_matrix_btru_inf.rds"))
saveRDS(trait_mat_shae, file.path(out_dir, "trait_matrix_shae.rds"))
