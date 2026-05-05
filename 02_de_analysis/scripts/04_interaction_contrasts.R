library(tidyverse)
library(DESeq2)
library(apeglm)

source("de_functions.R")

# ---- Paths ----

base_dir <- here::here()
btru_dir <- file.path(base_dir, "02_de_analysis/output/btru")
shae_dir <- file.path(base_dir, "02_de_analysis/output/shae")
interact_dir <- file.path(base_dir, "02_de_analysis/output/interactions")
dir.create(interact_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load DESeq2 Objects ----

btru_dds <- readRDS(file.path(btru_dir, "btru_dds.rds"))
shae_dds <- readRDS(file.path(shae_dir, "shae_dds.rds"))


# ---- Btru Interaction Contrasts ----


btru_tempXinf_16 <- results(btru_dds, name = "temp_C16.infectioninfected", alpha = 0.05)
btru_tempXinf_32 <- results(btru_dds, name = "temp_C32.infectioninfected", alpha = 0.05)

btru_tempXniclo_16 <- results(btru_dds, name = "temp_C16.niclo_ppm0.05", alpha = 0.05)
btru_tempXniclo_32 <- results(btru_dds, name = "temp_C32.niclo_ppm0.05", alpha = 0.05)

btru_infXniclo <- results(btru_dds, name = "infectioninfected.niclo_ppm0.05", alpha = 0.05)

btru_tempXinf_16_shr <- lfcShrink(btru_dds, coef = "temp_C16.infectioninfected", type = "apeglm")
btru_tempXinf_32_shr <- lfcShrink(btru_dds, coef = "temp_C32.infectioninfected", type = "apeglm")
btru_tempXniclo_16_shr <- lfcShrink(btru_dds, coef = "temp_C16.niclo_ppm0.05", type = "apeglm")
btru_tempXniclo_32_shr <- lfcShrink(btru_dds, coef = "temp_C32.niclo_ppm0.05", type = "apeglm")
btru_infXniclo_shr <- lfcShrink(btru_dds, coef = "infectioninfected.niclo_ppm0.05", type = "apeglm")

summarize_res(btru_tempXinf_16, "Temp(16):Infection")
summarize_res(btru_tempXinf_32, "Temp(32):Infection")
summarize_res(btru_tempXniclo_16, "Temp(16):Niclosamide")
summarize_res(btru_tempXniclo_32, "Temp(32):Niclosamide")
summarize_res(btru_infXniclo, "Infection:Niclosamide")

btru_temp_32v24 <- results(btru_dds, contrast = c("temp_C", "32", "24"), alpha = 0.05)
summarize_res(btru_temp_32v24, "Temperature (32C vs 24C)")

# ---- Shae Interaction Contrasts ----

shae_tempXniclo_16 <- results(shae_dds, name = "temp_C16.niclo_ppm0.05", alpha = 0.05)
shae_tempXniclo_32 <- results(shae_dds, name = "temp_C32.niclo_ppm0.05", alpha = 0.05)

shae_tempXniclo_16_shr <- lfcShrink(shae_dds, coef = "temp_C16.niclo_ppm0.05", type = "apeglm")
shae_tempXniclo_32_shr <- lfcShrink(shae_dds, coef = "temp_C32.niclo_ppm0.05", type = "apeglm")

summarize_res(shae_tempXniclo_16, "Temp(16):Niclosamide")
summarize_res(shae_tempXniclo_32, "Temp(32):Niclosamide")

shae_temp_32v24 <- results(shae_dds, contrast = c("temp_C", "32", "24"), alpha = 0.05)
summarize_res(shae_temp_32v24, "Temperature (32C vs 24C)")

# ---- Save Results ----

write_results(btru_tempXinf_16, file.path(interact_dir, "btru_tempXinf_16.tsv"))
write_results(btru_tempXinf_32, file.path(interact_dir, "btru_tempXinf_32.tsv"))
write_results(btru_tempXniclo_16, file.path(interact_dir, "btru_tempXniclo_16.tsv"))
write_results(btru_tempXniclo_32, file.path(interact_dir, "btru_tempXniclo_32.tsv"))
write_results(btru_infXniclo, file.path(interact_dir, "btru_infXniclo.tsv"))
write_results(btru_temp_32v24, file.path(interact_dir, "btru_temp32v24.tsv"))

write_results(btru_tempXinf_16_shr, file.path(interact_dir, "btru_tempXinf_16_shrunk.tsv"))
write_results(btru_tempXinf_32_shr, file.path(interact_dir, "btru_tempXinf_32_shrunk.tsv"))
write_results(btru_tempXniclo_16_shr, file.path(interact_dir, "btru_tempXniclo_16_shrunk.tsv"))
write_results(btru_tempXniclo_32_shr, file.path(interact_dir, "btru_tempXniclo_32_shrunk.tsv"))
write_results(btru_infXniclo_shr, file.path(interact_dir, "btru_infXniclo_shrunk.tsv"))

write_results(shae_tempXniclo_16, file.path(interact_dir, "shae_tempXniclo_16.tsv"))
write_results(shae_tempXniclo_32, file.path(interact_dir, "shae_tempXniclo_32.tsv"))
write_results(shae_temp_32v24, file.path(interact_dir, "shae_temp32v24.tsv"))

write_results(shae_tempXniclo_16_shr, file.path(interact_dir, "shae_tempXniclo_16_shrunk.tsv"))
write_results(shae_tempXniclo_32_shr, file.path(interact_dir, "shae_tempXniclo_32_shrunk.tsv"))

# ---- Volcano Plots ----

pdf(file.path(interact_dir, "btru_volcano_tempXinf_32.pdf"), width = 10, height = 8)
print(make_volcano(btru_tempXinf_32, "Btru: Temp(32C) x Infection interaction"))
dev.off()

pdf(file.path(interact_dir, "btru_volcano_tempXniclo_32.pdf"), width = 10, height = 8)
print(make_volcano(btru_tempXniclo_32, "Btru: Temp(32C) x Niclosamide interaction"))
dev.off()

pdf(file.path(interact_dir, "btru_volcano_infXniclo.pdf"), width = 10, height = 8)
print(make_volcano(btru_infXniclo, "Btru: Infection x Niclosamide interaction"))
dev.off()

pdf(file.path(interact_dir, "shae_volcano_tempXniclo_32.pdf"), width = 10, height = 8)
print(make_volcano(shae_tempXniclo_32, "Shae: Temp(32C) x Niclosamide interaction"))
dev.off()

# Combined multi-page PDF
pdf(file.path(interact_dir, "interaction_volcanos.pdf"), width = 10, height = 8)
print(make_volcano(btru_tempXinf_32, "Btru: Temp(32C) x Infection interaction"))
print(make_volcano(btru_tempXniclo_32, "Btru: Temp(32C) x Niclosamide interaction"))
print(make_volcano(btru_infXniclo, "Btru: Infection x Niclosamide interaction"))
print(make_volcano(shae_tempXniclo_32, "Shae: Temp(32C) x Niclosamide interaction"))
dev.off()
