library(dplyr)
library(readr)
library(clusterProfiler)
library(SummarizedExperiment)

# ---- Paths ----

base_dir <- here::here()
annot_dir <- file.path(base_dir, "00_setup/03_GO_annotation/output")
out_dir <- file.path(base_dir, "07_trait_integration/output/trait_gsea")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Stuff ----

btru_g2go <- readRDS(file.path(annot_dir, "step05_go_objects/btru_gene2go.Rds"))
btru_go2name <- readRDS(file.path(annot_dir, "step05_go_objects/btru_go2name.Rds"))
btru_kegg <- readRDS(file.path(annot_dir, "step06_kegg_objects/btru_kegg_term2gene.Rds"))

shae_g2go <- readRDS(file.path(annot_dir, "step05_go_objects/shae_gene2go.Rds"))
shae_go2name <- readRDS(file.path(annot_dir, "step05_go_objects/shae_go2name.Rds"))
shae_kegg <- readRDS(file.path(annot_dir, "step06_kegg_objects/shae_kegg_term2gene.Rds"))

btru_vsd <- readRDS("02_de_analysis/output/btru/btru_vsd.rds")
shae_vsd <- readRDS("02_de_analysis/output/shae/shae_vsd.rds")
btru_expr <- assay(btru_vsd)
shae_expr <- assay(shae_vsd)

meta <- readRDS("02_de_analysis/output/counts/metadata.rds")
meta_uninf <- meta %>% filter(infection == "uninfected")
meta_inf <- meta %>% filter(infection == "infected")

trait_uninf <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_btru_uninf.rds")
trait_inf <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_btru_inf.rds")
trait_shae <- readRDS("06_phenotype/output/trait_summaries/trait_matrix_shae.rds")

# ---- Trait-Ranked GSEA Function ----
# Ranks genes by signed -log10(p) of Spearman correlation with each trait

run_trait_gsea <- function(expr_mat, trait_mat, gene2go, kegg_t2g, species, context, out_dir) {
  
  shared <- intersect(colnames(expr_mat), rownames(trait_mat))
  expr <- expr_mat[, shared]
  traits <- trait_mat[shared, ]
  
  go_t2g <- select(gene2go, GO, gene_id)
  kegg_t2g_df <- select(kegg_t2g, KO, gene_id)
  
  summary_list <- list()
  
  for (trait in colnames(traits)) {
    
    rho <- cor(t(expr), traits[, trait], method = "spearman")[, 1]
    t_stat <- rho * sqrt((length(shared) - 2) / (1 - rho^2))
    pval <- 2 * pt(-abs(t_stat), df = length(shared) - 2)
    
    ranked <- sort(sign(rho) * -log10(pval), decreasing = TRUE)
    
    go_res <- GSEA(ranked, TERM2GENE = go_t2g, pvalueCutoff = 0.05)
    kegg_res <- GSEA(ranked, TERM2GENE = kegg_t2g_df, pvalueCutoff = 0.05)
    
    prefix <- file.path(out_dir, paste(species, context, trait, sep = "_"))
    write_tsv(as.data.frame(go_res), paste0(prefix, "_go_gsea.tsv"))
    write_tsv(as.data.frame(kegg_res), paste0(prefix, "_kegg_gsea.tsv"))
    
    summary_list[[trait]] <- data.frame(
      species = species, context = context, trait = trait,
      go_sig = sum(go_res@result$p.adjust < 0.05),
      kegg_sig = sum(kegg_res@result$p.adjust < 0.05)
    )
  }
  
  return(bind_rows(summary_list))
}

# ---- Run Trait GSEA ----

res_bu <- run_trait_gsea(btru_expr, trait_uninf, btru_g2go, btru_go2name,
                          btru_kegg, "btru", "uninf", out_dir)

res_bi <- run_trait_gsea(btru_expr, trait_inf, btru_g2go, btru_go2name,
                          btru_kegg, "btru", "inf", out_dir)

res_sh <- run_trait_gsea(shae_expr, trait_shae, shae_g2go, shae_go2name,
                          shae_kegg, "shae", "inf", out_dir)

# ---- Summary ----

all_res <- bind_rows(res_bu, res_bi, res_sh)
write_tsv(all_res, file.path(out_dir, "trait_gsea_summary.tsv"))
