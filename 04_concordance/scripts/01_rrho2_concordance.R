# RRHO2 concordance analysis b/w species
# Compare shared KO numbers (fold change and direction)

library(dplyr)
library(readr)
library(ggplot2)
library(RRHO2)

# ---- Paths ----

base_dir <- here::here()
btru_dir <- file.path(base_dir, "02_de_analysis/output/btru")
shae_dir <- file.path(base_dir, "02_de_analysis/output/shae")
kegg_dir <- file.path(base_dir, "00_setup/03_GO_annotation/output/step06_kegg_objects")
out_dir <- file.path(base_dir, "04_concordance/output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load KEGG Ortholog Mappings ----
btru_ko <- read_tsv(file.path(kegg_dir, "btru_kegg_term2gene.tsv"),
                    show_col_types = FALSE)
shae_ko <- read_tsv(file.path(kegg_dir, "shae_kegg_term2gene.tsv"),
                    show_col_types = FALSE)
shared_kos <- intersect(btru_ko$KO, shae_ko$KO)

# ---- Load Shrunk DE Results ----
# Temp refernce 24 C
load_de <- function(filepath) {
    df <- read_tsv(filepath, show_col_types = FALSE)
    df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
    df$rank_metric <- -log10(df$pvalue) * sign(df$log2FoldChange)
    df
}

btru_t16 <- load_de(file.path(btru_dir, "btru_res_temp16v24_shrunk.tsv"))
btru_t32 <- load_de(file.path(btru_dir, "btru_res_temp32v24_shrunk.tsv"))
btru_niclo <- load_de(file.path(btru_dir, "btru_res_niclosamide_shrunk.tsv"))

shae_t16 <- load_de(file.path(shae_dir, "shae_res_temp16v24_shrunk.tsv"))
shae_t32 <- load_de(file.path(shae_dir, "shae_res_temp32v24_shrunk.tsv"))
shae_niclo <- load_de(file.path(shae_dir, "shae_res_niclosamide_shrunk.tsv"))

# ---- KO-based RRHO2 ----

make_ko_ranked <- function(de_df, ko_map) {
    merged <- inner_join(de_df, ko_map, by = "gene_id")
    merged <- merged[merged$KO %in% shared_kos, ]
    agg <- merged %>%
        group_by(KO) %>%
        summarise(rank_metric = median(rank_metric), .groups = "drop") %>%
        arrange(desc(rank_metric))
    data.frame(gene = agg$KO, rank = agg$rank_metric)
}

run_rrho_ko <- function(btru_df, shae_df, label) {

    btru_ranked <- make_ko_ranked(btru_df, btru_ko)
    shae_ranked <- make_ko_ranked(shae_df, shae_ko)

    # Get identical gene names
    common <- intersect(btru_ranked$gene, shae_ranked$gene)
    btru_ranked <- btru_ranked[btru_ranked$gene %in% common, ]
    shae_ranked <- shae_ranked[shae_ranked$gene %in% common, ]

    n <- nrow(btru_ranked)

    if (n < 50) {
        cat("  Too few matched KOs, skipping RRHO2\n")
        return(NULL)
    }

    step <- max(2, round(sqrt(n) / 4))

    rrho_out <- tryCatch({
        RRHO2_initialize(btru_ranked, shae_ranked,
                         labels = c("Btru", "Shae"),
                         stepsize = step,
                         log10.ind = TRUE)
    }, error = function(e) {
        return(NULL)
    })

    if (is.null(rrho_out)) return(NULL)

    # Save RRHO object
    saveRDS(rrho_out, file.path(out_dir, paste0("rrho_ko_", label, ".rds")))

    # Heatmap
    pdf(file.path(out_dir, paste0("rrho_ko_", label, ".pdf")),
        width = 7, height = 6)
    RRHO2_heatmap(rrho_out,
                  main = paste("RRHO2 (KEGG orthologs):", label))
    dev.off()

    rrho_out
}

rrho_t16 <- run_rrho_ko(btru_t16, shae_t16, "temp_16v24")
rrho_t32 <- run_rrho_ko(btru_t32, shae_t32, "temp_32v24")
rrho_niclo <- run_rrho_ko(btru_niclo, shae_niclo, "niclosamide")

# ---- KO-level Fold-Change Concordance ----

# Compare log2fc between Btru and Shae
ko_fc_comparison <- function(btru_df, shae_df, label) {
    btru_ko_fc <- btru_df %>%
        inner_join(btru_ko, by = "gene_id") %>%
        filter(KO %in% shared_kos) %>%
        group_by(KO) %>%
        summarise(btru_lfc = median(log2FoldChange),
                  btru_padj = min(padj), .groups = "drop")

    shae_ko_fc <- shae_df %>%
        inner_join(shae_ko, by = "gene_id") %>%
        filter(KO %in% shared_kos) %>%
        group_by(KO) %>%
        summarise(shae_lfc = median(log2FoldChange),
                  shae_padj = min(padj), .groups = "drop")

    merged <- inner_join(btru_ko_fc, shae_ko_fc, by = "KO")
    merged$contrast <- label

    if (nrow(merged) < 10) return(list(merged = merged, rho = NA, p = NA))

    test <- cor.test(merged$btru_lfc, merged$shae_lfc, method = "spearman")
    list(merged = merged, rho = test$estimate, p = test$p.value)
}

contrasts <- list(
    "Temp_16v24" = list(btru_t16, shae_t16),
    "Temp_32v24" = list(btru_t32, shae_t32),
    "Niclosamide" = list(btru_niclo, shae_niclo)
)

fc_results <- list()
fc_tables <- list()
for (cname in names(contrasts)) {
    res <- ko_fc_comparison(contrasts[[cname]][[1]],
                            contrasts[[cname]][[2]], cname)
    fc_results[[cname]] <- data.frame(
        contrast = cname,
        rho = round(res$rho, 3),
        p = signif(res$p, 3),
        n_kos = nrow(res$merged)
    )
    fc_tables[[cname]] <- res$merged
}

fc_cor_summary <- bind_rows(fc_results)
print(fc_cor_summary)
write_tsv(fc_cor_summary, file.path(out_dir, "ko_lfc_correlation.tsv"))

# KO FC scatter plots
fc_all <- bind_rows(fc_tables)
write_tsv(fc_all, file.path(out_dir, "ko_lfc_merged.tsv"))

pdf(file.path(out_dir, "ko_fc_scatter.pdf"), width = 12, height = 4)
print(
    ggplot(fc_all, aes(x = btru_lfc, y = shae_lfc)) +
        geom_point(alpha = 0.3, size = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        geom_smooth(method = "lm", se = FALSE, color = "#e34a33",
                    linewidth = 0.5) +
        facet_wrap(~contrast, scales = "free") +
        labs(x = "Btru log2FC (median per KO)",
             y = "Shae log2FC (median per KO)",
             title = "Cross-species fold-change concordance (KEGG orthologs)") +
        theme_minimal()
)
dev.off()

# Direction agreement
for (cname in names(fc_tables)) {
    sig <- fc_tables[[cname]]
    sig <- sig[sig$btru_padj < 0.05 | sig$shae_padj < 0.05, ]
    if (nrow(sig) > 0) {
        agree <- sum(sign(sig$btru_lfc) == sign(sig$shae_lfc))
        cat(cname, ": direction agreement =", agree, "/", nrow(sig),
            "(", round(agree / nrow(sig) * 100, 1), "%)\n")
    }
}
