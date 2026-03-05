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
wgcna_btru <- file.path(base_dir, "03_wgcna/output/btru")
wgcna_shae <- file.path(base_dir, "03_wgcna/output/shae")
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

# ---- DE Landscape Comparison ----
# Compare proportion of DE genes in each direction per species

classify_de <- function(de_df, label) {
    de_df$direction <- "NS"
    de_df$direction[de_df$padj < 0.05 & de_df$log2FoldChange > 0] <- "Up"
    de_df$direction[de_df$padj < 0.05 & de_df$log2FoldChange < 0] <- "Down"
    counts <- as.data.frame(table(direction = de_df$direction))
    counts$pct <- round(counts$Freq / sum(counts$Freq) * 100, 1)
    counts$species <- label
    counts
}

contrasts <- list(
  "Temp_16v24" = list(btru_t16, shae_t16),
  "Temp_32v24" = list(btru_t32, shae_t32),
  "Niclosamide" = list(btru_niclo, shae_niclo)
)

landscape_list <- list()
for (cname in names(contrasts)) {
    b <- classify_de(contrasts[[cname]][[1]], "Btru")
    s <- classify_de(contrasts[[cname]][[2]], "Shae")
    both <- rbind(b, s)
    both$contrast <- cname
    landscape_list[[cname]] <- both
}
landscape <- bind_rows(landscape_list)
write_tsv(landscape, file.path(out_dir, "de_landscape.tsv"))

# Stacked barplot
landscape$direction <- factor(landscape$direction,
                              levels = c("Down", "NS", "Up"))

pdf(file.path(out_dir, "de_landscape.pdf"), width = 10, height = 5)
print(
    ggplot(landscape, aes(x = species, y = pct, fill = direction)) +
        geom_col(position = "stack") +
        facet_wrap(~contrast) +
        scale_fill_manual(values = c("Down" = "#2166ac",
                                     "NS" = "grey80",
                                     "Up" = "#b2182b")) +
        labs(y = "% of tested genes", x = NULL, fill = "DE status",
             title = "DE landscape: Btru vs Shae") +
        theme_minimal()
)
dev.off()

# ---- KO-level Fold-Change Concordance ----
# For shared KOs, compare log2FC between Btru and Shae

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

# ----WGCNA Module Eigengene Cross-Correlation ----

btru_me_file <- file.path(wgcna_btru, "btru_MEs.rds")
shae_me_file <- file.path(wgcna_shae, "shae_MEs.rds")

if (file.exists(btru_me_file) && file.exists(shae_me_file)) {
    btru_MEs <- readRDS(btru_me_file)
    shae_MEs <- readRDS(shae_me_file)

    # MEs already have "ME" prefix stripped (from WGCNA script)
    # Add prefix back for clarity in cross-species comparison
    colnames(btru_MEs) <- paste0("Btru_", colnames(btru_MEs))
    colnames(shae_MEs) <- paste0("Shae_", colnames(shae_MEs))

    shared_samples <- intersect(rownames(btru_MEs), rownames(shae_MEs))

    if (length(shared_samples) >= 5) {
        btru_sub <- btru_MEs[shared_samples, , drop = FALSE]
        shae_sub <- shae_MEs[shared_samples, , drop = FALSE]

        # Drop grey modules
        btru_sub <- btru_sub[, !grepl("grey$", colnames(btru_sub)),
                             drop = FALSE]
        shae_sub <- shae_sub[, !grepl("grey$", colnames(shae_sub)),
                             drop = FALSE]

        if (ncol(btru_sub) >= 1 && ncol(shae_sub) >= 1) {
            cross_cor <- cor(btru_sub, shae_sub,
                             use = "pairwise.complete.obs")

            clust_rows <- nrow(cross_cor) > 1
            clust_cols <- ncol(cross_cor) > 1

            pdf(file.path(out_dir, "eigengene_cross_cor.pdf"),
                width = 8, height = 6)
            pheatmap::pheatmap(
                cross_cor,
                main = "Btru vs Shae module eigengene correlation (18 infected)",
                color = colorRampPalette(
                    c("#2166ac", "white", "#b2182b"))(100),
                breaks = seq(-1, 1, length.out = 101),
                display_numbers = TRUE,
                number_format = "%.2f",
                fontsize_number = 7,
                cluster_rows = clust_rows,
                cluster_cols = clust_cols
            )
            dev.off()

            # P-values for cross-correlation
            n_samp <- length(shared_samples)
            cross_pval <- matrix(NA, nrow(cross_cor), ncol(cross_cor),
                                 dimnames = dimnames(cross_cor))
            for (i in seq_len(nrow(cross_cor))) {
                for (j in seq_len(ncol(cross_cor))) {
                    test <- cor.test(btru_sub[, i], shae_sub[, j])
                    cross_pval[i, j] <- test$p.value
                }
            }

            write.csv(cross_cor,
                      file.path(out_dir, "eigengene_cross_cor.csv"))
            write.csv(cross_pval,
                      file.path(out_dir, "eigengene_cross_pval.csv"))
            cat("Eigengene cross-correlation saved\n")
        }
    }
}

# ---- Summary Concordance Table ----

summary_rows <- list()
for (cname in names(contrasts)) {
    btru_df <- contrasts[[cname]][[1]]
    shae_df <- contrasts[[cname]][[2]]

    btru_de <- sum(btru_df$padj < 0.05, na.rm = TRUE)
    shae_de <- sum(shae_df$padj < 0.05, na.rm = TRUE)

    # Direction agreement from fc_tables
    sig <- fc_tables[[cname]]
    sig <- sig[sig$btru_padj < 0.05 | sig$shae_padj < 0.05, ]
    agree <- if (nrow(sig) > 0) {
        sum(sign(sig$btru_lfc) == sign(sig$shae_lfc))
    } else { NA }
    agree_pct <- if (!is.na(agree) && nrow(sig) > 0) {
        round(agree / nrow(sig) * 100, 1)
    } else { NA }

    rho_val <- fc_cor_summary$rho[fc_cor_summary$contrast == cname]

    summary_rows[[cname]] <- data.frame(
        contrast = cname,
        btru_de_genes = btru_de,
        shae_de_genes = shae_de,
        shared_kos = fc_cor_summary$n_kos[fc_cor_summary$contrast == cname],
        spearman_rho = rho_val,
        direction_agree_n = agree,
        direction_agree_total = nrow(sig),
        direction_agree_pct = agree_pct
    )
}

summary_tbl <- bind_rows(summary_rows)
write_tsv(summary_tbl, file.path(out_dir, "concordance_summary.tsv"))

print(summary_tbl)
