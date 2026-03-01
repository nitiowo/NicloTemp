# Shared functions for DE analysis

# DESeq results summaries
summarize_res <- function(res, label) {
  n_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
  cat(sprintf("  %-40s: %d up, %d down (padj < 0.05)\n", label, n_up, n_down))
}

# Write DESeq2 results to TSV
write_results <- function(res, filepath) {
  res_df <- as.data.frame(res)
  res_df <- rownames_to_column(res_df, "gene_id")
  res_df <- arrange(res_df, padj)
  write_tsv(res_df, filepath)
}

# Volcano plot from DESeq2 results object
make_volcano <- function(res, title, fc_thresh = 1, padj_thresh = 0.05) {
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  df$sig <- "NS"
  df$sig[df$padj < padj_thresh & df$log2FoldChange > fc_thresh] <- "Up"
  df$sig[df$padj < padj_thresh & df$log2FoldChange < -fc_thresh] <- "Down"

  ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Up" = "#e34a33", "Down" = "#2b8cbe", "NS" = "grey70")) +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", alpha = 0.5) +
    labs(title = title, x = "log2 fold change", y = "-log10(p-value)") +
    theme_minimal()
}
