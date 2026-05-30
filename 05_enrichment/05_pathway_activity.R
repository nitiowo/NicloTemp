library(tidyverse)
library(DESeq2)
library(pheatmap)

annot_dir <- "00_setup/03_GO_annotation/output/step05_go_objects"
de_dir    <- "02_de_analysis/output"
out_dir   <- "05_enrichment/output/pathway_activity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Rank-based enrichment score per sample per gene set (basaed on Barbie et al. 2009)
ssgsea_score <- function(expr_mat, gene_sets) {
  gene_sets <- Filter(function(gs) sum(gs %in% rownames(expr_mat)) >= 15, gene_sets)
  if (length(gene_sets) == 0) return(NULL)

  n <- nrow(expr_mat)
  scores <- matrix(NA, length(gene_sets), ncol(expr_mat),
                   dimnames = list(names(gene_sets), colnames(expr_mat)))

  for (samp in colnames(expr_mat)) {
    rnks <- rank(expr_mat[, samp])
    for (gs in names(gene_sets)) {
      r <- sort(rnks[intersect(gene_sets[[gs]], rownames(expr_mat))])
      k <- length(r)
      scores[gs, samp] <- sum(r / n - seq_len(k) / k)
    }
  }
  scores
}

for (species in c("btru", "shae")) {
  vsd     <- readRDS(file.path(de_dir, species, paste0(species, "_vsd.rds")))
  g2go    <- readRDS(file.path(annot_dir, paste0(species, "_gene2go.Rds")))
  go2name <- readRDS(file.path(annot_dir, paste0(species, "_go2name.Rds")))

  # Build GO gene sets (size 15-300), labelled as "GO:XXXX: term name"
  valid   <- g2go %>% dplyr::count(GO) %>% filter(n >= 15, n <= 300) %>% pull(GO)
  g2go_f  <- filter(g2go, GO %in% valid)
  sets    <- split(g2go_f$gene_id, g2go_f$GO)
  nm      <- setNames(go2name$name, go2name$GO)
  names(sets) <- ifelse(!is.na(nm[names(sets)]),
                        paste0(names(sets), ": ", nm[names(sets)]),
                        names(sets))

  scores <- ssgsea_score(assay(vsd), sets)
  if (is.null(scores)) next

  write.csv(scores, file.path(out_dir, paste0(species, "_pathway_scores.csv")))

  # Top 50 most variable pathways
  top  <- scores[order(apply(scores, 1, var), decreasing = TRUE)[1:min(50, nrow(scores))], , drop = FALSE]
  scld <- t(scale(t(top)))
  rownames(scld) <- str_trunc(rownames(scld), 60)

  cd      <- as.data.frame(colData(vsd))
  flds    <- if (species == "btru") c("temperature", "infection", "niclosamide") else c("temperature", "niclosamide")
  ann_col <- cd[, intersect(flds, colnames(cd)), drop = FALSE]

  pdf(file.path(out_dir, paste0(species, "_pathway_heatmap.pdf")),
      width = 12, height = max(8, nrow(scld) * 0.25 + 3))
  pheatmap(scld, annotation_col = ann_col, show_colnames = FALSE, fontsize_row = 6,
           clustering_method = "ward.D2",
           main = paste(toupper(species), "- Top Variable Pathway Activity Scores"),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  invisible(dev.off())

  pca     <- prcomp(t(top), scale. = TRUE)
  pca_df  <- data.frame(pca$x[, 1:2], cd[rownames(pca$x), ], check.names = FALSE)
  var_pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  p <- if (species == "btru") {
    ggplot(pca_df, aes(PC1, PC2, color = temperature, shape = infection))
  } else {
    ggplot(pca_df, aes(PC1, PC2, color = temperature))
  }
  p <- p + geom_point(size = 3) +
    labs(title = paste(toupper(species), "- PCA of Pathway Scores"),
         x = paste0("PC1 (", var_pct[1], "%)"),
         y = paste0("PC2 (", var_pct[2], "%)")) +
    theme_bw(base_size = 12)

  invisible(ggsave(file.path(out_dir, paste0(species, "_pathway_pca.pdf")),
                   p, width = 8, height = 6))
}
