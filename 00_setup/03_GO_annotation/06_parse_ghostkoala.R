# Parse GhostKOALA KEGG KO annotation output for Btru and Shae

library(dplyr)
library(readr)
library(stringr)
library(tibble)

# ---- Paths ----
repo <- here::here()
step01_dir   <- file.path(repo, "00_setup", "03_GO_annotation", "output",
                          "step01_ncbi")
step03_dir   <- file.path(repo, "00_setup", "03_GO_annotation", "output",
                          "step03_parsed")
out_dir      <- file.path(repo, "00_setup", "03_GO_annotation", "output",
                          "step06_kegg_objects")
btru_gk_file <- file.path(repo, "Data", "Btru_GhostKoala", "user_ko_btru.txt")
shae_gk_file <- file.path(repo, "Data", "Shae_GhostKoala", "user_ko_shae.txt")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Protein - Gene ID Maps ----
btru_p2g <- read_tsv(file.path(step01_dir, "btru_protid_to_geneid.tsv"),
                     col_names = c("protein_id", "gene_id"),
                     show_col_types = FALSE)

shae_p2g <- read_tsv(file.path(step01_dir, "shae_protid_to_geneid.tsv"),
                     col_names = c("protein_id", "gene_id"),
                     show_col_types = FALSE)

# ---- Parse GhostKOALA Output ----
parse_ghostkoala <- function(filepath, p2g_map, species) {
    if (!file.exists(filepath)) {
        return(NULL)
    }

    # Read raw lines and split
    raw_lines <- read_lines(filepath)
    split_mat <- str_split_fixed(raw_lines, "\t", 2)
    raw <- tibble(
        protein_id = split_mat[, 1],
        KO         = na_if(split_mat[, 2], "")
    )

    assigned <- raw %>%
        filter(!is.na(KO), KO != "") %>%
        inner_join(p2g_map, by = "protein_id") %>%
        select(gene_id, KO) %>%
        distinct()

    n_genes <- length(unique(assigned$gene_id))
    n_ko    <- length(unique(assigned$KO))
    total   <- length(unique(p2g_map$gene_id))

    assigned
}

btru_gk <- parse_ghostkoala(btru_gk_file, btru_p2g, "Btru")
shae_gk <- parse_ghostkoala(shae_gk_file, shae_p2g, "Shae")

# ---- Compare to eggNOG KEGG ----
compare_kegg <- function(ghostkoala_df, emapper_file, species) {
    if (is.null(ghostkoala_df)) return(invisible(NULL))

    if (!file.exists(emapper_file)) {
        return(invisible(NULL))
    }

    emap <- read_tsv(emapper_file, col_names = c("gene_id", "KO"),
                     show_col_types = FALSE)

    #Gene-level overlap
    gk_genes   <- unique(ghostkoala_df$gene_id)
    emap_genes <- unique(emap$gene_id)
    both_genes <- intersect(gk_genes, emap_genes)
    only_gk    <- setdiff(gk_genes, emap_genes)
    only_emap  <- setdiff(emap_genes, gk_genes)

    # KO agreement among shared genes
    shared <- inner_join(ghostkoala_df, emap, by = "gene_id",
                         suffix = c("_gk", "_emap")) %>%
        distinct()

    agreed <- shared %>%
        filter(KO_gk == KO_emap) %>%
        distinct(gene_id, KO_gk)

    disagreed <- shared %>%
        filter(KO_gk != KO_emap) %>%
        distinct(gene_id)

    # Save disagreement table
    if (nrow(disagreed) > 0) {
        disag_detail <- inner_join(ghostkoala_df, emap, by = "gene_id",
                                   suffix = c("_gk", "_emap")) %>%
            filter(KO_gk != KO_emap) %>%
            arrange(gene_id) %>%
            rename(ghostkoala_KO = KO_gk, emapper_KO = KO_emap)

        disag_file <- file.path(out_dir,
                                paste0(tolower(species), "_kegg_conflicts.tsv"))
        write_tsv(disag_detail, disag_file)
    }
}

compare_kegg(btru_gk, file.path(step03_dir, "btru_emapper_kegg.tsv"), "Btru")
compare_kegg(shae_gk, file.path(step03_dir, "shae_emapper_kegg.tsv"), "Shae")

# ---- Build Merged KEGG TERM2GENE Tables ----
# Prefer eggNOG when available, supplement with GhostKOALA for missed genes
build_kegg <- function(ghostkoala_df, emapper_file, species) {
    if (is.null(ghostkoala_df)) return(invisible(NULL))

    emap_ko <- if (file.exists(emapper_file)) {
        read_tsv(emapper_file, col_names = c("gene_id", "KO"),
                 show_col_types = FALSE)
    } else {
        NULL
    }

    if (!is.null(emap_ko)) {
        # All eggNOG KOs + GhostKOALA KOs for genes eggNOG missed
        emap_genes <- unique(emap_ko$gene_id)
        gk_only    <- ghostkoala_df %>% filter(!gene_id %in% emap_genes)
        merged     <- bind_rows(emap_ko, gk_only) %>% distinct()
        src_label  <- "eggNOG + GhostKOALA supplement"
    } else {
        merged    <- ghostkoala_df %>% distinct()
        src_label <- "GhostKOALA only (interim)"
    }

    merged
}

btru_kegg <- build_kegg(
    btru_gk,
    file.path(step03_dir, "btru_emapper_kegg.tsv"),
    "Btru"
)

shae_kegg <- build_kegg(
    shae_gk,
    file.path(step03_dir, "shae_emapper_kegg.tsv"),
    "Shae"
)

# ---- Save ----
if (!is.null(btru_gk)) {
    write_tsv(btru_gk, file.path(out_dir, "btru_ghostkoala_kegg.tsv"),
              col_names = FALSE)
}

if (!is.null(shae_gk)) {
    write_tsv(shae_gk, file.path(out_dir, "shae_ghostkoala_kegg.tsv"),
              col_names = FALSE)
}

if (!is.null(btru_kegg)) {
    # TERM2GENE format: KO first (check clusterprofiler requirements)
    btru_t2g <- btru_kegg[, c("KO", "gene_id")]
    saveRDS(btru_t2g, file.path(out_dir, "btru_kegg_term2gene.Rds"))
    write_tsv(btru_t2g, file.path(out_dir, "btru_kegg_term2gene.tsv"),
              col_names = TRUE)
}

if (!is.null(shae_kegg)) {
    shae_t2g <- shae_kegg[, c("KO", "gene_id")]
    saveRDS(shae_t2g, file.path(out_dir, "shae_kegg_term2gene.Rds"))
    write_tsv(shae_t2g, file.path(out_dir, "shae_kegg_term2gene.tsv"),
              col_names = TRUE)
}
