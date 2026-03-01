# Build gene2go data frames for clusterProfiler (enricher(), GSEA)

library(readr)
library(GO.db)
library(AnnotationDbi)

# ---- Paths ----
repo <- here::here()
step04_dir <- file.path(repo, "00_setup", "03_GO_annotation", "output",
                        "step04_merged_go")
step05_dir <- file.path(repo, "00_setup", "03_GO_annotation", "output",
                        "step05_go_objects")
dir.create(step05_dir, showWarnings = FALSE, recursive = TRUE)

btru_file <- file.path(step04_dir, "btru_gene2go.tsv")
shae_file <- file.path(step04_dir, "shae_gene2go.tsv")

# ---- Load Gene2GO Tables ----
btru_g2go <- read_tsv(btru_file, col_names = c("gene_id", "GO"),
                      show_col_types = FALSE)

shae_g2go <- read_tsv(shae_file, col_names = c("gene_id", "GO"),
                      show_col_types = FALSE)


# ---- Build GO Term name tables ----
# Gather all unique GO terms
all_go <- unique(c(btru_g2go$GO, shae_g2go$GO))

go_info <- AnnotationDbi::select(GO.db,
                                 keys = all_go,
                                 columns = c("GOID", "TERM", "ONTOLOGY"),
                                 keytype = "GOID")

# Drop NAs, split and filter\
go_info <- go_info[!is.na(go_info$TERM), ]
colnames(go_info) <- c("GO", "name", "ontology")

btru_go2name <- go_info[go_info$GO %in% btru_g2go$GO, ]
shae_go2name <- go_info[go_info$GO %in% shae_g2go$GO, ]

btru_g2go <- btru_g2go[btru_g2go$GO %in% go_info$GO, ]
shae_g2go <- shae_g2go[shae_g2go$GO %in% go_info$GO, ]

# ---- Save ----
saveRDS(btru_g2go, file.path(step05_dir, "btru_gene2go.Rds"))
saveRDS(shae_g2go, file.path(step05_dir, "shae_gene2go.Rds"))
saveRDS(btru_go2name, file.path(step05_dir, "btru_go2name.Rds"))
saveRDS(shae_go2name, file.path(step05_dir, "shae_go2name.Rds"))
