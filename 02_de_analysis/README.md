# 02_de_analysis

DESeq2 differential expression analysis for both species.

## Scripts

| Script | Description |
|--------|-------------|
| `01_qc_and_species_split.R` | Loads nf-core salmon counts, splits by species prefix, generates QC plots and species proportion summaries |
| `02_deseq2_btru.R` | DESeq2 model for Btru (all 36 samples). Main effects extracted using averaged contrasts over the model matrix. Shrinkage via ashr |
| `03_deseq2_shae.R` | DESeq2 model for Shae (18 infected samples). Two-factor design (temperature × niclosamide). Shrinkage via apeglm |
| `04_interaction_contrasts.R` | Extracts and shrinks 2-way interaction coefficients for both species |
| `05_de_visualization.R` | MA plots, top-gene heatmaps, dispersion plots, and UpSet contrast overlap plots |
| `06_variance_partition.R` | Variance partitioning of VST expression across treatment factors using variancePartition |
| `de_functions.R` | Shared  functions: `summarize_res()`, `write_results()`, `make_volcano()` |

## Outputs

- `output/counts/` — species-split raw count matrices and metadata Rds
- `output/btru/` — Btru DESeq2 objects, result TSVs, shrunk TSVs, PCA and volcano PDFs
- `output/shae/` — Shae DESeq2 objects, result TSVs, shrunk TSVs, PCA PDFs
- `output/interactions/` — Interaction contrast TSVs (shrunk and unshrunk)
- `output/visualization/` — MA plots, heatmaps, dispersion plots, UpSet plots
- `output/visualization/variance_partition/` — Variance partition violin PDFs and summary TSVs
