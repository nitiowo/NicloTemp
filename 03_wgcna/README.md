# 03_wgcna

Weighted Gene Co-expression Network Analysis for both species.

## Scripts

| Script | Description |
|--------|-------------|
| `01_wgcna_btru.R` | Constructs a signed WGCNA network for Btru (all 36 samples). Selects soft threshold, runs blockwiseModules, correlates eigengenes with experimental traits |
| `02_wgcna_shae.R` | Constructs a signed WGCNA network for Shae (18 infected samples). Same pipeline as Btru |
| `03_module_preservation.R` | Creates sub-networks per-contrast for Btru, tests cross-group module preservation using WGCNA::modulePreservation(). Compares infection, temperature, and niclosamide subsets |
| `04_variance_partition_eigengenes.R` | Partitions variance in module eigengenes across treatment factors using variancePartition |

## Outputs

- `output/btru/` — Btru network object, module assignments, eigengenes, trait correlations
- `output/shae/` — Shae network object, module assignments, eigengenes, trait correlations
- `output/preservation/` — Per-contrast RDS results, Z-summary PDFs, combined summary TSV
- `output/variance_partition/` — Eigengene variance partition TSVs and PDFs
