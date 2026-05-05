# 03_wgcna

Weighted Gene Co-expression Network Analysis for both species.

## Scripts

| Script | Description |
|--------|-------------|
| `01_wgcna_btru.R` | Constructs a signed WGCNA network for Btru (all 36 samples). Selects soft threshold, runs blockwiseModules, correlates eigengenes with experimental traits (main effects + pairwise interactions - no 3-way interaction) |
| `02_wgcna_shae.R` | Constructs a signed WGCNA network for Shae (18 infected samples). Same pipeline as Btru |
| `03_module_preservation.R` | Tests cross-group module preservation using WGCNA::modulePreservation(). Compares seven contrasts: infection, temperature (32v24, 16v24), niclosamide, infection at 16/32C, and niclosamide in infected samples |
| `04_variance_partition_eigengenes.R` | Partitions variance in module eigengenes across treatment factors using variancePartition |
| `05_btru_subset_networks.R` | Builds seven Btru subgroup networks (infected, uninfected, niclo-treated, niclo-control, 16C, 24C, 32C) |
| `06_btru_network_rewiring.R` | Computes pairwise eigengene-network correlations across subgroup network pairs |

## Outputs

- `output/btru/` — Btru network object, module assignments, eigengenes, trait correlations, top_genes/soft_power for subnet scripts
- `output/shae/` — Shae network object, module assignments, eigengenes, trait correlations
- `output/preservation/` — Per-contrast RDS results, Z-summary PDFs, combined summary TSV
- `output/variance_partition/` — Eigengene variance partition TSVs and PDFs
- `output/btru_<group>/` — Per-subgroup network objects, module assignments, eigengenes
- `output/rewiring/` — Delta-correlation tables, rewiring summary TSV, barplot and heatmaps
