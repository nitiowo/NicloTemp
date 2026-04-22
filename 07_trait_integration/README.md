# 07_trait_integration

Integrates WGCNA module eigengenes and gene expression with phenotype traits.

## Scripts

| Script | Description |
|--------|-------------|
| `01_module_trait_cor.R` | Pearson correlation between module eigengenes and phenotype traits, split by infection status.  BH-adjusted p-value matrices, heatmaps, top hits table,
| `02_gene_significance.R` | Spearman correlation between individual gene expression and trait values (gene significance). Also computes Pearson module membership. Generates GS vs MM scatter plots for top module-trait pairs |
| `03_trait_gsea.R` | Pre-ranked GSEA using trait-correlated gene statistics. Runs GO and KEGG GSEA for each trait, split by infection |

## Inputs

- `03_wgcna/output/` — Module eigengenes and gene module assignments
- `06_phenotype/output/trait_summaries/` — Trait matrices (RDS format, sample × trait)
- `02_de_analysis/output/` — VST expression matrices and metadata

## Outputs

- `output/module_trait/` — Correlation and BH-adjusted p-value matrices (CSV), top hits TSV, parasite burden correlation
- `output/gene_significance/` — Gene significance TSVs, module membership TSVs, GS-MM scatter PDFs
- `output/trait_gsea/` — Per-trait GSEA result TSVs, summary table
- `output/figures/` — Module-trait heatmaps, pathway-trait heatmaps
