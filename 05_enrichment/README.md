# 05_enrichment

Functional enrichment analysis for DE results and WGCNA modules.

## Scripts

| Script | Description |
|--------|-------------|
| `01_ora.R` | GO and KEGG over-representation analysis for all DE contrasts (up- and down-regulated gene sets separately) |
| `02_gsea.R` | GO and KEGG pre-ranked GSEA using Wald statistics from DESeq2 |
| `03_module_enrichment.R` | GO and KEGG ORA per WGCNA module for both species |
| `04_cross_species.R` | Compares enriched GO/KEGG terms between Btru and Shae for shared contrasts; GSEA NES scatter plots |
| `gene_table.R` | Master gene annotation table combining DE results, WGCNA module assignments, GO terms, and KEGG orthologs for both species |

## Outputs

- `output/ora/` — Per-contrast ORA results (TSV and dot plots)
- `output/gsea/` — Per-contrast GSEA results, dot plots, ridge plots, highlights table
- `output/module_enrichment/` — Per-module ORA results, tile heatmaps, summary tables
- `output/cross_species/` — Shared GO terms across species, NES scatter plots
