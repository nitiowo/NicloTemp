# 04_concordance

Cross-species transcriptomic concordance analysis between *B. truncatus* and *S. haematobium*.

Shared KEGG orthologs (KO numbers) are used to compre the two genomes. For each contrast (temperature 16 vs 24C, 32 vs 24C, niclosamide), genes are mapped to KOs, aggregated by median log2FC per KO, and RRHO2 run on matched KO lists.

## Outputs

- `output/rrho_ko_<contrast>.rds` — RRHO2 objects per contrast
- `output/rrho_ko_<contrast>.pdf` — RRHO2 heatmaps
- `output/de_landscape.tsv` — DE proportion summary per species/contrast
- `output/ko_lfc_correlation.tsv` — Spearman rho and p-values for KO-level FC agreement
- `output/ko_lfc_merged.tsv` — Per-KO fold-change table for both species
