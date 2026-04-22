# NicloTemp

Dual-species RNAseq analysis of *Schistosoma haematobium* and *Bulinus truncatus* under temperature and niclosamide stress. 36 whole-snail RNA-seq samples (3 temperatures × 2 niclosamide doses × 2 infection states × 3 replicates).

## Directory structure

| Directory | Contents |
|-----------|----------|
| `00_setup/` | Pipeline setup and reference genome GO and KEGG annotation through GhostKoala and eggNOG |
| `01_alignment/` | QC and alignment using nfcore-rnaseq Nextflow pipeline |
| `02_de_analysis/` | DESeq2 models, interaction contrasts, DE visualization, and variance partitioning |
| `03_wgcna/` | WGCNA co-expression network construction for Btru and Shae, module preservation, and eigengene variance partitioning |
| `04_concordance/` | RRHO2 cross-species concordance of DE responses |
| `05_enrichment/` | ORA, GSEA, module enrichment, cross-species enrichment, pathway activity, dose-response enrichment, and gene table |
| `06_phenotype/` | Host phenotype data summarization and overview |
| `07_trait_integration/` | Module-trait correlations, gene significance, trait GSEA, and pathway-trait heatmaps |

## Workflow

- **Genome:** Concatenated *B. trunculus* (GCA_021962125.1) + *S. haematobium* (GCF_000699445.3), species-prefixed (`Btru_` / `Shae_`)
- **Aligner:** STAR-Salmon via nf-core/rnaseq v3.21.0
- **DE analysis:** DESeq2 with temperature × niclosamide × infection interaction models
- **Co-expression:** WGCNA consensus networks per species, alongside per-contrast networks for module preservation analysis
- **Enrichment:** GO ORA and GSEA via clusterProfiler; KEGG pathway activity via GSVA
- **Trait integration:** Module-trait correlation (Pearson) and trait-specific gene-set enrichment
