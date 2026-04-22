# 00_setup

Reference genome preparation and GO/KEGG annotation.

## Scripts

| Script | Description |
|--------|-------------|
| `01_read_combine.sh` | Merges per-lane FASTQ files (run 52 + 53) into single per-sample files |
| `02_sampleSheet_gen.sh` | Generates the nf-core/rnaseq samplesheet CSV from merged FASTQs |
| `04_prefix_combine.sh` | Adds `Btru_`/`Shae_` prefixes to FASTA headers and GTF seqnames, then concatenates into a single dual-species reference |

## 03_GO_annotation/

Pipeline to build gene-to-GO and gene-to-KO tables for downstream enrichment analysis.

| Script | Description |
|--------|-------------|
| `01_extract_ncbi_go.sh` | Extracts GO annotations from NCBI GTF (Btru) and GAF (Shae); builds protein-to-gene ID maps for eggNOG |
| `02_run_emapper.job` | SGE job to run eggNOG-mapper on Btru and Shae protein FASTAs |
| `03_parse_emapper.sh` | Parses eggNOG-mapper output into gene-to-GO and gene-to-KO TSV files |
| `04_merge_go.sh` | Merges NCBI and eggNOG GO annotations into a unified gene2go table per species |
| `05_build_gene2go.R` | Builds clusterProfiler gene2go and GO-to-name Rds objects |
| `06_parse_ghostkoala.R` | Parses GhostKOALA KO output and merges with eggNOG KO coverage for KEGG enrichment |

## Outputs

- `Data/fastq_data/raw_combined/` — merged FASTQs
- `00_setup/samplesheet_raw.csv` — nf-core samplesheet
- `Data/reference_data/Combined/` — prefixed dual-species FASTA and GTF
- `00_setup/03_GO_annotation/output/step05_go_objects/` — Rds gene2go objects
- `00_setup/03_GO_annotation/output/step06_kegg_objects/` — Rds KEGG term2gene objects
