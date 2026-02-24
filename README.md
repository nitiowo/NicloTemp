# NicloTemp

Dual-species RNAseq analysis of *Schistosoma haematobium* and *Bulinus truncatus* under temperature and niclosamide stress. 36 whole-snail RNA-seq samples (3 temperatures × 2 niclosamide doses × 2 infection states × 3 replicates).

## Directory structure

| Directory | Contents |
|-----------|----------|
| `metadata/` | Sample metadata and nf-core samplesheet |
| `00_reference/` | Reference genome info (genomes stored in Data/ — not tracked by git) |
| `01_QC/` | Initial FastQC/MultiQC output; lane-merging and samplesheet scripts |
| `02_alignment/` | nf-core/rnaseq STAR-Salmon alignment scripts and output |
| `testing_grounds/` | Archived old scripts and output from earlier pipeline iterations |

## Key decisions

- **Aligner:** STAR-Salmon via nf-core/rnaseq v3.21.0
- **Genome:** Concatenated *B. truncatus* (GCA_021962125.1) + *S. haematobium* (GCF_000699445.3), species-prefixed (`Btru_` / `Shae_`)
- **BBsplit:** Disabled — reads compete for best alignment across both genomes simultaneously
- **Strandedness:** Reverse (RF), confirmed by nf-core auto-detection
- **Raw reads:** Two lanes (run_52, run_53) merged per sample before alignment

## Workflow order

1. `01_QC/read_combine.sh` — merge raw lanes per sample
2. `01_QC/sampleSheet_gen.sh` — generate nf-core samplesheet
3. `02_alignment/prefix_combine.sh` — build concatenated reference genome
4. `02_alignment/nfcore_rnaseq.job` — submit alignment pipeline