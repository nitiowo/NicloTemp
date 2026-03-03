# NicloTemp

Dual-species RNAseq analysis of *Schistosoma haematobium* and *Bulinus truncatus* under temperature and niclosamide stress. 36 whole-snail RNA-seq samples (3 temperatures × 2 niclosamide doses × 2 infection states × 3 replicates).

## Directory structure

| Directory | Contents |
|-----------|----------|
| `00_setup/` | Pipeline setup and reference genome GO and KEGG annotation through GhostKoala and eggNOG |
| `01_alignment/` | QC and alignment using nfcore-rnaseq Nextflow pipeline |

## Workflow

- **Genome:** Concatenated *B. truncatus* (GCA_021962125.1) + *S. haematobium* (GCF_000699445.3), species-prefixed (`Btru_` / `Shae_`)
- **Aligner:** STAR-Salmon via nf-core/rnaseq v3.21.0
- **DE analysis** DESeq2
