# 01_QC

Pre-alignment quality control and read preparation.

## Scripts

| Script | Purpose |
|--------|---------|
| `read_combine.sh` | Merge raw FASTQ files from run_52 and run_53 per sample (run first) |
| `initQC.job` | Run FastQC + MultiQC on raw reads from each lane (SGE job) |
| `sampleSheet_gen.sh` | Generate nf-core samplesheet from merged reads; writes to `metadata/samplesheet.csv` |

## Outputs

- `initial_qc/` — FastQC HTML reports and MultiQC summaries per lane
- Lane-merged reads are written to `Data/fastq_data/raw_combined/` (not tracked by git)
