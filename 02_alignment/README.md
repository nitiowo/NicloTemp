# 02_alignment

STAR-Salmon alignment via nf-core/rnaseq v3.21.0, mapping to the concatenated *Btru* + *Shae* genome.

## Scripts

| Script | Purpose |
|--------|---------|
| `prefix_combine.sh` | Add `Btru_`/`Shae_` prefixes to scaffold names and build combined FASTA + GTF |
| `nf-params.yml` | nf-core/rnaseq parameters (paths, aligner, QC modules) |
| `nfcore_rnaseq.job` | SGE job to submit the nf-core pipeline (run from this directory) |

## Outputs

- `nfcore_output/` — Full nf-core pipeline output including Salmon quant dirs, MultiQC, trimgalore, etc. (not tracked by git)
- `job_logs/` — SGE stdout/stderr logs

## Notes

- Run `prefix_combine.sh` once to build `Data/reference_data/Combined/` before submitting the job
- Requires `nfcore-env` conda environment with Nextflow installed
- Uses Singularity containers via `NXF_SINGULARITY_CACHEDIR`
- Split count matrix by `Btru_` vs `Shae_` gene ID prefix in R after pipeline completes
