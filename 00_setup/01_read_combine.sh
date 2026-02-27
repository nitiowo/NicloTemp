#!/bin/bash
set -euo pipefail

# Merge raw (untrimmed) FASTQ files from two sequencing lanes per sample
# nf-core handles trimming internally so this just combines lanes
# Run from: NicloTemp/00_setup/

REPO="$(cd "$(dirname "$0")/.." && pwd)"

LANE52_DIR="${REPO}/Data/fastq_data/run_52"
LANE53_DIR="${REPO}/Data/fastq_data/run_53"

OUTDIR="${REPO}/Data/fastq_data/raw_combined"
mkdir -p "$OUTDIR"

# Files are named like: L52_<SampleID>_R1_001.fastq.gz
# Strip the lane prefix (L52_) and suffix (_R1_001.fastq.gz) to get sample ID
for f in ${LANE52_DIR}/L52*_R1_001.fastq.gz; do
    sample=$(basename "$f" | sed 's/^L52_//' | sed 's/_R1_001.fastq.gz//')
    echo "Merging sample: $sample"

    # Concatenate R1 from both lanes
    cat "${LANE52_DIR}/L52_${sample}_R1_001.fastq.gz" \
        "${LANE53_DIR}/L53_${sample}_R1_001.fastq.gz" \
        > "${OUTDIR}/${sample}_R1_combined.fastq.gz"

    # Concatenate R2 from both lanes
    cat "${LANE52_DIR}/L52_${sample}_R2_001.fastq.gz" \
        "${LANE53_DIR}/L53_${sample}_R2_001.fastq.gz" \
        > "${OUTDIR}/${sample}_R2_combined.fastq.gz"
done
echo "Done — merged reads in ${OUTDIR}"
