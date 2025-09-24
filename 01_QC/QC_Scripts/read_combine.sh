#!/bin/bash
set -euo pipefail

# directories containing your trimmed lane outputs
LANE52_DIR="../trimmed/trim_52"
LANE53_DIR="../trimmed/trim_53"

# output directory for combined files
OUTDIR="../trimmed/combined"
mkdir -p "$OUTDIR"

# loop over samples (assuming files are named like L52_sampleX_R1_trimmed.fastq.gz)
for f in ${LANE52_DIR}/L52*_R1_trimmed.fastq.gz; do
    # extract sample name (remove dir, remove lane+R1 part)
    sample=$(basename "$f" | sed 's/^L52_//' | sed 's/_R1_trimmed.fastq.gz//')

    echo "Merging sample: $sample"

    # forward reads
    cat "${LANE52_DIR}/L52_${sample}_R1_trimmed.fastq.gz" \
        "${LANE53_DIR}/L53_${sample}_R1_trimmed.fastq.gz" \
        > "${OUTDIR}/${sample}_R1_combined.fastq.gz"

    # reverse reads
    cat "${LANE52_DIR}/L52_${sample}_R2_trimmed.fastq.gz" \
        "${LANE53_DIR}/L53_${sample}_R2_trimmed.fastq.gz" \
        > "${OUTDIR}/${sample}_R2_combined.fastq.gz"
done
