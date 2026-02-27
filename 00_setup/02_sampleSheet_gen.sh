#!/bin/bash
set -euo pipefail

# Generate nf-core samplesheet from raw lane-merged FASTQ files
# Run from: NicloTemp/00_setup/  (after 01_read_combine.sh)

REPO="$(cd "$(dirname "$0")/.." && pwd)"

FILEDIR="${REPO}/Data/fastq_data/raw_combined"
OUTFILE="${REPO}/00_setup/samplesheet_raw.csv"

# Header required by nf-core/rnaseq
echo "sample,fastq_1,fastq_2,strandedness" > "$OUTFILE"

for file in ${FILEDIR}/*_R1_combined.fastq.gz; do
    base=$(basename "$file" | sed 's/_R1_combined.fastq.gz//')
    fwd="$file"
    rev="${FILEDIR}/${base}_R2_combined.fastq.gz"
    printf "%s,%s,%s,%s\n" "$base" "$fwd" "$rev" "reverse" >> "$OUTFILE"
done
echo "Samplesheet written to ${OUTFILE}"
