#!/bin/bash
set -euo pipefail

# Generate nf-core samplesheet from lane-merged FASTQ files
# Run after read_combine.sh. Writes samplesheet to metadata/samplesheet.csv
# Run from: NicloTemp/01_QC/

REPO="$(cd "$(dirname "$0")/.." && pwd)"

FILEDIR="${REPO}/Data/fastq_data/raw_combined"
OUTFILE="${REPO}/metadata/samplesheet.csv"

echo "sample,fastq_1,fastq_2,strandedness" > "$OUTFILE"

for file in ${FILEDIR}/*_R1_combined.fastq.gz; do
    base=$(basename "$file" | sed 's/_R1_combined.fastq.gz//')
    fwd="$file"
    rev="${FILEDIR}/${base}_R2_combined.fastq.gz"
    printf "%s,%s,%s,%s\n" "$base" "$fwd" "$rev" "reverse" >> "$OUTFILE"
done

echo "Samplesheet written to ${OUTFILE}"
