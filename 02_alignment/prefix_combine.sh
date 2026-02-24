#!/bin/bash
set -euo pipefail

# Add species prefixes to Btru and Shae reference sequences and annotations,
# then concatenate into a single genome for dual-species STAR-Salmon alignment
# Run from: NicloTemp/02_alignment/
# Output: Data/reference_data/Combined/combined.genome.{fna,gtf}

REPO="$(cd "$(dirname "$0")/.." && pwd)"

B_fna="${REPO}/Data/reference_data/Btru/GCA_021962125.1/GCA_021962125.1_Btru.v1_genomic.fna"
B_gtf="${REPO}/Data/reference_data/Btru/GCA_021962125.1/Btru_genomic.gtf"

S_fna="${REPO}/Data/reference_data/Shae/GCF_000699445.3/GCF_000699445.3_UoM_Shae.V3_genomic.fna"
S_gtf="${REPO}/Data/reference_data/Shae/GCF_000699445.3/Shae_genomic.gtf"

COMBDIR="${REPO}/Data/reference_data/Combined"
mkdir -p "$COMBDIR"

# Add Btru_ prefix to genome sequences and GTF seqnames
sed '/^>/ s/^>/>Btru_/' "$B_fna" > "${COMBDIR}/Btru_prefixed.fna"

sed '/^>/ s/^>/>Shae_/' "$S_fna" > "${COMBDIR}/Shae_prefixed.fna"

# Use -F"\t" so awk does not split the attribute column on spaces,
# which would break the required 9-column GTF format
awk -F"\t" 'BEGIN{OFS="\t"}
    /^#/ { print; next }
    NF == 0 { next }
    $1 !~ /^Btru_/ { $1 = "Btru_" $1 }
    { print }
' "$B_gtf" > "${COMBDIR}/Btru_prefixed.gtf"

awk -F"\t" 'BEGIN{OFS="\t"}
    /^#/ { print; next }
    NF == 0 { next }
    $1 !~ /^Shae_/ { $1 = "Shae_" $1 }
    { print }
' "$S_gtf" > "${COMBDIR}/Shae_prefixed.gtf"

# Concatenate FASTA and GTF
cat "${COMBDIR}/Btru_prefixed.fna" "${COMBDIR}/Shae_prefixed.fna" \
    > "${COMBDIR}/combined.genome.fna"

awk -F"\t" 'NF == 9 && !/^#/' "${COMBDIR}/Btru_prefixed.gtf"  > "${COMBDIR}/combined.genome.gtf"
awk -F"\t" 'NF == 9 && !/^#/' "${COMBDIR}/Shae_prefixed.gtf" >> "${COMBDIR}/combined.genome.gtf"

echo "Combined FASTA: $(grep -c '^>' ${COMBDIR}/combined.genome.fna) scaffolds"
echo "Combined GTF:   $(awk '$3=="gene"' ${COMBDIR}/combined.genome.gtf | wc -l) gene records"
echo "Output: ${COMBDIR}/"
