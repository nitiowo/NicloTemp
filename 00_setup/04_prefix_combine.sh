#!/bin/bash

# Add species prefixes to Btru and Shae reference sequences and annotations,
# then concatenate into a single genome for dual-species STAR-Salmon alignment
# Run from: NicloTemp/00_setup/

REPO="$(cd "$(dirname "$0")/.." && pwd)"

B_fna="${REPO}/Data/reference_data/Btru/GCA_021962125.1/GCA_021962125.1_Btru.v1_genomic.fna"
B_gtf="${REPO}/Data/reference_data/Btru/GCA_021962125.1/Btru_genomic.gtf"

S_fna="${REPO}/Data/reference_data/Shae/GCF_000699445.3/GCF_000699445.3_UoM_Shae.V3_genomic.fna"
S_gtf="${REPO}/Data/reference_data/Shae/GCF_000699445.3/Shae_genomic.gtf"

COMBDIR="${REPO}/Data/reference_data/Combined"
mkdir -p $COMBDIR

# ---- Prefix FASTA headers ----

# Prepend "Btru_" to each scaffold name in the Btru genome
sed '/^>/ s/^>/>Btru_/' $B_fna > ${COMBDIR}/Btru_prefixed.fna

# Prepend "Shae_" to each scaffold name in the Shae genome
sed '/^>/ s/^>/>Shae_/' $S_fna > ${COMBDIR}/Shae_prefixed.fna

# ---- Prefix GTF seqnames ----

# Must use -F"\t" so awk treats tabs as the only field separator
# Without this, awk splits the 9th column (attributes) on spaces and
# re-joins with tabs, breaking the 9-column GTF format
awk -F"\t" 'BEGIN{OFS="\t"}
    /^#/ { print; next }
    NF == 0 { next }
    $1 !~ /^Btru_/ { $1 = "Btru_" $1 }
    { print }
' $B_gtf > ${COMBDIR}/Btru_prefixed.gtf

awk -F"\t" 'BEGIN{OFS="\t"}
    /^#/ { print; next }
    NF == 0 { next }
    $1 !~ /^Shae_/ { $1 = "Shae_" $1 }
    { print }
' $S_gtf > ${COMBDIR}/Shae_prefixed.gtf

# ---- Concatenate references ----

cat ${COMBDIR}/Btru_prefixed.fna ${COMBDIR}/Shae_prefixed.fna > ${COMBDIR}/combined.genome.fna

# Strip comment and blank lines to avoid header conflicts
awk -F"\t" 'NF == 9 && !/^#/' ${COMBDIR}/Btru_prefixed.gtf >  ${COMBDIR}/combined.genome.gtf
awk -F"\t" 'NF == 9 && !/^#/' ${COMBDIR}/Shae_prefixed.gtf >> ${COMBDIR}/combined.genome.gtf

# ---- Verify ----

echo "Combined FASTA: $(grep -c '^>' ${COMBDIR}/combined.genome.fna) scaffolds"
echo "Combined GTF:   $(awk '$3=="gene"' ${COMBDIR}/combined.genome.gtf | wc -l) gene records"
echo "Output: ${COMBDIR}/"
