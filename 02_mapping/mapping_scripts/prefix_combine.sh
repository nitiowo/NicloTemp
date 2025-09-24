#!/bin/bash

B_fna="../../Data/reference_data/Btru/ncbi_dataset/data/GCA_021962125.1/GCA_021962125.1_Btru.v1_genomic.fna"
B_gtf="../../Data/reference_data/Btru/ncbi_dataset/data/GCA_021962125.1/Btru_genomic.gtf"

S_fna="../../Data/reference_data/Shae/ncbi_dataset/data/GCF_000699445.3/GCF_000699445.3_UoM_Shae.V3_genomic.fna"
S_gtf="../../Data/reference_data/Shae/ncbi_dataset/data/GCF_000699445.3/Shae_genomic.gtf"

COMBDIR="../../Data/reference_data/Combined"
mkdir -p $COMBDIR

# add Btru prefix to genome
sed '/^>/ s/^>/>Btru_/' $B_fna > ${COMBDIR}/Btru_prefixed.fna

# add Shae prefix to genome
sed '/^>/ s/^>/>Shae_/' $S_fna> ${COMBDIR}/Shae_prefixed.fna

# prefix GTF seqnames
awk 'BEGIN{OFS="\t"}
    /^#/ { print; next }          # leave comment lines unchanged
    $1 !~ /^Btru_/ { $1 = "Btru_" $1 }  # add prefix only if not already present
    { print }
' $B_gtf > ${COMBDIR}/Btru_prefixed.gtf

awk 'BEGIN{OFS="\t"}
    /^#/ { print; next }          # leave comment lines unchanged
    $1 !~ /^Shae_/ { $1 = "Shae_" $1 }  # add prefix only if not already present
    { print }
' $S_gtf > ${COMBDIR}/Shae_prefixed.gtf


# concatenate files
cat ${COMBDIR}/Btru_prefixed.fna ${COMBDIR}/Shae_prefixed.fna > ${COMBDIR}/combined.genome.fna

grep -v "#" ${COMBDIR}/Btru_prefixed.gtf > ${COMBDIR}/combined.genome.gtf
grep -v "!#" ${COMBDIR}/Shae_prefixed.gtf | grep -v "#g" >> ${COMBDIR}/combined.genome.gtf
