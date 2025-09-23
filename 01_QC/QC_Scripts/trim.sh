#!/bin/bash --login

#Usage: bash trim.sh <sample_dir> <output_dir> <#_of_cores>

SAMPDIR=$1
OUTDIR=$2
ncores=$3

mkdir -p ${OUTDIR}

for r1 in "${SAMPDIR}"/*_R1_001.fastq.gz # for each sample

do
    echo "CURRENT INPUT FILE IS $f"
    case "$r1" in
        *Undetermined*) continue ;;  # skip
    esac
    
    f=${r1##*/}; base=${f%.fastq.gz}
    echo "BASENAME IS $base"
    
    IFS=_ read -r -a p <<< "$base"; new="${p[0]}_${p[1]}_${p[2]}"
    echo "Processing $new"
    
    n=${r1%%_R1_001.fastq.gz}
    trimmomatic PE -threads ${ncores} ${n}_R1_001.fastq.gz  ${n}_R2_001.fastq.gz \
    ${OUTDIR}/${new}_R1_trimmed.fastq.gz ${OUTDIR}/${new}_R1_unpaired.fastq.gz \
    ${OUTDIR}/${new}_R2_trimmed.fastq.gz ${OUTDIR}/${new}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done
