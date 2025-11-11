#!/bin/bash/

# Script to generate samplesheet for nf-core rnaseq pipeline from directory of input files

# usage: sampleSheet_gen.sh

filedir="/temp180/mpfrende/nvincen2/NicloTemp/01_QC/trimmed/combined_lanes"
outdir="/temp180/mpfrende/nvincen2/NicloTemp/01_QC/trimmed/combined_lanes"

outfile=${outdir}/samplesheet.csv

# Build header line and save to outfile
line1="sample,fastq_1,fastq_2,strandedness"
echo $line1 > $outfile

# For each forward read file
for file in ${filedir}/*_R1_combined.fastq.gz;
do
	# Grab sample name from filename
	base=$(basename $file | sed 's/_R1_combined.fastq.gz//')

	# Forward and reverse filenames with full path
	fwd=$file
	rev=${filedir}/${base}_R2_combined.fastq.gz

	# Build a line for each file pair and save to outfile
	line="$base,$fwd,$rev,auto"
	echo $line
	printf "%s\n" $line >> $outfile
done
