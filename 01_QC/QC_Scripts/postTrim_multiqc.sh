#!/bin/bash --login

fastqc52="../post_trim_qc/post_trim_52_fastqc"
fastqc53="../post_trim_qc/post_trim_53_fastqc"

unp52="${fastqc52}/unpaired"
unp53="${fastqc53}/unpaired"

mkdir -p $unp52 $unp53

mv ${fastqc52}/*_unpaired* $unp52
mv ${fastqc52}/*_unpaired* $unp53

mkdir -p ../trimmed/trim_52/unpaired ../trimmed/trim_53/unpaired

mv ../trimmed/trim_52/*_unpaired* ../trimmed/trim_52/unpaired
mv ../trimmed/trim_52/*_unpaired* ../trimmed/trim_53/unpaired

eval "$(conda shell.bash hook)"
conda activate multiqc-env

multiqc ${fastqc52} -o ../post_trim_qc/postTrim_multiqc_out_52/ -n multiqc_report_52_noUnp.html
multiqc ${fastqc53} -o ../post_trim_qc/postTrim_multiqc_out_53/ -n multiqc_report_53_noUnp.htmly

