#!/bin/bash
set -euo pipefail

FASTA="hg38.fa"

samples=('PS00718' 'PS00756' 'PS00757' 'PS00758' 'PS00867' 'PS00868' 'PS00869' 'PS00870' 'PS00871' 'PS00872' 'PS00873' 'PS00874')

for s in "${samples[@]}";
do
    LSC="./consensus_seqs/${s}_clipped_*.fa";
    OUT_BAM="consensus_bams/${s}_consensus_BothStrands_HG38.bam";
    minimap2 -t 30 --MD -Y -y -a -x map-pb $FASTA <(cat $LSC) | samtools sort > $OUT_BAM;
    samtools index $OUT_BAM;
    python ../../General/process_DddA_bam.py -b $OUT_BAM;
    samtools index ${OUT_BAM/.bam/_corrected.bam};
done

