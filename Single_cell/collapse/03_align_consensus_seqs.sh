#!/bin/bash
set -euo pipefail

FASTA="/mmfs1/gscratch/stergachislab/assemblies/simple-names/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

samples=('PS00718' 'PS00756' 'PS00757' 'PS00758' 'PS00867' 'PS00868' 'PS00869' 'PS00870' 'PS00871' 'PS00872' 'PS00873' 'PS00874')

for s in "${samples[@]}";
do
    LSC="./consensus_seqs/${s}_clipped_*.fa";
    OUT_BAM="consensus_bams/${s}_consensus_BothStrands_HG38.bam";
    minimap2 -t 30 --MD -Y -y -a -x map-pb $FASTA <(cat $LSC) | samtools sort > $OUT_BAM;
    samtools index $OUT_BAM;
    python /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/process_DddA_bam.py -b $OUT_BAM;
    samtools index ${OUT_BAM/.bam/_corrected.bam};
done

