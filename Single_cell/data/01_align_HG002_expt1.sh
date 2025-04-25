#!/bin/bash

hgFASTA="hg38.fa"

BAMS=$(ls ./*.bam)

# ALIGN TO hg38
for b in $BAMS; do 
minimap2 -t 40 --MD -Y -y -a -x map-pb $hgFASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/bam/map-pb_HG38.bam}
samtools index ${b/bam/map-pb_HG38.bam};
done
