#!/bin/bash

FASTA="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

BAMS=('../data/PS00685_GM12878_UBA1_m84046_240802_231812_s1.hifi_reads.bc2093.bam')

for b in "${BAMS[@]}"; do 
minimap2 -t 40 --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/bam/map-pb.bam}
samtools index ${b/bam/map-pb.bam};
done

# replace Da with IUPAC codes
for b in "${BAMS[@]}"; do 
python /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/process_DddA_bam.py -b ${b/bam/map-pb.bam};
done


COR_BAMS=('../data/PS00685_GM12878_UBA1_m84046_240802_231812_s1.hifi_reads.bc2093.map-pb_corrected.bam')

for b in "${COR_BAMS[@]}"; do
samtools index $b;
minimap2 -t 40 --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/corrected.bam/corrected_realigned.bam}
samtools index ${b/corrected.bam/corrected_realigned.bam};
done
