#!/bin/bash

FASTA="/mmfs1/gscratch/stergachislab/assemblies/simple-names/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

BAMS=('PS00626.m84046_240619_124816_s1.bc2072.ft.bam' 'PS00627.m84046_240619_124816_s1.bc2073.ft.bam' 'PS00628.m84046_240619_124816_s1.bc2074.ft.bam' 'PS00629.m84046_240619_124816_s1.bc2075.ft.bam' 'PS00630.m84046_240619_124816_s1.bc2076.ft.bam' 'PS00631.m84046_240619_124816_s1.bc2077.ft.bam')

for b in "${BAMS[@]}"; do 
minimap2 -t 40 --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/bam/map-pb.bam}
samtools index ${b/bam/map-pb.bam};
done

# replace Da with IUPAC codes
for b in "${BAMS[@]}"; do 
python /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/process_DddA_bam.py -b ${b/bam/map-pb.bam};
done


COR_BAMS=('PS00626.m84046_240619_124816_s1.bc2072.ft.map-pb_corrected.bam' 'PS00627.m84046_240619_124816_s1.bc2073.ft.map-pb_corrected.bam' 'PS00628.m84046_240619_124816_s1.bc2074.ft.map-pb_corrected.bam' 'PS00629.m84046_240619_124816_s1.bc2075.ft.map-pb_corrected.bam' 'PS00630.m84046_240619_124816_s1.bc2076.ft.map-pb_corrected.bam' 'PS00631.m84046_240619_124816_s1.bc2077.ft.map-pb_corrected.bam')

for b in "${COR_BAMS[@]}"; do
samtools index $b;
minimap2 -t 40 --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/corrected.bam/corrected_realigned.bam}
samtools index ${b/corrected.bam/corrected_realigned.bam};
done
