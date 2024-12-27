#!/bin/bash

FASTA="/mmfs1/gscratch/stergachislab/assemblies/simple-names/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

BAMS=('PS00661_COLO_Dig10-4_m84046_240802_231812_s1.hifi_reads.bc2078.bam' 'PS00664_COLO_Dig20-4_m84046_240802_231812_s1.hifi_reads.bc2081.bam' 'PS00679_COLO_reAMP_m84046_240802_231812_s1.hifi_reads.bc2087.bam' 'PS00662_COLO_Dig10-1_m84046_240802_231812_s1.hifi_reads.bc2079.bam' 'PS00665_COLO_Dig20-1_m84046_240802_231812_s1.hifi_reads.bc2082.bam' 'PS00663_COLO_Dig10-25_m84046_240802_231812_s1.hifi_reads.bc2080.bam' 'PS00666_COLO_Dig20-25_m84046_240802_231812_s1.hifi_reads.bc2083.bam')

for b in "${BAMS[@]}"; do 
minimap2 -t 40 --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/bam/map-pb.bam}
samtools index ${b/bam/map-pb.bam};
done

# replace Da with IUPAC codes
for b in "${BAMS[@]}"; do 
python /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/process_DddA_bam.py -b ${b/bam/map-pb.bam};
done


COR_BAMS=('PS00661_COLO_Dig10-4_m84046_240802_231812_s1.hifi_reads.bc2078.map-pb_corrected.bam' 'PS00664_COLO_Dig20-4_m84046_240802_231812_s1.hifi_reads.bc2081.map-pb_corrected.bam' 'PS00679_COLO_reAMP_m84046_240802_231812_s1.hifi_reads.bc2087.map-pb_corrected.bam' 'PS00662_COLO_Dig10-1_m84046_240802_231812_s1.hifi_reads.bc2079.map-pb_corrected.bam' 'PS00665_COLO_Dig20-1_m84046_240802_231812_s1.hifi_reads.bc2082.map-pb_corrected.bam' 'PS00663_COLO_Dig10-25_m84046_240802_231812_s1.hifi_reads.bc2080.map-pb_corrected.bam' 'PS00666_COLO_Dig20-25_m84046_240802_231812_s1.hifi_reads.bc2083.map-pb_corrected.bam')

for b in "${COR_BAMS[@]}"; do
samtools index $b;
minimap2 -t 40 --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/corrected.bam/corrected_realigned.bam}
samtools index ${b/corrected.bam/corrected_realigned.bam};
done
