#!/bin/bash

FASTA="/mmfs1/gscratch/stergachislab/assemblies/simple-names/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# BAMS=('Colon_SLC39A4_PS00682_m84046_240802_231812_s1.hifi_reads.bc2086.bam' 'GM12878_SLC39A4_PS00686_m84046_240802_231812_s1.hifi_reads.bc2094.bam' 'Heart_SLC39A4_PS00681_m84046_240802_231812_s1.hifi_reads.bc2085.bam' 'HepG2_SLC39A4_PS00683_m84046_240802_231812_s1.hifi_reads.bc2088.bam' 'Liver_SLC39A4_PS00680_m84046_240802_231812_s1.hifi_reads.bc2084.bam')

BAMS=('PS00744_HepG2_SLC39A4_m84055_241001_184151_s4.hifi_reads.bc2096.bam')

for b in "${BAMS[@]}"; do 
minimap2 -t 40 --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/bam/map-pb.bam}
samtools index ${b/bam/map-pb.bam};
done

# replace Da with IUPAC codes
for b in "${BAMS[@]}"; do 
python /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/process_DddA_bam.py -b ${b/bam/map-pb.bam};
done

# COR_BAMS=('Colon_SLC39A4_PS00682_m84046_240802_231812_s1.hifi_reads.bc2086.map-pb_corrected.bam' 'GM12878_SLC39A4_PS00686_m84046_240802_231812_s1.hifi_reads.bc2094.map-pb_corrected.bam' 'Heart_SLC39A4_PS00681_m84046_240802_231812_s1.hifi_reads.bc2085.map-pb_corrected.bam' 'HepG2_SLC39A4_PS00683_m84046_240802_231812_s1.hifi_reads.bc2088.map-pb_corrected.bam' 'Liver_SLC39A4_PS00680_m84046_240802_231812_s1.hifi_reads.bc2084.map-pb_corrected.bam')

COR_BAMS=('PS00744_HepG2_SLC39A4_m84055_241001_184151_s4.hifi_reads.bc2096.map-pb_corrected.bam')

for b in "${COR_BAMS[@]}"; do
samtools index $b;
minimap2 -t 40 --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/corrected.bam/corrected_realigned.bam}
samtools index ${b/corrected.bam/corrected_realigned.bam};
done
