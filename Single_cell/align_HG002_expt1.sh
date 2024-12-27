#!/bin/bash

FASTA="/gscratch/stergachislab/assemblies/hg002v1.0.1.fasta"
hgFASTA="/gscratch/stergachislab/assemblies/simple-names/hg38.fa"

# BAMS=('PS00716_UNTREATED_D1_HG002_m84055_240822_201326_s2.hifi_reads.bc2066.bam' 'PS00717_DddA_G1_LowYield_HG002_m84055_240822_201326_s2.hifi_reads.bc2067.bam' 'PS00718_DddA_A2_HighYield_HG002_m84055_240822_201326_s2.hifi_reads.bc2068.bam')

# # -I 8G for extra memory for diploid reference genome (-I NUM       split index for every ~NUM input bases [4G])
# for b in "${BAMS[@]}"; do 
# minimap2 -t 40 -I 8G --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/bam/map-pb.bam}
# samtools index ${b/bam/map-pb.bam};
# done


# # replace Da with IUPAC codes
# for b in "${BAMS[@]}"; do 
# python /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/process_DddA_bam.py -b ${b/bam/map-pb.bam};
# done


# COR_BAMS=('PS00717_DddA_G1_LowYield_HG002_m84055_240822_201326_s2.hifi_reads.bc2067.map-pb_corrected.bam' 'PS00718_DddA_A2_HighYield_HG002_m84055_240822_201326_s2.hifi_reads.bc2068.map-pb_corrected.bam')

# for b in "${COR_BAMS[@]}"; do
# # samtools index $b;
# minimap2 -t 40 -I 8G --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/corrected.bam/corrected_realigned.bam}
# samtools index ${b/corrected.bam/corrected_realigned.bam};
# done


# BAMS=('PS00718_DddA_A2_Seq2_m84046_241004_204600_s1.hifi_reads.bc2068.bam' 'PS00756_DddA_C2_Seq2_m84046_241004_224513_s2.hifi_reads.bc2001.bam' 'PS00757_DddA_E2_Seq2_m84046_241005_004435_s3.hifi_reads.bc2002.bam' 'PS00758_DddA_F2_Seq2_m84046_241005_024355_s4.hifi_reads.bc2003.bam')
# BAMS=('PS00758_DddA_F2_Seq2_m84046_241021_233822_s1.hifi_reads.bc2003.bam')


# ALIGN TO hg38
# -I 8G for extra memory for diploid reference genome (-I NUM split index for every ~NUM input bases [4G])
for b in "${BAMS[@]}"; do 
minimap2 -t 40 --MD -Y -y -a -x map-pb $hgFASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/bam/map-pb_HG38.bam}
samtools index ${b/bam/map-pb_HG38.bam};
done
