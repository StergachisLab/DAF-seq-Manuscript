#!/bin/bash

# Region-specific FASTAs

N_BAMS=('../data/PS00720_NAPA_CpG_DddA_m84055_240822_221245_s3.hifi_reads.bc2057.bam' '../data/PS00722_NAPA_DddA_m84055_240822_221245_s3.hifi_reads.bc2059.bam' '../data/PS00724_NAPA_CpG_m84055_240822_221245_s3.hifi_reads.bc2061.bam' '../data/PS00726_NAPA_NegCtrl_m84055_240822_221245_s3.hifi_reads.bc2063.bam')
N_FASTA="NAPA_hg38.fa"

for b in "${N_BAMS[@]}"; do 
minimap2 -t 40 --MD -Y -y -a -x map-pb $N_FASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/bam/map-pb_NAPA_Ref.bam}
samtools index ${b/bam/map-pb_NAPA_Ref.bam};
done


U_BAMS=('../data/PS00721_UBA1_CpG_DddA_m84055_240822_221245_s3.hifi_reads.bc2058.bam' '../data/PS00723_UBA1_DddA_m84055_240822_221245_s3.hifi_reads.bc2060.bam' '../data/PS00725_UBA1_CpG_m84055_240822_221245_s3.hifi_reads.bc2062.bam' '../data/PS00727_UBA1_NegCtrl_m84055_240822_221245_s3.hifi_reads.bc2065.bam')
U_FASTA="UBA1_hg38.fa"

for b in "${U_BAMS[@]}"; do 
minimap2 -t 40 --MD -Y -y -a -x map-pb $U_FASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/bam/map-pb_UBA1_Ref.bam}
samtools index ${b/bam/map-pb_UBA1_Ref.bam};
done

# replace Da with IUPAC codes
for b in "${N_BAMS[@]}"; do 
python ../../General/process_DddA_bam.py -b ${b/bam/map-pb_NAPA_Ref.bam}
done

for b in "${U_BAMS[@]}"; do 
python ../../General/process_DddA_bam.py -b ${b/bam/map-pb_UBA1_Ref.bam}
done


N_BAMS=('PS00720_NAPA_CpG_DddA_m84055_240822_221245_s3.hifi_reads.bc2057.map-pb_NAPA_Ref_corrected.bam' 'PS00722_NAPA_DddA_m84055_240822_221245_s3.hifi_reads.bc2059.map-pb_NAPA_Ref_corrected.bam')
U_BAMS=('PS00721_UBA1_CpG_DddA_m84055_240822_221245_s3.hifi_reads.bc2058.map-pb_UBA1_Ref_corrected.bam' 'PS00723_UBA1_DddA_m84055_240822_221245_s3.hifi_reads.bc2060.map-pb_UBA1_Ref_corrected.bam')

for b in "${N_BAMS[@]}"; do 
minimap2 -t 40 --MD -Y -y -a -x map-pb $N_FASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/corrected.bam/corrected_realigned.bam}
samtools index ${b/corrected.bam/corrected_realigned.bam};
done

for b in "${U_BAMS[@]}"; do 
minimap2 -t 40 --MD -Y -y -a -x map-pb $U_FASTA <(samtools fastq -T "*" "$b") | samtools sort > ${b/corrected.bam/corrected_realigned.bam}
samtools index ${b/corrected.bam/corrected_realigned.bam};
done

