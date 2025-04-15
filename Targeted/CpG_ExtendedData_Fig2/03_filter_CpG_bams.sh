#!/bin/bash
set -euo pipefail

# Filter BAMs for the reads used in CpG analysis

# NAPA
samtools view -b -N zmws_used_NAPA.txt PS00720_NAPA_CpG_DddA_m84055_240822_221245_s3.hifi_reads.bc2057.map-pb_NAPA_Ref_corrected_realigned.bam > PS00720_NAPA_CpG_DddA_filtered_reads.bam
samtools index PS00720_NAPA_CpG_DddA_filtered_reads.bam

# UBA1
samtools view -b -N zmws_used_UBA1.txt PS00721_UBA1_CpG_DddA_m84055_240822_221245_s3.hifi_reads.bc2058.map-pb_UBA1_Ref_corrected_realigned.bam > PS00721_UBA1_CpG_DddA_filtered_reads.bam
samtools index PS00721_UBA1_CpG_DddA_filtered_reads.bam

# downsample control BAMs for figure
samtools view -b -s 0.01 PS00724_NAPA_CpG_m84055_240822_221245_s3.hifi_reads.bc2061.map-pb_NAPA_Ref.bam > PS00724_NAPA_CpG_pos_cntl_01.bam
samtools index PS00724_NAPA_CpG_pos_cntl_01.bam
samtools view -b -s 0.01 PS00725_UBA1_CpG_m84055_240822_221245_s3.hifi_reads.bc2062.map-pb_UBA1_Ref.bam > PS00725_UBA1_CpG_pos_cntl_01.bam
samtools index PS00725_UBA1_CpG_pos_cntl_01.bam
samtools view -b -s 0.1 PS00726_NAPA_NegCtrl_m84055_240822_221245_s3.hifi_reads.bc2063.map-pb_NAPA_Ref.bam > PS00726_NAPA_NegCtrl_1.bam
samtools index PS00726_NAPA_NegCtrl_1.bam
samtools view -b -s 0.01 PS00727_UBA1_NegCtrl_m84055_240822_221245_s3.hifi_reads.bc2065.map-pb_UBA1_Ref.bam > PS00727_UBA1_NegCtrl_01.bam
samtools index PS00727_UBA1_NegCtrl_01.bam
