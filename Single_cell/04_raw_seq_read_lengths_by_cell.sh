#!/bin/bash
set -euo pipefail

# get the lengths of raw sequenced reads for each cell

BAM_DIR="/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/data"


# BAMS=('PS00718_DddA_A2_HighYield_HG002_m84055_240822_201326_s2.hifi_reads.bc2068.bam' 'PS00718_DddA_A2_Seq2_m84046_241004_204600_s1.hifi_reads.bc2068.bam' 'PS00756_DddA_C2_Seq2_m84046_241004_224513_s2.hifi_reads.bc2001.bam' 'PS00757_DddA_E2_Seq2_m84046_241005_004435_s3.hifi_reads.bc2002.bam' 'PS00758_DddA_F2_Seq2_m84046_241005_024355_s4.hifi_reads.bc2003.bam' 'PS00758_DddA_F2_Seq2_m84046_241014_163432_s2.hifi_reads.bc2003.bam' 'PS00758_DddA_F2_Seq2_m84046_241021_233822_s1.hifi_reads.bc2003.bam' 'PS00758_DddA_F2_Seq2_m84046_241104_232153_s1.hifi_reads.bc2003.bam')
BAMS=('PS00758_DddA_F2_Seq2_m84046_241005_024355_s4.hifi_reads.bc2003.bam' 'PS00758_DddA_F2_Seq2_m84046_241014_163432_s2.hifi_reads.bc2003.bam' 'PS00758_DddA_F2_Seq2_m84046_241021_233822_s1.hifi_reads.bc2003.bam' 'PS00758_DddA_F2_Seq2_m84046_241104_232153_s1.hifi_reads.bc2003.bam')
for B in "${BAMS[@]}"; do
    NAME=$(echo $B | cut -d'_' -f1)
    samtools view "${BAM_DIR}/${B}" | awk -v name=$NAME '{OFS="\t"}{ print length($10), name }' >> "${NAME}_raw_read_lengths.tsv"
    echo $NAME $B
done

bgzip -@ 20 ./*_raw_read_lengths.tsv
