#!/bin/bash
set -euo pipefail

FASTA="/gscratch/stergachislab/assemblies/hg002v1.0.1.fasta"

# MERGE PS00758 BAMS for HG002 ref alignment -------------------
MERGED_BAM="PS00758_DddA_F2_Seq2_unaligned_merged.bam"
samtools merge -@ 20 -o $MERGED_BAM PS00758_DddA_F2_Seq2_m84046_241005_024355_s4.hifi_reads.bc2003.bam PS00758_DddA_F2_Seq2_m84046_241014_163432_s2.hifi_reads.bc2003.bam PS00758_DddA_F2_Seq2_m84046_241021_233822_s1.hifi_reads.bc2003.bam

# -I 8G for extra memory for diploid reference genome (-I NUM split index for every ~NUM input bases [4G])
minimap2 -t 40 -I 8G --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" "$MERGED_BAM") | samtools sort > ${MERGED_BAM/bam/map-pb.bam}
samtools index ${MERGED_BAM/bam/map-pb.bam}

mv ${MERGED_BAM/bam/map-pb.bam} PS00758_DddA_F2_Seq2_HG002ref_merged.map-pb.bam
/mmfs1/gscratch/stergachislab/bin/ft clear-kinetics -t 40 PS00758_DddA_F2_Seq2_HG002ref_merged.map-pb.bam PS00758_DddA_F2_Seq2_HG002ref_merged.map-pb.no_kinetics.bam
samtools index PS00758_DddA_F2_Seq2_HG002ref_merged.map-pb.no_kinetics.bam
