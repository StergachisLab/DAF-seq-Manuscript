#!/bin/bash
set -euo pipefail


# Get Msss.I. CpG methylation stats

# https://github.com/PacificBiosciences/pb-CpG-tools


aligned_bam_to_cpg_scores --pileup-mode count --output-prefix CpG_tools_out/PS00724_NAPA_CpG --bam PS00724_NAPA_CpG_m84055_240822_221245_s3.hifi_reads.bc2061.map-pb_NAPA_Ref.bam
aligned_bam_to_cpg_scores --pileup-mode count --output-prefix CpG_tools_out/PS00725_UBA1_CpG --bam PS00725_UBA1_CpG_m84055_240822_221245_s3.hifi_reads.bc2062.map-pb_UBA1_Ref.bam
aligned_bam_to_cpg_scores --pileup-mode count --output-prefix CpG_tools_out/PS00726_NAPA_NegCtrl --bam PS00726_NAPA_NegCtrl_m84055_240822_221245_s3.hifi_reads.bc2063.map-pb_NAPA_Ref.bam
aligned_bam_to_cpg_scores --pileup-mode count --output-prefix CpG_tools_out/PS00727_UBA1_NegCtrl --bam PS00727_UBA1_NegCtrl_m84055_240822_221245_s3.hifi_reads.bc2065.map-pb_UBA1_Ref.bam
