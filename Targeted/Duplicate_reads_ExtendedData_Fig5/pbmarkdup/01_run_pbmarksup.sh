#!/bin/bash
set -euo pipefail

# Liver
samtools fasta ../../data/Liver_SLC39A4_PS00680_haplotype_corrected.bam | bgzip -@ 30 > liver_hap_corrected.fasta.gz
pbmarkdup -j 60 --log-level INFO --log-file liver_hap_corrected_info.log liver_hap_corrected.fasta.gz liver_pbdup.fa --dup-file liver_dups.fa

# Colon
samtools fasta ../../data/Colon_SLC39A4_PS00682_haplotype_corrected.bam | bgzip -@ 30 > colon_hap_corrected.fasta.gz
pbmarkdup -j 60 --log-level INFO --log-file colon_hap_corrected_info.log colon_hap_corrected.fasta.gz colon_pbdup.fa --dup-file colon_dups.fa

# Heart
samtools fasta ../../data/Heart_SLC39A4_PS00681_haplotype_corrected.bam | bgzip -@ 30 > heart_hap_corrected.fasta.gz
pbmarkdup -j 60 --log-level INFO --log-file heart_hap_corrected_info.log heart_hap_corrected.fasta.gz heart_pbdup.fa --dup-file heart_dups.fa

# GM12878 SLC39A4
samtools fasta ../../data/GM12878_SLC39A4_PS00686_haplotype_corrected.bam | bgzip -@ 30 > GM12878_SLC39A4_hap_corrected.fasta.gz
pbmarkdup -j 60 --log-level INFO --log-file GM12878_SLC39A4_hap_corrected_info.log GM12878_SLC39A4_hap_corrected.fasta.gz GM12878_SLC39A4_pbdup.fa --dup-file GM12878_SLC39A4_dups.fa

# GM12878 WASF1
samtools fasta ../../data/PS00626_map-pb_corrected_realigned_WASF1_region.bam | bgzip -@ 30 > WASF1_realigned.fasta.gz
pbmarkdup -j 60 --log-level INFO --log-file WASF1_realigned_info.log WASF1_realigned.fasta.gz WASF1_pbdup.fa --dup-file WASF1_dups.fa

# GM12878 NAPA
samtools fasta ../../data/PS00626_map-pb_corrected_realigned_NAPA_region.bam | bgzip -@ 30 > napa_realigned.fasta.gz
pbmarkdup -j 60 --log-level INFO --log-file napa_realigned_info.log napa_realigned.fasta.gz napa_pbdup.fa --dup-file napa_dups.fa

# GM12878 UBA1
samtools fasta ../../data/PS00685_GM12878_UBA1_map-pb_corrected_realigned_UBA1_region.bam | bgzip -@ 30 > UBA1_realigned.fasta.gz
pbmarkdup -j 60 --log-level INFO --log-file UBA1_realigned_info.log UBA1_realigned.fasta.gz UBA1_pbdup.fa --dup-file UBA1_dups.fa

# COLO829BLT
samtools fasta ../../data/PS00719_COLO_Region2_map-pb_corrected_realigned_COLO_region.bam | bgzip -@ 30 > COLO_realigned.fasta.gz
pbmarkdup -j 60 --log-level INFO --log-file COLO_realigned_info.log COLO_realigned.fasta.gz COLO_pbdup.fa --dup-file COLO_dups.fa

