#!/bin/bash

ft ddda-to-m6a ../UBA1_region_GA_phased.bam | ft add-nucleosomes -d 10 -n 60 -c 70 --min-distance-added 15 - GM12878_UBA1_PS00685_nuc.bam
samtools index GM12878_UBA1_PS00685_nuc.bam
samtools view -b -d HP:1 GM12878_UBA1_PS00685_nuc.bam > GM12878_UBA1_PS00685_nuc_H1.bam
samtools index GM12878_UBA1_PS00685_nuc_H1.bam
samtools view -b -d HP:2 GM12878_UBA1_PS00685_nuc.bam > GM12878_UBA1_PS00685_nuc_H2.bam
samtools index GM12878_UBA1_PS00685_nuc_H2.bam
ft extract GM12878_UBA1_PS00685_nuc_H1.bam --all GM12878_UBA1_PS00685_ft_extract_all_H1.bed
ft extract GM12878_UBA1_PS00685_nuc_H2.bam --all GM12878_UBA1_PS00685_ft_extract_all_H2.bed

cd ../
ft pileup --m6a -x "len(msp)>150" --out ft_pileup_GM12878_UBA1_H1.bed pileups/GM12878_UBA1_PS00685_nuc_H1.bam
ft pileup --m6a -x "len(msp)>150" --out ft_pileup_GM12878_UBA1_H2.bed pileups/GM12878_UBA1_PS00685_nuc_H2.bam

# MSP
sed 1d ft_pileup_GM12878_UBA1_H1.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_GM12878_UBA1_H1_msp150.bg
sed 1d ft_pileup_GM12878_UBA1_H2.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_GM12878_UBA1_H2_msp150.bg
bedGraphToBigWig ft_pileup_GM12878_UBA1_H1_msp150.bg ../hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H1_msp150.bw
bedGraphToBigWig ft_pileup_GM12878_UBA1_H2_msp150.bg ../hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H2_msp150.bw

# Nucleosomes
sed 1d ft_pileup_GM12878_UBA1_H1.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$7/$4 }' > ft_pileup_GM12878_UBA1_H1_nuc.bg
sed 1d ft_pileup_GM12878_UBA1_H2.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$7/$4 }' > ft_pileup_GM12878_UBA1_H2_nuc.bg
bedGraphToBigWig ft_pileup_GM12878_UBA1_H1_nuc.bg ../hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H1_nuc.bw
bedGraphToBigWig ft_pileup_GM12878_UBA1_H2_nuc.bg ../hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H2_nuc.bw

# Linkers
ft pileup --m6a -x "len(msp)<75" --out ft_pileup_GM12878_UBA1_H1_linker.bed GM12878_UBA1_PS00685_nuc_H1.bam
ft pileup --m6a -x "len(msp)<75" --out ft_pileup_GM12878_UBA1_H2_linker.bed GM12878_UBA1_PS00685_nuc_H2.bam
sed 1d ft_pileup_GM12878_UBA1_H1_linker.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_GM12878_UBA1_H1_msp75.bg
sed 1d ft_pileup_GM12878_UBA1_H2_linker.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_GM12878_UBA1_H2_msp75.bg
bedGraphToBigWig ft_pileup_GM12878_UBA1_H1_msp75.bg ../hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H1_msp75.bw
bedGraphToBigWig ft_pileup_GM12878_UBA1_H2_msp75.bg ../hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H2_msp75.bw

# C->T
sed 1d ft_pileup_GM12878_UBA1_H1.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$9/$4 }' > ft_pileup_GM12878_UBA1_H1_DA.bg
sed 1d ft_pileup_GM12878_UBA1_H2.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$9/$4 }' > ft_pileup_GM12878_UBA1_H2_DA.bg
bedGraphToBigWig ft_pileup_GM12878_UBA1_H1_DA.bg ../hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H1_DA.bw
bedGraphToBigWig ft_pileup_GM12878_UBA1_H2_DA.bg ../hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H2_DA.bw

