#!/bin/bash

NUC_LEN=60
ADD_D=10
ADD_C=70
MIN=15

ft ddda-to-m6a GM12878_SLC39A4_PS00686_haplotype_corrected.bam | ft add-nucleosomes -n $NUC_LEN -d $ADD_D -c $ADD_C --min-distance-added $MIN - GM12878_SLC39A4_PS00686_nuc.bam
samtools index GM12878_SLC39A4_PS00686_nuc.bam
samtools view -b -d ST:CT GM12878_SLC39A4_PS00686_nuc.bam | samtools view -b -d HP:1 - > GM12878_SLC39A4_PS00686_nuc_CT_H1.bam
samtools view -b -d ST:CT GM12878_SLC39A4_PS00686_nuc.bam | samtools view -b -d HP:2 - > GM12878_SLC39A4_PS00686_nuc_CT_H2.bam
samtools view -b -d ST:GA GM12878_SLC39A4_PS00686_nuc.bam | samtools view -b -d HP:1 - > GM12878_SLC39A4_PS00686_nuc_GA_H1.bam
samtools view -b -d ST:GA GM12878_SLC39A4_PS00686_nuc.bam | samtools view -b -d HP:2 - > GM12878_SLC39A4_PS00686_nuc_GA_H2.bam
samtools index GM12878_SLC39A4_PS00686_nuc_CT_H1.bam
samtools index GM12878_SLC39A4_PS00686_nuc_CT_H2.bam
samtools index GM12878_SLC39A4_PS00686_nuc_GA_H1.bam
samtools index GM12878_SLC39A4_PS00686_nuc_GA_H2.bam
ft extract -t 30 GM12878_SLC39A4_PS00686_nuc_CT_H1.bam --all GM12878_SLC39A4_PS00686_ft_extract_all_CT_H1.bed
ft extract -t 30 GM12878_SLC39A4_PS00686_nuc_CT_H2.bam --all GM12878_SLC39A4_PS00686_ft_extract_all_CT_H2.bed
ft extract -t 30 GM12878_SLC39A4_PS00686_nuc_GA_H1.bam --all GM12878_SLC39A4_PS00686_ft_extract_all_GA_H1.bed
ft extract -t 30 GM12878_SLC39A4_PS00686_nuc_GA_H1.bam --all GM12878_SLC39A4_PS00686_ft_extract_all_GA_H2.bed

ft ddda-to-m6a Liver_SLC39A4_PS00680_haplotype_corrected.bam | ft add-nucleosomes -n $NUC_LEN -d $ADD_D -c $ADD_C --min-distance-added $MIN - Liver_SLC39A4_PS00680_nuc.bam
samtools index Liver_SLC39A4_PS00680_nuc.bam
samtools view -b -d ST:CT Liver_SLC39A4_PS00680_nuc.bam | samtools view -b -d HP:1 - > Liver_SLC39A4_PS00680_nuc_CT_H1.bam
samtools view -b -d ST:CT Liver_SLC39A4_PS00680_nuc.bam | samtools view -b -d HP:2 - > Liver_SLC39A4_PS00680_nuc_CT_H2.bam
samtools view -b -d ST:GA Liver_SLC39A4_PS00680_nuc.bam | samtools view -b -d HP:1 - > Liver_SLC39A4_PS00680_nuc_GA_H1.bam
samtools view -b -d ST:GA Liver_SLC39A4_PS00680_nuc.bam | samtools view -b -d HP:2 - > Liver_SLC39A4_PS00680_nuc_GA_H2.bam
samtools index Liver_SLC39A4_PS00680_nuc_CT_H1.bam
samtools index Liver_SLC39A4_PS00680_nuc_CT_H2.bam
samtools index Liver_SLC39A4_PS00680_nuc_GA_H1.bam
samtools index Liver_SLC39A4_PS00680_nuc_GA_H2.bam
ft extract -t 30 Liver_SLC39A4_PS00680_nuc_CT_H1.bam --all Liver_SLC39A4_PS00680_ft_extract_all_CT_H1.bed
ft extract -t 30 Liver_SLC39A4_PS00680_nuc_CT_H2.bam --all Liver_SLC39A4_PS00680_ft_extract_all_CT_H2.bed
ft extract -t 30 Liver_SLC39A4_PS00680_nuc_GA_H1.bam --all Liver_SLC39A4_PS00680_ft_extract_all_GA_H1.bed
ft extract -t 30 Liver_SLC39A4_PS00680_nuc_GA_H2.bam --all Liver_SLC39A4_PS00680_ft_extract_all_GA_H2.bed


# ft pileup ------------------------------------------------------------
BW=bedGraphToBigWig
SIZES=hg38.analysisSet.chrom.sizes
REG=SLC39A4_region.bed


cd ft_pileup
ft pileup -t 30  --m6a -x "len(msp)>150" --out ft_pileup_Liver.bed Liver_SLC39A4_PS00680_nuc.bam
ft pileup -t 30  --m6a -x "len(msp)>150" --out ft_pileup_GM12878.bed GM12878_SLC39A4_PS00686_nuc.bam


# MSP
sed 1d ft_pileup_Liver.bed | bedtools intersect -a stdin -b $REG | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_Liver_msp150.bg
$BW ft_pileup_Liver_msp150.bg $SIZES ft_pileup_Liver_msp150.bw
sed 1d ft_pileup_GM12878.bed | bedtools intersect -a stdin -b $REG | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_GM12878_msp150.bg
$BW ft_pileup_GM12878_msp150.bg $SIZES ft_pileup_GM12878_msp150.bw


# Nucleosomes
sed 1d ft_pileup_Liver.bed | bedtools intersect -a stdin -b $REG | awk '{ print $1,$2,$3,$7/$4 }' > ft_pileup_Liver_nuc.bg
$BW ft_pileup_Liver_nuc.bg $SIZES ft_pileup_Liver_nuc.bw
sed 1d ft_pileup_GM12878.bed | bedtools intersect -a stdin -b $REG | awk '{ print $1,$2,$3,$7/$4 }' > ft_pileup_GM12878_nuc.bg
$BW ft_pileup_GM12878_nuc.bg $SIZES ft_pileup_GM12878_nuc.bw


# Linkers
$FT pileup -t 30 --m6a -x "len(msp)<75" --out ft_pileup_Liver_linker.bed ../Liver_SLC39A4_PS00680_nuc.bam
$FT pileup -t 30 --m6a -x "len(msp)<75" --out ft_pileup_GM12878_linker.bed ../GM12878_SLC39A4_PS00686_nuc.bam

sed 1d ft_pileup_Liver_linker.bed | bedtools intersect -a stdin -b $REG | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_Liver_msp75.bg
$BW ft_pileup_Liver_msp75.bg $SIZES ft_pileup_Liver_msp75.bw
sed 1d ft_pileup_GM12878_linker.bed | bedtools intersect -a stdin -b $REG | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_GM12878_msp75.bg
$BW ft_pileup_GM12878_msp75.bg $SIZES ft_pileup_GM12878_msp75.bw

# C->T
sed 1d ft_pileup_GM12878_UBA1_H1.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$9/$4 }' > ft_pileup_GM12878_UBA1_H1_DA.bg
sed 1d ft_pileup_GM12878_UBA1_H2.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$9/$4 }' > ft_pileup_GM12878_UBA1_H2_DA.bg
bedGraphToBigWig ft_pileup_GM12878_UBA1_H1_DA.bg $SIZES ft_pileup_GM12878_UBA1_H1_DA.bw
bedGraphToBigWig ft_pileup_GM12878_UBA1_H2_DA.bg $SIZES ft_pileup_GM12878_UBA1_H2_DA.bw




