#!/bin/bash

/mmfs1/gscratch/stergachislab/bin/ft ddda-to-m6a /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/UBA1/UBA1_region_GA_phased.bam | ft add-nucleosomes -d 10 -n 60 -c 70 --min-distance-added 15 - GM12878_UBA1_PS00685_nuc.bam
samtools index GM12878_UBA1_PS00685_nuc.bam
samtools view -b -d HP:1 GM12878_UBA1_PS00685_nuc.bam > GM12878_UBA1_PS00685_nuc_H1.bam
samtools index GM12878_UBA1_PS00685_nuc_H1.bam
samtools view -b -d HP:2 GM12878_UBA1_PS00685_nuc.bam > GM12878_UBA1_PS00685_nuc_H2.bam
samtools index GM12878_UBA1_PS00685_nuc_H2.bam
ft extract GM12878_UBA1_PS00685_nuc_H1.bam --all GM12878_UBA1_PS00685_ft_extract_all_H1.bed
ft extract GM12878_UBA1_PS00685_nuc_H2.bam --all GM12878_UBA1_PS00685_ft_extract_all_H2.bed

cd ../
/mmfs1/gscratch/stergachislab/bin/ft pileup --m6a -x "len(msp)>150" --out ft_pileup_GM12878_UBA1_H1.bed pileups/GM12878_UBA1_PS00685_nuc_H1.bam
/mmfs1/gscratch/stergachislab/bin/ft pileup --m6a -x "len(msp)>150" --out ft_pileup_GM12878_UBA1_H2.bed pileups/GM12878_UBA1_PS00685_nuc_H2.bam

# MSP
sed 1d ft_pileup_GM12878_UBA1_H1.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_GM12878_UBA1_H1_msp150.bg
sed 1d ft_pileup_GM12878_UBA1_H2.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_GM12878_UBA1_H2_msp150.bg
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_GM12878_UBA1_H1_msp150.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H1_msp150.bw
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_GM12878_UBA1_H2_msp150.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H2_msp150.bw

# Nucleosomes
sed 1d ft_pileup_GM12878_UBA1_H1.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$7/$4 }' > ft_pileup_GM12878_UBA1_H1_nuc.bg
sed 1d ft_pileup_GM12878_UBA1_H2.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$7/$4 }' > ft_pileup_GM12878_UBA1_H2_nuc.bg
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_GM12878_UBA1_H1_nuc.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H1_nuc.bw
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_GM12878_UBA1_H2_nuc.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H2_nuc.bw

# Linkers
/mmfs1/gscratch/stergachislab/bin/ft pileup --m6a -x "len(msp)<75" --out ft_pileup_GM12878_UBA1_H1_linker.bed pileups/GM12878_UBA1_PS00685_nuc_H1.bam
/mmfs1/gscratch/stergachislab/bin/ft pileup --m6a -x "len(msp)<75" --out ft_pileup_GM12878_UBA1_H2_linker.bed pileups/GM12878_UBA1_PS00685_nuc_H2.bam
sed 1d ft_pileup_GM12878_UBA1_H1_linker.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_GM12878_UBA1_H1_msp75.bg
sed 1d ft_pileup_GM12878_UBA1_H2_linker.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_GM12878_UBA1_H2_msp75.bg
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_GM12878_UBA1_H1_msp75.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H1_msp75.bw
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_GM12878_UBA1_H2_msp75.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H2_msp75.bw

# C->T
sed 1d ft_pileup_GM12878_UBA1_H1.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$9/$4 }' > ft_pileup_GM12878_UBA1_H1_DA.bg
sed 1d ft_pileup_GM12878_UBA1_H2.bed | bedtools intersect -a stdin -b region_UBA1.bed | awk '{ print $1,$2,$3,$9/$4 }' > ft_pileup_GM12878_UBA1_H2_DA.bg
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_GM12878_UBA1_H1_DA.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H1_DA.bw
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_GM12878_UBA1_H2_DA.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_GM12878_UBA1_H2_DA.bw

