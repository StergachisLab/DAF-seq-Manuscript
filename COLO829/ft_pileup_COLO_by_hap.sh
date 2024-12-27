#!/bin/bash

# samtools view hap2_T_only.bam | awk '{ print $1 }' > hap2_T_only_zmw.txt
# samtools view COLO_Region2_PS00719_GA_nuc.bam | awk '{ print $1 }' | grep -v -f hap2_T_only_zmw.txt > ref_hap_zmw.txt
# samtools view -b -N ref_hap_zmw.txt COLO_Region2_PS00719_GA_nuc.bam > ref_hap_COLO_Region2_PS00719_GA_nuc.bam
# samtools index ref_hap_COLO_Region2_PS00719_GA_nuc.bam

/mmfs1/gscratch/stergachislab/bin/ft ddda-to-m6a hap2_T_only.bam | /mmfs1/gscratch/stergachislab/bin/ft add-nucleosomes -d 10 -n 60 -c 80 --min-distance-added 10 - hap2_T_only_nuc.bam
/mmfs1/gscratch/stergachislab/bin/ft ddda-to-m6a ref_hap_COLO_Region2_PS00719_GA_nuc.bam | /mmfs1/gscratch/stergachislab/bin/ft add-nucleosomes -d 10 -n 60 -c 80 --min-distance-added 10 - ref_hap_only_nuc.bam

samtools index hap2_T_only_nuc.bam
samtools index ref_hap_only_nuc.bam
ft extract hap2_T_only_nuc.bam --all hap2_T_ft_extract_all.bed
ft extract ref_hap_only_nuc.bam --all ref_hap_ft_extract_all.bed


/mmfs1/gscratch/stergachislab/bin/ft pileup --m6a -x "len(msp)>150" --out ft_pileup_hap2_T.bed hap2_T_only_nuc.bam
/mmfs1/gscratch/stergachislab/bin/ft pileup --m6a -x "len(msp)>150" --out ft_pileup_ref_hap.bed ref_hap_only_nuc.bam


# MSP
sed 1d ft_pileup_hap2_T.bed | bedtools intersect -a stdin -b colo_region2.bed | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_hap2_T_msp150.bg
sed 1d ft_pileup_ref_hap.bed | bedtools intersect -a stdin -b colo_region2.bed | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_ref_hap_msp150.bg
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_hap2_T_msp150.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_hap2_T_msp150.bw
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_ref_hap_msp150.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_ref_hap_msp150.bw

# Nucleosomes
sed 1d ft_pileup_hap2_T.bed | bedtools intersect -a stdin -b colo_region2.bed | awk '{ print $1,$2,$3,$7/$4 }' > ft_pileup_hap2_T_nuc.bg
sed 1d ft_pileup_ref_hap.bed | bedtools intersect -a stdin -b colo_region2.bed | awk '{ print $1,$2,$3,$7/$4 }' > ft_pileup_ref_hap_nuc.bg
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_hap2_T_nuc.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_hap2_T_nuc.bw
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_ref_hap_nuc.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_ref_hap_nuc.bw

# Linkers
/mmfs1/gscratch/stergachislab/bin/ft pileup --m6a -x "len(msp)<75" --out ft_pileup_hap2_T_linker.bed hap2_T_only_nuc.bam
/mmfs1/gscratch/stergachislab/bin/ft pileup --m6a -x "len(msp)<75" --out ft_pileup_ref_hap_linker.bed ref_hap_only_nuc.bam
sed 1d ft_pileup_hap2_T_linker.bed | bedtools intersect -a stdin -b colo_region2.bed | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_hap2_T_msp75.bg
sed 1d ft_pileup_ref_hap_linker.bed | bedtools intersect -a stdin -b colo_region2.bed | awk '{ print $1,$2,$3,$8/$4 }' > ft_pileup_ref_hap_msp75.bg
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_hap2_T_msp75.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_hap2_T_msp75.bw
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_ref_hap_msp75.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_ref_hap_msp75.bw

# C->T
sed 1d ft_pileup_hap2_T.bed | bedtools intersect -a stdin -b colo_region2.bed | awk '{ print $1,$2,$3,$9/$4 }' > ft_pileup_hap2_T_DA.bg
sed 1d ft_pileup_ref_hap.bed | bedtools intersect -a stdin -b colo_region2.bed | awk '{ print $1,$2,$3,$9/$4 }' > ft_pileup_ref_hap_DA.bg
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_hap2_T_DA.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_hap2_T_DA.bw
/gscratch/stergachislab/install_dir/bedGraphToBigWig ft_pileup_ref_hap_DA.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes ft_pileup_ref_hap_DA.bw

