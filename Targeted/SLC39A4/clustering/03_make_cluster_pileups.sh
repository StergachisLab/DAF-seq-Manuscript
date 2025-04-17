#!/bin/bash
set -euo pipefail 

cd clust_pileups
FT=/mmfs1/gscratch/stergachislab/bin/ft
BW=/gscratch/stergachislab/install_dir/bedGraphToBigWig
SIZES=/mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes
REG=/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/SLC39A4_region.bed

GM_BAM="/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/pileups/GM12878_SLC39A4_PS00686_nuc.bam"
Liver_BAM="/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/pileups/Liver_SLC39A4_PS00680_nuc.bam"
MERGED="merged_GM12878_Liver_nuc.bam"

samtools merge -o $MERGED $GM_BAM $Liver_BAM
samtools index $MERGED

TEMP="temp_clust_subset.bam" # temp BAM file for subsetting reads
TEMP_IDX="${TEMP}.bai"

CL=$( ls ./*cluster_*_DA_GA_zmw.txt )

for c in $CL;
do
    samtools view -b -@ 20 -N $c $MERGED > $TEMP
    samtools index $TEMP
    NAME=${c%_DA_GA_zmw.txt};
    FT_OUT="${NAME}_ft_pileup_GA.bed";
    $FT pileup -t 30 --m6a -x "len(msp)>150" --out $FT_OUT $TEMP
    # MSP
    BG_NAME="${NAME}_ft_pileup_GA_msp150.bg"
    BW_NAME="${NAME}_ft_pileup_GA_msp150.bw"
    sed 1d $FT_OUT | bedtools intersect -a stdin -b $REG | awk '{ print $1,$2,$3,$8/$4 }' > $BG_NAME;
    $BW $BG_NAME $SIZES $BW_NAME;
    # Nucleosomes
    BG_NAME="${NAME}_ft_pileup_GA_nuc.bg"
    BW_NAME="${NAME}_ft_pileup_GA_nuc.bw"
    sed 1d $FT_OUT | bedtools intersect -a stdin -b $REG | awk '{ print $1,$2,$3,$7/$4 }' > $BG_NAME;
    $BW $BG_NAME $SIZES $BW_NAME;
    # Deamination
    BG_NAME="${NAME}_ft_pileup_GA_DA.bg"
    BW_NAME="${NAME}_ft_pileup_GA_DA.bw"
    sed 1d $FT_OUT | bedtools intersect -a stdin -b $REG | awk '{ print $1,$2,$3,$9/$4 }' > $BG_NAME;
    $BW $BG_NAME $SIZES $BW_NAME;
done

rm $TEMP
rm $TEMP_IDX


# MAP MSPs to modules ----------------------------------------------
bedmap --ec --fraction-ref 0.5 --echo --echo-map-id promoter_modules_sub_clust6.bed <( sort-bed clust6_Liver_SLC39A4_PS00680_H1_large_msp_positions.bed ) > clust6_large_msps_H1.txt
bedmap --ec --fraction-ref 0.5 --echo --echo-map-id promoter_modules_sub_clust6.bed <( sort-bed clust6_Liver_SLC39A4_PS00680_H2_large_msp_positions.bed ) > clust6_large_msps_H2.txt

