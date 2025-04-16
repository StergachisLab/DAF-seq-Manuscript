#!/bin/bash
set -euo pipefail


zcat FDR-FIRE-peaks.bed.gz | awk '$NF == "true"' | awk 'BEGIN {OFS="\t"}{ print $1,$2,$3,".",$6/$7,"."}' | bgzip > FIRE-peaks_HG002_pass_coverage.bed.gz
UNMAP="/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/unreliable_blacklist_highCov_regions_merged.bed"
AUTO="/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/hg38_autosomes.sizes.bed"
bedtools intersect -v -a FIRE-peaks_HG002_pass_coverage.bed.gz -b $UNMAP | bedtools intersect -a stdin -b $AUTO -u | sort-bed - | bgzip > FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz


# FIND FIBER NAMES & MSPs THAT OVERLAP FIRE PEAKS ---------------------------------
samples=('PS00718' 'PS00756' 'PS00757' 'PS00758' 'PS00867' 'PS00868' 'PS00869' 'PS00870' 'PS00871' 'PS00872' 'PS00873' 'PS00874')

# intersect collapsed reads with FIRE peaks (read must fully overlap)
for s in "${samples[@]}"; do
OUT_BED_FIB="fibertools_msp/${s}_collapsed_haplotagged_FIRE_read_intersect.bed";
bedtools intersect -bed -F 1 -abam "fibertools_msp/${s}_consensus_BothStrands_HG38_corrected.haplotagged.m6A_nuc.bam" -b FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz >> $OUT_BED_FIB;
bgzip -@ 10 $OUT_BED_FIB;
done

# Intersect MSPs with FIRE peaks (1bp overlap OK)
LENGTHS=(150)
for LEN in "${LENGTHS[@]}"; do
    for s in "${samples[@]}"; do
    OUT_BED_MSP="fibertools_msp/${s}_collapsed_haplotagged_FIRE_MSP_${LEN}_intersect.bed";     
        bedtools intersect -wb -f 0.5 -F 0.5 -e -a "fibertools_msp/${s}_consensus_BothStrands_HG38_corrected.haplotagged.MSP_BED6_single_${LEN}.bed.gz" -b FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz | awk -F "\t" '{OFS="\t"}{ print $7,$8,$9,$4,$5,$6 }' >> $OUT_BED_MSP;
    done
    for s in "${samples[@]}"; do
        OUT_BED_MSP="fibertools_msp/${s}_collapsed_haplotagged_FIRE_MSP_${LEN}_intersect.bed";
        bgzip -@ 10 $OUT_BED_MSP;
    done
done

# filtering genocde to Ensembl Canonical TSSs
cat gencode.v45.annotation_all_tss.bed | grep Ensembl_canonical | sort-bed - > gencodev45_Ensembl_canonical_TSS.bed

# Within 1 Kb -- TSS coordinates are 20 bp centered on predicted TSS
zcat FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz | sort-bed - | bedmap --ec --echo-ref-name --echo-map --skip-unmapped --range 990 - gencodev45_Ensembl_canonical_TSS.bed > FIRE-peak_HG002_noUnrealiableCov_autosome_v45_TSS_bedmap.txt
