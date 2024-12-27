#!/bin/bash
set -euo pipefail


# Find regions in each cell with >2x coverage for either CT or GA reads (misalignment)
# Find regions in each cell with >2x coverage for either haplotype (bad phasing)

samples=('PS00718' 'PS00756' 'PS00757' 'PS00758')
ST=('CT' 'GA')
HAP=('1' '2')

for s in "${samples[@]}"; do
OUT_BED="haplotagged/${s}_high_cov_filt_regions.bed"
if [[ -e $OUT_BED ]]; then rm $OUT_BED; fi
    for strand in "${ST[@]}"; do
        for h in "${HAP[@]}"; do
        samtools view -b -d "HP:${h}" "haplotagged/${s}_consensus_${strand}_HG38_corrected.haplotagged.bam" | samtools depth - | awk '{OFS="\t"}{ if($3>1) print $1,$2-1,$2 }' | sort-bed - | bedops --ec --merge - >> $OUT_BED;
        done
    done
bgzip -@ 20 $OUT_BED;
done


# MERGE ALL HIGH COV REGIONS TOGETHER
HIGH_COV_REG=$( ls haplotagged/PS*high_cov_filt_regions.bed.gz )
MERGED_HIGH_COV="haplotagged/High_Cov_regions_AllCells.bed.gz"

# Set merged high cov length cutoff at collapse overlap requirement of 11 bins (400 bp) (don't remove short overlaps, sometimes failed collapses).
zcat $HIGH_COV_REG | sort-bed - | bedops --ec --merge - | awk '{OFS="\t"}{ if(($3-$2) > 400) print $1,$2-1,$3 }' | bgzip -@ 10 > $MERGED_HIGH_COV



# IN_BAMS=$( ls haplotagged/PS00*_consensus_*_HG38_corrected.haplotagged.bam )
# for b in $IN_BAMS; do
# samtools view -@ 10 -b -L $MERGED_HIGH_COV -U ${b/bam/NoHighCov_filtered.bam} $b > ${b/bam/HighCovOnly.bam};
# samtools index ${b/bam/NoHighCov_filtered.bam};
# samtools index ${b/bam/HighCovOnly.bam};
# done


