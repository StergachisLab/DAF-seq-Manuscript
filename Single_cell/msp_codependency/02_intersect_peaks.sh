#!/bin/bash
set -euo pipefail


# Map FIRE peaks to p/q arms
FIRE_PEAKS=../msp_analysis/FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz

bedtools intersect -a $FIRE_PEAKS -b chrom_arm_coords_filtered.bed -wo | bgzip - > FIRE-peaks_HG002_noUnreliableCoverage.chrom_arm_map.bed.gz


# A/B compartments
bigWigToBedGraph Hi-C/GM12878_in_situ_Hi-C_compartments_4DNES3JX38V5_preferred_protocol_4DNFILYQ1PAY.bw Hi-C/GM12878_in_situ_Hi-C_compartments_4DNES3JX38V5_preferred_protocol_4DNFILYQ1PAY.bg

cat Hi-C/GM12878_in_situ_Hi-C_compartments_4DNES3JX38V5_preferred_protocol_4DNFILYQ1PAY.bg | awk '$4 > 0 && $4 != "nan"' | sort -k1,1 -k2,2n | bedtools merge -i stdin > Hi-C/GM12878_in_situ_Hi-C_A_compartments.bed
cat Hi-C/GM12878_in_situ_Hi-C_compartments_4DNES3JX38V5_preferred_protocol_4DNFILYQ1PAY.bg | awk '$4 < 0 && $4 != "nan"' | sort -k1,1 -k2,2n | bedtools merge -i stdin > Hi-C/GM12878_in_situ_Hi-C_B_compartments.bed

# FIRE peaks must be fully within a compartment region
bedtools intersect -f 1 -a $FIRE_PEAKS -b ../Hi-C/GM12878_in_situ_Hi-C_A_compartments.bed -wo | bgzip - > FIRE-peaks_HG002_noUnreliableCoverage.A_compartments.bed.gz
bedtools intersect -f 1 -a $FIRE_PEAKS -b ../Hi-C/GM12878_in_situ_Hi-C_B_compartments.bed -wo | bgzip - > FIRE-peaks_HG002_noUnreliableCoverage.B_compartments.bed.gz

# extract +/- 100 kb from compartment boundaries (Only use start coordinates of each comaprtment as the boundary coordinate to avoid double counting)
cat Hi-C/GM12878_in_situ_Hi-C_A_compartments.bed | awk '{OFS="\t"}{ print $1,($2-100000),($2+100000) }' | awk '$2 > 0' | sort-bed - > compartment_A_boundaries_100kb.bed
cat Hi-C/GM12878_in_situ_Hi-C_B_compartments.bed | awk '{OFS="\t"}{ print $1,($2-100000),($2+100000) }' | awk '$2 > 0' | sort-bed - > compartment_B_boundaries_100kb.bed

# FIRE peaks must be fully within a region
bedtools intersect -f 1 -a $FIRE_PEAKS -b compartment_A_boundaries_100kb.bed -wo | bgzip - > compartment_A_boundaries_100kb_FIRE-peaks.bed.gz
bedtools intersect -f 1 -a $FIRE_PEAKS -b compartment_B_boundaries_100kb.bed -wo | bgzip - > compartment_B_boundaries_100kb_FIRE-peaks.bed.gz
