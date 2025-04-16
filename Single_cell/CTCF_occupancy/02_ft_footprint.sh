#!/bin/bash
set -euo pipefail


# Footprint CTCF sites within each cell
FT=/mmfs1/gscratch/stergachislab/bin/ft
BED=Merged_CTCF_motifs_with_ChIP_HG002_Peaks.bed
YAML=ctcf_modules.yaml
BAMS=$(ls /gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/fibertools_msp/*_consensus_BothStrands_HG38_corrected.haplotagged.m6A_nuc.bam)

for b in $BAMS; do
    bname=$(basename $b)
    sname="${bname%%_*}"
    $FT footprint --bed $BED --yaml $YAML $b > "ft_out/ft_CTCF_${sname}.bed"
done


# Footprint CTCF in HG002 Fiber-seq data
FS_BAM=/mmfs1/gscratch/stergachislab/mvollger/projects/FIREv2.0/results/HG002/fire/HG002.fire.cram
$FT footprint --ftx "len(msp)>150" --bed $BED --yaml $YAML $FS_BAM > "ft_out/ft_CTCF_FS_HG002.bed"

# Footprint CTCF in GM12878 Fiber-seq data
FS_BAM=/mmfs1/gscratch/stergachislab/mvollger/projects/FIREv2.0/results/GM12878/fire/GM12878.fire.cram
$FT footprint --ftx "len(msp)>150" --bed $BED --yaml $YAML $FS_BAM > "ft_out/ft_CTCF_FS_GM12878.bed"


