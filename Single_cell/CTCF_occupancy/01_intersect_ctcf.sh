#!/bin/bash
set -euo pipefail


MOTIFS=Merged_CTCF_motifs.bed

# CTCF Motif must fully intersect a ChIP peak
bedtools intersect -u -f 1 -a $MOTIFS -b ctcf_narrowPeak_merged.bed > Merged_CTCF_motifs_with_ChIP.bed

# Filter motifs by HG002 FIRE peaks with good coverage (same peaks as used in figure 7)
PEAKS="/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz"
bedtools intersect -u -f 1 -a Merged_CTCF_motifs_with_ChIP.bed -b $PEAKS > Merged_CTCF_motifs_with_ChIP_HG002_Peaks.bed


# Map CTCF sites to ENCODE CTCF ChIA-PET Loops ----------------------------------
zcat ENCODE/GM12878_CTCF_ChIA-PET_Loops_ENCFF780PGS.bedpe.gz | awk '{OFS="\t"}{ print $1,$2,$3 }' | sort-bed - > ENCODE/GM12878_CTCF_ChIA-PET_Loops_ENCFF780PGS.first_regions.bed
zcat ENCODE/GM12878_CTCF_ChIA-PET_Loops_ENCFF780PGS.bedpe.gz | awk '{OFS="\t"}{ print $4,$5,$6 }' | sort-bed - > ENCODE/GM12878_CTCF_ChIA-PET_Loops_ENCFF780PGS.second_regions.bed

bedmap --ec --echo --echo-map --skip-unmapped --fraction-map 1 ENCODE/GM12878_CTCF_ChIA-PET_Loops_ENCFF780PGS.first_regions.bed Merged_CTCF_motifs_with_ChIP_HG002_Peaks.bed | bgzip - > Loops_first_regions_map.txt.gz
bedmap --ec --echo --echo-map --skip-unmapped --fraction-map 1 ENCODE/GM12878_CTCF_ChIA-PET_Loops_ENCFF780PGS.second_regions.bed Merged_CTCF_motifs_with_ChIP_HG002_Peaks.bed | bgzip - > Loops_second_regions_map.txt.gz


# shuffle ChIA-PET regions for a NULL distribution ----------------------------------

# Limit within the same chromosome
SIZES=/mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes
bedtools shuffle -chrom -seed 1 -i ENCODE/GM12878_CTCF_ChIA-PET_Loops_ENCFF780PGS.first_regions.bed -g $SIZES > ENCODE/SHUFFLED_GM12878_CTCF_ChIA-PET_Loops_ENCFF780PGS.first_regions.bed
bedtools shuffle -chrom -seed 2 -i ENCODE/GM12878_CTCF_ChIA-PET_Loops_ENCFF780PGS.second_regions.bed -g $SIZES > ENCODE/SHUFFLED_GM12878_CTCF_ChIA-PET_Loops_ENCFF780PGS.second_regions.bed

bedmap --ec --echo --echo-map --skip-unmapped --fraction-map 1 <(sort-bed ENCODE/SHUFFLED_GM12878_CTCF_ChIA-PET_Loops_ENCFF780PGS.first_regions.bed) Merged_CTCF_motifs_with_ChIP_HG002_Peaks.bed | bgzip - > Loops_first_regions_SHUFFLED_map.txt.gz
bedmap --ec --echo --echo-map --skip-unmapped --fraction-map 1 <(sort-bed ENCODE/SHUFFLED_GM12878_CTCF_ChIA-PET_Loops_ENCFF780PGS.second_regions.bed) Merged_CTCF_motifs_with_ChIP_HG002_Peaks.bed | bgzip - > Loops_second_regions_SHUFFLED_map.txt.gz
