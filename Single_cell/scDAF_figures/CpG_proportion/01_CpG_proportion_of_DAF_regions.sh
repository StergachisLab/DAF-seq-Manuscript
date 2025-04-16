#!/bin/bash
set -euo pipefail


bedops --ec --difference ../../../phasing/hg38_autosomes.sizes.bed ../../../phasing/unreliable_blacklist_highCov_STRAND_ONLY_regions_merged.bed > daf_mappable_regions.bed
bedtools getfasta -fi /gscratch/stergachislab/assemblies/simple-names/hg38.fa -bed daf_mappable_regions.bed -fo daf_mappable_regions.fa

# Stratefy by FIRE peak vs non-FIRE peak
bedops --difference daf_mappable_regions.bed <(zcat ../../../msp_analysis/FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz) > daf_mappable_regions_NO_FIRE.bed
bedtools getfasta -fi /gscratch/stergachislab/assemblies/simple-names/hg38.fa -bed daf_mappable_regions_NO_FIRE.bed -fo daf_mappable_regions_NO_FIRE.fa
bedops --intersect daf_mappable_regions.bed <(zcat ../../../msp_analysis/FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz) > daf_mappable_regions_FIRE_ONLY.bed
bedtools getfasta -fi /gscratch/stergachislab/assemblies/simple-names/hg38.fa -bed daf_mappable_regions_FIRE_ONLY.bed -fo daf_mappable_regions_FIRE_ONLY.fa


