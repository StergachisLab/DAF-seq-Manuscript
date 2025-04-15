#!/bin/bash
set -euo pipefail

SIZES=/mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes

bedtools genomecov -bg -i <(sortBed -i run_of_3_da_occurences_all.bed) -g $SIZES | grep chr19 > run_of_3_da_occurences_all.bg
/gscratch/stergachislab/install_dir/bedGraphToBigWig run_of_3_da_occurences_all.bg $SIZES run_of_3_da_occurences_all.bw
bedtools genomecov -dz -i <(sortBed -i run_of_3_da_occurences_all.bed) -g $SIZES | grep chr19 > run_of_3_da_occurences_all.tsv

bedtools genomecov -bg -i <(sortBed -i run_of_3_da_occurences_TOP.bed) -g $SIZES | grep chr19 > run_of_3_da_occurences_TOP.bg
/gscratch/stergachislab/install_dir/bedGraphToBigWig run_of_3_da_occurences_TOP.bg $SIZES run_of_3_da_occurences_TOP.bw
bedtools genomecov -dz -i <(sortBed -i run_of_3_da_occurences_TOP.bed) -g $SIZES | grep chr19 > run_of_3_da_occurences_TOP.tsv

bedtools genomecov -bg -i <(sortBed -i run_of_3_da_occurences_BOTTOM.bed) -g $SIZES | grep chr19 > run_of_3_da_occurences_BOTTOM.bg
/gscratch/stergachislab/install_dir/bedGraphToBigWig run_of_3_da_occurences_BOTTOM.bg $SIZES run_of_3_da_occurences_BOTTOM.bw
bedtools genomecov -dz -i <(sortBed -i run_of_3_da_occurences_BOTTOM.bed) -g $SIZES | grep chr19 > run_of_3_da_occurences_BOTTOM.tsv
