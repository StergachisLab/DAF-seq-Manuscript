#!/bin/bash

bedtools intersect -a nuc_positions.bed -b ../footprinting/napa_region.bed | bedtools genomecov -bg -i stdin -g ../footprinting/hg38.analysisSet.chrom.sizes > nuc_positions.bg
bedtools intersect -a small_msp_positions.bed -b ../footprinting/napa_region.bed | bedtools genomecov -bg -i stdin -g ../footprinting/hg38.analysisSet.chrom.sizes > small_msp_positions.bg
bedtools intersect -a large_msp_positions.bed -b ../footprinting/napa_region.bed | bedtools genomecov -bg -i stdin -g ../footprinting/hg38.analysisSet.chrom.sizes > large_msp_positions.bg
bedGraphToBigWig nuc_positions.bg ../footprinting/hg38.analysisSet.chrom.sizes nuc_positions.bw
bedGraphToBigWig small_msp_positions.bg ../footprinting/hg38.analysisSet.chrom.sizes small_msp_positions.bw
bedGraphToBigWig large_msp_positions.bg ../footprinting/hg38.analysisSet.chrom.sizes large_msp_positions.bw
