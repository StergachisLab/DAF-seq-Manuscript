#!/bin/bash

# conda activate fibertools

ft ddda-to-m6a hap2_T_only.bam | ft add-nucleosomes - hap2_T_only_nuc.bam
samtools index hap2_T_only_nuc.bam
ft extract hap2_T_only_nuc.bam --all hap2_T_only_ft_extract_all.bed
