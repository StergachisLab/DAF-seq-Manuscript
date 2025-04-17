#!/bin/bash

conda activate memesuite

# FIMO motifs
fimo --parse-genomic-coord ../../data/JASPAR2022_combined_meme_matrices_homo_sapiens.txt hg38_UBA1.fa

cat fimo_out/fimo.tsv | awk '$9 <= 0.05' | grep MA > motifs/filtered_fimo.tsv
cat motifs/filtered_fimo.tsv | awk '{OFS="\t"}{ print $3,$4,$5,$2,$7,$6 }' > motifs/filtered_fimo.bed
