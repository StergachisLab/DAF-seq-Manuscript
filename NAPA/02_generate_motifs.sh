#!/bin/bash

# Creating hg38 fasta for use in FIMO
bedtools getfasta -fi /gscratch/stergachislab/swansoe/ref_genomes/hg38.fa -bed napa_promoter.bed -fo hg38_NAPA_promoter.fa


# Sourcing env from Adriana
source /gscratch/stergachislab/asedeno/tools/miniconda3/bin/activate
conda activate memesuite

# FIMO motifs
fimo --parse-genomic-coord JASPAR2022_combined_meme_matrices_homo_sapiens.txt hg38_NAPA_promoter.fa


