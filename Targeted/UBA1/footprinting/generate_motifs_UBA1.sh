#!/bin/bash

# Creating hg38 fasta for use in FIMO
bedtools getfasta -fi /gscratch/stergachislab/swansoe/ref_genomes/hg38.fa -bed ../FIRE_narrowPeaks_GM12878.bed -fo hg38_UBA1.fa

# Sourcing env from Adriana
source /gscratch/stergachislab/asedeno/tools/miniconda3/bin/activate
conda activate memesuite

# FIMO motifs
fimo --parse-genomic-coord /mmfs1/gscratch/stergachislab/swansoe/projects/GM12878_scATAC/initial_RR/Fiber-seq/TF/JASPAR2022_combined_meme_matrices_homo_sapiens.txt hg38_UBA1.fa

