#!/bin/bash

alpha_ft=/mmfs1/gscratch/stergachislab/mvollger/projects/fibertools-rs/target/release/ft
mod_bam=/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/Napa_WASF1/haplotype_correct/NAPA_PS00626_haplotype_corrected.bam

$alpha_ft ddda-to-m6a $mod_bam NAPA_m6a.bam
samtools index NAPA_m6a.bam

ft add-nucleosomes -v NAPA_m6a.bam NAPA_nuc.bam
samtools index nuc.bam

samtools view -b -d ST:CT NAPA_nuc.bam > NAPA_nuc_CT.bam
samtools view -b -d ST:GA NAPA_nuc.bam > NAPA_nuc_GA.bam
samtools index NAPA_nuc_CT.bam
samtools index NAPA_nuc_GA.bam


# manually run ft footprint for CTCF using individual modules
ft footprint nuc_CT.bam motifs/ctcf_motif.bed modules/ctcf_modules.yaml > ft_out/ft_CT_ModCTCF.bed
ft footprint nuc_GA.bam motifs/ctcf_motif.bed modules/ctcf_modules.yaml > ft_out/ft_GA_ModCTCF.bed

