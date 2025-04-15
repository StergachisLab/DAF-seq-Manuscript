#!/bin/bash

FASTA="/mmfs1/gscratch/stergachislab/assemblies/simple-names/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
FASTA2="/gscratch/stergachislab/assemblies/hg002v1.0.1.fasta"


# K562 alignments
minimap2 -t 40 --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" PS00535_K562_WT.map-pb_corrected.bam) | samtools sort > PS00535_K562_WT.map-pb_corrected_realigned.bam
samtools index PS00535_K562_WT.map-pb_corrected_realigned.bam
minimap2 -t 40 --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" PS00536_K562_5.map-pb_corrected.bam) | samtools sort > PS00536_K562_5.map-pb_corrected_realigned.bam
samtools index PS00536_K562_5.map-pb_corrected_realigned.bam

# HG002 alignments
minimap2 -t 40 -I 8G --MD -Y -y -a -x map-pb $FASTA2 <(samtools fastq -T "*" PS00537_HG002_WT.map-pb_corrected.bam) | samtools sort > PS00537_HG002_WT.map-pb_corrected_realigned.bam
samtools index PS00537_HG002_WT.map-pb_corrected_realigned.bam
minimap2 -t 40 -I 8G --MD -Y -y -a -x map-pb $FASTA2 <(samtools fastq -T "*" PS00538_HG002_5_noWGA.map-pb_corrected.bam) | samtools sort > PS00538_HG002_5_noWGA.map-pb_corrected_realigned.bam
samtools index PS00538_HG002_5_noWGA.map-pb_corrected_realigned.bam
