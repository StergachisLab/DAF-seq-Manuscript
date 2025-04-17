#!/bin/bash

FASTA="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
FASTA2="hg002v1.0.1.fasta"

# PS00535	K562 nuclei (3.5 m) treated with 0.25 µM DddA WT at 25 ˚C for 10 min in Buffer A + 10 nM ZnCl2, then WGA
# PS00536	K562 nuclei (3.5 m) treated with 0.25 µM DddA 5 at 25 ˚C for 10 min in Buffer A + 10 nM ZnCl2, then WGA
# PS00537	HG002 nuclei (~1-2 m) treated with 0.25 µM DddA WT at 25 ˚C for 10 min in Buffer A + 10 nM ZnCl2, then WGA
# PS00538	HG002 nuclei (~1-2 m) treated with 0.25 µM DddA 5 at 25 ˚C for 10 min in Buffer A + 10 nM ZnCl2

# K562 alignments
minimap2 -t 40 --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" PS00535.phased.bam) | samtools sort > PS00535_K562_WT.map-pb.bam
samtools index PS00535_K562_WT.map-pb.bam
minimap2 -t 40 --MD -Y -y -a -x map-pb $FASTA <(samtools fastq -T "*" PS00536.phased.bam) | samtools sort > PS00536_K562_5.map-pb.bam
samtools index PS00536_K562_5.map-pb.bam

# HG002 alignments
minimap2 -t 40 -I 8G --MD -Y -y -a -x map-pb $FASTA2 <(samtools fastq -T "*" PS00537.phased.bam) | samtools sort > PS00537_HG002_WT.map-pb.bam
samtools index PS00537_HG002_WT.map-pb.bam
minimap2 -t 40 -I 8G --MD -Y -y -a -x map-pb $FASTA2 <(samtools fastq -T "*" PS00538.phased.bam) | samtools sort > PS00538_HG002_5_noWGA.map-pb.bam
samtools index PS00538_HG002_5_noWGA.map-pb.bam


