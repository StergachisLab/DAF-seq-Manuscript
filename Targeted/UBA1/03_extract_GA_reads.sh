#!/bin/bash

samtools view -b -d ST:GA UBA1_PS00685_GM12878_haplotype_corrected.bam > UBA1_region_GA.bam
samtools index UBA1_region_GA.bam
