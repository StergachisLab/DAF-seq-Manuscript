#!/bin/bash
set -euo pipefail


# conda activate isoseq

pbmm2 align --preset ISOSEQ --sort giab_na24385.hifi_reads.lima.0--0.lima.IsoSeqX_bc02_5p--IsoSeqX_3p.refined.bam /gscratch/stergachislab/assemblies/simple-names/hg38.fa giab_hg002_isoseq_mapped.bam

COLLAPSED="giab_hg002_isoseq_mapped.collapsed.gff"
isoseq collapse --do-not-collapse-extra-5exons giab_hg002_isoseq_mapped.bam $COLLAPSED

COL_SORTED="GM12878_mas_collapsed.sorted.gff"
cat $COLLAPSED | grep "#" > $COL_SORTED
cat $COLLAPSED | grep -v "#" | sort -k1,1 >> $COL_SORTED

