#!/bin/bash
set -euo pipefail

FASTA="/mmfs1/gscratch/stergachislab/assemblies/simple-names/hg38.fa"

# phasing with single VCF
VCF="analysis-HG002_WGS-4142-deepvariant.phased.vcf.gz"

BAMS=$(ls /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/collapse/consensus_bams/PS*G38_corrected.bam)

for b in $BAMS; do
bname=$(basename $b)
whatshap haplotag -o "haplotagged/${bname/bam/haplotagged.bam}" --reference $FASTA --ignore-read-groups --output-haplotag-list "haplotagged/${bname/bam/haplotag_list.txt.gz}" $VCF $b;
samtools index "haplotagged/${bname/bam/haplotagged.bam}";
done

