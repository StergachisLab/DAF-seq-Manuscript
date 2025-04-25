#!/bin/bash
set -euo pipefail

FASTA="hg38.fa"

# phasing with single VCF
VCF="analysis-HG002_WGS-4142-deepvariant.phased.vcf.gz"

BAMS=$(ls ../collapse/consensus_bams/PS*G38_corrected.bam)

for b in $BAMS; do
bname=$(basename $b)
whatshap haplotag -o "haplotagged/${bname/bam/haplotagged.bam}" --reference $FASTA --ignore-read-groups --output-haplotag-list "haplotagged/${bname/bam/haplotag_list.txt.gz}" $VCF $b;
samtools index "haplotagged/${bname/bam/haplotagged.bam}";
done

