#!/bin/bash
set -euo pipefail

FASTA="/mmfs1/gscratch/stergachislab/assemblies/simple-names/hg38.fa"

# phasing with single VCF
VCF="analysis-HG002_WGS-4142-deepvariant.phased.vcf.gz"

BAM_DIR="/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/collapse/consensus_bams"
CT_BAMS=('PS00718_consensus_CT_HG38_corrected.bam' 'PS00756_consensus_CT_HG38_corrected.bam' 'PS00757_consensus_CT_HG38_corrected.bam' 'PS00758_consensus_CT_HG38_corrected.bam')
GA_BAMS=('PS00718_consensus_GA_HG38_corrected.bam' 'PS00756_consensus_GA_HG38_corrected.bam' 'PS00757_consensus_GA_HG38_corrected.bam' 'PS00758_consensus_GA_HG38_corrected.bam')

for b in "${CT_BAMS[@]}"; do
whatshap haplotag -o "haplotagged/${b/bam/haplotagged.bam}" --reference $FASTA --ignore-read-groups --output-haplotag-list "haplotagged/${b/bam/haplotag_list.txt.gz}" $VCF $BAM_DIR/$b;
samtools index "haplotagged/${b/bam/haplotagged.bam}";
done

for b in "${GA_BAMS[@]}"; do
whatshap haplotag -o "haplotagged/${b/bam/haplotagged.bam}" --reference $FASTA --ignore-read-groups --output-haplotag-list "haplotagged/${b/bam/haplotag_list.txt.gz}" $VCF $BAM_DIR/$b;
samtools index "haplotagged/${b/bam/haplotagged.bam}";
done


# # Filter the collapsed, haplotagged reads by HG002 Fiber-seq unreliable coverage regions & sex chromosomes

# UNR_BED="/mmfs1/gscratch/stergachislab/mvollger/projects/FIREv2.0/results/HG002/coverage/unreliable-coverage-regions.bed.gz"

# HAP_BAMS=$( ls haplotagged/*HG38_corrected.haplotagged.bam)

# for b in $HAP_BAMS; do
# samtools view -e 'rname != "chrX" && rname != "chrY"' -@ 15 -b -L $UNR_BED -U ${b/bam/NoUnreliable_filtered.bam} $b > ${b/bam/UnreliableOnly.bam};
# samtools index ${b/bam/NoUnreliable_filtered.bam};
# samtools index ${b/bam/UnreliableOnly.bam};
# done

