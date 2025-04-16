#!/bin/bash
set -euo pipefail

MERGED="merged_GM12878_Liver_corrected_realigned_SLC39A4_region.bam"
samtools merge -o $MERGED ../processed_bams/Liver_SLC39A4_PS00680_map-pb_corrected_realigned_SLC39A4_region.bam ../processed_bams/GM12878_SLC39A4_PS00686_map-pb_corrected_realigned_SLC39A4_region.bam
samtools index $MERGED

# Extract fibers from each cluster for Supplemental Figure
zlists=$(ls clust_pileups/cluster*_zmw.txt)
for cl in $zlists; do
    BNAME=$(basename $cl)
    samtools view -b -N $cl $MERGED > ${BNAME/_zmw.txt/.bam}
    samtools index ${BNAME/_zmw.txt/.bam}
done


