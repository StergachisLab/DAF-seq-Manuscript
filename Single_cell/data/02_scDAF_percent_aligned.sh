#!/bin/bash
set -euo pipefail


echo "Sample" "Mapped" "Unmapped" "BAM" | awk '{OFS="\t"}{ print $1,$2,$3,$4 }' > percent_mapped.tsv

BAMS=$(ls PS*.map-pb_HG38.bam)
for b in $BAMS; do
    SAMPLE=$(basename $b | cut -d '_' -f1)
    MAPPED=$(samtools view -F 2308 -c $b)
    UNMAPPED=$(samtools view -f 4 -c $b)
    echo $SAMPLE $MAPPED $UNMAPPED $b | awk '{OFS="\t"}{ print $1,$2,$3,$4 }' >> percent_mapped.tsv
done
