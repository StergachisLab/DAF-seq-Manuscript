#!/bin/bash
set -euo pipefail

# get the lengths of raw sequenced reads for each cell

OUT_FILE="raw_read_lengths.tsv"
if [[ -e $OUT_FILE ]]; then rm $OUT_FILE; fi

BAMS=$(ls ../data/PS*.bam)

for B in $BAMS; do
    NAME=$(basename $B | cut -d'_' -f1)
    samtools view ${B} | awk -v name=$NAME '{OFS="\t"}{ print length($10), name }' >> $OUT_FILE
done

bgzip -@ 20 $OUT_FILE
