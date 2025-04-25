#!/bin/bash
set -euo pipefail


# Phased "contigs", unphased "contigs", and total "contigs"

FYAK=HG003.father.yak
MYAK=HG004.mother.yak

BAMS=$(ls ../haplotagged/PS00*.bam)

for FILE in $BAMS; do
    B=$(basename $FILE)
    HAP_LIST=${FILE/haplotagged.bam/haplotag_list.txt.gz}
    HL_BASE=$(basename $HAP_LIST)
    PHASED=${HL_BASE/txt.gz/phased_reads.txt}
    UN=${HL_BASE/txt.gz/UNphased_reads.txt}
    UN_FA=${B/bam/UNPHASED.fa}
    PHASED_FA=${B/bam/PHASED.fa}
    ALL_FA=${B/bam/ALL.fa}
    sed 1d <(zcat $HAP_LIST) | awk -F "\t" '{ if($2 != "none") print $1 }' > $PHASED
    sed 1d <(zcat $HAP_LIST) | awk -F "\t" '{ if($2 == "none") print $1 }' > $UN
    samtools view -b -N $UN $FILE | samtools fasta > $UN_FA
    samtools view -b -N $PHASED $FILE | samtools fasta > $PHASED_FA
    samtools view -b $FILE | samtools fasta > $ALL_FA
    yak trioeval -t32 $FYAK $MYAK $ALL_FA > ${ALL_FA/.fa/.trioeval.tbl} &
    yak trioeval -t32 $FYAK $MYAK $PHASED_FA > ${PHASED_FA/.fa/.trioeval.tbl} &
    yak trioeval -t32 $FYAK $MYAK $UN_FA > ${UN_FA/.fa/.trioeval.tbl} &
done
