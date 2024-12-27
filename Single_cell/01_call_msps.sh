#!/bin/bash
set -euo pipefail


FT=/mmfs1/gscratch/stergachislab/bin/ft

# Original ----------------
NUC_LEN=60
ADD_C=70
MIN=10

# input is haplotagged BAM
HAP_BAMS=$(ls /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/haplotagged/*HG38_corrected.haplotagged.bam)

LEN=300

for B in $HAP_BAMS;
do
BAM_BASE=$(basename "$B");
NUC_BAM="fibertools_msp/${BAM_BASE/bam/m6A_nuc.bam}";
MSP_BED="fibertools_msp/${BAM_BASE/bam/MSP_${LEN}.bed.gz}";
SINGLE_BED="fibertools_msp/${BAM_BASE/bam/MSP_BED6_single_${LEN}.bed.gz}";
# $FT ddda-to-m6a $B | $FT add-nucleosomes -n $NUC_LEN -c $ADD_C --min-distance-added $MIN - $NUC_BAM;
# samtools index $NUC_BAM;
$FT extract -r -t 20 -x "len(msp)>${LEN}" $NUC_BAM --msp $MSP_BED;
bed12ToBed6 -i $MSP_BED | awk '{OFS="\t"}{if(($3-$2)>100) print }' | bgzip -@ 20 > $SINGLE_BED; # (remove single base end features)
done

