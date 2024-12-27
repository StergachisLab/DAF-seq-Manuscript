#!/bin/bash
set -euo pipefail


GM_BAM="/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/pileups/GM12878_SLC39A4_PS00686_nuc.bam"
Liver_BAM="/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/pileups/Liver_SLC39A4_PS00680_nuc.bam"
MERGED="merged_GM12878_Liver_nuc.bam"

samtools merge -o $MERGED $GM_BAM $Liver_BAM
samtools index $MERGED


for cFILE in ./clust_pileups/cluster_*.txt; do
    echo $cFILE;
    samtools view -N $cFILE $MERGED
done



