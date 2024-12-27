#!/bin/bash
set -euo pipefail


# autocorrelation of deamination events on collapsed sequences, merged by strand

samples=('PS00718' 'PS00756' 'PS00757' 'PS00758' 'PS00758-4')

for s in "${samples[@]}"; do
    samtools merge -o - "fibertools_msp/${s}_consensus_CT_HG38_corrected.haplotagged.m6A_nuc.bam" \
    "fibertools_msp/${s}_consensus_CT_HG38_corrected.haplotagged.m6A_nuc.bam" | /mmfs1/gscratch/stergachislab/bin/ft qc \
    --acf > "${s}_ft_qc.tsv"
done

for s in "${samples[@]}"; do
    cat "${s}_ft_qc.tsv" | grep acf > "${s}_ft_acf.tsv"
done

