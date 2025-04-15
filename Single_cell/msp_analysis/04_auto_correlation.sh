#!/bin/bash
set -euo pipefail


# autocorrelation of deamination events on collapsed sequences, merged by strand
samples=('PS00718' 'PS00756' 'PS00757' 'PS00758' 'PS00867' 'PS00868' 'PS00869' 'PS00870' 'PS00871' 'PS00872' 'PS00873' 'PS00874')

for s in "${samples[@]}"; do
    /mmfs1/gscratch/stergachislab/bin/ft qc "fibertools_msp/${s}_consensus_BothStrands_HG38_corrected.haplotagged.m6A_nuc.bam" --acf > "${s}_ft_qc.tsv"
done

for s in "${samples[@]}"; do
    cat "${s}_ft_qc.tsv" | grep acf > "${s}_ft_acf.tsv"
done

