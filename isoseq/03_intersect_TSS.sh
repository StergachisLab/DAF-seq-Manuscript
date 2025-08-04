#!/bin/bash
set -euo pipefail

HG_PEAKS="/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz"
bedmap --ec --echo --echo-map-id --skip-unmapped <(zcat $HG_PEAKS | sort-bed -) <( sort-bed isoseq_counts_giab_hg002.bed ) | sort-bed - | bgzip -@ 10 > isoseq_counts_giab_hg002_Peaks_bedmap.bed.gz


