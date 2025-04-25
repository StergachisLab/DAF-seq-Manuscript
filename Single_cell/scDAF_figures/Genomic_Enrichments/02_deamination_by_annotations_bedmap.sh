#!/bin/bash
set -euo pipefail


dos2unix da_counts_by_pos_all_cells.tsv
gzip da_counts_by_pos_all_cells.tsv


RM_MAPPABLE="hg38.fa.out.repeatmasker.MAPPABLE.sort.bed"
UNMAP="../phasing/unreliable_blacklist_highCov_regions_merged.bed"
AUTO="../phasing/hg38_autosomes.sizes.bed"
PEAKS="../msp_analysis/FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz"

bedtools intersect -a hg38.fa.out.repeatmasker.sort.bed -b $AUTO -u | bedtools intersect -v -a stdin -b $UNMAP | sort-bed - > $RM_MAPPABLE

# Intersect TSV (as if BED) with annotations -----------------------------------------------------------------------
bedmap --ec --echo --echo-map --skip-unmapped $RM_MAPPABLE <(zcat da_counts_by_pos_all_cells.tsv.gz | sort-bed -) | bgzip -@ 10 > RepeatMasker_da_counts_by_pos_all_cells_BEDMAP.txt.gz

# Intersect TSV (as if BED) with HG002 FIRE PEAKS -----------------------------------------------------------------------
bedmap --ec --echo --echo-map --skip-unmapped <(zcat $PEAKS | sort-bed -) <(zcat da_counts_by_pos_all_cells.tsv.gz | sort-bed -) | bgzip -@ 10 > FIRE_Peaks_da_counts_by_pos_all_cells_BEDMAP.txt.gz

# Non RepeatMasker regions -------------------------------------------------------------
bedops --ec -n <(zcat da_counts_by_pos_all_cells.tsv.gz | sort-bed -) $RM_MAPPABLE $UNMAP | grep -v "chrX" | grep -v "chrY" | grep -v "chrEBV" | bgzip -@ 10 > NON_RepeatMasker_da_counts_by_pos_all_cells.bed.gz

