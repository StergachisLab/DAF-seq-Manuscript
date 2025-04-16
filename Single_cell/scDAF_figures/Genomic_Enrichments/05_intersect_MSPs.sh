#!/bin/bash
set -euo pipefail


# only include TSSs that are in HG002 FIRE peaks (active TSSs) with reliable coverage as determined in MSP analysis 
FIRE_PEAKS="/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz"
bedtools intersect -a gencodev45_Ensembl_canonical_TSS.bed -b $FIRE_PEAKS -u | sort-bed - > gencodev45_Ensembl_canonical_TSS_HG002_FIRE_noUnreliableCoverage.bed

# split by TSS strand
cat gencodev45_Ensembl_canonical_TSS_HG002_FIRE_noUnreliableCoverage.bed | awk '$6 == "+"' > gencodev45_Ensembl_canonical_TSS_HG002_FIRE_noUnreliableCoverage_PLUS_Strand.bed
cat gencodev45_Ensembl_canonical_TSS_HG002_FIRE_noUnreliableCoverage.bed | awk '$6 == "-"' > gencodev45_Ensembl_canonical_TSS_HG002_FIRE_noUnreliableCoverage_MINUS_Strand.bed


# merge all DAF MSP BEDs from all cells and all haplotype-strands
DAF_MSP_BED=DAF_all_MSPs_all_cells_and_strands.bed.gz
BEDS=$(ls /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/fibertools_msp/*.haplotagged.MSP_BED6_single_150.bed.gz)
zcat $BEDS | sort-bed - | bgzip -@ 10 > $DAF_MSP_BED


# filter BED for +/- 2 Kb from each TSS (overlap 1 bp)
bedops --ec --range 2000 --element-of 1 <(zcat $DAF_MSP_BED) gencodev45_Ensembl_canonical_TSS_HG002_FIRE_noUnreliableCoverage_PLUS_Strand.bed > DAF_tss_PLUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bed
bedtools genomecov -bg -i DAF_tss_PLUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bed -g /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes > DAF_tss_PLUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bg

bedops --ec --range 2000 --element-of 1 <(zcat $DAF_MSP_BED) gencodev45_Ensembl_canonical_TSS_HG002_FIRE_noUnreliableCoverage_MINUS_Strand.bed > DAF_tss_MINUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bed
bedtools genomecov -bg -i DAF_tss_MINUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bed -g /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes > DAF_tss_MINUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bg

# Repeat by FIRE Peaks
bedops --ec --range 2000 --element-of 1 <(zcat $DAF_MSP_BED) <(zcat $FIRE_PEAKS | sort-bed -) > DAF_MSPs_FIRE_peaks_2kb_all_cells_and_strands.bed
bedtools genomecov -bg -i DAF_MSPs_FIRE_peaks_2kb_all_cells_and_strands.bed -g /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes > DAF_MSPs_FIRE_peaks_2kb_all_cells_and_strands.bg


# Fiber-seq ---------------------------------------------------------
FT=/mmfs1/gscratch/stergachislab/bin/ft
FS_BAM=/mmfs1/gscratch/stergachislab/mvollger/projects/FIREv2.0/results/HG002/fire/HG002.fire.cram
MSP_BED=FS_FIRE_MSP_BED12.bed.gz
SINGLE_BED=FS_FIRE_MSP_BED6_single.bed.gz
$FT extract -r -t 20 -x "len(msp)>150" $FS_BAM --msp $MSP_BED
bed12ToBed6 -i $MSP_BED | awk '{OFS="\t"}{if(($3-$2)>100) print }' | sort-bed - | bgzip -@ 20 > $SINGLE_BED; # (remove single base end features)

bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) gencodev45_Ensembl_canonical_TSS_HG002_FIRE_noUnreliableCoverage_PLUS_Strand.bed > FS_tss_PLUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bed
bedtools genomecov -bg -i FS_tss_PLUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bed -g /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes > FS_tss_PLUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bg
bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) gencodev45_Ensembl_canonical_TSS_HG002_FIRE_noUnreliableCoverage_MINUS_Strand.bed > FS_tss_MINUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bed
bedtools genomecov -bg -i FS_tss_MINUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bed -g /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes > FS_tss_MINUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bg

bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) <(zcat $FIRE_PEAKS | sort-bed -) > FS_MSPs_FIRE_peaks_2kb_all_cells_and_strands.bed
bedtools genomecov -bg -i FS_MSPs_FIRE_peaks_2kb_all_cells_and_strands.bed -g /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes > FS_MSPs_FIRE_peaks_2kb_all_cells_and_strands.bg


# MSP intersect by RM categories (plus FIRE peaks and non-RM & non-FIRE) ---------------------------------------------------------

# REQUIRE >=50% overlap of MSP with region to avoid double counting MSPs !!
RM_MAPPABLE=hg38.fa.out.repeatmasker.MAPPABLE.sort.bed
UNMAP="/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/unreliable_blacklist_highCov_regions_merged.bed"
AUTO="/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/hg38_autosomes.sizes.bed"
NON_RM_REGIONS=NON_RM_regions.bed

bedops --ec -d $AUTO $RM_MAPPABLE > $NON_RM_REGIONS

bedmap --ec --echo --echo-map --skip-unmapped --fraction-either 0.5 $RM_MAPPABLE <(zcat $DAF_MSP_BED | sort-bed -) | gzip > DAF_RepeatMasker_MSP_BEDMAP.txt.gz
bedmap --ec --echo --echo-map --skip-unmapped --fraction-either 0.5 <(zcat $FIRE_PEAKS | sort-bed -) <(zcat $DAF_MSP_BED | sort-bed -) | gzip > DAF_FIRE_Peaks_MSP_BEDMAP.txt.gz
bedmap --ec --echo --echo-map --skip-unmapped --fraction-either 0.5 $RM_MAPPABLE <(zcat $SINGLE_BED | sort-bed -) | gzip > FS_RepeatMasker_MSP_BEDMAP.txt.gz
bedmap --ec --echo --echo-map --skip-unmapped --fraction-either 0.5 <(zcat $FIRE_PEAKS | sort-bed -) <(zcat $SINGLE_BED | sort-bed -) | gzip > FS_FIRE_Peaks_MSP_BEDMAP.txt.gz

# Non RepeatMasker regions
bedops --ec -n <(zcat $DAF_MSP_BED | sort-bed -) $RM_MAPPABLE $UNMAP | grep -v "chrX" | grep -v "chrY" | grep -v "chrEBV" | bgzip -@ 10 > DAF_NON_RepeatMasker_MSP.bed.gz
bedmap --ec --echo --echo-map --skip-unmapped --fraction-either 0.5 $NON_RM_REGIONS <(zcat DAF_NON_RepeatMasker_MSP.bed.gz | sort-bed -) | gzip > DAF_NON_RepeatMasker_MSP_BEDMAP.txt.gz
bedops --ec -n <(zcat $SINGLE_BED | sort-bed -) $RM_MAPPABLE $UNMAP | grep -v "chrX" | grep -v "chrY" | grep -v "chrEBV" | bgzip -@ 10 > FS_NON_RepeatMasker_MSP.bed.gz
bedmap --ec --echo --echo-map --skip-unmapped --fraction-either 0.5 $NON_RM_REGIONS <(zcat FS_NON_RepeatMasker_MSP.bed.gz | sort-bed -) | gzip > FS_NON_RepeatMasker_MSP_BEDMAP.txt.gz

# All mappable autosome regions (RM and FIRE included)
MAP_ALL="mappable_autosome_regions.bed"
# bedops --ec -d $AUTO $UNMAP > $MAP_ALL
bedmap --ec --echo --echo-map --skip-unmapped --fraction-either 0.5 $MAP_ALL <(zcat $DAF_MSP_BED | sort-bed -) | gzip > DAF_ALL_MSP_BEDMAP.txt.gz
bedmap --ec --echo --echo-map --skip-unmapped --fraction-either 0.5 $MAP_ALL <(zcat $SINGLE_BED | sort-bed -) | gzip > FS_ALL_MSP_BEDMAP.txt.gz
