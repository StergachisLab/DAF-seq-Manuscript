#!/bin/bash
set -euo pipefail


GENCODE=../../msp_analysis/gencodev45_Ensembl_canonical_TSS.bed
FT=ft
SIZES=../../msp_analysis/hg38.analysisSet.chrom.sizes
REF=hg38.analysisSet


# Fiber-seq ---------------------------------------------------------

# FIRE v0.1.1 ---------------------------------------------

# HG002
SNAME="HG002"
FS_BAM="../../data/HG002.fire.cram"
UNFILT_FIRE="../../data/FDR-FIRE-peaks.bed.gz"
FIRE_PEAKS="tss_enrichment_intermediate/${SNAME}_FIRE_peaks_filtered.bed.gz"
zcat $UNFILT_FIRE | awk '$NF == "true"' | bgzip -@ 10 > $FIRE_PEAKS
# only include TSSs that are in FIRE peaks (active TSSs) with reliable coverage as determined in MSP analysis 
ALL_TSS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE.bed"
bedtools intersect -a $GENCODE -b $FIRE_PEAKS -u | sort-bed - > $ALL_TSS
# split by TSS strand
TSS_PLUS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE_PLUS_Strand.bed"
TSS_MINUS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE_MINUS_Strand.bed"
cat  $ALL_TSS | awk '$6 == "+"' > $TSS_PLUS
cat  $ALL_TSS | awk '$6 == "-"' > $TSS_MINUS
MSP_BED="tss_enrichment_intermediate/${SNAME}_FS_FIRE_MSP_BED12.bed.gz"
SINGLE_BED="tss_enrichment_intermediate/${SNAME}_FS_FIRE_MSP_BED6_single.bed.gz"
$FT extract -r -t 20 -x "len(msp)>150" $FS_BAM --msp $MSP_BED
bed12ToBed6 -i $MSP_BED | awk '{OFS="\t"}{if(($3-$2)>100) print }' | sort-bed - | bgzip -@ 20 > $SINGLE_BED; # (remove single base end features)
PLUS_BED="tss_enrichment_intermediate/${SNAME}_FS_tss_PLUS_2kb_MSPs_FIRE_peaks.bed"
MINUS_BED="tss_enrichment_intermediate/${SNAME}_FS_tss_MINUS_2kb_MSPs_FIRE_peaks.bed"
PLUS_BG="tss_enrichment_intermediate/${SNAME}_FS_tss_PLUS_2kb_MSPs_FIRE_peaks.bg"
MINUS_BG="tss_enrichment_intermediate/${SNAME}_FS_tss_MINUS_2kb_MSPs_FIRE_peaks.bg"
bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) $TSS_PLUS > $PLUS_BED
bedtools genomecov -bg -i $PLUS_BED -g $SIZES > $PLUS_BG
bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) $TSS_MINUS > $MINUS_BED
bedtools genomecov -bg -i $MINUS_BED -g $SIZES > $MINUS_BG
MSP_BED="tss_enrichment_intermediate/${SNAME}_FS_MSPs_FIRE_peaks_2kb.bed"
MSP_BG="tss_enrichment_intermediate/${SNAME}_FS_MSPs_FIRE_peaks_2kb.bg"
bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) <(zcat $FIRE_PEAKS | sort-bed -) > $MSP_BED
bedtools genomecov -bg -i $MSP_BED -g $SIZES > $MSP_BG
$FT extract -r -t 30 $FS_BAM --all - | sed 1d | datamash sum 14 sum 15 | awk '{OFS="\t"}{ print }' > "tss_enrichment_intermediate/${SNAME}_m6A_prop.tsv"


SAMPLES=("GM12878")
for SNAME in "${SAMPLES[@]}"; do
    FS_BAM="../../data/${SNAME}/${SNAME}-fire-v0.1.1-filtered.cram"
    UNFILT_FIRE="../../data/${SNAME}/${SNAME}-fire-v0.1.1-peaks.bed.gz"
    FIRE_PEAKS="tss_enrichment_intermediate/${SNAME}_FIRE_peaks_filtered.bed.gz"
    zcat $UNFILT_FIRE | awk '$NF == "true"' | bgzip -@ 10 > $FIRE_PEAKS
    # only include TSSs that are in FIRE peaks (active TSSs) with reliable coverage as determined in MSP analysis 
    ALL_TSS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE.bed"
    bedtools intersect -a $GENCODE -b $FIRE_PEAKS -u | sort-bed - > $ALL_TSS
    # split by TSS strand
    TSS_PLUS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE_PLUS_Strand.bed"
    TSS_MINUS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE_MINUS_Strand.bed"
    cat  $ALL_TSS | awk '$6 == "+"' > $TSS_PLUS
    cat  $ALL_TSS | awk '$6 == "-"' > $TSS_MINUS
    MSP_BED="tss_enrichment_intermediate/${SNAME}_FS_FIRE_MSP_BED12.bed.gz"
    SINGLE_BED="tss_enrichment_intermediate/${SNAME}_FS_FIRE_MSP_BED6_single.bed.gz"
    $FT extract -r -t 20 -x "len(msp)>150" $FS_BAM --msp $MSP_BED
    bed12ToBed6 -i $MSP_BED | awk '{OFS="\t"}{if(($3-$2)>100) print }' | sort-bed - | bgzip -@ 20 > $SINGLE_BED; # (remove single base end features)
    PLUS_BED="tss_enrichment_intermediate/${SNAME}_FS_tss_PLUS_2kb_MSPs_FIRE_peaks.bed"
    MINUS_BED="tss_enrichment_intermediate/${SNAME}_FS_tss_MINUS_2kb_MSPs_FIRE_peaks.bed"
    PLUS_BG="tss_enrichment_intermediate/${SNAME}_FS_tss_PLUS_2kb_MSPs_FIRE_peaks.bg"
    MINUS_BG="tss_enrichment_intermediate/${SNAME}_FS_tss_MINUS_2kb_MSPs_FIRE_peaks.bg"
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) $TSS_PLUS > $PLUS_BED
    bedtools genomecov -bg -i $PLUS_BED -g $SIZES > $PLUS_BG
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) $TSS_MINUS > $MINUS_BED
    bedtools genomecov -bg -i $MINUS_BED -g $SIZES > $MINUS_BG
    MSP_BED="tss_enrichment_intermediate/${SNAME}_FS_MSPs_FIRE_peaks_2kb.bed"
    MSP_BG="tss_enrichment_intermediate/${SNAME}_FS_MSPs_FIRE_peaks_2kb.bg"
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) <(zcat $FIRE_PEAKS | sort-bed -) > $MSP_BED
    bedtools genomecov -bg -i $MSP_BED -g $SIZES > $MSP_BG
    $FT extract -r -t 30 $FS_BAM --all - | sed 1d | datamash sum 14 sum 15 | awk '{OFS="\t"}{ print }' > "tss_enrichment_intermediate/${SNAME}_m6A_prop.tsv"
done


# FIRE v0.1 ---------------------------------------------
SAMPLES=("CHM13")
for SNAME in "${SAMPLES[@]}"; do
    FS_BAM="../../data/${SNAME}/${SNAME}-fire-v0.1-filtered.cram"
    UNFILT_FIRE="../../data/${SNAME}/${SNAME}-fire-v0.1-peaks.bed.gz"
    FIRE_PEAKS="tss_enrichment_intermediate/${SNAME}_FIRE_peaks_filtered.bed.gz"
    zcat $UNFILT_FIRE | awk '$NF == "true"' | bgzip -@ 10 > $FIRE_PEAKS
    # only include TSSs that are in FIRE peaks (active TSSs) with reliable coverage as determined in MSP analysis 
    ALL_TSS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE.bed"
    bedtools intersect -a $GENCODE -b $FIRE_PEAKS -u | sort-bed - > $ALL_TSS
    # split by TSS strand
    TSS_PLUS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE_PLUS_Strand.bed"
    TSS_MINUS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE_MINUS_Strand.bed"
    cat  $ALL_TSS | awk '$6 == "+"' > $TSS_PLUS
    cat  $ALL_TSS | awk '$6 == "-"' > $TSS_MINUS
    MSP_BED="tss_enrichment_intermediate/${SNAME}_FS_FIRE_MSP_BED12.bed.gz"
    SINGLE_BED="tss_enrichment_intermediate/${SNAME}_FS_FIRE_MSP_BED6_single.bed.gz"
    $FT extract -r -t 20 -x "len(msp)>150" $FS_BAM --msp $MSP_BED
    bed12ToBed6 -i $MSP_BED | awk '{OFS="\t"}{if(($3-$2)>100) print }' | sort-bed - | bgzip -@ 20 > $SINGLE_BED; # (remove single base end features)
    PLUS_BED="tss_enrichment_intermediate/${SNAME}_FS_tss_PLUS_2kb_MSPs_FIRE_peaks.bed"
    MINUS_BED="tss_enrichment_intermediate/${SNAME}_FS_tss_MINUS_2kb_MSPs_FIRE_peaks.bed"
    PLUS_BG="tss_enrichment_intermediate/${SNAME}_FS_tss_PLUS_2kb_MSPs_FIRE_peaks.bg"
    MINUS_BG="tss_enrichment_intermediate/${SNAME}_FS_tss_MINUS_2kb_MSPs_FIRE_peaks.bg"
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) $TSS_PLUS > $PLUS_BED
    bedtools genomecov -bg -i $PLUS_BED -g $SIZES > $PLUS_BG
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) $TSS_MINUS > $MINUS_BED
    bedtools genomecov -bg -i $MINUS_BED -g $SIZES > $MINUS_BG
    MSP_BED="tss_enrichment_intermediate/${SNAME}_FS_MSPs_FIRE_peaks_2kb.bed"
    MSP_BG="tss_enrichment_intermediate/${SNAME}_FS_MSPs_FIRE_peaks_2kb.bg"
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) <(zcat $FIRE_PEAKS | sort-bed -) > $MSP_BED
    bedtools genomecov -bg -i $MSP_BED -g $SIZES > $MSP_BG
    $FT extract -r -t 30 $FS_BAM --all - | sed 1d | datamash sum 14 sum 15 | awk '{OFS="\t"}{ print }' > "tss_enrichment_intermediate/${SNAME}_m6A_prop.tsv"
done

# From FIRE publication ---------------------------------------------

SAMPLES=("PS00272" "PS00321" "PS00327" "PS00381" "PS00382" "PS00383" "PS00384" "ST001-liver" "ST001-lung" "PS00338_COLO829BL_1" "PS00356_COLO829BL_2" "COLO_T_2_PS00_418_451_488")
for SNAME in "${SAMPLES[@]}"; do
    FS_BAM="../../data/${SNAME}/fire/${SNAME}.fire.cram"
    UNFILT_FIRE="../../data/results/${SNAME}/FDR-peaks/FDR-FIRE-peaks.bed.gz"
    FIRE_PEAKS="tss_enrichment_intermediate/${SNAME}_FIRE_peaks_filtered.bed.gz"
    zcat $UNFILT_FIRE | awk '$NF == "true"' | bgzip -@ 10 > $FIRE_PEAKS
    # only include TSSs that are in FIRE peaks (active TSSs) with reliable coverage as determined in MSP analysis 
    ALL_TSS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE.bed"
    bedtools intersect -a $GENCODE -b $FIRE_PEAKS -u | sort-bed - > $ALL_TSS
    # split by TSS strand
    TSS_PLUS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE_PLUS_Strand.bed"
    TSS_MINUS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE_MINUS_Strand.bed"
    cat  $ALL_TSS | awk '$6 == "+"' > $TSS_PLUS
    cat  $ALL_TSS | awk '$6 == "-"' > $TSS_MINUS
    MSP_BED="tss_enrichment_intermediate/${SNAME}_FS_FIRE_MSP_BED12.bed.gz"
    SINGLE_BED="tss_enrichment_intermediate/${SNAME}_FS_FIRE_MSP_BED6_single.bed.gz"
    $FT extract -r -t 20 -x "len(msp)>150" $FS_BAM --msp $MSP_BED
    bed12ToBed6 -i $MSP_BED | awk '{OFS="\t"}{if(($3-$2)>100) print }' | sort-bed - | bgzip -@ 20 > $SINGLE_BED; # (remove single base end features)
    PLUS_BED="tss_enrichment_intermediate/${SNAME}_FS_tss_PLUS_2kb_MSPs_FIRE_peaks.bed"
    MINUS_BED="tss_enrichment_intermediate/${SNAME}_FS_tss_MINUS_2kb_MSPs_FIRE_peaks.bed"
    PLUS_BG="tss_enrichment_intermediate/${SNAME}_FS_tss_PLUS_2kb_MSPs_FIRE_peaks.bg"
    MINUS_BG="tss_enrichment_intermediate/${SNAME}_FS_tss_MINUS_2kb_MSPs_FIRE_peaks.bg"
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) $TSS_PLUS > $PLUS_BED
    bedtools genomecov -bg -i $PLUS_BED -g $SIZES > $PLUS_BG
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) $TSS_MINUS > $MINUS_BED
    bedtools genomecov -bg -i $MINUS_BED -g $SIZES > $MINUS_BG
    MSP_BED="tss_enrichment_intermediate/${SNAME}_FS_MSPs_FIRE_peaks_2kb.bed"
    MSP_BG="tss_enrichment_intermediate/${SNAME}_FS_MSPs_FIRE_peaks_2kb.bg"
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) <(zcat $FIRE_PEAKS | sort-bed -) > $MSP_BED
    bedtools genomecov -bg -i $MSP_BED -g $SIZES > $MSP_BG
    $FT extract -r -t 30 $FS_BAM --all - | sed 1d | datamash sum 14 sum 15 | awk '{OFS="\t"}{ print }' > "tss_enrichment_intermediate/${SNAME}_m6A_prop.tsv"
done


# Colon descending (SMaHT) ----------------------------------------------------------------------------
SNAME="PS30743"
FS_BAM="../../data/PS30743/PS30743-fire-v0.1.2-filtered.cram"
UNFILT_FIRE="../../data/results/pipelines/fire-run/results/PS30743/PS30743-fire-v0.1.2-peaks.bed.gz"
FIRE_PEAKS="tss_enrichment_intermediate/${SNAME}_FIRE_peaks_filtered.bed.gz"
zcat $UNFILT_FIRE | awk '$NF == "true"' | bgzip -@ 10 > $FIRE_PEAKS
# only include TSSs that are in FIRE peaks (active TSSs) with reliable coverage as determined in MSP analysis 
ALL_TSS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE.bed"
bedtools intersect -a $GENCODE -b $FIRE_PEAKS -u | sort-bed - > $ALL_TSS
# split by TSS strand
TSS_PLUS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE_PLUS_Strand.bed"
TSS_MINUS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE_MINUS_Strand.bed"
cat  $ALL_TSS | awk '$6 == "+"' > $TSS_PLUS
cat  $ALL_TSS | awk '$6 == "-"' > $TSS_MINUS
MSP_BED="tss_enrichment_intermediate/${SNAME}_FS_FIRE_MSP_BED12.bed.gz"
SINGLE_BED="tss_enrichment_intermediate/${SNAME}_FS_FIRE_MSP_BED6_single.bed.gz"
$FT extract -r -t 20 -x "len(msp)>150" $FS_BAM --msp $MSP_BED
bed12ToBed6 -i $MSP_BED | awk '{OFS="\t"}{if(($3-$2)>100) print }' | sort-bed - | bgzip -@ 20 > $SINGLE_BED; # (remove single base end features)
PLUS_BED="tss_enrichment_intermediate/${SNAME}_FS_tss_PLUS_2kb_MSPs_FIRE_peaks.bed"
MINUS_BED="tss_enrichment_intermediate/${SNAME}_FS_tss_MINUS_2kb_MSPs_FIRE_peaks.bed"
PLUS_BG="tss_enrichment_intermediate/${SNAME}_FS_tss_PLUS_2kb_MSPs_FIRE_peaks.bg"
MINUS_BG="tss_enrichment_intermediate/${SNAME}_FS_tss_MINUS_2kb_MSPs_FIRE_peaks.bg"
bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) $TSS_PLUS > $PLUS_BED
bedtools genomecov -bg -i $PLUS_BED -g $SIZES > $PLUS_BG
bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) $TSS_MINUS > $MINUS_BED
bedtools genomecov -bg -i $MINUS_BED -g $SIZES > $MINUS_BG
MSP_BED="tss_enrichment_intermediate/${SNAME}_FS_MSPs_FIRE_peaks_2kb.bed"
MSP_BG="tss_enrichment_intermediate/${SNAME}_FS_MSPs_FIRE_peaks_2kb.bg"
bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) <(zcat $FIRE_PEAKS | sort-bed -) > $MSP_BED
bedtools genomecov -bg -i $MSP_BED -g $SIZES > $MSP_BG
$FT extract -r -t 30 $FS_BAM --all - | sed 1d | datamash sum 14 sum 15 | awk '{OFS="\t"}{ print }' > "tss_enrichment_intermediate/${SNAME}_m6A_prop.tsv"


FS_PROP=tss_enrichment_intermediate/FS_m6A_prop.tsv
for TSV in $(ls tss_enrichment_intermediate/*.tsv); do
    BNAME=$(basename $TSV)
    SNAME=${BNAME/_m6A_prop.tsv/}
    P=$(cat $TSV | awk '{ print $2/$1 }')
    printf '%s\t%s\n' $SNAME $P >> $FS_PROP
done




# Single-cell DAF-seq ----------------------------------------------------------------------------

# only include TSSs that are in HG002 peaks (active TSSs) with reliable coverage as determined in MSP analysis 
FIRE_PEAKS="../../msp_analysis/FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz"
ALL_TSS="gencodev45_Ensembl_canonical_TSS_HG002_FIRE_noUnreliableCoverage.bed"
bedtools intersect -a $GENCODE -b $FIRE_PEAKS -u | sort-bed - > $ALL_TSS

SAMPLES=('PS00718' 'PS00756' 'PS00757' 'PS00758' 'PS00867' 'PS00868' 'PS00869' 'PS00870' 'PS00871' 'PS00872' 'PS00873' 'PS00874')
for SNAME in "${SAMPLES[@]}"; do
    TSS_PLUS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE_PLUS_Strand.bed"
    TSS_MINUS="tss_enrichment_intermediate/${SNAME}_gencodev45_Ensembl_canonical_TSS_FIRE_MINUS_Strand.bed"
    cat  $ALL_TSS | awk '$6 == "+"' > $TSS_PLUS
    cat  $ALL_TSS | awk '$6 == "-"' > $TSS_MINUS
    SINGLE_BED="../../msp_analysis/fibertools_msp/${SNAME}_consensus_BothStrands_HG38_corrected.haplotagged.MSP_BED6_single_150.bed.gz"
    PLUS_BED="tss_enrichment_intermediate/${SNAME}_DAF_tss_PLUS_2kb_MSPs_FIRE_peaks.bed"
    MINUS_BED="tss_enrichment_intermediate/${SNAME}_DAF_tss_MINUS_2kb_MSPs_FIRE_peaks.bed"
    PLUS_BG="tss_enrichment_intermediate/${SNAME}_DAF_tss_PLUS_2kb_MSPs_FIRE_peaks.bg"
    MINUS_BG="tss_enrichment_intermediate/${SNAME}_DAF_tss_MINUS_2kb_MSPs_FIRE_peaks.bg"
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED | sort-bed -) $TSS_PLUS > $PLUS_BED
    bedtools genomecov -bg -i $PLUS_BED -g $SIZES > $PLUS_BG
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED | sort-bed -) $TSS_MINUS > $MINUS_BED
    bedtools genomecov -bg -i $MINUS_BED -g $SIZES > $MINUS_BG
    MSP_BED="tss_enrichment_intermediate/${SNAME}_DAF_MSPs_FIRE_peaks_2kb.bed"
    MSP_BG="tss_enrichment_intermediate/${SNAME}_DAF_MSPs_FIRE_peaks_2kb.bg"
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED | sort-bed -) <(zcat $FIRE_PEAKS | sort-bed -) > $MSP_BED
    bedtools genomecov -bg -i $MSP_BED -g $SIZES > $MSP_BG
done


# scATAC-seq GM12878 --------------------------------------------------------
FRAGMENTS=../../data/fragments.tsv.gz
SINGLE_BED=../../data/passing_cell_fragments.tsv.gz
WHITELIST=../../data/10X_aggr_cell_barcodes.txt
zcat $FRAGMENTS | grep -f $WHITELIST | sort-bed - | bgzip -@ 20 > $SINGLE_BED
# TSS Enrichment Using ATAC-seq peaks
SNAME="scATAC_GM12878"
PEAKS="../../data/10X_GM12878_aggr_peaks.bed"
PEAKS_FILT="../../data/10X_GM12878_aggr_peaks_FILT.bed"
cat $PEAKS | grep -v KI | grep -v GL | grep -v "#" > $PEAKS_FILT
# only include TSSs that are in peaks (active TSSs)
ALL_TSS="../../data/${SNAME}_gencodev45_Ensembl_canonical_TSS_ATAC_Peaks.bed"
bedtools intersect -a $GENCODE -b $PEAKS_FILT -u | sort-bed - > $ALL_TSS
# split by TSS strand
TSS_PLUS="../../data/${SNAME}_gencodev45_Ensembl_canonical_TSS_PLUS_Strand.bed"
TSS_MINUS="../../data/${SNAME}_gencodev45_Ensembl_canonical_TSS_MINUS_Strand.bed"
cat  $ALL_TSS | awk '$6 == "+"' > $TSS_PLUS
cat  $ALL_TSS | awk '$6 == "-"' > $TSS_MINUS
PLUS_BED="../../data/${SNAME}_tss_PLUS_2kb_MSPs_ATAC_peaks.bed"
MINUS_BED="../../data/${SNAME}_tss_MINUS_2kb_MSPs_ATAC_peaks.bed"
PLUS_BG="../../data/${SNAME}_tss_PLUS_2kb_MSPs_ATAC_peaks.bg"
MINUS_BG="../../data/${SNAME}_tss_MINUS_2kb_MSPs_ATAC_peaks.bg"
bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) $TSS_PLUS > $PLUS_BED
bedtools genomecov -bg -i $PLUS_BED -g $SIZES > $PLUS_BG
bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED) $TSS_MINUS > $MINUS_BED
bedtools genomecov -bg -i $MINUS_BED -g $SIZES > $MINUS_BG


# ENCODE scATAC-seq (Pseudobulked) ------------------------------------
FRAGS=$(ls ../../data/scATAC_*_fragments.tsv.gz)
for SINGLE_BED in $FRAGS; do
    BNAME=$(basename $SINGLE_BED | cut -d'_' -f2)
    SNAME="scATAC_${BNAME}"
    PEAKS="../../data/${SNAME}_peaks.narrowPeak"
    PEAKS_FILT="../../data/${SNAME}_peaks_FILT.narrowPeak"
    # Call Peaks
    macs2 callpeak -f BED -t <(zcat $SINGLE_BED | sort-bed -) -g hs -n $SNAME --nomodel --shift -100 --extsize 200 --outdir ../../data
    cat $PEAKS | grep -v KI | grep -v GL | grep -v "#" > $PEAKS_FILT
done

for SINGLE_BED in $FRAGS; do
    BNAME=$(basename $SINGLE_BED | cut -d'_' -f2)
    SNAME="scATAC_${BNAME}"
    PEAKS_FILT="../../data/${SNAME}_peaks_FILT.narrowPeak"
    # only include TSSs that are in peaks (active TSSs)
    ALL_TSS="../../data/${SNAME}_gencodev45_Ensembl_canonical_TSS_ATAC_Peaks.bed"
    bedtools intersect -a $GENCODE -b $PEAKS_FILT -u | sort-bed - > $ALL_TSS
    # split by TSS strand
    TSS_PLUS="../../data/${SNAME}_gencodev45_Ensembl_canonical_TSS_PLUS_Strand.bed"
    TSS_MINUS="../../data/${SNAME}_gencodev45_Ensembl_canonical_TSS_MINUS_Strand.bed"
    cat  $ALL_TSS | awk '$6 == "+"' > $TSS_PLUS
    cat  $ALL_TSS | awk '$6 == "-"' > $TSS_MINUS
    PLUS_BED="../../data/${SNAME}_tss_PLUS_2kb_MSPs_ATAC_peaks.bed"
    MINUS_BED="../../data/${SNAME}_tss_MINUS_2kb_MSPs_ATAC_peaks.bed"
    PLUS_BG="../../data/${SNAME}_tss_PLUS_2kb_MSPs_ATAC_peaks.bg"
    MINUS_BG="../../data/${SNAME}_tss_MINUS_2kb_MSPs_ATAC_peaks.bg"
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED | sort-bed -) $TSS_PLUS > $PLUS_BED
    bedtools genomecov -bg -i $PLUS_BED -g $SIZES > $PLUS_BG
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED | sort-bed -) $TSS_MINUS > $MINUS_BED
    bedtools genomecov -bg -i $MINUS_BED -g $SIZES > $MINUS_BG
done


# OMNI-ATAC Bulk --------------------------
FQ=$(ls ../../data/*OMNI-ATAC_*.fastq.gz)
for F in $FQ; do
    OMNI_BAM=${F/.fastq.gz/.bam};
    DEDUP=${OMNI_BAM/.bam/.DEDUP.bam};
    METRICS=${F/.fastq.gz/.METRICS.txt};
    bowtie2 --very-sensitive --no-discordant --no-mixed --phred33 --threads 32 -X 1000 -x $REF -U $F | samtools sort - | samtools view -bS - > $OMNI_BAM;
    samtools index $OMNI_BAM;
    gatk MarkDuplicates --INPUT $OMNI_BAM --OUTPUT $DEDUP --METRICS_FILE $METRICS --REMOVE_DUPLICATES;
    samtools index $DEDUP;
    SINGLE_BED=${OMNI_BAM/.bam/.bed.gz};
    bedtools bamtobed -i $DEDUP | bgzip -@ 10 > $SINGLE_BED;
done

BAMS=$(ls ../../data/*OMNI-ATAC_*.DEDUP.bam)
for B in $BAMS; do
    # TSS Enrichment Using ATAC-seq peaks
    SNAME=$(basename $B | cut -d'_' -f2)
    PEAKS="../../data/OMNI-ATAC_${SNAME}_peaks.narrowPeak"
    PEAKS_FILT="../../data/OMNI-ATAC_${SNAME}_peaks_FILT.narrowPeak"
    SINGLE_BED=${B/.DEDUP.bam/.bed.gz}
    # Call Peaks
    macs2 callpeak -f BED -t <(zcat $SINGLE_BED | sort-bed -) -g hs -n "OMNI-ATAC_${SNAME}" --nomodel --shift -100 --extsize 200 --outdir ../../data/
    cat $PEAKS | grep -v KI | grep -v GL | grep -v "#" > $PEAKS_FILT
    # only include TSSs that are in peaks (active TSSs)
    ALL_TSS="../../data/OMNI-ATAC_${SNAME}_gencodev45_Ensembl_canonical_TSS_ATAC_Peaks.bed"
    bedtools intersect -a $GENCODE -b $PEAKS_FILT -u | sort-bed - > $ALL_TSS
    # split by TSS strand
    TSS_PLUS="../../data/OMNI-ATAC_${SNAME}_gencodev45_Ensembl_canonical_TSS_PLUS_Strand.bed"
    TSS_MINUS="../../data/OMNI-ATAC_${SNAME}_gencodev45_Ensembl_canonical_TSS_MINUS_Strand.bed"
    cat  $ALL_TSS | awk '$6 == "+"' > $TSS_PLUS
    cat  $ALL_TSS | awk '$6 == "-"' > $TSS_MINUS
    PLUS_BED="../../data/OMNI-ATAC_${SNAME}_tss_PLUS_2kb_MSPs_ATAC_peaks.bed"
    MINUS_BED="../../data/OMNI-ATAC_${SNAME}_tss_MINUS_2kb_MSPs_ATAC_peaks.bed"
    PLUS_BG="../../data/OMNI-ATAC_${SNAME}_tss_PLUS_2kb_MSPs_ATAC_peaks.bg"
    MINUS_BG="../../data/OMNI-ATAC_${SNAME}_tss_MINUS_2kb_MSPs_ATAC_peaks.bg"
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED | sort-bed -) $TSS_PLUS > $PLUS_BED
    bedtools genomecov -bg -i $PLUS_BED -g $SIZES > $PLUS_BG
    bedops --ec --range 2000 --element-of 1 <(zcat $SINGLE_BED | sort-bed -) $TSS_MINUS > $MINUS_BED
    bedtools genomecov -bg -i $MINUS_BED -g $SIZES > $MINUS_BG
done

