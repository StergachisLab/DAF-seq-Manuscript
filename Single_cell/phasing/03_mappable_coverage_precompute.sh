#!/bin/bash
set -euo pipefail

MERGED_HIGH_COV="haplotagged/High_Cov_regions_AllCells.bed.gz"
MERGED_HIGH_COV_ST="haplotagged/High_STRAND_Cov_regions_AllCells.bed.gz"
UNMAP="unreliable_blacklist_highCov_regions_merged.bed"
UNMAP_ST="unreliable_blacklist_highCov_STRAND_ONLY_regions_merged.bed"


# Merged set of unmappable positions to ignore in coverage calculations
zcat unreliable-coverage-regions.bed.gz blacklist_hg38_ENCODE_2020_ENCFF356LFX.bed.gz $MERGED_HIGH_COV_ST | sort-bed - | bedops --ec --merge - > $UNMAP_ST

# Merged set of unmappable positions to ignore in FIRE analyses
zcat unreliable-coverage-regions.bed.gz blacklist_hg38_ENCODE_2020_ENCFF356LFX.bed.gz $MERGED_HIGH_COV | sort-bed - | bedops --ec --merge - > $UNMAP


# Compute coverage against the mappable autosomal hg38 genome
AUTO_BED="hg38_autosomes.sizes.bed"
cat hg38_autosomes.sizes | awk '{OFS="\t"}{ print $1,0,$2 }' | sort-bed - > $AUTO_BED
echo "Mappable hg38 regions" > mappable_hg38_precompute.txt
bedops --ec --difference hg38_autosomes.sizes.bed $UNMAP_ST | awk -F "\t" '{ print $3-$2 }' | datamash sum 1 >> mappable_hg38_precompute.txt

samples=('PS00718' 'PS00756' 'PS00757' 'PS00758' 'PS00867' 'PS00868' 'PS00869' 'PS00870' 'PS00871' 'PS00872' 'PS00873' 'PS00874')
ST=('CT' 'GA')
HAP=('1' '2')


COV_BED="mappable_coverage_precompute.tsv"
if [[ -e $COV_BED ]]; then rm $COV_BED; fi

# all haps
for s in "${samples[@]}"; do
OUT_BED="${s}_ALL_mappable_regions.bed"
if [[ -e $OUT_BED ]]; then rm $OUT_BED; fi
    bedtools bamtobed -i "haplotagged/${s}_consensus_BothStrands_HG38_corrected.haplotagged.bam" | sort-bed - | bedtools intersect -a stdin -b $AUTO_BED | sort-bed - | bedops --ec --merge - >> $OUT_BED;
    TOT_BP=$(cat $OUT_BED | sort-bed - | bedops --ec --merge - | bedops --ec --difference - $UNMAP_ST | awk -F "\t" '{ print $3-$2 }' | datamash sum 1 );
    echo $s "ALL" $TOT_BP >> $COV_BED;
done

# H1 & H2
for s in "${samples[@]}"; do
    for h in "${HAP[@]}"; do
        OUT_BED="${s}_H${h}_mappable_regions.bed"
        if [[ -e $OUT_BED ]]; then rm $OUT_BED; fi
        samtools view -b -d "HP:${h}" "haplotagged/${s}_consensus_BothStrands_HG38_corrected.haplotagged.bam" | bedtools bamtobed -i stdin | sort-bed - | bedtools intersect -a stdin -b $AUTO_BED | sort-bed - | bedops --ec --merge - >> $OUT_BED;
        TOT_BP=$(cat $OUT_BED | sort-bed - | bedops --ec --merge - | bedops --ec --difference - $UNMAP_ST | awk -F "\t" '{ print $3-$2 }' | datamash sum 1 );
        echo $s "H${h}" $TOT_BP >> $COV_BED;
    done
done


# Filter uncollapsed and collapsed PS00758 chr2 reads for unmappable positions for KaryoplotR

# Collapsed reads
for BAM in $(ls haplotagged/*_consensus_BothStrands_HG38_corrected.haplotagged.bam); do
    outfile=$(basename "haplotagged/${BAM/_consensus_BothStrands_HG38_corrected.haplotagged.bam/_no_high_STRAND_coverage_subtract.bed.gz}")
    bedtools bamtobed -i $bam | sort-bed - | bedtools subtract -a stdin -b $UNMAP_ST | sort-bed - | bgzip -@ 20 > $outfile;
done


# Raw clipped read
cat /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/collapse/clipped_seqs/PS00756_clipped_CT.tsv | awk '{OFS="\t"}{ print $2,$3,$4,$1 }' | sort-bed - | bgzip -@ 10 > PS00756_clipped_CT.bed.gz
cat /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/collapse/clipped_seqs/PS00756_clipped_GA.tsv | awk '{OFS="\t"}{ print $2,$3,$4,$1 }' | sort-bed - | bgzip -@ 10 > PS00756_clipped_GA.bed.gz
zcat PS00756_clipped_CT.bed.gz PS00758_clipped_GA.bed.gz | sort-bed - | bedtools subtract -a stdin -b <(zcat $UNMAP_GZ $MERGED_HIGH_COV_ST | sort-bed -) | sort-bed - | bgzip -@ 20 > PS00756_clipped_FILTERED.bed.gz
