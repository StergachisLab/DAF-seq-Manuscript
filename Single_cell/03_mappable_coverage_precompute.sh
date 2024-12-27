#!/bin/bash
set -euo pipefail

MERGED_HIGH_COV="haplotagged/High_Cov_regions_AllCells.bed.gz"
UNMAP="unreliable_blacklist_highCov_regions_merged.bed"

# Merged set of unmappable positions to ignore in coverage calculations, FIRE analyses
zcat unreliable-coverage-regions.bed.gz blacklist_hg38_ENCODE_2020_ENCFF356LFX.bed.gz $MERGED_HIGH_COV | sort-bed - | bedops --ec --merge - > $UNMAP

# Compute coverage against the mappable autosomal hg38 genome
AUTO_BED="hg38_autosomes.sizes.bed"
cat hg38_autosomes.sizes | awk '{OFS="\t"}{ print $1,0,$2 }' | sort-bed - > $AUTO_BED
echo "Mappable hg38 regions" > mappable_hg38_precompute.txt
bedops --ec --difference hg38_autosomes.sizes.bed $UNMAP | awk -F "\t" '{ print $3-$2 }' | datamash sum 1 >> mappable_hg38_precompute.txt

samples=('PS00718' 'PS00756' 'PS00757' 'PS00758')
ST=('CT' 'GA')
HAP=('1' '2')


COV_BED="mappable_coverage_precompute.tsv"
if [[ -e $COV_BED ]]; then rm $COV_BED; fi

# all haps
for s in "${samples[@]}"; do
OUT_BED="${s}_ALL_mappable_regions.bed"
if [[ -e $OUT_BED ]]; then rm $OUT_BED; fi
    for strand in "${ST[@]}"; do
        bedtools bamtobed -i "haplotagged/${s}_consensus_${strand}_HG38_corrected.haplotagged.bam" | sort-bed - | bedtools intersect -a stdin -b $AUTO_BED | sort-bed - | bedops --ec --merge - >> $OUT_BED;
    done
    TOT_BP=$(cat $OUT_BED | sort-bed - | bedops --ec --merge - | bedops --ec --difference - $UNMAP | awk -F "\t" '{ print $3-$2 }' | datamash sum 1 );
    echo $s "ALL" $TOT_BP >> $COV_BED;
done

# H1 & H2
for s in "${samples[@]}"; do
    for h in "${HAP[@]}"; do
        OUT_BED="${s}_H${h}_mappable_regions.bed"
        if [[ -e $OUT_BED ]]; then rm $OUT_BED; fi
        for strand in "${ST[@]}"; do
            samtools view -b -d "HP:${h}" "haplotagged/${s}_consensus_${strand}_HG38_corrected.haplotagged.bam" | bedtools bamtobed -i stdin | sort-bed - | bedtools intersect -a stdin -b $AUTO_BED | sort-bed - | bedops --ec --merge - >> $OUT_BED;
        done
        TOT_BP=$(cat $OUT_BED | sort-bed - | bedops --ec --merge - | bedops --ec --difference - $UNMAP | awk -F "\t" '{ print $3-$2 }' | datamash sum 1 );
        echo $s "H${h}" $TOT_BP >> $COV_BED;
    done
done


# Filter uncollapsed and collapsed PS00758 chr2 reads for unmappable positions for KaryoplotR

# # Collapsed consensus seqs
# for strand in "${ST[@]}"; do
#     IN_BAM="haplotagged/PS00758_consensus_${strand}_HG38_corrected.haplotagged.bam";
#     BAM_PRE=$(basename $IN_BAM);
#     samtools view -@ 10 -b -L $UNMAP -U "${BAM_PRE/bam/NoHighCov_filtered_chr2_${strand}.bam}" $IN_BAM chr2 > "${BAM_PRE/bam/HighCovOnly_chr2_${strand}.bam}";
#     samtools index "${BAM_PRE/bam/NoHighCov_filtered_chr2_${strand}.bam}";
#     samtools index "${BAM_PRE/bam/HighCovOnly_chr2_${strand}.bam}";
# done

# # Uncollapsed reads
# IN_BAM=/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/data/PS00758_DddA_F2_map-pb.HG38.merged.no-kinetics.bam
# for strand in "${ST[@]}"; do
#     BAM_PRE=$(basename $IN_BAM);
#     samtools view -@ 20 -d "ST:${strand}" -b -L $UNMAP -U "${BAM_PRE/bam/NoHighCov_filtered_chr2_${strand}.bam}" $IN_BAM chr2 > "${BAM_PRE/bam/HighCovOnly_chr2_${strand}.bam}";
#     samtools index "${BAM_PRE/bam/NoHighCov_filtered_chr2_${strand}.bam}";
#     samtools index "${BAM_PRE/bam/HighCovOnly_chr2_${strand}.bam}";
# done


# FILTER OUT COVERAGE of >3x per strand if the overlap is <= 425 bp (failure to collapse)
# pad regions by 500bp for high coverage region merging

UNMAP_GZ=${UNMAP/.bed/.bed.gz}
cat $UNMAP | bgzip -@ 10 > $UNMAP_GZ

for s in "${samples[@]}"; do
    for strand in "${ST[@]}"; do
        OUT_BED="${s}_${strand}_WIDE_high_cov_filt_regions.bed.gz"
        FILT_BED="${s}_${strand}_WIDE_high_cov_FILTERED_consensus_reads.bed.gz"
        if [[ -e $OUT_BED ]]; then rm $OUT_BED; fi
        if [[ -e $FILT_BED ]]; then rm $FILT_BED; fi
        IN_BAM="haplotagged/${s}_consensus_${strand}_HG38_corrected.haplotagged.bam";
        samtools depth $IN_BAM | awk '{OFS="\t"}{ if($3>2) print $1,$2-1,$2 }' | sort-bed - | bedops --ec --merge --range 500 - | bedops --ec --merge --range -500 - | awk '($3-$2) < 425' | sort-bed - | bgzip -@ 20 > $OUT_BED;
        bedtools bamtobed -i $IN_BAM | sort-bed - | bedtools subtract -a stdin -b <(zcat $UNMAP_GZ $OUT_BED | sort-bed -) | sort-bed - | bgzip -@ 20 > $FILT_BED;
    done
done


# Raw clipped read
cat /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/collapse/clipped_seqs/PS00758_clipped_CT.tsv | awk '{OFS="\t"}{ print $2,$3,$4,$1 }' | sort-bed - | bgzip -@ 10 > PS00758_clipped_CT.bed.gz &
cat /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/collapse/clipped_seqs/PS00758_clipped_GA.tsv | awk '{OFS="\t"}{ print $2,$3,$4,$1 }' | sort-bed - | bgzip -@ 10 > PS00758_clipped_GA.bed.gz &

# Uncollapsed reads
zcat PS00758_clipped_CT.bed.gz PS00758_clipped_GA.bed.gz | sort-bed - | bedtools subtract -a stdin -b <(zcat $UNMAP_GZ $MERGED_HIGH_COV | sort-bed -) | sort-bed - | bgzip -@ 20 > PS00758_clipped_FILTERED.bed.gz

