#!/bin/bash
set -euo pipefail


# Phased "contigs", unphased "contigs", and total "contigs"

FYAK=/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/yak_switching/HG003.father.yak
MYAK=/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/yak_switching/HG004.mother.yak

# BAMS=$(ls ../haplotagged/PS007*.bam)

BAMS=$(ls ../haplotagged/PS00758-4*.bam)

for FILE in $BAMS; do
    B=$(basename $FILE)
    HAP_LIST=${FILE/haplotagged.bam/haplotag_list.txt.gz}
    HL_BASE=$(basename $HAP_LIST)
    PHASED=${HL_BASE/txt.gz/phased_reads.txt}
    UN=${HL_BASE/txt.gz/UNphased_reads.txt}
    UN_FA=${B/bam/UNPHASED.fa}
    PHASED_FA=${B/bam/PHASED.fa}
    ALL_FA=${B/bam/ALL.fa}
    sed 1d <(zcat $HAP_LIST) | awk -F "\t" '{ if($2 != "none") print $1 }' > $PHASED
    sed 1d <(zcat $HAP_LIST) | awk -F "\t" '{ if($2 == "none") print $1 }' > $UN
    samtools view -b -N $UN $FILE | samtools fasta > $UN_FA
    samtools view -b -N $PHASED $FILE | samtools fasta > $PHASED_FA
    samtools view -b $FILE | samtools fasta > $ALL_FA
    yak trioeval -t32 $FYAK $MYAK $ALL_FA > ${ALL_FA/.fa/.trioeval.tbl} &
    yak trioeval -t32 $FYAK $MYAK $PHASED_FA > ${PHASED_FA/.fa/.trioeval.tbl} &
    yak trioeval -t32 $FYAK $MYAK $UN_FA > ${UN_FA/.fa/.trioeval.tbl} &
done

# # Run on the clipped Reads ------------------------------------------------------

# # converted clipped sequences from TSV to FASTA using: clipped_read_tsv_to_fasta.py

# FASTA="/mmfs1/gscratch/stergachislab/assemblies/simple-names/hg38.fa"

# # Re-align the clipped fastas
# CLIP_FA=$(ls clipped_seqs/*clipped*.fa)
# for F in $CLIP_FA; do 
#     minimap2 -t 40 --MD -Y -y -a -x map-pb $FASTA $F | samtools sort > ${F/fa/map-pb_HG38.bam};
#     samtools index ${F/fa/map-pb_HG38.bam};
# done

# # Haplotag the clipped reads
# VCF="/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/analysis-HG002_WGS-4142-deepvariant.phased.vcf.gz"
# BAMS=$(realpath clipped_seqs/*clipped*HG38.bam)
# for B in $BAMS; do
#     whatshap haplotag -o ${B/bam/haplotagged.bam} --reference $FASTA --ignore-read-groups --output-haplotag-list ${B/bam/haplotag_list.txt.gz} $VCF $B;
#     samtools index ${B/bam/haplotagged.bam};
# done

# BAMS=$(ls clipped_seqs/PS007*.haplotagged.bam)
# for B in $BAMS; do
#     HAP_LIST=${B/haplotagged.bam/haplotag_list.txt.gz}
#     PHASED=${HAP_LIST/txt.gz/phased_reads.txt}
#     UN=${HAP_LIST/txt.gz/UNphased_reads.txt}
#     UN_FA=${B/bam/UNPHASED.fa}
#     PHASED_FA=${B/bam/PHASED.fa}
#     ALL_FA=${B/bam/ALL.fa}
#     sed 1d <(zcat $HAP_LIST) | awk -F "\t" '{ if($2 != "none") print $1 }' > $PHASED
#     sed 1d <(zcat $HAP_LIST) | awk -F "\t" '{ if($2 == "none") print $1 }' > $UN
#     samtools view -b -N $UN $B | samtools fasta > $UN_FA
#     samtools view -b -N $PHASED $B | samtools fasta > $PHASED_FA
#     samtools view -b $B | samtools fasta > $ALL_FA
#     yak trioeval -t32 $FYAK $MYAK $ALL_FA > ${ALL_FA/.fa/.trioeval.tbl} &
#     yak trioeval -t32 $FYAK $MYAK $PHASED_FA > ${PHASED_FA/.fa/.trioeval.tbl} &
#     yak trioeval -t32 $FYAK $MYAK $UN_FA > ${UN_FA/.fa/.trioeval.tbl} &
# done


# # Run Yak on consensus reads >= 10 Kb -----------------------------------------------------------------------------

# # filter BAM files for reads >= 10 Kb
# BAMS=$(ls ../haplotagged/PS007*.bam)
# for FILE in $BAMS; do
#     echo $FILE
#     samtools view $FILE | awk '{ if(length($10)>=10000) print $1 }' >> "length_cutoff/tenKb_consensus_reads.txt"
# done

# for FILE in $BAMS; do
#     B=$(basename $FILE)
#     samtools view -b -N length_cutoff/tenKb_consensus_reads.txt $FILE > "length_cutoff/${B/.bam/_10Kb.bam}";
#     samtools index "length_cutoff/${B/.bam/_10Kb.bam}";
# done

# # Haplotag the clipped reads
# VCF="/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/analysis-HG002_WGS-4142-deepvariant.phased.vcf.gz"
# BAMS=$(ls length_cutoff/PS007*.bam)
# for FILE in $BAMS; do
#     B=$(basename $FILE)
#     whatshap haplotag -o "length_cutoff/${B/bam/haplotagged.bam}" --reference $FASTA --ignore-read-groups --output-haplotag-list "length_cutoff/${B/bam/haplotag_list.txt.gz}" $VCF $FILE;
#     samtools index "length_cutoff/${B/bam/haplotagged.bam}";
# done

# # Phased "contigs", unphased "contigs", and total "contigs"
# FYAK=/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/yak_switching/HG003.father.yak
# MYAK=/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/yak_switching/HG004.mother.yak

# BAMS=$(ls length_cutoff/*haplotagged.bam)

# for FILE in $BAMS; do
#     HAP_LIST=${FILE/.haplotagged.bam/.haplotag_list.txt.gz}
#     PHASED=${HAP_LIST/txt.gz/phased_reads.txt}
#     UN=${HAP_LIST/txt.gz/UNphased_reads.txt}
#     UN_FA=${FILE/bam/UNPHASED.fa}
#     PHASED_FA=${FILE/bam/PHASED.fa}
#     ALL_FA=${FILE/bam/ALL.fa}
#     sed 1d <(zcat $HAP_LIST) | awk -F "\t" '{ if($2 != "none") print $1 }' > $PHASED
#     sed 1d <(zcat $HAP_LIST) | awk -F "\t" '{ if($2 == "none") print $1 }' > $UN
#     samtools view -b -N $UN $FILE | samtools fasta > $UN_FA
#     samtools view -b -N $PHASED $FILE | samtools fasta > $PHASED_FA
#     samtools view -b $FILE | samtools fasta > $ALL_FA
#     yak trioeval -t32 $FYAK $MYAK $ALL_FA > ${ALL_FA/.fa/.trioeval.tbl} &
#     yak trioeval -t32 $FYAK $MYAK $PHASED_FA > ${PHASED_FA/.fa/.trioeval.tbl} &
#     yak trioeval -t32 $FYAK $MYAK $UN_FA > ${UN_FA/.fa/.trioeval.tbl} &
# done


