import pysam
import csv

out_lines = []
out_lines.append(['Region','DAF_reads','FS_reads','DAF_prop','FS_prop','DAF_enrichment'])


bam_name = '/mmfs1/gscratch/stergachislab/mvollger/projects/FIREv2.0/results/GM12878/fire/GM12878.fire.bam'
bam = pysam.AlignmentFile(bam_name, "rb")
FS_total_reads = 0
for read in bam.fetch():
    if read.is_secondary == False and read.is_supplementary == False:
        FS_total_reads += 1


# NAPA & WASF1 BAM
bam_name = '/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/Napa_WASF1/PS00626.m84046_240619_124816_s1.bc2072.ft.map-pb_corrected_realigned.bam'
bam = pysam.AlignmentFile(bam_name, "rb")
total_reads = 0
for read in bam.fetch():
    if read.is_secondary == False and read.is_supplementary == False:
        total_reads += 1

# chr19:47,514,458-47,519,061 NAPA
reg_chrom = 'chr19'
reg_start = 47514458
reg_end = 47519061
napa_reg_reads = 0
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        napa_reg_reads += 1

# chr6:110,176,758-110,181,579 WASF1
reg_chrom = 'chr6'
reg_start = 110176758
reg_end = 110181579
wasf1_reg_reads = 0
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        wasf1_reg_reads += 1

daf_target = napa_reg_reads+wasf1_reg_reads
prop_DAF = daf_target/total_reads


# GM12878 Fiber-seq BAM
bam_name = '/mmfs1/gscratch/stergachislab/mvollger/projects/FIREv2.0/results/GM12878/fire/GM12878.fire.bam'
bam = pysam.AlignmentFile(bam_name, "rb")

# chr19:47,514,458-47,519,061 NAPA
reg_chrom = 'chr19'
reg_start = 47514458
reg_end = 47519061
napa_reg_reads = 0
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        napa_reg_reads += 1

# chr6:110,176,758-110,181,579 WASF1
reg_chrom = 'chr6'
reg_start = 110176758
reg_end = 110181579
wasf1_reg_reads = 0
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        wasf1_reg_reads += 1
FS_target = napa_reg_reads+wasf1_reg_reads
prop_FS = FS_target/FS_total_reads

out_lines.append(['NAPA_WASF1', daf_target, FS_target, prop_DAF, prop_FS, prop_DAF/prop_FS])



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# UBA1
bam_name = 'UBA1/PS00685_GM12878_UBA1_m84046_240802_231812_s1.hifi_reads.bc2093.map-pb_corrected_realigned.bam'
bam = pysam.AlignmentFile(bam_name, "rb")
total_reads = 0
for read in bam.fetch():
    if read.is_secondary == False and read.is_supplementary == False:
        total_reads += 1

# chrX:47190561-47194939 UBA1
reg_chrom = 'chrX'
reg_start = 47190561
reg_end = 47194939
target_reg_reads = 0
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        target_reg_reads += 1
daf_target = target_reg_reads
prop_DAF = daf_target/total_reads


# GM12878 Fiber-seq BAM
bam_name = '/mmfs1/gscratch/stergachislab/mvollger/projects/FIREv2.0/results/GM12878/fire/GM12878.fire.bam'
bam = pysam.AlignmentFile(bam_name, "rb")
reg_chrom = 'chrX'
reg_start = 47190561
reg_end = 47194939
target_reg_reads = 0
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        target_reg_reads += 1
FS_target = target_reg_reads
prop_FS = FS_target/FS_total_reads

out_lines.append(['UBA1', daf_target, FS_target, prop_DAF, prop_FS, prop_DAF/prop_FS])


# SLC39A4
bam_name = 'SLC39A4/processed_bams/GM12878_SLC39A4_PS00686_m84046_240802_231812_s1.hifi_reads.bc2094.map-pb_corrected_realigned.bam'
bam = pysam.AlignmentFile(bam_name, "rb")
total_reads = 0
for read in bam.fetch():
    if read.is_secondary == False and read.is_supplementary == False:
        total_reads += 1

# chr8:144,415,793-144,417,939 SLC39A4
reg_chrom = 'chr8'
reg_start = 144415793
reg_end = 144417939
target_reg_reads = 0
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        target_reg_reads += 1
daf_target = target_reg_reads
prop_DAF = daf_target/total_reads


# GM12878 Fiber-seq BAM
bam_name = '/mmfs1/gscratch/stergachislab/mvollger/projects/FIREv2.0/results/GM12878/fire/GM12878.fire.bam'
bam = pysam.AlignmentFile(bam_name, "rb")
reg_chrom = 'chr8'
reg_start = 144415793
reg_end = 144417939
target_reg_reads = 0
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        target_reg_reads += 1
FS_target = target_reg_reads
prop_FS = FS_target/FS_total_reads

out_lines.append(['SLC39A4', daf_target, FS_target, prop_DAF, prop_FS, prop_DAF/prop_FS])


with open('coverage_enrighment_DAF_vs_FS.csv','w') as fout:
    writer = csv.writer(fout)
    for row in out_lines:
        writer.writerow(row)

