import pysam
import csv
import os


BAMS = ["/mmfs1/gscratch/stergachislab/FIRE/results/GM12878/GM12878-fire-v0.1.1-filtered.cram",
        "/mmfs1/gscratch/stergachislab/FIRE/results/CHM13/CHM13-fire-v0.1-filtered.cram",
        "/mmfs1/gscratch/stergachislab/mvollger/projects/FIREv2.0/results/HG002/fire/HG002.fire.cram",
        "/mmfs1/gscratch/stergachislab/userprod/production/fibertools-pipeline/PacBio-Fiber-seq/PS30743/hg38/results/pipelines/fire-run/results/PS30743/PS30743-fire-v0.1.2-filtered.cram"]

FS_SAMPLES=["PS00272","PS00321","PS00327","PS00381","PS00382","PS00383","PS00384","ST001-liver","ST001-lung","PS00338_COLO829BL_1","PS00356_COLO829BL_2","COLO_T_2_PS00_418_451_488"]
for sname in FS_SAMPLES:
    BAMS.append(f"/mmfs1/gscratch/stergachislab/mvollger/projects/FIREv2.0/results/{sname}/fire/{sname}.fire.cram")


# NAPA
reg_chrom = 'chr19'
reg_start = 47515063
reg_end = 47515660

# NAPA promoter footprint regions for excluding bases
ft_bam = "/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/Napa_WASF1/footprinting/NAPA/merged_both_strands.bed"
ft_pos = set()
with open(ft_bam) as fr:
    for line in fr:
        s = line.strip().split('\t')[1]
        e = line.strip().split('\t')[2]
        for i in range(int(s),int(e)+1):
            ft_pos.add(i)

with open('napa_footprint_positions.txt','w') as fw:
    for p in sorted(ft_pos):
        fw.write(f'{p}\n')


bmods = [('A', 1, 'a'), ('T', 0, 'a'), ('A', 0, 'a'), ('T', 1, 'a')]
out_lines=[['sample','m6a_prop']]
pos_out_lines = [['ref_coord','sample','m6a_prop']]

for bam_name in BAMS:
    if "consensus" in bam_name:
        sname = os.path.basename(bam_name).split('_')[0]
    elif "filtered" in bam_name:
        sname = os.path.basename(bam_name).split('-')[0]
    else:
        sname = os.path.basename(bam_name).split('.')[0]
    bam = pysam.AlignmentFile(bam_name, "rb")
    read_props = []
    pos_m6a = dict()
    for read in bam.fetch(reg_chrom, reg_start, reg_end):
        if read.is_secondary == False and read.is_supplementary == False:
            ref_AT = set() # most precise per-read m6a prop to account for indel bases
            seq = read.seq
            pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
            mod_pos = set()
            for bm in bmods:
                if bm in read.modified_bases.keys():
                    mod_pos = mod_pos.union({p[0] for p in read.modified_bases[bm]})
            read_m6a = 0
            for pos in pair:
                mcoord = pos[0]
                ref_coord = pos[1]
                if ref_coord != None:
                    ref_coord += 1
                if mcoord != None and pos[1] != None and ref_coord >= reg_start and ref_coord <= reg_end and ref_coord not in ft_pos:
                    if pos[2].upper() in ['A','T']:
                        if ref_coord not in pos_m6a.keys():
                            pos_m6a[ref_coord] = [0,0] # cov, m6a
                        pos_m6a[ref_coord][0] += 1
                        ref_AT.add(ref_coord)
                        if mcoord in mod_pos:
                            pos_m6a[ref_coord][1] += 1
                            read_m6a += 1
            read_props.append(read_m6a/len(ref_AT))
    for prop in read_props:
        out_lines.append([sname, prop])
    for k,v in pos_m6a.items():
        pos_out_lines.append([k, sname, v[1]/v[0]])
        
with open('saturation_Fiber-seq_NAPA_m6a_read_props.csv','w') as fw:
    writer = csv.writer(fw)
    for line in out_lines:
        writer.writerow(line)

with open('saturation_Fiber-seq_NAPA_m6a_ref_coord_props.csv','w') as fw:
    writer = csv.writer(fw)
    for line in pos_out_lines:
        writer.writerow(line)

