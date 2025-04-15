import pysam
import csv

bam = pysam.AlignmentFile('SMHTCOLO829BLT50_COLO_reg2_merged.bam', "rb")
vars = [19447245, 19447246]

# sets to remove duplicate cluster IDs (R1 and R2 from the same fragment both overlap sites)
ref = set()
v = set()
for read in bam.fetch('chr17',19447231,19447270):
    if read.is_secondary == False and read.is_supplementary == False:
        seq = read.seq
        pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
        p1 = None
        p2 = None
        for pos in pair:
            if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                mcoord = pos[0]
                ref_coord = pos[1] + 1
                if ref_coord == vars[0]:
                    p1 = seq[mcoord]
                elif ref_coord == vars[1]:
                    p2 = seq[mcoord]
        if p1 == 'C' and p2 == 'C':
            ref.add(read.qname)
        elif p1 == 'T' and p2 == 'T':
            v.add(read.qname)

nref = len(ref)
nv = len(v)

with open('COLO_var_stats.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['ref','var','prop'])
    writer.writerow([nref, nv, nv/(nref+nv)])
