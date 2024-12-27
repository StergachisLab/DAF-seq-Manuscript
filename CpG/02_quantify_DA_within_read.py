""" Quantify CpG vs non-CpG cytosine deamination within individual reads. """

import pysam
import csv

NAPA_DddA_bam = 'PS00720_NAPA_CpG_DddA_m84055_240822_221245_s3.hifi_reads.bc2057.map-pb_NAPA_Ref_corrected_realigned.bam'
UBA1_DddA_bam = 'PS00721_UBA1_CpG_DddA_m84055_240822_221245_s3.hifi_reads.bc2058.map-pb_UBA1_Ref_corrected_realigned.bam'

# reads must spand the entire region minus this buffer
reg_buffer = 100
da_cutoff = 0.1
primer_buff = 28 # ignore priming sites

# Process NAPA
chrom = 'chr19:47514458-47519061'
start = 30
end = 4527

# read in reference C & CpG positions
with open('NAPA_hg38.fa') as fr:
    header = next(fr)
    ref_seq = next(fr).strip().upper()
NAPA_cpg_pos = []
NAPA_C_pos = []
for i in range(start-1, end):
    if ref_seq[i] == 'C':
        if ref_seq[i+1] == 'G':
            NAPA_cpg_pos.append(i)
        else:
            NAPA_C_pos.append(i)
    elif ref_seq[i] == 'G':
        if ref_seq[i-1] == 'C':
            NAPA_cpg_pos.append(i)
        else:
            NAPA_C_pos.append(i)

# Quantify deamination at each position -------------------------------------------------------------
NAPA_zmws = set()
NAPA_pos_counts = dict()
bam = pysam.AlignmentFile(NAPA_DddA_bam, "rb")
for read in bam.fetch(chrom, start, end):
    if read.is_secondary == False and read.is_supplementary == False:
        if (read.reference_start <= start+reg_buffer) and (read.reference_end >= end-reg_buffer): # read spans the entire region
            seq = read.query_sequence
            try:
                strand = read.get_tag('ST')
            except:
                strand = None
            if strand == 'CT':
                nY = seq.count('Y')
                nC = seq.count('C')
                if (nY / (nY+nC)) > da_cutoff:
                    NAPA_zmws.add(read.qname)
                    pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
                    for pos in pair:
                        qi = pos[0]
                        ref_pos = pos[1]
                        if qi == None: # deletion, ignore
                            pass
                        elif ref_pos == None: # insertion, use seq base
                            pass
                        else:
                            # count the basecalls
                            ref_base = pos[2].upper()
                            if ref_base == 'C' and ref_pos >= (start+primer_buff) and ref_pos <= (end-primer_buff):
                                if ref_pos not in NAPA_pos_counts.keys():
                                    NAPA_pos_counts[ref_pos] = {'C':0, 'Y':0}
                                if seq[qi] != ref_pos:
                                    if seq[qi] == 'Y':
                                        NAPA_pos_counts[ref_pos]['Y'] += 1
                                    elif seq[qi] == 'C':
                                        NAPA_pos_counts[ref_pos]['C'] += 1
            elif strand == 'GA':
                nR = seq.count('R')
                nG = seq.count('G')
                if (nR / (nR+nG)) > da_cutoff:
                    NAPA_zmws.add(read.qname)
                    pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
                    for pos in pair:
                        qi = pos[0]
                        ref_pos = pos[1]
                        if qi == None: # deletion, ignore
                            pass
                        elif ref_pos == None: # insertion, use seq base
                            pass
                        else:
                            # count the basecalls
                            ref_base = pos[2].upper()
                            if ref_base == 'G' and ref_pos >= (start+primer_buff) and ref_pos <= (end-primer_buff):
                                if ref_pos not in NAPA_pos_counts.keys():
                                    NAPA_pos_counts[ref_pos] = {'G':0, 'R':0}
                                if seq[qi] != ref_pos:
                                    if seq[qi] == 'R':
                                        NAPA_pos_counts[ref_pos]['R'] += 1
                                    elif seq[qi] == 'G':
                                        NAPA_pos_counts[ref_pos]['G'] += 1

# aggregate % DA at each class of C base
NAPA_CpG_DA_tot = [0, 0] # Ref, deaminated
NAPA_nonCpG_DA_tot = [0, 0] # Ref, deaminated
for k,v in NAPA_pos_counts.items():
    # CpG
    if k in NAPA_cpg_pos:
        if 'C' in v.keys():
            NAPA_CpG_DA_tot[0] += v['C']
            NAPA_CpG_DA_tot[1] += v['Y']
        elif 'G' in v.keys():
            NAPA_CpG_DA_tot[0] += v['G']
            NAPA_CpG_DA_tot[1] += v['R']
    elif k in NAPA_C_pos:
        if 'C' in v.keys():
            NAPA_nonCpG_DA_tot[0] += v['C']
            NAPA_nonCpG_DA_tot[1] += v['Y']
        elif 'G' in v.keys():
            NAPA_nonCpG_DA_tot[0] += v['G']
            NAPA_nonCpG_DA_tot[1] += v['R']


# individual positions for violin
NAPA_CpG_DA_sepBP = []
NAPA_nonCpG_DA_sepBP = []
for k,v in NAPA_pos_counts.items():
    # CpG
    if k in NAPA_cpg_pos:
        if 'C' in v.keys():
            NAPA_CpG_DA_sepBP.append(v['Y']/sum(v.values()))
        elif 'G' in v.keys():
            NAPA_CpG_DA_sepBP.append(v['R']/sum(v.values()))
    elif k in NAPA_C_pos:
        if 'C' in v.keys():
            NAPA_nonCpG_DA_sepBP.append(v['Y']/sum(v.values()))
        elif 'G' in v.keys():
            NAPA_nonCpG_DA_sepBP.append(v['R']/sum(v.values()))



# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------


# Process UBA1
chrom = 'chrX:47190561-47194939'
start = 1
end = 4378

# read in reference C & CpG positions
with open('UBA1_hg38.fa') as fr:
    header = next(fr)
    ref_seq = next(fr).strip().upper()
UBA1_cpg_pos = []
UBA1_C_pos = []
for i in range(start-1, end):
    if ref_seq[i] == 'C':
        if ref_seq[i+1] == 'G':
            UBA1_cpg_pos.append(i)
        else:
            UBA1_C_pos.append(i)
    elif ref_seq[i] == 'G':
        if ref_seq[i-1] == 'C':
            UBA1_cpg_pos.append(i)
        else:
            UBA1_C_pos.append(i)

# Quantify deamination at each position -------------------------------------------------------------
UBA1_zmws = set()
UBA1_pos_counts = dict()
bam = pysam.AlignmentFile(UBA1_DddA_bam, "rb")
for read in bam.fetch(chrom, start, end):
    if read.is_secondary == False and read.is_supplementary == False:
        if (read.reference_start <= start+reg_buffer) and (read.reference_end >= end-reg_buffer): # read spans the entire region
            seq = read.query_sequence
            try:
                strand = read.get_tag('ST')
            except:
                strand = None
            if strand == 'CT':
                nY = seq.count('Y')
                nC = seq.count('C')
                if (nY / (nY+nC)) > da_cutoff:
                    UBA1_zmws.add(read.qname)
                    pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
                    for pos in pair:
                        qi = pos[0]
                        ref_pos = pos[1]
                        if qi == None: # deletion, ignore
                            pass
                        elif ref_pos == None: # insertion, use seq base
                            pass
                        else:
                            # count the basecalls
                            ref_base = pos[2].upper()
                            if ref_base == 'C' and ref_pos >= (start+primer_buff) and ref_pos <= (end-primer_buff):
                                if ref_pos not in UBA1_pos_counts.keys():
                                    UBA1_pos_counts[ref_pos] = {'C':0, 'Y':0}
                                if seq[qi] != ref_pos:
                                    if seq[qi] == 'Y':
                                        UBA1_pos_counts[ref_pos]['Y'] += 1
                                    elif seq[qi] == 'C':
                                        UBA1_pos_counts[ref_pos]['C'] += 1
            elif strand == 'GA':
                nR = seq.count('R')
                nG = seq.count('G')
                if (nR / (nR+nG)) > da_cutoff:
                    UBA1_zmws.add(read.qname)
                    pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
                    for pos in pair:
                        qi = pos[0]
                        ref_pos = pos[1]
                        if qi == None: # deletion, ignore
                            pass
                        elif ref_pos == None: # insertion, use seq base
                            pass
                        else:
                            # count the basecalls
                            ref_base = pos[2].upper()
                            if ref_base == 'G' and ref_pos >= (start+primer_buff) and ref_pos <= (end-primer_buff):
                                if ref_pos not in UBA1_pos_counts.keys():
                                    UBA1_pos_counts[ref_pos] = {'G':0, 'R':0}
                                if seq[qi] != ref_pos:
                                    if seq[qi] == 'R':
                                        UBA1_pos_counts[ref_pos]['R'] += 1
                                    elif seq[qi] == 'G':
                                        UBA1_pos_counts[ref_pos]['G'] += 1

# aggregate % DA at each class of C base
UBA1_CpG_DA_tot = [0, 0] # Ref, deaminated
UBA1_nonCpG_DA_tot = [0, 0] # Ref, deaminated
for k,v in UBA1_pos_counts.items():
    # CpG
    if k in UBA1_cpg_pos:
        if 'C' in v.keys():
            UBA1_CpG_DA_tot[0] += v['C']
            UBA1_CpG_DA_tot[1] += v['Y']
        elif 'G' in v.keys():
            UBA1_CpG_DA_tot[0] += v['G']
            UBA1_CpG_DA_tot[1] += v['R']
    elif k in UBA1_C_pos:
        if 'C' in v.keys():
            UBA1_nonCpG_DA_tot[0] += v['C']
            UBA1_nonCpG_DA_tot[1] += v['Y']
        elif 'G' in v.keys():
            UBA1_nonCpG_DA_tot[0] += v['G']
            UBA1_nonCpG_DA_tot[1] += v['R']


# individual positions for violin
UBA1_CpG_DA_sepBP = []
UBA1_nonCpG_DA_sepBP = []
for k,v in UBA1_pos_counts.items():
    # CpG
    if k in UBA1_cpg_pos:
        if 'C' in v.keys():
            UBA1_CpG_DA_sepBP.append(v['Y']/sum(v.values()))
        elif 'G' in v.keys():
            UBA1_CpG_DA_sepBP.append(v['R']/sum(v.values()))
    elif k in UBA1_C_pos:
        if 'C' in v.keys():
            UBA1_nonCpG_DA_sepBP.append(v['Y']/sum(v.values()))
        elif 'G' in v.keys():
            UBA1_nonCpG_DA_sepBP.append(v['R']/sum(v.values()))


# output summary dataframe
header = ['Region','group','numPos','ref','deaminated','prop_da']
out_rows = []
out_rows.append(['NAPA','CpG', len(NAPA_cpg_pos), NAPA_CpG_DA_tot[0], NAPA_CpG_DA_tot[1], NAPA_CpG_DA_tot[1]/sum(NAPA_CpG_DA_tot)])
out_rows.append(['NAPA','nonCpG', len(NAPA_C_pos), NAPA_nonCpG_DA_tot[0], NAPA_nonCpG_DA_tot[1], NAPA_nonCpG_DA_tot[1]/sum(NAPA_nonCpG_DA_tot)])
out_rows.append(['UBA1','CpG', len(UBA1_cpg_pos), UBA1_CpG_DA_tot[0], UBA1_CpG_DA_tot[1], UBA1_CpG_DA_tot[1]/sum(UBA1_CpG_DA_tot)])
out_rows.append(['UBA1','nonCpG', len(UBA1_C_pos), UBA1_nonCpG_DA_tot[0], UBA1_nonCpG_DA_tot[1], UBA1_nonCpG_DA_tot[1]/sum(UBA1_nonCpG_DA_tot)])

with open('deamination_summary.csv','w') as fout:
    writer = csv.writer(fout)
    writer.writerow(header)
    for r in out_rows:
        writer.writerow(r)


# write ZMWs used in analysis for BAM filtering
with open('zmws_used_NAPA.txt','w')as fout:
    for z in NAPA_zmws:
        fout.write(f'{z}\n')
with open('zmws_used_UBA1.txt','w')as fout:
    for z in UBA1_zmws:
        fout.write(f'{z}\n')

