import pysam
import csv


# Identify positions of X straight deamiantions within the NAPA promoter, output as BED file --------------------------------------------

napa_bam = "/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/Napa_WASF1/haplotype_correct_NAPA/NAPA_PS00626_haplotype_corrected.bam"

prom_start = 47515063
prom_end = 47515660

hap_basecalls = {47514700: ['T', 'C'],
 47515000: ['T', 'G'],
 47515603: ['G', 'G'],
 47516117: ['A', 'T'],
 47517160: ['C', 'T'],
 47517464: ['T', 'C'],
 47518032: ['T', 'T'],
 47518383: ['C', 'G'],
 47518977: ['C', 'C']}

run_length = 3

def find_runs(da_bools, run_len):
    indices = []
    for i in range(len(da_bools)+1 - run_len):
        if all(da_bools[i+j] == False for j in range(run_len)):
            indices.append(i)
    return(indices)

run_lines_all = []
run_lines_CT = []
run_lines_GA = []
bam = pysam.AlignmentFile(napa_bam, "rb")
for read in bam.fetch('chr19', prom_start, prom_end):
    seq = read.query_sequence
    pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
    strand = read.get_tag('ST')
    try:
        hap = read.get_tag('HP')
    except:
        hap = 'UNK'
    if hap in [1,2]:
        da_call_bool = []
        da_call_coord = []
        for pos in pair:
            if pos[0] == None or pos[1] == None: # indel, ignore
                pass
            else:
                da = None
                qi = pos[0]
                ref_coord = pos[1] + 1
                base = seq[qi]
                if ref_coord >= prom_start and ref_coord <= prom_end: # limit to promoter FIRE peak footprinting region
                    if strand == 'CT':
                        if ref_coord in hap_basecalls.keys():
                            if hap_basecalls[ref_coord][hap-1] in ['C']:
                                if base == 'Y':
                                    da = True
                                elif base == 'C':
                                    da = False
                        else:
                            if base == 'Y':
                                da = True
                            elif base == 'C':
                                da = False
                    elif strand == 'GA':
                        if ref_coord in hap_basecalls.keys():
                            if hap_basecalls[ref_coord][hap-1] in ['G']:
                                if base == 'R':
                                    da = True
                                elif base == 'G':
                                    da = False
                        else:
                            if base == 'R':
                                da = True
                            elif base == 'G':
                                da = False
                    if da != None:
                        da_call_bool.append(da)
                        da_call_coord.append(ref_coord-1)
        # find runs within a read's deamination calls
        run_idx = find_runs(da_call_bool, run_length)
        if strand == 'CT':
            for idx in run_idx:
                run_lines_all.append(['chr19', da_call_coord[idx], da_call_coord[idx+run_length-1]+1])
                run_lines_CT.append(['chr19', da_call_coord[idx], da_call_coord[idx+run_length-1]+1])
        elif strand == 'GA':
            for idx in run_idx:
                run_lines_all.append(['chr19', da_call_coord[idx], da_call_coord[idx+run_length-1]+1])
                run_lines_GA.append(['chr19', da_call_coord[idx], da_call_coord[idx+run_length-1]+1])


with open(f'run_of_{run_length}_da_occurences_all.bed','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    for line in run_lines_all:
        writer.writerow(line)

with open(f'run_of_{run_length}_da_occurences_TOP.bed','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    for line in run_lines_CT:
        writer.writerow(line)

with open(f'run_of_{run_length}_da_occurences_BOTTOM.bed','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    for line in run_lines_GA:
        writer.writerow(line)
