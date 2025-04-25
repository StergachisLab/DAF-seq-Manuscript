import csv
import gzip
from math import floor
import pandas as pd


window_len = 2000


# TSSs within HG002 FIRE peaks split by strand (+/-)
tss_pos = []
with open('gencodev45_Ensembl_canonical_TSS_HG002_FIRE_noUnreliableCoverage_PLUS_Strand.bed') as fr:
    reader = csv.reader(fr, delimiter="\t")
    for line in reader:
        tss_pos.append(':'.join([line[0],line[1]]))
tss_neg = []
with open('gencodev45_Ensembl_canonical_TSS_HG002_FIRE_noUnreliableCoverage_MINUS_Strand.bed') as fr:
    reader = csv.reader(fr, delimiter="\t")
    for line in reader:
        tss_neg.append(':'.join([line[0],line[1]]))

cov_dict_pos = {}
with open('DAF_tss_PLUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bg') as fr:
    reader = csv.reader(fr, delimiter="\t")
    for line in reader:
        chrom = line[0]
        cov = int(line[3])
        for i in range(int(line[1]), int(line[2])):
            key = ':'.join([chrom,str(i)])
            cov_dict_pos[key] = cov
cov_dict_neg = {}
with open('DAF_tss_MINUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bg') as fr:
    reader = csv.reader(fr, delimiter="\t")
    for line in reader:
        chrom = line[0]
        cov = int(line[3])
        for i in range(int(line[1]), int(line[2])):
            key = ':'.join([chrom,str(i)])
            cov_dict_neg[key] = cov

centered_counts_stranded = [0 for i in range((window_len*2)+1)]

for t in tss_pos:
    chrom,coord = t.split(':')
    coord = int(coord)
    for i in range(coord-window_len, coord+window_len+1):
        key = ':'.join([chrom,str(i)])
        if key in cov_dict_pos.keys():
            centered_counts_stranded[i-coord+window_len] += cov_dict_pos[key]
for t in tss_neg:
    chrom,coord = t.split(':')
    coord = int(coord)
    for i in range(coord-window_len, coord+window_len+1):
        key = ':'.join([chrom,str(i)])
        if key in cov_dict_neg.keys():
            centered_counts_stranded[-(i-coord)+window_len] += cov_dict_neg[key]

with open('DAF_2kb_tss_counts_FIRE_peaks_Stranded.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['position','count'])
    for i in range(len(centered_counts_stranded)):
        writer.writerow([i-window_len, centered_counts_stranded[i]])


# Fiber-seq -------------------------------------------------------------------
cov_dict_pos = {}
with open('FS_tss_PLUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bg') as fr:
    reader = csv.reader(fr, delimiter="\t")
    for line in reader:
        chrom = line[0]
        cov = int(line[3])
        for i in range(int(line[1]), int(line[2])):
            key = ':'.join([chrom,str(i)])
            cov_dict_pos[key] = cov
cov_dict_neg = {}
with open('FS_tss_MINUS_2kb_MSPs_FIRE_peaks_all_cells_and_strands.bg') as fr:
    reader = csv.reader(fr, delimiter="\t")
    for line in reader:
        chrom = line[0]
        cov = int(line[3])
        for i in range(int(line[1]), int(line[2])):
            key = ':'.join([chrom,str(i)])
            cov_dict_neg[key] = cov

centered_counts_stranded = [0 for i in range((window_len*2)+1)]
for t in tss_pos:
    chrom,coord = t.split(':')
    coord = int(coord)
    for i in range(coord-window_len, coord+window_len+1):
        key = ':'.join([chrom,str(i)])
        if key in cov_dict_pos.keys():
            centered_counts_stranded[i-coord+window_len] += cov_dict_pos[key]
for t in tss_neg:
    chrom,coord = t.split(':')
    coord = int(coord)
    for i in range(coord-window_len, coord+window_len+1):
        key = ':'.join([chrom,str(i)])
        if key in cov_dict_neg.keys():
            centered_counts_stranded[-(i-coord)+window_len] += cov_dict_neg[key]

with open('FS_2kb_tss_counts_FIRE_peaks_Stranded.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['position','count'])
    for i in range(len(centered_counts_stranded)):
        writer.writerow([i-window_len, centered_counts_stranded[i]])



# --------------------------------------------------------------------------------------------------------------

# Enrichment over all FIRE peaks (not Just TSSs)
fire_peaks = "../msp_analysis/FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz"

fire_pos = dict()
with gzip.open(fire_peaks, 'rt') as fr:
    reader = csv.reader(fr, delimiter="\t")
    for line in reader:
        center = int(line[1])+round((int(line[2])-int(line[1])-1)/2)
        fire_pos[':'.join([line[0],str(center)])] = float(line[4])

cov_dict_FIRE = {}
with open('DAF_MSPs_FIRE_peaks_2kb_all_cells_and_strands.bg') as fr:
    reader = csv.reader(fr, delimiter="\t")
    for line in reader:
        chrom = line[0]
        cov = int(line[3])
        for i in range(int(line[1]), int(line[2])):
            key = ':'.join([chrom,str(i)])
            cov_dict_FIRE[key] = cov

centered_counts_FIRE = [[0 for i in range((window_len*2)+1)] for i in range(1,10)]

for pos,prop in fire_pos.items():
    bin = floor((prop*100)//10) # 10-20% actuation bin is bin 1
    if bin == 10: # 100% with 90-99%
        bin = 9
    chrom,coord = pos.split(':')
    coord = int(coord)
    for i in range(coord-window_len, coord+window_len+1):
        key = ':'.join([chrom,str(i)])
        if key in cov_dict_FIRE.keys():
            centered_counts_FIRE[bin-1][i-coord+window_len] += cov_dict_FIRE[key]

with open('DAF_MSPs_2kb_counts_binned_FIRE_peaks.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['position','count','bin'])
    for i in range((window_len*2)+1):
        for b in range(1,10):
            writer.writerow([i-window_len, centered_counts_FIRE[b-1][i], b])


# Fiber-seq -------------------------------------------------------------------
cov_dict_FIRE = {}
with open('FS_MSPs_FIRE_peaks_2kb_all_cells_and_strands.bg') as fr:
    reader = csv.reader(fr, delimiter="\t")
    for line in reader:
        chrom = line[0]
        cov = int(line[3])
        for i in range(int(line[1]), int(line[2])):
            key = ':'.join([chrom,str(i)])
            cov_dict_FIRE[key] = cov

centered_counts_FIRE = [[0 for i in range((window_len*2)+1)] for i in range(1,10)]

for pos,prop in fire_pos.items():
    bin = floor((prop*100)//10) # 10-20% actuation bin is bin 1
    if bin == 10: # 100% with 90-99%
        bin = 9
    chrom,coord = pos.split(':')
    coord = int(coord)
    for i in range(coord-window_len, coord+window_len+1):
        key = ':'.join([chrom,str(i)])
        if key in cov_dict_FIRE.keys():
            centered_counts_FIRE[bin-1][i-coord+window_len] += cov_dict_FIRE[key]

with open('FS_MSPs_2kb_counts_binned_FIRE_peaks.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['position','count','bin'])
    for i in range((window_len*2)+1):
        for b in range(1,10):
            writer.writerow([i-window_len, centered_counts_FIRE[b-1][i], b])


# RM enrichment of MSPs ------------------------------------------------------------------------------
samples = ['DAF','FS']

for samp in samples:
    col_names = ['chrom','start','end','repeat_name','repeat_length','strand','repeat_class','region_type','UNK','repeat_num']
    rm_annot = pd.read_csv('hg38.fa.out.repeatmasker.MAPPABLE.sort.bed', sep="\t", names=col_names, header=None)
    keys = rm_annot['chrom'] + ':' + rm_annot['start'].astype(str) + ':' + rm_annot['end'].astype(str)
    df_idx = {key:index for index,key in keys.items()}
    rm_annot['rm_bases'] = rm_annot['end'] - rm_annot['start']
    map_file = f"{samp}_RepeatMasker_MSP_BEDMAP.txt.gz"
    msp_list = [0 for i in range(len(rm_annot))]
    with gzip.open(map_file,'rt') as fr:
        for line in fr:
            element,overlaps = line.strip().split('|')
            element = element.split("\t")
            overlaps = overlaps.split(';')
            ec = element[0]
            key = ':'.join([ec, element[1], element[2]])
            msp_list[df_idx[key]] += len(overlaps)
    rm_annot['n_msp'] = msp_list
    rm_annot_grouped = []
    for g in rm_annot.groupby('repeat_class'):
        rm_annot_grouped.append([g[0], sum(g[1]['rm_bases']), sum(g[1]['n_msp'])])
    # add FIRE peak counts
    map_file = f"{samp}_FIRE_Peaks_MSP_BEDMAP.txt.gz"
    bases_fire = 0
    n_msp_fire = 0
    with gzip.open(map_file,'rt') as fr:
        for line in fr:
            element,overlaps = line.strip().split('|')
            element = element.split("\t")
            overlaps = overlaps.split(';')
            bases_fire += int(element[2]) - int(element[1])
            n_msp_fire += len(overlaps)
    rm_annot_grouped.append(['FIRE_peaks', bases_fire, n_msp_fire])
    # add non_RM counts
    map_file = f"{samp}_NON_RepeatMasker_MSP_BEDMAP.txt.gz"
    bases_nonRM = 0
    n_msp_nonRM = 0
    with gzip.open(map_file,'rt') as fr:
        for line in fr:
            element,overlaps = line.strip().split('|')
            element = element.split("\t")
            overlaps = overlaps.split(';')
            bases_nonRM += int(element[2]) - int(element[1])
            n_msp_nonRM += len(overlaps)
    rm_annot_grouped.append(['non_RM', bases_nonRM, n_msp_nonRM])
    # All regions
    map_file = f"{samp}_ALL_MSP_BEDMAP.txt.gz"
    bases_ALL = 0
    n_msp_ALL = 0
    with gzip.open(map_file,'rt') as fr:
        for line in fr:
            element,overlaps = line.strip().split('|')
            element = element.split("\t")
            overlaps = overlaps.split(';')
            bases_ALL += int(element[2]) - int(element[1])
            n_msp_ALL += len(overlaps)
    rm_annot_grouped.append(['ALL', bases_ALL, n_msp_ALL])
    with open(f"rm_grouped_MSP_counts_{samp}.csv",'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['repeat_class','rm_bases','n_msp'])
        for sg in rm_annot_grouped:
            writer.writerow(sg)

