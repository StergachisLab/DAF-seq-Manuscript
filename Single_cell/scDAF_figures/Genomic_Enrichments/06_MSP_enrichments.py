import csv
from statistics import mean
import pandas as pd


window_len = 2000
flank_len = 100

header = ['sample','assay','position','count','enrichment']
out_lines = []

# TSSs ----------------------
FS_samples = ['HG002','GM12878','CHM13',"PS00272","PS00321","PS00327","PS00381","PS00382","PS00383","PS00384","ST001-liver","ST001-lung","PS00338_COLO829BL_1","PS00356_COLO829BL_2","COLO_T_2_PS00_418_451_488",'PS30743']
for sname in FS_samples:
    tss_pos = []
    with open(f"tss_enrichment_intermediate/{sname}_gencodev45_Ensembl_canonical_TSS_FIRE_PLUS_Strand.bed") as fr:
        reader = csv.reader(fr, delimiter="\t")
        for line in reader:
            tss_pos.append(':'.join([line[0],line[1]]))
    tss_neg = []
    with open(f"tss_enrichment_intermediate/{sname}_gencodev45_Ensembl_canonical_TSS_FIRE_MINUS_Strand.bed") as fr:
        reader = csv.reader(fr, delimiter="\t")
        for line in reader:
            tss_neg.append(':'.join([line[0],line[1]]))
    cov_dict_pos = {}
    with open(f"tss_enrichment_intermediate/{sname}_FS_tss_PLUS_2kb_MSPs_FIRE_peaks.bg") as fr:
        reader = csv.reader(fr, delimiter="\t")
        for line in reader:
            chrom = line[0]
            cov = int(line[3])
            for i in range(int(line[1]), int(line[2])):
                key = ':'.join([chrom,str(i)])
                cov_dict_pos[key] = cov
    cov_dict_neg = {}
    with open(f"tss_enrichment_intermediate/{sname}_FS_tss_MINUS_2kb_MSPs_FIRE_peaks.bg") as fr:
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
    flank_counts = [centered_counts_stranded[c] for c in range(0,flank_len)] + [centered_counts_stranded[c] for c in range(len(centered_counts_stranded)-flank_len,len(centered_counts_stranded))]
    norm = mean(flank_counts)
    scores = [centered_counts_stranded[i]/norm for i in range(len(centered_counts_stranded))]
    for i in range(len(centered_counts_stranded)):
        out_lines.append([sname, 'Fiber-seq', i-window_len, centered_counts_stranded[i], scores[i]])
    print(f'COMPLETED: {sname}')


DAF_samples = ['PS00718','PS00756','PS00757','PS00758','PS00867','PS00868','PS00869','PS00870','PS00871','PS00872','PS00873','PS00874']
for sname in DAF_samples:
    tss_pos = []
    with open(f"tss_enrichment_intermediate/{sname}_gencodev45_Ensembl_canonical_TSS_FIRE_PLUS_Strand.bed") as fr:
        reader = csv.reader(fr, delimiter="\t")
        for line in reader:
            tss_pos.append(':'.join([line[0],line[1]]))
    tss_neg = []
    with open(f"tss_enrichment_intermediate/{sname}_gencodev45_Ensembl_canonical_TSS_FIRE_MINUS_Strand.bed") as fr:
        reader = csv.reader(fr, delimiter="\t")
        for line in reader:
            tss_neg.append(':'.join([line[0],line[1]]))
    cov_dict_pos = {}
    with open(f"tss_enrichment_intermediate/{sname}_DAF_tss_PLUS_2kb_MSPs_FIRE_peaks.bg") as fr:
        reader = csv.reader(fr, delimiter="\t")
        for line in reader:
            chrom = line[0]
            cov = int(line[3])
            for i in range(int(line[1]), int(line[2])):
                key = ':'.join([chrom,str(i)])
                cov_dict_pos[key] = cov
    cov_dict_neg = {}
    with open(f"tss_enrichment_intermediate/{sname}_DAF_tss_MINUS_2kb_MSPs_FIRE_peaks.bg") as fr:
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
    flank_counts = [centered_counts_stranded[c] for c in range(0,flank_len)] + [centered_counts_stranded[c] for c in range(len(centered_counts_stranded)-flank_len,len(centered_counts_stranded))]
    norm = mean(flank_counts)
    scores = [centered_counts_stranded[i]/norm for i in range(len(centered_counts_stranded))]
    for i in range(len(centered_counts_stranded)):
        out_lines.append([sname, 'scDAF-seq', i-window_len, centered_counts_stranded[i], scores[i]])
    print(f'COMPLETED: {sname}')


ATAC_samples = ['scATAC_GM12878','scATAC_K562','scATAC_Liver','scATAC_Colon','scATAC_Lung','scATAC_Thyroid','scATAC_CD8+','OMNI-ATAC_K562','OMNI-ATAC_GM12878-repA']
for sname in ATAC_samples:
    tss_pos = []
    with open(f"ATAC_REV2/{sname}_gencodev45_Ensembl_canonical_TSS_PLUS_Strand.bed") as fr:
        reader = csv.reader(fr, delimiter="\t")
        for line in reader:
            tss_pos.append(':'.join([line[0],line[1]]))
    tss_neg = []
    with open(f"ATAC_REV2/{sname}_gencodev45_Ensembl_canonical_TSS_MINUS_Strand.bed") as fr:
        reader = csv.reader(fr, delimiter="\t")
        for line in reader:
            tss_neg.append(':'.join([line[0],line[1]]))
    cov_dict_pos = {}
    with open(f"ATAC_REV2/{sname}_tss_PLUS_2kb_MSPs_ATAC_peaks.bg") as fr:
        reader = csv.reader(fr, delimiter="\t")
        for line in reader:
            chrom = line[0]
            cov = int(line[3])
            for i in range(int(line[1]), int(line[2])):
                key = ':'.join([chrom,str(i)])
                cov_dict_pos[key] = cov
    cov_dict_neg = {}
    with open(f"ATAC_REV2/{sname}_tss_MINUS_2kb_MSPs_ATAC_peaks.bg") as fr:
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
    flank_counts = [centered_counts_stranded[c] for c in range(0,flank_len)] + [centered_counts_stranded[c] for c in range(len(centered_counts_stranded)-flank_len,len(centered_counts_stranded))]
    norm = mean(flank_counts)
    scores = [centered_counts_stranded[i]/norm for i in range(len(centered_counts_stranded))]
    for i in range(len(centered_counts_stranded)):
        out_lines.append([sname, 'ATAC-seq', i-window_len, centered_counts_stranded[i], scores[i]])
    print(f'COMPLETED: {sname}')

with open("ALL_2kb_tss_counts_FIRE_peaks_Stranded.csv",'w') as fw:
    writer = csv.writer(fw)
    writer.writerow(header)
    for l in out_lines:
        writer.writerow(l)

