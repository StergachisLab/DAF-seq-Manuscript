import csv


# splitting bivalent promoter into two FIRE peaks, with core of well-positioned nucleosome in between

# chr6:110180013-110180227 Peak 1 WASF1 promoter
# chr6:110180356-110180508 Peak 2 CD40 promoter


# promoter coordinates
p_start = int(110180013)
p_end = int(110180508)


# deamination props
props = []
with open('../qc_plots/prop_da_titration_WASF1.tsv') as fr:
    reader = csv.DictReader(fr, delimiter="\t")
    for row in reader:
        if int(row['position']) > p_start and int(row['position']) <= p_end:
            props.append(row)

# classify positions as 'TC' or 'non-TC'
ft_pos = [int(pos['position']) for pos in props]
tc_calls = {pos:'' for pos in ft_pos}

with open('../footprinting/WASF1/hg38_WASF1_promoter.fa') as fr:
    header = next(fr).strip()
    seq = next(fr).strip().upper()

start_coord = 110176758+1

for i in range(1, len(seq)-1):
    coord = start_coord+i
    if coord in ft_pos:
        if seq[i] == 'C':
            if seq[i-1] == 'T':
                tc_calls[coord] = 'TC'
            else:
                tc_calls[coord] = 'non-TC'
        elif seq[i] == 'G':
            if seq[i+1] == 'A':
                tc_calls[coord] = 'TC'
            else:
                tc_calls[coord] = 'non-TC'


# extract props within each region
peak1_props = []
for p in props:
    pos = int(p['position'])
    if pos >= 110180013 and pos <= 110180227:
        peak1_props.append(list(p.values()))

peak2_props = []
for p in props:
    pos = int(p['position'])
    if pos >= 110180356 and pos <= 110180508:
        peak2_props.append(list(p.values()))


# taking the center 20 bp of the well-positioned nucleosome
nuc_s = 110180228
nuc_e = 110180355
nuc_center = round((nuc_e-nuc_s)/2) + nuc_s + 1
nuc_buff = 10

nuc_props = []
for p in props:
    pos = int(p['position'])
    if pos >= (nuc_center-nuc_buff) and pos <= (nuc_center+nuc_buff):
        nuc_props.append(list(p.values()))


# write out data
header = ['position', 'PS00626', 'PS00627', 'PS00628', 'PS00629', 'PS00630', 'PS00631','CT']
with open('non_ft_da_props_WASF1_peak1.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    writer.writerow(header)
    for p in peak1_props:
        ct_call = tc_calls[int(p[0])]
        p.append(ct_call)
        writer.writerow(p)

header = ['position', 'PS00626', 'PS00627', 'PS00628', 'PS00629', 'PS00630', 'PS00631','CT']
with open('non_ft_da_props_WASF1_peak2.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    writer.writerow(header)
    for p in peak2_props:
        ct_call = tc_calls[int(p[0])]
        p.append(ct_call)
        writer.writerow(p)

header = ['position', 'PS00626', 'PS00627', 'PS00628', 'PS00629', 'PS00630', 'PS00631','CT']
with open('non_ft_da_props_WASF1_nuc.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    writer.writerow(header)
    for p in nuc_props:
        ct_call = tc_calls[int(p[0])]
        p.append(ct_call)
        writer.writerow(p)

