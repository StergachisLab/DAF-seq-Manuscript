import csv

# Filter NAPA deamination proportions for bases within the promoter region but between footprint regions:
    # Deamination rate at defined non-TF element regions within NAPA promoter 
    # Deamination rates surrounding element 1, 4, 7, and 11 at the various concentrations 

# promoter coordinates
with open('../footprinting/napa_promoter_region.bed') as fr:
    chrom, p_start, p_end = next(fr).strip().split('\t')
    p_start = int(p_start)
    p_end = int(p_end)


# deamination props
props = []
with open('../qc_plots/prop_da_titration_NAPA.tsv') as fr:
    reader = csv.DictReader(fr, delimiter="\t")
    for row in reader:
        if int(row['position']) > p_start and int(row['position']) <= p_end:
            props.append(row)


# footprint regions
regions = dict()
with open('../footprinting/merged_ft_on_both_strands.bed') as fr:
    reader = csv.reader(fr, delimiter="\t")
    i = 1
    for line in reader:
        regions[i] = {'start':int(line[1]), 'end':int(line[2])}
        i += 1
# add fake region to handle positions following region 11
regions[12] = {'start':p_end, 'end':p_end+1000}


# extract props between footprint regions
non_reg_props = []
curr_reg = 1
last_filt_reg = 2 # begin after region 2
last_include_region = 11
for p in props:
    pos = int(p['position'])
    if pos > regions[curr_reg]['end'] and pos <= regions[curr_reg+1]['start'] and curr_reg >= last_filt_reg and curr_reg < last_include_region:
        data = list(p.values())
        data.append(curr_reg)
        data.append(curr_reg+1)
        non_reg_props.append(data)
    if pos >= regions[curr_reg+1]['start'] and pos < p_end:
        curr_reg += 1

# extract props within footprint 11 CTCF module 2 (same as used for footprinting) 'GTCGGCC' chr19:47,515,590-47,515,596
ctcf_props = []
for p in props:
    pos = int(p['position'])
    if pos >= 47515590 and pos <= 47515596:
        ctcf_props.append(list(p.values()))


# classify positions as 'TC' or 'non-TC'
ft_pos = [int(pos[0]) for pos in non_reg_props]
tc_calls = {pos:'' for pos in ft_pos}

with open('../footprinting/NAPA/hg38_NAPA_promoter.fa') as fr:
    header = next(fr).strip()
    seq = next(fr).strip().upper()

start_coord = 47514458+1

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


# classify positions as 'CpG' or 'non-CpG'
cpg_calls = {pos:'' for pos in ft_pos}

for i in range(1, len(seq)-1):
    coord = start_coord+i
    if coord in ft_pos:
        if seq[i] == 'C':
            if seq[i+1] == 'G':
                cpg_calls[coord] = 'CpG'
            else:
                cpg_calls[coord] = 'non-CpG'
        elif seq[i] == 'G':
            if seq[i-1] == 'C':
                cpg_calls[coord] = 'CpG'
            else:
                cpg_calls[coord] = 'non-CpG'


# write out data
header = ['position', 'PS00626', 'PS00627', 'PS00628', 'PS00629', 'PS00630', 'PS00631','left_region','right_region','TC','CpG']
with open('non_ft_da_props_NAPA.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    writer.writerow(header)
    for p in non_reg_props:
        ct_call = tc_calls[int(p[0])]
        p.append(ct_call)
        cpg_call = cpg_calls[int(p[0])]
        p.append(cpg_call)
        writer.writerow(p)


header = ['position', 'PS00626', 'PS00627', 'PS00628', 'PS00629', 'PS00630', 'PS00631']
with open('ctcf_mod2_da_props_NAPA.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    writer.writerow(header)
    for p in ctcf_props:
        writer.writerow(p)

