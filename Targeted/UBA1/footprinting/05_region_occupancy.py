import csv


# Export footprinting as DataFrame
c = 1
regions = dict()
with open('merged_UBA1_GM12878_footprinted_motifs_H1_sorted.bed') as fr:
    reader = csv.reader(fr, delimiter="\t")
    for line in reader:
        regions[c] = line
        c += 1
df = pd.DataFrame(motif_dict).T

df_reg_numbers = []
df_regions = []
for i in range(len(df)):
    r = None
    st = df.iloc[i,:]['start']
    end = df.iloc[i,:]['end']
    max_overlap = 0
    max_reg = 'None'
    for n,reg in regions.items():
        if st >= int(reg[1]) and end <= int(reg[2]):
            overlap = end - int(reg[1])
            if overlap > max_overlap:
                max_overlap = overlap
                r = n
                max_reg = f'{reg[0]}:{reg[1]}-{reg[2]}'
    if r == None:
        df_reg_numbers.append(0)
        df_regions.append('None')
    else:
        df_reg_numbers.append(r)
        df_regions.append(max_reg)
df['reg_num'] = df_reg_numbers
df['region'] = df_regions

df.to_csv('footprint_stats.csv', index_label='motif')



# DF of zmws by region footprint status (not overlapped by MSP -> -1, MSP & unbound -> 0, MSP & bound -> 1)
reg_zmw_dict_H1 = {r:{'msp':set(),'bound':set()} for r in regions.keys()}
for i in range(len(df)):
    key = df.index[i]
    r = df.iloc[i,:]['reg_num']
    if key in ft_H1:
        reg_zmw_dict_H1[r]['msp'] = reg_zmw_dict_H1[r]['msp'].union(motif_zmws[key]['H1'][0])
        reg_zmw_dict_H1[r]['bound'] = reg_zmw_dict_H1[r]['bound'].union(motif_zmws[key]['H1'][1])

reg_zmw_dict_H2 = {r:{'msp':set(),'bound':set()} for r in regions.keys()}
for i in range(len(df)):
    key = df.index[i]
    r = df.iloc[i,:]['reg_num']
    if key in ft_H1:      
        reg_zmw_dict_H2[r]['msp'] = reg_zmw_dict_H2[r]['msp'].union(motif_zmws[key]['H2'][0])
        reg_zmw_dict_H2[r]['bound'] = reg_zmw_dict_H2[r]['bound'].union(motif_zmws[key]['H2'][1])

H1_msp_zmws = set()
for k,v in reg_zmw_dict_H1.items():
    H1_msp_zmws = H1_msp_zmws.union(reg_zmw_dict_H1[k]['msp'])
out_rows = []
for z in H1_msp_zmws:
    z_row = [z]
    for k in reg_zmw_dict_H1.keys():
        if z in reg_zmw_dict_H1[k]['msp']:
            if z in reg_zmw_dict_H1[k]['bound']:
                z_row.append(1)
            else:
                z_row.append(0)
        else:
            z_row.append(-1)
    out_rows.append(z_row)
z_df = pd.DataFrame(out_rows)
z_df = z_df.set_index([0])
z_df.to_csv('zmw_footprint_regions_UBA1_H1.csv', index_label='zmw')

# prop of fibers with footprint by region
reg_bed = []
for k,v in regions.items():
    row = []+v
    row.append(k)
    counts = Counter(z_df[k])
    row.append(counts[1]/(counts[0]+counts[1]))
    row.append('+')
    reg_bed.append(row)
with open('region_footprint_props_UBA1_H1.bed','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    for row in reg_bed:
        writer.writerow(row)

H2_msp_zmws = set()
for k,v in reg_zmw_dict_H2.items():
    H2_msp_zmws = H2_msp_zmws.union(reg_zmw_dict_H2[k]['msp'])
out_rows = []
for z in H2_msp_zmws:
    z_row = [z]
    for k in reg_zmw_dict_H2.keys():
        if z in reg_zmw_dict_H2[k]['msp']:
            if z in reg_zmw_dict_H2[k]['bound']:
                z_row.append(1)
            else:
                z_row.append(0)
        else:
            z_row.append(-1)
    out_rows.append(z_row)
z_df = pd.DataFrame(out_rows)
z_df = z_df.set_index([0])
z_df.to_csv('zmw_footprint_regions_UBA1_H2.csv', index_label='zmw')

# prop of fibers with footprint by region
reg_bed = []
for k,v in regions.items():
    row = []+v
    row.append(k)
    counts = Counter(z_df[k])
    row.append(counts[1]/(counts[0]+counts[1]))
    row.append('+')
    reg_bed.append(row)
with open('region_footprint_props_UBA1_H2.bed','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    for row in reg_bed:
        writer.writerow(row)
