import pandas as pd
import csv
import gzip


col_names = ['chrom','start','end','repeat_name','repeat_length','strand','repeat_class','region_type','UNK','repeat_num']
rm_annot = pd.read_csv('hg38.fa.out.repeatmasker.sort.bed', sep="\t", names=col_names, header=None)


# 5) Read in intersected counts and aggregate by category ------------------------------------------------------------------------------------

# lookup dict of DF indices
keys = rm_annot['chrom'] + ':' + rm_annot['start'].astype(str) + ':' + rm_annot['end'].astype(str)
df_idx = {key:index for index,key in keys.iteritems()}

rm_annot['unmod'] = 0
rm_annot['deam'] = 0

map_file = "RepeatMasker_da_counts_by_pos_all_cells_BEDMAP.txt.gz"

unmod_list = [0 for i in range(len(rm_annot))]
deam_list = [0 for i in range(len(rm_annot))]

with gzip.open(map_file,'rt') as fr:
    for line in fr:
        element,overlaps = line.strip().split('|')
        element = element.split("\t")
        overlaps = overlaps.split(';')
        ec = element[0]
        es = int(element[1])
        ee = int(element[2])
        key = ':'.join([ec, element[1], element[2]])
        tot_unmod = 0
        tot_deam = 0
        for over in overlaps:
            unmod,deam = over.split("\t")[3:5]
            tot_unmod += int(unmod)
            tot_deam += int(deam)
        unmod_list[df_idx[key]] += tot_unmod
        deam_list[df_idx[key]] += tot_deam

rm_annot['unmod'] = unmod_list
rm_annot['deam'] = deam_list


sig_only = rm_annot[(rm_annot['unmod'] > 0) | (rm_annot['deam'] > 0)]


# all data ----------------------------------------
sig_only.to_csv("rm_deamination_counts.csv", index = False)


sig_only = pd.read_csv("rm_deamination_counts.csv")

# grouped by repeat class ----------------------------------------
sig_grouped = []
for g in sig_only.groupby('repeat_class'):
    sig_grouped.append([g[0], sum(g[1]['unmod']), sum(g[1]['deam'])])


# add non_RM counts
non_RM_counts = pd.read_csv("NON_RepeatMasker_da_counts_by_pos_all_cells.bed.gz", sep="\t", names=['chrom','start','end','unmod','deam'])
sig_grouped.append(['non_RM', non_RM_counts['unmod'].sum(), non_RM_counts['deam'].sum()])

# add FIRE peak counts
map_file = "FIRE_Peaks_da_counts_by_pos_all_cells_BEDMAP.txt.gz"

unmod_fire = 0
deam_fire = 0

with gzip.open(map_file,'rt') as fr:
    for line in fr:
        element,overlaps = line.strip().split('|')
        overlaps = overlaps.split(';')
        tot_unmod = 0
        tot_deam = 0
        for over in overlaps:
            unmod,deam = over.split("\t")[3:5]
            tot_unmod += int(unmod)
            tot_deam += int(deam)
        unmod_fire += tot_unmod
        deam_fire += tot_deam

sig_grouped.append(['FIRE_peaks', unmod_fire, deam_fire])


with open("rm_grouped_deamination_counts.csv",'w') as fout:
    writer = csv.writer(fout)
    writer.writerow(['repeat_class','unmod','deam'])
    for sg in sig_grouped:
        writer.writerow(sg)

