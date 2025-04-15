import csv
import statistics as st
import random
from glob import glob


files = glob('*_DAF_duplicate_reads.tsv')
depth_list=[10000, 25000, 100000, 250000, 1000000]

curr=1000
depth_list_high_range=[]
while curr <= 1200000:
    depth_list_high_range.append(curr)
    curr *= 2


# tsv_file:TSV duplicte read groups produced by "duplicate_ID.py", depth_list: number of ZMWs to downsample to.
def prop_dup_by_n_reads(tsv_file, depth_list, suffix):
    z_dict = {}
    with open(tsv_file) as fr:
        reader = csv.DictReader(fr, delimiter="\t")
        for row in reader:
            zmws = row['zmws'].split(';')
            for z in zmws:
                if z not in z_dict.keys():
                    z_dict[z] = set()
                not_z = {nz for nz in zmws if nz != z}
                z_dict[z] = z_dict[z].union(not_z)
        # compute downsampled dup counts
        out_stats = []
        for d in depth_list:
            if d <= len(z_dict):
                down_groups = []
                down_zmw = set(random.sample(sorted(z_dict.keys()), d))
                counted = set()
                for dz in down_zmw:
                    if dz not in counted:
                        counted.add(dz)
                        group = [dz]
                        for z in z_dict[dz]:
                            if z in down_zmw:
                                counted.add(z)
                                group.append(z)
                        down_groups.append(group)
                # get dup stats for downsampling
                dup = d - len(down_groups)
                n_copies = [len(g) for g in down_groups]
                out_stats.append([d, dup, (dup/d)*100, st.mean(n_copies), st.median(n_copies), max(n_copies)])
        # write deduplication stats
        out_stats_name = tsv_file.replace('duplicate_reads.tsv', suffix)
        with open(out_stats_name,'w') as fw:
            writer = csv.writer(fw, delimiter="\t")
            header = ['total_reads','duplicate_reads','percent_dup', 'mean_copies','median_copies','max_copies']
            writer.writerow(header)
            for row in out_stats:
                writer.writerow(row)


for f in files:
    prop_dup_by_n_reads(f, depth_list_high_range, 'DeDuplicated_stats_Downsample_highRange.tsv') # include 25k downsample

