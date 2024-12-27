from glob import glob
import csv
import os

collapsed_yak = glob('*HG38_corrected.haplotagged.ALL.trioeval.tbl')
clipped_yak = glob('clipped_seqs/*HG38.haplotagged.ALL.trioeval.tbl')

out_rows = []
out_rows.append(['sample','strand','read_type','yak_type','rate'])

for file in sorted(collapsed_yak):
    name = os.path.basename(file).split('_')[0]
    if 'CT' in os.path.basename(file).split('_')[2]:
        strand = 'CT'
    elif 'GA' in os.path.basename(file).split('_')[2]:
        strand = 'GA'
    with open(file) as fr:
        reader = csv.reader(fr, delimiter="\t")
        for row in reader:
            if row[0] == 'W':
                out_rows.append([name,strand,'collapsed',row[0],row[3]])
            elif row[0] == 'H':
                out_rows.append([name,strand,'collapsed',row[0],row[3]])

for file in sorted(clipped_yak):
    name = os.path.basename(file).split('_')[0]
    if 'CT' in os.path.basename(file).split('_')[2]:
        strand = 'CT'
    elif 'GA' in os.path.basename(file).split('_')[2]:
        strand = 'GA'
    with open(file) as fr:
        reader = csv.reader(fr, delimiter="\t")
        for row in reader:
            if row[0] == 'W':
                out_rows.append([name,strand,'clipped',row[0],row[3]])
            elif row[0] == 'H':
                out_rows.append([name,strand,'clipped',row[0],row[3]])

with open('yak_error_rate_summary.csv','w') as fout:
    writer = csv.writer(fout)
    for row in out_rows:
        writer.writerow(row)

