import csv
from glob import glob
import os

seq_dir = '/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/collapse/consensus_seqs'
fastas = glob(f'{seq_dir}/PS*.fa')

out_rows = []
for file in fastas:
    name = os.path.basename(file).split('_')[0]
    with open(file) as fr:
        for line in fr:
            if not line.startswith('>'):
                out_rows.append((len(line.strip()), name))

with open('consensus_read_lengths.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    for row in out_rows:
        writer.writerow(row)
