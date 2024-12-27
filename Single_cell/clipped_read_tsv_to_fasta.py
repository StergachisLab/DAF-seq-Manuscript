import csv
from glob import glob
import os

clipped_files = glob('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/collapse/clipped_seqs/*_clipped_*.tsv')

for f in clipped_files:
    out_fasta = f"clipped_seqs/{os.path.basename(f).replace('.tsv','.fa')}"
    with open(f) as clipped_f:
        reader = csv.reader(clipped_f, delimiter="\t")
        with open(out_fasta,'w') as fout:
            for row in reader:
                coord = f'{row[1]}:{row[2]}-{row[3]}'
                name = f'>{row[0]}'
                seq = row[4]
                fout.write(f'{name}\n')
                fout.write(f'{seq}\n')

