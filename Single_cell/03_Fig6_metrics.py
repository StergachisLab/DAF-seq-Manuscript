import pandas as pd
import numpy as np
from glob import glob


# average phasing rate from the four cells

haplotagged_dir = '/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/haplotagged'

tag_lists = glob(f'{haplotagged_dir}/*HG38_corrected.haplotag_list.txt.gz')

total = 0
phased = 0
for tl in tag_lists:
    df = pd.read_csv(tl, sep="\t")
    df_c = df['haplotype'].value_counts()
    total += sum(df_c)
    phased += df_c['H1'] + df_c['H2']

with open('avg_phased_read_prop.txt','w') as fout:
    fout.write(f'Avg Collapsed Read Phasing Proportion (all 4 cells): {phased/total}\n')
