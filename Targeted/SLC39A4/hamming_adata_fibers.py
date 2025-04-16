import pandas as pd
import numpy as np

import scanpy as sc
import anndata as ad

adata = sc.read('AnnData_GA_unfiltered_Liver_GM12878.h5ad')
with open('adata_zmw_all.txt','w') as fw:
    for sz in adata.obs.index:
        fw.write(f'{sz}\n')

sub_zmw_file = '/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/clustering/adata_zmw_all.txt'
sub_zmw = set()
with open(sub_zmw_file) as fr:
    for line in fr:
        sub_zmw.add(line.strip())

df_GA = df_GA[df_GA.index.isin(sub_zmw)]


# Filter out exact matches
dedup_dir='/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/duplicate_reads/exact_matches'
dz_files = ['Liver_SLC39A4_PS00680_map-pb_corrected_realigned_SLC39A4_region_DAF_DeDuplicated_zmws.txt', 'GM12878_SLC39A4_PS00686_map-pb_corrected_realigned_SLC39A4_region_DAF_DeDuplicated_zmws.txt']
dedup_zmws = set()
for dzf in dz_files:
    with open(f'{dedup_dir}/{dzf}') as fr:
        for line in fr:
            dedup_zmws.add(line.strip())

# GA reads
df_GA = pd.read_csv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/clustering/tissue_DA_features_GA.csv.gz', index_col='zmw')
df_GA = df_GA[df_GA.index.isin(dedup_zmws)]

df_h1_l = df_GA[(df_GA['hap'] == 1) & (df_GA['tissue'] == 'Liver')]
df_h2_l = df_GA[(df_GA['hap'] == 2) & (df_GA['tissue'] == 'Liver')]
df_h1_g = df_GA[(df_GA['hap'] == 1) & (df_GA['tissue'] == 'GM12878')]
df_h2_g = df_GA[(df_GA['hap'] == 2) & (df_GA['tissue'] == 'GM12878')]

# sample 25,000 zmws from each sample / Hap
n_sample = 25000
df_h1_l = df_h1_l.sample(n_sample)
df_h2_l = df_h2_l.sample(n_sample)
df_h1_g = df_h1_g.sample(n_sample)
df_h2_g = df_h2_g.sample(n_sample)


def distance_matrix(A, B):
    """
    Compute all pairwise distances between vectors in A and B.

    Parameters
    ----------
    A : np.array
        shape should be (M, K)
    B : np.array
        shape should be (N, K)

    Returns
    -------
    D : np.array
        A matrix D of shape (M, N).  Each entry in D i,j represnets the
        distance between row i in A and row j in B.

    See also
    --------
    A more generalized version of the distance matrix is available from
    scipy (https://www.scipy.org) using scipy.spatial.distance_matrix,
    which also gives a choice for p-norm.
    """
    M = A.shape[0]
    N = B.shape[0]
    assert (A.shape[1] == B.shape[1]), f"{A.shape[1]} does not match that of B {B.shape[1]}!"
    A_dots = (A * A).sum(axis=1).reshape((M, 1)) * np.ones(shape=(1, N))
    B_dots = (B * B).sum(axis=1) * np.ones(shape=(M, 1))
    D_squared = A_dots + B_dots - 2 * A.dot(B.T)
    return(D_squared)

def par_dist_matrix(a1, a2, keep, adj_set, chunk_n):
    dist = distance_matrix(a1, a2)
    # add to keep to an adjacency list
    kx,ky = np.where(dist <= keep)
    for idx in zip(kx,ky+i):
        adj_set.add(tuple(sorted(idx)))
    print(chunk_n)
    return(True)

max_dist = 3

# samples = {'Liver_H1':df_h1_l, 'Liver_H2':df_h2_l, 'GM12878_H1':df_h1_g, 'GM12878_H2':df_h2_g}

# samples = {'Liver_H1':df_h1_l}
# samples = {'Liver_H2':df_h2_l}
# samples = {'GM12878_H1':df_h1_g}
samples = {'GM12878_H2':df_h2_g}


for s_name,sam_df in samples.items():
    print(f'Starting {s_name}....')
    array_sub = sam_df.drop(['tissue','hap'], axis=1).to_numpy()
    dist = distance_matrix(array_sub, array_sub)
    kx,ky = np.where(dist <= max_dist)
    adj_set = {tuple(sorted(idx)) for idx in zip(kx,ky) if idx[0] != idx[1]}
    with open(f'{s_name}_adj_set.tsv','w') as fw:
        for adj in adj_set:
            fw.write(f'{adj[0]}\t{adj[1]}\n')
    # group fibers with low hamming dist
    groups = {i:set() for i in range(len(array_sub))}
    for f1,f2 in adj_set:
        groups[f1].add(f2)
        groups[f2].add(f1)
    # count single match
    counted = set()
    solo = set()
    for k,v in groups.items():
        if len(v) == 0:
            solo.add(k)
            counted.add(k)
        else:
            if len(v.intersection(counted)) == 0: # none of the matches have been tracked yet
                solo.add(list(v)[0])
            counted = counted.union(v)
            counted.add(k)
    solo_zmw = sam_df.iloc[list(solo)].index
    with open(f'{s_name}_unique_adata_zmw.txt','w') as fw:
        for z in solo_zmw:
            fw.write(f'{z}\n')
    print(f'Completed {s_name}')
    print(f'{len(solo)} solo reads')
    print(f'{len(solo) / len(counted)} prop')
    print()


dedup_dir = '/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/duplicate_reads'
dup_files = ['Liver_H1_unique_adata_zmw.txt', 'Liver_H2_unique_adata_zmw.txt', 'GM12878_H1_unique_adata_zmw.txt', 'GM12878_H2_unique_adata_zmw.txt',]
keep_zmw = []
for f in dup_files:
    file_zmw = set()
    with open(f'{dedup_dir}/{f}') as fr:
        for line in fr:
            file_zmw.add(int(line.strip()))
    keep_zmw.append(file_zmw)

LH1 = df_h1_l.iloc[list(keep_zmw[0])].index
LH2 = df_h2_l.iloc[list(keep_zmw[1])].index
GMH1 = df_h1_g.iloc[list(keep_zmw[2])].index
GMH2 = df_h2_g.iloc[list(keep_zmw[3])].index

all_zmw = []
for zmw in LH1:
    all_zmw.append(zmw)
for zmw in LH2:
    all_zmw.append(zmw)
for zmw in GMH1:
    all_zmw.append(zmw)
for zmw in GMH2:
    all_zmw.append(zmw)

with open('keep_zmw_all.txt','w') as fw:
    for z in all_zmw:
        fw.write(f'{z}\n')

