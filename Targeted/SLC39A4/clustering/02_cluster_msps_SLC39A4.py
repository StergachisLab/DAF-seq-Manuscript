import pandas as pd
import numpy as np
import itertools
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import csv

# Cluster GA reads only
df_DA_GA = pd.read_csv('tissue_DA_features_GA.csv.gz', low_memory=False, index_col='zmw')

# Filter DF to deduplicated ZMWs
dedup_dir='/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/duplicate_reads'
dz_files = ['GM12878_H1_unique_adata_zmw.txt','GM12878_H2_unique_adata_zmw.txt','Liver_H1_unique_adata_zmw.txt','Liver_H2_unique_adata_zmw.txt']
dedup_zmws = set()
for dzf in dz_files:
    with open(f'{dedup_dir}/{dzf}') as fr:
        for line in fr:
            dedup_zmws.add(line.strip())
df_DA_GA = df_DA_GA[df_DA_GA.index.isin(dedup_zmws)]


# Downsample by Hap (even fibers by sample & haplotype)
sample_n = 5000
df_H1 = df_DA_GA[df_DA_GA['hap'] == 1].groupby('tissue').sample(sample_n)
df_H2 = df_DA_GA[df_DA_GA['hap'] == 2].groupby('tissue').sample(sample_n)
df_s = pd.concat([df_H1, df_H2])
df_s_tissue = df_s['tissue']
df_s_hap = df_s['hap']
df_s = df_s.drop(['tissue','hap'], axis = 1)

# load AnnData
adata = sc.AnnData(df_s)
adata.obs["tissue"] = pd.Categorical(df_s_tissue) # add tissue as metadata
adata.obs["hap"] = pd.Categorical(df_s_hap) # add hap as metadata

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save='_DA_GA.pdf')
sc.pl.pca(adata, color=["tissue", "tissue", "hap", "hap"], dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2, size=2, save='_DA_GA.pdf')

# UMAP
sc.pp.neighbors(adata, n_pcs=0, n_neighbors=200)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["tissue","hap"],
    wspace=0.2, # increase horizontal space between panels
    save='_Tissue_Hap_DA_GA.pdf')

# Leiden Clustering
sc.tl.leiden(adata, flavor="igraph")
sc.pl.umap(adata, color=["leiden"], save='_Leiden_DA_GA.pdf')
cluster_counts = adata.obs['leiden'].value_counts()
cluster_counts.to_csv('cluster_counts_leiden_unfiltered_GA.csv')

# filter out fibers from clusters with < N ZMWs
N = 1000
save_clust = cluster_counts[cluster_counts > N].index
adata_2 = adata[adata.obs['leiden'].isin(save_clust)]

# change cluster 0 to 1
adata_2.obs['leiden'].replace('0','1', inplace=True)

cluster_counts = adata_2.obs['leiden'].value_counts()
cluster_counts.to_csv('cluster_counts_leiden_FINAL_GA.csv')


# Final UMAP with Leiden clusters
sc.pl.umap(adata_2, color=["leiden"], wspace=0.2, save='_FINAL_Leiden_DA_GA.pdf', palette="Set1")
# replace "#a65628" (brown) with "#4daf4a" (green)
adata_2.uns['leiden_colors'][4] = "#4daf4a"
sc.pl.umap(adata_2, color=["leiden"], wspace=0.2, save='_FINAL_Leiden_DA_GA.pdf')


# Quantify tissue contributions to each cluster
res="leiden"
header = ['cluster','num_Liver','num_GM12878']
cl_list = []
for cl in adata_2.obs[res].value_counts().index:
    cl_counts = adata_2[adata_2.obs[res] == cl].obs['tissue'].value_counts()
    try:
        liv_z = cl_counts['Liver']
    except:
        liv_z = 0
    try:
        gm_z = cl_counts['GM12878']
    except:
        gm_z = 0
    cl_list.append([cl, liv_z, gm_z])
with open('cluster_tissue_counts_GA_FINAL.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(header)
    for c in cl_list:
        writer.writerow(c)

# By Hap
for h in [1,2]:
    cl_list = []
    for cl in adata_2.obs[res].value_counts().index:
        cl_counts = adata_2[(adata_2.obs['hap'] == h) & (adata_2.obs[res] == cl)].obs['tissue'].value_counts()
        try:
            liv_z = cl_counts['Liver']
        except:
            liv_z = 0
        try:
            gm_z = cl_counts['GM12878']
        except:
            gm_z = 0
        cl_list.append([cl, liv_z, gm_z])
    header = ['cluster',f'Liver_H{h}',f'GM12878_H{h}']
    with open(f'cluster_tissue_counts_GA_FINAL_H{h}.csv','w') as fw:
        writer = csv.writer(fw)
        writer.writerow(header)
        for c in cl_list:
            writer.writerow(c)


# EXPORT CLUSTERS ------------------------------------------------------------------------------
cl_ids = sorted(adata_2.obs['leiden'].value_counts().index)
for c in cl_ids:
    z_list = list(adata_2[adata_2.obs['leiden'] == c].obs.index)
    with open(f'clust_pileups/cluster_{c}_DA_GA_zmw.txt','w') as fw:
        for z in z_list:
            fw.write(f'{z}\n')


# Recluster # --------------------------------------------

clust_6 = adata_2[adata_2.obs['leiden'] == '6']

sc.tl.pca(clust_6)
sc.pp.neighbors(clust_6, n_pcs=0, n_neighbors=50)
sc.tl.umap(clust_6)
sc.tl.leiden(clust_6, flavor="igraph")

sc.pl.umap(clust_6, color=["leiden"], save='_Clust6_Leiden_DA_GA.pdf')
cluster_counts = clust_6.obs['leiden'].value_counts()
cluster_counts.to_csv('cluster6_counts_leiden_unfiltered_GA.csv')

N = 200
save_clust = cluster_counts[cluster_counts > N].index
clust_6_2 = clust_6[clust_6.obs['leiden'].isin(save_clust)]

clust_6_2.obs['leiden'].replace('0','5', inplace=True)
clust_6_2.obs['leiden'].replace('8','6', inplace=True)

sc.pl.umap(clust_6_2, color=["leiden"], save='_Clust6_FINAL_Leiden_DA_GA.pdf')
cluster_counts = clust_6_2.obs['leiden'].value_counts()
cluster_counts.to_csv('cluster6_counts_leiden_FINAL_GA.csv')

# Sub Clust6 By Hap
for h in [1,2]:
    cl_list = []
    for cl in clust_6_2.obs['leiden'].value_counts().index:
        cl_counts = clust_6_2[(clust_6_2.obs['hap'] == h) & (clust_6_2.obs['leiden'] == cl)].obs['tissue'].value_counts()
        try:
            liv_z = cl_counts['Liver']
        except:
            liv_z = 0
        try:
            gm_z = cl_counts['GM12878']
        except:
            gm_z = 0
        cl_list.append([cl, liv_z, gm_z])
    header = ['cluster',f'Liver_H{h}',f'GM12878_H{h}']
    with open(f'cluster6_tissue_counts_GA_FINAL_H{h}.csv','w') as fw:
        writer = csv.writer(fw)
        writer.writerow(header)
        for c in cl_list:
            writer.writerow(c)



# EXPORT CLUSTERS ------------------------------------------------------------------------------
cl_ids = sorted(clust_6_2.obs['leiden'].value_counts().index)
for c in cl_ids:
    z_list = list(clust_6_2[clust_6_2.obs['leiden'] == c].obs.index)
    with open(f'clust_pileups/sub_Clust6_cluster_{c}_DA_GA_zmw.txt','w') as fw:
        for z in z_list:
            fw.write(f'{z}\n')

# Write AnnData objects to disk ----------------------------------------------------------
adata.write('AnnData_GA_unfiltered_Liver_GM12878.h5ad')
adata_2.write('AnnData_GA_FINAL_Liver_GM12878.h5ad')
clust_6.write('AnnData_GA_Cluster6_unfiltered_Liver_GM12878.h5ad')
clust_6_2.write('AnnData_GA_Cluster6_FINAL_Liver_GM12878.h5ad')


# Highlight Cluster 6 SubClusters on the FULL UMAP --------------------------------------------------------

adata_2 = sc.read('AnnData_GA_FINAL_Liver_GM12878.h5ad')
clust_6_2 = sc.read('AnnData_GA_Cluster6_FINAL_Liver_GM12878.h5ad')

sub_labels = []
sub_zmw = set(clust_6_2.obs.index)
all_zmw = list(adata_2.obs.index)
for i in range(len(all_zmw)):
    if all_zmw[i] in sub_zmw:
        sub_labels.append(clust_6_2.obs.loc[all_zmw[i]]['leiden'])
    else:
        sub_labels.append('-1')
        
adata_2.obs['sub_clust'] = sub_labels

sc.pl.umap(adata_2, color="sub_clust", save='GA_FINAL_SubCLust6_HIGHLIGHT.pdf')

