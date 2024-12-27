import pandas as pd


fire_df = pd.read_csv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/scDAF_FIRE_actuation_MSP150.tsv', sep="\t")

df1 = fire_df[['chrom','start','end','prop_acc','bin','PS00718_H1','PS00718_H2']]
df1 = df1.rename(columns={'PS00718_H1':'H1', 'PS00718_H2':'H2'})
df2 = fire_df[['chrom','start','end','prop_acc','bin','PS00756_H1','PS00756_H2']]
df2 = df2.rename(columns={'PS00756_H1':'H1', 'PS00756_H2':'H2'})
df3 = fire_df[['chrom','start','end','prop_acc','bin','PS00757_H1','PS00757_H2']]
df3 = df3.rename(columns={'PS00757_H1':'H1', 'PS00757_H2':'H2'})
df4 = fire_df[['chrom','start','end','prop_acc','bin','PS00758_H1','PS00758_H2']]
df4 = df4.rename(columns={'PS00758_H1':'H1', 'PS00758_H2':'H2'})

def check_columns(row):
    if row['H1'] == 0 and row['H2'] == 0:
        return('none')
    elif row['H1'] == 1 and row['H2'] == 0:
        return('H1')
    elif row['H1'] == 0 and row['H2'] == 1:
        return('H2')
    else:
        return('both')

filt1 = df1[(df1['H1']>=0) & (df1['H2']>=0)]
filt1['result'] = filt1.apply(check_columns, axis=1)
filt2 = df2[(df2['H1']>=0) & (df2['H2']>=0)]
filt2['result'] = filt2.apply(check_columns, axis=1)
filt3 = df3[(df3['H1']>=0) & (df3['H2']>=0)]
filt3['result'] = filt3.apply(check_columns, axis=1)
filt4 = df4[(df4['H1']>=0) & (df4['H2']>=0)]
filt4['result'] = filt4.apply(check_columns, axis=1)

freq_df = df1[['chrom','start','end','prop_acc','bin']]
freq_df['PS00718'] = filt1['result']
freq_df['PS00756'] = filt2['result']
freq_df['PS00757'] = filt3['result']
freq_df['PS00758'] = filt4['result']
freq_df['neither'] = freq_df.apply(lambda row: (row == 'none').sum(), axis=1)
freq_df['H1_only'] = freq_df.apply(lambda row: (row == 'H1').sum(), axis=1)
freq_df['H2_only'] = freq_df.apply(lambda row: (row == 'H2').sum(), axis=1)
freq_df['both'] = freq_df.apply(lambda row: (row == 'both').sum(), axis=1)

freq_df.to_csv('hap_actuation_bins.csv', index=False)

