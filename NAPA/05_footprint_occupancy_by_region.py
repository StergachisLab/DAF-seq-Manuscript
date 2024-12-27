import pandas as pd


stats = pd.read_csv('zmw_footprint_regions.csv')

stats = stats.loc[:,'1':'11']

df_prop = pd.DataFrame((stats==1).sum(), columns=['n_bound'])
df_prop['n_cov'] = (stats==1).sum() + (stats==0).sum()
df_prop['prop_bound'] = df_prop['n_bound'] / df_prop['n_cov']

df_prop.to_csv('prop_bound_by_region.csv', index_label='region')
