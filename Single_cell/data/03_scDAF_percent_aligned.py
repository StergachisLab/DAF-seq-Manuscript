import pandas as pd

df = pd.read_csv('percent_mapped.tsv', sep = "\t")

grouped = df[['Sample','Mapped','Unmapped']].groupby('Sample').sum().reset_index()
grouped['Total'] = grouped['Mapped'] + grouped['Unmapped']
grouped['Percent_Mapped'] = (grouped['Mapped']/grouped['Total'])*100

grouped.to_csv('percent_mapped_aggregated.tsv',sep ="\t", index=False)
