import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency


codep_scores = pd.read_csv('codep_scores_pairs.csv')
footprint_counts = pd.read_csv('zmw_footprint_regions.csv')

# count raw footprint occurences per pairwise comparison
neither = []
r1_only = []
r2_only = []
both = []

for i in range(len(codep_scores)):
    pair = codep_scores.iloc[i]
    r1 = str(int(pair['reg1']))
    r2 = str(int(pair['reg2']))
    neither.append(len(footprint_counts[(footprint_counts[r1] == 0) & (footprint_counts[r2] == 0)]))
    r1_only.append(len(footprint_counts[(footprint_counts[r1] == 1) & (footprint_counts[r2] == 0)]))
    r2_only.append(len(footprint_counts[(footprint_counts[r1] == 0) & (footprint_counts[r2] == 1)]))
    both.append(len(footprint_counts[(footprint_counts[r1] == 1) & (footprint_counts[r2] == 1)]))


codep_scores['neither'] = neither
codep_scores['r1_only'] = r1_only
codep_scores['r2_only'] = r2_only
codep_scores['both'] = both
codep_scores['tot_obs'] = codep_scores[['neither','r1_only','r2_only','both']].sum(axis=1)

# Min observations is 21,175
# min(codep_scores['tot_obs'])


# Chi Square Test of independence -------------------------------------
chi2_stats = []
p_values = []
for i in range(len(codep_scores)):
    pair = codep_scores.iloc[i]
    fp_obs = [[pair['neither'],pair['r1_only']],[pair['r2_only'],pair['both']]]
    chi2_stat,p_value,dof, expected = chi2_contingency(fp_obs)
    chi2_stats.append(chi2_stat)
    p_values.append(p_value)

codep_scores['chi2_stat'] = chi2_stats
codep_scores['p_value'] = p_values

# Bonferroni Correction
codep_scores['adj_p_value'] = codep_scores['p_value'] * len(codep_scores)

codep_scores.to_csv("napa_TF_footprinting_codep_chi_square.csv", index=False)
