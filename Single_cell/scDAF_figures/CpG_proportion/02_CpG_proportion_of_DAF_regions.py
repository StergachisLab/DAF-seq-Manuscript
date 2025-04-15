import csv

out_rows = []

seqs = []
with open('daf_mappable_regions.fa') as fr:
    for line in fr:
        if not line.startswith('>'):
            seqs.append(line.strip())

n_C = 0
n_CpG = 0
for s in seqs:
    for i in range(len(s)):
        if s[i] == 'C' and i < (len(s)-1):
            if s[i+1] == 'G':
                n_CpG += 1
            else:
                n_C += 1
prop_cpg = n_CpG / (n_CpG+n_C)

out_rows.append(['all','non-CpG',n_C,n_C/(n_CpG+n_C)])
out_rows.append(['all','CpG',n_CpG,n_CpG/(n_CpG+n_C)])


# Stratefy by FIRE peak vs non-FIRE peak
seqs_NO_FIRE = []
with open('daf_mappable_regions_NO_FIRE.fa') as fr:
    for line in fr:
        if not line.startswith('>'):
            seqs_NO_FIRE.append(line.strip())

n_C_NF = 0
n_CpG_NF = 0
for s in seqs_NO_FIRE:
    for i in range(len(s)):
        if s[i] == 'C' and i < (len(s)-1):
            if s[i+1] == 'G':
                n_CpG_NF += 1
            else:
                n_C_NF += 1
prop_cpg_NF = n_CpG_NF / (n_CpG_NF+n_C_NF)

out_rows.append(['non-FIRE','non-CpG',n_C_NF,n_C_NF/(n_CpG_NF+n_C_NF)])
out_rows.append(['non-FIRE','CpG',n_CpG_NF,n_CpG_NF/(n_CpG_NF+n_C_NF)])


seqs_FIRE_ONLY = []
with open('daf_mappable_regions_FIRE_ONLY.fa') as fr:
    for line in fr:
        if not line.startswith('>'):
            seqs_FIRE_ONLY.append(line.strip())

n_C_FIRE = 0
n_CpG_FIRE = 0
for s in seqs_FIRE_ONLY:
    for i in range(len(s)):
        if s[i] == 'C' and i < (len(s)-1):
            if s[i+1] == 'G':
                n_CpG_FIRE += 1
            else:
                n_C_FIRE += 1
prop_cpg_FIRE = n_CpG_FIRE / (n_CpG_FIRE+n_C_FIRE)

out_rows.append(['FIRE','non-CpG',n_C_FIRE,n_C_FIRE/(n_CpG_FIRE+n_C_FIRE)])
out_rows.append(['FIRE','CpG',n_CpG_FIRE,n_CpG_FIRE/(n_CpG_FIRE+n_C_FIRE)])


with open('cpg_prop_of_c_in_DAF_regions.csv','w') as fout:
    header = ['category','CpG','count','prop']
    writer = csv.writer(fout)
    writer.writerow(header)
    for row in out_rows:
        writer.writerow(row)
    
