from glob import glob
import csv
import gzip


raw_read_files = sorted(glob("*.fasta.gz"))
unique_read_files = sorted(glob("*_pbdup.fa"))

def raw_count(fasta):
    c = 0
    with gzip.open(fasta, 'rt') as fr:
        for line in fr:
            if not line.startswith('>'):
                c += 1
    return(c)

def uniq_count(fasta):
    c = 0
    with open(fasta, 'rt') as fr:
        for line in fr:
            if not line.startswith('>'):
                c += 1
    return(c)

out_data = []
for i in range(len(raw_read_files)):
    sample = raw_read_files[i].split('_')[0]
    rc = raw_count(raw_read_files[i])
    uc = uniq_count(unique_read_files[i])
    dc = rc-uc
    rate = dc/rc
    out_data.append([sample, rc, uc, dc, rate])

with open('duplicate_rates_pbmarkdub.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['sample','raw_reads','unique_reads','duplicate_reads','dup_rate'])
    for row in out_data:
        writer.writerow(row)
