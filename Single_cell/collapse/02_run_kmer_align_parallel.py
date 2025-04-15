import subprocess
import concurrent.futures
import sys


clipped_file = sys.argv[1]

chroms=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']

def launch_kmer_align(file, chrom):
    shell_command = f'python kmer_align.py -r {file} -c {chrom}'
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    return(True)

with concurrent.futures.ThreadPoolExecutor() as executor:
    for c in chroms:
        result = executor.submit(launch_kmer_align, clipped_file, c)

