#!/usr/bin/env python
# coding: utf-8

# In[2]:


import argparse

parser = argparse.ArgumentParser(description='The script searches for given T-cell epitopes in given microbial proteomes')

parser.add_argument('-microbiota_names', type=str, help='Path to the file containing names of the .fasta files with microbial proteomes')
parser.add_argument('-microbiota_path', type=str, help='Path to the directory with microbial proteomes')
parser.add_argument('-epitopes', type=str, help='Path to the file with epitopes')
parser.add_argument('-Hamming_distance', type=int, help='Hamming distance for search')
parser.add_argument('-results', type=str, help='Path to the file for resulting table')
args = parser.parse_args()

if args.microbiota_names is None or args.microbiota_path is None or args.epitopes is None or args.Hamming_distance is None or args.results is None:
    print('Missing one or more required arguments.')
else:
    print('microbiota_names:', args.microbiota_names)
    print('microbiota_path:', args.microbiota_path)
    print('epitopes:', args.epitopes)
    print('Hamming_distance:', args.Hamming_distance)
    print('results:', args.results)
    
from Bio import SeqIO
import rapidfuzz
from rapidfuzz.string_metric import hamming
import multiprocessing as mp
import pandas as pd
    
names = str(args.microbiota_names)
microbiota_path = str(args.microbiota_path)
search_epitopes = str(args.epitopes)
hamming = int(args.Hamming_distance)
resulting_table = str(args.results)

## List of proteome identificators
bacteria = []
with open("%s" %names) as f:
    for line in f:
        bacteria.append(line.strip())
        
## List of epitopes
epitopes = []
with open("%s" %search_epitopes) as f:
    for line in f:
        epitopes.append(line.strip())

def aligner(ep):
    '''Return a dataframe with [file, ids[i], epitopes[x], original, # of subs]'''
    global peptides
    global file
    df = rapidfuzz.process.extract(ep, list(peptides.keys()), scorer = hamming, limit = len(peptides), score_cutoff = ham)
    df = pd.DataFrame(df, columns =['original', '# of subs', 'index'])
    df = df[['original', '# of subs']]
    df['Description'] = ep
    df['Gene'] = ''
    df['Species'] = file
    for i in range(len(df)):
        df.iloc[i, 3] = ', '.join(peptides[df.iloc[i, 0]])
    result = df[['Species', 'Gene', 'Description', 'original', '# of subs']]
    return result

def main():
    global results
    pool = mp.Pool(15)
    result = pool.map(aligner, epitopes)
    results = results.append(result)
    pool.close()
    pool.join()
    return results

for file in bacteria:
    results = pd.DataFrame(columns = ['Species', 'Gene', 'Description', 'original', '# of subs'])
    records = []
    ids = []
    ## Path to the file with bacteria proteomes names according to the list of names
    with open("%s%s.fasta" %(microbiota_path, file)) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            records.append(record.seq)
            ids.append(record.id)
    peptides = dict()
    for x in range(len(records)):
        for s in (str(records[x])[i:i+n] for i in range(1 + len(str(records[x]))-n)):
            if len(s) != n:
                print('no!!!!')
            if s not in peptides.keys():
                peptides[s] = []
                peptides[s].append(ids[x])
            else:
                peptides[s].append(ids[x])
    if __name__ == "__main__":
        main()
    print(bacteria.index(file), 'out of', len(bacteria))
    ## Where the resulting table will be stored
    file_results = pd.read_csv('%s' %resulting_table, sep = "\t")
    file_results = pd.concat([file_results, results])
    file_results.to_csv('%s' %resulting_table, sep = '\t', index = None)

