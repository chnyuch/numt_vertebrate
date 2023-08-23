### Find true numt originated from gene flow
### numt_chk.3.py
### v3 2022/06/30
### include clustering
### numt_chk.2.py
### v2 2022/06/07
### numt_chk.1.py
### v1 2022/06/07
### Find numt originted from gene flow
### @author:Chen Yu Chi

import os,sys
import pandas as pd
import re
import numpy
from collections import defaultdict

# Path and file
# Mammal
dir_pair_out = '/vol/storage/MitoMammal/mito.01.11/dnds_all/pairwise/out_file/'
file_kaks = '/vol/storage/MitoMammal/4NUMTS/run_pairwise/summary_kaks2.csv'
file_best_hit = '/vol/storage/MitoMammal/mito.01.09/best_hits_len.tsv'
file_out = '/vol/storage/MitoMammal/mito.01.11/dnds_all/pairwise/'

# Avian
dir_pair_out = '/vol/storage/MitoAvian/mito.03.14/dnds_all/pairwise/out_file/'
file_kaks = '/vol/storage/MitoAvian/4NUMTS/run_pairwise/summary_kaks2.csv'
file_best_hit = '/vol/storage/MitoAvian/mito.03.11/best_hits_len.tsv'
file_out = '/vol/storage/MitoAvian/mito.03.14/dnds_all/pairwise/'

# def for replace the last nth occurance
def nth_repl(in_str, sub, repl, n):
    s = in_str[::-1]
    find = s.find(sub)
    # If find is not -1 we have found at least one match for the substring
    i = find != -1
    # loop util we find the nth or we find no match
    while find != -1 and i != n:
        # find + 1 means we start searching from after the last match
        find = s.find(sub, find + 1)
        i += 1
    # If i is equal to n we found nth match so replace
    if i == n:
        repl_str = s[:find] + repl + s[find+len(sub):]
        return repl_str[::-1]
    return in_str

# Prepare input
os.chdir(dir_pair_out)

# NUMT pairs form codeml runmode=-2
lst_cds = [cds for cds in os.listdir(dir_pair_out) if cds.endswith('out')]
lst_cds = sorted([cds.strip('.out') for cds in lst_cds])

dict_cds = {}
for cds in lst_cds:
    dict_cds['%s' % cds] = []

for cds in lst_cds:
    with open(dir_pair_out + cds + '.out', 'r') as pair_numt:
        for line in pair_numt:
            #print(line)
            if 'pairwise comparison, codon frequencies:' in line:
                for line in pair_numt:
                    lines = line.rstrip()
                    if lines:
                        dict_cds[cds].append(lines)

df_cds_3 = pd.DataFrame()
for key, lst in dict_cds.items():
    pd.options.mode.chained_assignment = None
    df_cds_1 = pd.DataFrame(numpy.array(lst).reshape(-1, 4)).add_prefix('col')
    df_cds_1['t'] = df_cds_1['col3'].str.replace(r'  S=.*', '').str.replace('t= ', '')
    #print(df_cds_1['t'])
    df_cds_1['S'] = df_cds_1['col3'].str.replace(r'  N=.*', '').str.replace(r't= [^S]*S=', '')
    #print(df_cds_1['S'])
    df_cds_1['dS'] = df_cds_1['col3'].str.replace(r'.*dS =', '')
    #print(df_cds_1[['S', 'dS']])
    df_cds_1['col0'] = df_cds_1['col0'].str.replace('(', '').str.replace(')', '')
    df_cds_1[['numt_1_no', 'numt_1_id', 'delim', 'numt_2_no', 'numt_2_id']] = df_cds_1['col0'].str.split(pat = ' ', expand = True)
    df_cds_2 = df_cds_1[['numt_1_id', 'numt_2_id', 't', 'S', 'dS']]
    df_cds_2['t'] = df_cds_2['t'].apply(pd.to_numeric, errors='coerce')
    df_cds_2 = df_cds_2.loc[df_cds_2['t'] < 0.1].reset_index(drop=True)
    df_cds_2.to_csv(dir_pair_out + key + '_lowt.tsv', sep = '\t', header = True, index = False)
    df_cds_3 = df_cds_3.append(df_cds_2, ignore_index = True)

df_cds_3['numt_1_id'] = df_cds_3['numt_1_id'].apply(lambda x: nth_repl(nth_repl(x, '_', ' ', 3), '_', ':', 2))
df_cds_3['numt_2_id'] = df_cds_3['numt_2_id'].apply(lambda x: nth_repl(nth_repl(x, '_', ' ', 3), '_', ':', 2))
df_cds_3[['species_1', 'numt_1_id']] = df_cds_3['numt_1_id'].str.split(pat = ' ', expand = True)
df_cds_3[['species_2', 'numt_2_id']] = df_cds_3['numt_2_id'].str.split(pat = ' ', expand = True)
# df_cds_3.to_csv(dir_pair_out + 'TEST', sep = '\t', header = True, index = False)

# KaKs between mito cds and numt
kaks = pd.read_csv(file_kaks, sep = "\t") 
info_ks = kaks[['Species', 'fragment_gene', 'ks']]
info_ks['fragment_gene'] = info_ks['fragment_gene'].str.replace(r'  S=.*', '').str.replace('t= ', '')
# merge dataframe
df_cds_ks = pd.merge(df_cds_3, info_ks[['ks', 'fragment_gene']], left_on = 'numt_1_id', right_on = 'fragment_gene', how='left')
df_cds_ks = pd.merge(df_cds_ks, info_ks[['ks', 'fragment_gene']], left_on = 'numt_2_id', right_on = 'fragment_gene', how='left', suffixes=('_1', '_2'))
df_cds_ks = df_cds_ks[['species_1', 'numt_1_id', 'species_2', 'numt_2_id', 'S', 'dS', 't', 'ks_1', 'ks_2']]
df_cds_ks['ks_1'] = df_cds_ks['ks_1'].apply(lambda x: round(x, 9))
df_cds_ks['ks_2'] = df_cds_ks['ks_2'].apply(lambda x: round(x, 9))
# df_cds_ks.to_csv(dir_pair_out + 'TEST', sep = '\t', header = True, index = False)

# Blast result
best_hit = pd.read_csv(file_best_hit, sep = "\t", index_col = 0) 
best_hit['qseqid'] = best_hit['qseqid'].str.replace(' ', '')
# merge dataframe
df_cds_ks_rm = df_cds_ks
df_cds_ks_rm['numt_1_id_mg'] = df_cds_ks_rm['numt_1_id'].str.replace(r'_.*\.fasta.*\.mafft', '')
df_cds_ks_rm['numt_2_id_mg'] = df_cds_ks_rm['numt_2_id'].str.replace(r'_.*\.fasta.*\.mafft', '')
df_cds_all = pd.merge(df_cds_ks_rm, best_hit[['qseqid', 'qspecies', 'sseqid', 'sspecies', 'pident', 'bitscore']], left_on = 'numt_1_id_mg', right_on = 'qseqid', how = 'left')
df_cds_all = pd.merge(df_cds_all, best_hit[['qseqid', 'qspecies', 'sseqid', 'sspecies', 'pident', 'bitscore']], left_on = 'numt_2_id_mg', right_on = 'qseqid', how='left', suffixes=('_1', '_2'))
df_cds_all['pident_1'] = df_cds_all['pident_1'].apply(lambda x: round(x, 9))
df_cds_all['bitscore_1'] = df_cds_all['bitscore_1'].apply(lambda x: round(x, 9))
df_cds_all['pident_2'] = df_cds_all['pident_2'].apply(lambda x: round(x, 9))
df_cds_all['bitscore_2'] = df_cds_all['bitscore_2'].apply(lambda x: round(x, 9))
df_cds_all.drop(['numt_1_id_mg', 'numt_2_id_mg'], axis=1, inplace=True)
df_cds_all.to_csv(file_out + 'gene_flow.1.tsv', sep = ',', header = True, index = False)

