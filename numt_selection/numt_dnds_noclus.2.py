### calculate dnds of non-cluster, single numt
### numt_dnds_noclus.1.py
### v1 2022/10/31
### calculate dnds of each numt 
### @author:Chen Yu Chi

import os,sys
from Bio.Seq import Seq
from Bio import SeqIO
from ete3 import PhyloTree, Tree
import re
import glob
import pandas as pd

# dir and list
lst_cds_mt = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
# avian
dir_work = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/'
# mammal
dir_work = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/'

def chunkify(total_item, item_in_lst):
    return [total_item[i::item_in_lst] for i in list(range(item_in_lst))]

# NUMT with mt cds
dir_clus = dir_work + 'clusters/'
fa_numt = dir_work + 'clusters/all_cds_numt.fa'
# avian
numt_ds01 = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.04_blastprep/list/numt.ds01.gene.1.lst'
# mammal
numt_ds01 = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.04/list/numt.ds01.gene.1.lst'

df_numt_ds01 = pd.read_csv(numt_ds01, header = None, sep = ' ')
lst_numt_ds01 = df_numt_ds01[1].tolist()

numt_clus = []
for files in os.listdir(dir_clus):
    if not files.endswith('.txt'): continue
    print(files)
    with open(dir_clus + files) as file_clus:
        for line in file_clus:
            id_numt = line.rstrip('\n')
            numt_clus.append(id_numt)

for cds in lst_cds_mt:
    os.system('mkdir -p ' + dir_work + cds + '/mt_tre_numt/sngl_numt/seq/')
    os.system('mkdir -p ' + dir_work + cds + '/mt_tre_numt/sngl_numt/aln/')

for numt_seq in SeqIO.parse(fa_numt, 'fasta'):
    header = numt_seq.name
    if header in numt_clus: continue 
    if header not in lst_numt_ds01: continue
    #print(header)
    mt_gene = re.sub('^.*_' ,'', header)
    mt_gene = re.sub('\..*$', '', mt_gene)
    #print(mt_gene)
    seq_to_file = ('>%s\n%s' % (numt_seq.name, numt_seq.seq + '\n'))
    with open(dir_work + mt_gene + '/mt_tre_numt/sngl_numt/seq/numt_seq.fa', 'a+') as out_fa:
        out_fa.write(seq_to_file)

# 3134 non-cluster numts
# 2244 non-cluster numts with ds > 0.1

# numt aln with cds
lst_sngl_numt = glob.glob(dir_work + '*/mt_tre_numt/sngl_numt/seq/numt_seq.fa')

for file_seq in lst_sngl_numt:
    mt_gene = re.sub('^.*\.01\.12/' ,'', file_seq)
    mt_gene = re.sub('/mt_tre.*$', '', mt_gene)
    i = 0
    for numt_seq in SeqIO.parse(file_seq, 'fasta'):
        #print(numt_seq)
        aa_seq = numt_seq.translate(table = "Vertebrate Mitochondrial")
        aa_seq.name = numt_seq.name
        seq_to_file = ('>%s\n%s' % (aa_seq.name, aa_seq.seq + '\n'))
        # aa files
        with open(dir_work + mt_gene + '/mt_tre_numt/seq/' + mt_gene + '_prank.best.fas', 'r') as in_prank, open(dir_work + mt_gene + '/mt_tre_numt/sngl_numt/aln/aa_numt_' + str(i) + '.fasta' ,'a+') as out_aa:
            out_aa.write(seq_to_file)
            for line in in_prank:
                out_aa.write(line)
        i += 1

# Alignment using prank
# Error message without pre-alignment with mafft
# alignmet with mafft
cmd_no = 100
cmd_exe = 'mafft_exe.sh'
# avian
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/sh/mafft.3/'
# mammal
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/mafft.3/'


lst_sngl_numt = glob.glob(dir_work + '*/mt_tre_numt/sngl_numt/aln/*.fasta')
chunk_lst = list(chunkify(lst_sngl_numt, cmd_no))

os.system('mkdir -p ' + dir_sh)
os.chdir(dir_sh)
i = 0
for index, line in enumerate(chunk_lst):
    with open(dir_sh + 'job_mafft{}.sh'.format(index), 'w') as bash_out:
        bash_out.write('#!/bin/bash\nexport PATH=/prj/ycc-backup/miniconda3/bin:$PATH\nsource activate mafft\n')
        for chunk in chunk_lst[i]:
            dir_fas = re.sub('aa_numt.*$' ,'', chunk)
            to_dir = 'cd ' + dir_fas + '\n'
            mafft = 'ginsi --thread 1 ' + chunk + ' > ' + chunk[:-6] + '_mafft.aln \n'
            bash_out.write(to_dir + mafft)
        i += 1

with open(dir_sh + cmd_exe, 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    #exe_out.write('for job in ' + dir_sh + 'job_iqtre*.sh; do nice -n 1 ${job} & done;')
    for i in range(0, cmd_no):
        exe_out.write('qsub -P fair_share -l idle=1 -l vf=8G -o ' + dir_sh + 'job_mafft{}.out -e '.format(i) + dir_sh + 'job_mafft{}.err -cwd '.format(i) + dir_sh + 'job_mafft{}.sh'.format(i) + '\n')

# nohup bash /prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/mafft.3/mafft_exe.sh

# alignment with prank
cmd_no = 100
cmd_exe = 'prank_exe.sh'
# avian
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/sh/prank.3/'
# mammal
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/prank.3/'

lst_sngl_numt = glob.glob(dir_work + '*/mt_tre_numt/sngl_numt/aln/*_mafft.aln')
chunk_lst = list(chunkify(lst_sngl_numt, cmd_no))

os.system('mkdir -p ' + dir_sh)
os.chdir(dir_sh)
i = 0
for index, line in enumerate(chunk_lst):
    with open(dir_sh + 'job_prank{}.sh'.format(index), 'w') as bash_out:
        bash_out.write('#!/bin/bash\nexport PATH=/prj/ycc-backup/miniconda3/bin:$PATH\nsource activate\n')
        for chunk in chunk_lst[i]:
            dir_fas = re.sub('aa_numt.*$' ,'', chunk)
            to_dir = 'cd ' + dir_fas + '\n'
            prank = 'prank -d=' + chunk[:-9] + 'mafft.aln -njtree -o=' + chunk[:-9] + 'prank.aln -showtree +F \n'
            bash_out.write(to_dir + prank)
        i += 1

with open(dir_sh + cmd_exe, 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    #exe_out.write('for job in ' + dir_sh + 'job_iqtre*.sh; do nice -n 1 ${job} & done;')
    for i in range(0, cmd_no):
        exe_out.write('qsub -P fair_share -l idle=1 -l vf=8G -o ' + dir_sh + 'job_prank{}.out -e '.format(i) + dir_sh + 'job_prank{}.err -cwd '.format(i) + dir_sh + 'job_prank{}.sh'.format(i) + '\n')

os.system('chmod -R 755 ' + dir_work)

# cd /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/sh/prank.3/
# nohup bash /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/sh/prank.3/prank_exe.sh
# cd /prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/ND6/mt_tre_numt/sngl_numt/aln
# for f in *__prank.aln.best.fas; do mv "$f" "$(echo "$f" | sed s/__/_/)"; done
# for f in *__prank.aln.best.dnd; do mv "$f" "$(echo "$f" | sed s/__/_/)"; done

# prepare dna sequences
lst_sngl_numt = glob.glob(dir_work + '*/mt_tre_numt/sngl_numt/seq/numt_seq.fa')

for file_seq in lst_sngl_numt:
    mt_gene = re.sub('^.*\.01\.12/' ,'', file_seq)
    mt_gene = re.sub('/mt_tre.*$', '', mt_gene)
    #print(mt_gene)
    i = 0
    for numt_seq in SeqIO.parse(file_seq, 'fasta'):
        seq_to_file = ('>%s\n%s' % (numt_seq.name, numt_seq.seq + '\n'))
        # dna files
        with open(dir_work + mt_gene + '/mt_tre_numt/seq/' + mt_gene + '.fasta', 'r') as in_prank, open(dir_work + mt_gene + '/mt_tre_numt/sngl_numt/aln/aa_numt_' + str(i) + '_dna.fasta' ,'a+') as out_nt:
            out_nt.write(seq_to_file)
            for line in in_prank:
                out_nt.write(line)
                #print(line)
        i += 1

os.system('sed -i \'s/:/_/g\' ' + dir_work + '*/mt_tre_numt/sngl_numt/aln/aa_numt_*_dna.fasta')

#find */*/clus_aln/*/*.best.fas -type f | wc -l
#find */*/clus_aln/*/tmp* -type d

# alignment aa with trimal and produce codon file for paml
lst_sngl_numt = glob.glob(dir_work + '*/mt_tre_numt/sngl_numt/aln/aa_numt_[0-9]*dna.fasta')

for sngl_numt in lst_sngl_numt:
    file_numt_pre = re.sub('_dna.fasta', '', sngl_numt)
    #print(file_numt_pre)
    numt_no = re.sub('^.*/aln/', '', file_numt_pre)
    #print(numt_no)
    lst_seq = []
    for aa_seq in SeqIO.parse(file_numt_pre + '_prank.aln.best.fas' , 'fasta'):
        lst_seq.append(aa_seq.name)
    with open(file_numt_pre + '_dna_sel.fasta', 'w') as out_nt:
        for nt_seq in SeqIO.parse(file_numt_pre + '_dna.fasta', 'fasta'):
            if nt_seq.name not in lst_seq: continue
            SeqIO.write(nt_seq, out_nt, 'fasta')
    os.system('seqkit sort -n ' + file_numt_pre + '_dna_sel.fasta -o ' + file_numt_pre + '_dna_sel_nm.fasta')
    os.system('seqkit sort -n ' + file_numt_pre + '_prank.aln.best.fas -o ' + file_numt_pre + '_prank.best.nm.fas')
    os.system('pal2nal.pl ' + file_numt_pre + '_prank.best.nm.fas ' +  file_numt_pre + '_dna_sel_nm.fasta -output paml -codontable 2 -nogap 1>' + file_numt_pre + '.phylip 2>>' + dir_work + 'pal2nal.log')
    os.system('pal2nal.pl ' + file_numt_pre + '_prank.best.nm.fas ' +  file_numt_pre + '_dna_sel_nm.fasta -output fasta -codontable 2 -nogap 1>' + file_numt_pre + '_codon.fasta 2>>' + dir_work + 'pal2nal.log')
    os.system('rm ' + file_numt_pre + '_dna_sel_nm.fasta ' + file_numt_pre + '_prank.best.nm.fas')

# check error message in log file
# numt_25 in ND5 too short

# FastTree using DNA sequences
# transfer tree running to Cebitec
# nohup scp -r -P 30183 -i /prj/ycc-backup/denbi.vs.txt ubuntu@129.70.51.6:/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16 . &
cmd_no = 100
cmd_exe = 'fastTree_exe.sh'
# avian
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/sh/fastTree.2/'
# mammal
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/fastTree.2/'

lst_sngl_numt = glob.glob(dir_work + '*/mt_tre_numt/sngl_numt/aln/aa_numt_[0-9]*_codon.fasta')
chunk_lst = list(chunkify(lst_sngl_numt, cmd_no))

os.system('mkdir -p ' + dir_sh)
os.chdir(dir_sh)

i =0
for index, line in enumerate(chunk_lst):
    with open(dir_sh + 'fastTree{}.sh'.format(index), 'w') as bash_out:
        bash_out.write('#!/bin/bash\nexport PATH=/prj/ycc-backup/miniconda3/bin:$PATH\nsource activate\n')
        for seq_codon in chunk_lst[i]:
            dir_numt = re.sub('aa_numt.*', '', seq_codon)
            file_numt_pre = re.sub('_codon.fasta', '', seq_codon)
            #print(dir_numt)
            numt_no = re.sub('^.*/aln/', '', file_numt_pre)
            #print(numt_no)
            to_dir = 'cd ' + dir_numt + '\n'
            fastTree = 'FastTree -gtr -nt ' + seq_codon + ' > ' + dir_numt + numt_no + '_dna_FT.tree \n'
            bash_out.write(to_dir + fastTree)
        i += 1

with open(dir_sh + cmd_exe, 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    for i in range(0, cmd_no):
        exe_out.write('qsub -P fair_share -l idle=1 -l vf=8G -o ' + dir_sh + 'fastTree{}.out -e '.format(i) + dir_sh + 'fastTree{}.err -cwd '.format(i) + dir_sh + 'fastTree{}.sh'.format(i) + '\n')

#os.system('nohup sh ' + dir_sh + cmd_exe)
# cd /prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/fastTree.2/
# nohup bash /prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/fastTree.2/fastTree_exe.sh

# extract subtree with 3 more hierarchical outgroup
# find sister lineages three times
# deroot

len(lst_sngl_numt)
len([item for item in lst_sngl_numt if 'ATP6' not in item])
lst_sngl_numt = glob.glob(dir_work + '*/mt_tre_numt/sngl_numt/aln/aa_numt_[0-9]*_codon.fasta')
#lst_sngl_numt = [item for item in lst_sngl_numt if 'ATP6' not in item]

for trees in lst_sngl_numt:
    dir_numt = re.sub('aa_.*', '', trees)
    os.chdir(dir_numt)
    numt_no = re.sub('^.*/aa_', '', trees)
    numt_no = re.sub('_codon.*$', '', numt_no)
    #print(dir_numt, numt_no)
    try:
        os.system('nw_reroot aa_' + numt_no + '_dna_FT.tree | nw_clade -c 3 -r - \'.*\\.fasta.*mafft\' | nw_reroot - | nw_reroot -d - | cat > aa_' + numt_no + '_paml_FT.treefile')
        tree_all = Tree('aa_' + numt_no + '_paml_FT.treefile')
        tree_all.write(format=8, outfile= numt_no + "_paml_FT.tre")
    except:
        os.system('nw_reroot aa_' + numt_no + '_dna_FT.tree | nw_clade -c 3 -r - \'.*\\.fasta.*mafft\' | nw_reroot - | cat > aa_' + numt_no + '_paml_FT.treefile')
        with open('/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/fail_deroot.lst', 'a+') as out_deroot:
            out_deroot.write(trees + '\n')
        tree_all = Tree('aa_' + numt_no + '_paml_FT.treefile')
        tree_all.write(format=8, outfile= numt_no + "_paml_FT.tre")

os.system('sed -i \'s/NoName//g\' ' + dir_work + '*/*/sngl_numt/*/*_paml_FT.tre')

# sequence files for subtree dNdS calculation
lst_sngl_numt = glob.glob(dir_work + '*/mt_tre_numt/sngl_numt/aln/numt_[0-9]*_paml_FT.tre')

for trees in lst_sngl_numt:
    dir_numt = re.sub('numt_[0-9]*_paml_FT.tre', '', trees)
    #print(dir_numt)
    os.chdir(dir_numt)
    clus_no = re.sub('^.*/aln/', '', trees)
    clus_no = re.sub('_paml_FT.tre$', '', clus_no)
    #print(clus_no)
    tree_all = Tree(trees)
    lst_seq = tree_all.get_leaf_names()
    #print(lst_seq)
    for aa_seq in SeqIO.parse('aa_' + clus_no + '_prank.aln.best.fas', 'fasta'):
        if aa_seq.name not in lst_seq: continue
        #print(clus_no, aa_seq.name)
        with open(clus_no + '_prank_best_sub.fasta', 'a+') as out_aa:
            SeqIO.write(aa_seq, out_aa, 'fasta')
    for nt_seq in SeqIO.parse('aa_' + clus_no + '_dna.fasta', 'fasta'):
        if nt_seq.name not in lst_seq: continue
        #print(clus_no, nt_seq.name)
        with open(clus_no + '_dna_sub.fasta', 'a+') as out_nt:
            SeqIO.write(nt_seq, out_nt, 'fasta')
    os.system('seqkit sort -n ' + clus_no + '_prank_best_sub.fasta -o ' + clus_no + '_prank_best_sub_nm.fas')
    os.system('seqkit sort -n ' + clus_no + '_dna_sub.fasta -o ' + clus_no + '_dna_sub_nm.fasta')
    os.system('pal2nal.pl ' + clus_no + '_prank_best_sub_nm.fas ' + clus_no + '_dna_sub_nm.fasta -output paml -codontable 2 -nogap 1>' + clus_no + '.phylip 2>>' + dir_work + 'pal2nal.log')
    os.system('rm ' + clus_no + '_prank_best_sub_nm.fas ' + clus_no + '_dna_sub_nm.fasta')
    os.chdir(dir_work)

#find /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/*/mt_tre_numt/sngl_numt/aln/*_prank_best_sub_nm.fas -delete
#find /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/*/mt_tre_numt/sngl_numt/aln/*_dna_sub_nm.fasta -delete
os.system('sed -i \' s/mafft/mafft #1/g\' ' + dir_work + '*/*/sngl_numt/*/*_paml_FT.tre')

# paml qsub version
lst_sngl_numt = glob.glob(dir_work + '*/mt_tre_numt/sngl_numt/aln/numt_[0-9]*_paml_FT.tre')
for file_seq in lst_sngl_numt:
    mt_gene = re.sub('^.*\.03\.16/' ,'', file_seq)
    mt_gene = re.sub('/mt_tre.*$', '', mt_gene)
    numt_no = re.sub('^.*aln/', '', file_seq)
    numt_no = re.sub('_paml.*$', '', numt_no)
    #print(numt_no)
    os.system('cp ' + file_seq + ' /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_sngl/tree/' + mt_gene + '_' + numt_no + '_paml.tre')
    os.system('cp ' + dir_work + mt_gene + '/mt_tre_numt/sngl_numt/aln/' + numt_no + '.phylip /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_sngl/seq/' + mt_gene + '_' + numt_no + '.phylip')

# for paml setting file, check mito.03.02_avian.docx
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_sngl/m0/sh/'
cmd_exe = 'm0_exe.sh'
cmd_no = 51
with open(dir_sh + cmd_exe, 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    #exe_out.write('for job in ' + dir_sh + 'job_iqtre*.sh; do nice -n 1 ${job} & done;')
    for i in range(1, cmd_no):
        exe_out.write('qsub -P fair_share -l idle=1 -l vf=8G -o ' + dir_sh + 'job{}.out -e '.format(i) + dir_sh + 'job{}.err -cwd '.format(i) + dir_sh + 'job{}.sh'.format(i) + '\n')

# cd /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_sngl/m0/sh
# nohup bash /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_sngl/m0/sh/m0_exe.sh &

dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_sngl/md2/sh/'
cmd_exe = 'md2_exe.sh'
cmd_no = 51
with open(dir_sh + cmd_exe, 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    #exe_out.write('for job in ' + dir_sh + 'job_iqtre*.sh; do nice -n 1 ${job} & done;')
    for i in range(1, cmd_no):
        exe_out.write('qsub -P fair_share -l idle=1 -l vf=8G -o ' + dir_sh + 'job{}.out -e '.format(i) + dir_sh + 'job{}.err -cwd '.format(i) + dir_sh + 'job{}.sh'.format(i) + '\n')

# cd /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_sngl/md2/sh
# nohup bash /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_sngl/md2/sh/md2_exe.sh &

dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_sngl/md2_fix/sh/'
cmd_exe = 'md2_fix_exe.sh'
cmd_no = 51
with open(dir_sh + cmd_exe, 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    #exe_out.write('for job in ' + dir_sh + 'job_iqtre*.sh; do nice -n 1 ${job} & done;')
    for i in range(1, cmd_no):
        exe_out.write('qsub -P fair_share -l idle=1 -l vf=8G -o ' + dir_sh + 'job{}.out -e '.format(i) + dir_sh + 'job{}.err -cwd '.format(i) + dir_sh + 'job{}.sh'.format(i) + '\n')

# cd /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_sngl/md2_fix/sh
# nohup bash /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_sngl/md2_fix/sh/md2_fix_exe.sh &

# qstat -f | grep chnyuch | wc -l
# grep 'Time used' */*/*/*.out | wc -l

# Extract likelihood
dir_dnds = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_sngl/'
set_dnds = ['m0', 'md2', 'md2_fix']

for dnds in set_dnds:
    os.system('mkdir -p ' + dir_dnds + 'out/' + dnds)
    os.system('ln -snf ' + dir_dnds + dnds + '/out/*/*.out ' + dir_dnds + 'out/' + dnds)

for dnds in set_dnds:
    os.chdir(dir_dnds + 'out/' + dnds)
    for file_out in os.listdir(dir_dnds + 'out/' + dnds):
        #print(file_out)
        with open(file_out, 'r') as in_lnL, open(dir_dnds + 'out/' + dnds + '_lnL.out', 'a+') as out_lnL:
            for line in in_lnL:
                if line.startswith('lnL'):
                    line = re.sub('^.*:', '', line)
                    line = re.sub('\+.*$', '', line)
                    clus_no = re.sub('\.out', '', file_out)
                    #print(clus_no)
                    out_line = clus_no + '\t' + line
                    out_lnL.write(out_line)

df_m0 = pd.read_csv(dir_dnds + 'out/m0_lnL.out', sep = '\t', names = ['clus_no', 'm0'], header = None)
df_md2 = pd.read_csv(dir_dnds + 'out/md2_lnL.out', sep = '\t', names = ['clus_no', 'md2'], header = None)
df_md2_fix = pd.read_csv(dir_dnds + 'out/md2_fix_lnL.out', sep = '\t', names = ['clus_no', 'md2_fix'], header = None)

df_merge_1 = pd.merge(df_m0, df_md2 , on = 'clus_no')
df_merge_2 = pd.merge(df_merge_1, df_md2_fix , on = 'clus_no')
df_merge_2.to_csv(dir_dnds + 'out/lnL.out', sep='\t', index = False)

for file_dnds in os.listdir(dir_dnds + 'out/md2/'):
    with open(dir_dnds + 'out/md2/' + file_dnds, 'r') as in_lnL, open(dir_dnds + 'out/md2_dnds_tree.out', 'a+') as out_lnL:
        #print(file_dnds)
        for line in in_lnL:
            if line.startswith('w ratios as labels for TreeView:'):
                w_tree = next(in_lnL, '').strip()
                #print(w_tree)
                out_line = file_dnds[:-4] + '\t' + w_tree + '\n'
                out_lnL.write(out_line)



