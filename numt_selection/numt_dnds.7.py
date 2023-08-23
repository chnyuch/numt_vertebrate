### calculate dnds of cluster óf numt
### numt_dnds.7.py
### v7 2023/05/09
### Tree marking after unroot, exclude the marking the anstral branch
### numt_dnds.6.py
### v6 2022/01/18
### Change iqtree to FastTree
### v5 2022/12/06
### Mammal
### v4 2022/12/06
### Cebitec
### v3 2022/11/09
### change working serial
### v2 2022/11/07
### modifying pipeline
### v1 2022/10/31
### calculate dnds of each numt cluster 
### @author:Chen Yu Chi

import os,sys
from Bio.Seq import Seq
from Bio import SeqIO
from ete3 import PhyloTree
import re
import glob
from ete3 import Tree
from Bio.Phylo.PhyloXML import Phylogeny
import pandas as pd

def chunkify(total_item, item_in_lst):
    return [total_item[i::item_in_lst] for i in list(range(item_in_lst))]

# dir and list
list_cds_mt = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
# avian
dir_work = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/'
dir_mito = '/prj/ycc-backup/Mito_dnds/MitoAvian/4NUMTS/coding_mito/'
# mammal
dir_work = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/'
dir_mito = '/prj/ycc-backup/Mito_dnds/MitoMammal/4NUMTS/coding_mito/'

# In avian, ND3 in Centrocercus_minimus and ND6 Sterna_hirundo have early stop codon. Both sequences were unannotated and then annotated by mitos.
# Since this might be the sequencing error. There's nothing I can do. Just remove it.   

# Prepare folder and copy sequences
for numt_cds in list_cds_mt: 
    dir_mito_seq = dir_work + numt_cds + '/mt_tre_numt/seq/'
    #os.system('mkdir -p ' + dir_work + numt_cds + '/mt_tre_all/seq/')
    os.system('mkdir -p ' + dir_mito_seq)
    # move mito genes from speceis with both mito genomes and genomes
    for species in os.listdir(dir_mito):
        if not os.path.isdir(dir_mito + '/' + species):continue
        for files in os.listdir(dir_mito + species):
            if files[:-6] == numt_cds: 
                os.system('cp ' + dir_mito + species + '/' + numt_cds + '.fasta ' + dir_mito_seq + species + '.fasta')
        # order sequences by length if multifasta in one file
        os.system('seqkit sort -l -r ' + dir_mito_seq + species + '.fasta -o ' + dir_mito_seq + species + '_sorted.fasta')
        os.system('rm ' + dir_mito_seq + species + '.fasta')

for numt_cds in list_cds_mt:
    dir_mito_seq = dir_work + numt_cds + '/mt_tre_numt/seq/'
    for files in os.listdir(dir_mito_seq):
        if not files.endswith('.fasta'): continue
        species = files[:-13]
        with open(dir_mito_seq + files ,'r') as seq_in, open(dir_mito_seq + species + '_nm.fa', 'w') as seq_out:
            i = 0
            for nt_fa in SeqIO.parse(seq_in, 'fasta'):
                header = str(nt_fa.name)
                seq_to_file = ('>%s\n%s' % (species + '_' + str(i), nt_fa.seq + '\n'))
                seq_out.write(seq_to_file)
                i += 1
    os.system('rm ' + dir_mito_seq + '*.fasta')
    os.chdir(dir_mito_seq)
    os.system('cat *_nm.fa > ' + numt_cds + '.fasta')
    os.system('rm ' + dir_mito_seq + '*_nm.fa')

#for numt_cds in list_cds_mt:
#    dir_mito_seq = dir_work + numt_cds + '/mt_tre_numt/seq/'
#    os.system('grep -c \'_1\' ' + dir_mito_seq + numt_cds + '.fasta') 

# Create aa sequence file
for numt_cds in list_cds_mt:
    dir_mito_seq = dir_work + numt_cds + '/mt_tre_numt/seq/'
    #os.system('seqkit translate -x -T 2 ' + dir_mito_seq + numt_cds + '.fasta > ' + dir_mito_seq + numt_cds + '_aa.fasta')
    #print(dir_mito_seq)
    with open(dir_mito_seq + numt_cds + '_aa.fasta' ,'w') as out_aa:
        for dna_seq in SeqIO.parse(dir_mito_seq + numt_cds + '.fasta', 'fasta'):
            aa_seq = dna_seq.translate(table = "Vertebrate Mitochondrial")
            #print(aa_seq)
            aa_seq.id = dna_seq.name
            aa_seq.description = ""
            SeqIO.write(aa_seq, out_aa, 'fasta')

# Balaeniceps rex atp6 and atp8 reverse 
# need to re-do the analysis      

# Remove short length sequences
for numt_cds in list_cds_mt:
    seq_long = []
    dir_mito_seq = dir_work + numt_cds + '/mt_tre_numt/seq/'
    file_in = open(dir_mito_seq + numt_cds + '_aa.fasta', 'r')
    with open(dir_mito_seq + numt_cds + '_aa_rm.fa', 'w') as file_out:
        if numt_cds == 'ATP6':
            len_cut = 227 # nt 681
        elif numt_cds == 'ATP8':
            len_cut = 56 # nt 168
        elif numt_cds == 'COX1':
            len_cut = 516 # nt 1551
        elif numt_cds == 'COX2':
            len_cut = 227 # nt 684    
        elif numt_cds == 'COX3':
            len_cut = 261 # nt 784
        elif numt_cds == 'CYTB':
            len_cut = 380 # nt 1143
        elif numt_cds == 'ND1':
            len_cut = 325 # nt 978
        elif numt_cds == 'ND2':
            len_cut = 346 # nt 1039
        elif numt_cds == 'ND3':
            len_cut = 116 # nt 351
        elif numt_cds == 'ND4':
            len_cut = 455 # nt 1368
        elif numt_cds == 'ND4L':
            len_cut = 98 # nt 297
        elif numt_cds == 'ND5':
            len_cut = 605 # nt 1818
        elif numt_cds == 'ND6':
            len_cut = 173 # nt 522
        for seqs in SeqIO.parse(file_in, 'fasta'):
            header = str(seqs.name)
            if len(seqs.seq) >= (len_cut*0.8):
                if not "XXXXX" in seqs.seq:
                    seq_to_file = ('>%s\n%s' % (seqs.name, seqs.seq + '\n'))
                    file_out.write(seq_to_file)

# Remove duplicate sequences
for numt_cds in list_cds_mt:
    dir_mito_seq = dir_work + numt_cds + '/mt_tre_numt/seq/'
    os.system('grep \'.*_1\' ' + dir_mito_seq + numt_cds + '_aa_rm.fa >> ' + dir_work + 'repeat_seq.lst')
    with open (dir_mito_seq + numt_cds + '_aa_rm.fa', 'r') as file_in, open(dir_mito_seq + numt_cds + '_aa_rm.1.fa', 'w') as file_out:
        for seqs in SeqIO.parse(file_in, 'fasta'):
            header = str(seqs.name)
            if header.endswith('_0'):
                seq_to_file = ('>%s\n%s' % (seqs.name, seqs.seq + '\n'))
                file_out.write(seq_to_file)
    os.system('rm ' + dir_mito_seq + numt_cds + '_aa_rm.fa')

#for numt_cds in list_cds_mt:
#    dir_mito_seq = dir_work + numt_cds + '/mt_tre_numt/seq/'
#    os.system('grep -c \'>\' ' + dir_mito_seq + numt_cds + '_aa_rm.1.fa') 

# Alignment using mafft -chk
# avian
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/sh/mafft.1/'
# mammal 
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/mafft.1/'
os.system('mkdir -p ' + dir_sh)
os.chdir(dir_sh)

with open(dir_sh + 'job_mafft.sh', 'a+') as bash_out:
    bash_out.write('#!/bin/bash\nexport PATH=/prj/ycc-backup/miniconda3/bin:$PATH\nsource activate mafft\n')
    for numt_cds in list_cds_mt:
        dir_mito_seq = dir_work + numt_cds + '/mt_tre_numt/seq/'
        to_dir = 'cd ' + dir_work + numt_cds + '/mt_tre_numt/seq/' + '\n'
        mafft = 'ginsi --thread 12 ' + dir_mito_seq + numt_cds + '_aa_rm.1.fa > ' + dir_mito_seq + numt_cds + '_aa_rm.mafft.aln \n'
        # recommended for <200 sequences with global homology; https://mafft.cbrc.jp/alignment/server/index.html 
        bash_out.write(to_dir + mafft)

os.system('nohup bash ' + dir_sh + 'job_mafft.sh')
# manually check if there's any large gap

# Alignment using prank
cmd_no = 13
cmd_exe = 'prank_exe.sh'
# avian
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/sh/prank.1/'
# Mammal
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/prank.1/'

chunk_lst = list(chunkify(list_cds_mt, cmd_no))

os.system('mkdir -p ' + dir_sh)
os.chdir(dir_work)
i = 0
for index, line in enumerate(chunk_lst):
    with open(dir_sh + 'job_prank{}.sh'.format(index), 'w') as bash_out:
        bash_out.write('#!/bin/bash\n')
        for chunk in chunk_lst[i]:
            to_dir = 'cd ' + dir_work + chunk + '/mt_tre_numt/seq/' + '\n'
            prank = 'prank -d=' + chunk + '_aa_rm.mafft.aln -njtree -o=' + chunk + '_prank -showtree +F \n'
            #prank = 'prank -d=' + chunk + '_aa_rm.mafft.aln -t=genome_mt_species_namerp_prn.nwk -o=' + chunk + '_prank -F -once \n'
            bash_out.write(to_dir + prank)
        i += 1

with open(dir_sh + cmd_exe, 'w') as exe_out:
     exe_out.write('#!/bin/bash\n')
     exe_out.write('for job in ' + dir_sh + 'job_prank*.sh; do nice -n 1 ${job} & done;')

os.system('chmod -R 755 ' + dir_work)
#os.system('nohup bash ' + dir_sh + cmd_exe)

# Insert numt
# avian
dir_clus = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/clusters/'
fa_numt = '/prj/ycc-backup/Mito_dnds/mito.03.16/clusters/all_cds_numt.fa'
# mammal
dir_clus = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.11/clusters/'
fa_numt = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.11/clusters/all_cds_numt.fa'

for files in os.listdir(dir_clus):
    if not files.endswith('.txt'): continue
    no_clus = files[:-4]
    for numt_seq in SeqIO.parse(fa_numt, 'fasta'):
        header = numt_seq.name
        #print(header)
        with open(dir_clus + files ,'r') as lst_in:
            lst_numt = lst_in.read().splitlines()
            #print(lst_numt)
            mt_gene_1 = re.sub('^.*_' ,'', lst_numt[0])
            mt_gene_2 = re.sub('\..*$', '', mt_gene_1)
            #print(mt_gene_2)   
            for numt in lst_numt:
                if numt in header:
                    seq_to_file = ('>%s\n%s' % (numt_seq.name, numt_seq.seq + '\n'))
                    os.system('mkdir -p ' + dir_work + mt_gene_2 + '/mt_tre_numt/clus_aln/' + no_clus + '/')
                    with open(dir_work + mt_gene_2 + '/mt_tre_numt/clus_aln/' + no_clus + '/' + 'clus_seq.fa', 'a+') as out_fa:
                        out_fa.write(seq_to_file)
                    os.system('cp ' + dir_clus + files + ' ' + dir_work + mt_gene_2 + '/mt_tre_numt/clus_aln/' + no_clus + '/')

lst_clus = glob.glob(dir_work + '/*/mt_tre_numt/clus_aln/*/')

# Realign mt gene alignment with numt
# translate and merge files
for dir_clus in lst_clus:
    #print(dir_clus)
    #mt_gene_1 = re.sub('^.*\.03\.16/' ,'', dir_clus)
    mt_gene_1 = re.sub('^.*\.01\.12/' ,'', dir_clus)
    mt_gene_2 = re.sub('/mt_tre.*$', '', mt_gene_1)
    #print(mt_gene_2)
    #os.system('seqkit translate -x -T 2 ' + dir_clus + 'clus_seq.fa > ' + dir_clus + 'clus_seq_aa.fasta')
    #os.system('cat ' + dir_clus + 'clus_seq_aa.fasta ' + dir_work + mt_gene_2 + '/mt_tre_numt/seq/' + mt_gene_2 + '_prank.best.fas > ' + dir_clus + 'all_aa.fasta')
    with open(dir_work + mt_gene_2 + '/mt_tre_numt/seq/' + mt_gene_2 + '_prank.best.fas', 'r') as in_prank, open(dir_clus + 'all_aa.fasta' ,'a+') as out_aa:
        for dna_seq in SeqIO.parse(dir_clus + 'clus_seq.fa' , 'fasta'):
            aa_seq = dna_seq.translate(table = "Vertebrate Mitochondrial")
            aa_seq.id = dna_seq.name
            aa_seq.description = ""
            SeqIO.write(aa_seq, out_aa, 'fasta')
        for line in in_prank:
            out_aa.write(line)

# Alignment using prank
# Error message without pre-alignment with mafft
# alignmet with mafft
# avian
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/sh/mafft.2/'
# mammal
cmd_no = 50
cmd_exe = 'mafft_exe.sh'
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/mafft.2/'

chunk_lst = list(chunkify(lst_clus, cmd_no))

os.system('mkdir -p ' + dir_sh)
os.chdir(dir_sh)
i = 0
for index, line in enumerate(chunk_lst):
    with open(dir_sh + 'job_mafft{}.sh'.format(index), 'a+') as bash_out:
        bash_out.write('#!/bin/bash\nexport PATH=/prj/ycc-backup/miniconda3/bin:$PATH\nsource activate mafft\n')
        for chunk in chunk_lst[i]:
            to_dir = 'cd ' + chunk + '\n'
            mafft = 'ginsi --thread 12 ' + chunk + 'all_aa.fasta > ' + chunk + 'all.mafft.aln \n'
            bash_out.write(to_dir + mafft)
        i += 1

with open(dir_sh + cmd_exe, 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    for i in range(0, cmd_no):
        exe_out.write('qsub -P fair_share -l idle=1 -l vf=10G -o ' + dir_sh + 'job_mafft{}.out -e '.format(i) + dir_sh + 'job_mafft{}.err -cwd '.format(i) + dir_sh + 'job_mafft{}.sh'.format(i) + '\n')

os.system('chmod -R 755 ' + dir_work)
# cd /prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/mafft.2/
# nohup sh /prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/mafft.2/mafft_exe.sh
#os.system('nohup sh ' + dir_sh + cmd_exe)

# alignment with prank
cmd_no = 50
cmd_exe = 'prank_exe.sh'
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/prank.2/'

chunk_lst = list(chunkify(lst_clus, cmd_no))

os.system('mkdir -p ' + dir_sh)
os.chdir(dir_sh)
i = 0
for index, line in enumerate(chunk_lst):
    with open(dir_sh + 'job_prank{}.sh'.format(index), 'w') as bash_out:
        bash_out.write('#!/bin/bash\nexport PATH=/prj/ycc-backup/miniconda3/bin:$PATH\nsource activate mafft\n')
        for chunk in chunk_lst[i]:
            to_dir = 'cd ' + chunk + '\n'
            prank = 'prank -d=' + chunk + 'all.mafft.aln -njtree -o=' + chunk + 'all.prank -showtree +F \n'
            bash_out.write(to_dir + prank)
        i += 1

with open(dir_sh + cmd_exe, 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    for i in range(0, cmd_no):
        exe_out.write('qsub -P fair_share -l idle=1 -l vf=10G -o ' + dir_sh + 'job_prank{}.out -e '.format(i) + dir_sh + 'job_prank{}.err -cwd '.format(i) + dir_sh + 'job_prank{}.sh'.format(i) + '\n')

os.system('chmod -R 755 ' + dir_work)
#os.system('nohup bash ' + dir_sh + cmd_exe)
# cd /prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/prank.2/
# nohup sh /prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/prank.2/prank_exe.sh

# prepare dna sequences
for dir_clus in lst_clus:
    os.chdir(dir_clus)
    #mt_gene_1 = re.sub('^.*\.03\.16/' ,'', dir_clus)
    mt_gene_1 = re.sub('^.*\.01\.12/' ,'', dir_clus)
    mt_gene_2 = re.sub('/mt_tre.*$', '', mt_gene_1)
    #print(mt_gene_2)
    dir_mito_seq = dir_work + mt_gene_2 + '/mt_tre_numt/seq/'
    with open('clus_seq.fa', 'r') as in_numt,open(dir_mito_seq + mt_gene_2 + '.fasta', 'r') as in_nt, open('all_dna.fasta' ,'a+') as out_nt:
        for line in in_numt:
            out_nt.write(line)
        for line in in_nt:
            out_nt.write(line)

os.system('sed -i \'s/:/_/g\' ' + dir_work + '*/mt_tre_numt/clus_aln/*/all_dna.fasta')

#find */*/clus_aln/*/*.best.fas -type f | wc -l
#find */*/clus_aln/*/tmp* -type d

# alignment aa with trimal and produce codon file for paml
for clus_aln in lst_clus:
    #print(clus_aln)
    os.chdir(clus_aln)
    clus_no = re.sub('^.*/clus_aln/', '', clus_aln)
    clus_no = re.sub('/', '', clus_no)
    lst_seq = []
    for aa_seq in SeqIO.parse('all.prank.best.fas' , 'fasta'):
        lst_seq.append(aa_seq.name)
    #print(lst_seq)
    with open('all_dna_sel.fasta', 'w') as out_nt:
        for nt_seq in SeqIO.parse('all_dna.fasta', 'fasta'):
            if nt_seq.name not in lst_seq: continue
            SeqIO.write(nt_seq, out_nt, 'fasta')
    os.system('seqkit sort -n all_dna_sel.fasta -o all_dna_sel_nm.fasta')
    os.system('seqkit sort -n all.prank.best.fas -o all.prank.best.nm.fas')
    os.system('pal2nal.pl all.prank.best.nm.fas all_dna_sel_nm.fasta -output fasta -codontable 2 -nogap 1> all_codon.fas 2>>' + dir_work + 'pal2nal.log')
    # use codontable 2 but the program still ran with codontable 1
    os.system('rm all_dna_sel.fasta')
    os.chdir(dir_work)

# check error message in log file

# FastTree using DNA sequences
# transfer tree running to Cebitec
# nohup scp -r -P 30183 -i /prj/ycc-backup/denbi.vs.txt ubuntu@129.70.51.6:/vol/storage/MitoAvian/mito.03.16 . &
cmd_no = 50
cmd_exe = 'fastTree_exe.sh'
# avian
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/sh/fastTree.1/'
# mammal
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/sh/fastTree.1/'

chunk_lst = list(chunkify(lst_clus, cmd_no))

os.system('mkdir -p ' + dir_sh)
os.chdir(dir_sh)
i = 0
for index, line in enumerate(chunk_lst):
    with open(dir_sh + 'job_fastTree{}.sh'.format(index), 'w') as bash_out:
        bash_out.write('#!/bin/bash\nexport PATH=/prj/ycc-backup/miniconda3/bin:$PATH\nsource activate\n')
        for chunk in chunk_lst[i]:
            to_dir = 'cd ' + chunk + '\n'
            fastTree = 'FastTree -gtr -nt all_codon.fas > all_codon_FT.tree \n'
            bash_out.write(to_dir + fastTree)
        i += 1

with open(dir_sh + cmd_exe, 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    for i in range(0, cmd_no):
        exe_out.write('qsub -P fair_share -l idle=1 -l vf=8G -o ' + dir_sh + 'job_fastTree{}.out -e '.format(i) + dir_sh + 'job_fastTree{}.err -cwd '.format(i) + dir_sh + 'job_fastTree{}.sh'.format(i) + '\n')

os.system('chmod -R 755 ' + dir_work)
#os.system('nohup bash ' + dir_sh + cmd_exe)
# cd /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/sh/fastTree.1/
# nohup sh /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/sh/fastTree.1/fastTree_exe.sh &
# cd /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16
# find . -mindepth 4 -maxdepth 4 -type d '!' -exec test -e "{}/all_codon.ckp.gz" ';' -print

#os.chdir('/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/ATP8/mt_tre_numt/clus_aln/clus_18')
#os.system('blastn -db /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/ATP8/mt_tre_numt/seq/ATP8_prank_nt -query clus_seq.fa -task megablast -outfmt 7 -word_size 20 -num_threads 1 -out test.out')

# Check if all numt are in one cluster
for trees in lst_clus:
    os.chdir(trees)
    #print(trees)
    clus_no = re.sub('^.*/clus_aln/', '', trees)
    clus_no = re.sub('/', '', clus_no)
    #print(clus_no)
    os.system('nw_reroot all_codon_FT.tree | nw_clade -r - \'.*\\.fasta.*mafft\' | cat > ' + clus_no + '_numt_FT.treefile')
    os.system('nw_reroot all_codon_FT.tree | nw_clade -r - \'.*\\.fasta.*mafft\' | nw_display - > ' + clus_no + '_numt_graph_FT.treefile')
    tree_numt = Tree(clus_no + '_numt_FT.treefile')
    leaf_lst = tree_numt.get_leaf_names()
    #print(species)
    with open(clus_no + '.txt') as in_numt, open(dir_work + 'numt_0out_FT.lst', 'a+') as out0_lst, open(dir_work + 'numt_1out_FT.lst', 'a+') as out1_lst, open(dir_work + 'numt_2out_FT.lst', 'a+') as out2_lst:
        numt_lst = [line.rstrip('\n') for line in in_numt]
        numt_lst = [numt.replace(':', '_') for numt in numt_lst]
        #print(numt)
        i = 0
        for leaf in leaf_lst:
            if leaf in numt_lst: continue
            i += 1
        if i == 0:
            out0_lst.write(clus_no + '\n')
        elif i == 1:
            out1_lst.write(clus_no + '\n')
        else:
            out2_lst.write(clus_no + '\n')

os.system('mkdir -p ' + dir_work + 'numt_tre')
os.chdir(dir_work + 'numt_tre')
os.system('cp ' + dir_work + '*/*/clus_aln/*/*numt_graph_FT.treefile ' + dir_work + 'numt_tre')
os.system('cp ' + dir_work + '*/*/clus_aln/*/*numt_FT.treefile ' + dir_work + 'numt_tre')
# manually checked the graph

# extract subtree with 3 more hierarchical outgroup
# find sister lineages three times
for trees in lst_clus:
    #print(trees)
    os.chdir(trees)
    print(trees)
    os.system('nw_reroot all_codon_FT.tree | nw_clade -c 3 -r - \'.*\\.fasta.*mafft\' | cat > numt_paml_FT.1.tree')
    clus_no = re.sub('^.*/clus_aln/', '', trees)
    clus_no = re.sub('/', '', clus_no)
    tree_all = Tree('numt_paml_FT.1.tree')
    tree_numt = Tree(clus_no + '_numt_FT.treefile')
    leaf_lst = tree_numt.get_leaf_names()
    mrca = tree_all.get_common_ancestor(leaf_lst)
    for node in tree_all.traverse():
        if node.is_leaf(): continue
        if node == mrca:
            node.name = '$1'
        else:
            node.name = ''
    #print(tree_all.write(format=8))
    tree_all.write(format=1, outfile= 'numt_paml_FT.2.tree')
    os.system('gotree unroot -i numt_paml_FT.2.tree -o numt_paml_FT.3.tree')
    os.system('sed -i \'s/:0\.[0-9]*//g; s/:[0-9]*e-[0-9*]//g\' numt_paml_FT.3.tree; rm numt_paml_FT.1.tree numt_paml_FT.2.tree')
    info = open('numt_paml_FT.3.tree').read()
    outwrite = open(clus_no + '_FT_paml.tre', 'w')
    parts=info.split(')$1')
    #print(trees + 'numt_paml_FT.2.tree', len(parts))
    p1=parts[0]
    p1+='$1)' # first branch	
    sec=0
    first=0
    count = -1
    sep = []
    for p in range(len(p1))[::-1]:
        if p1[p]==')':			
            count+=1
            #print(count,')')
        if p1[p]=='(':
            count-=1
			#print(count,'(')
        if count==0 and p1[p]==',': 
            sep.append(p)
        #print(count)
	#print(sec)
    sec=sep[0]
    newp=p1[:sec]+'$1'+p1[sec:]
    outwrite.write(newp+parts[1])
    outwrite.close()
    os.system('rm numt_paml_FT.3.tree ; sed -i \'s/mafft$1/mafft #1/g\' ' + clus_no + '_FT_paml.tre')


# sequence files for subtree dNdS calculation
for clus_aln in lst_clus:
    #print(clus_aln)
    os.chdir(clus_aln)
    clus_no = re.sub('^.*/clus_aln/', '', clus_aln)
    clus_no = re.sub('/', '', clus_no)
    tree_all = Tree('numt_paml_FT.tree')
    lst_seq = tree_all.get_leaf_names()
    #print(lst_seq)
    for aa_seq in SeqIO.parse('all.prank.best.fas', 'fasta'):
        if aa_seq.name not in lst_seq: continue
        #print(clus_no, aa_seq.name)
        with open('all_prank_best_sub_FT.fasta', 'a+') as out_aa:
            SeqIO.write(aa_seq, out_aa, 'fasta')
    for nt_seq in SeqIO.parse('all_dna.fasta', 'fasta'):
        if nt_seq.name not in lst_seq: continue
        #print(clus_no, nt_seq.name)
        with open('all_dna_sub_FT.fasta', 'a+') as out_nt:
            SeqIO.write(nt_seq, out_nt, 'fasta')
    os.system('seqkit sort -n all_prank_best_sub_FT.fasta -o all_prank_best_sub_nm_FT.fas')
    os.system('seqkit sort -n all_dna_sub_FT.fasta -o all_dna_sub_nm_FT.fas')
    os.system('pal2nal.pl all_prank_best_sub_nm_FT.fas all_dna_sub_nm_FT.fas -output paml -codontable 2 -nogap 1>' + clus_no + '_FT.phylip 2>' + dir_work + 'pal2nal.log')
    os.system('rm all_dna_sub_nm_FT.fas all_prank_best_sub_nm_FT.fas numt_paml_FT.tree')
    os.chdir(dir_work)

# paml qsub version
# for paml setting file, check mito.03.02_avian.docx
# avian
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_FT/m0/sh/'
# mammal
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/dnds_FT/m0/sh/'

cmd_exe = 'm0_exe.sh'
cmd_no = 51
with open(dir_sh + cmd_exe, 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    #exe_out.write('for job in ' + dir_sh + 'job_iqtre*.sh; do nice -n 1 ${job} & done;')
    for i in range(1, cmd_no):
        exe_out.write('qsub -P fair_share -l idle=1 -l vf=8G -o ' + dir_sh + 'job{}.out -e '.format(i) + dir_sh + 'job{}.err -cwd '.format(i) + dir_sh + 'job{}.sh'.format(i) + '\n')

# cd /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_FT/m0/sh
# nohup bash /prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/dnds_FT/m0/sh/m0_exe.sh &

# avian
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_FT/md2.1/sh/'
# mammal
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/dnds_FT/md2.1/sh/'
cmd_exe = 'md2_exe.sh'
cmd_no = 51
with open(dir_sh + cmd_exe, 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    #exe_out.write('for job in ' + dir_sh + 'job_iqtre*.sh; do nice -n 1 ${job} & done;')
    for i in range(1, cmd_no):
        exe_out.write('qsub -P fair_share -l idle=1 -l vf=8G -o ' + dir_sh + 'job{}.out -e '.format(i) + dir_sh + 'job{}.err -cwd '.format(i) + dir_sh + 'job{}.sh'.format(i) + '\n')

# cd /prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/dnds_FT/md2.1/sh
# nohup bash /prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/dnds_FT/md2.1/sh/md2_exe.sh &

# avian
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_FT/md2_fix/sh/'
# mammal
dir_sh = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/dnds_FT/md2_fix/sh/'

cmd_exe = 'md2_fix_exe.sh'
cmd_no = 51
with open(dir_sh + cmd_exe, 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    #exe_out.write('for job in ' + dir_sh + 'job_iqtre*.sh; do nice -n 1 ${job} & done;')
    for i in range(1, cmd_no):
        exe_out.write('qsub -P fair_share -l idle=1 -l vf=8G -o ' + dir_sh + 'job{}.out -e '.format(i) + dir_sh + 'job{}.err -cwd '.format(i) + dir_sh + 'job{}.sh'.format(i) + '\n')

# cd /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_FT/md2_fix/sh
# nohup bash /prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/dnds_FT/md2_fix/sh/md2_fix_exe.sh &

dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds/mdC/sh/'
cmd_exe = 'mdC_exe.sh'
cmd_no = 51
with open(dir_sh + cmd_exe, 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    #exe_out.write('for job in ' + dir_sh + 'job_iqtre*.sh; do nice -n 1 ${job} & done;')
    for i in range(1, cmd_no):
        exe_out.write('qsub -P fair_share -l idle=1 -l vf=15G -o ' + dir_sh + 'job{}.out -e '.format(i) + dir_sh + 'job{}.err -cwd '.format(i) + dir_sh + 'job{}.sh'.format(i) + '\n')

# cd /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds/mdC/sh
# nohup bash /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds/mdC/sh/mdC_exe.sh &

dir_sh = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds/m1a/sh/'
cmd_exe = 'm1a_exe.sh'
cmd_no = 51
with open(dir_sh + cmd_exe, 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    #exe_out.write('for job in ' + dir_sh + 'job_iqtre*.sh; do nice -n 1 ${job} & done;')
    for i in range(1, cmd_no):
        exe_out.write('qsub -P fair_share -l idle=1 -l vf=15G -o ' + dir_sh + 'job{}.out -e '.format(i) + dir_sh + 'job{}.err -cwd '.format(i) + dir_sh + 'job{}.sh'.format(i) + '\n')

# cd /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds/m1a/sh
# nohup bash /prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds/m1a/sh/m1a_exe.sh &


# qstat -f | grep chnyuch | wc -l
# grep 'Time used' */*/*/*.out | wc -l

# Extract likelihood
# mammal
dir_dnds = '/prj/ycc-backup/Mito_dnds/MitoMammal/mito.01.12/dnds_FT/'
# avian
dir_dnds = '/prj/ycc-backup/Mito_dnds/MitoAvian/mito.03.16/dnds_FT/'
set_dnds = ['m0', 'md2.1', 'md2_fix']

for dnds in set_dnds:
    os.system('mkdir -p ' + dir_dnds + 'out/' + dnds)
    os.system('ln -snf ' + dir_dnds + dnds + '/out/*/*.out ' + dir_dnds + 'out/' + dnds)

for dnds in set_dnds:
    os.chdir(dir_dnds + 'out/' + dnds)
    for file_out in os.listdir(dir_dnds + 'out/' + dnds):
        print(file_out)
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
df_md2 = pd.read_csv(dir_dnds + 'out/md2.1_lnL.out', sep = '\t', names = ['clus_no', 'md2'], header = None)
df_md2_fix = pd.read_csv(dir_dnds + 'out/md2_fix_lnL.out', sep = '\t', names = ['clus_no', 'md2_fix'], header = None)

df_merge_1 = pd.merge(df_m0, df_md2 , on = 'clus_no')
df_merge_2 = pd.merge(df_merge_1, df_md2_fix , on = 'clus_no')
df_merge_2.to_csv(dir_dnds + 'out/lnL.out', sep='\t', index = False)

for file_dnds in os.listdir(dir_dnds + 'out/md2.1/'):
    with open(dir_dnds + 'out/md2.1/' + file_dnds, 'r') as in_lnL, open(dir_dnds + 'out/md2_dnds_tree.out', 'a+') as out_lnL:
        #print(file_dnds)
        for line in in_lnL:
            if line.startswith('w ratios as labels for TreeView:'):
                w_tree = next(in_lnL, '').strip()
                #print(w_tree)
                out_line = file_dnds[:-4] + '\t' + w_tree + '\n'
                out_lnL.write(out_line)




