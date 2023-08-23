import os,sys
import time

## MAKE DATABASE
# make dirs
os.system('mkdir -p Results')
os.system('mkdir -p ResultsClean')

# list mito genome dir
indir1='mitos/'
alle=[mito_genome for mito_genome in os.listdir(indir1) if mito_genome.endswith('.fasta')]

# for species list
noending=[]

# make database for mito genomes
for a in alle:
	os.system('makeblastdb -in '+indir1+a+' -parse_seqids -dbtype nucl')
		## print out command line: makeblastdb -in dir/filename -parse_seqids -dbtype nucl
	noending.append(a.split('.fasta')[0])
		## store species name to the noending list

# list genome dir
indir2='genomes/'
balle=[genome for genome in os.listdir(indir2) if genome.endswith('.fna')]

for a in noending:
	name=a.split('-')[0]
	os.system('mkdir -p '+indir2+a+'/')
	os.system('mv '+indir2+a+'.fna '+indir2+a+'/')

all_sh=[]
for a in balle:
	if a.split('.fna')[0] not in noending: continue
	name=a.split('-')[0]
	#print(name)
	all_sh.append('cd '+indir2+a[:-4]+'/ && '+'makeblastdb -in '+a+' -parse_seqids -dbtype nucl && cd ..')
#print(all_sh)

sh_num=12
split_sh=[all_sh[x:x+sh_num] for x in range(0, len(all_sh), sh_num)]
print(len(split_sh))
##54

pro_num=0
while pro_num < len(split_sh):
	for sh in split_sh[pro_num]:
		os.system(sh+"&")
	time.sleep(90)
	## if running python from outside: nohup python runNUMTidentification2.1.cyc.py &
	## extend the time here cause we need to wait this step finish to start the next
	pro_num+=1

## make database for genomes
#for a in balle:
#	if a.split('.fna')[0] not in noending: continue
#	os.system('makeblastdb -in '+indir2+a+' -parse_seqids -dbtype nucl')
#		## only make database for the genomes with mito genomes
#		## too slow, have to make this parallel

# DO BLAST
# make the full command list
for a in noending:
	name=a.split('-')[0]
	os.system('mkdir -p Results/'+name+'/')

all_sh=[]
for a in alle:
	if not a.endswith('.fasta'): continue
	name=a.split('-')[0]
	#print(name)
	all_sh.append('cd Results/'+name[:-6]+'/ && blastn -db /vol/storage/MitoAvian/4NUMTS/'+indir2+a[:-6]+'/'+a[:-4]+'na -query /vol/storage/MitoAvian/4NUMTS/'+indir1+a+' -outfmt 7 -word_size 20 -num_threads 1 -out '+name+'.blast.out && cd ..')
#print(all_sh)

sh_num=12
split_sh=[all_sh[x:x+sh_num] for x in range(0, len(all_sh), sh_num)]
print(len(split_sh))
##39

pro_num=0
while pro_num < len(split_sh):
	for sh in split_sh[pro_num]:
		os.system(sh+"&")
	time.sleep(10)
	## ask a better way
	pro_num+=1

# CHECK FOR FAILED BLAST
os.chdir('Results/')
os.system("find */*blast.out -empty > empty_files.txt")
os.system("grep -v 'empty_files.txt' empty_files.txt > Failed_BLAST_searches.txt")
os.system('rm empty_files.txt')
os.system('cat Failed_BLAST_searches.txt | while read i; do rm $i; done')
os.system('mv Failed_BLAST_searches.txt ../')

os.chdir('..')

# CLEAN BALST
for a in alle:
	if not a.endswith('.fasta'): continue
	name=a.split('-')[0]
	os.system("grep -v '#' Results/"+name[:-6]+'/'+name+".blast.out | cut -d '\t' -f2,3,4,9,10 > ResultsClean/"+name+".blast.clean.out")

os.chdir('ResultsClean/')
os.system("find *blast.clean.out -empty > empty_files_blastcln.txt")
#os.system("grep -v 'empty_files.txt' empty_files.txt > Species_with_no_numts.txt")
#os.system("rm empty_files.txt")
os.system("cat empty_files_blastcln.txt | while read i; do rm $i; done")
os.system("mv empty_files_blastcln.txt ../")

os.chdir('..')

os.system('ls -1 ResultsClean/*.clean.out > files_list.txt')
os.system('Rscript get_bed2.R')
	# only get NUMTs longer than 200 bp

os.system('mkdir -p Beds')
os.system('mv ResultsClean/*.bed Beds/')

os.system('mkdir -p FASTA_sep')
#os.system('mkdir -p FASTA_merged')

for a in noending:
	name=a.split('-')[0]
	os.system('mkdir -p FASTA_sep/'+name+'/')

#Get fasta
all_sh=[]
for a in alle:
	if not a.endswith('.fasta'): continue
	name=a.split('-')[0]
	#print(name)
	all_sh.append('cd FASTA_sep/'+name[:-6]+'/ && bedtools getfasta -fi /vol/storage/MitoAvian/4NUMTS/'+indir2+a[:-6]+'/'+a[:-4]+'na -bed /vol/storage/MitoAvian/4NUMTS/Beds/'+name+'.blast.clean.out_AllFrags.bed -fo '+name+'.sep.fa && cd ../..')

sh_num=12
split_sh=[all_sh[x:x+sh_num] for x in range(0, len(all_sh), sh_num)]

pro_num=0
while pro_num < len(split_sh):
	for sh in split_sh[pro_num]:
		os.system(sh+"&")
	time.sleep(70)
	## ask a better way
	pro_num+=1

os.system("find FASTA_sep/*/* -empty -print > no_sep_fasta.txt")

#for a in alle:
#	if not a.endswith('.fasta'): continue
#	name=a.split('-')[0]
#	print(name)
#	print('bedtools getfasta -fi '+indir2+a[:-4]+'na -bed Beds/'+name+'.blast.clean.out_AllFrags.bed -fo FASTA_sep/'+name+'.sep.fa')
#	os.system('bedtools getfasta -fi '+indir2+a[:-4]+'na -bed Beds/'+name+'.blast.clean.out_AllFrags.bed -fo FASTA_sep/'+name+'.sep.fa')
#	#os.system('bedtools getfasta -fi '+indir2+a[:-1]+'na -bed Beds/'+name+'.blast.clean.out_AllFrags_merged10k.bed  -fo FASTA_merged/'+name+'.merged.fa')
#	## too slow, have to make this parallel
