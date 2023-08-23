import os, sys
import time

indir='../FASTA_sep/' # NUMT
indir2='../coding_mito/' # coding

alle=os.listdir(indir)
alle2=os.listdir(indir2)

os.system('mv '+indir+'*/* '+indir)
os.system('find '+indir+' -empty -delete')

def run_alignment(species):
	# copy files
	os.system('mkdir '+species)
	os.system('mkdir '+species+'/alignmentsNEW')
	#os.system('mkdir '+a+'/pairs')
	os.system('cp '+indir+'/'+species+'.fasta.sep.fa '+species+'/fragments.fa') 
	os.system('cp '+indir2+species+'/*.fasta '+species+'/') 
	
	os.system('cd '+species+'/')
	k=open(species+'/fragments.fa').read()
	
	# every numt fragment
	for x in k.split('>')[1:]:
		name=x.split('\n')[0]
		# evey gene
		for gene in os.listdir(species):
			if not gene.endswith('.fasta'): continue
			os.system('cat '+species+'/'+gene+' > '+species+'/pair.fa')
			
			# Reverse sequence
			g=open(species+'/pair.fa','a')
			g.write('>'+x)
			g.close()
			
			g=open(species+'/pairTMP.fa','w')
			g.write('>'+x)
			g.close()
			
			os.system('revseq -sequence '+species+'/pairTMP.fa -outseq '+species+'/pairREV.fa')
			os.system('cat '+species+'/'+gene+' >> '+species+'/pairREV.fa')
			# run alignment software
			for k in os.listdir(species+'/alignmentsNEW'):
				if os.stat(species+'/alignmentsNEW/'+k).st_size<1: os.system('rm '+species+'/alignmentsNEW/'+k)
			if name+'_'+gene+'.mafft' not in os.listdir(species+'/alignmentsNEW/'): 
				os.system('mafft --thread 8 '+species+'/pair.fa > '+species+'/alignmentsNEW/'+name+'_'+gene+'.mafft')
			if name+'_'+gene+'.REV.mafft' not in os.listdir(species+'/alignmentsNEW/'): 
				os.system('mafft --thread 8 '+species+'/pairREV.fa > '+species+'/alignmentsNEW/'+name+'_'+gene+'.REV.mafft')
	os.system('cd ..')

sh_num=12
all_sh=[]
for a in alle2:
	if a+'.fasta.sep.fa' not in os.listdir(indir): continue
	all_sh.append(a)

split_sh=[all_sh[x:x+sh_num] for x in range(0, len(all_sh), sh_num)]

from multiprocessing import Pool
pro_num=0
while pro_num < len(split_sh):
	pool = Pool()
	pool.map(run_alignment, split_sh[pro_num])
	## ask a better way
	pro_num+=1

"""
pro_num=0
while pro_num < len(split_sh):
	for a in split_sh[pro_num]:
		print(a)
		if a+'.fasta.sep.fa' not in os.listdir(indir):
			print('nohit')
			continue
		run_alignment(a)
	## ask a better way
	pro_num+=1
"""
"""
#####original_code#####

for a in alle2:
	#for a2 in alle2:
	print(a)
	if a+'.fasta.sep.fa' not in os.listdir(indir):
		print('nohit') 
		continue
	#if a=='ant1' :continue
	
	
	#os.system('needle -asequence '+indir+a+' -bsequence mitogenome/sequence.fasta -gapopen 10 -gapextend 0.5 -outfile pairs/'+a+'.needle')
	#os.system('cat '+indir2+a2+' '+indir+a+' > pair.fa')
	
	# copy files
	os.system('mkdir '+a)
	os.system('mkdir '+a+'/alignmentsNEW')
	#os.system('mkdir '+a+'/pairs')
	os.system('cp '+indir+'/'+a+'.fasta.sep.fa '+a+'/fragments.fa') 
	os.system('cp '+indir2+a+'/*.fasta '+a+'/') 
	
	
	k=open(a+'/fragments.fa').read()
	
	# every numt fragment
	for x in k.split('>')[1:]:
		name=x.split('\n')[0]
		# evey gene
		for gene in os.listdir(a):
			if not gene.endswith('.fasta'): continue
			os.system('cat '+a+'/'+gene+' > '+a+'/pair.fa')
			
			# Reverse sequence
			
			g=open(a+'/pair.fa','a')
			g.write('>'+x)
			g.close()
			
			g=open(a+'/pairTMP.fa','w')
			g.write('>'+x)
			g.close()
			
			os.system('revseq -sequence '+a+'/pairTMP.fa -outseq '+a+'/pairREV.fa')
			os.system('cat '+a+'/'+gene+' >> '+a+'/pairREV.fa')
			# run alignment software
			for k in os.listdir(a+'/alignmentsNEW'):
				if os.stat(a+'/alignmentsNEW/'+k).st_size<1: os.system('rm '+a+'/alignmentsNEW/'+k)
			if name+'_'+gene+'.mafft' not in os.listdir(a+'/alignmentsNEW/'): 
				#os.system('cp pair.fa pairs/'+name+'_'+gene+'.pair')
				os.system('mafft --thread 8 '+a+'/pair.fa > '+a+'/alignmentsNEW/'+name+'_'+gene+'.mafft')
			if name+'_'+gene+'.REV.mafft' not in os.listdir(a+'/alignmentsNEW/'): 
				#os.system('cp pair.fa pairs/'+name+'_'+gene+'.pair')
				os.system('mafft --thread 8 '+a+'/pairREV.fa > '+a+'/alignmentsNEW/'+name+'_'+gene+'.REV.mafft')
			#os.system('prank '+a+'/pair.fa')
			#os.system('mv output.best.fas '+a+'/alignments/'+name+'_'+gene+'.prank')
				
	#
	#os.system('mafft --thread 2 pair.fa >  gene_pairwise/'+a2+'_'+a+'.mafft')
	#os.system('prank pair.fa')
	#os.system('mv output.best.fas gene_pairwise/'+a2+'_'+a+'.prank')
"""