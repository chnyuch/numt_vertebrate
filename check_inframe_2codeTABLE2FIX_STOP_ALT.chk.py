import os, sys
from Bio import SeqIO
from Bio.Seq import Seq
#>>> my_seq = Seq("AGTACACTGGT")

def compare(a1,a2):
	diff=0
	for l in range(len(a1)):
		if a1[l]!=a2[l]: diff+=1
	return diff


def create_diffs(a1,a2, name):
	#open('pairwisePOS/'+name+'.txt')
	diff=0
	##print a1
	##print a2
	for l in range(len(a1)):
		if a1[l]!=a2[l]: 
			diff+=1
			#print name, l+1,  a1[l],  a2[l]
	#return diff

kaksout=open('summary_kaks4.gyn.csv','w')
kaksout.write('Species\tfragment_gene\tka\tks\tkaks\tgene_length\tNo_Substitutions\n')


for u in os.listdir('./'):
	if u.find('.')>-1: continue
	if u=='KaKs_Calculator': continue
	#if u!='seal1': continue
	indir=u+'/alignmentsNEW/'
	alle=os.listdir(indir)
	#record_dict = SeqIO.index("example.fasta", "fasta")
	#print(record_dict["gi:12345678"])  # use any record ID
	liste=[]
	ausgabe=open(u+'/allseqs2.axt','w')
	for a in alle:
			print(a)
			if not a.endswith('mafft'): continue
			fasta = SeqIO.index(indir+a, "fasta")
			#print a
			for r in fasta:
				#print r
				if r.find('-')>-1: 
					seqO= str(fasta[r].seq).lstrip('-').rstrip('-')
					num= seqO.count('-')
					seq1=str(fasta[r].seq)
					#if num==0: print a
				# fragment
				if not r.find('-')>-1: 
					seqO= str(fasta[r].seq).lstrip('-').rstrip('-')
					num1= seqO.count('-')
					#if num==0: print a
					seq2=str(fasta[r].seq)
			if num1==0 and num==0: 
				liste.append(a)
				
				#print seq2
				l1=seq2.rstrip('-').count('-')
				l2=seq2.lstrip('-').count('-')
				
				if l2!=0 and l1!=0:
					s1= seq1[l1:-l2]
					s2= seq2[l1:-l2]
				elif l2==0:
					s1= seq1[l1:]
					s2= seq2[l1:]
				elif l1==0:
					s1= seq1[:-l2]
					s2= seq2[:-l2]
					
				#if s2.find('-')>-1: continue
				# case gene is not completeley covered
				if s1.find('-')>-1: continue
				
				#print l1, l2, a, compare(s1.upper(),s2.upper())
				#print "------------------------"
				#print s1
				#print s2
				#print l1, l2, a, compare(s1.upper(),s2.upper())
				
				# TEST WHETHER CODING!!!!!!!!!!!!!!
				
				if len(s1)%3==1:
					s1=s1[:-2]
					s2=s2[:-2]
				if len(s1)%3==2:
					s1=s1[:-1]
					s2=s2[:-1]
				p1= Seq(s1).translate(table=2,to_stop=True)
				p2= Seq(s2).translate(table=2,to_stop=True)
				
				#print (a,p1,p2)
				if len(str(p1))!=len(str(p2)): continue
				#if len(str(p1))!=len(str(s1))/3: continue
				
				
				# For some of the genes A or AA are added later to form the stop codon
				if len(s1)%3==0:
					ausgabe.write(a+'\n'+s1.upper()+'\n'+s2.upper()+'\n\n')
				elif len(s1)%3==1:
					ausgabe.write(a+'\n'+s1.upper()[:-1]+'\n'+s2.upper()[:-1]+'\n\n')

				elif len(s1)%3==2:
					ausgabe.write(a+'\n'+s1.upper()[:-2]+'\n'+s2.upper()[:-2]+'\n\n')

				p1= Seq(s1).translate(table=2,to_stop=True) # insect 
				p2= Seq(s2).translate(table=2,to_stop=True) # insect
				#p2x= Seq(s2).translate(table=1,to_stop=True)
				#print p1
				#print p2
				#print l1, l2, a, len(s1), len(p1)*3, len(p2)*3, len(p2x)*3,compare(p1.upper(),p2.upper()), compare(s1.upper(),s2.upper())
				#if compare(s1.upper(),s2.upper())> 3:	
				#	print ">"+a.split('.fa_')[1]+' '+a.split('.fa_')[0]+"\n",s1.upper()
				#	print p2x
				#print p2.upper()
				create_diffs(s1.upper(),s2.upper(), a)

	liste.sort()

	#for l in liste:
	#	if l.endswith('mafft'): continue
	#	#print l
	#	os.system('cp '+indir+l+' ADD_analysis/'+l)

	ausgabe.close()
	command=('./KaKs_Calculator -i '+u+'/allseqs2.axt -c 2 -o '+u+'/allseqs4.gyn.csv -m LWL')#	//use Gamma-MYN method
	print(command)
	os.system(command)
	info=open(u+'/allseqs4.gyn.csv').readlines()
	if len(info)<2:continue
	for i in info[1:]:
		zeile=i.split('\t')
		kaksout.write(u+'\t'+zeile[0]+'\t'+zeile[2]+'\t'+zeile[3]+'\t'+zeile[4]+'\t'+zeile[6]+'\t'+zeile[10]+'\n')
		
		
	#print liste

