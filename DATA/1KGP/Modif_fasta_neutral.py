# TRANSFORM fasta sequence after neutral mutation
# /Memoire/Neutral_mutation

# coding=utf-8
from numpy import *
import petl as etl
import operator
import pandas as pd
import re
import csv

#-------------------------------------------------------------------------------
# Creation of file with fasta sequence mutant
#-------------------------------------------------------------------------------
# dida_vars_1kg.csv = 78'428 lines
#normalement 8'669 unique variants
f1 = csv.reader(open('dida_vars_1kg.csv', 'rb'),delimiter='\t')
writer = csv.writer(open('neutral_variants.csv', 'wb'))
hgvsprot = set()
#SORT BY dbsnp_id
for row in f1:
	if row[33] not in hgvsprot:
		writer.writerow(row)
		hgvsprot.add(row[33]) #8'669 lines

df = pd.read_csv('neutral_variants.csv',',') #176 columns
df2 = df[['id','hgvs_protein','gene_symbol','protein_pos','protein_length','gene_ensembl','transcript_ensembl','transcript_uniprot_id','dbsnp_id','snpeff_effect','Protein_position_Biomart','Protein_ref_allele','Protein_alt_allele']]
df2.to_csv('neutral_variants_modif.csv', index=False)

#-------------------------------------------------------------------------------
V = open('neutral_variants_modif.csv','r')
c1 = csv.reader(V,delimiter=',')
variants = list(c1) #8'670
variants.pop(0)

file_fasta_wt = 'DIDAgenes_wt.txt' #136 sequences
file_fasta_mut = 'DIDAgenes_neutralmut.txt' # OUTPUT file
fasta = open(file_fasta_wt,'r')
mutant = open(file_fasta_mut,'w')
lines = fasta.readlines()

sequences={} #all fasta sequences without the 1st line >
listID={} #to keep the first line of fasta sequence with all informations
for i in range(0,len(lines)):
    if lines[i].startswith('>'):
        splitline = lines[i].split('|') 
        accessorID = splitline[1]
        listID[accessorID] = lines[i].split(' ')[0] 
        sequences[accessorID] = ''
    else:
        sequences[accessorID] = sequences[accessorID] + lines[i].rstrip('\n').rstrip('*')
        
#++++++++++++++++++++++++++++++++++++++++
# MODIFICATION of fasta sequence after mutation
#++++++++++++++++++++++++++++++++++++++++
compteur = 0
pb = open('ProblemsFile.txt','w')
for v in variants:
    for ID in sequences.keys():              
        if ID.split('-')[0] == v[7]:
            #---------------------------------------------------------------
            if v[9]=='deletion' and v[10]!='0':
                aa1 = v[11]
                longueur = len(aa1)
                position1 = int(v[10])
                new_aa1 = v[12]
                s1 = list(sequences[ID])
                if new_aa1=='-0':
                    if ''.join(s1[(position1 - 1):(position1-1+longueur)])==aa1 :
                        s1[(position1 - 1):(position1-1+longueur)] = ''            
                        seq1=''.join(s1)
                        mutant.write(listID[ID] + '|variant_' + v[8] + '\n')
                        for i in range(0,len(seq1),60):
                            mutant.write(seq1[i:min(len(seq1),i+60)] + '\n')
                    else:
                        pb.write('Position problem in '+ ID +' (' + v[8] + ') position ' + str(position1) + ' with effect deletion'+ '\n')
                else:
                    if ''.join(s1[(position1 - 1):(position1-1+longueur)])==aa1 :
                        s1[(position1 - 1):(position1-1+longueur)] = new_aa1            
                        seq1=''.join(s1)
                        mutant.write(listID[ID] + '|variant_' + v[8] + '\n')
                        for i in range(0,len(seq1),60):
                            mutant.write(seq1[i:min(len(seq1),i+60)] + '\n')
                    else:
                        pb.write('Position problem in '+ ID +' (' + v[8] + ') position ' + str(position1) + ' with effect deletion'+ '\n')
            #---------------------------------------------------------------    
            elif v[9]=='insertion' and v[10]!='0':
                aa3 = v[11]
                position3 = int(v[10])
                new_aa3= v[12]
                s3 = list(sequences[ID])
                if s3[position3 - 1]==aa3 :
                    s3[position3-1] =  new_aa3 + s3[position3-1]
                    seq3=''.join(s3)
                    mutant.write(listID[ID] + '|variant_' + v[8] + '\n')
                    for i in range(0,len(seq3),60):
                        mutant.write(seq3[i:min(len(seq3),i+60)] + '\n')
                else:
                    pb.write('Position problem in '+ ID +' (' + v[8] + ') position '+ str(position3) + ' with effect insertion'+ '\n')
            #-----------------------------------------------------------------
            elif v[9]=='missense' and v[10]!='0' and v[12]!='*':
                #-------------
                if v[1].startswith('p.Ter'):
                    position5 = int(v[10])
                    new_aa5=v[12]
                    s5 = list(sequences[ID])
                    s5.append(new_aa5)
                    seq5=''.join(s5)
                    mutant.write(listID[ID] + '|variant_' + v[8] + '\n')
                    for i in range(0,len(seq5),60):
                        mutant.write(seq5[i:min(len(seq5),i+60)] + '\n')
                #-------------
                else:
                    aa4 = v[11]
                    position4 = int(v[10])
                    new_aa4=v[12]
                    s4 = list(sequences[ID])
                    if position4 < len(s4):
                        if s4[position4 - 1]==aa4:
                            s4[position4 - 1] =  new_aa4
                            seq4=''.join(s4)
                            mutant.write(listID[ID] + '|variant_' + v[8] + '\n')
                            for i in range(0,len(seq4),60):
                                mutant.write(seq4[i:min(len(seq4),i+60)] + '\n')
                        else:
                            pb.write('Position problem in '+ ID +' (' + v[8] + ') position '+ str(position4) + ' with effect missense'+ '\n')
                    else:
                        pb.write('Isoform problem in '+ ID +' (' + v[8] + ') position '+ str(position4) + ' with effect missense'+ '\n')
            #---------------------------------------------------------------                     
            elif v[9]=='synonymous' and v[10]!='0' : 
                if v[12]!='*':
                    aa6 = v[11]
                    position6 = int(v[10])
                    s6 = list(sequences[ID])
                    if position6 < len(s6):
                        if s6[position6 - 1]==aa6 :
                            seq6 = sequences[ID] 
                            mutant.write(listID[ID] + '|variant_' + v[8] + '\n')
                            for i in range(0,len(seq6),60):
                                mutant.write(seq6[i:min(len(seq6),i+60)] + '\n')
                            compteur = compteur + 1
                        else:
                            pb.write('Position problem in '+ ID +' (' + v[8] + ') position '+ str(position6) + ' with effect synonymous'+ '\n')               
                    else:
                        pb.write('Isoform problem in '+ ID +' (' + v[8] + ') position '+ str(position6) + ' with effect synonymous'+ '\n')
                elif v[11] == '*' and v[12]=='*':
                    aa6 = v[11]
                    position6 = int(v[10])
                    s6 = list(sequences[ID])
                    if position6 < len(s6):
                        if s6[position6 - 1]==aa6 :
                            seq6 = sequences[ID] 
                            mutant.write(listID[ID] + '|variant_' + v[8] + '\n')
                            for i in range(0,len(seq6),60):
                                mutant.write(seq6[i:min(len(seq6),i+60)] + '\n')
                            compteur = compteur + 1
                        else:
                            pb.write('Position problem in '+ ID +' (' + v[8] + ') position '+ str(position6) + ' with effect synonymous'+ '\n')               
                    else:
                        pb.write('Isoform problem in '+ ID +' (' + v[8] + ') position '+ str(position6) + ' with effect synonymous'+ '\n')
fasta.close()
V.close()
mutant.close()   #6'348 sequences 
pb.close() #2'262 lines

# 49 variants with no protein change -> Protein_position_Biomart = 0
# 11 variants with Protein_alt_allele = *

#TO KNOW NUMBER OF SEQUENCES 
file=open('DIDAgenes_neutralmut.txt','r')  
lines = file.readlines()        
liste=[]
for i in range(0,len(lines)):
    if lines[i].startswith('>'):
        liste.append(lines[i])            
file.close()

file=open('ProblemsFile.txt','r')  
lines = file.readlines()        
liste=[]
for i in range(0,len(lines)):
    liste.append(lines[i])            
file.close()
