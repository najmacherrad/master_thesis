# TRANSFORM fasta sequence after mutation
# /Memoire/DIDAfasta

# coding=utf-8
from numpy import *
import petl as etl
import operator
import pandas as pd
import re
import csv

#-------------------------------------------------------------------------------
# Download and modify files from DIDA website
#-------------------------------------------------------------------------------

DIDAgenes = 'DIDA_Genes_a000.csv' #Gene name + Uniprot ACC  (delimiter='\t')

table1 = etl.fromcsv('DIDA_Genes_a000.csv', delimiter='\t')
table2 = etl.rename(table1, {'Uniprot ACC': 'Uniprot_ACC','Gene name': 'Gene_name' })
table3 = etl.cut(table2, 'Uniprot_ACC', 'Gene_name')
etl.totsv(table3,'_didagenes.csv') #137 lignes

table4 = etl.fromcsv('DIDA_Variants_108c0.csv', delimiter='\t')
table5 = etl.split(table4, 'Protein change', '\)', ['Protein change',''])
table6 = etl.cutout(table5,'') 
table7 = etl.split(table6, 'Protein change', '\(', ['','Protein change'])
table8 = etl.cutout(table7,'') 
table9 = etl.rename(table8, {'Gene name': 'Gene_name','Variant effect':'Variant_effect','Protein change':'Protein_change' })
etl.totsv(table9,'_didavariants.csv')

#-------------------------------------------------------------------------------
# Creation of file with fasta sequence mutant
#-------------------------------------------------------------------------------
file_variants = '_didavariants.csv' #364 variants
file_genes = '_didagenes.csv' #136 genes
a = pd.read_csv(file_variants,'\t')
b = pd.read_csv(file_genes,'\t')
merged = a.merge(b, on='Gene_name')
merged.to_csv('_didavariantsgenes.csv', index=False)

V = open('_didavariantsgenes.csv','r')
c1 = csv.reader(V,delimiter=',')
variants = list(c1)
variants.pop(0)

file_fasta_wt = 'DIDAgenes_wt.txt' #136 sequences
file_fasta_mut = 'DIDAgenes_mut.txt' # OUTPUT file
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

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LENGTH OF PROTEINS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fasta = open('DIDAgenes_wt.txt','r')
lines = fasta.readlines()
sequences={}
for i in range(0,len(lines)):
    if lines[i].startswith('>'):
        splitline = lines[i].split('|')
        accessorID = splitline[1]
        sequences[accessorID] = ''
    else:
        sequences[accessorID] = sequences[accessorID] + lines[i].rstrip('\n').rstrip('*')

csv_file = open('DIDAlength.csv','w')
cL = csv.writer(csv_file)
head_row = ['Uniprot_ACC','length']
cL.writerow(head_row)
length_list=[]
for key, value in sequences.items():
    length = len(value)
    length_list.append(len(value))
    cL.writerow([key,length])
csv_file.close()

# BARCHART statistique
import numpy as np
import matplotlib.pyplot as plt
df = pd.read_csv('DIDAlength2.csv','\t')
#df = df.set_index('Uniprot_ACC')
df.plot(kind='bar',legend=False) #,title='Prelude conformation differences for DIDA variants')
plt.ylabel('Uniprot ACC')
plt.savefig('barplotDIDA_length.png')

# HISTOGRAMME
from scipy import stats
from pylab import plot, show, savefig, xlim, figure, \
    hold, ylim, legend, boxplot, setp, axes
fig = figure()
mu2, std2 = stats.norm.fit(log(length_list))
bins = np.linspace(0, 99, 30)
plt.hist(log(length_list),normed=True,bins=15,alpha=0.8)
xmin2, xmax2 = plt.xlim()
x2 = np.linspace(xmin2, xmax2, 100)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x2, p2, 'k--',linewidth=2)
plt.xlabel('log(protein length)')
plt.ylabel('Frequence')
#plt.xlim(-4,4)
#plt.ylim(0,0.8)
plt.title('fit results: mu=%.2f, std=%.2f'%(mu2, std2))
fig.savefig('histo_DIDA_length.png')

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MODIFICATION of fasta sequence after mutation
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PROBLEME =[]
compteur = 0
for v in variants:
    for ID in sequences.keys():              
        if ID == v[4]:
            #---------------------------------------------------------------
            if v[2]=='deletion':
                l1=list(v[3])
                aa1 = l1[0]
                objet1 = re.search(r'[0-9]+',v[3])
                position1 = int(objet1.group())
                s1 = list(sequences[ID])
                if s1[position1 - 1]==aa1 :
                    s1[position1 - 1] = ''
                    seq1=''.join(s1)
                    mutant.write(listID[ID] + '|variant_' + v[0] + '\n')
                    for i in range(0,len(seq1),60):
                        mutant.write(seq1[i:min(len(seq1),i+60)] + '\n')
                else:
                    str1='PROBLEME in '+ ID + ' position ' + str(position1)
                    PROBLEME.append(str1)
                compteur = compteur + 1
            #---------------------------------------------------------------    
            elif v[2]=='insertion':
                l3=list(v[3])
                aa3 = l3[0]
                objet3 = re.search(r'[0-9]+',v[3])
                position3 = int(objet3.group())
                #new AA
                objet3bis = re.search(r'[A-Z]+$',v[3] )
                new_aa3=objet3bis.group()
                s3 = list(sequences[ID])
                if s3[position3 - 1]==aa3 :
                    s3[position3] =  new_aa3 + s3[position3]
                    seq3=''.join(s3)
                    mutant.write(listID[ID] + '|variant_' + v[0] + '\n')
                    for i in range(0,len(seq3),60):
                        mutant.write(seq3[i:min(len(seq3),i+60)] + '\n')
                else:
                    str3 = 'PROBLEME in '+ ID + ' position '+ str(position3)
                    PROBLEME.append(str3)
                compteur = compteur + 1
            #-----------------------------------------------------------------
            elif v[2]=='missense':
                l4=list(v[3])
                aa4 = l4[0]
                objet4 = re.search(r'[0-9]+',v[3])
                position4 = int(objet4.group())
                #new AA
                new_aa4=l4[-1]
                s4 = list(sequences[ID])
                if s4[position4 - 1]==aa4 :
                    s4[position4 - 1] =  new_aa4
                    seq4=''.join(s4)
                    mutant.write(listID[ID] + '|variant_' + v[0] + '\n')
                    for i in range(0,len(seq4),60):
                        mutant.write(seq4[i:min(len(seq4),i+60)] + '\n')
                else:
                    str4 = 'PROBLEME in '+ ID + ' position '+ str(position4)
                    PROBLEME.append(str4)
                compteur = compteur + 1
            #---------------------------------------------------------------NO CHANGE                     
            elif v[2]=='silent': 
                seq6 = sequences[ID] 
                mutant.write(listID[ID] + '|variant_' + v[0] + '\n')
                for i in range(0,len(seq6),60):
                    mutant.write(seq6[i:min(len(seq6),i+60)] + '\n')
                compteur = compteur + 1
            #elif v[2]=='intronic' or v[2]=='splicing' or v[2]=='frameshift' or v[2]=='nonsense': 
                #On les prend pas en compte
                                    
#['PROBLEME in Q14032 position 66', 'PROBLEME in A7E2Y1 position 890', 'PROBLEME in Q9C004 position 241']   
            
fasta.close()
V.close()
mutant.close()     

#TO KNOW NUMBER OF SEQUENCES 
file=open('DIDAgenes_mutALL.txt','r')  
lines = file.readlines()        
liste=[]
for i in range(0,len(lines)):
    if lines[i].startswith('>'):
        liste.append(lines[i])            
file.close()
