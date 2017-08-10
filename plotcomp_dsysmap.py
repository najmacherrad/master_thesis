# DynaMine
# /Memoire/dynamine

# coding=utf-8
from numpy import *
import petl as etl
from re import *
import operator
import glob 
import pandas as pd
import re
import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

#dsysmap
b = pd.read_csv('all_mutations.dat','\t') #28'695 rows x 11 columns
b['Protein_change'] = ''
b['Protein_change'] = b.RES_ORIG + b.RES_NUM.map(str) + b.RES_MUT
b = b.rename(columns = {'UNIPROT_AC':'Uniprot_ACC'})

#DIDA
a = pd.read_csv('didavariantsgenes.csv',',')
d = a.merge(b, on=['Protein_change','Uniprot_ACC'])
d = d.drop(['Protein_change','Uniprot_ACC','GENE_NAME','RES_NUM','RES_ORIG','RES_MUT'],1)
d['ID'].value_counts().to_csv('nb_DIDA.csv')#118 variants
d.to_csv('dsysmap_results_DIDA.csv', index=False) #132 rows x 9 columns

d_fil = d.drop_duplicates('ID', take_last=True)
d_fil.to_csv('dsysmap_results_DIDA_filtered.csv', index=False)

#Neutral
c = pd.read_csv('neutral_variants_final.csv','\t')
c['Protein_change'] = ''
c['Protein_change'] = c.Protein_ref_allele  + c.Protein_position_Biomart.map(str) + c.Protein_alt_allele
c=c.rename(columns = {'transcript_uniprot_id':'Uniprot_ACC','gene_symbol':'Gene_name','dbsnp_id':'ID','snpeff_effect':'Variant_effect'})
e = c.merge(b, on=['Protein_change','Uniprot_ACC'])
e = e.drop(['Protein_position_Biomart','Protein_ref_allele','Protein_alt_allele','GENE_NAME','RES_NUM','RES_ORIG','RES_MUT','id','hgvs_protein','protein_pos','protein_length','gene_ensembl','transcript_ensembl','Uniprot_ACC','Protein_change'],1)
e = e[['ID','Gene_name','Variant_effect','DISEASE_GROUP','DISEASE','MIM','SWISSVAR_ID','STRUCT_CLASS','INTERFACE']]
e['ID'].value_counts().to_csv('nb_1kgp.csv')#241 variants
e.to_csv('dsysmap_results_neutral.csv', index=False) #286 rows x 9 columns

e_fil = e.drop_duplicates('ID', take_last=True)
e_fil.to_csv('dsysmap_results_neutral_filtered.csv', index=False)


####################################################################################
# BARCHARTS en fonction de SURF/BURIED/NA=Not_classified/INTERFACE
####################################################################################
#DIDA -> 241
D = open ('dsysmap_results_DIDA_filtered.csv','r')
c1 = csv.reader(D, delimiter=',')
buried1,surface1,interface1,notclassified1 = 0,0,0,0
for line in c1:
    if line[7] == 'BURIED':
        buried1 = buried1 + 1
    elif line[7]== 'SURF' and line[8]=='NO':
        surface1 = surface1 + 1
    elif line[7]== 'SURF' and line[8]=='YES':
        interface1 = interface1 + 1
    elif line[7]=='': # Not classified
        notclassified1 = notclassified1 + 1

D.close()

#1KGP -> 5'851
N = open ('dsysmap_results_neutral_filtered.csv','r')
c2=csv.reader(N, delimiter=',')
buried2,surface2,interface2,notclassified2 = 0,0,0,0
for line in c2:
    if line[7] == 'BURIED':
        buried2 = buried2 + 1
    elif line[7]== 'SURF' and line[8]=='NO':
        surface2 = surface2 + 1
    elif line[7]== 'SURF' and line[8]=='YES':
        interface2 = interface2 + 1
    elif line[7]=='': # Not classified
        notclassified2 = notclassified2 + 1

N.close()

# BAR PLOT
N = 5
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
DIDA = (float(buried1)/241,float(surface1)/241,float(interface1)/241,float(notclassified1)/241,float(123)/241)
rects1 = ax.bar(ind, DIDA, width, color='red')
neutral = (float(buried2)/5851,float(surface2)/5851,float(interface2)/5851,float(notclassified2)/5851,float(5610)/5851)
rects2 = ax.bar(ind + width, neutral, width, color='blue')
ax.set_ylabel('Frequency')
ax.set_xticks(ind + width)
ax.set_xticklabels(('Buried', 'Surface','Interface','Not classified','NA'))
ax.legend((rects1[0], rects2[0]), ('Deleterious DIDA mutants', 'Neutral 1KGP mutants'),loc='upper left')
#plt.ylim(0,0.35)
fig.savefig('barplot_dsysmap.png')

stats.chi2_contingency(np.column_stack((DIDA,neutral))) #(0.51988356450737627, 0.97153666825577745, 4, array([[ 0.03246648,  0.03233785],[ 0.06887172,  0.06859884],[ 0.02025051,  0.02017028],[ 0.14651059,  0.14593011],[ 0.73605008,  0.73313383]]))
